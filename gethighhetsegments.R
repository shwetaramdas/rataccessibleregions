########################
26 March 2018
Author: Shweta Ramdas
From a raw file, calculates the heterozygosity in windows across the genome (part 1), and then combines adjacent high-het windows into segments (part 2)
Length of each window set by NUM_MARKERS
Threshold to define 'high-heterizygosity': PERCENT_HET_THRESHOLD

########################

library(data.table)

NUM_MARKERS = 1000
PERCENT_HET_THRESHOLD = 0.1

#read in the genotype files, and calculate the number of hets per 1000-marker window in matrix N1s, and HETS
NON_MISSING = list()
for(CHR in 20:1){
	print(CHR)
	F = fread(paste("/scratch/junzli_flux/sramdas/haplotyping/founders/rn6/chr",CHR,"raw.raw",sep=""),skip=1,stringsAsFactors=FALSE)
	F = as.data.frame(F[,-c(1:6)])
	N1s = matrix(0, nrow=8, ncol=round(ncol(F)/NUM_MARKERS))
	nonmissing = matrix(NUM_MARKERS, nrow=8, ncol = round(ncol(F)/NUM_MARKERS))
	
	for(i in 0:(ncol(N1s)-2)){
		s = (i*NUM_MARKERS)+1
		e = (i*NUM_MARKERS)+NUM_MARKERS
		ns = apply(F[,s:e], 1, function(x){length(which(x==1))})
		N1s[,i+1] = ns
		nonmissing[,i+1] = apply(F[,s:e],1, function(x){NUM_MARKERS - length(which(is.na(x)));})
	}
	nonmissing[1,which(nonmissing[1,]==0)] = 1
	nonmissing[2,which(nonmissing[2,]==0)] = 1
	nonmissing[3,which(nonmissing[3,]==0)] = 1
	nonmissing[4,which(nonmissing[4,]==0)] = 1
	nonmissing[5,which(nonmissing[5,]==0)] = 1
	nonmissing[6,which(nonmissing[6,]==0)] = 1
	nonmissing[7,which(nonmissing[7,]==0)] = 1
	nonmissing[8,which(nonmissing[8,]==0)] = 1
	
	NON_MISSING[[CHR]] = nonmissing
	HETS[[CHR]] = N1s
}

SEGMENTS = list()
BED = list()

SEGMENT_starts = list()
SEGMENT_ends = list()
segmentlengths = list()
segmentlengths[[1]] = c(NA)
segmentlengths[[2]] = c(NA)
segmentlengths[[3]] = c(NA)
segmentlengths[[4]] = c(NA)
segmentlengths[[5]] = c(NA)
segmentlengths[[6]] = c(NA)
segmentlengths[[7]] = c(NA)
segmentlengths[[8]] = c(NA)

#strain to calculate hets for (1 to 8)
for(j in 1:8){
	towrite_bed = data.frame()
	print(j)
	for(CHR in 1:20){
		print(CHR)
		chr = CHR
		S = -9
		E = -9
		map = fread(paste("/scratch/junzli_flux/sramdas/haplotyping/founders/rn6/chr",chr,".map",sep=""),stringsAsFactors=FALSE, data.table=FALSE)
		positions_in_raw = fread(paste("/scratch/junzli_flux/sramdas/haplotyping/founders/rn6/chr",CHR,"raw.raw",sep=""),nrow=1,stringsAsFactors=FALSE,data.table=FALSE,header=FALSE)
		positions_in_raw = positions_in_raw[,-c(1:6)]
		positions_in_raw = positions_in_raw[1,]
		positions_in_raw = gsub("_.*","",positions_in_raw)
		map = merge(positions_in_raw, map, by.x=1,by.y=2,sort=FALSE)
		
		segmentstarts = c()
		segmentends = c()
		s_lengths = c()
		N1s = HETS[[CHR]]
		
		starts = c()
		ends = c()
		indices_s = c()
		indices_e = c()

		for(i in 1:(ncol(N1s)-1)){
	#        print(i)
			if(N1s[j,i]/NON_MISSING[[CHR]][j,i] > PERCENT_HET_THRESHOLD){
				if(S == -9){
					S = i
				}else{
					E = i
				}
			}
			if(N1s[j,i]/NON_MISSING[[CHR]][j,i] < PERCENT_HET_THRESHOLD){
				if(S != -9){
					E = i
					starts = c(starts, S)
					ends = c(ends, E-1)
					indices_s = c(indices_s, (S*1000)+1)
					indices_e = c(indices_e, (E*1000)+1000)
					S = -9
					E = -9
				}else{
					donothing = 1
				}
			}
		}
		positions_s = map[indices_s,4]
		positions_e = map[indices_e,4]
		if(indices_e[length(indices_e)] > nrow(map)){
			positions_e[length(positions_e)] = map[nrow(map),4]
		}
		segmentlengths[[j]] = c(segmentlengths[[j]], positions_e - positions_s)

		segmentstarts = c(segmentstarts, starts)
		segmentends = c(segmentends, ends)

		SEGMENT_starts[[chr]] = positions_s
		SEGMENT_ends[[chr]] = positions_e
		towrite_bed = rbind(towrite_bed, cbind(rep(CHR, length(positions_s)), positions_s, positions_e));
	}
	BED[[j]] = towrite_bed
	SEGMENTS[[j]] = list()
	SEGMENTS[[j]][['starts']] = SEGMENT_starts
	SEGMENTS[[j]][['ends']] = SEGMENT_ends
	write.table(towrite_bed, file=paste("strain_",strains[j],"highhetsegments.bed",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
}

