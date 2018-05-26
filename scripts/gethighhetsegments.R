names = c("AC","BN","BU","F3","M5","MR","WK","WN")
strains = names

#library(readr)
library(data.table)
HETS = list()

#

#read in the genotype files, and calculate the number of hets per 1000-marker window in matrix N1s, and HETS
NON_MISSING = list()
for(CHR in 21:1){
	print(CHR)
	F = fread(paste("/scratch/junzli_flux/sramdas/haplotyping/founders/rn6/unifiedgenotyper/chr",CHR,"raw.raw",sep=""),skip=1,stringsAsFactors=FALSE)
    F = as.data.frame(F[,-c(1:6)])
    N1s = matrix(0, nrow=8, ncol=round(ncol(F)/1000))
	nonmissing = matrix(1000, nrow=8, ncol = round(ncol(F)/1000))
	
    for(i in 0:(ncol(N1s)-2)){
            s = (i*1000)+1
            e = (i*1000)+1000
            ns = apply(F[,s:e], 1, function(x){length(which(x==1))})
            N1s[,i+1] = ns
			nonmissing[,i+1] = apply(F[,s:e],1, function(x){1000 - length(which(is.na(x)));})
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

#strain to calculate hets for
for(j in 1:8){
	towrite_bed = data.frame()
	print(j)
	for(CHR in 1:21){
		print(CHR)
		chr = CHR
		S = -9
		E = -9

		F = fread(paste("/scratch/junzli_flux/sramdas/haplotyping/founders/rn6/unifiedgenotyper/chr",CHR,"raw.raw",sep=""),skip=1,stringsAsFactors=FALSE)
		F = as.data.frame(F[,-c(1:6)])		
		
		positions_in_raw = fread(paste("chr",CHR,"raw.raw",sep=""),nrow=1,stringsAsFactors=FALSE,data.table=FALSE,header=FALSE)
		positions_in_raw = positions_in_raw[,-c(1:6)]
		positions_in_raw = positions_in_raw[1,]
		positions_in_raw = gsub("_.*","",positions_in_raw)
		
		toremove_raw = c()
		for(i in length(positions_in_raw):2){if(positions_in_raw[i] == positions_in_raw[i-1]){toremove_raw = c(toremove_raw, i);}}
		if(length(toremove_raw) > 0){
			positions_in_raw = positions_in_raw[-toremove_raw]
			F = F[,-toremove_raw]

		}

		map = fread(paste0("chr",chr,".map"),stringsAsFactors=FALSE,data.table=FALSE)
		map = map[!duplicated(map[,4]),]
		map = merge(positions_in_raw, map, by.x=1,by.y=2,sort=FALSE)
			
		N1s = matrix(0, nrow=8, ncol=round(ncol(F)/1000))
		nonmissing = matrix(1000, nrow=8, ncol = round(ncol(F)/1000))

		for(i in 0:(ncol(N1s)-2)){
				s = (i*1000)+1
				e = (i*1000)+1000
				ns = apply(F[,s:e], 1, function(x){length(which(x==1))})
				N1s[,i+1] = ns
				nonmissing[,i+1] = apply(F[,s:e],1, function(x){1000 - length(which(is.na(x)));})
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


		segmentstarts = c()
		segmentends = c()
		s_lengths = c()
		N1s = HETS[[CHR]]
		
		starts = c()
		ends = c()
		indices_s = c()
		indices_e = c()

		for(i in 1:(ncol(N1s)-1)){
			if(N1s[j,i]/NON_MISSING[[CHR]][j,i] > 0.25){
				if(S == -9){
					S = i
				}else{
					E = i
				}
			}
			if(N1s[j,i]/NON_MISSING[[CHR]][j,i] < 0.25){
#				print(i)
#				print(E)
				if(S != -9){
					E = i - 1
#					if(E == -9){
#						E = S
#					}
					starts = c(starts, S)
					ends = c(ends, E-1)
					indices_s = c(indices_s, ((S-1)*1000)+1)
					indices_e = c(indices_e, ((E-1)*1000)+1000)
					S = -9
					E = -9
				}else{
					donothing = 1
				}
#				break

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

	write.table(towrite_bed, file=paste("strain_",strains[j],"_0.25_highhetsegments_withX.bed",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
}

