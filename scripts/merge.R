#Script to merge adjoining neighbouring high-het segments if they are separated by 1 marker

names = c("AC","BN","BU","F3","M5","MR","WK","WN")
strains = names

POSITIONS = list()
MAP = list()
for(chr in 1:21){
	positions_in_raw = fread(paste("../chr",CHR,"raw.raw",sep=""),nrow=1,stringsAsFactors=FALSE,data.table=FALSE,header=FALSE)
	positions_in_raw = positions_in_raw[,-c(1:6)]
	positions_in_raw = positions_in_raw[1,]
	positions_in_raw = gsub("_.*","",positions_in_raw)
	
	toremove_raw = c()
	for(i in length(positions_in_raw):2){if(positions_in_raw[i] == positions_in_raw[i-1]){toremove_raw = c(toremove_raw, i);}}
	if(length(toremove_raw) > 0){
		positions_in_raw = positions_in_raw[-toremove_raw]
	}

	map = fread(paste0("../chr",chr,".map"),stringsAsFactors=FALSE,data.table=FALSE)
	map = map[!duplicated(map[,4]),]
	map = merge(positions_in_raw, map, by.x=1,by.y=2,sort=FALSE)

	POSITIONS[[chr]] = positions_in_raw
	MAP[[chr]] = map
}

#library(readr)
library(data.table)

strain_num = 0
for(strain in strains){
	print(strain)
	strain_num = strain_num + 1
	towrite_bed = data.frame()
	bed = fread(paste0("../RETRY_strain_",strain,"_0.25_highhetsegments_withX.bed"),stringsAsFactors=FALSE,data.table=FALSE)
	for(chr in 1:21){
		print(chr)
		CHR = chr
		
		map = MAP[[chr]]
		
		bedchr = bed[which(bed[,1] == chr),]
		for(i in nrow(bedchr):2){
			s_i = bedchr[i,2]
			e_iminus1 = bedchr[i-1,3]

			p_s = which(map[,4] == e_iminus1)
			p_e = which(map[,4] == s_i)

			w_s = (p_s) / 1000
			w_e = (p_e-1)/1000 + 1

			print(w_e - w_s)
			if(w_e - w_s == 2){
				if(HETS[[CHR]][strain_num,w_s+1] > 175){
					#merge
					bedchr[i-1,3] = bedchr[i,3]
					bedchr = bedchr[-i,]
				}
			}
			
		}
		towrite_bed = rbind(towrite_bed, bedchr)
	}
	write.table(towrite_bed, file=paste0("strain_",strain,"_merged.bed"),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
}
