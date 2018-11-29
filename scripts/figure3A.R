library(data.table)
library(ggplot2)
library(gridExtra)

library(data.table)
library(reshape2)
library(ggplot2)
library(stringr)

names = c("AC","BN","BU","F3","M5","MR","WK","WN")
#geneannotations = read.table("geneannotation_rn6.txt",skip=1,stringsAsFactors=F,sep=" ")

num_windows = matrix(0, 8, 21)

#the read depth file is available from: https://www.dropbox.com/s/lqnsyjfaoqny84e/readdepths.txt.gz?dl=0
depths = fread("readdepths.txt",stringsAsFactors=FALSE,sep="\t",data.table=FALSE)
#depths[,1] = as.numeric(gsub("chr","",depths[,1]))
HETs = list()

N1S = list()
DEP = list()


for(chr in 20:1){
	CHR = chr
	toplot = c()
	print(chr)

	if(CHR > 1){
		raw = fread(paste("chr",chr,"raw.raw",sep=""),stringsAsFactors=FALSE,skip=1)
	}else{
		raw = fread('temp',stringsAsFactors=FALSE,data.table=FALSE,header=FALSE)
	}
	raw = as.data.frame(raw[,-c(1:6)])
	raw = as.matrix(raw)
	
	map = fread(paste("chr",chr,".map",sep=""),stringsAsFactors=FALSE, data.table=FALSE)
	if(CHR > 1){
		positions_in_raw = fread(paste("chr",CHR,"raw.raw",sep=""),nrow=1,stringsAsFactors=FALSE,data.table=FALSE,header=FALSE)
	}else{
		positions_in_raw = fread("tempheader",nrow=1,stringsAsFactors=FALSE,data.table=FALSE,header=FALSE)
	}
	positions_in_raw = positions_in_raw[,-c(1:6)]
	positions_in_raw = positions_in_raw[1,]
	positions_in_raw = gsub("_.*","",positions_in_raw)
	map = merge(positions_in_raw, map, by.x=1,by.y=2,sort=FALSE)		

	chrdepth =depths[which(depths[,1] == paste0("chr",chr)),]
	depthschr = merge(map, chrdepth, by.x=4,by.y=2,sort=FALSE)
	
#	par(mfrow=c(2,1))
	
	N1s = matrix(0, nrow=8, ncol=round(ncol(raw)/1000))
	dep = c()
	GCs = c()
	positions = c()
	
	for(i in 0:(ncol(N1s)-2)){
		s = (i*1000)+1
		e = (i*1000)+1000
		ns = apply(raw[,s:e], 1, function(x){length(which(x==1))})
		N1s[,i+1] = ns
		positions = c(positions, map[s,4])
		dep = c(dep,mean(as.numeric(depthschr[s:e,11]),))
   }
	
	DEP[[CHR]] = dep
	N1S[[CHR]] = N1s
	
}

num = c()
d = c()

for(chr in 1:20){
	for(i in 1:ncol(N1S[[chr]])){
		num = c(num, length(which(N1S[[chr]][,i] > 250)))
		d = c(d, DEP[[chr]][i])
	}
}



ggplot(toplot, aes(x=as.factor(num),y=d)) + geom_boxplot() + geom_point(aes(size=0.75)) + labs(x="Num. founders with heterozygous calls", y="Read Depth")
ggsave("figure3a.pdf",width=3.5,height=3.5,device="pdf")
