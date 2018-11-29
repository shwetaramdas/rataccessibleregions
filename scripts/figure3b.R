library(data.table)
library(ggplot2)
library(gridExtra)

library(data.table)
library(reshape2)
library(ggplot2)
library(stringr)


names = c("AC","BN","BU","F3","M5","MR","WK","WN")
strains = names

num_windows = matrix(0, 8, 20)
depths = fread("readdepths.txt",stringsAsFactors=FALSE,sep="\t",data.table=FALSE)
HETs = list()

chr=1
CHR = chr
toplot = c()
print(chr)

raw = fread(paste("chr",chr,"raw.raw",sep=""),stringsAsFactors=FALSE,skip=1)
raw = as.data.frame(raw[,-c(1:6)])
raw = as.matrix(raw)

map = fread(paste("chr",chr,".map",sep=""),stringsAsFactors=FALSE, data.table=FALSE)
positions_in_raw = fread(paste("chr",CHR,"raw.raw",sep=""),nrow=1,stringsAsFactors=FALSE,data.table=FALSE,header=FALSE)
positions_in_raw = positions_in_raw[,-c(1:6)]
positions_in_raw = positions_in_raw[1,]
positions_in_raw = gsub("_.*","",positions_in_raw)

#remove duplicate positions in raw file
toremove = c()
for(i in length(positions_in_raw):2){
	if(positions_in_raw[i] == positions_in_raw[i-1]){toremove = c(toremove, i); positions_in_raw = positions_in_raw[-i];}
}
raw = raw[,-toremove]

map = unique(map)
map = merge(positions_in_raw, map, by.x=1,by.y=2,sort=FALSE)		
chrdepth =depths[which(depths[,1] == paste0("chr",chr)),]
chrdepth = chrdepth[!duplicated(chrdepth[,2]),]

depthschr = merge(map, chrdepth, by.x=4,by.y=2,sort=FALSE)

N1s = matrix(0, nrow=8, ncol=round(ncol(raw)/1000))
positions = c()
for(i in 0:(ncol(N1s)-2)){
	s = (i*1000)+1
	e = (i*1000)+1000
	ns = apply(raw[,s:e], 1, function(x){length(which(x==1))})
	N1s[,i+1] = ns
	positions = c(positions, map[s,4])
}
toplot = data.frame(pos=depthschr[,1],depth = as.numeric(depthschr$V8), geno= as.numeric(raw[1,]))
toplot1 = toplot[which(toplot$geno == 1),]
toplot2 = toplot[union(which(toplot$geno == 0),which(toplot$geno==2)),]


p1 = ggplot(toplot1,aes(x=pos, y=log10(depth))) + geom_point(size=0.5,alpha=1/5) + theme(legend.position="none",panel.background=element_blank()) + xlim(1.46e+08,1.51e+08)
p2 = ggplot(toplot2,aes(x=pos, y=log10(depth))) + geom_point(size=0.5,alpha=1/5) + theme(legend.position="none",panel.background=element_blank()) + xlim(1.46e+08,1.51e+08)
grid.arrange(p1,p2, ncol=1)
dev.off()
