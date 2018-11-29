library(data.table)
library(ggplot2)
library(reshape2)

CHR = 1
print(CHR)
F = fread(paste("chr",CHR,"raw.raw",sep=""),skip=1,stringsAsFactors=FALSE)
F = as.data.frame(F[,-c(1:6)])
N1s = matrix(0, nrow=8, ncol=round(ncol(F)/1000))

for(i in 0:(ncol(N1s)-2)){
		s = (i*1000)+1
		e = (i*1000)+1000
		ns = apply(F[,s:e], 1, function(x){length(which(x==1))})
		N1s[,i+1] = ns
}

names = c("ACI","BN","BUF","F344","M520","MR","WKY","WN")
toplot = data.frame(names, N1s)
df = melt(toplot, id=c("names"))
df[,2] = as.numeric(gsub("X","",df[,2]))
colnames(df)[colnames(df) == "value"] = 'Val'
colnames(df)[colnames(df) == 'variable'] = paste("Chr",CHR,sep="")

ggplot(df, aes_string(x=colnames(df)[2], y=colnames(df)[3],group=colnames(df)[1])) + geom_line() + facet_grid(names ~ .) + theme_bw() + labs(x="1000-marker windows along Chromosome 1",y="Fraction of the heterozygosity in window")
ggsave(paste("/scratch/junzli_flux/sramdas/haplotyping/founders/rn6/unifiedgenotyper/chr",CHR,"_1000markers.pdf",sep=""),width=7,height=8,device="pdf")
