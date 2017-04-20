#!/usr/bin/env Rscript
library(LSD)

infile      = commandArgs(TRUE)[1]
outfile     = commandArgs(TRUE)[2]

d = read.table(infile, header=F)
kmers = as.character(d[,1])
mat   = t(as.matrix(d[,2:ncol(d)]))[,nrow(d):1]

klength = length(strsplit(kmers[1],'')[[1]])

if(klength == 2){
    limit = 2.5
}
if(klength == 3){
    limit = 2.5
}
if(klength == 4){
    limit = 3.0
}
if(klength == 5){
    limit = 4.0
}
if(klength > 5){
    limit = 4.0
}

mat[mat>limit]    = limit
mat[mat<(-limit)] = (-limit)

pdf(file=outfile, width=nrow(mat)*0.15, height=ncol(mat)*0.15)
    image(x=1:nrow(mat), y=1:ncol(mat), z=mat, zlim=c(-limit, limit), xlab='PAR-CLIP bins', sub='1000 sites per bin', ylab='', yaxt='n', xaxt='n', col=colorpalette('piyg',100))
    axis(s=2, at=length(kmers):1, labels=kmers, las=2, cex.axis=0.6)
dev.off()



