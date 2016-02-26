#!/usr/bin/env Rscript
library(LSD)
m = as.matrix(read.table(commandArgs(TRUE)[1],header=F))
zlimit  = as.numeric(commandArgs(TRUE)[2])
outfile = commandArgs(TRUE)[3]
fnameA  = commandArgs(TRUE)[4]
fnameB  = commandArgs(TRUE)[5]

m[m>zlimit] = zlimit
pdf(file=outfile, width=6, height=6)
    image(x=1:ncol(m), y=1:nrow(m), z=t(m)[,nrow(m):1], zlim=c(0,zlimit), col=rev(colorpalette('greens', 100)), xaxt='n', yaxt='n', xlab=fnameB, ylab=fnameA)
    axis(s=1, at=1:ncol(m), labels=(ncol(m):1)*10)
    axis(s=3, at=1:ncol(m), labels=(ncol(m):1)*10)
    axis(s=2, at=1:nrow(m), labels=(1:nrow(m))*10, las=2)
    axis(s=4, at=1:nrow(m), labels=(1:nrow(m))*10, las=2)
dev.off()
