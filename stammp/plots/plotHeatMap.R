library(LSD)

input_sense    = commandArgs(TRUE)[1]
input_asense   = commandArgs(TRUE)[2]
outfile_sense  = commandArgs(TRUE)[3]
outfile_asense = commandArgs(TRUE)[4]
qvalue         = as.numeric(commandArgs(TRUE)[5])
start          = as.numeric(commandArgs(TRUE)[6])
ysize          = as.numeric(commandArgs(TRUE)[7])
xsize          = as.numeric(commandArgs(TRUE)[8])

sense = as.matrix(read.table(input_sense,header=F))
asense  = as.matrix(read.table(input_asense,header=F))
pA = sense[, 1]

sense = sense[,2:ncol(sense)]
sense = sense - min(sense)

asense = asense[,2:ncol(asense)]
asense = asense - min(asense)

q = quantile(cbind(sense,asense), qvalue)

sense[sense>q] = q
sense = sense/q

asense[asense>q] = q
asense = asense/q

png(file=outfile_sense, width=xsize, height=ysize, units="px", res=300)
    par(mar=c(0,0,0,0))
    image(x=1:ncol(sense), y=1:nrow(sense), z=t(sense)[,nrow(sense):1], zlim=c(0,1), col=rev(colorpalette('blues',100)), xlab='', ylab='', xaxt='n', yaxt='n')
    lines(x=pA,y=length(pA):1, lwd=0.5)
    abline(v=start, lwd=0.5)
dev.off()

png(file=outfile_asense, width=xsize, height=ysize, units="px", res=300)
    par(mar=c(0,0,0,0))
    image(x=1:ncol(asense), y=1:nrow(asense), z=t(asense)[,nrow(asense):1], zlim=c(0,1), col=rev(colorpalette('greens',100)), xlab='', ylab='', xaxt='n', yaxt='n')
    lines(x=pA,y=length(pA):1, lwd=0.5)
    abline(v=start, lwd=0.5)
dev.off()
