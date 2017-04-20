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

sense  = sense[,2:ncol(sense)]
asense = asense[,2:ncol(asense)]

sense[sense==-1]   = NA
asense[asense==-1] = NA

sense[sense<0]=0
asense[asense<0]=0

q = quantile(cbind(sense,asense), qvalue, na.rm=T)
if(q == 0){
    cat('Warning: data set seems very sparse. The 0.98 quantile of the matrices is 0. Quantile is set to the maximum of all values!\n')
    q = max(c(sense,asense), na.rm=T)
}

sense[sense>q] = q
sense = sense/q

asense[asense>q] = q
asense = asense/q

sense[is.na(sense)] = 0
asense[is.na(asense)] = 0

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
