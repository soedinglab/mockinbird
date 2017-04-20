#!/usr/bin/env Rscript
library(LSD)

s = as.matrix(read.table(commandArgs(TRUE)[1],header=F))
a = as.matrix(read.table(commandArgs(TRUE)[2],header=F))
outputfile  = commandArgs(TRUE)[3]
fname       = commandArgs(TRUE)[4]
rStart      = as.numeric(commandArgs(TRUE)[5])
rStop       = as.numeric(commandArgs(TRUE)[6])
body        = as.numeric(commandArgs(TRUE)[7])
hw          = as.numeric(commandArgs(TRUE)[8])
centerName1 = commandArgs(TRUE)[9]
midName     = commandArgs(TRUE)[10]
centerName2 = commandArgs(TRUE)[11]


rmean2 = function(x, hw=10){
    smx = numeric(length(x))
    smx[1:(hw-1)]                   = rep(mean(x[1:hw]), length(1:(hw-1)))
    smx[(length(x)-hw+1):length(x)] = rep(mean(x[(length(x)-hw+1):length(x)]), length((length(x)-hw+1):length(x)))
    for(i in hw:(length(x)-hw)){
        smx[i] = mean(x[(i-hw):(i+hw)])
    }
    return(smx)
}

plotBothEnds = function(sense, asense, relativeSTART, relativeSTOP, bodysize, centerNameA, midName, centerNameB, main=''){
    lwd  = 3
    gradient_sense   = rev(colorpalette('blues', 100))
    gradient_asense  = rev(colorpalette('greens', 100))
    maxy = max(c(sense,asense))
    shift = 50
    plot(x=10000, y=10000, xlim=c(1,length(sense)+shift), ylim=c(0,maxy), ylab='Averaged Occupancy', xlab='Genomic position', xaxt='n', main=main)
        sel = 1:(relativeSTART+bodysize)
        lines(x=sel, sense[sel],  lwd=lwd, col=gradient_sense[95])
        lines(x=sel, asense[sel], lwd=lwd, col=gradient_asense[95])
        #abline(v=c(relativeSTART, (relativeSTART+bodysize)), lwd=2, lty = 2 , col='black')
        abline(v=relativeSTART, lwd=1, lty = 1 , col='black')
        abline(v=relativeSTART+bodysize, lwd=1, lty = 2 , col='darkgrey')
        abline(v=relativeSTART+bodysize+shift, lwd=1, lty = 2 , col='darkgrey')
        sel = (relativeSTART+bodysize):length(sense)
        lines(x=sel+shift, sense[sel],  lwd=lwd, col=gradient_sense[95])
        lines(x=sel+shift, asense[sel], lwd=lwd, col=gradient_asense[95])
        abline(v=length(sense)+shift-relativeSTOP, lwd=1, lty = 1 , col='black')
        axis(side=1, at=c(1,round(relativeSTART/2), relativeSTART, relativeSTART+round(bodysize/2), relativeSTART+bodysize+round(shift/2), 
        relativeSTART+bodysize+shift+round(bodysize/2), length(sense)+shift-relativeSTOP,length(sense)+shift-round(relativeSTOP/2),length(sense)+shift),
        labels=c(paste('-',relativeSTART), paste('-',round(relativeSTART/2)), centerNameA, paste('+',round(bodysize/2)), 
        midName, paste('-',round(bodysize/2)), centerNameB, paste('+',round(relativeSTOP/2)), relativeSTOP))
}

plotAsMat = function(sense, asense, relativeSTART, relativeSTOP, bodysize, centerNameA, midName, centerNameB, main=''){
    mat = rbind(sense, asense)
    image(x=1:ncol(mat), y=1:nrow(mat), z=t(mat)[,nrow(mat):1], col=rev(colorpalette('rdbu', 100)), xaxt='n', yaxt='n', ylab='', xlab='Genomic position', main=main)
    abline(v=relativeSTART+bodysize, lwd=3)
    axis(side=1, at=c(1,round(relativeSTART/2), relativeSTART, relativeSTART+round(bodysize/2), relativeSTART+bodysize, 
    relativeSTART+bodysize+round(bodysize/2), length(sense)-relativeSTOP,length(sense)-round(relativeSTOP/2),length(sense)),
    labels=c(paste('-',relativeSTART), paste('-',round(relativeSTART/2)), centerNameA, paste('+',round(bodysize/2)), 
    midName, paste('-',round(bodysize/2)), centerNameB, paste('+',round(relativeSTOP/2)), relativeSTOP))
    axis(side=2, at=c(1,2), labels=c('asense', 'sense'))
}

css = colSums(s)
csa = colSums(a)

css = css/nrow(s)
csa = csa/nrow(a)


gradient_sense   = rev(colorpalette('blues', 100))
gradient_asense  = rev(colorpalette('greens', 100))
lwd  = 3
maxy = max(c(css,csa))

#raw
pdf(file=outputfile, width=11, height=12.5)
par(mfcol=c(4,1))
plotBothEnds(css, csa, rStart, rStop, body, centerName1, midName, centerName2, paste('Raw signal', fname))

#smooth
css_smooth = rmean2(css, hw=hw)
csa_smooth = rmean2(csa, hw=hw)

maxy = max(c(css_smooth,csa_smooth))

plotBothEnds(css_smooth, csa_smooth, rStart, rStop, body, centerName1, midName, centerName2, paste('Smoothed Signal [window size=',(2*hw),']', fname))

#rescale to unit interval
miny = min(css_smooth, csa_smooth)
css_smooth_r = css_smooth - miny
csa_smooth_r = csa_smooth - miny

maxy = max(c(css_smooth_r, csa_smooth_r))

if (maxy > 0) {
	css_smooth_r = css_smooth_r / maxy
	csa_smooth_r = csa_smooth_r / maxy
}

plotBothEnds(css_smooth_r, csa_smooth_r, rStart, rStop, body, centerName1, midName, centerName2, paste('Smoothed signal is rescaled to unit interval', fname))

plotAsMat(css_smooth_r, csa_smooth_r, rStart, rStop, body, centerName1, midName, centerName2, paste('Smoothed and rescaled signal as heatplot', fname))
dev.off()













