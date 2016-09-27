#!/usr/bin/env Rscript
library(LSD)

sense = as.matrix(read.table(commandArgs(TRUE)[1],header=F))
sense_bs = as.matrix(read.table(commandArgs(TRUE)[2],header=F))
asense = as.matrix(read.table(commandArgs(TRUE)[3],header=F))
asense_bs = as.matrix(read.table(commandArgs(TRUE)[4],header=F))

outputfile  = commandArgs(TRUE)[5]
fname       = commandArgs(TRUE)[6]
rStart      = as.numeric(commandArgs(TRUE)[7])
rStop       = as.numeric(commandArgs(TRUE)[8])
body        = as.numeric(commandArgs(TRUE)[9])
hw          = as.numeric(commandArgs(TRUE)[10])
centerName1 = commandArgs(TRUE)[11]
midName     = commandArgs(TRUE)[12]
centerName2 = commandArgs(TRUE)[13]


plotBothEnds = function(sense, sense_bs, asense, asense_bs, relativeSTART, relativeSTOP, bodysize, centerNameA, midName, centerNameB, main=''){
    lwd  = 3

    bss_quantile <- apply(sense_bs, 2, quantile, probs = c(0.025, 0.975),  na.rm = TRUE)
    bsa_quantile <- apply(asense_bs, 2, quantile, probs = c(0.025, 0.975),  na.rm = TRUE)
    gradient_sense   = rev(colorpalette('blues', 100))
    gradient_asense  = rev(colorpalette('greens', 100))
    s_col = gradient_sense[95]
    s_ci_col = paste(s_col, '40', sep='')
    a_col = gradient_asense[95]
    a_ci_col = paste(a_col, '40', sep='')

    max_sbs = max(bss_quantile)
    max_abs = max(bsa_quantile)
    maxy = max(c(sense, asense, max_sbs, max_abs))
    shift = 50
    plot(x=10000, y=10000, xlim=c(1,length(sense)+shift), ylim=c(0,maxy), ylab='Averaged Occupancy', xlab='Genomic position', xaxt='n', main=main)
        sel = 1:(relativeSTART+bodysize+1)
        lines(x=sel, sense[sel],  lwd=lwd, col=s_col)
        lines(x=sel, asense[sel], lwd=lwd, col=a_col)
        polygon(c(sel, rev(sel)), c(bss_quantile[1,][sel], rev(bss_quantile[2,][sel])), col=s_ci_col, border=NA)
        polygon(c(sel, rev(sel)), c(bsa_quantile[1,][sel], rev(bsa_quantile[2,][sel])), col=a_ci_col, border=NA)

        #abline(v=c(relativeSTART, (relativeSTART+bodysize)), lwd=2, lty = 2 , col='black')
        abline(v=relativeSTART, lwd=1, lty = 1 , col='black')
        abline(v=relativeSTART+bodysize, lwd=1, lty = 2 , col='darkgrey')
        abline(v=relativeSTART+bodysize+shift, lwd=1, lty = 2 , col='darkgrey')
        sel = (relativeSTART+bodysize+2):length(sense)
        lines(x=sel+shift, sense[sel],  lwd=lwd, col=s_col)
        lines(x=sel+shift, asense[sel], lwd=lwd, col=a_col)
        polygon(c(sel + shift, rev(sel + shift)), c(bss_quantile[1,][sel], rev(bss_quantile[2,][sel])), col=s_ci_col, border=NA)
        polygon(c(sel + shift, rev(sel + shift)), c(bsa_quantile[1,][sel], rev(bsa_quantile[2,][sel])), col=a_ci_col, border=NA)
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

pdf(file=outputfile, width=11, height=12.5)
par(mfcol=c(2,1))

plotBothEnds(sense, sense_bs, asense, asense_bs, rStart, rStop, body, centerName1, midName, centerName2, paste('Signal', fname))

plotAsMat(sense, asense, rStart, rStop, body, centerName1, midName, centerName2, paste('Signal as heatplot', fname))
dev.off()













