#!/usr/bin/env Rscript
library(LSD)

rmean2 = function(x, hw=10){
    smx = numeric(length(x))
    smx[1:(hw-1)]                   = rep(mean(x[1:hw]), length(1:(hw-1)))
    smx[(length(x)-hw+1):length(x)] = rep(mean(x[(length(x)-hw+1):length(x)]), length((length(x)-hw+1):length(x)))
    for(i in hw:(length(x)-hw)){
        smx[i] = mean(x[(i-hw):(i+hw)])
    }
    return(smx)
}

plotKmers = function(inputfile, main, motifs=character(), smooth=10, lwd=3, cex.axis=1.2, r=50, mid=101, colpal='spectral', plotQuantiles=0){
    kmer_pos       = read.table(inputfile, header=F)
    kmer_mat       = as.matrix(kmer_pos[,2:ncol(kmer_pos)])
    
    kmers          = as.character(kmer_pos[,1])
    cols           = colorpalette(colpal, length(kmers))
    if(length(motifs)>1){
        cols           = colorpalette(colpal, length(motifs))
    }
    
    savemax        = numeric(length(kmers))
    names(savemax) = kmers
    
    supersel = numeric()
    ylim     = c(0,0)
    plotmat = matrix(0, nr=nrow(kmer_mat), nc=ncol(kmer_mat))
    for(i in 1:nrow(kmer_mat)){
        x = kmer_mat[i,]
        if(smooth > 0){ x = rmean2(x, hw=smooth) }
        x[x==-Inf] = 0
        plotmat[i, ] = x
        savemax[i] = max(x)
        if(max(x) > ylim[2]){ ylim[2] = max(x)}
    }
    if(length(motifs)!=0){
        plotmat = matrix(0, nr=length(motifs), nc=ncol(kmer_mat))
        for(i in 1:length(motifs)){
            sel = which(kmers==motifs[i])
            supersel = c(supersel, sel)
            x = kmer_mat[sel,]
            if(smooth > 0){ x = rmean2(x, hw=smooth) }
            x[x==-Inf] = 0
            plotmat[i, ] = x
            savemax[i] = max(x)
        }
    }
    
    sel = (mid-r):(mid+r)
    plotmat=plotmat[, sel]
    plot(x=-100, y=-100, xlim=c(1,length(sel)), ylim=ylim, ylab="absolute counts", xlab="genomic position", xaxt="n", main=main, cex.axis=cex.axis, cex.lab=cex.axis)
    
    for(i in 1:nrow(plotmat)){
        points(plotmat[i,], type="l", col=cols[i], lwd=lwd)
    }
    if(plotQuantiles != 0){
        baseline_75 = sapply(1:ncol(plotmat), function(z){ return(quantile(plotmat[,z], 0.75)) })
        baseline_50 = sapply(1:ncol(plotmat), function(z){ return(quantile(plotmat[,z], 0.50)) })
        baseline_25 = sapply(1:ncol(plotmat), function(z){ return(quantile(plotmat[,z], 0.25)) })
        points(baseline_75, type="l", col="lightgrey", lwd=lwd-1)
        points(baseline_50, type="l", col="darkgrey", lwd=lwd-1)
        points(baseline_25, type="l", col="lightgrey", lwd=lwd-1)
    }
    size = ncol(plotmat)
    
    #abline(v=(size/2), lty=2, col="grey")
    abline(v=r+1, lty=2, col="grey")
    axis(side=1, at=c(1,r+1,ncol(plotmat)), labels=c(paste("-",r, sep=""), 0, paste("+",r, sep="")), cex.axis=cex.axis)
    #axis(side=1, at=c(1,(size/2),size), labels=c(paste("-",round(size/2), sep=""), 0, paste("+",round(size/2), sep="")), cex.axis=cex.axis)
    box()
    bla = rev(sort(savemax))
    cat(main,'\n')
    legendnames= character()
    legendcols = character()
    if(length(motifs)==0){
        for(i in 1:10){
            legendnames = c(legendnames, names(bla)[i])
            sel = which(kmers == names(bla)[i])
            legendcols = c(legendcols, cols[sel])
            cat('\t',names(bla)[i],'\t',bla[i],'\n')
        }
    }else{
        for(i in 1:length(motifs)){
            legendnames = c(legendnames, motifs[i])
            legendcols  = c(legendcols, cols[i])
        }
    }
    legend("topright", lwd=rep(lwd,10), lty=rep(1,10), legend=legendnames, col=legendcols)
}

inputfile    = commandArgs(TRUE)[1]
outputfile   = commandArgs(TRUE)[2]
width        = as.numeric(commandArgs(TRUE)[3])
smooth       = as.numeric(commandArgs(TRUE)[4])
mid          = as.numeric(commandArgs(TRUE)[5])
#pQuantiles   = as.numeric(commandArgs(TRUE)[6])
#motifs       = commandArgs(TRUE)[7]


pdf(file=outputfile, width=13,height=5.5)
plotKmers(inputfile, main="", motifs=character(), smooth=smooth, r=width, mid=mid, plotQuantiles=1)
dev.off()






