library(LSD)
#rows = reference
#cols = mutations
#        m u t 
#       A C G T
# r  A         
# e  C
# f  G
#    T
makeJointProb = function(mat){
    mat = rbind(mat, colSums(mat))
    mat = cbind(mat, rowSums(mat))
    mat[1:4,1:4] = mat[1:4,1:4]/mat[5,5]
    mat[5,1:4] = mat[5,1:4]/mat[5,5]
    mat[1:4,5] = mat[1:4,5]/mat[5,5]
    mat[5,5] = 1
    return(mat)
}

plotConditionals = function(jointmat, maxy=0, maxToZero=TRUE, main=''){
    col_A = '#008000'
    col_C = '#0000FF'
    col_G = '#FFA500'
    col_T = '#FF0000'
    mut_given_A = jointmat[1,1:4]/jointmat[1,5]
    mut_given_C = jointmat[2,1:4]/jointmat[2,5]
    mut_given_G = jointmat[3,1:4]/jointmat[3,5]
    mut_given_T = jointmat[4,1:4]/jointmat[4,5]
    if(maxToZero){
        mut_given_A[1] = 0
        mut_given_C[2] = 0
        mut_given_G[3] = 0
        mut_given_T[4] = 0
    }
    if(maxy == 0){
        maxy = max(c(mut_given_A, mut_given_C, mut_given_G, mut_given_T))
    }
    plot(x=-100, y=-100, xlim=c(1,19), ylim=c(0, maxy),
         ylab='P(mutation|reference)', xlab='', xaxt='n', main=main)
     rect(1,0,1.75,mut_given_A[1], col=col_A)
     rect(2,0,2.75,mut_given_A[2], col=col_A)
     rect(3,0,3.75,mut_given_A[3], col=col_A)
     rect(4,0,4.75,mut_given_A[4], col=col_A)
    
     rect(6,0,6.75,mut_given_C[1], col=col_C)
     rect(7,0,7.75,mut_given_C[2], col=col_C)
     rect(8,0,8.75,mut_given_C[3], col=col_C)
     rect(9,0,9.75,mut_given_C[4], col=col_C)
    
     rect(11,0,11.75,mut_given_G[1], col=col_G)
     rect(12,0,12.75,mut_given_G[2], col=col_G)
     rect(13,0,13.75,mut_given_G[3], col=col_G)
     rect(14,0,14.75,mut_given_G[4], col=col_G)
    
     rect(16,0,16.75,mut_given_T[1], col=col_T)
     rect(17,0,17.75,mut_given_T[2], col=col_T)
     rect(18,0,18.75,mut_given_T[3], col=col_T)
     rect(19,0,19.75,mut_given_T[4], col=col_T)
     axis(s=1, at=c(1,2,3,4,6,7,8,9,11,12,13,14,16,17,18,19)+0.375,
          labels=c('P(A|A)','P(C|A)','P(G|A)','P(T|A)','P(A|C)','P(C|C)',
                   'P(G|C)','P(T|C)','P(A|G)','P(C|G)','P(G|G)','P(T|G)',
                   'P(A|T)','P(C|T)','P(G|T)','P(T|T)'), las=2)
    abline(v=c(5.375,10.375,15.375))
}

infile    = commandArgs(TRUE)[1]
outfile   = commandArgs(TRUE)[2]
limit     = as.numeric(commandArgs(TRUE)[3])

mat = as.matrix(read.table(infile, header=F))

pdf(file=outfile, width=7, height=7)
    plotConditionals(makeJointProb(mat), limit)
dev.off()










