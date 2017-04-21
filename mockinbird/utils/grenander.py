import numpy as np

from rpy2.robjects.packages import importr, STAP
from rpy2 import robjects as ro
import rpy2.robjects.numpy2ri as rpyn


PVAL_ECDF_STR = '''
ecdf_pval <- function (x)
{
    # compute empirical CDF as usual
    x = sort(x)
    n = length(x)
    if (n < 1)
        stop("'x' must have 1 or more non-missing values")
    vals = sort(unique(x))
    F.raw = cumsum(tabulate(match(x, vals)))/n

    # if necessary add an atom at 1 to make it a proper CDF
    if (vals[length(vals)] != 1)
    {
       F.raw = c(F.raw, 1)
       vals = c(vals, 1)
    }

    # if necessary also add an atom at 0 with weight zero to get support [0,1]
    if (vals[1] != 0)
    {
       F.raw = c(0, F.raw)
       vals = c(0, vals)
    }

    rval <- approxfun(vals, F.raw,
        method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
    class(rval) = c("ecdf", "stepfun", class(rval))
    attr(rval, "call") <- sys.call()
    rval
}

require(fdrtool)
plot_ecdf_lcm <- function(pvals, outfile, title) {
    require(fdrtool)
    x_uniq <- sort(unique(pvals))
    pval_ecdf <- ecdf(pvals)
    y <- pval_ecdf(x_uniq)
    res <- fdrtool::gcmlcm(x_uniq, y, type="lcm")
    pdf(outfile)
    plot(pval_ecdf, main=title, xlab='p-value', ylab='cumul. density', lwd=2)
    lines(res$x.knots, res$y.knots, col='red', lwd=2)
    dev.off()
}
'''

ecdf_pkg = STAP(PVAL_ECDF_STR, 'ecdf_pkg')


def pval_grenander_fit(pvals):
    fdrtool = importr('fdrtool')
    pval_vec = ro.FloatVector(pvals)
    pval_ecdf = ecdf_pkg.ecdf_pval(pval_vec)
    gre_fit = fdrtool.grenander(pval_ecdf, type='decreasing')
    x_knots = rpyn.ri2py(gre_fit.rx2('x.knots'))
    f_knots = rpyn.ri2py(gre_fit.rx2('f.knots'))
    if len(f_knots) == 0:
        x_knots = np.array([0, 1])
        f_knots = np.array([1, 1])
    assert len(f_knots) == len(x_knots)
    return x_knots, f_knots


def plot_cumul_density(pvals, plot_basename, title):
    pval_vec = ro.FloatVector(pvals)
    ecdf_pkg.plot_ecdf_lcm(pval_vec, plot_basename + '.pdf', title)
