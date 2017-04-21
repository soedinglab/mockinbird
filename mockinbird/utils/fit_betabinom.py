from functools import lru_cache
import warnings

import numpy as np

from rpy2 import robjects as ro
from rpy2.robjects.packages import STAP
import rpy2.robjects.numpy2ri as rpyn

from scipy.special import betaln


FIT_RFUN_STR = """
options(warn=-1)
library(VGAM)

fit_betabinom <- function(n, k) {
    fit <- vglm(cbind(k, n-k) ~ 1, betabinomialff)
    return(coef(fit, matrix = TRUE))
}

fit_betabinom_w <- function(n, k, w) {
    fit <- vglm(cbind(k, n-k) ~ 1, betabinomialff, weights=w)
    return(coef(fit, matrix = TRUE))
}

eval_betabinom <- function(n, k, a, b) {
    p <- dbetabinom.ab(k, size=n, shape1=a, shape2=b)
    return(p)
}
"""

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    bbfit = STAP(FIT_RFUN_STR, 'bbfit')


def fit_betabinom_ab(n, k, weights=None):
    assert np.all((k >= 0) & (k <= n))
    n_r = ro.FloatVector(n)
    k_r = ro.FloatVector(k)
    if weights is not None:
        assert len(weights) == len(k)
        assert len(weights) == len(n)
        weights_r = ro.FloatVector(weights)
        result_r = bbfit.fit_betabinom_w(n_r, k_r, weights_r)
    else:
        result_r = bbfit.fit_betabinom(n_r, k_r)
    result = rpyn.ri2py(result_r)
    log_a, log_b = result.flatten()
    return np.exp(log_a), np.exp(log_b)


def p_kn_wrapper(a, b):
    @lru_cache(maxsize=2**15)
    def p_kn(k, n):
        #p = binom(n, k) * beta(a+k,b+n-k) / beta(a, b)
        p = np.exp(
            - betaln(1 + n - k, 1 + k)
            - np.log(n + 1)
            + betaln(a + k, b + n - k)
            - betaln(a, b)
        )
        assert 0 <= p <= 1, "(n=%s, k=%s) -> p=%s" % (n, k, p)
        return p
    return p_kn
