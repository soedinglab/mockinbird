import argparse
import json
import os
import pickle
import logging

import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mockinbird.utils.fit_betabinom import fit_betabinom_ab


logger = logging.getLogger()


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('transition_table')
    parser.add_argument('max_mixture_components', type=int, default=3)
    parser.add_argument('bam_stat_json')
    parser.add_argument('out_dir')
    parser.add_argument('--no_global_fit', action='store_true')
    parser.add_argument('--n_iterations', default=100, type=int)
    return parser


def plot_fit(k_vals, w, a, title, file_name, out_dir):

    def geom_distr(k, w, a):
        return np.sum(w * (a / (1 + a) ** (k + 1)))

    fit = [geom_distr(k, w, a) for k in range(np.max(k_vals) + 1)]
    fit = np.array(fit)
    fig, ax = plt.subplots()
    ax.hist(k_vals, log=True, bins=50, normed=True, alpha=0.5)
    ax.plot(fit, color='red')
    fig.suptitle(title, fontsize=18)
    fname = os.path.join(out_dir, file_name)
    plt.savefig(fname)


def geom_mm_fit_fast(x, n_components, n_iter=100):
    # unfortunately this time I parameterized with p = (1-q). Shame on me.
    weights = [0.01, 0.99]
    x_mean = np.mean(x)
    q_mean = x_mean / (1 + x_mean)
    x_max = np.max(x)
    q_max = x_max / (1 + x_max)

    p_mean = 1 - q_mean
    p_max = 1 - q_max

    p_initial = [p_mean, p_mean + 0.05 * p_mean]

    for m in range(3, n_components):
        p_initial.append(p_max)
        weights = np.hstack((weights, [0.01]))
        weights /= np.sum(weights)
        ps, weights = fast_geom_mm_fit(x, p_initial, weights, n_iter=5)
    p_initial.append(p_max)
    weights = np.hstack((weights, [0.01]))
    weights /= np.sum(weights)

    ps, weights = fast_geom_mm_fit(x, p_initial, weights, n_iter=n_iter)
    w_m = np.array(weights)
    g_m = np.array([1 - p for p in ps])
    g_m = (1 - g_m) / g_m
    return w_m, g_m


def main():
    parser = create_parser()
    args = parser.parse_args()
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    with open(args.bam_stat_json) as bam_json:
        bam_stat = json.load(bam_json)

    mock_table = pd.read_table(args.transition_table,
                               usecols=['n_mock', 'k_mock'])
    k_vals = mock_table.k_mock
    n_vals = mock_table.n_mock
    coverage = bam_stat['total_coverage']

    parameters = {}
    parameters['coverage'] = coverage
    parameters['n_mixture_components'] = args.max_mixture_components

    # fit p(k|z=0)
    logger.info('Fitting p(k)')
    pw_m, g_m = geom_mm_fit_fast(k_vals, args.max_mixture_components, args.n_iterations)
    plot_fit(k_vals, pw_m, g_m, 'p(k)', 'pk_fit', args.out_dir)
    pg_m = coverage * g_m
    parameters['pk_params'] = pw_m, pg_m

    # fit p(n)
    w, g = geom_mm_fit_fast(n_vals, args.max_mixture_components, args.n_iterations)
    plot_fit(n_vals, w, g, 'p(n)', 'pn_fit', args.out_dir)
    pg = coverage * g
    parameters['pn_params'] = w, pg

    if not args.no_global_fit:
        tab = mock_table[['k_mock', 'n_mock']]
        agg_counts = tab.groupby(['k_mock', 'n_mock']).agg(len)
        weights = np.array(agg_counts)
        weights_nk = agg_counts.reset_index()
        k = np.array(weights_nk.k_mock)
        n = np.array(weights_nk.n_mock)
        assert len(k) == len(n)
        assert len(weights) == len(k)
        alpha, beta = fit_betabinom_ab(n, k, weights=weights)
    else:
        alpha_estim = []
        beta_estim = []
        for n in 3, 4, 5:
            k = mock_table.k_mock[mock_table.n_mock == n]
            a, b = fit_betabinom_ab(n * np.ones(len(k)), k)
            alpha_estim.append(a)
            beta_estim.append(b)
        alpha = np.mean(alpha_estim)
        beta = np.mean(beta_estim)

    parameters['pkn_params'] = alpha, beta

    model_file = os.path.join(args.out_dir, 'model.pkl')
    with open(model_file, 'wb') as pkl:
        pickle.dump(parameters, pkl)


def fast_geom_mm_fit(x, p_init, pi_init, n_iter=250):
    N = len(x)
    p = np.array(p_init)
    pi = np.array(pi_init)
    assert len(p) == len(pi)

    # calculate the responsibility for each distinct x value just once
    # and weight by the number of occurences
    x_agg = np.bincount(x + 1)[1:]
    x_new = np.arange(len(x_agg)) + 1

    n, d = len(x_new), len(p)
    r_matrix = np.zeros((d, n))
    for i in range(n_iter):
        # calculate the responsibilities
        for k in range(d):
            r_matrix[k, :] = pi[k] * (1 - p[k]) ** (x_new - 1) * p[k]
        r_matrix /= r_matrix.sum(axis=0)

        # optimize the parameters
        for k in range(d):
            Nk = np.sum(r_matrix[k, :] * x_agg)
            p[k] = Nk / np.dot(x_agg * x_new, r_matrix[k, :])
            pi[k] = Nk / N

    assert np.isclose(np.sum(pi), 1)
    return p, pi


if __name__ == '__main__':
    main()
