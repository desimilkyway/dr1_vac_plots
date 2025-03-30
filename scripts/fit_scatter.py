import numpy as np
import dynesty
from dataclasses import dataclass


def _trans_uni(x, a, b):
    return a + (b - a) * x


def _trans_loguni(x, a, b):
    return a * np.exp(x * np.log(b / a))


@dataclass
class Prior:
    max_bg_frac: float = 0.3
    min_obj_sig: float = 0.02
    max_obj_sig: float = 0.7
    min_bg_sig: float = 0.7
    max_bg_sig: float = 1.
    min_bg_cen: float = -2
    max_bg_cen: float = 0
    min_obj_cen: float = -4
    max_obj_cen: float = 0.5

    def __call__(self, p):
        """
        Prior transformation function.
        ----------
        p : array-like of shape (5,)
            p[0] in [0,1]: fraction -> scale to [0, max_bg_frac]
            p[1] in [0,1]: object center -> scale to [min_obj_cen, max_obj_cen]
            p[2] in [0,1]: object sigma -> log-scale to [min_obj_sig, max_obj_sig]
            p[3] in [0,1]: background center -> scale to [min_bg_cen, max_bg_cen]
            p[4] in [0,1]: background sigma -> log-scale to [min_bg_sig, max_bg_sig]
        """
        # Create an output array filled with NaN (same shape as p).
        p1 = np.full_like(p, np.nan)

        # 1) Fraction in [0, max_bg_frac]
        p1[0] = _trans_uni(p[0], 0, self.max_bg_frac)

        # 2) Object center in [min_obj_cen, max_obj_cen]
        p1[1] = _trans_uni(p[1], self.min_obj_cen, self.max_obj_cen)
        # 3) Object sigma in [min_obj_sig, max_obj_sig]
        p1[2] = _trans_uni(p[2], self.min_obj_sig, self.max_obj_sig)

        # 4) Background center in [min_bg_cen, max_bg_cen]
        p1[3] = _trans_uni(p[3], self.min_bg_cen, self.max_bg_cen)
        # 5) Background sigma in [min_bg_sig, max_bg_sig]
        p1[4] = _trans_uni(p[4], self.min_bg_sig, self.max_bg_sig)
        return p1


def like(p, args):
    v, ev = args
    bg_frac, obj_cen, obj_sig, bg_cen, bg_sig = p

    # Calculate the log-likelihood for the given parameters
    # This is a placeholder implementation; replace with actual calculation
    obj_sig1 = np.sqrt(obj_sig**2 + ev**2)
    bg_sig1 = np.sqrt(bg_sig**2 + ev**2)
    logp_obj = -0.5 * (((v - obj_cen) / obj_sig1)**2) - np.log(obj_sig1)
    logp_bg = -0.5 * (((v - bg_cen) / bg_sig1)**2) - np.log(bg_sig1)
    logp = np.logaddexp(
        np.log(bg_frac) + logp_bg,
        np.log(1 - bg_frac) + logp_obj)
    return logp.sum()


def fitter(v, ev, seed=43632, max_bg_frac=0.3):
    ndim = 5
    args = (v, ev)
    rstate = np.random.default_rng(seed=seed)
    transf = Prior(max_bg_frac=max_bg_frac)
    dns = dynesty.DynamicNestedSampler(
        like,
        transf,
        ndim,
        sample='rslice',
        logl_args=(args, ),
        rstate=rstate,
    )
    dns.run_nested(n_effective=1000, print_progress=False)
    samps = dns.results.samples_equal(rstate=rstate)
    samps_d = {}
    # Convert the samples to a dictionary with meaningful keys
    key_names = ['bg_frac', 'obj_cen', 'obj_sig', 'bg_cen', 'bg_sig']
    for i, key in enumerate(key_names):
        samps_d[key] = samps[:, i]
    return samps_d


def get_hist(x, min, max, nbins=10):
    """
    Create a histogram of the data.
    ----------
    x : array-like
        Data to be histogrammed.
    min : float
        Minimum value for the histogram.
    max : float
        Maximum value for the histogram.
    nbins : int, optional
        Number of bins (default is 10).
    ----------
    Returns:
        hist : array-like
            Histogram counts.
    """
    hist = np.histogram(x, bins=nbins, range=(min, max))
    return hist[0]


def get_scatter(v, ev, max_bg_frac=0.3, seed=43632):
    """
    Fit the scatter of the data using a nested sampling approach.
    ----------
    v : array-like
        Velocities of the data points.
    ev : array-like
        Errors of the data points.
    seed : int, optional
        Random seed for reproducibility (default is 43632).
    ----------
    Returns:
        
    """
    samps_d = fitter(v, ev, seed=seed, max_bg_frac=max_bg_frac)
    # here we first want to flag the case where the background fraction is close to the
    # maximum value, which indicates a bad fit
    bg_frac_hist = get_hist(samps_d['bg_frac'], 0, max_bg_frac, nbins=10)
    warn_list = []
    if bg_frac_hist[:-1].max() < bg_frac_hist.max():
        warn_list.append('bg_frac near the right edge')
    if len(warn_list) > 0:
        warn_flag = True
    else:
        warn_flag = False
    obj_sig_hist = get_hist(samps_d['obj_sig'], 0, max_bg_frac, nbins=100)
    if obj_sig_hist[1:].max() < obj_sig_hist.max():
        # limit
        ret = [
            -np.inf,
            np.percentile(samps_d['obj_sig'], 50),
            np.percentile(samps_d['obj_sig'], 84)
        ]
    else:
        ret = [
            np.percentile(samps_d['obj_sig'], 16),
            np.percentile(samps_d['obj_sig'], 50),
            np.percentile(samps_d['obj_sig'], 84),
        ]
    ret_mean = [
        np.percentile(samps_d['obj_cen'], 16),
        np.percentile(samps_d['obj_cen'], 50),
        np.percentile(samps_d['obj_cen'], 84),
    ]

    return ret_mean, ret, warn_flag, warn_list
