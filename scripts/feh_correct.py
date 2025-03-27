import numpy as np


def betw(x, x1, x2):
    return (x >= x1) & (x < x2)


def calibrate(feh_orig, teff, logg, pipeline='RV', release='DR1'):
    assert release == 'DR1'
    logteff_scale = 0.1
    teff_ref = 5000
    bounds = {
        'RVS': [[-0.5, 5.1], [4200, 6600]],
        'SP': [[0, 5.4], [4000, 6700]]
    }
    cuts = {'RVS': 4.4, 'SP': 3.6}
    coeffs_low = {'RVS': [0.197, -0.785, .417], 'SP': [0.112, 0.195, -.165]}
    coeffs_high = {'RVS': [0.038, -0.281, .084], 'SP': [0.041, -0.095, -.117]}

    coeff_l = coeffs_low[pipeline]
    coeff_h = coeffs_high[pipeline]
    cut = cuts[pipeline]
    bound = bounds[pipeline]
    if pipeline not in ["RVS", "SP"]:
        raise RuntimeError('oops, unknown pipeline')

    subset = betw(logg, bound[0][0], bound[0][1]) & betw(
        teff, bound[1][0], bound[1][1])
    ret = feh_orig * 0 + np.nan
    subset_l = subset & (logg < cut)
    subset_h = subset & (logg >= cut)
    for cur_sub, cur_coeff in zip([subset_l, coeff_l], [subset_h, coeff_h]):
        ret[cur_sub] = (feh_orig[cur_sub] - np.poly1d(cur_coeff)(
            np.log10(teff[cur_sub] / teff_ref) / logteff_scale))
    return ret
