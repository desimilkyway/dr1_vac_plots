import minimint
import numpy as np
import imf


def fancy_sampler(iso, imf, Ntot, nbatch=100):
    pass


def doit(hh_frac, thresh=10, feh=-2, nsamp=1e7):
    iso_interp = minimint.Interpolator(['DECam_g', 'DECam_r'])
    imf_distr = imf.chabrier2005.distr
    mass = imf_distr.rvs(int(nsamp))
    ii = iso_interp(mass, np.log10(10e9), feh)
    xind = np.isfinite(ii['DECam_r'])
    print(xind.sum())
    g, r = ii['DECam_g'][xind], ii['DECam_r'][xind]
    totlum = -2.5 * np.log10((10**(-g / 2.5)).sum())
    dists = np.exp(np.linspace(np.log(10), np.log(250), 250))
    res = []
    for curd in dists:
        curr = r + 5 * np.log10(curd * 1e3) - 5
        curhh = np.histogram2d(g - r,
                               curr,
                               range=[[-.3, 1.8], [16, 19]],
                               bins=[30, 30])[0]
        curn = (curhh * hh_frac).sum()
        xrat = curn / thresh
        res.append(totlum + 2.5 * np.log10(xrat))
    return np.asarray(dists), np.asarray(res)
