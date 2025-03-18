import numpy as np
import scipy.optimize

arr = np.loadtxt('compar.dat')
rv_logg, rv_teff, rv_feh, sp_logg, sp_teff, sp_feh, feh_ref = arr.T


def func(p, data):
    logg, teff, feh, feh_ref, cut, get_corr = data
    X = (teff - 5500) / 2000
    xind = logg < cut
    p1 = p[:3]
    p2 = p[3:]
    V1 = np.poly1d(p1)(X)
    V2 = np.poly1d(p2)(X)
    V1[xind] = V2[xind]
    if get_corr:
        return V1
    ret = np.abs(feh - feh_ref - V1).sum()
    print(ret, p)
    return ret


cuts = np.linspace(0, 6, 100)
funcs = []
for cut in cuts:
    RR = scipy.optimize.minimize(
        func, np.zeros(6), ((sp_logg, sp_teff, sp_feh, feh_ref, cut, False), ))
    funcs.append(RR.fun)
