import matplotlib.pyplot as plt
import plot_preamb as pp
import numpy as np
import hashlib
import os
import h5py
from config import cache_dir as CACHE_DIR

pp.run()


def _hash_args(args, kwargs):
    hasher = hashlib.sha256()
    hasher.update((repr(args) + repr(kwargs)).encode("utf-8"))
    return hasher.hexdigest()


def get_maglims(frac, feh, nsamp):
    args = (frac, )
    kw = dict(feh=feh, nsamp=nsamp)
    _hash = _hash_args(args, kw)
    fname = CACHE_DIR + '/' + _hash + '.hdf5'
    if not os.path.exists(fname):
        import maglim_calculator as maglim
        A, B = maglim.doit(*args, **kw)
        with h5py.File(fname, 'w') as R:
            R.create_dataset('A', data=A)
            R.create_dataset('B', data=B)
    with h5py.File(fname, 'r') as R:
        A, B = np.array(R['A']), np.array(R['B'])
    return A, B


frac = 0.2  # completeness in desi DR1
frac2 = 0.38  # in DESI DR2
# frac = 0.2  # in dr1
nsamp = 1e7
A1, B1 = get_maglims(frac, -2.5, nsamp)
A2, B2 = get_maglims(frac, -2., nsamp)
A3, B3 = get_maglims(frac, -1.5, nsamp)
plt.figure(figsize=(3.37, 2.5))
shift = 2.5 * np.log(frac2 / frac)  # dr1 vs dr2
plt.plot(A1, B1, label='[Fe/H]=$-2.5$', color='blue')
plt.gca().set_xscale('log')

plt.xlabel('Distance [kpc]', )
plt.ylabel('M$_V$ [mag]', )
plt.plot(A2, B2, label='[Fe/H]=$-2.0$', color='green')
plt.plot(A3, B3, label='[Fe/H]=$-1.5$', color='red')
plt.plot(A2, B2 + shift, label='DESI DR2', color='grey')
plt.legend()
plt.tight_layout()
plt.fill_between(A1,
                 np.minimum(np.minimum(B1, B2), B3),
                 B3 * 0 - 55,
                 zorder=2,
                 fc='grey',
                 alpha=0.2)
plt.ylim(0, -13)
plt.xlim(10, 250)
plt.text(35, -8, '>10 stars in DESI')
plt.savefig('plots//maglim10.pdf')
