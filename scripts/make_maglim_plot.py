import maglim_calculator as cc
import matplotlib.pyplot as plt
import plot_preamb as pp
import numpy as np

pp.run()

frac = 0.2  # completeness in desi DR1
frac2 = 0.38  # in DESI DR2
# frac = 0.2  # in dr1
nsamp = 1e7
A1, B1 = cc.doit(frac, feh=-2.5, nsamp=nsamp)
A2, B2 = cc.doit(frac, feh=-2, nsamp=nsamp)
A3, B3 = cc.doit(frac, feh=-1.5, nsamp=nsamp)
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
