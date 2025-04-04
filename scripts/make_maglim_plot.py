import maglim_calculator as cc
from idlplotInd import plot, oplot
import matplotlib.pyplot as plt
import plot_preamb as pp
import numpy as np

pp.run()

frac = 0.17  # completeness in desi DR1
frac2 = 0.38  # in DESI DR2
# frac = 0.2  # in dr1
nsamp = 1e7
A1, B1 = cc.doit(frac, feh=-2.5, nsamp=nsamp)
A2, B2 = cc.doit(frac, feh=-2, nsamp=nsamp)
A3, B3 = cc.doit(frac, feh=-1.5, nsamp=nsamp)
plt.figure(figsize=(3.37, 2.5))
shift = 2.5 * np.log(frac2 / frac)  # dr1 vs dr2
plot(A1,
     B1,
     xlog=True,
     xtitle='Distance [kpc]',
     ytitle='M$_V$ [mag]',
     label='[Fe/H]=$-2.5$',
     color='blue')
oplot(A2, B2, label='[Fe/H]=$-2.0$', color='green')
oplot(A3, B3, label='[Fe/H]=$-1.5$', color='red')
oplot(A2, B2 + shift, label='DESI DR2', color='grey')
plt.legend()
plt.tight_layout()
plt.fill_between(A1,
                 np.minimum(np.minimum(B1, B2), B3),
                 B3 * 0 - 55,
                 zorder=2,
                 fc='grey',
                 alpha=0.2)
plt.ylim(0, -13)

plt.text(35, -8, '>10 stars in DESI')
plt.savefig('plots//maglim10.pdf')
