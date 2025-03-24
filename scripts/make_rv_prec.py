import astropy.table as atpy
import numpy as np
import matplotlib.pyplot as plt
import plot_preamb as pp
import scipy.optimize
import matplotlib.colors as maco
import matplotlib.gridspec as gridspec


def like(p, args):
    zmag, lerv, feh = args
    lpred = p[0] + p[1] * (zmag - 16) + p[2] * feh
    ret = np.abs(lpred - lerv).sum()
    print(ret, p)
    return ret


def fitter(zmag, rv_err, feh):
    lerv = np.log10(rv_err)
    R = scipy.optimize.minimize(like, [0, 0, 0],
                                args=((zmag, lerv, feh), ),
                                method='Nelder-Mead')
    return R.x


pp.run()

RV_T = atpy.Table().read('../data/mwsall-pix-iron.fits',
                         'RVTAB',
                         mask_invalid=False)
SP_T = atpy.Table().read('../data/mwsall-pix-iron.fits',
                         'SPTAB',
                         mask_invalid=False)
FM_T = atpy.Table().read('../data/mwsall-pix-iron.fits',
                         'FIBERMAP',
                         mask_invalid=False)
G_T = atpy.Table().read('../data/mwsall-pix-iron.fits',
                        'GAIA',
                        mask_invalid=False)
SC_T = atpy.Table().read('../data/mwsall-pix-iron.fits',
                         'SCORES',
                         mask_invalid=False)

main_sel = (RV_T['RVS_WARN'] == 0) & (RV_T['RR_SPECTYPE'] == 'STAR')


def betw(x, x1, x2):
    return (x >= x1) & (x < x2)


plt.clf()
fig = plt.figure(figsize=(3.37 * 1, 3.37 * .6))
cnt = 0
cnt = 0
cur_sel = (main_sel & (RV_T['SURVEY'] == 'main') &
           (RV_T['PROGRAM'] == 'bright') & (RV_T['VSINI'] < 30)
           & betw(RV_T['TEFF'], 4500, 7000)
           & betw(SC_T['TSNR2_LRG'] * 12.15, 180, 220))
zmag = 22.5 - 2.5 * np.log10(FM_T['FLUX_Z'])

zgrid = np.linspace(15.5, 19, 100)
fit_sub = ((RV_T['RR_SPECTYPE'] == 'STAR') & (RV_T['RVS_WARN'] == 0) &
           (RV_T['SURVEY'] == 'main') & (RV_T['PROGRAM'] == 'bright') &
           (RV_T['VSINI'] < 30)
           & betw(RV_T['TEFF'], 4500, 7000)
           & betw(SC_T['TSNR2_LRG'] * 12.15, 180, 220) & betw(zmag, 15.5, 19))
coeffs = fitter(zmag[fit_sub], RV_T['VRAD_ERR'][fit_sub], RV_T['FEH'][fit_sub])
coeffs = np.round(coeffs, 2)
mag_mult = coeffs[1]
feh_mult = coeffs[2]
zpt = coeffs[0]
gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 0.03])
for i in range(2):
    ax = fig.add_subplot(gs[0, i])
    # plt.subplot(1, 2, i + 1)
    feh = {0: 0, 1: -2}[i]
    cur_sel1 = cur_sel & betw(RV_T['FEH'] - feh, -0.05, 0.05)
    im = ax.hist2d(zmag[cur_sel1],
                   np.log10(RV_T['VRAD_ERR'][cur_sel1]),
                   range=[[15.5, 19.5], [-.8, 1.2]],
                   bins=[60, 60],
                   norm=maco.PowerNorm(gamma=.5))[-1]
    plt.plot(
        zgrid,
        zpt + feh * feh_mult + (zgrid - 16) * mag_mult,
        color='red',
        linestyle='--',
        label=(r'$y= %.2f \cdot (z-16)$' + '\n' + '$ %g\cdot [Fe/H] + %g$') %
        (mag_mult, feh_mult, zpt))
    # norm=maco.LogNorm()

    plt.xlabel('z [mag]')
    plt.text(17, 1, '[Fe/H]=$%d$' % feh, color='white')
    plt.xlim(15.5, 19.5)

    im.set_rasterized(True)
    if i == 0:
        plt.ylabel(r'$\log_{10} [\sigma_{RV}/(1 {\rm km/s})]$ ')
    else:
        plt.legend(numpoints=1, loc='lower right', fontsize='x-small')
        plt.gca().yaxis.set_major_formatter(plt.NullFormatter())
# plt.colorbar()
cax = fig.add_subplot(gs[0, 2])
fig.colorbar(im, cax=cax)

# plt.subplot(122)
# plt.title('[Fe/H]=-2')
# cur_sel2 = cur_sel & betw(RV_T['FEH'], -2, -1.9)
# plt.xlabel('z [mag]')
plt.tight_layout()
plt.subplots_adjust(wspace=0.03, left=.15, top=.99)

plt.savefig('plots/rv_prec.pdf')
