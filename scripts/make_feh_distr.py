import astropy.table as atpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as maco
import plot_preamb as pp

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

main_sel = (RV_T['RVS_WARN'] == 0) & (RV_T['RR_SPECTYPE'] == 'STAR')

plt.clf()
cnt = 0
cur_sel0 = main_sel & (RV_T['SURVEY'] == 'main') & (
    RV_T['PROGRAM'] == 'bright') & (RV_T['SN_R'] > 10)


def betw(x, x1, x2):
    return (x >= x1) & (x < x2)


fig = plt.figure(figsize=(3.37 * 1, 3.37 * .75))

cur_sel = cur_sel0 & betw(RV_T['TEFF'], 4500, 7000) & (RV_T['VSINI'] < 30)
opt = dict(alpha=.5, bins=100, range=[-5, .5], histtype='step')
plt.hist(RV_T['FEH'][cur_sel], label='RVS', **opt)
print('RV', (cur_sel & (SP_T['FEH'] < -3)).sum())

cur_sel = cur_sel0 & (SP_T['BESTGRID'] != 's_rdesi1') & betw(
    SP_T['TEFF'], 4500, 7000)
plt.hist(SP_T['FEH'][cur_sel], label='SP', **opt)
print('SP', (cur_sel & (SP_T['FEH'] < -3)).sum())

plt.gca().set_yscale('log')
plt.xlabel('[Fe/H]')
plt.legend()
plt.tight_layout()
plt.savefig('plots/feh_distr.pdf')
