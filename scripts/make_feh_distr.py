import astropy.table as atpy
import matplotlib.pyplot as plt
# import numpy as np
# import matplotlib.colors as maco
import plot_preamb as pp
from config import main_file, data_path, external_path

fname = data_path + '/' + main_file

pp.run()


def betw(x, x1, x2):
    return (x >= x1) & (x < x2)


RV_T = atpy.Table().read(fname, 'RVTAB', mask_invalid=False)
SP_T = atpy.Table().read(fname, 'SPTAB', mask_invalid=False)
FM_T = atpy.Table().read(fname, 'FIBERMAP', mask_invalid=False)
G_T = atpy.Table().read(fname, 'GAIA', mask_invalid=False)

main_sel = (RV_T['RVS_WARN'] == 0) & (RV_T['RR_SPECTYPE']
                                      == 'STAR') & (RV_T['SN_R'] > 10)

plt.clf()
cnt = 0
rv_sel = betw(RV_T['TEFF'], 4500, 7000) & (RV_T['VSINI'] < 30)
sp_sel = (SP_T['BESTGRID'] != 's_rdesi1') & betw(SP_T['TEFF'], 4500, 7000)

cur_sel0 = main_sel & RV_T['PRIMARY'] & (RV_T['SURVEY'] == 'main')
#& (RV_T['PROGRAM']
#                                                    == 'bright')

fig = plt.figure(figsize=(3.37 * 1, 3.37 * .75))

cur_sel = cur_sel0 & rv_sel
opt = dict(alpha=.5, bins=55 * 2, range=[-5, .5], histtype='step')
plt.hist(RV_T['FEH'][cur_sel], label='RVS', **opt)
cur_sel = cur_sel0 & sp_sel
plt.hist(SP_T['FEH'][cur_sel], label='SP', **opt)

plt.plot([-1.4, -1.4], [100, 3e3], linestyle='--', color='grey')
plt.plot([-2.1, -2.1], [100, 3e3], linestyle='--', color='grey')
plt.text(-1.6, 50, 'GSE', color='grey')
plt.gca().set_yscale('log')
plt.xlabel('[Fe/H]')
plt.ylabel('stars/bin')
plt.legend()
plt.tight_layout()
plt.savefig('plots/feh_distr.pdf')

print('main/bright, then all')
for cur_sel0 in [
        main_sel & (RV_T['SURVEY'] == 'main') & (RV_T['PROGRAM'] == 'bright'),
        main_sel & RV_T['PRIMARY']
]:
    cur_sel = cur_sel0 & rv_sel
    print('RV', (cur_sel & (RV_T['FEH'] < -2)).sum())
    print('RV', (cur_sel & (RV_T['FEH'] < -3)).sum())

    cur_sel = cur_sel0 & sp_sel
    plt.hist(SP_T['FEH'][cur_sel], label='SP', **opt)
    print('SP', (cur_sel & (SP_T['FEH'] < -2)).sum())
    print('SP', (cur_sel & (SP_T['FEH'] < -3)).sum())
