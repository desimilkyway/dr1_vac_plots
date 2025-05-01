import astropy.table as atpy
import numpy as np
import matplotlib.pyplot as plt
import plot_preamb as pp
import scipy.optimize
import scipy.stats
import healpy
from config import main_file, data_path
import sky_plotter

fname = data_path + '/' + main_file

pp.run()

RV_T = atpy.Table().read(fname, 'RVTAB', mask_invalid=False)
SP_T = atpy.Table().read(fname, 'SPTAB', mask_invalid=False)
FM_T = atpy.Table().read(fname, 'FIBERMAP', mask_invalid=False)
G_T = atpy.Table().read(fname, 'GAIA', mask_invalid=False)

main_sel = (RV_T['RVS_WARN'] == 0) & (RV_T['RR_SPECTYPE'] == 'STAR')
cur_sel0 = main_sel & (RV_T['SURVEY'] == 'main') & (RV_T['PROGRAM']
                                                    == 'backup')


def betw(x, x1, x2):
    return (x >= x1) & (x < x2)


ra, dec = RV_T['TARGET_RA'], RV_T['TARGET_DEC']

plt.clf()
fig = plt.figure(figsize=(3.37 * 2., 3.37 * 1))
delt = RV_T['VRAD'] - G_T['RADIAL_VELOCITY']
cur_sel = np.isfinite(delt) & cur_sel0

xr = 0, 360
yr = -32, 90
bins = [360 // 2, 120 // 2]
S = scipy.stats.binned_statistic_2d(RV_T['TARGET_RA'][cur_sel],
                                    RV_T['TARGET_DEC'][cur_sel],
                                    delt[cur_sel],
                                    'median',
                                    range=[xr, yr],
                                    bins=bins)
SC = scipy.stats.binned_statistic_2d(RV_T['TARGET_RA'][cur_sel],
                                     RV_T['TARGET_DEC'][cur_sel],
                                     delt[cur_sel],
                                     'count',
                                     range=[xr, yr],
                                     bins=bins)
minval = 5
stat = S.statistic
stat[SC.statistic < minval] = np.nan

nside = 32

hpx = healpy.ang2pix(nside,
                     RV_T['TARGET_RA'],
                     RV_T['TARGET_DEC'],
                     lonlat=True,
                     nest=True)

S = scipy.stats.binned_statistic(hpx[cur_sel],
                                 delt[cur_sel],
                                 'median',
                                 range=[-.5, 12 * nside**2 - .5],
                                 bins=12 * nside**2)
SC = scipy.stats.binned_statistic(hpx[cur_sel],
                                  delt[cur_sel],
                                  'count',
                                  range=[-.5, 12 * nside**2 - 0.5],
                                  bins=12 * nside**2)
minval = 5
stat = S.statistic
stat[SC.statistic < minval] = np.nan

print('Total', np.isfinite(stat).sum())
print('<20', (np.abs(stat) < 20).sum())

R = sky_plotter.hpx_show(
    stat,
    ra_shift=180,
    vmin=-20,
    vmax=20,
    cmap='turbo',
    dra_label=60,
)
R.set_rasterized(True)
plt.colorbar(R, label=r'$\delta V_{rad}$ [km/s]', shrink=.8)

if False:
    plt.imshow(stat.T,
               extent=list(xr) + list(yr),
               origin='lower',
               cmap='turbo',
               aspect='auto',
               vmax=25,
               vmin=-25)
plt.subplots_adjust(right=.999, left=0.04, top=.98)
plt.savefig('plots/backup_delt.pdf')
