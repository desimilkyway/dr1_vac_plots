# import astropy.table as atpy
import astropy.io.fits as pyfits
import sqlutilpy as sqlutil
import matplotlib.pyplot as plt
import plot_preamb as pp
from matplotlib.colors import TABLEAU_COLORS

pp.run()
fname = '../data/mwsall-pix-iron.fits'
T = pyfits.getdata('../data/mwsall-pix-iron.fits', 'RVTAB')
SP_T = pyfits.getdata('../data/mwsall-pix-iron.fits', 'SPTAB')
lam_t = sqlutil.get('select gaia_g_mean_mag from lamost_dr9.lrs_stellar')
gaia_t = sqlutil.get('select phot_g_mean_mag from gaia_dr3_aux.gaia_source_rv')
GT = pyfits.getdata('../data/mwsall-pix-iron.fits', 'GAIA')
T2 = pyfits.getdata('../../loa/rvpix-loa.fits', 'RVTAB')
GT2 = pyfits.getdata('../../loa/rvpix-loa.fits', 'GAIA')
sdss_t = sqlutil.get(
    '''select  phot_g_mean_mag from gaia_dr3.gaia_source as g ,
gaia_edr3_aux.sdssdr13bestneighbour as x,
sdssdr14.sppparams as sp, sdssdr14.specphotoall as s
where x.original_ext_source_id =s.objid and
    g.source_id=x.source_id and (elodiervfinal between -1000 and 1000) and
s.scienceprimary=1 and
s.specobjid=sp.specobjid and fehadop > -6 and fehadop <10 ;''')
sub = (T['RVS_WARN'] == 0) & (T['RR_SPECTYPE'] == 'STAR') & (T['PRIMARY']) & (
    T['SURVEY'] == 'main') & (T['PROGRAM'] == 'bright')
sub2 = (T2['RVS_WARN'] == 0) & (T2['RR_SPECTYPE'] == 'STAR') & (
    T2['PRIMARY']) & (T2['SURVEY'] == 'main') & (T2['PROGRAM'] == 'bright')
kw = dict(histtype='step', range=[15, 21], bins=60)
plt.figure(figsize=(3.37, 2.5))

plt.hist(GT2['PHOT_G_MEAN_MAG'][sub2],
         color='lightgrey',
         label='DESI DR2',
         **kw)
colors = list(TABLEAU_COLORS.values())
plt.hist(GT['PHOT_G_MEAN_MAG'][sub],
         color='black',
         linewidth=2,
         label='DESI DR1',
         **kw)
plt.hist(lam_t, color=colors[1], label='LAMOST LRS DR9', **kw)
plt.hist(sdss_t, color=colors[2], label='SDSS DR14', **kw)
plt.hist(gaia_t, color=colors[3], label='Gaia DR3 RVS', **kw)
plt.text(18, 3e5, 'DESI DR2', color='lightgrey')
plt.text(18, .4e5, 'DESI DR1')
plt.text(18.4, 70, 'LAMOST DR9', color=colors[1])
plt.text(18, 3.9e3, 'SDSS DR14', color=colors[2])
plt.text(16.7, 100, 'Gaia DR3 RVS', color=colors[3])

plt.xlabel('G [mag]')
plt.gca().set_yscale('log')
plt.xlim(15.8, 20.5)
plt.ylim(10, 7e5)
# plt.legend()
plt.tight_layout()
plt.savefig('plots/survey_comparison.pdf')
