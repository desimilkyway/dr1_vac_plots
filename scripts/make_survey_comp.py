# import astropy.table as atpy
import astropy.io.fits as pyfits
import sqlutil_cache as sqlutil
import matplotlib.pyplot as plt
import plot_preamb as pp
from matplotlib.colors import TABLEAU_COLORS
from config import main_file, data_path
from matplotlib.colors import to_rgba

fname = data_path + '/' + main_file

pp.run()
T = pyfits.getdata(fname, 'RVTAB')
SP_T = pyfits.getdata(fname, 'SPTAB')
lam_t = sqlutil.get('select gaia_g_mean_mag from lamost_dr9.lrs_stellar')
gaia_t = sqlutil.get('select phot_g_mean_mag from gaia_dr3_aux.gaia_source_rv')
GT = pyfits.getdata(fname, 'GAIA')
# T2 = pyfits.getdata('../../loa/rvpix-loa.fits', 'RVTAB')
# GT2 = pyfits.getdata('../../loa/rvpix-loa.fits', 'GAIA')
sdss_t = sqlutil.get(
    '''select  phot_g_mean_mag from gaia_dr3.gaia_source as g ,
gaia_edr3_aux.sdssdr13bestneighbour as x,
sdssdr14.sppparams as sp, sdssdr14.specphotoall as s
where x.original_ext_source_id =s.objid and
    g.source_id=x.source_id and (elodiervfinal between -1000 and 1000) and
s.scienceprimary=1 and
s.specobjid=sp.specobjid and fehadop > -6 and fehadop <10 ;''')

min_mag = 10
max_mag = 22
bins = (max_mag - min_mag) * 10
shift = 1 / 7.
kw = dict(histtype='stepfilled', range=[min_mag, max_mag], bins=bins)
plt.figure(figsize=(3.37 * 2, 2.5))
colors = list(TABLEAU_COLORS.values())

cnt = 0
progs = ['dark', 'bright', 'backup']
color_map = {'dark': colors[0], 'bright': colors[1], 'backup': colors[2]}
for program in progs:
    sub = (T['RVS_WARN'] == 0) & (T['RR_SPECTYPE'] == 'STAR') & (
        T['PRIMARY']) & (T['SURVEY'] == 'main') & (T['PROGRAM'] == program)
    # sub2 = (T2['RVS_WARN'] == 0) & (T2['RR_SPECTYPE'] == 'STAR') & (
    #    T2['PRIMARY']) & (T2['SURVEY'] == 'main') & (T2['PROGRAM'] == program)
    # if False:
    # plt.hist(GT2['PHOT_G_MEAN_MAG'][sub2],
    #         color='lightgrey',
    #         label='DESI DR2',
    #         linestyle=ls,
    #         **kw)
    plt.hist(
        GT['PHOT_G_MEAN_MAG'][sub],
        edgecolor=color_map[program],  # ('black'),
        #             linewidth=2,
        # linestyle=ls,
        #label='DESI DR1',
        color=to_rgba(colors[cnt], alpha=.1),
        **kw)
    kw['range'] = [kw['range'][0] + shift, kw['range'][1] + shift]
    cnt += 1

plt.text(13, 3000, 'DESI DR1 backup', rotation=10, color=color_map['backup'])
plt.text(18, 1.3e5, 'DESI DR1 bright', color=color_map['bright'])
plt.text(16.5, 2000, 'DESI DR1 dark', rotation=7, color=color_map['dark'])

kw['histtype'] = 'step'
ls = {'LAMOST': '--', 'Gaia': None, 'SDSS': ':'}
plt.hist(lam_t,
         color='black',
         label='LAMOST LRS DR9',
         linestyle=ls['LAMOST'],
         **kw)
kw['range'] = [kw['range'][0] + shift, kw['range'][1] + shift]
plt.hist(sdss_t, color='black', label='SDSS DR14', linestyle=ls['SDSS'], **kw)
kw['range'] = [kw['range'][0] + shift, kw['range'][1] + shift]
plt.hist(gaia_t,
         color='black',
         label='Gaia DR3 RVS',
         linestyle=ls['Gaia'],
         **kw)
kw['range'] = [kw['range'][0] + shift, kw['range'][1] + shift]
plt.xlim(11, 21)
plt.ylim(10, 2e6)
# plt.text(18, 3e5, 'DESI DR2', color='lightgrey')
plt.text(12, 3e4, 'LAMOST DR9', color='black', rotation=10)
plt.text(14.5, 2e3, 'SDSS DR14', color='black', rotation=20)
plt.text(12, 1.5e5, 'Gaia DR3 RVS', color='black', rotation=10)
plt.xlabel('G [mag]')
plt.gca().set_yscale('log')
# plt.xlim(15.8, 20.5)
plt.ylabel('N$_{stars}$/bin')
# plt.legend()
plt.tight_layout()
plt.savefig('plots/survey_comparison.pdf')
