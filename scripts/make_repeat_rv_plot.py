import glob
import astropy.table as atpy
import duckdb
import sqlutilpy as sqlutil
import scipy.stats
import matplotlib.pyplot as plt
import numpy as np
import plot_preamb as pp
from config import main_file, data_path, rvexp_path

fname = data_path + '/' + main_file


def func(p, args):
    xerr, yerr = args
    err_calib = np.sqrt(p[0]**2 * xerr**2 + p[1]**2)
    ret = np.sum(np.abs(np.log10(yerr) - np.log10(err_calib)))
    # print(ret, p)
    return ret


def fitter(xerr, yerr):
    args = ((xerr, yerr), )
    R = scipy.optimize.minimize(func, [1, 1], args=args, method='Nelder-Mead')
    R = scipy.optimize.minimize(func, R.x, args=args)
    return R.x


pp.run()
fs = glob.glob('/*exp*fits')
mask = rvexp_path + '/*exp*fits'
fs = glob.glob(mask)
if len(fs) == 0:
    raise Exception(f'failed to find files in the mask {mask}')

tabs = []
for f in fs:
    T = atpy.Table().read(f, 'RVTAB', mask_invalid=False)
    T1 = atpy.Table().read(f, 'FIBERMAP', mask_invalid=False)
    T2 = atpy.Table().read(f, 'GAIA', mask_invalid=False)
    T.convert_bytestring_to_unicode()
    T1.convert_bytestring_to_unicode()
    T2.convert_bytestring_to_unicode()
    D = {}
    XF = f.replace('.fits', '').split('-')
    survey = XF[-2]
    program = XF[-1]
    D['survey'] = [survey] * len(T)
    D['program'] = [program] * len(T)
    for k in [
            'TARGETID', 'VRAD', 'VRAD_ERR', 'RVS_WARN', 'FEH', 'LOGG', 'TEFF',
            'TARGET_RA', 'TARGET_DEC', 'RR_SPECTYPE', 'SN_R'
    ]:
        D[k] = T[k]
    D['MJD'] = T1['MJD']
    D['TILEID'] = T1['TILEID']
    tabs.append(atpy.Table(D))
tabs = atpy.vstack(tabs)
tab_p = tabs.to_pandas()
conn1 = duckdb.connect(':memory:')

Tstack = atpy.Table().read(fname, mask_invalid=False)
Tstack2 = atpy.Table().read(fname, 'FIBERMAP', mask_invalid=False)
Tstack3 = atpy.Table().read(fname, 'GAIA', mask_invalid=False)
Tstack.convert_bytestring_to_unicode()
Tstack2.convert_bytestring_to_unicode()
Tstack3.convert_bytestring_to_unicode()
for k in [
        'SCND_TARGET', 'MWS_TARGET', 'SV1_MWS_TARGET', 'SV2_MWS_TARGET',
        'SV3_MWS_TARGET', 'SCND_TARGET', 'SV1_SCND_TARGET', 'SV1_SCND_TARGET',
        'SV1_SCND_TARGET', 'FLUX_G', 'FLUX_R'
]:
    try:
        Tstack[k] = Tstack2[k].filled(0)
    except AttributeError:
        Tstack[k] = Tstack2[k]
Tstack['source_id'] = Tstack3['SOURCE_ID'].filled(-1)
for k in [
        'PMRA', 'PMDEC', 'PARALLAX', 'PARALLAX_ERROR', 'RUWE',
        'PHOT_G_MEAN_MAG', 'PHOT_BP_MEAN_MAG', 'PHOT_RP_MEAN_MAG', 'DEC_ERROR',
        'ASTROMETRIC_N_GOOD_OBS_AL'
]:
    Tstack[k] = Tstack3[k]
Tstack = Tstack.to_pandas()

conn1.register('repeats0', tab_p)
conn1.register('stack0', Tstack)

conn1.execute('create table repeats as select * from repeats0')
conn1.execute('create table stack as select * from stack0')
conn1.execute('create index qq1 on stack(source_id)')
conn1.execute('create index qq2 on stack(targetid)')

PAIRS = sqlutil.get('''
select r1.mjd as mjd1,
r2.mjd as mjd2,
r1.vrad as vrad1,
r2.vrad as vrad2,
r1.vrad_err as vrad_err1,
r2.vrad_err as vrad_err2,
r1.survey as survey,
r1.program as program,
r1.rvs_warn as rvs_warn1,
r2.rvs_warn as rvs_warn2,
from
repeats as r1,
repeats as r2
where r1.targetid = r2.targetid
and r1.survey=r2.survey
and r1.program=r2.program
and r1.mjd<r2.mjd
and r1.rr_spectype='STAR'
''',
                    conn=conn1,
                    driver='duckdb',
                    asDict=True)

comb_err, delt = ((PAIRS['vrad_err1']**2 + PAIRS['vrad_err2']**2) /
                  2)**.5, (PAIRS['vrad1'] - PAIRS['vrad2']) / np.sqrt(2)

xgrid = np.linspace(-1, 1.5, 100)


def funcer(x):
    return (scipy.stats.scoreatpercentile(x, 84) -
            scipy.stats.scoreatpercentile(x, 16)) / 2.


plt.clf()
survey = 'main'
bins = 20
fig = plt.figure(figsize=(3.37 * 1, 3.37 * .9))
xres2 = {}
xres1 = {}
floor_dict = {}
for i, prog in enumerate(['bright', 'backup', 'dark']):
    sel1 = (PAIRS['program'] == prog) & (PAIRS['survey'] == survey) & (
        PAIRS['rvs_warn1'] == 0) & (PAIRS['rvs_warn2'] == 0)
    erange = [-1, 1.5]
    sel2 = sel1 & (np.abs(PAIRS['mjd1'] - PAIRS['mjd2']) > 1)
    SS1 = scipy.stats.binned_statistic(np.log10(comb_err[sel1]),
                                       delt[sel1],
                                       funcer,
                                       range=erange,
                                       bins=bins)
    SC1 = scipy.stats.binned_statistic(np.log10(comb_err[sel1]),
                                       delt[sel1],
                                       'count',
                                       range=erange,
                                       bins=bins)
    SS2 = scipy.stats.binned_statistic(np.log10(comb_err[sel2]),
                                       delt[sel2],
                                       funcer,
                                       range=erange,
                                       bins=bins)
    SC2 = scipy.stats.binned_statistic(np.log10(comb_err[sel2]),
                                       delt[sel2],
                                       'count',
                                       range=erange,
                                       bins=bins)

    plt.subplot(3, 1, i + 1)
    xsub1 = SC1.statistic > 100
    xsub2 = SC2.statistic > 100
    A1, B1 = (10**(SS1.bin_edges[:-1] + .5 * np.diff(SS1.bin_edges)),
              SS1.statistic)
    A2, B2 = (10**(SS2.bin_edges[:-1] + .5 * np.diff(SS2.bin_edges)),
              SS2.statistic)
    xres1[prog] = A1[xsub1], B1[xsub1]
    xres2[prog] = A2[xsub2], B2[xsub2]
    lab1, lab2 = None, None
    if prog == 'bright':
        lab2 = r'Data, |$\Delta$ t$_{obs}$|>1 d'
        lab1 = 'Data, all'
    plt.plot(A1[xsub1], B1[xsub1], '.', color='grey', label=lab1)
    plt.plot(A2[xsub2], B2[xsub2], '.', color='black', label=lab2)
    plt.ylim(.5, 30)
    plt.xlim(.05, 40)
    plt.gca().set_yscale('log')
    plt.gca().set_xscale('log')

    if i == 1:
        plt.ylabel(r'$\frac{1}{\sqrt{2}}$ StdDev($V_{i}-V_{j}$) [km/s]')
    # else:
    # plt.gca().yaxis.set_major_formatter(plt.NullFormatter())
    if i == 2:
        plt.xlabel(r'$\sqrt{\frac{\sigma_{i}^2+\sigma_{j}^2}{2}}$ [km/s]')
    # plt.text(.2, 10, f'{survey},{prog}')
    plt.text(1.2, 15, f'{survey},{prog}')
    coeffs1 = fitter(*xres1[prog])
    print(prog, 'all', np.round(coeffs1, 2))
    coeffs2 = fitter(*xres2[prog])
    print(prog, '>1d', np.round(coeffs2, 2))
    floor_dict[prog] = np.round(coeffs2[-1], 1)
    floor = floor_dict[prog]
    plt.plot(10**xgrid,
             np.sqrt(10**(2 * xgrid) + floor**2),
             label='Model, floor %.1f km/s' % floor,
             linewidth=1,
             color='black')
    plt.gca().tick_params(top=False, which='both')
    plt.legend(loc='upper left', handlelength=1)
plt.tight_layout()
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig('plots/repeats_accuracy.pdf')

plt.clf()
fig = plt.figure(figsize=(3.37 * 1, 3.37 * .8))
for i, prog in enumerate(['bright', 'backup', 'dark']):
    sel1 = (PAIRS['program'] == prog) & (PAIRS['survey'] == survey) & (
        PAIRS['rvs_warn1'] == 0) & (PAIRS['rvs_warn2'] == 0) & (
            np.abs(PAIRS['mjd1'] - PAIRS['mjd2']) > 1)

    plt.subplot(3, 1, i + 1)
    floor = floor_dict[prog]
    xdelt = delt / (comb_err**2 + floor**2)**.5
    xr = [-6, 6]
    nbins = 100
    plt.hist(xdelt[sel1], range=xr, bins=nbins, histtype='step', color='black')

    xgrid = np.linspace(-6, 6, 1000)
    plt.plot(xgrid,
             scipy.stats.norm(0, 1).pdf(xgrid) * sel1.sum() * (xr[1] - xr[0]) /
             nbins,
             color='red',
             label=r'${\mathcal N}(0,1)$')
    plt.text(.7, .9, f'{survey},{prog}', transform=plt.gca().transAxes)

    if i == 0:
        plt.legend()
    plt.ylabel(r'N/bin')
    if i < 2:
        plt.gca().xaxis.set_major_formatter(plt.NullFormatter())
    if i == 2:
        plt.xlabel(r'$\delta_{RV}/\sigma_{RV,calib}$')
    plt.gca().tick_params(top=False, which='both')
    # plt.legend()
plt.tight_layout()
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig('plots/repeats_delta.pdf')
