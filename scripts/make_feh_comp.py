import astropy.table as atpy
import astropy.units as auni
import astropy.coordinates as acoo
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as maco
import plot_preamb as pp
import crossmatcher_cache as crossmatcher
from matplotlib.colors import TABLEAU_COLORS
import scipy.optimize
from config import main_file, data_path, external_path

fname = data_path + '/' + main_file

teff_ref = 5000
logteff_scale = 0.1

# minteff, maxteff = 4500, 7000
minteff, maxteff = 500, 37000


def betw(x, x1, x2):
    return (x >= x1) & (x < x2)


def combiner(*args):
    # fill nansby us
    ret = args[0] * 1
    for x in args[1:]:
        ind = ~np.isfinite(ret)
        ret[ind] = x[ind]
    return ret


def get_saga(ra, dec):
    # Read SAGA data
    SAGAT = atpy.Table().read(external_path + '/saga_cleaned_catalog.tsv',
                              format='ascii')

    # Build Astropy SkyCoord objects
    c_input = acoo.SkyCoord(ra=ra * auni.deg, dec=dec * auni.deg)
    c_saga = acoo.SkyCoord(ra=SAGAT['RAdeg'] * auni.deg,
                           dec=SAGAT['DECdeg'] * auni.deg)

    # Match input coordinates to SAGA catalog
    idx, sep2d, _ = c_input.match_to_catalog_sky(c_saga)

    # We only accept matches within 1 arcsecond (1")
    # 1" = 1/3600 degrees
    match_mask = (sep2d.arcsec < 1.0)

    # Prepare output. Initialize fe_h with NaN.
    SAGA_R = {'fe_h': np.full(len(ra), np.nan, dtype=float)}

    matched_fe_h = SAGAT['[M/H]'].filled(np.nan)[idx]
    SAGA_R['fe_h'][match_mask] = matched_fe_h[match_mask]

    return SAGA_R


def get_ges(ra, dec):
    GEST = atpy.Table().read(external_path + '/gaia_eso_cat_dr4_esoarch.fits',
                             mask_invalid=False)
    sub = (((GEST['SFLAGS'] == '                                       ') |
            (GEST['SFLAGS'] == 'NIA                                    ')) &
           (GEST['E_FEH'] < .1))
    GEST = GEST[sub]
    c_input = acoo.SkyCoord(ra=ra * auni.deg, dec=dec * auni.deg)
    c_ges = acoo.SkyCoord(ra=GEST['RA'], dec=GEST['DECLINATION'])

    # Match input coordinates to SAGA catalog
    idx, sep2d, _ = c_input.match_to_catalog_sky(c_ges)

    # We only accept matches within 1 arcsecond (1")
    # 1" = 1/3600 degrees
    match_mask = (sep2d.arcsec < 1.0)

    # Prepare output. Initialize fe_h with NaN.
    GES_R = {'fe_h': np.full(len(ra), np.nan, dtype=float)}

    matched_fe_h = GEST['FEH'][idx]
    GES_R['fe_h'][match_mask] = matched_fe_h[match_mask]

    return GES_R


def func(p, X, Y):
    return np.mean(
        np.abs(Y - np.poly1d(p)((X - np.log10(teff_ref)) / logteff_scale)))


def fitter(teff, feh_ref, feh_obs):
    X, Y = np.log10(teff), feh_obs - feh_ref
    aind = main_sel & np.isfinite(X + Y) & betw(teff, minteff, maxteff)
    X, Y = X[aind], Y[aind]
    R1 = scipy.optimize.minimize(func, [0, 0, 0], args=(X, Y))
    return R1.x


def func_cut(p, X, Y, Z, cut):
    P_low = np.poly1d(p[:3])
    P_high = np.poly1d(p[3:])
    xscale = (X - np.log10(teff_ref)) / logteff_scale
    pred_low = P_low(xscale)
    pred_high = P_high(xscale)
    pred = pred_low * (Y < cut) + pred_high * (Y >= cut)
    return np.mean(np.abs(Z - pred))


def fitter_cut(teff, logg, feh_ref, feh_obs, sel0=None):
    X, Y, Z = np.log10(teff), logg, feh_obs - feh_ref
    aind = sel0 & np.isfinite(X + Y + Z) & betw(teff, minteff, maxteff)
    X, Y, Z = np.array(X[aind]), np.array(Y[aind]), np.array(Z[aind])
    print('Number of stars', aind.sum())
    perc_val = 0.2, 99.8  # 10 stars of 5000
    print('Percentiles')
    for q in ['teff', 'logg']:
        print(q, end=' ')
        arr = {'teff': teff, 'logg': logg}[q][aind]
        for curp in perc_val:
            print(np.round(np.percentile(arr, curp), 1), end=' ')
        print(' ')
    funcs = []
    cuts = np.arange(0, 6, 0.1)
    for cut in cuts:
        R1 = scipy.optimize.minimize(func_cut,
                                     np.zeros(6),
                                     args=(X, Y, Z, cut),
                                     method='Nelder-Mead')
        R1 = scipy.optimize.minimize(func_cut, R1.x, args=(X, Y, Z, cut))
        funcs.append(R1.fun)
    best_cut = cuts[np.argmin(funcs)]
    R1 = scipy.optimize.minimize(func_cut,
                                 np.zeros(6),
                                 args=(X, Y, Z, best_cut))
    return best_cut, R1.x


pp.run()

RV_T = atpy.Table().read(fname, 'RVTAB', mask_invalid=False)
SP_T = atpy.Table().read(fname, 'SPTAB', mask_invalid=False)
FM_T = atpy.Table().read(fname, 'FIBERMAP', mask_invalid=False)
G_T = atpy.Table().read(fname, 'GAIA', mask_invalid=False)

main_sel0 = (RV_T['RVS_WARN'] == 0) & (RV_T['RR_SPECTYPE'] == 'STAR')
main_sel = main_sel0 & (RV_T['SURVEY'] == 'main') & (RV_T['SN_R'] > 10)
rv_sel = (RV_T['FEH_ERR'] < 0.1) & (RV_T['VSINI'] < 30)
sp_sel = ((SP_T['BESTGRID'] != 's_rdesi1') & (SP_T['COVAR'][:, 0, 0]**.5 < .1))

ra, dec = RV_T['TARGET_RA'], RV_T['TARGET_DEC']

#  fetch surveys
D_SAGA = get_saga(ra, dec)

D_GES = get_ges(RV_T['TARGET_RA'], RV_T['TARGET_DEC'])
D_GA = crossmatcher.doit_by_key(
    'galah_dr4.allstar',
    G_T['SOURCE_ID'],
    'fe_h,teff,logg,mg_fe,ca_fe,c_fe,flag_fe_h,flag_sp',
    db='wsdb',
    asDict=True,
    key_col='gaiadr3_source_id')

D_AP = crossmatcher.doit_by_key(
    'apogee_dr17.allstar',
    G_T['SOURCE_ID'],
    '''alpha_m,fe_h,c_fe,n_fe,o_fe,na_fe,mg_fe,si_fe,ca_fe,ti_fe,mn_fe,ni_fe,ce_fe,vhelio_avg,logg,teff,teff_spec,logg_spec,
            aspcapflag,starflag, fe_h_flag''',
    key_col='gaiaedr3_source_id',
    db='wsdb',
    asDict=True)

D_GAIA = crossmatcher.doit_by_key(
    'gaia_dr3.astrophysical_parameters',
    G_T['SOURCE_ID'],
    '''mh_gspspec,teff_gspspec,fem_gspspec,logg_gspspec,
        coalesce(flags_gspspec like '0000000000000%', false) as good_flag ''',
    key_col='source_id',
    db='wsdb',
    asDict=True)

D_GAIA['fe_h'] = D_GAIA['mh_gspspec']
D_GAIA['fe_h'][~D_GAIA['good_flag']] = np.nan
D_AP['fe_h'][(D_AP['fe_h_flag'] != 0) | (D_AP['starflag'] != 0) |
             (D_AP['aspcapflag'] != 0)] = np.nan
D_GA['fe_h'][(D_GA['flag_fe_h'] != 0) | (D_GA['flag_sp'] != 0)] = np.nan

print('---------')
print('SP')
comb_feh = combiner(D_GA["fe_h"], D_AP['fe_h'])

cut_sp, coeff_sp = fitter_cut(SP_T['TEFF'],
                              SP_T['LOGG'],
                              comb_feh,
                              SP_T['FEH'],
                              sel0=main_sel & sp_sel)
print('---------')
print('RV')
cut_rv, coeff_rv = fitter_cut(RV_T['TEFF'],
                              RV_T['LOGG'],
                              comb_feh,
                              RV_T['FEH'],
                              sel0=main_sel & rv_sel)
print('---------')
coeff_sp = np.round(coeff_sp, 3)
coeff_rv = np.round(coeff_rv, 3)

coeff_rv_low = coeff_rv[:3][::-1]
coeff_rv_high = coeff_rv[3:][::-1]
coeff_sp_low = coeff_sp[:3][::-1]
coeff_sp_high = coeff_sp[3:][::-1]
print('RV', coeff_rv_low, '\n', coeff_rv_high, '\n cut', cut_rv, '\n', 'SP',
      coeff_sp_low, '\n', coeff_sp_high, '\n cut', cut_sp)

fp = open('output/feh_coeffs.txt', 'w')
print('RV',
      'low', ("%s %s %s %s") % ((cut_rv, ) + tuple(coeff_rv_low)),
      file=fp)
print('RV',
      'high', ("%s %s %s %s") % ((cut_rv, ) + tuple(coeff_rv_high)),
      file=fp)

print('SP',
      'low', ("%s %s %s %s") % ((cut_sp, ) + tuple(coeff_sp_low)),
      file=fp)
print('SP',
      'high', ("%s %s %s %s") % ((cut_sp, ) + tuple(coeff_sp_high)),
      file=fp)
fp.close()
STR = f'''
    cuts = {{'RVS': {cut_rv}, 'SP': {cut_sp}}}
    coeffs_low = {{'RVS': {np.array2string(coeff_rv_low, separator=',')},
'SP':  {np.array2string(coeff_sp_low, separator=',')} }}
    coeffs_high = {{'RVS': {np.array2string(coeff_rv_high, separator=',')},
 'SP':  {np.array2string(coeff_sp_high, separator=',')} }}
'''
print(STR)
print('-------------')

# RVS, low
print(f"    RVS & $\\log g <{cut_rv:.1f}$  "
      f"& ${coeff_rv_low[0]:.3f}$ "
      f"& ${coeff_rv_low[1]:.3f}$ "
      f"& ${coeff_rv_low[2]:.3f}$ \\\\")

# RVS, high
print(f"    RVS & $\\log g \\geq {cut_rv:.1f}$ "
      f"& ${coeff_rv_high[0]:.3f}$ "
      f"& ${coeff_rv_high[1]:.3f}$ "
      f"& ${coeff_rv_high[2]:.3f}$ \\\\")

# SP, low
print(f"    SP  & $\\log g <{cut_sp:.1f}$ "
      f"& ${coeff_sp_low[0]:.3f}$ "
      f"& ${coeff_sp_low[1]:.3f}$ "
      f"& ${coeff_sp_low[2]:.3f}$ \\\\")

# SP, high
print(f"    SP  & $\\log g \\geq {cut_sp:.1f}$ "
      f"& ${coeff_sp_high[0]:.3f}$ "
      f"& ${coeff_sp_high[1]:.3f}$ "
      f"& ${coeff_sp_high[2]:.3f}$ \\\\")
print('-----')
import feh_correct  # noqa

SP_T['FEH_CALIB2'] = (SP_T['FEH'] - np.poly1d(coeff_sp_low[::-1])
                      (np.log10(SP_T['TEFF'] / teff_ref) / logteff_scale) *
                      (SP_T['LOGG'] < cut_sp) - np.poly1d(coeff_sp_high[::-1])
                      (np.log10(SP_T['TEFF'] / teff_ref) / logteff_scale) *
                      (SP_T['LOGG'] >= cut_sp))

RV_T['FEH_CALIB2'] = (RV_T['FEH'] - np.poly1d(coeff_rv_low[::-1])
                      (np.log10(RV_T['TEFF'] / teff_ref) / logteff_scale) *
                      (RV_T['LOGG'] < cut_rv) - np.poly1d(coeff_rv_high[::-1])
                      (np.log10(RV_T['TEFF'] / teff_ref) / logteff_scale) *
                      (RV_T['LOGG'] >= cut_rv))

xf1 = feh_correct.calibrate(RV_T['FEH'],
                            RV_T['TEFF'],
                            RV_T['LOGG'],
                            pipeline='RVS')
xf2 = feh_correct.calibrate(SP_T['FEH'],
                            SP_T['TEFF'],
                            SP_T['LOGG'],
                            pipeline='SP')


def checker(x1, x2):
    xind = np.isfinite(x1 + x2)
    assert (np.allclose(x1[xind], x2[xind]))


if True:
    checker(xf1, RV_T['FEH_CALIB2'])
    checker(xf2, SP_T['FEH_CALIB2'])
else:
    print('!! NOT CHECKING feh_correct !!! ')

plt.clf()
fig = plt.figure(figsize=(3.37 * 1, 3.37 * 1.4))
pad = 10
gs = plt.GridSpec(300 + pad, 200)
cnt = 0

for x1, pipe in enumerate(['RVS', 'SP']):
    for x2 in range(3):
        curT = {'RVS': RV_T, 'SP': SP_T}[pipe]
        titl = pipe
        if pipe == 'RVS':
            cur_sel = main_sel & rv_sel
        else:
            cur_sel = main_sel & sp_sel
        cur_sel = cur_sel & betw(curT['TEFF'], minteff, maxteff)
        # plt.subplot(3, 2, 1 + 2 * x2 + x1)
        plt.subplot(gs[100 * x2 + pad * (x2 == 2):100 * x2 + 100 + pad *
                       (x2 == 2), 100 * x1:100 * x1 + 100])
        COMP = [D_GA, D_AP, D_SAGA][x2]
        curfeh = curT["FEH"]
        if x2 < 2:
            plt.hist2d(COMP['fe_h'][cur_sel],
                       curfeh[cur_sel],
                       range=[[-2.5, .6], [-2.5, .6]],
                       bins=[50, 50],
                       norm=maco.PowerNorm(gamma=0.5, vmax=100))
            for xx in 'xy':
                plt.gca().tick_params(axis=xx,
                                      which='both',
                                      labelcolor='black',
                                      color='white')
            for xx in plt.gca().spines.keys():
                plt.gca().spines[xx].set_color("white")
            plt.gca().tick_params(top=False, which='both')
            plt.ylim(-2.49, .6)
            plt.gci().set_rasterized(True)
            if x1 == 0:
                plt.text(-2., 0.2, ['GALAH', 'APOGEE'][x2], color='white')
        else:
            plt.plot(COMP['fe_h'][cur_sel], curfeh[cur_sel], '.')
            plt.xlim(-4, -.01)
            plt.ylim(-4, -.01)
            if x1 == 0:
                plt.text(-3, -.5, 'SAGA')
        plt.plot([-4, 1], [-4, 1], color='red', linestyle='--', dashes=(3, 10))

        plt.xlabel('[Fe/H]$_{survey}$ [dex]')

        # plt.title(titl)
        if x2 == 0:
            plt.annotate(titl, (0.5, .93),
                         xycoords='axes fraction',
                         color='white',
                         horizontalalignment='center')
        if x1 == 0:
            plt.ylabel(r'[Fe/H]$_{DESI}$ [dex]')
        else:
            plt.gca().yaxis.set_major_formatter(plt.NullFormatter())
        cnt += 1
plt.subplots_adjust(top=.99, right=.99, left=.13, bottom=.065)
plt.savefig('plots/feh_compar.pdf')

for var_name in ['FEH', 'FEH_CALIB2']:
    fig = plt.figure(figsize=(3.37 * 1, 3.37 * 1))
    ncnt = 4
    for cnt in range(ncnt):
        plt.subplot(ncnt, 1, cnt + 1)
        COMP = [D_GA, D_AP, D_GAIA, D_SAGA][(cnt)]
        tit = ['GALAH', 'APOGEE', 'Gaia', 'SAGA'][cnt]
        bins = [100, 100, 100, 10][cnt]
        for i, (curT, label) in enumerate(zip([RV_T, SP_T], ['RVS', 'SP'])):
            if i == 0:
                cur_sel = main_sel & rv_sel
            else:
                cur_sel = main_sel & sp_sel
            curfeh = curT[var_name]

            cur_sel = cur_sel & betw(curT['TEFF'], minteff,
                                     maxteff)  # & (T['FEH_ERR'] < .1)
            delt = curfeh[cur_sel] - COMP['fe_h'][cur_sel]
            delt = delt[np.isfinite(delt)]
            percs = [np.percentile(delt, _) for _ in [16, 50, 84]]
            med = percs[1]

            plt.hist(
                delt,
                range=[-1, 1],
                # linestyle=[':', '--'][i],
                histtype='step',
                bins=bins,
                label=label,
            )
            plt.annotate(r'$%.2f_{%.2f}^{+%.2f}$' %
                         (percs[1], percs[0] - percs[1], percs[2] - percs[1]),
                         (.8, .8 - .3 * i),
                         color=list(TABLEAU_COLORS.values())[i],
                         xycoords='axes fraction')
        plt.annotate(tit, (.5, .86),
                     xycoords='axes fraction',
                     horizontalalignment='center')
        postfix = {'FEH': '', 'FEH_CALIB2': r'$_{\rm calibrated}$'}[var_name]
        if cnt == ncnt - 1:
            plt.xlabel(r'$\delta$ [Fe/H]' + postfix + ' [dex]')
        else:
            plt.gca().xaxis.set_major_formatter(plt.NullFormatter())
        plt.ylim(.1, plt.ylim()[1] * 1.15)
        plt.gca().tick_params(top=False, which='both')

        if cnt == 0:
            plt.legend()
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.0)
    plt.savefig({
        'FEH': 'plots/feh_compar_delta.pdf',
        'FEH_CALIB2': 'plots/feh_compar_delta_calibrated.pdf'
    }[var_name])

if False:
    # disable
    plt.clf()
    print('temporary validation plot')
    # this is internal validation
    from idlplotInd import tvhist2d  # noqa
    tvhist2d(SP_T['LOGG'],
             comb_feh - SP_T['FEH_CALIB2'],
             -1,
             6,
             -1,
             1,
             normx='sum',
             vmax=.2,
             bins=[40, 40],
             subplot=221,
             ind=sp_sel,
             title='SP',
             ytitle='delta feh')
    tvhist2d(RV_T['LOGG'],
             comb_feh - RV_T['FEH_CALIB2'],
             -1,
             6,
             -1,
             1,
             normx='sum',
             vmax=.2,
             bins=[40, 40],
             subplot=223,
             ind=rv_sel,
             title='RVS',
             xtitle='logg')

    tvhist2d(SP_T['TEFF'],
             comb_feh - SP_T['FEH_CALIB2'],
             3500,
             8000,
             -1,
             1,
             normx='sum',
             vmax=.2,
             ind=sp_sel,
             subplot=222,
             bins=[40, 40])
    tvhist2d(RV_T['TEFF'],
             comb_feh - RV_T['FEH_CALIB2'],
             3500,
             8000,
             -1,
             1,
             normx='sum',
             vmax=.2,
             subplot=224,
             bins=[40, 40],
             ind=rv_sel,
             xtitle='teff')
    plt.savefig('tmp_plots/valid.pdf')
