import astropy.table as atpy
import astropy.io.fits as pyfits
import match_lists
import crossmatcher
import numpy as np

GT = pyfits.getdata('../data/mwsall-pix-iron.fits', 'GAIA')
FT = pyfits.getdata('../data/mwsall-pix-iron.fits', 'FIBERMAP')
TT = atpy.Table().read('../data/mwsall-pix-iron.fits', 'RVTAB')
xind = (TT['PROGRAM'] == 'backup') & (TT['SURVEY'] == 'main')
TT = TT[xind]
GT = GT[xind]
tiles = atpy.Table().read('tiles-iron.fits')

tiles = tiles[(tiles['PROGRAM'] == 'backup') & (tiles['SURVEY'] == 'main')]
DD, xind = match_lists.match_lists(TT['TARGET_RA'], TT['TARGET_DEC'],
                                   tiles['TILERA'], tiles['TILEDEC'], 100)
tileid = tiles[xind]['TILEID']

D_AP = crossmatcher.doit('apogee_dr17.allstar',
                         TT['TARGET_RA'],
                         TT['TARGET_DEC'],
                         '''vhelio_avg,
    aspcapflag,starflag, fe_h_flag''',
                         db='wsdb',
                         asDict=True)
D_GA = crossmatcher.doit('galah_dr4.allstar',
                         TT['TARGET_RA'],
                         TT['TARGET_DEC'],
                         'fe_h,teff,logg,fe_h,mg_fe,ca_fe,rv_comp_1',
                         host='cappc127.ast.cam.ac.uk',
                         db='wsdb',
                         asDict=True)
corr = TT['VRAD'] * 0
res = {}
for t in np.unique(tileid):
    xind = tileid == t
    res[t] = np.nanmedian((TT['VRAD'] - GT['RADIAL_VELOCITY'])[xind])
    corr[xind] = res[t] + np.zeros(xind.sum())
    # TT['VRAD'][xind] - res[t]
atpy.Table({
    'TARGETID': TT['TARGETID'],
    'VRAD_BIAS': corr,
}).write('backup_correction.fits', overwrite=True)
