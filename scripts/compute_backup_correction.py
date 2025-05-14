import astropy.table as atpy
import astropy.io.fits as pyfits
import crossmatcher_cache as crossmatcher
import numpy as np
import astropy.coordinates as acoo
import astropy.units as u

from config import main_file, data_path, external_path

fname = data_path + '/' + main_file
GT = pyfits.getdata(fname, 'GAIA')
FT = pyfits.getdata(fname, 'FIBERMAP')
TT = atpy.Table().read(fname, 'RVTAB')
xind = (TT['PROGRAM'] == 'backup') & (TT['SURVEY'] == 'main')
TT = TT[xind]
GT = GT[xind]
tiles = atpy.Table().read(external_path + '/tiles-iron.fits')

tiles = tiles[(tiles['PROGRAM'] == 'backup') & (tiles['SURVEY'] == 'main')]
coord1 = acoo.SkyCoord(ra=TT['TARGET_RA'] * u.deg,
                       dec=TT['TARGET_DEC'] * u.deg)
coord2 = acoo.SkyCoord(ra=tiles['TILERA'] * u.deg,
                       dec=tiles['TILEDEC'] * u.deg)

# Match each source in coord1 to nearest source in coord2
idx, sep2d, _ = coord1.match_to_catalog_sky(coord2)
tileid = tiles[idx]['TILEID']

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
                         'teff,logg,fe_h,mg_fe,ca_fe,rv_comp_1',
                         db='wsdb',
                         asDict=True)
corr = TT['VRAD'] * 0
res = {}
for t in np.unique(tileid):
    xind = tileid == t
    res[t] = np.nanmedian((TT['VRAD'] - GT['RADIAL_VELOCITY'])[xind])
    corr[xind] = res[t] + np.zeros(xind.sum())
    # TT['VRAD'][xind] - res[t]
print(np.percentile(corr, [1, 16, 50, 84, 99]))
atpy.Table({
    'TARGETID': TT['TARGETID'],
    'VRAD_BIAS': corr,
}).write('output/backup_correction.fits', overwrite=True)
