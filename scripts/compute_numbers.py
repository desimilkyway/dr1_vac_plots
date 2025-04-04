import astropy.table as atpy
import numpy as np
from config import main_file, data_path, external_path

fname = data_path + '/' + main_file

RV_T = atpy.Table().read(fname, 'RVTAB', mask_invalid=False)
# Colour-magnitude distribution by survey/program
surveys = np.unique(RV_T['SURVEY'][:])
mysurveys = ['cmx', 'special', 'sv1', 'sv2', 'sv3', 'main']
assert (tuple(sorted(mysurveys)) == tuple(surveys))
for s in mysurveys:
    xind = RV_T['SURVEY'] == s
    for p in np.unique(RV_T['PROGRAM'][xind]):
        xind1 = (RV_T['SURVEY'] == s) & (RV_T['PROGRAM'] == p)
        xind2 = (RV_T['SURVEY'] == s) & (RV_T['PROGRAM'] == p) & (
            RV_T['RR_SPECTYPE'] == 'STAR')
        xind3 = (RV_T['SURVEY'] == s) & (RV_T['PROGRAM']
                                         == p) & (RV_T['RVS_WARN'] == 0)
        print('%s & %s &  %s & %s & %s \\\\' %
              (s, p, format(xind1.sum(), ','), format(
                  xind2.sum(), ','), format(xind3.sum(), ',')))

print(np.median(RV_T['VRAD_ERR']))
sub = (RV_T['SURVEY'] == 'main') & (RV_T['PROGRAM'] == 'bright')
print(np.median(RV_T['VRAD_ERR'][sub]))
sub = ((RV_T['SURVEY'] == 'main') & (RV_T['PROGRAM'] == 'bright') &
       (RV_T['RR_SPECTYPE'] == 'STAR') & (RV_T['RVS_WARN'] == 0)
       & RV_T['PRIMARY'])
SP_T = atpy.Table().read(fname, 'SPTAB', mask_invalid=False)
sub = ((RV_T['RR_SPECTYPE'] == 'STAR') & (RV_T['RVS_WARN'] == 0)
       & RV_T['PRIMARY'])
print('RV MP', (sub & (RV_T['FEH'] < -2)).sum())
print('SP MP', (sub & (SP_T['FEH'] < -2)
                & (SP_T['BESTGRID'] != 's_rdesi1')).sum())
print('RV MP', (sub & (RV_T['FEH'] < -3)).sum())
print('SP MP', (sub & (SP_T['FEH'] < -3)
                & (SP_T['BESTGRID'] != 's_rdesi1')).sum())
print('RV MP', (sub & (RV_T['FEH'] < -4)).sum())
print('SP MP', (sub & (SP_T['FEH'] < -4)
                & (SP_T['BESTGRID'] != 's_rdesi1')).sum())

sub = ((RV_T['RR_SPECTYPE'] == 'STAR') & (RV_T['RVS_WARN'] == 0))
rv = RV_T['VRAD_ERR'][sub]
print([np.percentile(rv, _) for _ in [1, 50, 99]])
sub = ((RV_T['RR_SPECTYPE'] == 'STAR') & (RV_T['RVS_WARN'] == 0) &
       (RV_T['SURVEY'] == 'main') & (RV_T['PROGRAM'] == 'backup'))
rv = RV_T['VRAD_ERR'][sub]
print([np.percentile(rv, _) for _ in [1, 50, 99]])
