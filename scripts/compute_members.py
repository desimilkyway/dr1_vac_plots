import astropy.table as atpy
import numpy as np
import sqlutilpy as sqlutil
import collections

PL = atpy.Table().read(
    '/home/skoposov/science/desi/spphot_fit/validation/paceli/allmembs.fits')
PL = PL[PL['mem_fixed'] > .9]
BAU = atpy.Table().read(
    '/home/skoposov/science/datasets/vasiliev_baumgardt2021/cluster_members.fits'
)
BAU = BAU[BAU['memberprob'] > .9]
CA = atpy.Table().read(
    '/home/skoposov/science/datasets/cantat_gaugin_ocmem_2020.fits')
CA = CA[CA['proba'] > .9]
rid = np.arange(len(CA))
xid = CA['GaiaDR2']
CA_dr3_id, = sqlutil.local_join(
    '''
with x as (select rid, dr3_source_id, row_number() over (
    partition by rid order by angular_distance asc ) from gaia_edr3.dr2_neighbourhood as g, m
where m.source_id = g.dr2_source_id )
select dr3_source_id from x where row_number =1 order by rid
''', 'm', (
        rid,
        xid,
    ), (
        'rid',
        'source_id',
    ))

IB = atpy.Table().read(
    '/home/skoposov/science/datasets/Streamfinder_Ibata24/members.fits')
IB2 = atpy.Table().read(
    '/home/skoposov/science/datasets/Streamfinder_Ibata24/streams.fits')
stream_map = {}
for i in range(len(IB2)):
    stream_map[IB2['s_ID'][i]] = IB2['Name'][i]

RV_T = atpy.Table().read('../data/mwsall-pix-iron.fits',
                         'RVTAB',
                         mask_invalid=False)
G_T = atpy.Table().read('../data/mwsall-pix-iron.fits',
                        'GAIA',
                        mask_invalid=False)

sub = np.isin(PL['source_id'], G_T['SOURCE_ID'])
PL = PL[sub]
for k, v in collections.Counter(PL['key']).items():
    if v > 10:
        print('dSph', k, v)

sub = np.isin(BAU['source_id'], G_T['SOURCE_ID'])
BAU = BAU[sub]
for k, v in collections.Counter(BAU['name']).items():
    if v > 10:
        print('GC', k, v)

sub = np.isin(CA_dr3_id, G_T['SOURCE_ID'])
CA = CA[sub]
for k, v in collections.Counter(CA['Cluster']).items():
    if v > 10:
        print('OCL', k, v)

sub = np.isin(IB['Gaia'], G_T['SOURCE_ID'])
IB = IB[sub]
for k, v in collections.Counter(IB['Stream']).items():
    if v > 10:
        print('Stream', stream_map[k], v)
