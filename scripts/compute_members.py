import astropy.table as atpy
import numpy as np
import sqlutilpy as sqlutil
import collections
import feh_correct
import fit_scatter
import multiprocessing as mp


def get_lists():
    PL = atpy.Table().read(
        '/home/skoposov/science/desi/spphot_fit/validation/paceli/allmembs.fits'
    )
    PL = PL[PL['mem_fixed'] > .9]
    PL['name'] = PL['key']
    # PL['source_id'] is already there
    BAU = atpy.Table().read(
        '/home/skoposov/science/datasets/vasiliev_baumgardt2021/cluster_members.fits'
    )
    BAU = BAU[BAU['memberprob'] > .9]
    # name and source_id are already there

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
    CA['source_id'] = CA_dr3_id
    CA['name'] = CA['Cluster']
    IB = atpy.Table().read(
        '/home/skoposov/science/datasets/Streamfinder_Ibata24/members.fits')
    IB2 = atpy.Table().read(
        '/home/skoposov/science/datasets/Streamfinder_Ibata24/streams.fits')
    stream_map = {}
    for i in range(len(IB2)):
        stream_map[IB2['s_ID'][i]] = IB2['Name'][i]
    IB['name'] = ['                                '] * len(IB)
    for i in range(len(IB)):
        IB['name'][i] = stream_map[IB['Stream'][i]]
    IB['source_id'] = IB['Gaia']
    return {'dSph': PL, 'GC': BAU, 'OCL': CA, 'Stream': IB}


if __name__ == '__main__':
    RV_T = atpy.Table().read('../data/mwsall-pix-iron.fits',
                             'RVTAB',
                             mask_invalid=False)
    G_T = atpy.Table().read('../data/mwsall-pix-iron.fits',
                            'GAIA',
                            mask_invalid=False)
    feh_calib = feh_correct.calibrate(np.asarray(RV_T["FEH"]),
                                      np.asarray(RV_T["TEFF"]),
                                      np.asarray(RV_T['LOGG']),
                                      pipeline='RVS',
                                      release='DR1')
    feh_calib = np.asarray(feh_calib)
    feh_err = np.asarray(RV_T['FEH_ERR'])
    lists = get_lists()
    to_print = []
    queue = []
    with mp.Pool(36) as poo:
        for typ, tab in lists.items():
            sub = np.isin(tab['source_id'], G_T['SOURCE_ID'])
            curtab = tab[sub]
            for k, v in collections.Counter(curtab['name']).items():
                cursid = curtab['source_id'][curtab['name'] == k]
                if v > 10:
                    # print(typ, k, v)
                    to_print.append((typ, k, v))
                    sub2 = (np.isin(G_T['SOURCE_ID'], cursid)
                            & np.isfinite(feh_calib + feh_err) &
                            (RV_T['RVS_WARN'] == 0)
                            & (RV_T['PRIMARY']) &
                            (RV_T['RR_SPECTYPE'] == 'STAR'))
                    queue.append(
                        poo.apply_async(fit_scatter.get_scatter,
                                        (feh_calib[sub2], feh_err[sub2])))
                    #ret_mean, ret, warn, _ = fit_scatter.get_scatter(
                    #    feh_calib[sub2], np.asarray(RV_T['FEH_ERR'][sub2]))
        for tp, it in zip(to_print, queue):
            ret_mean, ret, warn, _ = it.get()
            typ, k, v = tp
            print(typ, k, v)
            print('%.2f %.2f %.2f' % tuple(ret_mean),
                  '%.2f %.2f %.2f' % tuple(ret))
