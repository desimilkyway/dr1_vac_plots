import astropy.table as atpy
import numpy as np
import sqlutilpy as sqlutil
import collections
import feh_correct
import fit_scatter
import multiprocessing as mp
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt


def get_lists():
    PL = atpy.Table().read(
        '/home/skoposov/science/desi/spphot_fit/validation/paceli/'
        'allmembs.fits')
    PL = PL[PL['mem_fixed'] > .9]
    PL['name'] = PL['key']
    # PL['source_id'] is already there
    BAU = atpy.Table().read(
        '/home/skoposov/science/datasets/vasiliev_baumgardt2021/'
        'cluster_members.fits')
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
        partition by rid order by angular_distance asc )
        from gaia_edr3.dr2_neighbourhood as g, m
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
    fname = '../data/mwsall-pix-iron.fits'
    RV_T = atpy.Table().read(fname, 'RVTAB', mask_invalid=False)
    SP_T = atpy.Table().read(fname, 'SPTAB', mask_invalid=False)
    G_T = atpy.Table().read(fname, 'GAIA', mask_invalid=False)

    calibrated = True
    lists = get_lists()

    to_print = []
    queue = []
    out_tab = {}
    with mp.Pool(36) as poo:
        for pipeline in ['RVS', 'SP']:
            if pipeline == 'SP':
                XTAB = SP_T
                bad = (SP_T['BESTGRID'] == 's_rdesi1')
                SP_T['FEH'][bad] = np.nan
                feh_err = np.asarray(SP_T['COVAR'][:, 0, 0]**.5)
            elif pipeline == 'RVS':
                XTAB = RV_T
                feh_err = np.asarray(RV_T['FEH_ERR'])
            else:
                raise Exception('oops')
            if calibrated:
                feh = feh_correct.calibrate(np.asarray(XTAB["FEH"]),
                                            np.asarray(XTAB["TEFF"]),
                                            np.asarray(XTAB['LOGG']),
                                            pipeline=pipeline,
                                            release='DR1')
            else:
                feh = np.asarray(XTAB['FEH'])

            feh = np.asarray(feh)

            for typ, tab in lists.items():
                sub = np.isin(tab['source_id'], G_T['SOURCE_ID'])
                curtab = tab[sub]
                for k, v in collections.Counter(curtab['name']).items():
                    cursid = curtab['source_id'][curtab['name'] == k]
                    if v > 10:
                        # print(typ, k, v)
                        to_print.append((typ, k, v, pipeline))
                        sub2 = (np.isin(G_T['SOURCE_ID'], cursid)
                                & np.isfinite(feh + feh_err) &
                                (RV_T['RVS_WARN'] == 0)
                                & (RV_T['PRIMARY']) &
                                (RV_T['RR_SPECTYPE'] == 'STAR'))
                        queue.append(
                            (poo.apply_async(fit_scatter.get_scatter,
                                             (feh[sub2], feh_err[sub2])),
                             feh[sub2], feh_err[sub2]))

        pdf_file = PdfPages('tmp_plots/multipage_plots.pdf')
        for tp, it in zip(to_print, queue):
            (ret_mean, ret_sig, warn,
             _), cur_feh, cur_efeh = it[0].get(), it[1], it[2]
            # manually kill LMS-1
            typ, k, v, pipeline = tp
            if len(cur_feh) < 10 or k == 'LMS-1' or k == 'Leiptr':
                ret_mean = [np.nan] * 3
                ret_sig = [np.nan] * 3
            if k not in out_tab:
                out_tab[k] = (typ, k, v, {'RVS': [ret_mean, ret_sig]})
            else:
                out_tab[k][3]['SP'] = [ret_mean, ret_sig]
            fig = plt.figure()
            plt.hist(cur_feh, range=[-4, .5], bins=45, histtype='step')
            plt.axvline(ret_mean[1], color='red')
            plt.title(k + ' ' + pipeline)
            plt.xlabel('feh')
            pdf_file.savefig()
            plt.close()
        pdf_file.close()
        with open('output/objs.txt', 'w') as fp:
            print('type name count', end=' ', file=fp)
            for pp in ['rvs', 'sp']:
                end = {'rvs': ' ', 'sp': None}[pp]
                print(
                    f'feh_{pp}_1 feh_{pp}_2 feh_{pp}_3 '
                    f'sfeh_{pp}_1 sfeh_{pp}_2 sfeh_{pp}_3',
                    end=end,
                    file=fp)
            for k in out_tab.keys():
                typ, k, v, param = out_tab[k]
                print(typ, k, v, end=" ", file=fp)
                ret_mean, ret_sig = param['RVS']
                print('%.2f %.2f %.2f' % tuple(ret_mean),
                      '%.2f %.2f %.2f' % tuple(ret_sig),
                      end=" ",
                      file=fp)
                ret_mean, ret_sig = param['SP']
                print('%.2f %.2f %.2f' % tuple(ret_mean),
                      '%.2f %.2f %.2f' % tuple(ret_sig),
                      file=fp)
