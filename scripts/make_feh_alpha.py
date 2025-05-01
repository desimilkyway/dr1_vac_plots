import astropy.table as atpy
import matplotlib.pyplot as plt
import matplotlib.colors as maco
import matplotlib.gridspec as gridspec
import matplotlib.patches as mapa
import plot_preamb as pp
import numpy as np
from config import main_file, data_path

fname = data_path + '/' + main_file

pp.run()

RV_T = atpy.Table().read(fname, 'RVTAB', mask_invalid=False)
SP_T = atpy.Table().read(fname, 'SPTAB', mask_invalid=False)
FM_T = atpy.Table().read(fname, 'FIBERMAP', mask_invalid=False)
G_T = atpy.Table().read(fname, 'GAIA', mask_invalid=False)

main_sel = (RV_T['RVS_WARN'] == 0) & (RV_T['RR_SPECTYPE'] == 'STAR')

plt.clf()
fig = plt.figure(figsize=(3.37 * 1, 3.37 * .85))
cnt = 0
cnt = 0
cur_sel0 = (
    main_sel & (RV_T['SURVEY'] == 'main') &
    #            (RV_T['PROGRAM'] == 'bright') &
    (RV_T['SN_R'] > 10))


def betw(x, x1, x2):
    return (x >= x1) & (x < x2)


gs = gridspec.GridSpec(nrows=2,
                       ncols=3,
                       width_ratios=[1, 1, 0.05],
                       wspace=0.05,
                       hspace=.13,
                       bottom=.1)
for xcnt in range(2):
    for cnt in range(2):
        T = [RV_T, SP_T][cnt]
        if cnt == 0:
            cur_sel = cur_sel0 & (RV_T['VSINI'] < 30)
        else:
            cur_sel = cur_sel0 & (SP_T['BESTGRID'] != 's_rdesi1')
        cur_sel = cur_sel & betw(T['TEFF'], 4500, 7000)
        # plt.subplot(1, 2, cnt + 1)
        plt.subplot(gs[xcnt, cnt])
        range_dict = {0: [[-5, 1], [-.5, 1.3]], 1: [[-2, .5], [-.2, .7]]}
        rr = range_dict[xcnt]
        if xcnt == 0:
            rr1 = np.array(range_dict[1])
        vmax = {0: 10000, 1: 3000}[xcnt]
        im = plt.hist2d((T['FEH'][cur_sel]),
                        T['ALPHAFE'][cur_sel],
                        range=rr,
                        bins=[100, 100],
                        norm=maco.PowerNorm(gamma=.5, vmax=vmax))
        if xcnt == 0:

            rect = mapa.Rectangle((rr1[0, 0], rr1[1, 0]),
                                  rr1[0, 1] - rr1[0, 0],
                                  rr1[1, 1] - rr1[1, 0],
                                  linewidth=1,
                                  edgecolor='grey',
                                  linestyle=':',
                                  facecolor='none')
            plt.gca().add_patch(rect)
        plt.gci().set_rasterized(True)
        if xcnt == 1:
            plt.xlabel('[Fe/H] ')
        for xx in 'xy':
            plt.gca().tick_params(axis=xx,
                                  which='both',
                                  labelcolor='black',
                                  color='white')
        for xx in plt.gca().spines.keys():
            plt.gca().spines[xx].set_color("white")
        plt.gca().tick_params(top=False, which='both')

        # plt.title(['RVS', 'SP'][cnt])
        if xcnt == 0:
            plt.annotate(['RVS', 'SP'][cnt], (.5, .9),
                         color='white',
                         xycoords='axes fraction',
                         horizontalalignment='center')

        if cnt == 0:
            plt.ylabel(r'[$\alpha$/Fe]')
        else:
            plt.gca().yaxis.set_major_formatter(plt.NullFormatter())
    ax = plt.subplot(gs[xcnt, 2])
    plt.colorbar(im[-1], cax=ax)

# plt.tight_layout()
plt.subplots_adjust(wspace=0., hspace=0.001, top=.97, bottom=.16)

plt.savefig('plots/feh_alpha.pdf')
