import astropy.table as atpy
import numpy as np
import matplotlib.pyplot as plt
# import matplotlib.colors as maco
import matplotlib.colors as colors
import plot_preamb as pp
from config import main_file, data_path, external_path

fname = data_path + '/' + main_file

pp.run()

RV_T = atpy.Table().read(fname, 'RVTAB', mask_invalid=False)
FM_T = atpy.Table().read(fname, 'FIBERMAP', mask_invalid=False)
G_T = atpy.Table().read(fname, 'GAIA', mask_invalid=False)

main_sel = (RV_T['RVS_WARN'] == 0) & (RV_T['RR_SPECTYPE'] == 'STAR')

fig = plt.figure(figsize=(3.37 * 1, 3.37 * 1.))
cnt = 0
bitmasks = [
    ('MWS_MAIN_BLUE', 8),
    ('MWS_MAIN_RED', 11),
    ('MWS_BROAD', 0),
    ('MWS_WD', 1),
    ('MWS_NEARBY', 2),
    ('MWS_BHB', 6),
]

survey, program = ('main', 'bright')
for curt, bit in bitmasks:
    objtype_sel = (FM_T['MWS_TARGET'] & (2**bit)) > 0
    cur_sel = objtype_sel & main_sel & (RV_T['SURVEY'] == survey) & (
        RV_T['PROGRAM'] == program)
    plt.subplot(2, 3, cnt + 1)
    vmax = {
        'MWS_MAIN_BLUE': 2000,
        'MWS_MAIN_RED': 2000,
        'MWS_BROAD': 2000,
        'MWS_WD': 50,
        'MWS_NEARBY': 50,
        'MWS_BHB': 50
    }
    plt.hist2d(-2.5 *
               np.log10(FM_T['FLUX_G'][cur_sel] / FM_T['FLUX_R'][cur_sel]),
               22.5 - 2.5 * np.log10(FM_T['FLUX_R'][cur_sel]),
               bins=[100, 100],
               range=[[-0.49, 2.1], [15.31, 20.5]],
               norm=colors.PowerNorm(gamma=0.5, vmax=vmax[curt]))
    plt.gci().set_rasterized(True)
    if cnt % 3 == 0:
        plt.ylabel('r [mag]')
    else:
        plt.gca().yaxis.set_major_formatter(plt.NullFormatter())
    if cnt % 3 == 2:
        plt.colorbar(shrink=.9)
    cnt += 1
    # plt.title(f'survey, program: {survey},{program}')
    plt.ylim(20.5, 15.31)
    plt.text(-.3, 16., f'{curt}\n bitmask {2**bit}', color='white')
    if cnt > 3:
        plt.xlabel('g-r [mag]')
        plt.xticks([0, 1, 2])
    else:
        plt.gca().xaxis.set_major_formatter(plt.NullFormatter())
plt.tight_layout()
plt.subplots_adjust(wspace=0., hspace=0.01, top=.99, right=.94, left=.096)

plt.savefig('plots/targ_classes.pdf')
