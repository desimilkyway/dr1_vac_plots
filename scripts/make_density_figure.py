import astropy.table as atpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as maco
import plot_preamb as pp
import astropy.units as auni
import astropy.coordinates as acoo
from config import main_file, data_path

fname = data_path + '/' + main_file

pp.run()

RV_T = atpy.Table().read(fname, 'RVTAB', mask_invalid=False)

main_sel = (RV_T['RVS_WARN'] == 0) & (RV_T['RR_SPECTYPE'] == 'STAR')

plt.clf()
fig = plt.figure(figsize=(3.37 * 2, 3.37 * .65))
cnt = 0
shift = 60
lgrid = np.linspace(0, 360, 1000)

icrs_coord = acoo.SkyCoord(l=lgrid * auni.deg,
                           b=lgrid * 0 * auni.deg,
                           frame='galactic')
xra = icrs_coord.icrs.ra.deg
xdec = icrs_coord.icrs.dec.deg
xra = (xra + shift) % 360 - shift
xind = np.argsort(xra)
xra, xdec = xra[xind], xdec[xind]
print(xra[0], xra[-1])
for survey, program in [('main', 'bright'), ('main', 'backup'),
                        ('main', 'dark')]:
    cur_sel = main_sel & (RV_T['SURVEY'] == survey) & (RV_T['PROGRAM']
                                                       == program)
    plt.subplot(1, 3, cnt + 1)
    plt.hist2d((RV_T['TARGET_RA'][cur_sel] + shift) % 360 - shift,
               RV_T['TARGET_DEC'][cur_sel],
               bins=[360, 123],
               range=[[-shift, 360 - shift], [-30, 93]],
               weights=1. / np.cos(np.deg2rad(RV_T['TARGET_DEC'][cur_sel])),
               vmax=700)
    plt.gci().set_rasterized(True)
    plt.plot(xra, xdec, linestyle='--', color='grey')
    if cnt > 0:
        plt.gca().yaxis.set_major_formatter(plt.NullFormatter())
    else:
        plt.ylabel(r'$\delta$ [deg]')
    # else:
    plt.xlabel(r'$\alpha$ [deg]')
    cnt += 1
    #plt.title(f'survey, program: {survey},{program}')
    plt.text(150, 86, f'{survey},{program}', color='white')
    plt.xlim(360 - shift, -shift)
    for xx in 'xy':
        plt.gca().tick_params(axis=xx,
                              which='both',
                              labelcolor='black',
                              color='white')
    for xx in plt.gca().spines.keys():
        plt.gca().spines[xx].set_color("white")
    plt.gca().tick_params(top=False, which='both')
    if cnt == 3:
        plt.colorbar()
plt.tight_layout()
plt.subplots_adjust(wspace=0.02, hspace=0.01, top=.98)

plt.savefig('plots/density.pdf')
