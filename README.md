This is the repository of codes used to make all the figures/tables in the Koposov+2025 paper 
on DESI DR1 stellar catalogue.


The scripts used to make figures are in scripts folder

The code also relies on cached results of some queries stored in
query_cache/
Available here: https://doi.org/10.5281/zenodo.15317970

and external datasets stored in
external/
Available here: https://doi.org/10.5281/zenodo.15317970

You also need to have DESI VAC tables.

You may need to change the paths in scripts/config.py to point at
the DESI data.

You will need the
mwsall-pix-iron.fits fits file
from https://data.desi.lbl.gov/public/dr1/vac/dr1/mws/iron/v1.0/
And the files for single epoch measurements:
rvpix_exp-sv1-bright.fits
rvpix_exp-cmx-other.fits       rvpix_exp-sv1-dark.fits
rvpix_exp-main-backup.fits     rvpix_exp-sv1-other.fits
rvpix_exp-main-bright.fits     rvpix_exp-sv2-backup.fits
rvpix_exp-main-dark.fits       rvpix_exp-sv2-bright.fits
rvpix_exp-special-backup.fits  rvpix_exp-sv2-dark.fits
rvpix_exp-special-bright.fits  rvpix_exp-sv3-backup.fits
rvpix_exp-special-dark.fits    rvpix_exp-sv3-bright.fits
rvpix_exp-special-other.fits   rvpix_exp-sv3-dark.fits
rvpix_exp-sv1-backup.fits


Required python modules :
 numpy
 scipy
 matplotlib
 astropy 
 dynesty
 duckdb
 healpy
 pandas

If on nersc you need this module
texlive/2024


After that you should be able to run all the .py
scripts in the scripts/ folder to make all the figures.
You only need to run the make_*scripts to create figures
