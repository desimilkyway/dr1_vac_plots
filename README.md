The scripts used to make figures are in scripts folder
you need to put DES DR1 files in the data folder.

You may need to change the paths in config.py to point at
the DESI data and several external tables

The code also relies on cached results of some queries stores in
scripts/cache


Required python modules :
 numpy scipy
 matplotlib
 astropy 
 dynesty
 duckdb
 sqlutilpy

If on nersc you need this module
texlive/2024