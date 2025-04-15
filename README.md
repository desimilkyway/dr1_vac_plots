The scripts used to make figures are in scripts folder
you need to put DES DR1 files in the data folder.

You may need to change the paths in scripts/config.py to point at
the DESI data and several external tables

The code also relies on cached results of some queries stored in
query_cache/
and external datasets stored in
external/


Required python modules :
 numpy scipy
 matplotlib
 astropy 
 dynesty
 duckdb
 sqlutilpy

If on nersc you need this module
texlive/2024


After that you should be able to run all the .py
scripts in the scripts/ folder to make all the figures