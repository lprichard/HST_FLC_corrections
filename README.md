# HST FLC corrections

Codes to download data from the Mikulski Archive for Space Telescopes (MAST) with astroquery and apply additional calibrations to reduced Hubble Space Telescope (HST) single-visit images (called FLCs). Different corrections are needed for the Wide Field Camera 3 (WFC3)/UV-Visible (UVIS) and Advanced Camera for Surveys (ACS) FLCs. See HST_FLC_corrections.ipynb for a full list of corrections with examples and a step-by-step guide to applying them to your data.

Written by Laura Prichard, May 2021. Includes codes developed Ben Sunnquist (on [GitHub](https://github.com/bsunnquist/uvis-skydarks)) & Marc Rafelski.

Please reference [Prichard et al. 2022](https://ui.adsabs.harvard.edu/abs/2022ApJ...924...14P/abstract) and codes by Ben Sunnquist if you use any of the corrections outlined here.

Notebook tested with: Python v3.7, astropy v4.0, astroquery v0.4, drizzlepac 3.1.8, photutils v1.0.2, and stwcs v1.6.1. These versions and higher are recommended.
