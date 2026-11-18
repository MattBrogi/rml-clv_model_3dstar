#!/usr/bin/env python
# encoding: utf-8

''' This is a list of system and analysis parameters to be imported in multiple
files '''

import numpy as np

b = 0.267       # Impact parameter
a = 13.9        # Semi-major axis in units of Rstar
rr = 0.0764     # Planet-to-star radius
gRSh = 0.68     # Gravitational red shift for the star (km/s)
ks = 0.00821    # Stellar RV semi-amplitude (km/s)
dph = 0.        # Extra phase shift (optional)
vsini = 1.2     # Stellar projected rotational velocity (from RStar and Prot)
lam_deg = 98.0
lam_rad = lam_deg/180*np.pi  # Obliquity angle
iStar_deg = 90.0
iStar_rad = iStar_deg/180*np.pi  # Stellar inclination angle
vSys = 24.6

''' Output files '''
fullSpecFile = 'integ_flux.npy'
fullContFile = 'integ_cont.npy'
plShadowFile = 'plnt_shadow_gj3470.npy'
finalModelFile = 'spec_series_gj3470.npy'