#!/usr/bin/env python3
# -*- coding: utf-8 -*-


''' Necessary imports '''
import sys
import numpy as np

''' Dynamic progress bar '''
def progress(count, total, status=''):
    bar_len = 50
    filled_len = int(round(bar_len * count / float(total)))
    percents = round(100. * count / float(total), 1)
    bar = '■' * filled_len + '-' * (bar_len - filled_len)
    fmt = '[%s] %s%s ...%s' % (bar, percents, '%', status)
    print('\b' * len(fmt), end='')  # clears the line
    sys.stdout.write(fmt)
    if count < total:
        sys.stdout.flush()
    else:
        sys.stdout.write('\n')
    return

''' Gaussian instrument profile calculation given target resolution 'res' and
wavelength vector 'wl'. Wavelengths MUST BE already spaced at constant steps in
velocity (or in delta(lambda)/lambda), that is constant resolving power R. 
The code will not check for the spacing, so be aware! '''
def get_gaussian_ip(res, wl):
    lamPix = 0.5*(wl[1:]+wl[0:-1])
    dlamPix = (wl[1:]-wl[0:-1])
    dvPix = np.mean(dlamPix / lamPix)    # Dropping a factor 'c' (in common with dvFWHM below)
    dvFWHM = 1.0 / res
    fwhmPix = dvFWHM / dvPix
    sigmaPix = 0.5 * fwhmPix / np.sqrt(2.0 * np.log(2.0))
    hker = int(4.0 * fwhmPix)
    nker = hker * 2 + 1
    xker = np.arange(nker) - hker
    yker = np.exp(-0.5*(xker/sigmaPix)**2)
    yker /= yker.sum()
    return yker, hker

''' Air to vacuum correction - wavelengths must be inputted in Angstrom '''
def air_to_vac(wAir):
    s2 = (1E4 / wAir)**2.0
    n = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s2) + 0.0001599740894897 / (38.92568793293 - s2)
    return wAir.copy() * n