# Modelling of time- and wavelength-resolved Rossiter-McLaughlin (RML) and Center-to-Limb Variations (CLVs)
Code and example files for the RML + CLV modelling described in [Chiavassa et al. (2026)](https://scixplorer.org/abs/2026arXiv261105740C/abstract). It uses 3D stellar models (provided separately) as input and it solves the transit geometry of a planet to produce a spectral sequence in the format (wavelength, time).

## Main features

- The code can take in input *any* stellar model resolved in $\mu = 1/\cos(\theta)$ where $\theta$ is the angle between the LOS and the stellar normal. In this version - the whole point of the paper - the model comes from post-processing 3D hydrodynamical simulations of stellar convection with the Stagger grid [Stagger grid](https://scixplorer.org/abs/2013A%26A...557A..26M/abstract, https://scixplorer.org/abs/2024A%26A...688A.212R/abstract), but it can be substituted by any spatially-resolved model.
- The code accounts for instrumental resolution (Gaussian profile) and stellar $v\sin\,i$, but not for the RV corrections related to the reference frame (BERV), systemic velocity, or stellar reflex motion around the center of mass. Note that gravitational redshift e stellar RV semi-amplitude are listed in `parameters.py` but so far unused.
- With respect to [Chiavassa & Brogi (2019)](https://scixplorer.org/abs/2019A%26A...631A.100C/abstract), this version takes into account explicitly the obliquity ($\lambda$) and the angle between the stellar spin and the plane of the sky ($i_\star$) thus allowing modelling of systems with non-zero obliquity.

## Content of the repository

This repository hosts the main `create_rml+clv_model.py` file that is responsible for the calculation of the model, and takes in input:
* A file containing the system parameters (`parameters.py` in the example provided)
* A file containing the $\mu$-resolved stellar spectrum (`3d_stellar_spec.sav` in the example provided)

Additional functions are kept separated for convenience in a `pipeline.py` file.

## Usage

Run `create_rml+clv_model.py`. If `doFullDisk` is set to `True` at line 16 it recomputes the full disk spectrum (slower - Just redo it if basic system parameters such as obliquity angles or radius ratios change). If `doPlShadow` is set to true at line 17 it recomputes the intersection between the planet and the disk of the star.

### Harcoded information in `create_rml+clv_model.py`:
* The exact values of $\mu$ (43 angles in the case of the stellar model provided in the example)
* The way the (3D) stellar models are loaded and read in (this will arguably vary on a case-by-case fashion)
* Resolution of the grid (200 pixel per stellar diameter, matching the precision tests explained in [Chiavassa et al. (2026)](https://scixplorer.org/abs/2026arXiv261105740C/abstract) => this can be increased if higher precision is required).
* Vector of orbital phases (in the example manually defined, in reality it should be the actual phases of the observations)

### Outputs of the code

* Disk-integrated stellar spectrum => 1D array (`integ_flux.npy` in the example provided)
* Disk-integrated continuum spectrum (used to normalize the previous spectrum) => 1D array (`cont_flux.npy` in the example provided)
* Sequence of the stellar strip occluded by the planet => 3D array (time, y coordinate, x coordinate), `plnt_shadow_gj3470.npy` in the example provided
* Final model of the combined RML + CLVs => 2D array (time, wavelength), `spec_series_gj3470.npy` in the example provided.
