import numpy as np
import matplotlib.pyplot as plt
from scipy.io import readsav
from scipy.interpolate import splrep, splev, interp1d
import sys
sys.path.insert(0,'/Users/mbrogi/work/spectra/pipeline/')
import pipeline as pipe

''' Complete revamp of the CLV+RME model based on Andrea Chiavassa's 3D models '''

thisFold = './'

''' *** PARAMETERS TO BE MODIFIED FOR EACH RUN *** '''

import parameters as pars  # Imports model_parameters.py containing the parameters for the system to simulate
doFullDisk = False   # If True, it recomputes the full-disk stellar spectrum
doPlShadow = True    # If True, it recomputes the planet shadow
instrRes = 1E5  # Instrumental resolving power
outFold = thisFold+'outputs/'   # Output folder for model products
vEq = pars.vsini / np.sin(pars.iStar_rad)

''' Loading in spatially-resolved spectra. This may sometimes be a .sav file,
sometimes a .dat file. PLEASE EDIT ACCORDINGLY to your needs. '''

modFold = 'stellar_model/'
fileName = modFold + 'gj3470_3d_43angles.sav'

data = readsav(fileName)
fMod = np.array(data['data']['flux'][0]).copy()   # Stellar fluxes
cMod = np.array(data['data']['cont'][0]).copy()   # Continuum fluxes
wMod, = pipe.air_to_vac(data['data']['lambda'].copy()) # Air => Vacuum conversion
wMod /= 10.0                                      # Convert wlen to nm
muv = np.array([0.1045, 0.1386, 0.1727, 0.2068, 0.2409, 0.2749, 0.3090, 0.3408,
                0.3727, 0.4045, 0.4363, 0.4682, 0.5000, 0.5282, 0.5564, 0.5846, 
                0.6128, 0.6409, 0.6691, 0.6924, 0.7158, 0.7391, 0.7624, 
                0.7857, 0.8090, 0.8264, 0.8439, 0.8613, 0.8787, 0.8961, 
                0.9135, 0.9243, 0.9351, 0.9458, 0.9566, 0.9674, 0.9781, 
                0.9818, 0.9854, 0.9891, 0.9927, 0.9964, 1.0000])

print('Wavelength range: {:4.2f}-{:4.2f} nm'.format(wMod[0],wMod[-1]))

''' Applying convolution with the Instrumental Profile and pre-computing the
spline interpolation of fluxes and continuum to save time. '''
modResFull = 0.5*(wMod[1:]+wMod[0:-1])/(wMod[1:]-wMod[0:-1])
modRes = np.mean(modResFull)
print('Model equivalent resolution: {:6.1f}'.format(modRes))
ker, hk = pipe.get_gaussian_ip(instrRes,wMod)
plt.plot(ker,'+')
plt.title('Convolution kernel for stellar model spectra')
plt.show()
for i in range(len(muv)):
    fMod[i,] = np.convolve(fMod[i,],ker,mode='same')
    cMod[i,] = np.convolve(cMod[i,],ker,mode='same')
fMod = fMod[:,hk:-hk]
cMod = cMod[:,hk:-hk]
wMod = wMod[hk:-hk]
nWlen = len(wMod)
np.save(outFold+'wlmodel.npy',wMod)
''' Precompute the linear interpolation of the model along the wlen axis '''
fModLinInt = interp1d(muv, fMod, axis=0, fill_value='extrapolate')
cModLinInt = interp1d(muv, cMod, axis=0, fill_value='extrapolate')

''' FIXED SIMULATION PARAMETERS - These come from the precision tests in Chiavassa et al. (2026) '''

nS = 100          # Stellar radius in pixels
nP = nS * pars.rr # Planet radius in pixels
fb = 1 / nS       # Full bin size in stellar radius units
hb = fb / 2.0     # Half bin size in stellar radius units

''' Creating the super-sampled star '''
n = 5000
xcen = 2499.5
ycen = 2499.5
radius = 2500.0
x_2 = (np.arange(n)-xcen)**2   # X distances squared
y_2 = (np.arange(n)-ycen)**2   # Y distances squared
distMatrix = np.zeros((n,n))
for i in range(n): distMatrix[i,] = np.sqrt(x_2 + y_2[i])
starFull = np.zeros((n,n))
starFull[distMatrix <= 2500.0] = 1.0

''' Building auxiliary vectors +ve pixel grid and full grid '''
halfPixGrid = np.arange(nS) / float(nS) + hb  # (vec) Positive half, centred on half pixel
fullPixGrid = np.append(-1*halfPixGrid[::-1],halfPixGrid)     # (vec2) Appending the negative half

''' Rebinning to desired output size and computing the velocity field '''
nOut = nS*2
assert (n / nOut) % 1 == 0.0 # Supersampled must be an exact multiple of the final size
scaleFact = int(n / nOut)
star = np.zeros((nOut,nOut))
vLOS = np.zeros((nOut,nOut))
muStar = np.zeros((nOut,nOut))
for j in range(nOut):
    yCoord = fullPixGrid[j]
    for i in range(nOut):
        xCoord = fullPixGrid[i]
        sqrtArg = 1.0 - distMatrix[j*scaleFact:(j+1)*scaleFact,i*scaleFact:(i+1)*scaleFact].flatten() / radius
        weightValue = np.mean(starFull[j*scaleFact:(j+1)*scaleFact,i*scaleFact:(i+1)*scaleFact])
        star[j,i] = weightValue
        if weightValue == 1:  # Fully illuminated => assign the LOS velocity of the pixel centre
            muStar[j,i] = np.mean(np.sqrt(sqrtArg))
            vLOS[j,i] = vEq * np.sin(pars.iStar_rad) * (np.sin(pars.lam_rad)*yCoord - np.cos(pars.lam_rad)*xCoord)
        elif weightValue > 0: # Partially illuminated => assign the LOS velocity of the edge of the slice ** CAN BE IMPROVED **
            muStar[j,i] = np.mean(np.sqrt(sqrtArg[sqrtArg >= 0.0]))
            xCoord = np.sqrt(1.0 - yCoord**2) * np.sign(fullPixGrid[i])
            vLOS[j,i] = vEq * np.sin(pars.iStar_rad) * (np.sin(pars.lam_rad)*yCoord - np.cos(pars.lam_rad)*xCoord)
        else:                 # Not illuminated => assign 'NaN'
            vLOS[j,i] = 'NaN'
rvStar = -1.0 * vLOS

''' Plotting the star and the velocity field '''

plt.figure(figsize=(10,10))
plt.imshow(rvStar,origin='lower',aspect=1.0,interpolation='None',extent=[-nS,nS,-nS,nS],vmin=-vEq,vmax=vEq, cmap='coolwarm')
plt.title('$veq={:2.1f}$ km/s, $\lambda={:3.0f}^\circ, i_\star={:3.0f}^\circ$'.format(vEq, pars.lam_deg, pars.iStar_deg))
plt.colorbar(label='RV (km/s)')
plt.show()

plt.figure(figsize=(10,10))
plt.imshow(star,origin='lower',aspect=1.0,interpolation='None',extent=[-nS,nS,-nS,nS])
plt.show()

''' Full-disk recomputation if the keyword is set (this is slow!) '''

if doFullDisk:
    fullSpec = np.zeros(nWlen)
    fullCont = np.zeros(nWlen)
    cond = star != 0.0
    nidx = cond.sum()  # Number of non-zero pixels
    idxStar = np.where(cond) # A list of indices for each dimension of substar
    ix = idxStar[1] # Indices of the horizontal dimension
    iy = idxStar[0] # Indices of the vertical dimension
    xCoord = fullPixGrid[ix]
    yCoord = fullPixGrid[iy]
    ''' All slices '''
    mu = muStar[cond]
    pixelWeight = star[cond]
    rvtot = rvStar[cond]
    print('- Computing interpolation for full-disk stellar spectrum:')
    for i in range(nidx):
        pipe.progress(i+1,nidx,str(nidx))
        fracPix = pixelWeight[i]
        spec = fModLinInt(mu[i]) * fracPix
        cont = cModLinInt(mu[i]) * fracPix
        wModShifted = wMod * (1.0 + rvtot[i]/2.998E5)  # Doppler shift
        ''' Uncomment / comment below to enable / disable SPLINE interpolation '''
        #fullSpec += splev(wMod,splrep(wModShifted,spec,s=0),der=0)
        #fullCont += splev(wMod,splrep(wModShifted,cont,s=0),der=0)
        ''' Uncomment / comment below to enable / disable LINEAR interpolation '''
        fullSpec += interp1d(wModShifted,spec,bounds_error=None, fill_value='extrapolate')(wMod)
        fullCont += interp1d(wModShifted,cont,bounds_error=None, fill_value='extrapolate')(wMod)
    plt.figure(figsize=(11,3))
    plt.plot(wMod, fullSpec / fullCont)
    plt.xlabel('Wavelength')
    plt.ylabel('Normalized stellar flux')
    plt.show()
    np.save(outFold+pars.fullSpecFile,fullSpec)
    np.save(outFold+pars.fullContFile,fullCont)
else:
    print('- Skipping computation of the full-disk spectrum, loading from '+pars.fullSpecFile)
    fullSpec = np.load(outFold+pars.fullSpecFile)
    fullCont = np.load(outFold+pars.fullContFile)

''' Computation of the stellar regions occulted by the exoplanet transit.
It only uses the strip effectively crossed by the exoplanet, therefore in this case
the code needs to correctly rotate the velocity field based on the obliquity. '''

thetaLim = np.arcsin( np.sqrt(((1+pars.rr)**2 - pars.b**2)/(pars.a**2-pars.b**2)) )  # 4th contact
bLow = pars.b * np.cos(thetaLim)
yMax = pars.b + pars.rr + hb
yMin = bLow - pars.rr - hb

iStrip = (fullPixGrid >= yMin) * (fullPixGrid <= yMax)  # Corresponding elements of +ve half coordinate vector
nStripY = iStrip.sum()         # Height in pixels of the planet slice crossing the star
yStrip = fullPixGrid[iStrip].copy()  # (vec3) Sliced coordinate vector

print('- Planet strip ranging from y = {:1.3f} to {:1.3f} ({:4.0f} pixels)'.format(yMin,yMax,nStripY))

''' Precomputing sub-grid for planet-sub-pixel approximation '''
nrand = 10000
xrand = np.random.uniform(size=nrand)*fb - hb  # fb = Full bin size in stellar radii
yrand = np.random.uniform(size=nrand)*fb - hb

''' Initialise the final stellar spectrum '''
phases = np.arange(-0.01,0.0101,0.002)   # Orbital phases to simulate
nPhases = len(phases)
fstar = np.zeros((nPhases,nWlen))

if doPlShadow:
    plnt = np.zeros((nPhases, nStripY, nS*2))  # Empty strip cube
    ''' Computing the (x, y) coordinates of planet centre '''
    xP = pars.a * np.sin(2.0 * np.pi * phases)
    yP = pars.b * np.cos(2.0 * np.pi * phases)
    dd = np.sqrt(xP**2 + yP**2)  # Distance of planet centre from stellar centre
    distPix2PlCenLim = pars.rr + hb * np.sqrt(2)
    print('- Computing planet disk shadow:')
    for j in range(nPhases):
        print(phases[j], xP[j], yP[j], dd[j])
        plnt[j,] = star[iStrip,].copy()
        if dd[j] < (1+pars.rr): # computing occlusion only if distance is less than 1+rr
            for iy in range(nStripY):
                yCoord = yStrip[iy] - yP[j]
                distPix2PlCen = np.sqrt(yCoord**2 + (xP[j]-fullPixGrid)**2)
                occludedPixels = distPix2PlCen < distPix2PlCenLim
                nOccluded = occludedPixels.sum()
                plntStrip = plnt[j,iy,occludedPixels]
                # print(distPix2PlCen[occludedPixels]<)
                for ix in range(nOccluded):
                    xCoord = fullPixGrid[occludedPixels][ix] - xP[j]
                    # Distances of pixel vertices from planet centre
                    dtr = np.sqrt((xCoord+hb)**2.+(yCoord+hb)**2.)
                    dtl = np.sqrt((xCoord-hb)**2.+(yCoord+hb)**2.)
                    dbl = np.sqrt((xCoord-hb)**2.+(yCoord-hb)**2.)
                    dbr = np.sqrt((xCoord+hb)**2.+(yCoord-hb)**2.)
                    # Counting vertices distant less than rr from pl cnt
                    dvec = np.array([dtr,dtl,dbl,dbr])
                    dtot = (dvec <= pars.rr).sum()
                    ''' If all vertices are inside 1 planet radius assign 0, 
                    otherwise subgrid and assign the fraction of pixel not occulted. '''
                    if dtot == 4: 
                        plntStrip[ix] = 0.0
                    elif dtot > 0:
                        term1 = (xrand+xCoord)**2.
                        term2 = (yrand+yCoord)**2.
                        subGridDist = term1 + term2
                        condition =  subGridDist > pars.rr**2.
                        frac = np.sum(condition)
                        plntStrip[ix] *= float(frac) / float(nrand)
                if iy == nStripY-1:
                    print(len(plntStrip)-plntStrip.sum())
                plnt[j,iy,occludedPixels] = plntStrip.copy()
        plt.figure(figsize=(0.1*nS*2,0.1*nStripY))
        plt.imshow(plnt[j,],aspect=1.0,interpolation='None',origin='lower')
        plt.yticks([])
        plt.xticks([])
        plt.show()
    np.save(outFold+pars.plShadowFile,plnt)
else:
    plnt = np.load(outFold+pars.plShadowFile)
        
''' For each timestep subtract the occulted spectrum '''
print('Subtracting occulted disk from each frame...')
starStrip = star[iStrip,].copy()
muStrip = muStar[iStrip,].copy()
rvStrip = rvStar[iStrip,].copy()
for j in range(nPhases):
    plnt0 = plnt[j,].copy()
    diff = starStrip - plnt0
    plt.imshow(diff,aspect=1.0,interpolation='None',origin='lower')
    plt.show()
    ipl = diff != 0
    ni = ipl.sum()
    plsp = np.zeros(nWlen)
    print('- Frame {0}/{1} has {2} occulted pixels'.format(j+1,nPhases,diff[ipl].sum()))
    if ni > 0:
        mu = muStrip[ipl]
        xmu = splev(mu, splrep(np.arange(len(muv)),muv,s=0),der=0)
        subDiff = diff[ipl]
        rvtot = rvStrip[ipl]
        for i in range(ni):
            wMod1 = wMod * (1.0 + rvtot[i]/2.998E5)
            spc = fModLinInt(xmu[i]) * subDiff[i]
            plsp += interp1d(wMod1,spc,bounds_error=None, fill_value='extrapolate')(wMod)       
    fstar[j,] = (fullSpec - plsp) / fullCont
plt.imshow(fstar,aspect='auto')
plt.show()
np.save(outFold+pars.finalModelFile,fstar)