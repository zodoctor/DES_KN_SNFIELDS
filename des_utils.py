import numpy as np
import des_io


def extract_colors(infiles):
    nobjects = len(infiles)

    triggers = np.zeros(nobjects, dtype='bool')
    colors = np.zeros(nobjects)
    ifluxes = np.zeros(nobjects)
    detections = np.zeros(nobjects,dtype='bool')
    SNIDset = set()
    zcutflag = 0
    allowmultitrig = True
    for i, infile in enumerate(infiles):
        # get a list with all the values in the data table
        (obs,headerdict)= des_io.parse_observations(infile)

        # skip files associated with objects w/ redshift > zmax
        if zcutflag and photoZcut(headerdict):
            continue
        # Separate deep and shallow fields
        deep_sel = deepfield(obs)
        deep = obs[deep_sel]
        shallow = obs[~deep_sel]

        # Look just at shallow fields for now
        # make a list of all the nites there were observations and make a nitedictlist
        nitelist = np.unique(shallow['MJD'].astype('int'))

        zobs, zMJD, zflux, zSNR = obsinband(shallow, 'z') # identify whether there was a z observation on a nite and get its MJD
        iobs, iMJD, iflux, iSNR = obsinband(shallow, 'i') # identify whether there was an i observation on a nite and get its MJD

        ztrig = zobs == 2
        itrig = iobs == 2
        zdet = zobs > 0
        idet = iobs > 0
    
        # get boolean selector list of common trigger nites for z and i observations       
        zsel, isel = common_nites(zMJD, iMJD)

        # if SNR is defined (only for sims) make sure at least one trigger obs has adequate SNR
        if zSNR.any():
            zSNRpass = zSNR >= 5
            iSNRpass = iSNR >= 5
            trig_flags = ztrig[zsel] & itrig[isel] & (zSNRpass[zsel] | iSNRpass[isel])
        else:
            trig_flags = ztrig[zsel] & itrig[isel] 
        trig = np.any(trig_flags)
        MJDtrig = zMJD[trig_flags]

        # cut out object if it has multiple triggers within maxtrignites
        if ((not allowmultitrig) and multitrig(MJDtrig)):
            continue

        # record whether the trigger has a follow up observation for a full detection
        detections[i] = followupdet(MJDtrig,zMJD,iMJD)
        if detections[i]:
            SNIDset.add(headerdict['SNID'])
        triggers[i] = trig
        if trig:
            zflux1 = zflux[ztrig & zsel][0] # this doesn't quite work if there are multiple triggers -- need fix
            iflux1 = iflux[itrig & isel][0]
            colors[i] = -2.5*(np.log10(iflux1)-np.log10(zflux1))
            if np.isnan(colors[i]):
                print zflux1,iflux1,headerdict['SNID']
            ifluxes[i] = iflux[iobs > 0][-1] - iflux[iobs > 0][0]
    return triggers, colors, ifluxes, np.sum(detections), SNIDset

def followupdet(MJDtrig,zMJD,iMJD,nitesepmin=7,nitesepmax=7,iandzfollowup = 1):
    detected = False
    for tnite in MJDtrig:
        nitesepz = zMJD-tnite
        nitesepi = iMJD-tnite
        detectedz = any(nitesepmin <= nitesep <= nitesepmax for nitesep in nitesepz)
        detectedi = any(nitesepmin <= nitesep <= nitesepmax for nitesep in nitesepi)
        if iandzfollowup:
            detected = detectedz and detectedi
        else:
            detected = detectedz or detectedi
        if detected:
            break
    return detected  

def multitrig(MJDtrig,maxtrignites=10):
    multitrigflag = any((((0-maxtrignites) < MJD1 - MJD2 < maxtrignites) and MJD1 != MJD2) for MJD1 in MJDtrig for MJD2 in MJDtrig)
    return multitrigflag

def deepfield(obs):
    sel = (obs['FIELD'] == 'X3') | (obs['FIELD'] == 'C3')
    return sel

def obsinband(obslist, band='i'):
    # returns list with 0's for no observations in that band on that nite,
    # 1's for if there is a near field observation in that band on that nite with all cuts/detects met, 2 if 
    # there is a deep field obs on that nite with cuts/detects met.  3 and 4 signal that there was a 
    # near or deep (resp.) observation that night but with no requirement on those observations passing
    # cuts and detects.  Also returns obsMJD which is a list of the MJD's of those observations in that band
    # -----------------------------------
    obs = obslist[obslist['FLT'] == band]

    nitelist = np.unique(obs['MJD'].astype('int'))
    nnites = len(nitelist)

    obsband = np.zeros(nnites, dtype='int')
    obsflux = np.empty(nnites)
    obsSNR = np.empty(nnites)

    for x, nite in enumerate(nitelist):
        nite_obs = obs[obs['MJD'].astype('int') == nite]
        det = detected(nite_obs)
        passed = within_cuts(nite_obs)
        if np.any(det & passed):
            sel = det & passed
            obsband[x] = 2
            obsflux[x] = np.mean(nite_obs[sel]['FLUXCAL'])
	    try:
	        obsSNR[x] = np.mean(nite_obs[sel]['SNR'])
	    except ValueError:
		pass
        else:
            obsband[x] = 1
            obsflux[x] = np.mean(nite_obs['FLUXCAL'])
	    try:
                obsSNR[x] = np.mean(nite_obs['SNR'])
            except ValueError:
                pass
    return (obsband, nitelist, obsflux, obsSNR)

def exist_deep_trigs(zobs, iobs, zMJD,iMJD):
    zcnites,icnites = common_nites(zMJD,iMJD)    
    ztrig = zobs == 2
    itrig = iobs == 2
    trig = any((zMJDtrig - iMJDtrig <=1) for zMJDtrig in zMJD[ztrig] for iMJDtrig in iMJD[itrig])
    return trig

def common_nites(nitelist1, nitelist2):
    cnites = np.array([nite for nite in nitelist1 if np.any(nitelist2 == nite)])

    sel1 = np.array([True if nite in cnites else False for nite in nitelist1], dtype='bool')
    sel2 = np.array([True if nite in cnites else False for nite in nitelist2], dtype='bool')
    return sel1, sel2
    
def detected(obs):
    # given an observation dictionary, check if that observation counts as a detection
    # ----------------------------------------------------------
    sel = obs['PHOTFLAG'] > 1

    try:
        sel = sel & (obs['PHOTPROB'] >= 0.5)
    except ValueError:
        pass

    return sel 

def within_cuts(obs, zp_lower=30.5, zp_upper=34.0, zp_fwhm_lower=-.5, zp_fwhm_upper=2.0):
    # given an observation dictionary, check if that observation passes the cuts defined here
    # ------------------------------------------------------
    zpname = 'ZPFLUX'
    try:
        obs[zpname]
    except ValueError:
        zpname = 'ZPT'

    zp_sel = (obs[zpname] > zp_lower) & (obs[zpname] < zp_upper)

    try:
        psf_sel = (obs['PSF'] > zp_fwhm_lower) & (obs['PSF'] < zp_fwhm_upper)
    except ValueError:
        psf_sel = np.ones(len(obs), dtype='bool')

    return zp_sel & psf_sel

def photoZcut(headerdict,zmax=.1): 
    if float(headerdict['REDSHIFT_HELIO']) > zmax:
        photoZpass = True
    else:
        photoZpass = False        
    return photoZpass
