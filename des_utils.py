import numpy as np
import des_io


def extract_colors(infiles):
    nobjects = len(infiles)

    triggers = np.empty(nobjects, dtype='bool')
    colors = np.zeros(nobjects)
    ifluxes = np.zeros(nobjects)
    for i, infile in enumerate(infiles):
        # get a list with all the values in the data table
        obs = des_io.parse_observations(infile)

        # Separate deep and shallow fields
        deep_sel = deepfield(obs)
        deep = obs[deep_sel]
        shallow = obs[~deep_sel]

        # Look just at shallow fields for now
        # make a list of all the nites there were observations and make a nitedictlist
        nitelist = np.unique(shallow['MJD'].astype('int'))

        zobs, zMJD, zflux, zSNR = obsinband(shallow, 'z') # identify whether there was a z observation on a nite and get its MJD
        iobs, iMJD, iflux, iSNR = obsinband(shallow, 'i') # identify whether there was an i observation on a nite and get its MJD

        zdet = zobs == 2
        idet = iobs == 2
        zsel, isel = common_nites(zMJD, iMJD)
        trig_flags = zdet[zsel] & idet[isel]
        trig = np.any(trig_flags)
        triggers[i] = trig
        if trig:
            zflux1 = zflux[zdet & zsel][0]
            iflux1 = iflux[idet & isel][0]
            colors[i] = -2.5*(np.log10(iflux1)-np.log10(zflux1))
            ifluxes[i] = iflux[iobs > 0][-1] - iflux[iobs > 0][0]
    return triggers, colors, ifluxes


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
