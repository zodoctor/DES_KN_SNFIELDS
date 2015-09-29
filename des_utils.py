import numpy as np
import des_io

def get_all_obs(infiles):
# parses data from each file into a structured array and adds it to a list.
# Output: a list of structured arrays and a list of header dictionaries
# corresponding to each object  
    all_obs_list = []
    all_headerdicts = []
    for infile in infiles:
        (obs,headerdict) = des_io.parse_observations(infile)
        all_obs_list.append(obs)
        all_headerdicts.append(headerdict)
    return all_obs_list, all_headerdicts

def hostgalAssociate(infile=None,headerlist=None,datatype='sim'):
    if (infile==None) and (headerlist==None):
        print 'you need either a file for random hostgal association or a headerlist'
    elif (headerlist!=None) and (infile!=None): 
        readfile = np.genfromtxt(infile,names=True,dtype=float)
        randomindices = np.random.randint(0,len(readfile),size=len(headerlist))
        DLR = readfile['DLR'][randomindices]
        SEP = readfile['SEP'][randomindices]
        ZGAL_MATCH = readfile['ZGAL_MATCH'][randomindices]
        return ZGAL_MATCH, DLR, SEP
    elif (infile==None) and (datatype=='sim'):
        HOSTGAL_PHOTOZ = np.array([float(eye['HOSTGAL_PHOTOZ']) for i,eye in enumerate(headerlist)])
        SEP = np.array([float(eye['HOSTGAL_SNSEP']) for i,eye in enumerate(headerlist)])
        return HOSTGAL_PHOTOZ,SEP
    elif (infile==None) and (datatype=='data'):
        HOSTGAL_PHOTOZ = np.array([float(eye['HOSTGAL_PHOTOZ']) for i,eye in enumerate(headerlist)])
        DLR = np.array([float(eye['PRIVATE(DES_hostgal_dlr)']) for i,eye in enumerate(headerlist)])
        SEP = np.array([float(eye['HOSTGAL_SNSEP']) for i,eye in enumerate(headerlist)])
        return HOSTGAL_PHOTOZ,DLR,SEP 
    else:
        print 'something went wrong...try again?'

def cut_on_photoZ(all_headerdicts,zmax=0.1):
# outputs a boolean selector list with true elements corresponding to objects
# with photoZ greater than zmax
    photoZcutsel = np.zeros(len(all_headerdicts),dtype='bool')
    for i,headerdict in enumerate(all_headerdicts):
        photoZcutsel[i] = photoZcut(headerdict,zmax)
    return photoZcutsel

def get_depth_lists(all_obs_list):
# takes in all_obs_list and sorts it based on deep and shallow fields
# outputs shallow and deep lists with same length as all_obs_list
    deeplist = np.empty(len(all_obs_list),dtype=object)
    shallowlist = np.empty(len(all_obs_list),dtype=object)
    for i,obs in enumerate(all_obs_list):
        deep_sel = deepfield(obs)
        deeplist[i] = obs[deep_sel]
        shallowlist[i] = obs[~deep_sel]
    return shallowlist, deeplist
 
def get_band_info(all_obs_list,band,zp_lower=30.5,zp_upper=34.0, zp_fwhm_lower=-.5, zp_fwhm_upper=2.0,photprobmin=0.5):
# outputs a list with length of all_obs_list which contains data about a given band
# for all objects in all_obs_list
    bandinfolist = np.empty(len(all_obs_list),dtype=object)
    for i,obs in enumerate(all_obs_list):
        bandinfolist[i] = obsinband(obs,band,zp_lower,zp_upper, zp_fwhm_lower, zp_fwhm_upper,photprobmin)
    return bandinfolist

def trigger_selector(bandinfolist):
# outputs a list of binary selector lists which have a true element for a trigger nite 
    bandtrig = np.zeros(len(bandinfolist),dtype=object)
    for i,info in enumerate(bandinfolist):
        bandtrig[i] = info[0] == 2
    return bandtrig

def existobs_selector(bandinfolist):
# don't think this function is needed....
    bandobs = np.zeros(len(bandinfolist),dtype=object)
    for i,info in enumerate(bandinfolist):
        bandobs[i] = info[0] > 0
    return bandobs

def common_trignite_selector(bandinfolist1,bandinfolist2,SNRsel1,SNRsel2,trigreq = 1,deep=0,sim=1,SNRand=0):
# returns a list of binary selector lists which give common trigger nites
    sel_list1 = np.empty(len(bandinfolist1),dtype=object)
    sel_list2 = np.empty(len(bandinfolist1),dtype=object)
    cMJDlist = np.empty(len(bandinfolist1),dtype=object)
    if trigreq==1 and sim==1:
        for i in range(0,len(bandinfolist1)):
            sel1, sel2,cnites = common_nites(bandinfolist1[i][1],bandinfolist2[i][1],bandinfolist1[i][0],bandinfolist2[i][0],SNRsel1[i],SNRsel2[i],deep,SNRand=SNRand)
            sel_list1[i] = sel1
            sel_list2[i] = sel2
            cMJDlist[i] = cnites
    elif trigreq==1 and sim==0:
        for i in range(0,len(bandinfolist1)):
            sel1, sel2,cnites = common_nites(bandinfolist1[i][1],bandinfolist2[i][1],bandinfolist1[i][0],bandinfolist2[i][0],None,None,deep,SNRand=SNRand)
            sel_list1[i] = sel1
            sel_list2[i] = sel2
            cMJDlist[i] = cnites
    else:
        for i in range(0,len(bandinfolist1)):
            sel1, sel2,cnites = common_nites(bandinfolist1[i][1],bandinfolist2[i][1])
            sel_list1[i] = sel1
            sel_list2[i] = sel2
            cMJDlist[i] = cnites
    return sel_list1, sel_list2,cMJDlist

def get_SNR_selector(bandinfolist,SNRmin=5.0):
# outputs a list of binary selector lists that have a true element for each nite which has
# an SNR that is > SNRmin in that band
    SNR_selector = np.zeros(len(bandinfolist),dtype=object)
    for i,info in enumerate(bandinfolist):
        if info[3].any():
            SNR_selector[i] = info[3] >= SNRmin
    return SNR_selector

def get_trig_flags_list(trig_sel1,trig_sel2,common_trig_sel1,common_trig_sel2,SNR_sel1=None,SNR_sel2=None,SNRand=0):
# outputs a trig_flags_list which is a list of binary selector lists with a true element for nites that have
# a trigger and an SNR greater than SNRmin in at least one of the trigger observations.  If SNRand flag is 
# set to 1, require both trigger observations have SNR greater than SNRmin 
    trig_flags_list = np.empty(len(trig_sel1),dtype=object)
    anytrigs = np.zeros(len(trig_sel1),dtype='bool')
    for i in range(0,len(trig_sel1)):
        if SNR_sel1 is None and SNR_sel2 is None:
            trig_flags_list[i] = trig_sel1[i][common_trig_sel1[i]] & trig_sel2[i][common_trig_sel2[i]]
        elif np.any(SNR_sel1[i]) and SNRand:
            trig_flags_list[i] = trig_sel1[i][common_trig_sel1[i]] & trig_sel2[i][common_trig_sel2[i]] & (SNR_sel1[i][common_trig_sel1[i]] & SNR_sel2[i][common_trig_sel2[i]])
        elif np.any(SNR_sel1[i]):
            trig_flags_list[i] = trig_sel1[i][common_trig_sel1[i]] & trig_sel2[i][common_trig_sel2[i]] & (SNR_sel1[i][common_trig_sel1[i]] | SNR_sel2[i][common_trig_sel2[i]])
        else:
            continue
            #trig_flags_list[i] = trig_sel1[i][common_trig_sel1[i]] & trig_sel2[i][common_trig_sel2[i]]
        anytrigs[i] = trig_flags_list[i].any()
    return trig_flags_list,anytrigs


def get_trig_MJD_list(band_info_list,trig_flags_list,trig_sel_list):
# gives a list of lists of MJDs which had triggers
    MJDtriglist = np.empty(len(trig_flags_list),dtype=object)
    for i,trig in enumerate(trig_flags_list):
        MJDtriglist[i] = band_info_list[i][1][trig_sel_list[i]][trig]
    return MJDtriglist


def get_detection_flags_list(MJDtriglist,bandinfolist1,bandinfolist2,nitesepmin=7,nitesepmax=7,iandzfollowup=1):
# gives a binary selector list for each object that had a detection
    detection_flags_list = np.zeros(len(MJDtriglist),dtype='bool')
    for i in range(0,len(MJDtriglist)):
        detection_flags_list[i] = followupdet(MJDtriglist[i],bandinfolist1[i][1],bandinfolist2[i][1],nitesepmin,nitesepmax,iandzfollowup)
    return detection_flags_list  

def get_detection_info(MJDtriglist,bandinfolist,nitesepmin=7,nitesepmax=7,iandzfollowup=1):
    detection_fluxes = np.empty(len(MJDtriglist),dtype='float')
    detection_fluxes[:] = np.NAN
    detection_MJD = np.empty(len(MJDtriglist),dtype='int')
    detection_MJD[:] = np.NAN
    detection_flags_list = np.zeros(len(MJDtriglist),dtype='bool')
    for i in range(0,len(MJDtriglist)):
        if MJDtriglist[i].size != 0:
            detection_nites_min_sel = nitesepmin <= (bandinfolist[i][1]-MJDtriglist[i][-1])
            detection_nites_max_sel = (bandinfolist[i][1]-MJDtriglist[i][-1]) <= nitesepmax
            detection_nites_sel = detection_nites_min_sel & detection_nites_max_sel
            if np.any(detection_nites_sel):
                detection_flags_list[i] = np.any(detection_nites_sel)
                detection_MJD[i] = bandinfolist[i][1][detection_nites_sel][0]
                detection_fluxes[i] = bandinfolist[i][2][detection_nites_sel][0]
    return detection_fluxes, detection_MJD,detection_flags_list

def extract_colors(infiles):
    nobjects = len(infiles)

    triggers = np.zeros(nobjects, dtype='bool')
    colors = np.zeros(nobjects)
    ifluxes = np.zeros(nobjects)
    detections = np.zeros(nobjects,dtype='bool')
    deltaT = []
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
        zsel, isel,cnites = common_nites(zMJD, iMJD, zobs, iobs)

        trig_flags = ztrig[zsel] & itrig[isel]
        # if SNR is defined (only for sims) make sure at least one trigger obs has adequate SNR
        if zSNR.any():
            zSNRpass = zSNR >= 5
            iSNRpass = iSNR >= 5
            trig_flags = trig_flags & (zSNRpass[zsel] | iSNRpass[isel])

        trig = np.any(trig_flags)
        MJDtrig = zMJD[zsel][trig_flags]
        if len(MJDtrig)>1: 
            deltaT.append(MJDtrig.max() - MJDtrig.min())
        # cut out object if it has multiple triggers within maxtrignites
        if ((not allowmultitrig) and multitrig(MJDtrig,maximumtrignites)):
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
    return triggers, colors, ifluxes, np.sum(detections), SNIDset,np.array(deltaT)

def followupdet(MJDtrig,zMJD,iMJD,nitesepmin=7,nitesepmax=7,iandzfollowup = 1):
    detected = False
    if MJDtrig.size > 0:
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
    if len(MJDtrig) == 0:
        multitrigflag = 1
    else:   
        trigdiff = MJDtrig.max() - MJDtrig.min()
        multitrigflag = trigdiff>maxtrignites
    return multitrigflag

def deepfield(obs):
    sel = (obs['FIELD'] == 'X3') | (obs['FIELD'] == 'C3')
    return sel

def obsinband(obslist, band='i',zp_lower=30.5, zp_upper=34.0, zp_fwhm_lower=-.5, zp_fwhm_upper=2.0,photprobmin=0.5):
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
    obsZP = np.empty(nnites)
    obsField = np.empty(nnites,dtype=object)
    obsfluxcalerr = np.empty(nnites)
    obsPSF = np.empty(nnites,dtype=float)
    obsPHOTFLAG = np.empty(nnites,dtype=bool)
    for x, nite in enumerate(nitelist):
        nite_obs = obs[obs['MJD'].astype('int') == nite]
        det = detected(nite_obs,photprobmin)
        passed = within_cuts(nite_obs,zp_lower, zp_upper, zp_fwhm_lower, zp_fwhm_upper)
        sel = det & passed
### When there are multiple triggers on the same night in the same band, we choose the one with the better SNR.  If the SNRs are the same, we pick the first one in the list to get the flux.
        if np.any(sel):
            try:
                SNRmax = np.max(nite_obs[sel]['SNR'])
                SNRsel = nite_obs[sel]['SNR'] >= SNRmax
            except ValueError:
                obsSNRlist = np.absolute(np.divide(nite_obs[sel]['FLUXCAL'],nite_obs[sel]['FLUXCALERR']))
                SNRmax = np.max(obsSNRlist)
                SNRsel = obsSNRlist >= SNRmax
            obsband[x] = 2
            obsflux[x] = nite_obs[sel]['FLUXCAL'][SNRsel][0]
            obsfluxcalerr[x] = nite_obs[sel]['FLUXCALERR'][SNRsel][0]
            obsSNR[x] = SNRmax
            obsField[x] = nite_obs[sel]['FIELD'][SNRsel][0]
            try:
                obsPSF[x] = nite_obs[sel]['PSF'][SNRsel][0]
            except ValueError:
                obsPSF[x] = 0.0
            zpname = 'ZPFLUX'
            try:
                nite_obs[sel][zpname]
            except ValueError:
                zpname = 'ZPT'
            obsZP[x] = nite_obs[sel][zpname][SNRsel][0]
            obsPHOTFLAG[x] = nite_obs[sel]['PHOTFLAG'][SNRsel][0] > 1000
        else:
            obsband[x] = 1
            try:
                SNRmax = np.max(nite_obs['SNR'])
                SNRsel = nite_obs['SNR'] >= SNRmax
            except:
                obsSNRlist = np.absolute(np.divide(nite_obs['FLUXCAL'],nite_obs['FLUXCALERR']))
                SNRmax = np.max(obsSNRlist)
                SNRsel = obsSNRlist >= SNRmax
            obsflux[x] = nite_obs['FLUXCAL'][SNRsel][0]
            obsfluxcalerr[x] = nite_obs['FLUXCALERR'][SNRsel][0]
            obsSNR[x] = SNRmax
            obsField[x] = nite_obs['FIELD'][SNRsel][0]
            try:
                obsPSF[x] = nite_obs['PSF'][SNRsel][0]
            except ValueError:
                obsPSF[x] = 0.0
            zpname = 'ZPFLUX'
            try:
                nite_obs[sel][zpname]
            except ValueError:
                zpname = 'ZPT'
            obsZP[x] = nite_obs[zpname][SNRsel][0]
            obsPHOTFLAG[x] = nite_obs['PHOTFLAG'][SNRsel][0] > 1000
    return (obsband, nitelist, obsflux, obsSNR,obsfluxcalerr,obsField,obsPSF,obsPHOTFLAG,obsZP)

def exist_deep_trigs(zobs, iobs, zMJD,iMJD):
    zcnites,icnites = common_nites(zMJD,iMJD)
    ztrig = zobs == 2
    itrig = iobs == 2
    trig = any((zMJDtrig - iMJDtrig <=1) for zMJDtrig in zMJD[ztrig] for iMJDtrig in iMJD[itrig])
    return trig

def common_nites(nitelist1, nitelist2, obs1=None, obs2=None,SNRsel1=None,SNRsel2=None,deep=0,SNRand=0):
    if deep == 0:
        if SNRsel1 == None:
            goodnites1 = nitelist1 if obs1 is None else nitelist1[obs1 == 2]
            goodnites2 = nitelist2 if obs2 is None else nitelist2[obs2 == 2]
            cnites = np.intersect1d(goodnites1, goodnites2)
        elif SNRand==0:
            goodnites1 = nitelist1 if obs1 is None else nitelist1[obs1 == 2]
            goodnites1SNR = nitelist1 if obs1 is None else (nitelist1[(obs1 == 2) & SNRsel1])
            goodnites2 = nitelist2 if obs2 is None else nitelist2[obs2 == 2]
            goodnites2SNR = nitelist2 if obs2 is None else (nitelist2[(obs2 == 2) & SNRsel2])
            cnitesSNR1 = np.intersect1d(goodnites1SNR,goodnites2)
            cnitesSNR2 = np.intersect1d(goodnites1,goodnites2SNR)
            cnites = np.union1d(cnitesSNR1,cnitesSNR2)
        elif SNRand==1:
            goodnites1 = nitelist1 if obs1 is None else (nitelist1[(obs1 == 2) & SNRsel1])
            goodnites2 = nitelist2 if obs2 is None else (nitelist2[(obs2 == 2) & SNRsel2])
            cnites = np.intersect1d(goodnites1,goodnites2)
        sel1 = np.in1d(nitelist1, cnites)
        sel2 = np.in1d(nitelist2, cnites)
    else:
        sel1, sel2,cnites = common_nites_deep(nitelist1,nitelist2,obs1,obs2,SNRsel1,SNRsel2,SNRand)
    return sel1, sel2,cnites

def common_nites_deep(nitelist1,nitelist2,obs1=None,obs2=None,SNRsel1=None,SNRsel2=None,SNRand=0):
    sel1 = np.zeros(len(nitelist1),dtype='bool')
    sel2 = np.zeros(len(nitelist2),dtype='bool')
    cnites = np.array([])
    if SNRsel1 == None:
        for lag in [-1,0,1]:
            goodnites1 = nitelist1 if obs1 is None else (nitelist1[obs1 == 2]+lag)
            goodnites2 = nitelist2 if obs2 is None else nitelist2[obs2 == 2]
            cnites_intersect = np.intersect1d(goodnites1,goodnites2)
            cnites = np.union1d(cnites,cnites_intersect)
            sel1 = sel1 | np.in1d(nitelist1,cnites-lag)
            sel2 = sel2 | np.in1d(nitelist2,cnites)
    elif SNRand == 0: 
        for lag in [-1,0,1]:
            goodnites1 = nitelist1 if obs1 is None else (nitelist1[obs1 == 2]+lag)
            goodnites1SNR = nitelist1 if obs1 is None else (nitelist1[(obs1 == 2) & SNRsel1]+lag)
            goodnites2 = nitelist2 if obs2 is None else nitelist2[obs2 == 2]
            goodnites2SNR = nitelist2 if obs2 is None else (nitelist2[(obs2 == 2) & SNRsel2])
            cnitesSNR1 = np.intersect1d(goodnites1SNR,goodnites2)
            cnitesSNR2 = np.intersect1d(goodnites1,goodnites2SNR)
            cnites_intersect = np.union1d(cnitesSNR1,cnitesSNR2)
            #cnites = np.union1d(cnites,cnites_intersect)
            cnites = np.union1d(np.union1d(cnites,cnites_intersect),np.union1d(cnites,cnites_intersect-lag))
            sel1 = (sel1 | np.in1d(nitelist1,cnites-lag)) & (obs1 == 2)
            sel2 = (sel2 | np.in1d(nitelist2,cnites)) & (obs2 == 2)
    elif SNRand == 1:
        for lag in [-1,0,1]:
            goodnites1SNR = nitelist1 if obs1 is None else (nitelist1[(obs1 == 2) & SNRsel1]+lag)
            goodnites2SNR = nitelist2 if obs2 is None else (nitelist2[(obs2 == 2) & SNRsel2])
            cnites_intersect = np.intersect1d(goodnites1SNR,goodnites2SNR)
            #cnites = np.union1d(cnites,cnites_intersect)
            cnites = np.union1d(np.union1d(cnites,cnites_intersect),np.union1d(cnites,cnites_intersect-lag))
            sel1 = (sel1 | np.in1d(nitelist1,cnites-lag)) & (obs1 == 2)
            sel2 = (sel2 | np.in1d(nitelist2,cnites)) & (obs2 == 2)
    return sel1, sel2,cnites
        
def detected(obs,photprobmin=0.5):
    # given an observation dictionary, check if that observation counts as a detection
    # ---------------------------------------------------------
    try:
        sel = (max(obs['PHOTFLAG']) > 1000) and ((max(obs['PHOTFLAG']) & 1016) == 0)
    except ValueError:

        print obs['PHOTFLAG'], (obs['PHOTFLAG'] & 1016)
    try:
        sel = sel & (obs['PHOTPROB'] >= photprobmin)
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

    #zp_sel = (obs[zpname] > zp_lower) & (obs[zpname] < zp_upper)
    zp_sel = np.ones(len(obs), dtype='bool')
    try:
        psf_sel = (obs['PSF'] > zp_fwhm_lower) & (obs['PSF'] < zp_fwhm_upper)
    except ValueError:
        psf_sel = np.ones(len(obs), dtype='bool')

    return zp_sel #& psf_sel

def photoZcut(headerdict,zmax=.1):
    if (float(headerdict['HOSTGAL_PHOTOZ']) > zmax):
        photoZpass = True
    else:
        photoZpass = False
    return photoZpass
