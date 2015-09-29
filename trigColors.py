import numpy as np
from glob import glob

import des_utils
import des_io

import matplotlib.pyplot as plt

def trigColors(path,datatype,fieldtype,useSNRflag=1,zmax=0.1,zp_lower=30.5,zp_upper=34.0, zp_fwhm_lower=-.5, zp_fwhm_upper=2.0,photprobmin=0.5,maxtrignite=20,SNRand=0):
    dict1 = dict()
    files = glob(path)
    thelist,theheaders = des_utils.get_all_obs(files)
    photoZcutsel = des_utils.cut_on_photoZ(theheaders,zmax)
    dict1['photoZcutsel'] = photoZcutsel
    shallow_list,deep_list = des_utils.get_depth_lists(thelist)
    if fieldtype == 'deep':
        depth_list = deep_list
    else:
        depth_list = shallow_list
    dict1['depthlist'] = depth_list
    zbandinfo = des_utils.get_band_info(depth_list,'z',zp_lower,zp_upper, zp_fwhm_lower, zp_fwhm_upper,photprobmin)
    ibandinfo = des_utils.get_band_info(depth_list,'i',zp_lower,zp_upper, zp_fwhm_lower, zp_fwhm_upper,photprobmin)
    gbandinfo = des_utils.get_band_info(depth_list,'g',zp_lower,zp_upper, zp_fwhm_lower, zp_fwhm_upper,photprobmin)
    rbandinfo = des_utils.get_band_info(depth_list,'r',zp_lower,zp_upper, zp_fwhm_lower, zp_fwhm_upper,photprobmin)
    dict1['zbandinfo']=zbandinfo
    dict1['ibandinfo']=ibandinfo
    dict1['gbandinfo']=gbandinfo
    dict1['rbandinfo']=rbandinfo
    if useSNRflag==1:
        zSNR_sel = des_utils.get_SNR_selector(zbandinfo)
        iSNR_sel = des_utils.get_SNR_selector(ibandinfo)
        if fieldtype == 'shallow':
            zsellist,isellist,cnites = des_utils.common_trignite_selector(zbandinfo,ibandinfo,zSNR_sel,iSNR_sel,1,0,1,SNRand=SNRand)
        elif fieldtype == 'deep':
            zsellist,isellist,cnites = des_utils.common_trignite_selector(zbandinfo,ibandinfo,zSNR_sel,iSNR_sel,1,1,1,SNRand=SNRand)
    else:
        if fieldtype == 'shallow':
            zsellist,isellist,cnites = des_utils.common_trignite_selector(zbandinfo,ibandinfo,None,None,1,0,0)
        elif fieldtype == 'deep':
            zsellist,isellist,cnites = des_utils.common_trignite_selector(zbandinfo,ibandinfo,None,None,1,1,0)
    dict1['zsel'] = zsellist
    dict1['isel'] = isellist
    dict1['cnites']=cnites
    colors = np.empty(len(zbandinfo))
    anytrigs = np.zeros(len(zbandinfo),dtype='bool')
    nitediff = np.zeros(len(zbandinfo))
    #multitrig = np.zeros(len(zbandinfo),dtype='bool')
    for i in range(0,len(zbandinfo)):
        nitediff[i] = -1
        if np.any(zsellist[i]):
            #imag = -2.5*np.log10(ibandinfo[i][2][isellist[i]][0])
            #zmag = -2.5*np.log10(zbandinfo[i][2][zsellist[i]][0])
            imag = -2.5*np.log10(ibandinfo[i][2][isellist[i]][0])
            zmag = -2.5*np.log10(zbandinfo[i][2][zsellist[i]][0])
            colors[i] = imag - zmag
            anytrigs[i] = np.any(zsellist[i]) and np.any(isellist[i])
            if cnites[i].size == 0:
                nitediff[i] = -1
            else:
                nitediff[i] = cnites[i].max() - cnites[i].min()
    multitrig = nitediff <= maxtrignite
    dict1['colors'] = colors
    dict1['trigs'] = anytrigs
    dict1['headers'] = theheaders
    dict1['multitrig'] = multitrig
    return dict1

def get_nitediff(zbandinfo,ibandinfo,rbandinfo,gbandinfo,trigseplim=20):
    nitediffsel = np.zeros(len(zbandinfo),dtype=bool)
    nitediffs = np.zeros(len(zbandinfo))
    for a,zband in enumerate(zbandinfo):
        ztrigsel = zband[0] == 2
        itrigsel = ibandinfo[a][0] == 2
        rtrigsel = rbandinfo[a][0] == 2
        gtrigsel = gbandinfo[a][0] == 2
        ztrignites = zband[1][ztrigsel]
        itrignites = ibandinfo[a][1][itrigsel]
        rtrignites = rbandinfo[a][1][rtrigsel]
        gtrignites = gbandinfo[a][1][gtrigsel]
        allnites = np.concatenate((ztrignites,itrignites,rtrignites,gtrignites))
        try:
            nitediff = np.amax(allnites)-np.amin(allnites)
        except ValueError:
            nitediff = 1000
        nitediffs[a] = nitediff
        nitediffsel[a] = nitediff < trigseplim
    return nitediffs,nitediffsel

def get_slopes(bandinfo,trignites,slopelim1,slopelim2):
    slopes1 = np.empty(len(bandinfo))
    slopes2 = np.empty(len(bandinfo))
    for r,zband in enumerate(bandinfo):
        try:
            ztrigsel = np.in1d(trignites[r],zband[1])
            ztrignite = trignites[r][ztrigsel][0]
            ztrigflux = zband[2][zband[1]==ztrignite]
            zfollownite = zband[1][zband[1]>ztrignite][0]
            zfollowflux = zband[2][zband[1]>ztrignite][0]
            zfollownite2 = zband[1][zband[1]>ztrignite][1]
            zfollowflux2 = zband[2][zband[1]>ztrignite][1]
        except IndexError:
            ztrignite = np.NAN
            ztrigflux = np.NAN
            zfollownite = np.NAN
            zfollowflux = np.NAN
            zfollownite2 = np.NAN
            zfollowflux2 = np.NAN
        try:
            slopes1[r] = (1-zfollowflux/ztrigflux)/(zfollownite-ztrignite)
        except ValueError:
            print zfollowflux,ztrigflux,zfollownite,ztrignite
        try:
            slopes2[r] = (1-zfollowflux2/ztrigflux)/(zfollownite2-ztrignite)
        except ValueError:
            print zfollowflux2,ztrigflux,zfollownite2,ztrignite
    slopesel1 = (slopes1 > slopelim1) | np.isnan(slopes1)
    slopesel2 = (slopes2 > slopelim2) | np.isnan(slopes2)
    return slopes1,slopes2,slopesel1,slopesel2

def get_zpcutsel(zbandinfo,ibandinfo,trignites,maxzp):
    zp = np.zeros(len(zbandinfo),dtype=float)
    zpsel = np.zeros(len(zbandinfo),dtype=bool)
    for i,info in enumerate(ibandinfo):
        isel = np.in1d(info[1],trignites[i])
        zsel = np.in1d(zbandinfo[i][1],trignites[i])
        if trignites[i].size != 0:
            izp = info[8][isel][0]
            zzp = zbandinfo[i][8][zsel][0]
            zp[i] = min(izp,zzp)
            zpsel[i] = (izp < maxzp) & (zzp < maxzp)
        else:
            zp[i] = 1000.
            zpsel[i] = False
    return zp,zpsel
        
def get_firsttrigSNRsel(ibandinfo,zbandinfo,trignites,field='shallow'):
    firsttrigSNRsel = np.zeros(len(ibandinfo),dtype=bool)
    for i,info in enumerate(ibandinfo):
        inites = info[1][info[0]==2]
        znites = zbandinfo[i][1][zbandinfo[i][0]==2]        
        if field == 'shallow':
            nites = np.intersect1d(inites,znites)
        else:
            trigs = np.intersect1d(inites,znites)
            ptrigs = np.intersect1d(inites,znites+1)
            ntrigs = np.intersect1d(inites,znites-1)
            ptrigs = np.union1d(ptrigs,ptrigs-1)
            ntrigs = np.union1d(ntrigs,ntrigs+1)
            pntrigs = np.union1d(ptrigs,ntrigs)
            nites = np.union1d(trigs,pntrigs)
        if len(nites) == 0:
            firsttrigSNRsel[i] = False
        elif len(trignites[i]) == 0:
            firsttrigSNRsel[i] = False
        elif not nites[0] in trignites[i]:
            firsttrigSNRsel[i] = False
        else: 
            firsttrigSNRsel[i] = True
    return firsttrigSNRsel

def get_gzcutsel(zbandinfo,gbandinfo,trignites,maxgzratio):
    gzcutsel = np.zeros(len(gbandinfo),dtype=bool)
    gzratio = np.zeros(len(gbandinfo),dtype=float)
    for i,info in enumerate(zbandinfo):
        gtrigsel = gbandinfo[i][0] == 2 
        ztrigsel = info[0] == 2
        ztrignites = info[1][ztrigsel]
        gtrignites = gbandinfo[i][1]
        try:
            gkntrigsel = np.in1d(trignites[i],gtrignites)
        except ValueError:
            print gtrignites,trignites
        zkntrigsel = np.in1d(trignites[i],ztrignites)
        try:
            kntrignite = trignites[i][gkntrigsel & zkntrigsel][0]
        except IndexError:
            kntrignite = trignites[i][gkntrigsel & zkntrigsel]
        newzsel = ztrignites == kntrignite
        newgsel = gtrignites == kntrignite
        zflux = info[2][ztrigsel][newzsel]
        gflux = gbandinfo[i][2][newgsel]
        try:
            gzratio[i] = gflux/zflux
            gzcutsel[i] = gzratio[i] < maxgzratio
        except ValueError:
            gzratio[i] = np.NAN
            gzcutsel[i] = True
    return gzratio,gzcutsel
        
def get_noRGtrigsel(rbandinfo,gbandinfo,trignites):
    noRGtrigsel = np.zeros(len(rbandinfo),dtype=bool)
    for i,info in enumerate(rbandinfo):
        rtrigsel = info[0] == 2
        gtrigsel = gbandinfo[i][0] == 2
        rtrignites = info[1][rtrigsel]
        gtrignites = gbandinfo[i][1][gtrigsel]
        rcnitesel = np.in1d(trignites[i],rtrignites)
        gcnitesel = np.in1d(trignites[i],gtrignites)
        noRGtrigsel[i] = ~np.any(gcnitesel)
        #noRGtrigsel[i] = ~(np.any(rcnitesel) or np.any(gcnitesel))
    return noRGtrigsel

def get_angsep(headers,angsepmax):
    angsep = np.zeros(len(headers))
    for h,head in enumerate(headers):
        try: 
            angsep[h] = float(head['PRIVATE(DES_angsep_trigger)'])
        except KeyError:
            angsep[h] = 0
    angsepsel = angsep < angsepmax
    return angsep, angsepsel

def get_laterobssel(zbandinfo,ibandinfo,rbandinfo,gbandinfo,trignites,minnites,maxnites):
# get a selector list which has a true element when its corresponding trigger has an
# observation in any band between minnites and maxnites nights after the trigger.
    bands = ['z','i','r','g']
    allbandinfos = [zbandinfo,ibandinfo,rbandinfo,gbandinfo]
    laterobssel = np.zeros(len(zbandinfo),dtype=bool)
    for b,bandinfo in enumerate(allbandinfos):
        for i,info in enumerate(bandinfo):
            try:
                firsttrignite = trignites[i][0]
                obsaftersel = (info[1] >= (firsttrignite+minnites)) & (info[1] <= (firsttrignite+maxnites))
                laterobssel[i] = np.any(obsaftersel) or laterobssel[i]
            except IndexError:
                laterobssel[i] = False
    return laterobssel

def get_laterdetsel(zbandinfo,ibandinfo,rbandinfo,gbandinfo,trignites,minnites,maxnites):
# get a selector list which has a true element when its corresponding trigger has no pipeline 
# trigger in any band between minnites and maxnites nights after the trigger.
    bands = ['z','i','r','g']
    allbandinfos = [zbandinfo,ibandinfo,rbandinfo,gbandinfo]
    laterdetsel = np.ones(len(zbandinfo),dtype=bool)
    for b,bandinfo in enumerate(allbandinfos):
        for i,info in enumerate(bandinfo):
            try:
                firsttrignite = trignites[i][0]
                detsel = info[7]
                detaftersel = (info[1][detsel] >= (firsttrignite+minnites)) & (info[1][detsel] <= (firsttrignite+maxnites))
                laterdetsel[i] = ~np.any(detaftersel) and laterdetsel[i]
            except IndexError:
                laterdetsel[i] = False
    return laterdetsel

def get_priorobssel(zbandinfo,ibandinfo,rbandinfo,gbandinfo,trignites,minnites,maxnites):
# get a selector list which has a true element when its corresponding trigger has an
# observation in any band between minnites and maxnites nights before the trigger.
    bands = ['z','i','r','g']
    allbandinfos = [zbandinfo,ibandinfo,rbandinfo,gbandinfo]
    priorobssel = np.zeros(len(zbandinfo),dtype=bool)
    for b,bandinfo in enumerate(allbandinfos):
        for i,info in enumerate(bandinfo):
            try:
                firsttrignite = trignites[i][0]
                obsbeforesel = (info[1] <= (firsttrignite-minnites)) & (info[1] >= (firsttrignite-maxnites))
                priorobssel[i] = np.any(obsbeforesel) or priorobssel[i]
            except IndexError:
                priorobssel[i] = False
    return priorobssel

def get_priordetsel(zbandinfo,ibandinfo,rbandinfo,gbandinfo,trignites,minnites,maxnites):
# get a selector list which has a true element when its corresponding KN trigger has no 
# pipeline trigger in any band between minnites and maxnites nights before the KN trigger.
    bands = ['z','i','r','g']
    allbandinfos = [zbandinfo,ibandinfo,rbandinfo,gbandinfo]
    priordetsel = np.ones(len(zbandinfo),dtype=bool)
    for b,bandinfo in enumerate(allbandinfos):
        for i,info in enumerate(bandinfo):
            try:
                firsttrignite = trignites[i][0]
                detsel = info[7]
                detbeforesel = (info[1][detsel] <= (firsttrignite-minnites)) & (info[1][detsel] >= (firsttrignite-maxnites))
                priordetsel[i] = ~np.any(detbeforesel) & priordetsel[i]
            except IndexError:
                priordetsel[i] = False
    return priordetsel

def get_PSFsel(zbandinfo,ibandinfo,trignites,PSFmax):
    PSFsel = np.zeros(len(zbandinfo),dtype=bool)
    PSFlist = np.empty(len(zbandinfo),dtype=float)
    for i,info in enumerate(zbandinfo):
        znitesel = np.in1d(info[1],trignites[i])
        initesel = np.in1d(ibandinfo[i][1],trignites[i])
        try:
            zPSF = info[6][znitesel][0]
            iPSF = ibandinfo[i][6][initesel][0]
            PSFsel[i] = ((zPSF < PSFmax) and (iPSF < PSFmax))
            PSFlist[i] = max(zPSF,iPSF)
        except IndexError:
            PSFsel[i] = True
            PSFlist[i] = np.NAN 
    return PSFlist,PSFsel
