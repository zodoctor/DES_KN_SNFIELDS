import numpy as np
from glob import glob

import des_utils
import des_io

import matplotlib.pyplot as plt

def trigColors(path,datatype,fieldtype,useSNRflag=1,zmax=0.1,zp_lower=30.5,zp_upper=34.0, zp_fwhm_lower=-.5, zp_fwhm_upper=2.0,photprobmin=0.5):
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
            zsellist,isellist,cnites = des_utils.common_trignite_selector(zbandinfo,ibandinfo,zSNR_sel,iSNR_sel,1,0,1)
        elif fieldtype == 'deep':
            zsellist,isellist,cnites = des_utils.common_trignite_selector(zbandinfo,ibandinfo,zSNR_sel,iSNR_sel,1,1,1)
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
    for i in range(0,len(zbandinfo)):
        if np.any(zsellist[i]):
            imag = -2.5*np.log10(ibandinfo[i][2][isellist[i]][0])
            zmag = -2.5*np.log10(zbandinfo[i][2][zsellist[i]][0])
            colors[i] = imag - zmag
            anytrigs[i] = np.any(zsellist[i]) and np.any(isellist[i])
    dict1['colors'] = colors
    dict1['trigs'] = anytrigs
    dict1['headers'] = theheaders
    return dict1
