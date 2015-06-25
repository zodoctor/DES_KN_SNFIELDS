import numpy as np
from itertools import ifilter

dtypes = {
    'CID':'|S15',
    'MJD':float,
    'Tobs':float,
    'FLUXCAL':float,
    'FLUXCAL_ERR':float,
    'FLUXCALERR':float,
    'DATAFLAG':int,
    'BAND':'|S1',
    'CHI2':float,
    'IFIT':float,
    'FLT':'|S1',
    'FIELD':'|S2',
    'PHOTFLAG':int,
    'PHOTPROB':float,
    'PSF':float,
    'ZPFLUX':float,
    'TEMPLERR':float,
    'ZPT':float,
    'SNR':float,
    'MAG':float,
    'MAGERR':float,
    'SIM_MAGOBS':float,
    'RA':float,
    'DEC':float,
    'HOST_ZPHOT':float,
    'N_NITE_OBS':int,
    'N_NITE_TRIG':int,
    'g1_MJD':float,
    'g1_FLUX':float,
    'g1_FLUXERR':float,
    'r1_MJD':float,
    'r1_FLUX':float,
    'r1_FLUXERR':float,
    'i1_MJD':float,
    'i1_FLUX':float,
    'i1_FLUXERR':float,
    'z1_MJD':float,
    'z1_FLUX':float,
    'z1_FLUXERR':float,
    'g2_MJD':float,
    'g2_FLUX':float,
    'g2_FLUXERR':float,
    'r2_MJD':float,
    'r2_FLUX':float,
    'r2_FLUXERR':float,
    'i2_MJD':float,
    'i2_FLUX':float,
    'i2_FLUXERR':float,
    'z2_MJD':float,
    'z2_FLUX':float,
    'z2_FLUXERR':float,
    'g3_MJD':float,
    'g3_FLUX':float,
    'g3_FLUXERR':float,
    'r3_MJD':float,
    'r3_FLUX':float,
    'r3_FLUXERR':float,
    'i3_MJD':float,
    'i3_FLUX':float,
    'i3_FLUXERR':float,
    'z3_MJD':float,
    'z3_FLUX':float,
    'z3_FLUXERR':float,
    'g4_MJD':float,
    'g4_FLUX':float,
    'g4_FLUXERR':float,
    'r4_MJD':float,
    'r4_FLUX':float,
    'r4_FLUXERR':float,
    'i4_MJD':float,
    'i4_FLUX':float,
    'i4_FLUXERR':float,
    'z4_MJD':float,
    'z4_FLUX':float,
    'z4_FLUXERR':float}

def parse_observations(infile):
    with open(infile, 'r') as inp:
        (header,headerdict) = parse_header(inp)
        # Prepare data type for struc. array
        while 'VARLIST:' not in header:
            header = inp.readline().split()
        params = header
        cols = [i for i in range(len(params)) if params[i] in dtypes.keys()]
        dtype = [(params[col], dtypes[params[col]]) for col in cols]

        # Filter out anything that isn't an observation
        filtered_inp = ifilter(lambda line: line.startswith('OBS:'), inp)
        obs = np.genfromtxt(filtered_inp, usecols=cols, dtype=dtype)

    return obs,headerdict

def parse_header(inp):
    # get the header field names
    headerdict = dict()
    header = inp.readline().split()
    while (1==1):
        header = inp.readline().split()
        if not header:
            pass
        elif header[0]=='#':
            break
        else:
            headerdict[header[0][:-1]] = header[1]
    return header,headerdict
