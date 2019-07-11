import os
from astropy.io import fits

direc = "/Users/zhn11tau/Documents/DATA/cont/"

def create_directory_dict( dir, on_ra ):
    individual_dicts = []


    for file in sorted( os.listdir( dir ) ):

        fileF = dir + file

        try:
            hdul = fits.open( fileF )
        except OSError:
            print( "File {} did not match ASCII signature required for a fits file".format( file ) )
            continue

        # Get the frequency band used in the observation.
        try:
            frontend = hdul[0].header[ 'FRONTEND' ]
        except OSError:
            print( "Could not find any frontend information in file {}".format( file ) )
            hdul.close()
            continue

        # Get the observation mode.
        try:
            obs = hdul[0].header[ 'OBS_MODE' ]
        except OSError:
            print( "Could not find any observation mode information in file {}".format( file ) )
            hdul.close()
            continue

        # Get the start MJD
        try:
            cal_mjd = hdul[0].header[ 'STT_IMJD' ]
        except OSError:
            print( "Could not find any MJD information in file {}".format( file ) )
            hdul.close()
            continue


        # Get the cal mode.
        try:
            ra = hdul[0].header[ 'RA' ]
        except OSError:
            print( "Could not find any RA information in file {}".format( file ) )
            hdul.close()
            continue

        hdul.close()

        if obs == 'CAL':

            if ra == on_ra:
                cal_mode = 'ON'
            else:
                cal_mode = 'OFF'

            individual_dicts.append({ "FILE": file, "MODE": cal_mode, "FE": frontend, "MJD": cal_mjd })


        else:
            continue


    return individual_dicts

def create_ordered_onoff_list( direc_dict ):

    l_bands = []
    b_430 = []
    on_off_list_l = []
    on_off_list_430 = []

    for elem in direc_dict:
        if elem['FE'] == "lbw":
            l_bands.append( elem )
        if elem['FE'] == "430":
            b_430.append( elem )

    for x in [ l_bands, b_430 ]:

        for a in x:
            for b in x:
                if a == b:
                    continue
                if a['MODE'] == "OFF":
                    continue
                if b['MODE'] == "ON":
                    continue
                if (a['MJD'] == b['MJD']) and (a['FILE'][30:] == b['FILE'][30:]):
                    if x == l_bands:
                        on_off_list_l.append({ "ON": a['FILE'], "OFF": b['FILE'], "FE": a['FE'], "MJD": a['MJD'] })
                    else:
                        on_off_list_430.append({ "ON": a['FILE'], "OFF": b['FILE'], "FE": a['FE'], "MJD": a['MJD'] })

    return on_off_list_l, on_off_list_430

def get_closest_calibrator( band, oool, psr_mjd ):

    b = []
    if band == "lbw":
        list = oool[0]
    elif band == "430":
        list = oool[1]

    for a in list:
        delta_mjd = abs( a['MJD'] - psr_mjd )
        if all( elem > delta_mjd for elem in b ):
            b.append( delta_mjd )
            on = a['ON']
            off = a['OFF']

    return on, off
