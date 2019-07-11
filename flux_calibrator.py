# Flux calibrator

"""
TODO:
Main calobration logic
Smart file name handling for both load and save
    Try to load
Only works for Format 2 at the moment...

On / Off dict --- [{ 'MJD' : mjd, 'ON' : on_file, 'OFF' : off_file, 'FE' : fe, 'NUM' : num }]

DO NOT MAKE 'C++' CODE AGAIN!!!!!!! THIS IS PYTHON
"""

# Imports
import numpy as np
import os
import math
import pickle
import matplotlib.pyplot as plt

from pypulse.archive import Archive

from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as astu

from utils.calculate_flux import find_source_params_f2, getFlux

file_root = os.path.dirname( os.path.abspath( __file__ ) )


class FluxCalibrator:

    """
    Class for flux calibration

    Parameters
    ----------
    psr_name      : str
        Name of PSR as given in the PSRFITS files
    cont_name     : str
        Name of continuum source (or alias) as given in fluxcal.cfg
    *dirs         : str, os.Path, [str, ..., str], [os.Path, ..., os.Path], optional
        Directories to look for PSRFITS files in for calibration (default is saveddata_dir)
    cont_dir      : str, os.Path, optional
        Directory containing all continuum source cal files you plan to use (must be one directory).
    saveddata_dir : str, os.Path, optional
        Location to save newly calibrated PSRFITS files (local to this file)
    verbose       : bool, optional
        Displays more information to the console
    """

    def __init__( self, psr_name, cont_name, *dirs, cont_dir = None, saveddata_dir = "data/", verbose = False ):

        self.psr_name = str( psr_name )
        self.cont_name = str( cont_name )
        self.saveddata_dir = self.psr_name + '/' + saveddata_dir

        # Directory setup
        self.check_directories()
        if len( dirs ) == 0:
            self.dirs = [os.path.join( file_root, self.saveddata_dir )]
        else:
            self.dirs = dirs
        if cont_dir is None:
            self.cont_dir = os.path.join( file_root, self.saveddata_dir )
        else:
            if not os.path.exists( cont_dir ):
                raise FileNotFoundError( "Could not find continuum source directory: {}".format( cont_dir ) )
            self.cont_dir = cont_dir

        self.verbose = verbose
        self.pkl_dir = os.path.join( file_root, self.psr_name, 'pickle_dumps/' )
        self.pklfile = os.path.join( self.pkl_dir, "{}_calibration_save.pkl".format( self.psr_name ) )
        self.dirs = ["/Users/zhn11tau/Documents/DATA/1829+2456_2017/"]

    def __repr__( self ):
        return "FluxCalibrator( psr_name = {}, cont_name = {} )".format( self.psr_name, self.cont_name )

    def __str__( self ):
        return self.psr_name + "_" + self.cont_name

    def check_directories( self ):
        if not os.path.exists( os.path.join( file_root, self.psr_name ) ):
            os.makedirs( os.path.join( file_root, self.psr_name ) )
        if not os.path.exists( os.path.join( file_root, self.saveddata_dir ) ):
            os.makedirs( os.path.join( file_root, self.saveddata_dir ) )
        if not os.path.exists( os.path.join( file_root, self.psr_name, 'pickle_dumps' ) ):
            os.makedirs( os.path.join( file_root, self.psr_name, 'pickle_dumps' ) )
        if not os.path.exists( os.path.join( file_root, self.psr_name, 'pickle_dumps', 'calibration' ) ):
            os.makedirs( os.path.join( file_root, self.psr_name, 'pickle_dumps', 'calibration' ) )


    def load_session( self ):

        filename = self.pklfile

        try:
            pickle_in = open( filename, "rb" )
            load_dict = pickle.load( pickle_in )
            pickle_in.close()
        except OSError:
            if self.verbose:
                print( "Could not retrieve fluxcal saved data." )
                print( "Attempting to create a new save file..." )
            raise Exception( "Nah" )

    def hdul_setup( self, dir, file, is_cal = True ):

        hdul = fits.open( os.path.join( dir, file ) )

        root, ext = os.path.splitext( file )
        obs_num = root[-4:]

        mjd = hdul[0].header[ 'STT_IMJD' ]
        fe = hdul[0].header[ 'FRONTEND' ]
        obs_mode = hdul[0].header[ 'OBS_MODE' ]
        if (is_cal and obs_mode == "PSR") or (not is_cal and obs_mode == "CAL"):
            hdul.close()
            raise OSError( "Incorrect observation mode." )

        return hdul, mjd, fe, obs_num, obs_mode

    def get_onoff_list( self, tolerance = 1 ):

        dict_file = "{}_onoff_list.pkl".format( self.cont_name )
        abs_dict_file = os.path.join( self.pkl_dir, "calibration", dict_file )

        if os.path.isfile( abs_dict_file ):
            if self.verbose:
                print( "Loading previously saved continuum data..." )
            pickle_in = open( abs_dict_file, "rb" )
            onoff_list = pickle.load( pickle_in )
            pickle_in.close()
        else:

            if self.verbose:
                print( "Making new continuum data list..." )

            onoff_list = []

            pos, params = find_source_params_f2( self.cont_name )
            m_coordinates = SkyCoord( "{0} {1}".format( pos[0], pos[1] ), unit = ( astu.hourangle, astu.degree ) )

            for file in sorted( os.listdir( self.cont_dir ) ):

                try:
                    hdul, mjd, fe, obs_num, obs_mode = self.hdul_setup( self.cont_dir, file )
                except OSError:
                    continue

                coords = SkyCoord( "{0} {1}".format( hdul[0].header[ 'RA' ], hdul[0].header[ 'DEC' ] ), unit = ( astu.hourangle, astu.degree ) )
                hdul.close()

                onoff_list.append( { 'MJD' : mjd, 'ON' : None, 'OFF' : None, 'FE' : fe, 'NUM' : obs_num } )

                if m_coordinates.separation( coords ) <= ( tolerance * astu.arcmin ):
                    mode = 'ON'
                else:
                    mode = 'OFF'

                for dict in onoff_list:
                    if (dict[ 'MJD' ] == mjd) and (dict[ 'FE' ] == fe) and (dict[ 'NUM' ] == obs_num) and (dict[ mode ] is None):
                        dict[ mode ] = file

            for dict in reversed( onoff_list ):
                if (dict[ 'ON' ] is None) or (dict[ 'OFF' ] is None):
                    onoff_list.remove( dict )

            if self.verbose:
                print( "Saving as {}".format( dict_file ) )

            pickle_out = open( abs_dict_file, "wb" )
            pickle.dump( onoff_list, pickle_out )
            pickle_out.close()

        return onoff_list

    def get_closest_contfile( self, mjd_tol = 50 ):

        onoff_list = self.get_onoff_list( tolerance = 1 )

        a = []

        for directory in self.dirs:
            for psr_file in sorted( os.listdir( directory ) ):

                mjd_list = []

                try:
                    hdul, psr_mjd, psr_fe, obs_num, obs_mode = self.hdul_setup( directory, psr_file )
                    if self.verbose:
                        print( "Opening {}".format( psr_file ) )
                except OSError:
                    if self.verbose:
                        print( "Couldn't open {}".format( psr_file ) )
                    continue

                psr_file = os.path.join( directory, psr_file )

                if obs_mode == 'CAL':

                    for dict in onoff_list:
                        if dict[ 'FE' ] == psr_fe:
                            delta_mjd = abs( dict[ 'MJD' ] - psr_mjd )
                            if all( elem > delta_mjd for elem in mjd_list ):
                                mjd_list.append( delta_mjd )
                                if delta_mjd < mjd_tol:
                                    on = os.path.join( self.cont_dir, dict[ 'ON' ] )
                                    off = os.path.join( self.cont_dir, dict[ 'OFF' ] )
                                else:
                                    on, off = [None, None]
                    a.append( [psr_file, on, off] )

        return a


    def calculate_Fcal_functions( self, cal_file_list ):

        if type( cal_file_list ) != np.ndarray:
            cal_file_list = np.array( cal_file_list )

        if cal_file_list.ndim != 1:
            raise ValueError( "Should be a vector" )

        archives = []
        for f in cal_file_list:
            archives.append( Archive( f, verbose = self.verbose ) )

        polarizations = []

        for t in archives:
            t.tscrunch()
            A, B, C, D = self.convert_subint_pol_state( t.getData(), t.subintheader[ 'POL_TYPE' ], "AABBCRCI", linear = t.header[ 'FD_POLN' ] )
            polarizations.append( [ A, B ] )

        polarizations = np.array( polarizations )
        psr_cal, on_cal, off_cal = archives



        return self

    def cal_to_count( self ):


        return


    def calibrate( self ):



        return


    def convert_subint_pol_state( self, subint, input, output, linear = 'LIN' ):

        if input == output:
            out_S = subint
        elif input == "AABBCRCI" and output == "IQUV": # Coherence -> Stokes
            A, B, C, D = subint
            if linear == 'LIN':
                I = A+B
                Q = A-B
                U = 2*C
                V = 2*D
            else:
                I = A+B
                Q = 2*C
                U = 2*D
                V = A-B
            out_S = [ I, Q, U, V ]
        elif input == "IQUV" and output == "AABBCRCI": # Stokes -> Coherence
            I, Q, U, V = subint
            if linear == 'LIN':
                pass
                A = (I+Q)/2.0
                B = (I-Q)/2.0
                C = U/2.0
                D = V/2.0
            else:
                A = (I+V)/2.0
                B = (I-V)/2.0
                C = Q/2.0
                D = U/2.0
            out_S = [ A, B, C, D ]
        if type( subint ) == np.ndarray:
            out_S = np.array( out_S )

        return out_S




if __name__ == "__main__":

    n = "J1829+2456"
    cn = "B1442"
    c_d = "/Users/zhn11tau/Documents/DATA/cont/"
    c = FluxCalibrator( n, cn, cont_dir = c_d, verbose = False )
    cont_match = c.get_closest_contfile()
    for e in cont_match:
        c.calculate_Fcal_functions( e )
