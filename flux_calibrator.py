# Flux calibrator
# Henryk T. Haniewicz, 2019

"""
TODO:
Main calobration logic
Smart file name handling for both load and save
    Try to load
Only works for Format 2 at the moment...

DO NOT MAKE 'C++' CODE AGAIN!!!!!!! THIS IS PYTHON

CAL_DICT_LIST --- [ [{ 'ARC' : arc, 'DATA' : aabb_data, 'S_DUTY' : arc.getValue( 'CAL_PHS' ) , 'DUTY' : arc.getValue( 'CAL_DCYC' ), 'BW' : arc.getBandwidth(), 'CTR_F' : arc.getCenterFrequency( weighted = True ) }, {  }, {  }], [{}{}{}] ]
"""

# Imports
import numpy as np
import os
import math
import pickle
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

from pypulse.archive import Archive

from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as astu

from utils.calculate_flux import find_source_params_f2, getFlux
from utils.saving import save_psrfits, save_session
from utils.loading import load_session as load

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

    def __init__( self, psr_name, cont_name, *dirs, cont_dir = None, saveddata_dir = "data", verbose = False ):

        self.psr_name = str( psr_name )
        self.cont_name = str( cont_name )
        self.saveddata_dir = os.path.join( self.psr_name, saveddata_dir )

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
        self.pkl_dir = os.path.join( file_root, self.psr_name, 'pickle_dumps' )
        self.pklfile = os.path.join( self.pkl_dir, "{}_calibration_save.pkl".format( self.psr_name ) )

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


    def load_session( self, bin_file ):
        return load( bin_file, mode = 'c' )

    def save_position( self, bin_file, *ex ):
        return save_session( bin_file, *ex, mode = 'c' )

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
            onoff_list = self.load_session( abs_dict_file )
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

            self.save_position( abs_dict_file, onoff_list )

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
                    a.append( [psr_file, on, off, psr_mjd] )

        return a


    def calculate_Jy_per_count( self, cal_file_list, gain = 11.0 ):

        """
        Input list: [ PSR_CAL, ON_CAL, OFF_CAL, CAL_MJD ]

        Returns:
        conversion_factor  :  np.ndarray
        """

        G = gain

        if type( cal_file_list ) != np.ndarray:
            cal_file_list = np.array( cal_file_list )

        if cal_file_list.ndim != 1:
            raise ValueError( "Should be a vector" )

        archives = []
        freqs = []
        for f in cal_file_list[:-1]:
            hdul = fits.open( f )
            freqs.append( hdul[3].data[ 'DAT_FREQ' ][0] )
            hdul.close()
            archives.append( Archive( f, verbose = self.verbose ) )


        aabb_list = []

        for i, arc in enumerate( archives ):
            arc.reset()
            arc.tscrunch()
            A, B, C, D = self.convert_subint_pol_state( arc.getData( weight = False ), arc.subintheader[ 'POL_TYPE' ], "AABBCRCI", linear = arc.header[ 'FD_POLN' ] )
            l = { 'ARC' : arc, 'DATA' : [ A, B ], 'FREQS' : freqs[i], 'S_DUTY' : arc.getValue( 'CAL_PHS' ) , 'DUTY' : arc.getValue( 'CAL_DCYC' ), 'BW' : arc.getBandwidth() }
            aabb_list.append( l )


        H, L, T0 = self._prepare_calibration( aabb_list )
        F_ON = ( H[1]/L[1] ) - 1
        F_OFF = ( H[2]/L[2] ) - 1

        C0 = T0[1:] / ( ( 1 / F_ON ) - ( 1 / F_OFF ) )
        T_sys = C0 / F_OFF
        F_cal = ( T_sys * F_OFF ) / G


        Fa, Fb = interp1d( aabb_list[1][ 'FREQS' ], F_cal[0][0], kind='cubic', fill_value = 'extrapolate' ), interp1d( aabb_list[2][ 'FREQS' ], F_cal[0][1], kind='cubic', fill_value = 'extrapolate' )

        conversion_factor = [ np.array(Fa( aabb_list[0][ 'FREQS' ] ) / ( H[0][0] - L[0][0] )), np.array( Fb( aabb_list[0][ 'FREQS' ] ) / ( H[0][1] - L[0][1] ) ) ]
        conversion_factor = np.array( conversion_factor )

        return conversion_factor


    def _prepare_calibration( self, archive_list, r_err = 8 ):

        H = []
        L = []
        T0 = []

        for dict in archive_list:
            all_high_means = []
            all_low_means = []
            T0_pol = []

            for pol in dict[ 'DATA' ]:
                high_means = []
                low_means = []
                T0_chan = []

                for i, channel in enumerate( pol ):

                    flux = getFlux( float( dict[ 'FREQS' ][i]/1000 ), self.cont_name, False )

                    start_bin = math.floor( len( channel ) * dict[ 'S_DUTY' ] )
                    mid_bin = math.floor( len( channel ) * ( dict[ 'S_DUTY' ] + dict[ 'DUTY' ] ) )
                    end_bin = mid_bin + ( math.floor( len( channel ) * dict[ 'DUTY' ] ) )
                    bin_params = [ start_bin, mid_bin, end_bin ]

                    low_mean = np.mean( channel[ bin_params[0] : bin_params[1] ] )
                    high_mean = np.mean( channel[ bin_params[1] : bin_params[2] ] )

                    low_mean = round( low_mean, r_err )
                    high_mean = round( high_mean, r_err )

                    high_means.append( high_mean )
                    low_means.append( low_mean )
                    T0_chan.append( flux )

                all_high_means.append( high_means )
                all_low_means.append( low_means )
                T0_pol.append( T0_chan )

            H.append( all_high_means )
            L.append( all_low_means )
            T0.append( T0_pol )

        H = np.array(H)
        L = np.array(L)
        T0 = np.array(T0)

        return H, L, T0


    def calibrate( self ):

        """
        Master calibration method
        """

        conv_file = "{}_{}_fluxcalibration_conversion_factors.pkl".format( self.psr_name, self.cont_name )
        cal_mjd_file = "{}_{}_fluxcalibration_cal_mjds.pkl".format( self.psr_name, self.cont_name )
        conv_abs_path, cal_abs_path = os.path.join( self.pkl_dir, 'calibration', conv_file ), os.path.join( self.pkl_dir, 'calibration', cal_mjd_file )

        if os.path.isfile( conv_abs_path ):
            if self.verbose:
                print( "Loading previously saved conversion factor data..." )
            conversion_factors = self.load_session( conv_abs_path )
            cal_mjds = self.load_session( cal_abs_path )
        else:
            if self.verbose:
                print( "Making new conversion factor list..." )

            conversion_factors = []
            cal_mjds = []
            for e in self.get_closest_contfile():
                conversion_factors.append( self.calculate_Jy_per_count( e ) )
                cal_mjds.append( e[3] )

            conversion_factors = np.array( conversion_factors )

            if self.verbose:
                print( "Saving as {}".format( conv_file ) )

            self.save_position( conv_abs_path, conversion_factors )
            self.save_position( cal_abs_path, cal_mjds )


        if type( conversion_factors ) != np.ndarray:
            conversion_factors = np.array( conversion_factors )
        if type( cal_mjds ) != np.ndarray:
            cal_mjds = np.array( cal_mjds )

        print(conversion_factors)

        counter = 0

        for directory in self.dirs:
            for psr_file in sorted( os.listdir( directory ) ):
                try:
                    hdul, psr_mjd, psr_fe, obs_num, obs_mode = self.hdul_setup( directory, psr_file, False )
                    if self.verbose:
                        print( "Opening {}".format( psr_file ) )
                except OSError:
                    if self.verbose:
                        try:
                            print( "Couldn't open {}".format( psr_file ) )
                        except UnicodeEncodeError:
                            print( "Couldn't open {}".format( psr_file.encode( "utf-8" ) ) )
                    continue

                if obs_mode != "PSR":
                    continue

                ar = Archive( os.path.join( directory, psr_file ), verbose = self.verbose )
                ar.reset()
                data = ar.getData( weight = False )
                new_data = []
                for sub in data:
                    A, B, C, D = self.convert_subint_pol_state( sub, ar.subintheader[ 'POL_TYPE' ], "AABBCRCI", linear = ar.header[ 'FD_POLN' ] )
                    new_data.append( [ A, B ] )

                new_data = np.array( new_data )


                while psr_mjd != cal_mjds[ counter ]:
                    print( psr_mjd, cal_mjds[ counter ] )
                    counter += 1
                    if counter >= len( conversion_factors ):
                        break
                else:
                    for sub in new_data:
                        sub = conversion_factors[ counter ] * sub
                        print(sub.shape)
                    counter = 0



                save_psrfits(  )

        return self


    def convert_subint_pol_state( self, subint, input, output, linear = 'LIN' ):

        """

        """

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
        else:
            print("WTF")


        out_S = np.array( out_S )

        return out_S




if __name__ == "__main__":

    n = "J1829+2456"
    cn = "B1442"
    d = [ "/Users/zhn11tau/Documents/DATA/J1829+2456/1829+2456_2017", "/Users/zhn11tau/Documents/DATA/J1829+2456/1829+2456_2018_0001" ]
    c_d = "/Users/zhn11tau/Documents/DATA/J1829+2456/cont/"
    c = FluxCalibrator( n, cn, *d, cont_dir = c_d, verbose = True )
    c.calibrate()
