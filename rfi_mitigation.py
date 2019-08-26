# RFI mitigation classes
# Henryk T. Haniewicz, 2019

"""
TODO:
Create four RFI mitigation classes based on a base class
    User enforced channel zapping
    Curve fit and sigma clipping (S)
    Bayesian mitigation (B)
    DLNN 1D image recognition (N)
Modular saving (with detailed ignore list)
    save_dict = { 'FILE' : filename, 'DATA' : self.data, 'POS' : [ s, c ], 'IG_LIST' : [ {}, ... ] }
"""

FIT = True

# Local imports
import utils.pulsarUtilities as pu
import utils.plotUtils as pltu
import utils.otherUtilities as u
import utils.mathUtils as mathu
from utils.saving import save_psrfits
from utils.saving import load_session as load
from custom_exceptions import TemplateLoadError
from template_builder import FD_Template

import matplotlib.pyplot as plt
import scipy.stats as spyst

# Imports
from pypulse.archive import Archive
import numpy as np
import pickle
from astropy.io import fits
import os

TOL = 1e-7

file_root = os.path.dirname( os.path.abspath( __file__ ) )

# Base RFI excision class
class RFIBlaster:

    """
    Base class for RFI mitigation in PulseBlast

    Initializing this base class and mitigating will return the data as input.
    """

    def __init__( self, psr_name, *dirs, iterations = 1, temp_dir = "templates", saveddata_dir = "data", epoch_avg = False, save_as_np = False, verbose = False ):

        self.psr_name = str( psr_name )
        self.temp_dir = os.path.join( file_root, self.psr_name, temp_dir )
        self.pkl_dir = os.path.join( file_root, self.psr_name, 'pickle_dumps' )
        self.saveddata_dir = os.path.join( file_root, self.psr_name, saveddata_dir )

        # Check that necessary directories exist (if not, make them)
        self.check_directories()

        self.verbose = verbose
        if len( dirs ) == 0:
            self.dirs = [self.saveddata_dir]
        else:
            self.dirs = dirs
        self.pklfile = os.path.join( self.pkl_dir, "{}_rfimitigation_save.pkl".format( self.psr_name ) )
        self.iterations = iterations
        self.epoch_average = epoch_avg
        self.save_as_np = save_as_np
        self.method = self.get_method()

    def __repr__( self ):
        return "RFIBlaster( psr_name = {} )".format( self.psr_name )

    def __str__( self ):
        return self.psr_name

    def check_directories( self ):
        if not os.path.exists( os.path.join( file_root, self.psr_name ) ):
            os.makedirs( os.path.join( file_root, self.psr_name ) )
        if not os.path.exists( self.temp_dir ):
            os.makedirs( self.temp_dir )
        if not os.path.exists( self.saveddata_dir ):
            os.makedirs( self.saveddata_dir )
        if not os.path.exists( self.pkl_dir ):
            os.makedirs( os.path.join( self.pkl_dir ) )
        return self


    def get_method( self ):
        return 'NO METHOD'


    def load_template( self, dir, filename ):

        try:
            template = np.load( os.path.join( dir, filename ) )
        except FileNotFoundError:
            raise TemplateLoadError( "n = 1 template not found in {}".format( dir ) )

        return template

    def load_session( self ):
        return load( self.pklfile, arg = 'r' )

    def save_position( self, file, data, position, ignore_list ):

        if len( position ) != 2:
            raise ValueError( "Position data corrupted. Should be length 2. Actual length: {}".format( len( position ) ) )

        save_dict = { 'FILE' : file, 'DATA': data, 'POS' : position, 'IG_LIST': ignore_list }
        pickle_out = open( self.pklfile, "wb" )
        pickle.dump( save_dict, pickle_out )
        pickle_out.close()

        return save_dict


    def prepare_file( self, file, do_fit = False ):

        """
        Prepares the PSRFITS file in the correct format for the program.
        """

        try:
            hdul = fits.open( file )
        except OSError:
            return -1

        name = hdul[0].header[ 'SRC_NAME' ]
        fe = hdul[0].header[ 'FRONTEND' ]
        mjd = hdul[0].header[ 'STT_IMJD' ]
        if hdul[0].header[ 'OBS_MODE' ] != "PSR" or name != self.psr_name:
            hdul.close()
            return -1
        hdul.close()

        tmp_fn = "{0}_{1}_nchan1_template.npy".format( self.psr_name, fe )
        try:
            template = self.load_template( self.temp_dir, tmp_fn )
        except TemplateLoadError:
            print( "Template not found" )
            reply = str( input( "Would you like to make a suitable one? (y / n)" ) ).lower().strip()
            if reply[0] == 'y':
                temp = FD_Template( self.psr_name, fe, 1, template_dir = "templates", verbose = self.verbose, *self.dirs )
                template = temp.make_template( gaussian_fit = do_fit )
            else:
                raise TemplateLoadError( "You can make a suitable template via the following command: python template_builder.py psr_name -b [frontend] -d [dirs]" )

        ar = Archive( file, verbose = self.verbose )
        if self.epoch_average:
            ar.tscrunch( nsubint = 1 )

        return ar, template, fe, mjd


    def zap( self, file, archive, p, ignore_list, index ):

        if isinstance( index, int ):
            if self.verbose:
                print( "Zapping channel {}".format( index ) )
            archive.setWeights( 0.0, f = index )
        elif isinstance( index, ( list, tuple, np.ndarray ) ):
            for elem in index:
                if p[1] > elem:
                    print( p[1], elem )
                    continue
                save = self.save_position( file, archive.getData(), [0, elem], ignore_list )
                if self.verbose:
                    print( "Zapping channel {}".format( elem ) )
                archive.setWeights( 0.0, f = elem )
        data = archive.getData()

        return data, archive


    def mitigate( self, file, p, template, archive, ignore_list ):

        archive = archive
        mu = 0
        sigma = 1
        data = None

        return archive, mu, sigma, data


    def mitigation_setup( self ):

        for directory in self.dirs:

            if not os.path.exists( directory ):
                continue

            for f in sorted( os.listdir( directory ) ):

                root, ext = os.path.splitext( f )
                obs_num = root[-4:]

                last_file, data, p, ignore_list = self.load_session()
                ig_dict = { 'FILE' : f, 'METHOD' : self.method }

                ig = False
                for dict in ignore_list:
                    if (f == dict[ 'FILE' ]) and (dict[ 'METHOD' ] == self.method ):
                        if self.verbose:
                            print( "{} has already been excised using method {}".format( f, self.method ) )
                        ig = True
                        break
                if ig:
                    continue

                if (last_file is not None) and (last_file != f):
                    if self.verbose:
                        print( "File {} does not match file in saved data. Skipping...".format( f ) )
                    continue

                prep = self.prepare_file( os.path.join( directory, f ), do_fit = FIT )
                if prep == -1:
                    if self.verbose:
                        try:
                            print( "Preparation of file {} failed. Skipping...".format( f ) )
                        except UnicodeEncodeError:
                            print( "Preparation of file {} failed. Skipping...".format( f.encode( 'utf-8' ) ) )
                    continue

                ar, template, fe, mjd = prep[0], prep[1], prep[2], prep[3]

                # Start mitigation

                for it in np.arange( self.iterations ):
                    ar, mu, sigma, data = self.mitigate( f, p, template, ar, ignore_list )
                    if self.verbose:
                        print( "Data loaded for iteration {}".format( it + 1 ) )
                    try:
                        if (abs( old_mu - mu ) < TOL) and (abs( old_s - sigma ) < TOL):
                            if self.verbose:
                                print( "Stopping after iteration {} as data is fully excised.".format( it + 1 ) )
                            break
                    except NameError:
                        pass
                    old_mu, old_s = mu, sigma

                # End mitigation

                save_fn = "{0}_{1}_{2}_{3}".format( self.psr_name, mjd, fe, obs_num )
                if self.save_as_np:
                    np.save( os.join( self.saveddata_dir, save_fn ), data )
                    save_fn += "_DATWTS"
                    np.save( os.join( self.saveddata_dir, save_fn ), ar.getWeights() )
                else:
                    save_fn += ext
                    save_psrfits( os.path.join( self.saveddata_dir, save_fn ), ar )
                ignore_list.append( ig_dict )
                save_dict = self.save_position( None, None, [0, 0], ignore_list )

                if self.verbose:
                    print( "{0} fully mitigated using method {1}.".format( f, self.method ) )

        return self


# Bayesian heirarchy class
class Bayesian_Mitigator( RFIBlaster ):

    """
    Class for Bayesian heirarchy mitigation

    Parameters
    ----------
    psr_name      : str
        Name of PSR as given in the PSRFITS files
    *dirs         : str, os.Path, [str, ..., str], [os.Path, ..., os.Path], optional
        Directories to look for files with which to mitigate RFI (default is saveddata_dir)
    iterations    : int, optional
        Number of mitigation iterations to conduct
    temp_dir      : str, os.Path, optional
        Location to save / load templates to / from (local to this file)
    saveddata_dir : str, os.Path, optional
        Location to save newly RFI excised PSRFITS files (local to this file)
    epoch_avg     : bool, optional
        Determines whether profiles should be day averaged before mitigation
    verbose       : bool, optional
        Displays more information to the console
    """

    def get_method( self ):
        return 'B'

    def mitigate( self, file, p, template, archive, ignore_list ):

        archive = archive
        mu = 0
        sigma = 1
        data = None
        raise TypeError( "Bayesian excision features are not currently implemented" )

        return archive, mu, sigma, data


# DLNN class
class NN_Mitigator( RFIBlaster ):

    """
    Class for neural network image recognition mitigation

    Parameters
    ----------
    psr_name      : str
        Name of PSR as given in the PSRFITS files
    *dirs         : str, os.Path, [str, ..., str], [os.Path, ..., os.Path]
        Directories to look for files with which to mitigate RFI
    iterations    : int, optional
        Number of mitigation iterations to conduct
    temp_dir      : str, os.Path, optional
        Location to save / load templates to / from (local to this file)
    saveddata_dir : str, os.Path, optional
        Location to save newly RFI excised PSRFITS files (local to this file)
    epoch_avg     : bool, optional
        Determines whether profiles should be day averaged before mitigation
    verbose       : bool, optional
        Displays more information to the console
    """

    def get_method( self ):
        return 'N'

    def mitigate( self, file, p, template, archive, ignore_list ):

        archive = archive
        mu = 0
        sigma = 1
        data = None
        raise TypeError( "Neural Network features are not currently implemented" )

        return archive, mu, sigma, data


# Sigma-clipping class
class SigmaClip_Mitigator( RFIBlaster ):

    """
    Class for sigma clipping mitigation

    Parameters
    ----------
    psr_name      : str
        Name of PSR as given in the PSRFITS files
    *dirs         : str, os.Path, [str, ..., str], [os.Path, ..., os.Path]
        Directories to look for files with which to mitigate RFI
    iterations    : int, optional
        Number of mitigation iterations to conduct
    temp_dir      : str, os.Path, optional
        Location to save / load templates to / from (local to this file)
    saveddata_dir : str, os.Path, optional
        Location to save newly RFI excised PSRFITS files (local to this file)
    epoch_avg     : bool, optional
        Determines whether profiles should be day averaged before mitigation
    verbose       : bool, optional
        Displays more information to the console
    """

    def get_method( self ):
        return 'S'

    def mitigate( self, file, p, template, archive, ignore_list ):

        data = archive.getData()

        templateMask = pu.get_1D_OPW_mask( template, windowsize = (archive.getNbin() - 150) )
        rmsArray, linearRmsArray, mu, sigma = u.getRMSArrayProperties( data, templateMask, 1.0 ) # Needs to input 2D array

        # Creates the histogram
        pltu.histogram_and_curves( linearRmsArray, mean = mu, std_dev = sigma, bins = (archive.getNchan() // 4), x_axis = 'Root Mean Squared', y_axis = 'Frequency Density', title = r'$\mu={},\ \sigma={}$'.format( mu, sigma ), show = False, curve_list = [spyst.norm.pdf] )

        rejectionCriterion = mathu.chauvenet( rmsArray, mu, sigma, 2 )

        archive.reset()

        if self.verbose:
            print( "Rejection criterion created." )

        data, archive = self.zap( file, archive, p, ignore_list, np.asarray( rejectionCriterion ).nonzero()[ rejectionCriterion.ndim - 1 ] )

        return archive, mu, sigma, data

if __name__ == "__main__":

    np.set_printoptions( threshold = np.inf )
    v = ["/Users/zhn11tau/Documents/DATA/J1829+2456/1829+2456_2017"]
    s = SigmaClip_Mitigator( "J1829+2456", *v, iterations = 1, epoch_avg = True, verbose = True )
    s.mitigation_setup()
