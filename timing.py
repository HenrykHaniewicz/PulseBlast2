# PulseBlast2 timing file
# Henryk T. Haniewicz, 2019

"""
TODO:
Do timing
"""

# Local imports
import utils.pulsarUtilities as pu
import utils.plotUtils as pltu
import utils.otherUtilities as u
import utils.mathUtils as mathu
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
class Timer:

    """
    Base class for timing in PulseBlast
    """

    def __init__( self, psr_name, *dirs, temp_dir = "templates", saveddata_dir = "data", toa_dir = 'toas', tim_ext = 'tim', jump_flags = "", epochs = 1, subbands = 1, verbose = False ):

        self.psr_name = str( psr_name )
        self.temp_dir = os.path.join( file_root, self.psr_name, temp_dir )
        self.pkl_dir = os.path.join( file_root, self.psr_name, 'pickle_dumps' )
        self.saveddata_dir = os.path.join( file_root, self.psr_name, saveddata_dir )
        self.toa_dir = os.path.join( file_root, self.psr_name, toa_dir )

        # Check that necessary directories exist (if not, make them)
        self.check_directories()

        self.verbose = verbose
        if len( dirs ) == 0:
            self.dirs = [self.saveddata_dir]
        else:
            self.dirs = dirs
        self.pklfile = os.path.join( self.pkl_dir, "{}_timing_save.pkl".format( self.psr_name ) )
        self.epochs = epochs
        self.subbands = subbands
        self.tim_ext = tim_ext
        self.jump_flags = str( jump_flags )

    def __repr__( self ):
        return "Timer( psr_name = {}, epochs = {}, subbands = {} )".format( self.psr_name, self.epochs, self.subbands )

    def __str__( self ):
        return self.psr_name


    def check_directories( self ):
        if not os.path.exists( os.path.join( file_root, self.psr_name ) ):
            os.makedirs( os.path.join( file_root, self.psr_name ) )
        if not os.path.exists( self.temp_dir ):
            os.makedirs( self.temp_dir )
        if not os.path.exists( self.saveddata_dir ):
            os.makedirs( self.saveddata_dir )
        if not os.path.exists( self.toa_dir ):
            os.makedirs( self.toa_dir )
        if not os.path.exists( self.pkl_dir ):
            os.makedirs( os.path.join( self.pkl_dir ) )
        return self


    def load_template( self, dir, filename ):

        try:
            template = np.load( os.path.join( dir, filename ) )
        except FileNotFoundError:
            raise TemplateLoadError( "n = {} template not found in {}".format( self.subbands, dir ) )

        return template


    def prepare_file( self, file ):

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

        tmp_fn = "{0}_{1}_nchan{2}_template.npy".format( self.psr_name, fe, self.subbands )
        try:
            template = self.load_template( self.temp_dir, tmp_fn )
        except TemplateLoadError:
            print( "Template not found" )
            reply = str( input( "Would you like to make a suitable one? (y / n)" ) ).lower().strip()
            if reply[0] == 'y':
                temp = FD_Template( self.psr_name, fe, self.subbands, template_dir = "templates", verbose = self.verbose, *self.dirs )
                template = temp.make_template()
            else:
                raise TemplateLoadError( "You can make a suitable template via the following command: python template_builder.py psr_name -b [frontend] -d [dirs]" )

        ar = Archive( file, verbose = self.verbose )
        ar.tscrunch( nsubint = self.epochs )
        ar.fscrunch( nchan = self.subbands )

        return ar, template, fe, mjd

    def time( self ):

        """
        """

        if not self.verbose:
            sys.stdout.write( '\n {0:<7s}  {1:<7s}\n'.format( 'Files', '% done' ) )

        for directory in self.dirs:
        # Cycle through each file in the stored directory
            for f in sorted( os.listdir( directory ) ):

                prep = self.prepare_file( os.path.join( directory, f ) )

                if prep == -1:
                    if self.verbose:
                        try:
                            print( "Preparation of file {} failed. Skipping...".format( f ) )
                        except UnicodeEncodeError:
                            print( "Preparation of file {} failed. Skipping...".format( f.encode( 'utf-8' ) ) )
                    continue

                ar, temp, fe, mjd = prep[0], prep[1], prep[2], prep[3]

                save_fn = "{0}_{1}.{2}".format( self.psr_name, fe, self.tim_ext )
                abs_save = os.path.join( self.toa_dir, save_fn )

                ar.time( temp, filename = abs_save, MJD = True, flags = self.jump_flags, appendto = True )


if __name__ == "__main__":

    v = ["/Volumes/Physics_Group/pulsar/data/J1851+00/fold"]
    t = Timer( "J1851+00", *v, jump_flags = "-f TEST_JUMP", epochs = 2, subbands = 1, verbose = True )
    t.time()
