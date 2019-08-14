# Template builder
# Henryk T. Haniewicz, 2019

# Imports
import os
import numpy as np
import pickle
from astropy.io import fits
from pypulse.archive import Archive

file_root = os.path.dirname( os.path.abspath( __file__ ) )

# Master template class
class Template:

    """
    Master class dedicated to creating high SNR profiles for use in pulsar timing.
    """

    def __init__( self, psr_name, frontend, subbands, *dirs, template_dir = "templates", verbose = False ):

        """
        Template class

        Parameters
        ----------
        psr_name            : str
            Pulsar name in PSRFITS files
        frontend            : str
            Frontend (frequency band) of observations
        subbands            : int
            Number of sub-bands (F / T) to compart data into
        *dirs                  : str, [str, str, ...]
            Directories to look for PSRFITS files
        template_dir        : str, optional
            Directory to save final template
        verbose             : bool, optional
            Prints information to the console
        """



        self.psr_name = str( psr_name )
        self.frontend = str( frontend )
        self.subbands = subbands
        self.temp_dir = os.path.join( file_root, self.psr_name, template_dir )
        self.pkl_dir = os.path.join( file_root, self.psr_name, 'pickle_dumps' )

        self.check_directories()
        if len( dirs ) == 0:
            raise ValueError( "Must provide at least one directory to search through!" )
        else:
            self.dirs = dirs
        self.dirs = dirs
        self.verbose = verbose
        self.pklfile = os.path.join( self.pkl_dir, "{0}_{1}_{2}.pkl".format( self.psr_name, self.frontend, self.subbands ) )
        self.savefile = "{0}_{1}_{2}_template.npy".format( self.psr_name, self.frontend, self.subbands )

    def __repr__( self ):
        return "Template( psr_name = {}, frontend = {}, sub-bands = {} )".format( self.psr_name, self.frontend, self.subbands )

    def __str__( self ):
        return self.psr_name + "_" + self.frontend + "_" + str( self.subbands )

    def check_directories( self ):
        if not os.path.exists( os.path.join( file_root, self.psr_name ) ):
            os.makedirs( os.path.join( file_root, self.psr_name ) )
        if not os.path.exists( self.temp_dir ):
            os.makedirs( self.temp_dir )
        if not os.path.exists( self.pkl_dir ):
            os.makedirs( os.path.join( self.pkl_dir ) )
        return self


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
        if hdul[0].header[ 'OBS_MODE' ] != "PSR" or name != self.psr_name or fe != self.frontend:
            hdul.close()
            return -1
        hdul.close()

        ar = Archive( file, verbose = self.verbose )
        ar.tscrunch( nsubint = 1 )
        nsubint = ar.getNsubint()
        ar.fscrunch( nchan = 1 )
        n = 1
        nbin = ar.getNbin()
        data = ar.getData()

        return np.copy( data ), n, nbin


    def load_template( self ):

        """
        Attempts to load template from save. Otherwise, creates a new template with conventional name.
        """

        filename = self.pklfile

        try:
            pickle_in = open( filename, "rb" )
            temp_dict = pickle.load( pickle_in )
            pickle_in.close()
        except OSError:
            temp_dict = { 'DATA': None, 'IG_LIST': [] }
            pickle_out = open( filename, "wb" )
            pickle.dump( temp_dict, pickle_out )
            pickle_out.close()

        return temp_dict[ 'DATA' ], temp_dict[ 'IG_LIST' ], filename


    def make_template( self ):

        """
        Script to create the template
        """

        for directory in self.dirs:
            for f in sorted( os.listdir( directory ) ):
                template, ignore_list, pkl_name = self.load_template()

                if f in ignore_list:
                    if self.verbose:
                        print( "{} has already been added to the template.".format( f ) )
                    continue

                prep = self.prepare_file( os.path.join( directory, f ) )
                if prep == -1:
                    if self.verbose:
                        print( "Preparation of file {} failed. Skipping file...".format( f ) )
                    continue
                else:
                    d, c, b = prep[0], prep[1], prep[2]


                if template is None and c > 1:
                    if self.verbose:
                        print( "First time pass for subband > 1" )
                    self.template = np.zeros( (c, b), dtype = float )
                elif template is None:
                    if self.verbose:
                        print( "First time pass for subband = 1" )
                    self.template = np.zeros( b, dtype = float )
                else:
                    if self.verbose:
                        print( "Continuing template..." )
                    self.template = template

                if c > 1:
                    for i in range(c):
                        self.template[i] += d[i]
                else:
                    self.template += d


                if self.verbose:
                    print( "Template appended with data from {}".format( f ) )

                ignore_list.append( f )

                temp_dict = { 'DATA': self.template, 'IG_LIST': ignore_list }
                pickle_out = open( pkl_name, "wb" )
                pickle.dump( temp_dict, pickle_out )
                pickle_out.close()


        save_file = self.savefile
        np.save( os.path.join( self.temp_dir, save_file ), self.template )
        print( "Template creation finished." )
        return self.template


# Frequency dependent template class
class FD_Template( Template ):

    """
    Class dedicated to creating high SNR, frequency dependent, profiles for use in pulsar timing.
    """

    def __init__( self, psr_name, frontend, subbands, *dirs, template_dir = "templates", verbose = False ):

        """
        FD_Template (Frequency-Dependent Template) class

        Parameters
        ----------
        psr_name            : str
            Pulsar name in PSRFITS file
        frontend            : str
            Frontend (frequency band) of observations
        subbands            : int
            Number of frequency sub-bands to compart data into
        *dirs               : str, [str, str, ...]
            Directories to look for PSRFITS files
        template_dir        : str, optional
            Directory to save final template
        verbose             : bool, optional
            Prints information to the console
        """

        self.psr_name = str( psr_name )
        self.frontend = str( frontend )
        self.subbands = subbands

        self.temp_dir = os.path.join( file_root, self.psr_name, template_dir )
        self.pkl_dir = os.path.join( file_root, self.psr_name, 'pickle_dumps' )

        self.check_directories()

        if len( dirs ) == 0:
            raise ValueError( "Must provide at least one directory to search through!" )
        else:
            self.dirs = dirs
        self.verbose = verbose
        self.pklfile = os.path.join( self.pkl_dir, "{0}_{1}_nchan{2}.pkl".format( self.psr_name, self.frontend, self.subbands ) )
        self.savefile = "{0}_{1}_nchan{2}_template.npy".format( self.psr_name, self.frontend, self.subbands )


    def prepare_file( self, file ):

        """
        Prepares the PSRFITS file in the correct format for the program to use.
        """

        try:
            hdul = fits.open( file )
        except OSError:
            return -1

        name = hdul[0].header[ 'SRC_NAME' ]
        fe = hdul[0].header[ 'FRONTEND' ]
        if hdul[0].header[ 'OBS_MODE' ] != "PSR" or name != self.psr_name or fe != self.frontend:
            hdul.close()
            return -1
        hdul.close()

        ar = Archive( file, verbose = self.verbose )
        ar.tscrunch( nsubint = 1 )
        nsubint = ar.getNsubint()
        ar.fscrunch( nchan = self.subbands )
        nchan = ar.getNchan()
        nbin = ar.getNbin()
        data = ar.getData()

        return np.copy( data ), nchan, nbin


# Epoch dependent template class
class TD_Template( Template ):

    """
    Class dedicated to creating high SNR, time dependent, profiles for use in pulsar timing.
    """

    def __init__( self, psr_name, frontend, subbands, *dirs, template_dir = "templates", verbose = False ):

        """
        TD_Template (Time-Dependent Template) class

        Parameters
        ----------
        psr_name            : str
            Pulsar name in PSRFITS file
        frontend            : str
            Frontend (frequency band) of observations
        subbands            : int
            Number of epoch sub-bands to compart data into
        *dirs               : str, [str, str, ...]
            Directories to look for PSRFITS files
        template_dir        : str, optional
            Directory to save final template
        verbose             : bool, optional
            Prints information to the console
        """

        self.psr_name = str( psr_name )
        self.frontend = str( frontend )
        self.subbands = subbands

        self.temp_dir = os.path.join( file_root, self.psr_name, template_dir )
        self.pkl_dir = os.path.join( file_root, self.psr_name, 'pickle_dumps' )

        self.check_directories()
        if len( dirs ) == 0:
            raise ValueError( "Must provide at least one directory to search through!" )
        else:
            self.dirs = dirs
        self.verbose = verbose
        self.pklfile = os.path.join( pkl_dir, "{0}_{1}_nsubint{2}.pkl".format( self.psr_name, self.frontend, self.subbands ) )
        self.savefile = "{0}_{1}_nsubint{2}_template.npy".format( self.psr_name, self.frontend, self.subbands )


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
        if hdul[0].header[ 'OBS_MODE' ] != "PSR" or name != self.psr_name or fe != self.frontend:
            hdul.close()
            return -1
        hdul.close()

        ar = Archive( file, verbose = self.verbose )
        ar.tscrunch( nsubint = self.subbands )
        nsubint = ar.getNsubint()
        ar.fscrunch( nchan = 1 )
        nsubint = ar.getNsubint()
        nbin = ar.getNbin()
        data = ar.getData()

        return np.copy( data ), nsubint, nbin



if __name__ == "__main__":

    from custom_exceptions import ArgumentError
    import argparse

    def main():

        parser = argparse.ArgumentParser( formatter_class = argparse.RawDescriptionHelpFormatter,
                        prog = 'template_builder.py', description = '''\
                             PulseBlast Smart Template Builder
                       -------------------------------------------
                           Creates a very high SNR profile to
                                   create TOAs from
                            ''' )

        parser.add_argument( dest = 'psrname', nargs = 1, help = 'Name of the pulsar as given in the PSRFITS files. (e.g. J0737-3039)' )
        parser.add_argument( '-b', dest = 'frontend', nargs = 1, default = None, help = 'Frequency band of observation as given in the PSRFITS files.' )
        parser.add_argument( '-n', dest = 'subbands', nargs = '?', type = int, default = 1, help = 'Number of sub-bands to separate data into.' )
        parser.add_argument( '-o', dest = 'output_directory', nargs = '?', default = "templates", help = 'Path to save templates to (relative or global).' )
        parser.add_argument( '-d', dest = 'directories', nargs = '*', default = None, help = 'Directories to search for PSRFITS files in.' )
        parser.add_argument( '-t', dest = 'TD', action = 'store_true', default = False, help = 'Makes an epoch dependent template instead of a frequency dependent one.' )
        parser.add_argument( '-v', dest = 'verbose', action = 'store_true', default = False, help = 'Prints information to the console.' )

        args = parser.parse_args()


        # Checks argument requirements for non-optional flags (as a double check)
        if ( not args.frontend ):
            raise ArgumentError( "Frontend data required to create template." )
        if ( not args.directories ):
            raise ArgumentError( "At least one directory containing PSRFITS files is required." )


        # Initialize the template class object as normal and run the template creation script
        if args.TD:
            temp = TD_Template( args.psrname[0], args.frontend[0], args.subbands, *args.directories, template_dir = args.output_directory, verbose = args.verbose )
        else:
            temp = FD_Template( args.psrname[0], args.frontend[0], args.subbands, *args.directories, template_dir = args.output_directory, verbose = args.verbose )
        template = temp.make_template()


    main()
