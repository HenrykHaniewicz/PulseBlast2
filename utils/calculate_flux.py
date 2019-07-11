# Calculates flux from a known continuum based on one of two formats

import os.path as osp
import mmap
import math
import argparse

config_file = 'fluxcal.cfg'
config_abs = osp.join( osp.dirname( osp.abspath( __file__ ) ), config_file )


def find_flux_f1( frequency, source ):

    """
    Returns the Format 1 (see .cfg file) flux as a float from a given source.
    """

    if not isinstance( source, str ):
        raise TypeError( "Source parsed in must be a string" )

    new_f = float( frequency )

    src = source.encode()
    amp_line = b""


    # Opens cfg file and creates read-only mmap. Then, searches for '%' (start of new source), before then searching for aliases.
    # Once the alias is found, the last recorded '%' line is parsed into memory and the file closes.

    with open( config_abs, 'rb', 0 ) as file, mmap.mmap( file.fileno(), 0, access = mmap.ACCESS_READ ) as s:
        for lin_num, line in enumerate( iter(s.readline, b"") ):
            if line.find( b'%' ) != -1:
                amp_line = line

            if src in line and not line.startswith( b"#" ):
                file.close()
                break

            if line.find( b'Format 2' ) != -1:
                file.close()
                break

    amp_line = amp_line.decode()
    if amp_line == "":
        raise ValueError( "No source matching name given was found" )
    name, ra, dec, freq, flux, spec = amp_line.split()

    freq, flux, spec = float( freq ), float( flux ), float( spec )
    freq /= 1000

    # Use spectral index to calculate flux at a different frequency

    new_flux = flux * (( new_f / freq )**spec)
    return new_flux


def find_source_params_f2( source ):

    """
    Returns a list of Format 2 (see .cfg file) flux parameters for a given continuum source as strings.
    """

    if not isinstance( source, str ):
        raise TypeError( "Source parsed in must be a string" )

    src = source.encode()
    amp_line = b""


    # Opens cfg file and creates read-only mmap. Then, searches for '&' (start of new source), before then searching for aliases.
    # Once the alias is found, the last recorded '&' line is parsed into memory and the file closes.

    with open( config_abs, 'rb', 0 ) as file, mmap.mmap( file.fileno(), 0, access = mmap.ACCESS_READ ) as s:
        for lin_num, line in enumerate( iter(s.readline, b"") ):
            if line.find( b'&' ) != -1:
                amp_line = line

            if src in line and not line.startswith( b"#" ):
                file.close()
                break

    amp_line = amp_line.decode()
    if amp_line == "":
        raise ValueError( "No source matching name given was found" )
    params = amp_line.split()
    return params[1:3], params[3:]


def calculate_flux_f2( frequency, params ):

    """
    Returns the Format 2 (see .cfg file) flux at a given frequency for a set of coefficients.
    """

    # Turn all inputs to floats if not already
    f = float( frequency )
    p = list( map( float, params ) )

    LogS = p[0]

    for i, elem in enumerate( p[1:] ):
        corr = elem*((math.log10(f))**(i+1))
        LogS += corr

    flux = 10**LogS
    return flux


def getFlux( frequency, source, format1 = False ):

    """
    Master method. Returns the flux of a given frequency, source and format.
    """

    if format1:
        flux = find_flux_f1( frequency, str( source ) )
    else:
        position, params = find_source_params_f2( str( source ) )
        flux = calculate_flux_f2( frequency, params )

    return flux


if __name__ == "__main__":

    def parser( progname ):

        """
        Initializes argparse and collects the command line arguments.
        Returns a list of input arguments.
        """

        # Initialize the parser
        parser = argparse.ArgumentParser( formatter_class = argparse.RawDescriptionHelpFormatter,
                    prog = progname,
                    description = '''\
                                Flux Calculator
                    -------------------------------------
                         Calculates continuum source
                            flux for calibration
                        ''' )

        # Arguments list
        parser.add_argument( '-f', '--freq', dest = 'frequency', required = True, nargs = 1, help = "Freqency in GHz" )
        parser.add_argument( '-s', '--src', dest = 'source', required = True, nargs = 1, help = "Provide source as a string e.g. \"J1445+0958\"" )
        parser.add_argument( '--format1', dest = 'format1', action = 'store_true', default = False, help = "If true, looks for flux parameters in Format 1 (otherwise, Format 2)" )
        parser.add_argument( '-c', '--compare', dest = 'compare', nargs = '*', default = None, help = "Compares flux at given frequencies with frequency given in -f" )

        args = parser.parse_args()

        return args


    args = parser( 'calculate_flux.py' )

    flux = getFlux( args.frequency[0], args.source[0], args.format1 )
    print( "Flux (Jy): ", flux )

    if args.compare:
        print( "Freq (GHz)", ", ", "Flux (Jy)", ", ", "% error" )
        for arg in args.compare:
            comp_flux = getFlux( arg, args.source[0], args.format1 )
            err = 100 * ( abs( comp_flux - flux ) / flux )
            print( arg, " ", comp_flux, " ", err, "%" )
