# Advanced Gaussian fitting

import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import scipy.optimize
from numpy import exp, pi, sqrt
from utils.pulsarUtilities import get_1D_OPW_mask, get_data_from_asc
from utils.mathUtils import multi_norm, norm
from custom_exceptions import DimensionError

def print_parameters( p0, p1, n ):
    for i in np.arange( 0, 3*n, 3 ):
        print( p0[i], " -> ", p1[i], "\t", p0[i + 1], " -> ", p1[i + 1], "\t", p0[i + 2], " -> ", p1[i + 2] )

def get_gaussians( p1, n ):
    ind_gaussians = []
    for i in np.arange( 0, 3*n, 3 ):
        ind_gaussians.append( [ p1[i], p1[i + 1], p1[i + 2] ] )
    ind_gaussians = np.array( ind_gaussians )
    return ind_gaussians

def save_all( abs_dir, file_prefix, n, ig, m ):
    np.save( os.path.join( abs_dir, file_prefix + "_{}gaussianfit_individualgaussians.npy".format( n ) ), ig )
    np.save( os.path.join( abs_dir, file_prefix + "_{}gaussianfit.npy".format( n ) ), m )

def get_best_gaussian_fit( x, y, m_gauss = 15, bp = 15, p_wid = 150, guess = [ 1024, 20, 6000 ], plot_chisq = True, save_dir = None, f_pre = "", verbose = True ):
    if not isinstance( guess, np.ndarray ):
        guess_shape = np.array( guess ).shape
    else:
        guess_shape = guess.shape

    n_gauss = 0
    params, c = [], []
    while ( len(c) < m_gauss ) and ( n_gauss < bp ):

        if len( guess_shape ) == 1:
            params.extend( guess )
        elif len( guess_shape ) == 2:
            try:
                params.extend( guess[ n_gauss ] )
            except IndexError:
                params.extend( guess[0] )
        else:
            raise DimensionError( "Initial guess parameters must be a Nx3 array. Current shape of array is: {}".format( guess_shape ) )

        n_gauss = len( params )//3

        try:
            fitted_params,_ = scipy.optimize.curve_fit( multi_norm, x, y, p0 = params )
            if ( n_gauss == bp ) and verbose:
                print( "Maximum number of tries reached ({})".format( n_gauss ) )
        except RuntimeError:
            if len( guess_shape ) == 1:
                fitted_params = np.append( fitted_params, guess )
            elif len( guess_shape ) == 2:
                try:
                    fitted_params = np.append( fitted_params, guess[ n_gauss ] )
                except IndexError:
                    fitted_params = np.append( fitted_params, guess[0] )
            else:
                raise DimensionError( "You definitely shouldn't be able to see this error message." )

            if verbose:
                print("No fit for {} gaussians".format( n_gauss ))
                if n_gauss == bp:
                    print( "Maximum number of tries reached ({})".format( n_gauss ) )
            continue

        m = multi_norm( x, *fitted_params )
        mask = get_1D_OPW_mask( m, windowsize = ( len(m) - p_wid )  )
        for i, elem in enumerate( m ):
            if mask[i] == False:
                m[i] = 0

        chi2, p = scipy.stats.chisquare( y[mask == 1], f_exp = m[mask == 1] )
        if verbose:
            print( "Chi-sq for {} gaussians: ".format( n_gauss ), chi2 )
        c.append( chi2 )

    if verbose:
        print_parameters( params, fitted_params, n_gauss )

    if plot_chisq:
        plt.plot( c[1:] )
        plt.show()
        plt.close()

    ind_gaussians = get_gaussians( fitted_params, n_gauss )

    if save_dir is not None:
        save_all( save_dir, f_pre, n_gauss, ind_gaussians, m )

    return m, c, ind_gaussians, mask


# TESTING
if __name__ == "__main__":

    asc = r'/Users/zhn11tau/Documents/Programs/Python/J1829+2456_lbw_nchan1_template.ascii'
    x, y = get_data_from_asc( asc, duty = 0.05 )
    m, c, ig, msk = get_best_gaussian_fit( x, y, m_gauss = 7, save_dir = "/Users/zhn11tau/Documents/Programs/Python", f_pre = "J1829+2456_lbw_nchan1_template" )

    #plt.plot( x, y, 'k' )
    for p in ig:
        plt.plot( x[ msk == 1 ], norm( x[ msk == 1 ], p[0], p[1], p[2] ) )
    plt.plot( x[ msk == 1 ], m[ msk == 1 ], 'k' )
    plt.show()
