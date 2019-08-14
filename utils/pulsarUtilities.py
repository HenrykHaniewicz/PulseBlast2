# Utilities for pulsar calculations (code by Robert Ferdman)

# Imports
import numpy as np
import math
import utils.otherUtilities as u
import utils.mathUtils as mu
from custom_exceptions import DimensionError

from pypulse.singlepulse import SinglePulse

# Functions

def get_data_from_asc( asc_file, duty = None ):
    x, y = u.get_data_from_asc( asc_file )
    if duty is not None:
        y = removeBase( y, duty )
    return x, y

def getBase( profData, duty ):
     # get profile mask to determine off-pulse bins
     mask = get_1D_OPW_mask( profData, windowsize = (len( profData ) - 100) )
     # select those with mask==0, i.e. baseline bins
     baseline = profData[mask == 0]
     # get mean and rms of baseline
     baseMean = np.mean( baseline )
     #baseRMS = np.std( baseline )
     baseRMS = mu.rootMeanSquare( baseline )

     # return tuple consisting of mean and rms of baseline
     return baseMean, baseRMS


# Returns the profile data minus the baseline
def removeBase( profData, duty ):

     baseline, baseRMS = getBase( profData, duty )

     # remove baseline mean from profile in place
     profData = profData - baseline

     return profData

def chan_to_freq( ctr_freq, bandwidth, nchan ):

    start_freq = ctr_freq - ( bandwidth / 2 )
    end_freq = ctr_freq + ( bandwidth / 2 )
    freq = np.linspace( start_freq, end_freq, num = nchan )
    return freq


def get_1D_OPW_mask( vector, **kwargs ):

    if vector.ndim is not 1:
        raise DimensionError( "Input data must be 1 dimensional to create an OPW mask." )

    mask = []
    sp_dat = SinglePulse( vector, **kwargs )

    while len( mask ) < len( vector ):
        if len( mask ) in sp_dat.opw:
            mask.append( False )
        else:
            mask.append( True )

    mask = np.asarray( mask )
    return mask


# Contour things

def loadContourArrays( fileprefix ):

    """
    Loads a numpy array of a set filename structure
    """

    load_array = np.load(fileprefix+'_params.npy')
    # for the purposes of this routine, only need the following
    # things in p_out
    p_out = {'m2':load_array[0],
             'mtot':load_array[1],
             'm1':load_array[2],
             'm1_prob':load_array[3]}
    p_out['norm_like'] = np.load(fileprefix+'_prob.npy')

    return p_out


def get_prob_2D_levels( z, prob_intervals, norm = False, n_steps = 64, weights = None ):

    """
    Calculates the confidence levels of a set of data
    """


# If weights are None, assign them to ones, with the same shape as the
# input z array:
    if weights is None:
        weights = np.ones_like( z, dtype = float )
# If weights are given, ensure they are the same shape as z:
    else:
        if weights.shape is not z.shape:
            print( 'Shape of weight array ', weights.shape, ' does not match the input data array ', z.shape )
            return None

# Do normalization
    if norm is True:
        z = (z*weights) / np.sum( z*weights )
# Find maximum value of normalized 2D probability function
    z_max = np.amax( z )
# Set up contour_levels array
    contour_level = np.zeros_like( prob_intervals ) # Ensures same dimenstions

# Initialize step size to half the maximum probability value, as well as
# intial pdf value
    step_size = z_max / 2
    z_level = z_max - step_size

# Now, for each of the given contour levels, determine z level that corresponds to
# an enclosed probability of that contour level.  Do this by stepping around the
# final value until we arrive at the step threshold.
    for i_prob in np.arange( len( prob_intervals ) ):
# Run through number of steps given by user, dividing in half each time
        for i_step in np.arange( n_steps ):
            test_ind = np.where( z >= z_level )
            test_prob = np.sum( z[test_ind]*weights[test_ind] )
            step_size = step_size / 2
            if test_prob > prob_intervals[i_prob]:
                z_level += step_size
            else:
                z_level -= step_size
# Now reset step_size to half the current prob interval (e.g. 0.683, 0.954, 0.9973)
        step_size = z_level / 2
# Now that we have gone down to desired step threshold, set current z_level to
# the level corresponding to desired current interval in loop
        contour_level[i_prob] = z_level


    return contour_level



def plot_contour_pdf( x_val, y_val, contour_data, n_steps = 64, linecolour = 'black', **kwargs ):

    """
    Plots confidence contour maps of one array compared to another
    """

    # Check for default kwargs
    u.check_kwarg( None, 'weights', 'canvassize', 'xlim', 'ylim', 'figtext', **kwargs )
    u.check_kwarg( True, 'xticks', 'yticks', 'xlabel', 'ylabel', **kwargs )
    u.check_kwarg( False, 'norm', 'hgrid', 'vgrid', **kwargs )
    u.check_kwarg( 16, 'figtextsize', **kwargs )
    u.check_kwarg( 18, 'ticklabelsize', 'axislabelsize', **kwargs )
    u.check_kwarg( 35, 'xstart', 'xend', 'ystart', 'yend', **kwargs )
    u.check_kwarg( 1, 'xscale', 'yscale', **kwargs )


    # Modify limits and scales for plots
    xlim = np.asarray(xlim)
    ylim = np.asarray(ylim)

    xlim_scale = xlim / xscale
    ylim_scale = ylim / yscale

    x_val_scale = x_val / xscale
    y_val_scale = y_val / yscale

# If weights are None, assign them to ones, with the same shape as the
# input z array:
    if weights is None:
        weights = np.ones_like( contour_data, dtype = float )
# If weights are given, ensure they are the same shape as z:
    else:
        if weights.shape is not contour_data.shape:
            print( 'Shape of weight array ', weights.shape, ' does not match the input data array ', contour_data.shape )
            return None

# Start by setting up lot limits.  Assuming 1-D array input:
    xmin = np.min( x_val )
    xmax = np.max( x_val )
    ymin = np.min( y_val )
    ymax = np.max( y_val )
    xspan = abs( xmax - xmin )
    yspan = abs( ymax - ymin )

# Set up the plot:
    fig = plt.figure( figsize = canvassize )
    ax = fig.add_axes( [xstart, ystart, xend, yend] )
    ax.xaxis.set_tick_params( labelsize = ticklabelsize, pad = 8 )
    ax.yaxis.set_tick_params( labelsize = ticklabelsize, pad = 8 )
    ax.ticklabel_format( axis = 'x', useOffset = False )
    ax.ticklabel_format( axis = 'y', useOffset = False )
    if xlim is None:
        ax.set_xlim( xmin - 0.01*xspan, xmax + 0.02*xspan )
    else:
        ax.set_xlim( xlim_scale )
    if ylim is None:
        ax.set_ylim( ymin - 0.01*yspan, ymax + 0.02*yspan )
    else:
        ax.set_ylim( ylim_scale )

    if xlabel is not None:
        ax.set_xlabel( xlabel, fontsize = axislabelsize, labelpad = 12 )
    if ylabel is not None:
        ax.set_ylabel( ylabel, fontsize = axislabelsize, labelpad = 12 )

    if not xticks:
        for tick in ax.xaxis.get_major_ticks():
            tick.label1On = False
            tick.label2On = False

    if not yticks:
        for tick in ax.yaxis.get_major_ticks():
            tick.label1On = False
            tick.label2On = False

    if hgrid:
        ax.yaxis.grid(linestyle = '--', color = 'black', linewidth = 0.4 )
    if vgrid:
        ax.xaxis.grid( linestyle = '--', color = 'black', linewidth = 0.4 )

    prob_intervals = np.array( [0.683, 0.954, 0.9973] )

# Create levels at which to plot contours at each of the above intervals.
# Will not assume going in that Z values are normalized to total volume of 1.
    contour_level = get_prob_2D_levels( contour_data, prob_intervals, n_steps = n_steps )
    print( "CONTOUR_LEVEL = ", contour_level )

    if norm is True:
        z_val = ( contour_data*weights ) / np.sum( contour_data*weights )
    else:
        z_val = contour_data

# Now plot the pdf data
    ax.contour( x_val_scale, y_val_scale, z_val, levels = np.flip( contour_level, axis = 0 ), colors = ( 'red', 'blue', 'green' ) )

    if figtext is not None:
        for txt in figtext:
            ax.text( txt[0], txt[1], txt[2], fontsize = figtextsize, horizontalalignment = 'center', verticalalignment = 'center' )
