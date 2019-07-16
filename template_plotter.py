# Latex template plotter

import numpy as np
import matplotlib.pyplot as plt

clr = 'b'


def plot_templates( BW, *filenames ):

    fig = plt.figure( figsize = ( 8, 8 ) )
    plt.rc( 'text', usetex = True )
    plt.rc( 'font', family = 'serif' )

    if len( filenames ) == 0:
        plt.close()
        raise ValueError( "Need at least one template to load!" )
    elif len( filenames ) == 1:
        temp = np.load( filenames[0] )
        lin = np.linspace( 0, 1, temp.shape[-1] )
        if temp.ndim == 1:
            ax = fig.add_subplot(111)
            ax.plot( lin, temp, color = clr )
            ax.set_xlabel(r'Phase')
            ax.set_ylabel(r'Flux (Jy)')
            ax.set_xlim( 0, 1 )
        elif temp.ndim == 2:
            n_chan = temp.shape[0]
            c_sqrt = np.ceil( np.sqrt( n_chan ) )
            ax = []
            for i, chan in enumerate( temp ):
                ax.append( fig.add_subplot( c_sqrt, c_sqrt, i+1 ) )
            for j, axis in enumerate( ax ):
                axis.plot( lin, temp[j], color = clr )
                axis.ticklabel_format( style = 'sci', axis = 'y', scilimits = (-2, 2) )
                #axis.set_xlabel(r'Phase')
                axis.set_ylabel(r'Flux (Jy)')
                axis.set_xlim( 0, 1 )
                nu = (j * (BW/n_chan)) + 1030
                axis.set_title( r'{} MHz'.format( nu ) )
    elif len( filenames ) == 2:
        temps, ax = [], []
        for f in filenames:
            t = np.load( f )
            if t.ndim == 2:
                raise IndexError( "Can't plot more than one template in more than one frequency..." )
            temps.append( t )
        for i, tmp in enumerate( temps ):
            lin = np.linspace( 0, 1, len(tmp) )
            ax.append( fig.add_subplot( 2, 1, i+1 ) )
        for j, axis in enumerate( ax ):
            axis.plot( lin, temps[j], color = clr )
            axis.set_xlabel(r'Phase')
            axis.set_ylabel(r'Flux (Jy)')
            axis.set_xlim( 0, 1 )
            if j == 0:
                axis.set_title( r'L-band average profile' )
            elif j == 1:
                axis.set_title( r'430-band average profile' )

    else:
        raise ValueError( "Do again!" )


    fig.tight_layout()
    plt.show()
    return fig

a = ["/Users/zhn11tau/Documents/PulseBlastV2/J1829+2456/templates/J1829+2456_lbw_nchan1_template.npy", "/Users/zhn11tau/Documents/PulseBlastV2/J1829+2456/templates/J1829+2456_430_nchan1_template.npy"]
b = ["/Users/zhn11tau/Documents/PulseBlastV2/J1829+2456/templates/J1829+2456_lbw_nchan1_template.npy"]
c = ["/Users/zhn11tau/Documents/PulseBlastV2/J1829+2456/templates/J1829+2456_lbw_nchan8_template.npy"]

plot_templates( 800, *c )
