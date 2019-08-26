# Saving utilities

import os
import pickle
import numpy as np
try:
    import astropy.io.fits as pyfits
except:
    import pyfits
from pypulse.archive import Archive

def save_psrfits( save_filename, ar ):

      primaryhdu = pyfits.PrimaryHDU(header=ar.header) #need to make alterations to header
      hdulist = pyfits.HDUList(primaryhdu)

      if ar.history is not None:
          cols = []
          for name in ar.history.namelist:
              fmt,unit,array = ar.history.dictionary[name]
              #print name,fmt,unit,array
              col = pyfits.Column(name=name,format=fmt,unit=unit,array=array)
              cols.append(col)
          historyhdr = pyfits.Header()
          for key in ar.history.headerlist:
              historyhdr[key] = ar.history.header[key]
          historyhdu = pyfits.BinTableHDU.from_columns(cols,name='HISTORY',header=historyhdr)
          hdulist.append(historyhdu)
          # Need to add in PyPulse changes into a HISTORY
      #else: #else start a HISTORY table


      if ar.params is not None:
          cols = [pyfits.Column(name='PSRPARAM',format='128A',array=ar.params.filename)] #PARAM and not PSRPARAM?
          paramhdr = pyfits.Header()
          for key in ar.paramheaderlist:
              paramhdr[key] = ar.paramheader[key]
          paramhdu = pyfits.BinTableHDU.from_columns(cols,name='PSRPARAM')
          hdulist.append(paramhdu)
          # Need to include mode for PSREPHEM




      if ar.polyco is not None:
          cols = []
          for name in ar.polyco.namelist:
              fmt,unit,array = ar.polyco.dictionary[name]
              #print name,fmt,unit,array
              col = pyfits.Column(name=name,format=fmt,unit=unit,array=array)
              cols.append(col)
          polycohdr = pyfits.Header()
          for key in ar.polyco.headerlist:
              polycohdr[key] = ar.polyco.header[key]
          polycohdu = pyfits.BinTableHDU.from_columns(cols,name='POLYCO',header=polycohdr)
          hdulist.append(polycohdu)




      if len(ar.tables) > 0:
          for table in ar.tables:
              hdulist.append(table)

      cols = []
      for name in ar.subintinfolist:
          fmt,unit,array = ar.subintinfo[name]
          col = pyfits.Column(name=name,format=fmt,unit=unit,array=array)
          cols.append(col)
          # finish writing out SUBINT!

      cols.append(pyfits.Column(name='DAT_FREQ',format='%iE'%np.shape(ar.freq)[1],unit='MHz',array=ar.freq)) #correct size? check units?
      cols.append(pyfits.Column(name='DAT_WTS',format='%iE'%np.shape(ar.weights)[1],array=ar.weights)) #call getWeights()

      nsubint = ar.getNsubint()
      npol = ar.getNpol()
      nchan = ar.getNchan()
      nbin = ar.getNbin()
      DAT_OFFS = ar.data_offset
      DAT_SCL = ar.data_scale
      saveDATA = ar.data_beforeanything

      cols.append(pyfits.Column(name='DAT_OFFS',format='%iE'%np.size(DAT_OFFS[0]),array=DAT_OFFS))
      cols.append(pyfits.Column(name='DAT_SCL',format='%iE'%np.size(DAT_SCL[0]),array=DAT_SCL))
      cols.append(pyfits.Column(name='DATA',format='%iI'%np.size(saveDATA[0]),array=saveDATA,unit='Jy',dim='(%s,%s,%s)'%(nbin,nchan,npol))) #replace the unit here

      subinthdr = pyfits.Header()
      for key in ar.subintheaderlist:
          subinthdr[key] = ar.subintheader[key]
      subinthdu = pyfits.BinTableHDU.from_columns(cols,name='SUBINT',header=subinthdr)
      hdulist.append(subinthdu)

      hdulist.writeto( save_filename, overwrite = True )


def load_session( bin_file, arg = None ):

    if len( arg ) is not 1:
        raise ValueError( "Please only try to load one session at a time..." )

    root, ext = os.path.splitext( bin_file )

    if ext is 'pkl':
        if arg == 'm':
            return load_session_template_pickle( bin_file )
        elif arg == 'r':
            return load_session_rfi_pickle( bin_file )
        elif arg == 'c':
            return load_session_calibration_pickle( bin_file )
        else:
            return 0
    elif ext is 'json':
        raise TypeError( "Need to figure out how JSON works..." )
    else:
        return 0

def load_session_template_pickle( bin_file ):

    try:
        pickle_in = open( bin_file, "rb" )
        temp_dict = pickle.load( pickle_in )
        pickle_in.close()
    except OSError:
        temp_dict = { 'DATA': None, 'IG_LIST': [] }
        pickle_out = open( bin_file, "wb" )
        pickle.dump( temp_dict, pickle_out )
        pickle_out.close()

    return temp_dict[ 'DATA' ], temp_dict[ 'IG_LIST' ], bin_file

def load_session_rfi_pickle( bin_file ):

    try:
        pickle_in = open( bin_file, "rb" )
        load_dict = pickle.load( pickle_in )
        pickle_in.close()
    except OSError:
        load_dict = { 'FILE' : None, 'DATA': None, 'POS' : [0, 0], 'IG_LIST': [] }
        pickle_out = open( bin_file, "wb" )
        pickle.dump( load_dict, pickle_out )
        pickle_out.close()

    return load_dict[ 'FILE' ], load_dict[ 'DATA' ], load_dict[ 'POS' ], load_dict[ 'IG_LIST' ]

def load_session_calibration_pickle( bin_file ):

    # Needs work...

    try:
        pickle_in = open( bin_file, "rb" )
        load_dict = pickle.load( pickle_in )
        pickle_in.close()
    except OSError:
        if self.verbose:
            print( "Could not retrieve fluxcal saved data." )
            print( "Attempting to create a new save file..." )
