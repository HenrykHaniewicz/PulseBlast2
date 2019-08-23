# Saving utilities

import os
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
