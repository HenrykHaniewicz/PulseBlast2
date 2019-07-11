# Equate multi-files of different overall length together

# Imports
import os
import numpy as np
from astropy.io import fits

def get_like_files( dir ):

    like_files = []
    file_list = []

    for file in sorted( os.listdir( os.path.dirname( dir ) ) ):

        root, ext = os.path.splitext( file )
        obs_num = root[-4:]

        if float( obs_num ) == 0001:
            if len( file_list ) != 0:
                like_files.append( file_list )
            file_list = []

        file_list.append( file )


    return like_files

def equate_file_lists( list1, list2 ):

    """
    Matches up files with each other by finding a normalized midpoint and comparing.
    """

    l1, l2 = np.array( list1 ), np.array( list2 )
