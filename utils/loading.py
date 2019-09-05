# Loading utilities

import os
import pickle
import numpy as np


def load_session( bin_file, mode = None ):

    if len( mode ) is not 1:
        raise ValueError( "Please only try to load one session at a time..." )

    root, ext = os.path.splitext( bin_file )

    if ext == '.pkl':
        if mode == 'm':
            return load_session_template_pickle( bin_file )
        elif mode == 'r':
            return load_session_rfi_pickle( bin_file )
        elif mode == 'c':
            return load_session_calibration_pickle( bin_file )
        else:
            return 0
    elif ext == '.json':
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

    return temp_dict[ 'DATA' ], temp_dict[ 'IG_LIST' ]

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
        load_dict = None

    return load_dict
