# Custom exceptions file
# Add exceptions as necessary

class ArgumentError( Exception ):

    "Any place where incorrect or invalid arguments might exist should utilize this class."

    def __init__( self, message ):
        self.message = message

class DimensionError( Exception ):

    "Raise when the number of dimensions in an array is not desired."

    def __init__( self, message ):
        self.message = message

class TemplateLoadError( Exception ):

    ""

    def __init__( self, message ):
        self.message = message
