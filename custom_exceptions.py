# Custom exceptions file
# Add exceptions as necessary

class ArgumentError( Exception ):

    "Arguments on CLI do not match known arguments"

    def __init__( self, message ):
        self.message = message

class DimensionError( Exception ):

    "Mismatch of array dimensions"

    def __init__( self, message ):
        self.message = message

class TemplateLoadError( Exception ):

    "There was a problem loading the template profile"

    def __init__( self, message ):
        self.message = message
        # This is a thing
