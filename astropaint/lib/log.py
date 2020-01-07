"""
library for logging functions and classes
"""

__author__ = "Siavash Yasini"
__email__ = "yasini@usc.edu"



class ParameterNotFound(KeyError):
    pass

class CMBAlreadyAdded(Exception):
    pass

class NoiseAlreadyAdded(Exception):
    pass
