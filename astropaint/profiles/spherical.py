"""
library containing simple [projected] spherical halo profiles
"""

__author__ = ["Siavash Yasini"]
__email__ = ["yasini@usc.edu"]

import numpy as np

from astropaint.lib.utils import interpolate, LOS_integrate
from astropaint.lib import transform



# -----------
# 3D profiles
# -----------

# Add profiles here!

# ---------------
# LOS projections
# ---------------

def solid_sphere_2D(R, M_200c, R_200c):
    """	
    projected mass density of uniform sphere	
    Parameters	
    ----------	
    r: [Mpc]	
        distance from the center	
    M_200c: [M_sun]	
        total mass of the sphere	
    R_200c: [Mpc]	
        total radius (edge) of the sphere	
    Returns	
    -------	
    Sigma = M_200c /2/pi * sqrt(R_tot**2 - r**2)/R_tot**3	
    """
    Sigma = M_200c / 2 / np.pi * np.sqrt(R_200c ** 2 - R ** 2) / R_200c ** 3

    return Sigma

def constant_density_2D(R, constant):

    """	
    return a constant value at every input r	
    Parameters	
    ----------	
    R: [Mpc]	
        distance from the center	
    constant:	
        multiplicative constant	
    Returns	
    -------	
    constant	
    """

    return constant


def linear_density_2D(R, intercept, slope):

    """	
    return a R*constant at every input R	
    Parameters	
    ----------	
    R: [Mpc]	
        distance from the center	
    intercept:	
       intercept of the line	
    slope:	
        slope of the line	
    Returns	
    -------	
    intercept	
    """

    return intercept + R * slope