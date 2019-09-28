"""
library containing [projected] halo profiles
"""

__author__ = "Siavash Yasini"
__email__ = "yasini@usc.edu"


import numpy as np

from astropy.constants import sigma_T
c = 299792. #km/s


#########################################################
#                Profiles
#########################################################

# ------------------------
#           3D
# ------------------------

def NFW_mass(r, rho_s, R_s):
    """
    Calculate the NFW profile #TODO: add reference Eq.

    Parameters
    ----------
    r:
        distance from the center
    rho_s:
        density at radius R_s
    R_s:
        characterisic radius R_200c/c_200c

    Returns
    -------
    rho = 4 * rho_s * R_s ** 3 / r / (r + R_s) ** 2
    """

    rho = 4 * rho_s * R_s ** 3 / r / (r + R_s) ** 2

    return rho
# ------------------------
#        Projected
# ------------------------


def projected_NFW_mass(r, rho_s, R_s):

    """
    projected NFW profile
    Eq. 7 in Bartlemann 1996: https://arxiv.org/abs/astro-ph/9602053

    Returns
    -------

    """
    x = np.asarray(r/R_s, dtype=np.complex)
    f = 1 - 2 / np.sqrt(1 - x ** 2) * np.arctanh(np.sqrt((1 - x) / (1 + x)))
    f = np.true_divide(f, x ** 2 - 1)

    return 8 * rho_s * R_s * f


def projected_solid_sphere(r, M_tot, R_tot):
    """
    projected mass density of uniform sphere

    Parameters
    ----------
    r: [Mpc]
        distance from the center
    M_tot: [M_sun]
        total mass of the sphere
    R_tot: [Mpc]
        total radius (edge) of the sphere

    Returns
    -------
    Sigma = M_tot /2/pi * sqrt(R_tot**2 - r**2)/R_tot**3

    """
    Sigma = M_tot / 2/np.pi * np.sqrt(R_tot ** 2 - r ** 2) / R_tot ** 3

    return Sigma

def constant_density(r, constant):

    """
    return a constant value at every input r

    Parameters
    ----------
    r: [Mpc]
        distance from the center
    constant:
        multiplicative constant

    Returns
    -------
    constant

    """

    return constant


def linear_density(r, intercept, slope):

    """
    return a r*constant at every input r

    Parameters
    ----------
    r: [Mpc]
        distance from the center
    intercept:
       intercept of the line
    slope:
        slope of the line

    Returns
    -------
    intercept

    """

    return intercept + r*slope


# tests:

def kSZ_T_solid_sphere(r, catalog_dataframe):


    M_200c = catalog_dataframe["M_200c"]
    R_200c = catalog_dataframe["R_200c"]
    v_r = catalog_dataframe["v_r"]

    Sigma = projected_solid_sphere(r, M_200c, R_200c)

    dT_T = Sigma * v_r
    return dT_T

def kSZ_T_NFW(r, rho_s, R_s, v_r):


    Sigma = projected_NFW_mass(r, rho_s, R_s)

    dT_T = -Sigma * v_r
    return dT_T
