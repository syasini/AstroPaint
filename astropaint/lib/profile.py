"""
library containing [projected] halo profiles
"""

__author__ = "Siavash Yasini"
__email__ = "yasini@usc.edu"


import numpy as np
from . import transform
from astropy.constants import sigma_T, m_p
from astropy.cosmology import Planck15 as cosmo
sigma_T = sigma_T.value # [m^2]
m_p = m_p.value # [kg]
f_b = cosmo.Ob0/cosmo.Om0
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


def NFW_mass_density_proj(r, rho_s, R_s):

    """
    projected NFW mass profile
    Eq. 7 in Bartlemann 1996: https://arxiv.org/abs/astro-ph/9602053

    Returns
    -------
    surface mass density: [M_sun/Mpc^2]
    """

    #FIXME: remove this
    r[r < 0.2] = 0.2  # flatten the core

    x = np.asarray(r/R_s, dtype=np.complex)
    f = 1 - 2 / np.sqrt(1 - x ** 2) * np.arctanh(np.sqrt((1 - x) / (1 + x)))
    f = np.true_divide(f, x ** 2 - 1)
    Sigma = 8 * rho_s * R_s * f
    return Sigma

def NFW_tau_density_proj(r, rho_s, R_s):

    """
    projected NFW tau profile
    Eq. 7 in Battaglia 2016 :

    Returns
    -------
    tau: [NA]
    """
    X_H = 0.76
    x_e = (X_H+1)/2*X_H
    f_s = 0.02
    mu = 4/(2*X_H+1+X_H*x_e)

    Sigma = NFW_mass_density_proj(r, rho_s, R_s)
    tau = sigma_T * x_e * X_H * (1-f_s) * f_b * Sigma / mu / m_p
    return tau

def solid_sphere_proj(r, M_200c, R_200c):
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
    Sigma = M_200c / 2 / np.pi * np.sqrt(R_200c ** 2 - r ** 2) / R_200c ** 3

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

def kSZ_T_solid_sphere(r, M_200c, R_200c, v_r):

    Sigma = projected_solid_sphere(r, M_200c, R_200c)
    tau = transform.M_to_tau(Sigma)
    dT_over_T = -tau * v_r

    return dT_over_T

def kSZ_T_NFW(r, rho_s, R_s, v_r):

    tau = NFW_tau_density_proj(r, rho_s, R_s)
    dT_over_T = -tau * v_r

    return dT_over_T
