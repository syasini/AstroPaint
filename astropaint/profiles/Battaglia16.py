"""
library containing [projected] gas profiles

E. Schaan:
Tau profile, from Battaglia 2016
watch typos in paper. This code is correct.
"""

__author__ = ["Siavash Yasini", "Emmanuel Schaan"]
__email__ = ["yasini@usc.edu", "eschaan@lbl.gov"]

import numpy as np

from astropy import units as u
from astropy.constants import sigma_T, m_p
from astropy.cosmology import Planck18_arXiv_v2 as cosmo

from astropaint.lib.utils import interpolate, LOS_integrate

# ---------------Caching----------------
# To cache templates use
#  the @memory.cache decorator
from joblib import Memory
cachedir = 'cache'
memory = Memory(cachedir, verbose=False)
# --------------------------------------

# ---------
# constants
# ---------

sigma_T = sigma_T.to(u.Mpc**2).value # [Mpc^2]
m_p = m_p.to(u.M_sun).value # [M_sun]
f_b = cosmo.Ob0/cosmo.Om0
c = 299792. #km/s
h = cosmo.h
T_cmb = 2.7251
Gcm2 = 4.785E-20 #(Mpc/M_sun)

mMin = 0.
mMax = np.inf
trunc = 2

# parameters from Battaglia 17
xc = 0.5
gamma = -0.2


def _rho0(M_200c, redshift):
    return 4.e3 * ((M_200c / h) / 1.e14) ** 0.29 * (1. + redshift) ** (-0.66)


def _alpha(M_200c, redshift):
    return 0.88 * ((M_200c / h) / 1.e14) ** (-0.03) * (1. + redshift) ** 0.19


def _beta(M_200c, redshift):
    return 3.83 * ((M_200c / h) / 1.e14) ** 0.04 * (1. + redshift) ** (-0.025)


# -----------
# 3D profiles
# -----------

def rho_gas_3D(r, R_200c, M_200c, redshift):
    """3D physical gas density profile in (Msun/h) / (Mpc/h)^3

    Parameters
    ----------
    r: arraylike, float
        3D distance from the center of the halo
    R_200c: float [Mpc/h]
        comoving radius of the halo
    M_200c: float [Msun/h]
        mass of the halo
    redshift:
        redshift of the halo
    """

    # calculate the dimensionless 3D scaled gas density profile
    # fitting function from Battaglia 17
    x = r / R_200c
    fit = 1. + (x / xc) ** _alpha(M_200c, redshift)
    fit **= -(_beta(M_200c, redshift) + gamma) / _alpha(M_200c, redshift)
    fit *= (x / xc) ** gamma
    fit *= _rho0(M_200c, redshift)

    # rescale with comoving critical density, in (Msun/h) / (Mpc/h)^3
    rho_3d = fit * cosmo.critical_density(redshift).to(u.Msun/u.Mpc**3).value

    # rescale with the baryon fraction
    rho_3d *= cosmo.Ob0 / cosmo.Om0

    return rho_3d


def ne_3D(r, R_200c, M_200c, redshift):
    """3D physical electron number density profile in 1/(Mpc/h)^3
    assuming fully ionized gas and primordial He abundance

    Parameters
    ----------
    r: arraylike, float
        3D distance from the center of the halo
    R_200c: float [Mpc/h]
        comoving radius of the halo
    M_200c: float [Msun/h]
        mass of the halo
    redshift:
        redshift of the halo
    """
    result = rho_gas_3D(r, R_200c, M_200c, redshift)  # (Msun/h) / (Mpc/h)^3

    # convert from baryon mass to electron number
    me = 9.10938291e-31  # electron mass (kg)
    mH = 1.67262178e-27  # proton mass (kg)
    mHe = 4. * mH  # helium nucleus mass (kg)
    xH = 0.76  # primordial hydrogen fraction by mass
    nH_ne = 2. * xH / (xH + 1.)
    nHe_ne = (1. - xH) / (2. * (1. + xH))
    factor = (me + nH_ne * mH + nHe_ne * mHe)  # in kg
    msun = 1.989e30  # solar mass (kg)
    factor /= msun  # in Msun
    factor *= h  # in Msun/h

    # get the physical 3D electron number density
    result /= factor  # in (Mpc/h)^(-3)

    return result


def tau_3D(r, R_200c, M_200c, redshift):
    """Thompson scattering optical depth 3D profile in 1/(Mpc/h) comoving
    ie you get the 2D tau profile by projecting this profile
    along the physical (not comoving) radial coordinate
    assuming fully ionized gas and primordial He abundance

    Parameters
    ----------
    r: arraylike, float
        3D distance from the center of the halo
    R_200c: float [Mpc/h]
        comoving radius of the halo
    M_200c: float [Msun/h]
        mass of the halo
    redshift:
        redshift of the halo
    """

    ne3d = ne_3D(r, R_200c, M_200c, redshift)

    # multiply by Thompson cross section (physical)
    tau = sigma_T * ne3d

    return tau


# ---------------
# LOS projections
# ---------------

@interpolate(min_frac=0.1, sampling_method="logspace")
@LOS_integrate
def rho_gas_2D_interp(R, R_200c, M_200c, redshift):
    """Interpolated version of rho_gas_2D
    2D physical gas density profile in (Msun/h) / (Mpc/h)^3

    Parameters
    ----------
    R: arraylike, float
        2D projected distance from the center of the halo
    R_200c: float [Mpc/h]
        comoving radius of the halo
    M_200c: float [Msun/h]
        mass of the halo
    redshift:
        redshift of the halo
    """

    return rho_gas_3D(R, R_200c, M_200c, redshift)


@LOS_integrate
def rho_gas_2D(R, R_200c, M_200c, redshift):
    """2D physical gas density profile in (Msun/h) / (Mpc/h)^3

    Parameters
    ----------
    R: arraylike, float
        2D projected distance from the center of the halo
    R_200c: float [Mpc/h]
        comoving radius of the halo
    M_200c: float [Msun/h]
        mass of the halo
    redshift:
        redshift of the halo
    """
    return rho_gas_3D(R, R_200c, M_200c, redshift)


@interpolate
@LOS_integrate
def ne_2D(R, R_200c, M_200c, redshift):
    """2D physical electron number density profile in 1/(Mpc/h)^3
        assuming fully ionized gas and primordial He abundance

    Parameters
    ----------
    R: arraylike, float
        2D projected distance from the center of the halo
    R_200c: float [Mpc/h]
        comoving radius of the halo
    M_200c: float [Msun/h]
        mass of the halo
    redshift:
        redshift of the halo
        """
    return ne_3D(R, R_200c, M_200c, redshift)


@LOS_integrate
def tau_2D(R, R_200c, M_200c, redshift):
    """Thompson scattering optical depth 2D projected profile in 1/(Mpc/h) comoving
        ie you get the 2D tau profile by projecting this profile
        along the physical (not comoving) radial coordinate
        assuming fully ionized gas and primordial He abundance

    Parameters
    ----------
    R: arraylike, float
        2D projected distance from the center of the halo
    R_200c: float [Mpc/h]
        comoving radius of the halo
    M_200c: float [Msun/h]
        mass of the halo
    redshift:
        redshift of the halo
        """
    return tau_3D(R, R_200c, M_200c, redshift)


@interpolate(n_samples=15)
@LOS_integrate
def tau_2D_interp(R, R_200c, M_200c, redshift):
    """ Interpolated version of tau_2D:
    Thompson scattering optical depth 2D projected profile in 1/(Mpc/h) comoving
    ie you get the 2D tau profile by projecting this profile
    along the physical (not comoving) radial coordinate
    assuming fully ionized gas and primordial He abundance

    Parameters
    ----------
    R: arraylike, float
        2D projected distance from the center of the halo
    R_200c: float [Mpc/h]
        comoving radius of the halo
    M_200c: float [Msun/h]
        mass of the halo
    redshift:
        redshift of the halo
        """
    return tau_3D(R, R_200c, M_200c, redshift)


def kSZ_T(R, R_200c, M_200c, v_r, redshift, *, T_cmb=T_cmb):
    """kinetic Sunyaev Zeldovich effect
    #TODO: add reference"""
    tau = tau_2D_interp(R, R_200c, M_200c, redshift)
    dT = -tau * v_r / c * T_cmb * (1+redshift)

    return dT