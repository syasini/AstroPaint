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

from astropaint.lib.utilities import interpolate, LOS_integrate

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

# M_200c is in Msun/h, redshift is redshift


def _rho0(M_200c, redshift):
    return 4.e3 * ((M_200c / h) / 1.e14) ** 0.29 * (1. + redshift) ** (-0.66)


def _alpha(M_200c, redshift):
    return 0.88 * ((M_200c / h) / 1.e14) ** (-0.03) * (1. + redshift) ** 0.19


def _beta(M_200c, redshift):
    return 3.83 * ((M_200c / h) / 1.e14) ** 0.04 * (1. + redshift) ** (-0.025)


def rho_gas_3D(r, R_200c, M_200c, redshift):
    """3d physical gas density profile
        in (Msun/h) / (Mpc/h)^3
        M_200c is in Msun/h
        R_200c is comoving radius in Mpc/h
        """

    # calculate the dimensionless 3d scaled gas density profile
    # fitting function from Battaglia 17
    x = r / R_200c
    fit = 1. + (x / xc) ** _alpha(M_200c, redshift)
    fit **= -(_beta(M_200c, redshift) + gamma) / _alpha(M_200c, redshift)
    fit *= (x / xc) ** gamma
    fit *= _rho0(M_200c, redshift)

    # rescale with comoving critical density, in (Msun/h) / (Mpc/h)^3
    rho_3d = fit * cosmo.critical_density(redshift).to(u.Msun/u.Mpc**3)

    # rescale with the baryon fraction
    rho_3d *= cosmo.Ob0 / cosmo.Om0

    return rho_3d


@interpolate(min_frac=0.1, sampling_method="logspace")
@LOS_integrate
def rho_gas_2D_interp(R, R_200c, M_200c, redshift):
    """2D scaled gas density profile
    dimensionless
    """
    return rho_gas_3D(R, R_200c, M_200c, redshift)


@LOS_integrate
def rho_gas_2D(R, R_200c, M_200c, redshift):
    """2D scaled gas density profile
    dimensionless
    """
    return rho_gas_3D(R, R_200c, M_200c, redshift)

def ne_3D(r, R_200c, M_200c, redshift):
    """3d physical electron number density profile
    in 1/(Mpc/h)^3
    assuming fully ionized gas and primordial He abundance
    m is mVir in Msun/h
    r is comoving radius in Mpc/h
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

    # get the physical 3d electron number density
    result /= factor  # in (Mpc/h)^(-3)

    return result


@interpolate
@LOS_integrate
def ne_2D(R, R_200c, M_200c, redshift):

    return ne_3D(R, R_200c, M_200c, redshift)


