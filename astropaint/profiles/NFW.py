"""
library containing [projected] halo profiles
"""

__author__ = ["Siavash Yasini"]
__email__ = ["yasini@usc.edu"]

import numpy as np

from astropy import units as u
from astropy.constants import sigma_T, m_p
from astropy.cosmology import Planck18_arXiv_v2 as cosmo

from astropaint.lib.utilities import interpolate, LOS_integrate
from astropaint.lib import transform

# ---------------Caching----------------
# To cache templates use
#  the @memory.cache decorator
from joblib import Memory
cachedir = 'cache'
memory = Memory(cachedir, verbose=False)
# --------------------------------------

sigma_T = sigma_T.to(u.Mpc**2).value # [Mpc^2]
m_p = m_p.to(u.M_sun).value # [M_sun]
f_b = cosmo.Ob0/cosmo.Om0
c = 299792. #km/s
h = cosmo.h
T_cmb = 2.7251
Gcm2 = 4.785E-20 # G/c^2 (Mpc/M_sun)


# -----------
# 3D profiles
# -----------

def rho_3D(r, rho_s, r_s):
    """
    Calculate the 3D NFW density profile #TODO: add reference Eq.

    Parameters
    ----------
    r:
        distance from the center
    rho_s:
        density at radius r_s
    r_s:
        characterisic radius R_200c/c_200c

    Returns
    -------
    rho = 4 * rho_s * r_s ** 3 / r / (r + r_s) ** 2
    """

    rho = 4 * rho_s * r_s ** 3 / r / (r + r_s) ** 2

    return rho


# ---------------
# LOS projections
# ---------------

def rho_2D_bartlemann(R, rho_s, R_s):
    """
    projected NFW mass profile
    Eq. 7 in Bartlemann 1996: https://arxiv.org/abs/astro-ph/9602053

    Returns
    -------
    surface mass density: [M_sun/Mpc^2]
    """

    x = np.asarray(R / R_s, dtype=np.complex)
    f = 1 - 2 / np.sqrt(1 - x ** 2) * np.arctanh(np.sqrt((1 - x) / (1 + x)))
    f = f.real
    f = np.true_divide(f, x ** 2 - 1)
    Sigma = 8 * rho_s * R_s * f
    return Sigma


@LOS_integrate
def rho_2D(R, rho_s, R_s):
    """
    3D NFW profile intgrated along the line of sight

    Returns
    -------
    surface mass density: [M_sun/Mpc^2]
    """

    return rho_3D(R, rho_s, R_s)


@interpolate(n_samples=20, sampling_method="logspace")
@LOS_integrate
def rho_2D_interp(R, rho_s, R_s):
    """
    3D NFW profile intgrated along a sampled number of line of sights
    and then interpolated

    Returns
    -------
    surface mass density: [M_sun/Mpc^2]
    """

    return rho_3D(R, rho_s, R_s)


def deflection_angle(R, c_200c, R_200c, M_200c, *, suppress=True, suppression_R=8):
    """
    calculate the deflection angle of a halo with NFW profile
    Using Eq 6 in Baxter et al 2015 (1412.7521)

    Parameters
    ----------
    R:
        distance from the center of halo [Mpc]
    c_200c:
        halo concentration parameter
    R_200c:
        halo 200c radius in [Mpc]
    M_200c:
        halo 200c mass of halo in M_sun


    Returns
    -------
        the deflection angle at distance R from the center of halo
    """

    A = M_200c * c_200c ** 2 / (np.log(1 + c_200c) - c_200c / (1 + c_200c)) / 4. / np.pi
    C = 16 * np.pi * Gcm2 * A / c_200c / R_200c

    R_s = R_200c / c_200c
    x = R / R_s
    x = x.astype(np.complex)

    f = np.true_divide(1, x) * (np.log(x / 2) + 2 / np.sqrt(1 - x ** 2) *
                                np.arctanh(np.sqrt(np.true_divide(1 - x, 1 + x))))

    alpha = C * f

    # suppress alpha at large radii
    if suppress:
        suppress_radius = suppression_R * R_200c
        alpha *= np.exp(-(R / suppress_radius) ** 3)

    return alpha.real


def tau_2D(R, rho_s, R_s):
    """
    projected NFW tau profile
    Eq. 7 in Battaglia 2016 :

    Returns
    -------
    tau: [NA]
    """
    X_H = 0.76
    x_e = (X_H + 1) / 2 * X_H
    f_s = 0.02
    mu = 4 / (2 * X_H + 1 + X_H * x_e)

    Sigma = rho_2D_bartlemann(R, rho_s, R_s)
    tau = sigma_T * x_e * X_H * (1 - f_s) * f_b * Sigma / mu / m_p
    return tau


def kSZ_T(R, rho_s, R_s, v_r, *, T_cmb=T_cmb):
    """kinetic Sunyaev Zeldovich effect
    #TODO: add reference"""
    tau = tau_2D(R, rho_s, R_s)
    dT = -tau * v_r / c * T_cmb

    return dT


def BG(R_vec, c_200c, R_200c, M_200c, theta, phi, v_th, v_ph, *, T_cmb=T_cmb):
    """Birkinshaw-Gull effect
    aka moving lens
    aka Rees-Sciama (moving gravitational potential)"""

    R = np.linalg.norm(R_vec, axis=-1)
    # R_hat = np.true_divide(R_vec, R[:, None])
    R_hat = np.true_divide(R_vec, np.expand_dims(R, axis=-1))

    alpha = deflection_angle(R, c_200c, R_200c, M_200c)
    v_vec = transform.convert_velocity_sph2cart(theta, phi, 0, v_th, v_ph)
    dT = -alpha * np.dot(R_hat, v_vec) / c * T_cmb
    return dT