#!/usr/bin/env python
# coding: utf-8

"""
library for transforming halo properties and coordinates
"""

__author__ = "Siavash Yasini"
__email__ = "yasini@usc.edu"

import numpy as np
from scipy.ndimage import rotate, zoom
import warnings
import logging
#logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()
logger.setLevel(logging.ERROR)


# ------------------------
#        constants
# ------------------------
from astropy.constants import sigma_T, m_p
from astropy.cosmology import z_at_value, Planck18_arXiv_v2 as cosmo
from astropy import units as u

T_0 = 2.725E6 #uK
c = 299792. #km/s
c2 = c**2 #(km/s)^2
#Gcm3= 1.597E-25 #(s/km) with masses in M_sol and lengths in Mpc
Gcm2 = 4.785E-20 #(Mpc/M_sun)

h = cosmo.H(0).value/100.
crit_dens_0 = cosmo.critical_density0.to(u.solMass/u.Mpc**3).value #M_sun/Mpc**3

sigma_T = sigma_T.to(u.Mpc**2).value # [Mpc^2]
m_p = m_p.to(u.M_sun).value # [M_sun]
#sigma_T_over_mp = sigma_T.to(value/m_p.to(u.M_sun).value # m^2
f_b = cosmo.Ob0/cosmo.Om0
c = 299792. #km/s
#########################################################
#                Halo Transformations
#########################################################

# --------------------
# Mass Transformations
# --------------------

def M_200c_to_R_200c(M_200c, redshift):
    """calculate R_200c from M_200c at the given redshift
    see Eq. 1 in Huang et al 1701.04001"""

    crit_dens = crit_dens_0 * (1 + redshift) ** 3
    R_200c = np.power((3*M_200c)/(800*np.pi*crit_dens), 1/3)

    return R_200c

def M_200c_to_c_200c(M_200c, redshift):
    """calculate the concentration parameter from M_200c at the given redshift
    use fitting formula in Eq 19 of Child et al 2018 (1804.10199)"""

    #TODO: Double check this
    #fitting parameters for stacked NFW (Table 2 of Child et al 2018)
    #A = 57.6
    #d = -0.376
    #m = -0.078

    #TODO: add a "model" argument

    #fitting parameters from Coe 2010 (1005.0411)
    A = 5.71
    d = -0.047
    m = -0.097
    M_0 = 2.E12/h #M_sun

    c_200c = A * (1 + redshift) ** d * (M_200c / M_0) ** m

    return c_200c


def M_200c_to_rho_s(M_200c, redshift, R_200c=None, c_200c=None):
    """calculate the NFW rho_s parameter from M_200c at the given redshift
    if R_200c and c_200c are not given, calculate them"""

    if R_200c is None:
        R_200c = M_200c_to_R_200c(M_200c, redshift)
    if c_200c is None:
        c_200c = M_200c_to_c_200c(M_200c, redshift)

    R_s = R_200c / c_200c
    A_c = np.log(1+c_200c) - c_200c/(1+c_200c)
    rho_s = M_200c / (16*np.pi*(R_s**3)*A_c)

    return rho_s

def M_to_tau(M):

    X_H = 0.76
    x_e = (X_H+1)/2*X_H
    f_s = 0.02
    mu = 4/(2*X_H+1+X_H*x_e)

    tau = sigma_T * x_e * X_H * (1-f_s) * f_b * M / mu / m_p

    return tau
# ------------------------
# Distance Transformations
# ------------------------

def D_c_to_D_a(D_c, redshift):
    """calculate the angular diameter distance (D_a) from comoving distance (D_c) and redshift (
    redshift)"""

    return D_c/(1 + redshift)

def D_c_to_redshift(D_c, units=u.Mpc):
    """calculate the redshift from comoving distance (D_c)"""

    return z_at_value(cosmo.comoving_distance, D_c * units)

def radius_to_angsize(radius, D_a, arcmin=True):
    """calculate  the angular radius (theta) of the halo from its radius and angular diameter
    distance (D_a).

    if arcmin == True: return the value in arcmin

    *NOTE: radius and D_a must have the same units"""

    ang_size = np.true_divide(radius, D_a)

    if arcmin: ang_size = rad2arcmin(ang_size)

    return ang_size



#########################################################
#           Coordinate Transformations
#########################################################

#TODO: add thph2lonlat conversion and vice versa

def rad2arcmin(angle):
    """convert radians to arcmins"""
    return np.rad2deg(angle)*60


def arcmin2rad(angle):
    """convert arcmins to radians"""
    return np.deg2rad(angle/60)


def get_cart2sph_jacobian(th, ph):
    """calculate the transformation matrix (jacobian) for spherical to cartesian coordinates at
    line of sight (th, ph) [radians]
        see https://en.wikipedia.org/wiki/Vector_fields_in_cylindrical_and_spherical_coordinates

        th is the polar angle with respect to z and ph is the azimuthal angle with respect to x

    example:

    th_rad = np.deg2rad(df['th'].values)
    ph_rad = np.deg2rad(df['ph'].values)

    v_cart = np.array([df['vx'],df['vy'],df['vz']])


    thph_grid = np.array([th_rad,ph_rad])

    J_cart2sph = cart2sph2(th_rad,ph_rad)
    v_cart2sph = np.einsum('ij...,i...->j...',J_cart2sph,v_cart)

    """

    row1 = np.stack((np.sin(th) * np.cos(ph), np.cos(th) * np.cos(ph), -np.sin(ph)))
    row2 = np.stack((np.sin(th) * np.sin(ph), np.cos(th) * np.sin(ph), np.cos(ph)))
    row3 = np.stack((np.cos(th), -np.sin(th), 0.0 * np.cos(th)))

    return np.stack((row1, row2, row3))


def get_sph2cart_jacobian(th, ph):
    """calculate the transformation matrix (jacobian) for spherical to cartesian coordinates at
    line of sight (th, ph) [radians]
        see https://en.wikipedia.org/wiki/Vector_fields_in_cylindrical_and_spherical_coordinates

        th is the polar angle with respect to z and ph is the azimuthal angle with respect to x

    example: see cart2sph2"""

    row1 = np.stack((np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)))
    row2 = np.stack((np.cos(th) * np.cos(ph), np.cos(th) * np.sin(ph), -np.sin(th)))
    row3 = np.stack((-np.sin(ph), np.cos(ph), 0.0 * np.cos(th)))

    return np.stack((row1, row2, row3))


def convert_velocity_sph2cart(th, ph, v_r, v_th, v_ph):
    """
    Calculate the cartesian velocity components from the spherical ones

    Parameters
    ----------

    Returns
    -------

    """
    #th = np.deg2rad(cat['co-lat'].values)
    #ph = np.deg2rad(cat['lon'].values)

    J = get_sph2cart_jacobian(th, ph)

    vij_sph = np.array([ v_r, v_th, v_ph])
    vij_cart = np.einsum('ij...,i...->j...', J, vij_sph)
    logger.info(f"v_ij.shape = {vij_cart.shape}")

    #v_x = vij_cart[0, :]
    #v_y = vij_cart[1, :]
    #v_z = vij_cart[2, :]
    v_x = vij_cart[0]
    v_y = vij_cart[1]
    v_z = vij_cart[2]

    return v_x, v_y, v_z


def convert_velocity_cart2sph(th, ph, v_x, v_y, v_z):

    # find the cartesian to spherical coords transformation matrix
    J_cart2sph = get_cart2sph_jacobian(th, ph)
    #J_sph2cart = transform.sph2cart(cat['co-lat'].values,cat['lon'].values)

    # transform the velocity field and define v_r (radial), v_th (co-latitude), v_ph (longitude)
    v_cart = np.array([v_x, v_y, v_z])
    v_r, v_th, v_ph = np.einsum('ij...,i...->j...', J_cart2sph, v_cart)
    # TODO: add columns or v_lat and v_lon

    return v_r, v_th, v_ph


#########################################################
#          Patch Transformations (for stacking)
#########################################################

def rotate_patch(patch, angle):
    """Rotate input patch by angle [deg]"""
    rotated_patch = rotate(patch, angle, reshape=False)

    return rotated_patch