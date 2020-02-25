import os
import time

import numpy as np

__author__ = "Siavash Yasini"
__email__ = "yasini@usc.edu"

from decorator import decorator

from scipy import integrate
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline

from .transform import arcmin2rad, fwhm2sigma
import yaml
from pprint import pprint
from contextlib import contextmanager

#########################################################
#           CMB and Noise Power Spectra
#########################################################

# ------------------
# CMB Power Spectrum
# ------------------

def load_Cl_Planck2018(lmin=0):
    """
    load Cl from camb generated Dl file

    Parameters
    ----------
    lmax: max ell number
    lmin: min ell number

    Returns
    -------
    Cls

    available keys in Cls : L, TT, EE, BB, TE
    """
    assert lmin == 0, "lmin=0 cannot be changed. It is only to indicate explicitly that the " \
                    "returned results will start from lmin=0.\n" \
                    "If you want to get the Cl in a custom ell range use utilities.get_CMB_Cl."

    Cl_fname = os.path.join(os.path.dirname(__file__),
                            "Cl_Planck2018_camb.npz")

    Cls = np.load(Cl_fname)
    # L = Cls['ell']
    # Cl_TT = Cls['tt']
    # Cl_EE = Cls['ee']
    # Cl_BB = Cls['bb']
    # Cl_TE = Cls['te']

    return Cls


def get_CMB_Cl(lmax, lmin=0, mode="TT", return_ell=False, uK=False):
    """
    load Cl from camb generated Dl file

    Parameters
    ----------
    lmax: int
        max ell number
    lmin: int
        min ell number
    mode: str
        CMB mode to return (e.g. "TT", "EE", etc)
    return_ell: bool
        if True, returns the corresponding ell array as well

    Returns
    -------
    Cl [K^2]
    or
    ell, Cl [K^2]

    available keys in Cls : L, TT, EE, BB, TE
    """

    Cl_fname = os.path.join(os.path.dirname(__file__),
                            "Cl_Planck2018_camb.npz")

    Cls = load_Cl_Planck2018()

    L = Cls['L'][lmin:lmax+1]
    Cl = Cls[mode][lmin:lmax+1]

    if uK:
        Cl *= 1E12
    if return_ell:
        return L, Cl
    else:
        return Cl


# --------------------
# Noise Power Spectrum
# --------------------

def _compute_Nl(sigma_n, lmax, lmin=0, fwhm=None, apply_beam=False, uK=False, return_ell=False):
    """
    compute the instrumental noise power spectrum (uK^2)

    Parameters
    ----------
    sigma_n [uK-arcmin]:
        noise level in uK-arcmin
        can a scalar or an array for multiple channels
        the length must match that of fwhm
    lmax:
        maximum ell mode in the power spectrum
    lmin:
        minimum ell mode in the power spectrum
    fwhm [arcmin]:
        beam fwhm in arcmins
        can be scalar or an array for multiple channels
    apply_beam:
        if True, deconvolves the noise with beam
    uK: bool
        if True, the returned Nl will be in uK^2 units
    return_ell: bool
        if True, returns the corresponding ell array as well

    Returns
    -------
    Nl [K^2]
    or
    ell, Nl [K^2]
    """

    # make sure input is an array
    if np.isscalar(sigma_n):
        sigma_n = [sigma_n]

    #convert sigma_n to radians
    sigma_n = arcmin2rad(np.asarray(sigma_n))

    # set up ell
    L = np.arange(lmin, lmax + 1)
    # determine number of channels
    n_channels = len(sigma_n)
    # calculate w^(-1) noise power spectrum prefactor
    w_inverse = sigma_n ** 2

    # calculate the noise power spectrum N_l (no beam) for each channel
    Nl_channel = np.array([w_inverse[channel] * np.ones(L.shape) for channel in range(n_channels)])

    # apply beam to the noise power spectrum
    if apply_beam:
        if fwhm:
            # convert scalar fwhm to list for consistency
            if np.isscalar(fwhm):
                fwhm = [fwhm]

            # convert fwhm to radians
            fwhm = arcmin2rad(np.asarray(fwhm))

            # check the length and dimensions of the input
            assert len(fwhm) == len(sigma_n), "fwhm and sigma_n must have the same length"
            assert np.ndim(fwhm) == np.ndim(sigma_n) == 1

            # convert fwhm to sigma
            sigma_theta = fwhm ** 2 / 8 / np.log(2)

            # calculate the beam power spectrum
            B2l_channel = np.array(
                [np.exp(-L * (L + 1) * sigma_theta[channel]) for channel in range(n_channels)])

            Nl_channel = np.array([Nl / Bl for Nl, Bl in zip(Nl_channel, B2l_channel)])

        else:
            raise ValueError("fwhm is not provided")

    if uK:
        Nl_channel *= 1E12

    if return_ell:
        return L, Nl_channel
    else:
        return Nl_channel


def get_experiment_Nl(lmax, lmin=0, name="Planck", frequency=[217], apply_beam=False, uK=False,
                      return_ell=False):
    """
    get temperature and polarization noise power spectra for various experiments

    Parameters
    ----------
    lmax: scalar
        maximum ell mode of the power spectrum
    lmin: scalar
        minimum ell mode of the power spectrum
    name:
        name of the experiment
        string in ["Planck" , "SO" , "S4"]
    apply_beam:
        if True, deconvolves the noise with beam
    return_ell: bool
        if True, returns the corresponding ell array as well

    Returns
    -------
    Nl_TT [K^2]
    or
    ell, Nl_TT [K^2]
    """

    # make sure input is an array
    if np.isscalar(frequency):
        frequency = [frequency]

    # load the noise configuration yaml file
    noise_yml = load_noise_yaml()

    try:
        # extract noise configuration (fwhm and sigma_n)
        # for each frequency channel from noise dictionary
        fwhm = [noise_yml[name][f]["fwhm"] for f in frequency]
        sigma_T = [noise_yml[name][f]["sigma_T"] for f in frequency]
        # sigma_P = [noise_yml[name][f]["sigma_P"] for f in frequency]
    except KeyError:
        print(
            "\n\nKey not found. Make sure both the experiment's 'name' and 'frequency' exist in "
            "the noise dictionary.\n")
        pprint(noise_yml)
        raise

    print(f"{name} noise @ {frequency} GHz:\n"
          f"fwhm    [arcmin]    = {fwhm}\n"
          f"sigma_T [uK-arcmin] = {sigma_T}\n"
          # f"sigma_P [uK-arcmin] = {sigma_P}"
          )

    # calculate the noise power spectrum
    Nl_TT = _compute_Nl(sigma_n=sigma_T, lmax=lmax, lmin=lmin, fwhm=fwhm, apply_beam=apply_beam,
                        uK=False)
    # Nl_EE = compute_Nl(sigma_n=sigma_P, lmax=lmax, lmin=lmin, fwhm=fwhm, apply_beam=apply_beam)

    # combine the frequency channels
    Nl_TT = combine_Nl(Nl_TT)
    # Nl_EE = combine_Nl(Nl_EE)
    L = np.arange(lmin, lmax+1)

    if not uK:
        Nl_TT *= 1E-12

    if return_ell:
        return L, Nl_TT
    else:
        return Nl_TT


def get_custom_Nl(sigma_n, lmax, fwhm=None, frequency=[217], lmin=0, apply_beam=False, uK=False,
                  return_ell=False):
    """
    get temperature and polarization noise power spectra for a custom experiment

    Parameters
    ----------
    sigma_n [uK-arcmin]:
        noise level in uK-arcmin
        can a scalar or an array for multiple channels
        the length must match that of fwhm
    lmax: scalar
        maximum ell mode of the power spectrum
    fwhm [arcmin]:
        beam fwhm in arcmins
        can be scalar or an array for multiple channels
    lmin: scalar
        minimum ell mode of the power spectrum
    apply_beam:
        if True, deconvolves the noise with beam
    return_ell: bool
        if True, returns the corresponding ell array as well

    Returns
    -------
    Nl [K^2]
    or
    ell, Nl [K^2]
    """

    # make sure input is an array
    if np.isscalar(frequency):
        frequency = [frequency]
    if np.isscalar(fwhm):
        fwhm = [fwhm]
    if np.isscalar(sigma_n):
        sigma_n = [sigma_n]

    print(f"Custom noise @ {frequency} GHz:\n"
          f"fwhm    [arcmin]    = {fwhm}\n"
          f"sigma_T [uK-arcmin] = {sigma_n}\n"
          # f"sigma_P [uK-arcmin] = {sigma_P}"
          )

    # calculate the noise power spectrum in K^2
    Nl_TT = _compute_Nl(sigma_n=sigma_n, lmax=lmax, lmin=lmin, fwhm=fwhm, apply_beam=apply_beam,
                        uK=False)
    # Nl_EE = compute_Nl(sigma_n=sigma_P, lmax=lmax, lmin=lmin, fwhm=fwhm, apply_beam=apply_beam)

    # combine the frequency channels
    Nl_TT = combine_Nl(Nl_TT)
    # Nl_EE = combine_Nl(Nl_EE)
    L = np.arange(lmin, lmax + 1)

    if not uK:
        Nl_TT *= 1E-12

    if return_ell:
        return L, Nl_TT
    else:
        return Nl_TT


def combine_Nl(Nls):
    """
    combine the input noise power spectra

    Parameters
    ----------
    Nls: list
        list of noise power spectra

    Returns
    -------
    combined noise power spectrum
    1/N_tot = 1/N1 + 1/N2 + ...
    """

    assert isinstance(Nls, (list, np.ndarray))
    assert np.ndim(Nls) == 2, "Input must be a list of 1D noise power spectra"

    N_tot_inv = np.array([np.true_divide(1, Nl) for Nl in Nls])

    return 1 / np.sum(N_tot_inv, axis=0)

def load_noise_yaml():
    # load the noise configuration yaml file
    lib_dir = os.path.dirname(os.path.abspath(__file__))
    noise_file = os.path.join(lib_dir, "noise.yml")
    with open(noise_file, "r") as file:
        noise_yml = yaml.load(file, Loader=yaml.FullLoader)

    return noise_yml


# --------------------
# Beam Power Spectrum
# --------------------
def get_custom_B2l(fwhm, lmax, lmin=0, arcmin=True, return_ell=False):
    """
    Compute the instrumental Beam power spectrum

    After smoothing the map with a beam of size fwhm, the power spectrum would be suppressed by a
    factor

    B2l= np.exp(-ell * (ell + 1) * sigma_b)

    where sigma_b = fwhm ** 2 / 8 / np.log(2)

    Parameters
    ----------
    fwhm [arcmin]:
        beam fwhm in arcmins (or radians if arcmin=False)
    lmax:
        maximum ell mode in the power spectrum
    lmin:
        minimum ell mode in the power spectrum
    arcmin: bool
        set to True if fwhm is in arcmin
    return_ell: bool
        if True, returns the corresponding ell array as well

    Returns
    -------
    Bl^2
    or
    ell, Bl^2
    """


    # set up ell
    L = np.arange(lmin, lmax + 1)

    sigma_b = fwhm2sigma(fwhm, arcmin=arcmin)


    # calculate the beam power spectrum
    B2l = np.exp(-L * (L + 1) * sigma_b)

    if return_ell:
        return L, B2l
    else:
        return B2l


#########################################################
#       sampling, projection, and interpolation
#########################################################

def sample_array(array, n_samples, method="linspace", eps=0.001):
    """sample an input array"""

    assert method in ["linspace", "logspace", "random"]

    array_min = array.min()
    array_max = array.max()

    if method == "linspace":
        samples = np.linspace(array_min, array_max, n_samples)
    elif method == "logspace":
        #TODO: fix the min max range issue
        if array_min < eps:
            array_min = eps
            print(array_min)
        samples = np.logspace(np.log10(array_min), np.log10(array_max), n_samples)
    elif method == "random":
        samples = np.random.uniform(array_min, array_max, size=n_samples)

    return samples


@decorator
def LOS_integrate(profile_3D, *args, **kwargs):
    """integrate along the line of sight for all 3D r's that correspond to the 2D R"""

    # extract R
    R = args[0]
    args = args[1:]
    #TODO: Add support for R_vec too

    # see if R is a scalar or a list (array)
    R_is_scalar = False
    if not hasattr(R, "__iter__"):
        R_is_scalar = True
        R = [R]

    # integrate along LOS for each R
    LOS_integrated = []
    for R_i in R:
        #TODO: Take f outside and profile the funtion
        f = lambda r: profile_3D(r, *args, **kwargs) * 2. * r / np.sqrt(r ** 2 - R_i ** 2)
        LOS_integrated.append(integrate.quad(f, R_i, np.inf, epsabs=0., epsrel=1.e-2)[0])

    # if R was a scalar, return a scalar
    if R_is_scalar:
        LOS_integrated = LOS_integrated[0]

    return LOS_integrated


@decorator
def interpolate(profile,
                n_samples=20,
                min_frac=None,
                sampling_method="linspace",
                k=3,
                interpolator=InterpolatedUnivariateSpline,
                *args, **kwargs):
    """
    interpolate the profile function instead of calculating it at every given R

    Parameters
    ----------
    profile:
        wrapped profile to be interpolated (implicit)

    n_samples: int
        number of points to sample in R

    min_frac: float
        fraction of points in R to sample, unless n_sample is larger
        if min_frac: n_samples = max(n_samples, min_frac * len(R))

        e.g. for n_sample=10, min_frac=0.1 if len(R)=200, 20 (0.1*200) points will be sampled,
        but if len(R)=50 10 points will be sampled

    sampling_method: str in ["linspace", "logspace", "random"]
        determines how the points are sampled

    k: int
        interpolation order

    interpolator: func
        interpolating function

    Returns
    -------
    Interpolated profile
    """

    assert n_samples > 1, "number of samples must be larger than 1"


    # extract R
    R = args[0]
    args = args[1:]
    # TODO: Add support for R_vec too

    if min_frac:
        assert 0 <= min_frac <= 1, "min_frac must be between 0 to 1"
        n_samples = max(n_samples, min_frac * len(R))

    print(n_samples)
    # if the input R vector is small, just calculate the profile directly
    if len(R) < n_samples:
        return profile(R, *args, **kwargs)

    else:
        # sample the input R and evaluate profile at those points
        R_samples = sample_array(R, n_samples, method=sampling_method)
        sample_values = np.array([profile(R_samp, *args, **kwargs) for R_samp in R_samples])

        # initialize the scipy interpolator
        profile_interp = interpolator(R_samples, sample_values, k=k)
        #print(k)
        return profile_interp(R)


#########################################################
#                   custom timer
#########################################################

@contextmanager
def timeit(process_name="Process"):
    """Time the code in mins"""

    time_stamp = time.strftime("%H:%M:%S %p")
    print("{:=>50}\n{} started at {}\n".format("", process_name, time_stamp))
    t_i = time.time()

    yield

    t_f = time.time()

    t = t_f-t_i

    print("{} was done in {:.1f} min.\n{:=>50}\n".format(process_name, t/60,""))


