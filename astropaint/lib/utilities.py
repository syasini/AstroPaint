import os
import numpy as np

__author__ = "Siavash Yasini"
__email__ = "yasini@usc.edu"

from .transform import arcmin2rad
import yaml
from pprint import pprint

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


# --------------------
# Noise Power Spectrum
# --------------------

def compute_Nl(sigma_n, lmax, lmin=0, fwhm=None, apply_beam=False):
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

    Returns
    -------
    N_ell [uK^2]
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


    return Nl_channel


def get_experiment_Nl(lmax, lmin=0, name="Planck", frequency=[217], apply_beam=False, uK=False):
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


    Returns
    -------
    Nl_TT [uK^2]
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
    Nl_TT = compute_Nl(sigma_n=sigma_T, lmax=lmax, lmin=lmin, fwhm=fwhm, apply_beam=apply_beam)
    # Nl_EE = compute_Nl(sigma_n=sigma_P, lmax=lmax, lmin=lmin, fwhm=fwhm, apply_beam=apply_beam)

    # combine the frequency channels
    Nl_TT = combine_Nl(Nl_TT)
    # Nl_EE = combine_Nl(Nl_EE)
    if uK==False:
        Nl_TT *= 1E-12
    return Nl_TT  # , Nl_EE

def get_custom_Nl(lmax, sigma_n, fwhm, frequency=[217], lmin=0, apply_beam=False, uK=False):
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


    Returns
    -------
    Nl_TT [uK^2]
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

    # calculate the noise power spectrum
    Nl_TT = compute_Nl(sigma_n=sigma_n, lmax=lmax, lmin=lmin, fwhm=fwhm, apply_beam=apply_beam)
    # Nl_EE = compute_Nl(sigma_n=sigma_P, lmax=lmax, lmin=lmin, fwhm=fwhm, apply_beam=apply_beam)

    # combine the frequency channels
    Nl_TT = combine_Nl(Nl_TT)
    # Nl_EE = combine_Nl(Nl_EE)
    if uK==False:
        Nl_TT *= 1E-12
    return Nl_TT  # , Nl_EE


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