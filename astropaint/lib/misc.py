import os
import numpy as np

__author__ = "Siavash Yasini"
__email__ = "yasini@usc.edu"

def load_Cl_Planck2018(lmax, lmin=0,):
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

    Cl_fname = os.path.join(os.path.dirname(__file__),
                            "Cl_Planck2018_camb.npz")

    Cls = np.load(Cl_fname)
    # L = Cls['ell']
    # Cl_TT = Cls['tt']
    # Cl_EE = Cls['ee']
    # Cl_BB = Cls['bb']
    # Cl_TE = Cls['te']

    return Cls