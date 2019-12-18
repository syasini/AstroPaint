import os
import numpy as np

__author__ = "Siavash Yasini"
__email__ = "yasini@usc.edu"

def load_Planck2018_Cl(lmax, lmin=0,):
    """
    load Cl from camb generated Dl file

    Parameters
    ----------
    lmax: max ell number
    lmin: min ell number

    Returns
    -------
    L, Cl

    Cl = np.array([Cl_TT, Cl_EE, Cl_BB, Cl_TE])
    """

    Cl_fname = os.path.join(os.path.dirname(__file__),
                            "Planck2018_Cl_raw.npy")

    Cls = np.load(Cl_fname, allow_pickle=True, encoding="latin1").item()

    L = Cls['ell']
    Cl_TT = Cls['tt']
    Cl_EE = Cls['ee']
    Cl_BB = Cls['bb']
    Cl_TE = Cls['te']

    return L, Cl_TT, Cl_EE, Cl_BB, Cl_TE