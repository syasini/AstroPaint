"""
library containing fun profiles!
"""

__author__ = ["Siavash Yasini"]
__email__ = ["yasini@usc.edu"]

import numpy as np

from astropaint.lib.utils import interpolate, LOS_integrate
from astropaint.lib import transform

# -----------
# 3D profiles
# -----------

# Add profiles here!

# ---------------
# LOS projections
# ---------------

def drops(R, R_200c, x, y, z):
    """profile to emulate water droplets...
    use with cm.Blues colormap"""

    droplet = np.sin(R / R_200c * np.random.rand() * 10)
    droplet *= np.exp(-(R / R_200c / 2) ** 2) * 2 * (x + y + z)
    return droplet


def bacteria(R, R_200c, M_200c, x):
    """Are these microbes?
    use with cm.Greys_r colormap
    """

    microbe = np.log(M_200c) * np.exp(-(2 * (R - 1) / R_200c) ** 2) * x
    microbe -= np.log(M_200c) ** 2 * np.exp(-(4 * R / R_200c) ** 2)
    return microbe


def twilight(R, R_200c, x, y, z):
    """nobody knows where this profile came from...
    use with cm.twilight profile
    """

    layer = np.exp(-(R / R_200c) / 4) * (x + y + z)
    layer += np.exp(-(R / R_200c * 1.5) ** 2) * (x - y - z)

    return layer



