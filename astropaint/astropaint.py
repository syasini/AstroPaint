"""
library for simulating semi-analytic mock maps of CMB secondary anisotropies
"""

__author__ = "Siavash Yasini"
__email__ = "yasini@usc.edu"

import numpy as np
import pandas as pd
from matplotlib import cm
from warnings import warn
import inspect

try:
    import healpy as hp

except ModuleNotFoundError:
    warn("Healpy is not installed. You cannot use the full sky canvas without it.")

from astropy.coordinates import cartesian_to_spherical
from astropaint.lib import transform


#########################################################
#                  Halo Catalog Object
#########################################################

class Catalog:

    """halo catalog containing halo masses, locations, velocities, and redshifts

    Units
    -----

    x, y, z: [Mpc]
    v_x, v_y, v_z: [km/s]
    M_200c: [M_sun]
    """
    def __init__(self,
                 data=None,
                 redshift=0,
                 build_dataframe=False):

        #TODO: define attribute dictionary with __slots__

        self.redshift = redshift
        self._data = data

        # if no input is provided generate a random catalog
        if self.data is None:
            self.data = self.generate_random_halos()

        # TODO: check input type/columns/etc
        self.size = len(self.data)

        # build the complete data frame
        # e.g. angular distances, radii, etc.
        if build_dataframe:
            self.build_dataframe()

    # ------------------------
    #       properties
    # ------------------------

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, val):
        self._data = val
        print("Input data has been modified. Rebuild the dataframe using catalog.build_dataframe "
              "to "
              "update the other parameters as well.")

    # ------------------------
    #         methods
    # ------------------------

    def build_dataframe(self):

        #TODO: add units documentation to the catalog for reference


        print("Building the dataframe...\n")

        # calculate the comoving distance and angular position (theta and phi in radians)
        self.data['D_c'], self.data['lat'], self.data['lon'] = cartesian_to_spherical(
                                                                    self.data['x'].values,
                                                                    self.data['y'].values,
                                                                    self.data['z'].values)

        # theta = pi/2 - lat , phi = lon
        self.data['theta'] = np.pi / 2 - self.data['lat']
        self.data['phi'] = self.data['lon']

        # convert lonlat coords to deg
        self.data['lon'], self.data['lat'] = np.rad2deg((self.data['lon'], self.data['lat']))

        # calculate angular diameter distance, virial radius and angular size
        self.data['D_a'] = transform.D_c_to_D_a(self.data['D_c'], self.redshift)
        self.data['R_200c'] = transform.M_200c_to_R_200c(self.data['M_200c'], self.redshift)
        self.data['c_200c'] = transform.M_200c_to_c_200c(self.data['M_200c'], self.redshift)
        self.data['R_th_200c'] = transform.radius_to_angsize(self.data['R_200c'],
                                                             self.data['D_a'], arcmin=True)



        # find the cartesian to spherical coords transformation matrix
        J_cart2sph = transform.get_cart2sph_jacobian(self.data['theta'].values,
                                                     self.data['phi'].values)
        # J_sph2cart = transform.sph2cart(self.data['co-lat'].values,self.data['lon'].values)

        # transform the velocity field and define v_r (radial), v_th (co-latitude), v_ph (longitude)
        v_cart = np.array([self.data['v_x'], self.data['v_y'], self.data['v_z']])
        self.data['v_r'], self.data['v_th'], self.data['v_ph'] = np.einsum('ij...,i...->j...',
                                                                           J_cart2sph, v_cart)

        self.data['v_lat'] = -self.data['v_th']
        self.data['v_lon'] = self.data['v_ph']



    @staticmethod
    def _initialize_catalog(n_tot):
        """initialize an empty catalog with x, y, z, v_x, v_y, v_z, M_200c columns"""

        dtype = {"names": ["x", "y", "z", "v_x", "v_y", "v_z", "M_200c"],
                 "formats": 7 * [np.float32]}

        catalog = np.zeros(n_tot, dtype)
        return catalog

    @staticmethod
    def generate_random_halos(box_size=100,
                              v_max=100,
                              mass_min=1E10,
                              mass_max=1E15,
                              n_tot=100,
                              put_on_shell=True):

        catalog = Catalog._initialize_catalog(n_tot)

        print("generating random catalog...\n")
        # generate random positions
        x, y, z = np.random.uniform(low=-box_size/2,
                                    high=box_size/2,
                                    size=(3, n_tot))

        if put_on_shell:
            (x, y, z) = box_size * np.true_divide((x, y, z), np.linalg.norm((x, y, z), axis=0))

        catalog["x"], catalog["y"], catalog["z"] = x, y, z

        # generate random velocities
        v_x, v_y, v_z = np.random.uniform(low=-v_max,
                                           high=v_max,
                                           size=(3, n_tot))

        catalog["v_x"], catalog["v_y"], catalog["v_z"] = v_x, v_y, v_z

        # generate random log uniform masses
        catalog["M_200c"] = np.exp(np.random.uniform(low=np.log(mass_min),
                                                     high=np.log(mass_max),
                                                     size=n_tot))

        return pd.DataFrame(catalog)  # convert catalog to pandas data frame



#########################################################
#                  Canvas Object
#########################################################

class Canvas:
    """healpy or flat-sky canvas with the location of the halos to paint the signal on"""

    def __init__(self,
                 catalog,
                 nside,
                 mode="healpy"):

        #TODO: define attribute dictionary with __slots__

        assert mode == "healpy", "currently only full sky is supported"

        self._nside = nside
        self._npix = hp.nside2npix(self.nside)
        self._cmap = cm.Greys_r

        self.pixels = np.zeros(self.npix)



        self.catalog = catalog
        self.centers_D_a = self.catalog.data.D_a


        assert isinstance(catalog, Catalog), "input catalog has to be an instance of " \
                                             "astroPaint.Catalog"

        self._proj_dict = {"mollweide": hp.mollview,
                           "mollview": hp.mollview,
                           "cartesian": hp.cartview,
                           "cartview": hp.cartview,
                           }

    # ------------------------
    #       properties
    # ------------------------

    # Immutables:

    @property
    def nside(self):
        return self._nside

    @property
    def npix(self):
        return self._npix

    # Mutables:

    @property
    def cmap(self):
        return self._cmap

    @cmap.setter
    def cmap(self, val):
        assert type(val) is type(cm.Greys), "cmap must be an instance of cm. \n" \
                                            "You can import it using:\n" \
                                            "from matplotlib import cm"
        self._cmap = val
        self._cmap.set_under("white")

    # ------------------------
    #         methods
    # ------------------------

    def clean(self):
        """
        Clean the canvas and set all pixels to zero

        Returns
        -------
        None
        """

        self.pixels = np.zeros(self.npix)

    def find_centers_indx(self):
        """
        Find the pixel indices of the halo centers

        Returns
        -------
        None
        Sets Canvas.centers_indx to array of pixels.
        Element [i] of the array points to the center of halo [i].

        """

        self.centers_indx = hp.ang2pix(self.nside,
                                       self.catalog.data.theta.to_list(),
                                       self.catalog.data.phi.to_list())

        print("Done! You can now get the center pixels using Canvas.centers_indx.")

    def find_discs_indx(self, k):
        """
        Find the pixel indices of discs of size k times R_200 around halo centers

        Parameters
        ----------
        k: int
            multiplicative factor indicating the extent of the queried disc in units of R_200

        Returns
        -------
        None

        Sets Canvas.discs_indx to a list of pixel index arrays. Element [i] of the list holds the
        pixel indices around halo [i].

        """

        #FIXME: list comprehension
        self.discs_indx = ([np.asarray(
                                    hp.query_disc(self.nside,
                                                  (self.catalog.data.x[halo],
                                    self.catalog.data.y[halo],
                                    self.catalog.data.z[halo]),
                                                  k * transform.arcmin2rad(self.catalog.data.R_th_200c[halo]))
                                       )
                                for halo in range(self.catalog.size)])

        print("Done! You can now get the discs using Canvas.discs_indx.")

    def find_centers_ang(self):
        """
        Store the theta and phi coordinates of the halos in Canvas.centers_ang

        Returns
        -------
        None
        """

        self.centers_ang = np.asarray([self.catalog.data.theta.to_list(),
                                     self.catalog.data.phi.to_list()])

        print("Done! You can now get the angular position of the discs using Canvas.centers_ang.")


    def find_discs_ang(self):
        """
        Find the angular coordinates of the disc pixels

        Returns
        -------
        None
        """
        try:
            self.discs_indx
        except AttributeError:
            print("Canvas.discs_indx is not defined. Use Canvas.find_discs_indx to set it up.")

        self.discs_ang = [np.asarray(
            hp.pix2ang(self.nside, indx)
            )
            for indx in self.discs_indx]

        print("Done! You can now get the angular position of the discs using Canvas.discs_ang.")

    def find_discs_2center_distance(self):
        """
        Find the angular distance [radians] of disc pixels to the halo center pixel

        Returns
        -------
        None
        """

        # squeeze the disc_ang arrays to remove redundant second dimensions
        # this is necessary at the moment to avoid a bug in healpy.rotator.angdist
        # when calculating the angdist o=for arrays of shape (2,1) and (2,)
        # the returned results is 3 dimensional instead of 1
        # squeezing the array will resolve the issue though
        # TODO: update this and post issue on healpy github

        #FIXME: list comprehension
        self.discs_ang = [np.squeeze(self.discs_ang[halo]) for halo in range(self.catalog.size)]

        #FIXME: list comprehension
        self.discs_2center_rad = [hp.rotator.angdist(self.discs_ang[halo],\
                                                    self.centers_ang[:, halo])
                                for halo in range(self.catalog.size)]

        #FIXME: list comprehension
        self.discs_2center_mpc = [self.centers_D_a[halo]*self.discs_2center_rad[halo]
                                  for halo in range(self.catalog.size)]

    # ------------------------
    #  visualization  methods
    # ------------------------
    #FIXME: remove graticule args

    def _viewer(self,
                map_,
                projection="mollview",
                #graticule=True,
                min=None,
                max=None,
                *args,
                **kwargs):
        """
        wrapper for healpy visualization functions

        Returns
        -------
        None
        """

        # select healpy projection type (e.g. mollview, cartview)
        hp_viewer = self._proj_dict[projection]

        #if graticule: hp.graticule()

        hp_viewer(map_,
                        cmap=self.cmap,
                        min=min,
                        max=max,
                        *args,
                        **kwargs,
                        )

    def show_halo_centers(self,
                          projection="mollweide",
                          graticule=True,
                          ):

        # TODO: implement quiver on the sphere


        try:
            # draw and empty map
            self._viewer(projection=projection,)
        except TypeError:
            # get around the bug in healpy
            map_ = np.zeros(self.npix) + np.inf
            cbar = False
            self._viewer(map_,
                         projection=projection,
                         cbar=cbar)

        if graticule: hp.graticule()

        hp.projscatter(self.catalog.data.theta,
                       self.catalog.data.phi,
                       s=np.log(self.catalog.data.M_200c),
                       )

    def show_discs(self,
                   projection="mollweide",
                   graticule=False,
                   *args,
                   **kwargs,
                   ):

        # select healpy projection type (e.g. mollview, cartview)
        #viewer = self._proj_dict[projection]
        junk_pixels = np.zeros(self.npix)

        def set_to_1(disc):
            junk_pixels[disc] = 1

        [set_to_1(disc) for disc in self.discs_indx]

        if graticule: hp.graticule()

        self._viewer(junk_pixels,
                     projection=projection,
                     *args,
                     **kwargs,
                     )

        del junk_pixels


    def show_map(self,
                 projection="mollweide",
                 graticule=True):

        self._viewer(self.pixels,
                     projection=projection,
                     #graticule=graticule,
                     )
        #TODO: add min max args



#########################################################
#                   Painter Object
#########################################################

class Painter:

    """
    Painter object sprays a signal over the canvas using a template
    """

    def __init__(self,
                 template,
                 ):

        self._template = template
        self._announce_template()

    # ------------------------
    #       properties
    # ------------------------

    # Mutables:

    @property
    def template(self):
        return self._template

    @template.setter
    def template(self, val):
        self._template = val
        self._announce_template()

    # ------------------------
    #         methods
    # ------------------------

    def spray(self, canvas, distance_units="Mpc", **template_kwargs):

        #TODO: check the arg list and if the parameter is not in the catalog add it there

        # check the units
        if distance_units.lower() == "mpc":
            r = canvas.discs_2center_mpc
        elif distance_units.lower() in ["radians", "rad", "rads"]:
            r = canvas.discs_2center_rad
        else:
            raise KeyError("distance_units must be either 'mpc' or 'radians'.")


        #FIXME: list comprehension
        [np.add.at(canvas.pixels, canvas.discs_indx[halo],
                   self.template(r[halo],
                                 canvas.catalog.data.loc[halo],
                                 ))
         for halo in range(canvas.catalog.size)]

        print("Your artwork is fininshed. Check it out with Canvas.show_map()")

    def _announce_template(self):
        """
        Get the template name and list of arguments

        Returns
        -------
        None
        """

        self.template_name = self.template.__name__
        self.template_args_list = inspect.getfullargspec(self.template)[0]

        message = f"The template '{self.template_name}' takes in the following arguments:\n" \
                  f"{self.template_args_list}"
        print(message)


if __name__ == "__main__":

    catalog = Catalog()

    canvas = Canvas(catalog=catalog, nside=32)
    canvas.show_halo_centers()

