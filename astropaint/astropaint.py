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
            self.data = self.generate_random_box()

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
        self.size = len(self._data)
        print("Input data has been modified. Rebuilding the dataframe using "
              "catalog.build_dataframe to update all the parameters...\n")
        self.build_dataframe()


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
        #TODO: change redshift to nonuniversal value
        self.data["rho_s"] = transform.M_200_to_rho_s(self.data["M_200c"],
                                                      self.redshift,
                                                      self.data["R_200c"],
                                                      self.data["c_200c"])

        self.data["R_s"] = np.true_divide(self.data["R_200c"], self.data["c_200c"])

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

        print("Done!")

    @staticmethod
    def _initialize_catalog(n_tot):
        """initialize an empty catalog with x, y, z, v_x, v_y, v_z, M_200c columns"""

        dtype = {"names": ["x", "y", "z", "v_x", "v_y", "v_z", "M_200c"],
                 "formats": 7 * [np.float32]}

        catalog = np.zeros(n_tot, dtype)
        return catalog

    @staticmethod
    def generate_random_box(box_size=100,
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
                 mode="healpy",
                 analyze=True,
                 R_times=1, # the discs will be found around R_times x virial radius
                 ):

        #TODO: define attribute dictionary with __slots__

        assert mode == "healpy", "currently only full sky is supported"

        self._nside = nside
        self._npix = hp.nside2npix(self.nside)
        self._cmap = cm.Greys_r
        self.R_times = R_times

        self.pixels = np.zeros(self.npix)

        self._catalog = catalog
        self.centers_D_a = self._catalog.data.D_a

        if analyze:
            self.analyze()

        #TODO: remove this
        #assert isinstance(catalog, Catalog), "input catalog has to be an instance of " \
        #                                     "astroPaint.Catalog"

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
    def catalog(self):
        return self._catalog

    @catalog.setter
    def catalog(self, val):
        self._catalog = val
        self.analyze()

    @property
    def cmap(self):
        return self._cmap

    @cmap.setter
    def cmap(self, val):
        #FIXME: find the parent class of cm
        assert type(val) is type(cm.Greys), "cmap must be an instance of cm. \n" \
                                            "You can import it using:\n" \
                                            "from matplotlib import cm"
        self._cmap = val
        self._cmap.set_under("white")

    # ------------------------
    #         methods
    # ------------------------

    def analyze(self,):
        """
        Analyze the catalog and find the relevant pixels on the canvas

        Returns
        -------
        None
        """

        self.centers_D_a = self.catalog.data.D_a

        # update the index and angular location of the center pixel
        self.find_centers_indx()
        self.find_centers_ang()

        self.find_discs_indx(self.R_times)
        self.find_discs_ang()
        self.find_discs_2center_distance()


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

    def find_discs_indx(self, R_times):
        """
        Find the pixel indices of discs of size k times R_200 around halo centers

        Parameters
        ----------
        R_times: int
            multiplicative factor indicating the extent of the queried disc in units of R_200

        Returns
        -------
        None

        Sets Canvas.discs_indx to a list of pixel index arrays. Element [i] of the list holds the
        pixel indices around halo [i].

        """

        #FIXME: list comprehension
        self.R_times = R_times
        self.discs_indx = ([np.asarray(
                                    hp.query_disc(self.nside,
                                                  (self.catalog.data.x[halo],
                                    self.catalog.data.y[halo],
                                    self.catalog.data.z[halo]),
                                                  R_times * transform.arcmin2rad(
                                                      self.catalog.data.R_th_200c[halo]))
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
                          marker="o",
                          color=None,
                          s=None,
                          *args,
                          **kwargs,
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
                         cbar=cbar,
                         *args,
                         **kwargs,
                         )

        if graticule: hp.graticule()

        if s is None: s=np.log(self.catalog.data.M_200c),
        hp.projscatter(self.catalog.data.theta,
                       self.catalog.data.phi,
                       color=color,
                       s=s,
                       marker=marker,
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
                 graticule=True,
                 *args,
                 **kwargs):

        self._viewer(self.pixels,
                     projection=projection,
                     #graticule=graticule,
                     *args,
                     **kwargs,
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

        self.template = template
        #self._analyze_template()

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
        self._analyze_template()
        #self._check_template()

    # ------------------------
    #         methods
    # ------------------------

    def spray(self,
              canvas,
              distance_units="Mpc",
              **template_kwargs):

        """
        #TODO: add example

        Parameters
        ----------
        canvas
        distance_units
        template_kwargs

        Returns
        -------

        """
        print("Painting the canvas...")
        #TODO: check the arg list and if the parameter is not in the catalog add it there

        # TODO: check the length and type of the extra_params

        # if it's a scalar dictionary extend it to the size of the catalog
        # also make sure the length matches the size of the catalog


        # convert the template_kwargs into a dataframe
        template_kwargs_df = self._check_template_kwargs(**template_kwargs)
        # use template args to grab the relevant columns from the catalog dataframe
        template_args_df = self._check_template_args(canvas.catalog)

        #TODO: remove this block
        # check the canvas catalog and make sure all the template arguments are already there
        # for parameter in self.template_args_list[1:]:
        #     try:
        #         canvas.catalog.data[parameter]
        #     except KeyError:
        #         try:
        #             template_kwargs[parameter]
        #         except KeyError:
        #             raise KeyError(f"Parameter {parameter} was not found either in the canvas.catalog.data "
        #                   f"or the extra_params.")

        # match the size of the args and kwargs dataframes
        # if template kwargs are scalars, extend then to the size of the catalog
        if template_kwargs_df is None:
            pass
        elif len(template_kwargs_df) == 1:
            template_kwargs_df = pd.concat([template_kwargs_df]*len(template_args_df),
                                           ignore_index=True)

        #TODO: check for other conditions (e.g. longer len, shorter, etc.)

        # concatenate the two dataframes together
        spray_df = pd.concat((template_args_df, template_kwargs_df), axis=1)
        print(f"spray_df.columns = {spray_df.columns}")

        # check the units
        if distance_units.lower() in ["mpc", "megaparsecs", "mega parsecs"]:
            r = canvas.discs_2center_mpc
        elif distance_units.lower() in ["radians", "rad", "rads"]:
            r = canvas.discs_2center_rad
        else:
            raise KeyError("distance_units must be either 'mpc' or 'radians'.")


        #TODO: this has been checked elsewhere... remove it
        # make sure r (distance) is in the argument list
        assert 'r' in self.template_args_list

        #TODO: think about how to redo this part
        if len(self.template_args_list) == 1:

            #FIXME: list comprehension
            [np.add.at(canvas.pixels,
                       canvas.discs_indx[halo],
                       self.template(r[halo]))
             for halo in range(canvas.catalog.size)]

        else:
            #FIXME: list comprehension
            [np.add.at(canvas.pixels,
                       canvas.discs_indx[halo],
                       self.template(r[halo],
                                     **spray_df.loc[halo]))
             for halo in range(canvas.catalog.size)]

        print("Your artwork is fininshed. Check it out with Canvas.show_map()")



    def _analyze_template(self):
        """
        Get the template name and list of arguments

        Returns
        -------
        None
        """

        self.template_name = self.template.__name__

        # get the list of args and keyword args
        self.template_args_list = inspect.getfullargspec(self.template).args
        self.template_kwargs_list = inspect.getfullargspec(self.template).kwonlyargs

        # print out the list of args and kwargs
        message = f"The template '{self.template_name}' takes in the following arguments:\n" \
                  f"{self.template_args_list}\n" \
                  f"and the following keyword-only arguments:\n" \
                  f"{self.template_kwargs_list}"

        # ensure the first argument of the profile template is 'r'
        assert self.template_args_list[0] == "r", "The first argument of the profile template " \
                                                  "must be 'r' (the distance from the center of " \
                                                  "the halo)."
        print(message)

    def _check_template_kwargs(self, **template_kwargs):
        """Ensure the template_kwargs is pandas compatible"""

        if template_kwargs:
            try:
                #TODO: find the type of input (scalar, array, DF)?
                for key, value in template_kwargs.items():
                    if not hasattr(value, "__len__"):
                        template_kwargs[key] = [value]

                template_kwargs_df = pd.DataFrame(template_kwargs)
                return template_kwargs_df

                #self.template_kwargs_data = pd.DataFrame(template_kwargs)
            except:
                raise
        else:
            #TODO: add warning if template has kwargs but no template_kwargs are provided
            print("No template_kwargs provided")

            return None


    def _check_template_args(self, catalog):
        """Check to see if the template profile function arguments exist in the catalog """

        # check the canvas catalog and make sure all the template arguments are already there
        params_not_found = []
        #params_not_found_anywhere = []
        for parameter in self.template_args_list[1:]:
            try:
                catalog.data[parameter]
            except KeyError:
                params_not_found.append(parameter)

        if len(params_not_found) > 0:
            print("The following parameters were not found in the canvas.catalog.data\n"
                  f"{params_not_found}\n"
                  "Make sure you pass them as kwargs (key=value), dictionary (**dict), or Pandas "
                  "DataFrame (**df) in the .spray method. Check the spray docstring"
                  "(.spray.__doc__) for examples. ")


        parameters = list(set(self.template_args_list[1:]) - set(params_not_found))

        template_args_df = catalog.data[parameters]
        return template_args_df

if __name__ == "__main__":

    catalog = Catalog()

    canvas = Canvas(catalog=catalog, nside=32)
    canvas.show_halo_centers()

