"""
library for simulating semi-analytic mock maps of CMB secondary anisotropies
"""

__author__ = "Siavash Yasini"
__email__ = "yasini@usc.edu"

import os
import numpy as np
import pandas as pd
from matplotlib import cm
from warnings import warn
import inspect
from itertools import product
import operator
from memory_profiler import profile

try:
    import healpy as hp

except ModuleNotFoundError:
    warn("Healpy is not installed. You cannot use the full sky canvas without it.")

#import sys
#print(sys.path)
from astropy.coordinates import cartesian_to_spherical
from .lib import transform

# find the package path; same as __path__
path_dir = os.path.dirname(os.path.abspath(__file__))

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
                 ):

        #TODO: define attribute dictionary with __slots__

        self.redshift = redshift
        # if no input is provided generate a random catalog
        if data is None:
            self.data = self.generate_random_box()
        elif isinstance(data, str):
            self.load_sample(data)
        else:
            #FIXME: check data type and columns
            self.data = data

        # .................
        # octant signatures
        # .................

        # (x,y,z) signatures for each octant e.g. (+,+,+) , (+,+,-) etc.
        self.octant_signature = self._get_octant_signatures(mode="user")

        # same thing but for use in calculations
        self._octant_shift_signature = self._get_octant_signatures(mode="shift")
        self._octant_mirror_signature = self._get_octant_signatures(mode="mirror")
        self._octant_rotate_signature = self._get_octant_signatures(mode="rotate")

        # TODO: check input type/columns/etc

    # ------------------------
    #       properties
    # ------------------------

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, val):
        self._data = val
        self._data = pd.DataFrame(self.data).reset_index(drop=True)

        self.size = len(self._data)
        self.box_size = self._get_box_size()
        print("Input data has been modified. Rebuilding the dataframe using "
              "catalog.build_dataframe to update all the parameters...\n")

        # build the complete data frame
        # e.g. angular distances, radii, etc.
        self.build_dataframe()

    # ------------------------
    #         sample data
    # ------------------------

    #TODO: support inputs other than csv
    def load_sample(self, sample_name="MICE"):
        """load sample data using the name of dataset"""
        fname = os.path.join(path_dir, "data", f"{sample_name}.csv")

        self.data = pd.read_csv(fname, index_col=0)

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

    def _get_box_size(self):
        """find the catalog box size from x, y, z coordinates"""

        Lx = self.data["x"].max() - self.data["x"].min()
        Ly = self.data["y"].max() - self.data["y"].min()
        Lz = self.data["z"].max() - self.data["z"].min()

        return Lx, Ly, Lz

    @staticmethod
    def _get_octant_signatures(mode="user"):
        """calculate the octant signatures to be used in replication function later"""

        # set up the octant signature with +1, -1 indicating the sign of each axis
        # e.g. (+1,+1,+1) is the first octant/ (-1,+1,+1) is the second octant etc.
        x_signs = np.sign(np.exp(1.j * (np.arange(8) * np.pi / 2 + np.pi / 4)).real).astype(int)
        y_signs = np.sign(np.exp(1.j * (np.arange(8) * np.pi / 2 + np.pi / 4)).imag).astype(int)
        z_signs = np.array(4 * [1] + 4 * [-1])

        # put them together as a reference dictionary
        oct_sign_dict = dict(enumerate(zip(x_signs, y_signs, z_signs)))

        if mode == "user":
            sign_dict = {1: "+",
                         -1: "-"}
            # (x,y,z) signatures for each octant e.g. (+,+,+) , (+,+,-) etc.
            octant_signature = [(sign_dict[i], sign_dict[j], sign_dict[k])
                                for (i, j, k) in oct_sign_dict.values()]

        elif mode == "shift":
            sign_dict = {1 : 0,
                         -1: -1}
            octant_signature = [(sign_dict[i], sign_dict[j], sign_dict[k])
                                for (i, j, k) in oct_sign_dict.values()]


        elif mode == "mirror":
            octant_signature = product((+1, -1), repeat=3)

        elif mode == "rotate":
            # octant signature for replication by rotation
            octant_signature = sorted(product([0, 1, 2, 3], [1, -1]),
                                                              key=operator.itemgetter(1),
                                                              reverse=True)
        else:
            raise KeyError("octant signature mode not defined")

        octant_signature = dict(enumerate(octant_signature))
        return octant_signature

    @staticmethod
    def _initialize_catalog(n_tot):
        """initialize an empty catalog with x, y, z, v_x, v_y, v_z, M_200c columns"""

        dtype = {"names": ["x", "y", "z", "v_x", "v_y", "v_z", "M_200c"],
                 "formats": 7 * [np.float32]}

        catalog = np.zeros(n_tot, dtype)
        return catalog

    @staticmethod
    def generate_random_box(box_size=50,
                            v_max=100,
                            mass_min=1E14,
                            mass_max=1E15,
                            n_tot=50000,
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

    @staticmethod
    def generate_test_box(configuration=["front"],
                          distance=100,
                          mass=1E15,
                          ):

        catalog = pd.DataFrame(Catalog._initialize_catalog(0))
        config_dict = {"front": (1, 0, 0),
                       "back": (-1, 0, 0),
                       "left": (0, 1, 0),
                       "right": (0, -1, 0),
                       "top": (0, 0, 1),
                       "bottom": (0, 0, -1),
                       }

        # set configuration for "all" keyword
        if "all" in configuration:
            configuration = config_dict.keys()

        for key in configuration:
            # get the coordinates from config_dic and load it in a dataframe
            x, y, z = config_dict[key]
            df = pd.DataFrame(Catalog._initialize_catalog(1))
            df["x"], df["y"], df["z"] = x, y, z
            df[["x", "y", "z"]] *= distance

            # set the mass
            df["M_200c"] = mass

            # append the test case to the catalog
            catalog = catalog.append(df, ignore_index=True)

        return pd.DataFrame(catalog)  # return the pandas dataframe

    @staticmethod
    def _set_octant(df, octant):
        """Affix an octant column to a copy of the data frame """
        df_copy = df.copy() #FIXME: Make sure shallow copy is safe
        df_copy["octant"] = octant
        return df_copy

    @staticmethod
    def _tile_by_shifting(coords, boxsize, move_signature):
        """tile a 3d box by shifting coordinates according to the move signatures in
        _octant_shift_signature
        e.g. (0,0,-1) shifts the box down one unit in the -z direction"""

        #TODO: assert 3d arrays
        #TODO: assert array move signature
        x, y, z = coords
        dx, dy, dz = move_signature * boxsize

        x += dx
        y += dy
        z += dz

        return pd.Series([x, y, z])

    @staticmethod
    def _tile_by_mirroring(coords, boxsize, move_signature):
        """tile a 3d box by reflecting coordinates according to the move signatures in
        _octant_mirror_signature
        e.g. (+1,+1,-1) reflects the box along the x-y plane (in the -z direction)"""

        # TODO: assert 3d arrays
        # TODO: assert array move signature
        x, y, z = coords
        # boxsize is redundant and is used for consistency
        dx, dy, dz = move_signature

        x *= dx
        y *= dy
        z *= dz

        return pd.Series([x, y, z])

    @staticmethod
    def _tile_by_rotating(coords, boxsize, move_signature):
        """tile a 3d box by rotating coordinates according to the move signatures in
        _octant_rotate_signature
        e.g. (1,1) rotates the box counter-clockwise around the z axis"""

        # TODO: assert 3d arrays
        # TODO: assert array move signature
        x, y, z = coords
        # boxsize is redundant and is used for consistency
        n, z_sign = move_signature

        # rotate the coordinates according to the transformation in Sehgah 2010 Eq 28
        xiy = (x + np.sign(z_sign)*1.j*y)*np.exp(1.j*(n-0.5*np.sign(z_sign)+0.5)*np.pi/2)
        x = xiy.real
        y = xiy.imag
        z *= z_sign

        return pd.Series([x, y, z])

    def replicate(self,
                  mode="rotate"):
        """
        Replicate an octant to get a whole-sky box

        Parameters
        ----------
        mode

        Returns
        -------

        """
        assert mode in ["shift", "rotate", "mirror"], "mode can be either 'shift', 'rotate', " \
                                                      "or 'mirror."

        # add the octant column to data
        # define local variable data to force catalog rebuilding at the end
        data = pd.concat([self._set_octant(f, i) for (i, f) in enumerate(8*[self.data])],
                              axis=0, ignore_index=True)

        # set the replication mode and signature based on the given kwarg
        if mode.lower() == "shift":
            tile = self._tile_by_shifting
            move_signature = self._octant_shift_signature
        elif mode.lower() == "mirror":
            tile = self._tile_by_mirroring
            move_signature = self._octant_mirror_signature
        elif mode.lower() == "rotate":
            tile = self._tile_by_rotating
            move_signature = self._octant_rotate_signature

        # replicate the octants using the tiling function set above
        data[["x", "y", "z"]] = \
            data.apply(lambda row:
                            tile(
                                row[["x", "y", "z"]],
                                self.box_size,
                                np.array(move_signature[row["octant"]])
                                ),
                            axis=1)

        # reset data and rebuild the dataframe
        self.data = data

    def move_to_box_center(self):
        """move the observer from (0,0,0) to the center of the box (Lx/2, Ly/2, Lz/2) to make
        coordinates symmetric
        *Not recommended for light-cone catalogs"""

        data = self.data  # trick for forcing catalog rebuilding at the end
        Lx, Ly, Lz = self.box_size

        data["x"] -= Lx / 2
        data["y"] -= Ly / 2
        data["z"] -= Lz / 2

        # reset data and rebuild the dataframe
        self.data = data

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
                 R_times=1,  # the discs will be found around R_times x virial radius,
                 inclusive=False,
                 ):

        #TODO: define attribute dictionary with __slots__

        assert mode == "healpy", "currently only full sky is supported"

        self._nside = nside
        self._npix = hp.nside2npix(self.nside)
        self._lmax = 3 * self.nside-1
        self._ell = np.arange(self.lmax+1)
        self._cmap = cm.Greys_r
        self.R_times = R_times
        self.inclusive = inclusive

        self._pixels = np.zeros(self.npix)
        self._Cl = np.zeros(self.lmax+1)
        self._Cl_is_outdated = False

        self._catalog = catalog
        self.centers_D_a = self._catalog.data.D_a

        self.generate_discs()

        if analyze:
            self.discs.analyze()

        #TODO: remove this
        #assert isinstance(catalog, Catalog), "input catalog has to be an instance of " \
        #                                     "astroPaint.Catalog"

        self._proj_dict = {"mollweide": hp.mollview,
                           "mollview": hp.mollview,
                           "cartesian": hp.cartview,
                           "cartview": hp.cartview,
                           }

        self.template_name = None

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

    @property
    def lmax(self):
        return self._lmax

    @property
    def ell(self):
        return self._ell

    @property
    def Cl(self):
        if self._Cl_is_outdated:
            self.get_Cl()
        return self._Cl

    @property
    def Dl(self):
        Dl = self.ell*(self.ell+1)*self.Cl/(2*np.pi)
        return Dl

    # Mutables:

    @property
    def catalog(self):
        return self._catalog

    @catalog.setter
    def catalog(self, val):
        self._catalog = val
        self.discs.analyze()

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

    @property
    def pixels(self):
        return self._pixels

    @pixels.setter
    def pixels(self, val):
        self._pixels = val
        self._Cl_is_outdated = True

    # ------------------------
    #     Disc inner class
    # ------------------------

    class Disc:
        __slots__ = ['catalog',
                     'nside',
                     'R_times',
                     'inclusive',
                     'center_D_a',
                     'center_index',
                     'center_ang',
                     'center_vec',
                     'pixel_index',
                     'pixel_ang',
                     'pix2cent_rad',
                     'pix2cent_mpc',
                     'pix2cent_vec',
                     ]

        def __init__(self,
                     catalog,
                     nside,
                     R_times,
                     inclusive,
                     ):



            #FIXME: the whole catalog does not need to be passed to disc here
            #Check if this affects performance
            self.catalog = catalog
            self.nside = nside
            self.R_times = R_times
            self.inclusive = inclusive

            self.center_D_a = self.catalog.data.D_a

        # ------------------------
        #       finder methods
        # ------------------------

        #FIXME: decide whether to keep this or discard it
        def analyze(self,):
            """
            Analyze the catalog and find the relevant pixels on the canvas

            Returns
            -------
            None
            """


            # update the index and angular location of the center pixel
            # self.find_centers_indx()
            # self.find_centers_ang()
            #
            # self.find_discs_indx(self.R_times)
            # self.find_discs_ang()
            # self.find_discs_2center_distance()


        def find_centers_indx(self):
            """
            Find the pixel indices of the halo centers

            Returns
            -------
            None
            Sets Canvas.centers_indx to array of pixels.
            Element [i] of the array points to the center of halo [i].

            """

            self.center_index = hp.ang2pix(self.nside,
                                           self.catalog.data.theta.to_list(),
                                           self.catalog.data.phi.to_list())

            print("Done! You can now get the center pixels using Canvas.centers_indx.")

        def find_centers_ang(self):
            """
            Store the theta and phi coordinates of the halos in
            Canvas.centers_ang

            Returns
            -------
            None
            """

            self.center_ang = np.asarray([self.catalog.data.theta.to_list(),
                                           self.catalog.data.phi.to_list()])

            print(
                "Done! You can now get the angular position of the discs using Canvas.centers_ang.")

        def find_centers_vec(self):
            """
            Find the unit vectors pointing to the halo centers

            Returns
            -------
            None
            Sets Canvas.centers_vec to array of pixels.
            Element [i] of the array points to the center of halo [i].

            """

            self.centers_vec = hp.ang2vec(self.catalog.data.theta.to_list(),
                                          self.catalog.data.phi.to_list())

            print("Done! You can now get the center pixel vectors using Canvas.centers_vec.")

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

            Sets Canvas.discs_indx to a list of pixel index arrays. Element [i] of the list holds
            the
            pixel indices around halo [i].

            """

            # FIXME: list comprehension
            self.R_times = R_times
            self.pixel_index = ([np.asarray(
                hp.query_disc(self.nside,
                              (self.catalog.data.x[halo],
                               self.catalog.data.y[halo],
                               self.catalog.data.z[halo]),
                              R_times * transform.arcmin2rad(
                                  self.catalog.data.R_th_200c[halo]),
                              inclusive=self.inclusive,
                              )
                )
                for halo in range(self.catalog.size)])

            print("Done! You can now get the discs using Canvas.discs_indx.")

        def find_discs_ang(self):
            """
            Find the angular coordinates of the disc pixels

            Returns
            -------
            None
            """
            try:
                self.pixel_index
            except AttributeError:
                print("Canvas.discs_indx is not defined. Use Canvas.find_discs_indx to set it up.")

            # FIXME: list comprehension
            self.pixel_ang = [np.asarray(
                hp.pix2ang(self.nside, indx)
                )
                for indx in self.pixel_index]

            print("Done! You can now get the angular position of the discs using Canvas.discs_ang.")

        def find_discs_vec(self):
            """
            Find the unit vectors pointing to the disc pixels

            Returns
            -------
            None
            """
            try:
                self.discs_indx
            except AttributeError:
                print("Canvas.discs_indx is not defined. Use Canvas"
                  ".find_discs_indx to set it up.")

            # FIXME: list comprehension
            self.discs_vec = [np.asarray(
                hp.pix2vec(self.nside, indx)
                ).T
                              for indx in self.discs_indx]

            print("Done! You can now get the vectots pointing to the disc pixels using "
                  "Canvas.discs_vec.")

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

            # FIXME: list comprehension
            self.pixel_ang = [np.squeeze(self.pixel_ang[
                                             halo]) for halo in range(self.catalog.size)]

            # FIXME: list comprehension
            self.pix2cent_rad = [hp.rotator.angdist(self.pixel_ang[halo],
                                                         self.center_ang[:, halo])
                                      for halo in range(self.catalog.size)]

            # FIXME: list comprehension
            self.pix2cent_mpc = [self.center_D_a[halo] * self.pix2cent_rad[halo]
                                      for halo in range(self.catalog.size)]

        def find_discs_2center_vec(self):
            """
                Find the 3D unit vector pointing from the disc pixels to the halo center pixel

                Returns
                -------
                None
                """

            # if discs_vec does not exist, find it
            try:
                self.discs_vec
            except AttributeError:
                self.find_discs_vec()

            # if centers_vec does not exist, find it
            try:
                self.centers_vec
            except AttributeError:
                self.find_centers_vec()

            #FIXME: list comprehension
            self.discs_2center_vec = [Canvas._normalize_vec(self.discs_vec[halo] -
                                                          self.centers_vec[halo],
                                          axis=-1)
                                          for halo in range(self.catalog.size)]

        # ------------------------
        #    generator methods
        # ------------------------

        def gen_center_index(self, halo_list="All"):
            if halo_list is "All":
                halo_list = range(self.catalog.size)

            assert hasattr(halo_list, '__iter__')

            for halo in halo_list:
                yield hp.ang2pix(self.nside,
                                 self.catalog.data.theta[halo],
                                 self.catalog.data.phi[halo])

        def gen_center_ang(self, halo_list="All"):
            if halo_list is "All":
                halo_list = range(self.catalog.size)

            assert hasattr(halo_list, '__iter__')

            #TODO: check if this is faster with pandas .iterrows or .itertuples
            for halo in halo_list:
                yield (self.catalog.data.theta[halo],
                       self.catalog.data.phi[halo])

        def gen_center_vec(self, halo_list="All"):
            if halo_list is "All":
                halo_list = range(self.catalog.size)

            assert hasattr(halo_list, '__iter__')

            for halo in halo_list:
                yield hp.ang2vec(self.catalog.data.theta[halo],
                                 self.catalog.data.phi[halo])

        def gen_pixel_index(self, halo_list="All"):
            if halo_list is "All":
                halo_list = range(self.catalog.size)

            assert hasattr(halo_list, '__iter__')

            for halo in halo_list:
                yield hp.query_disc(self.nside,
                              (self.catalog.data.x[halo],
                               self.catalog.data.y[halo],
                               self.catalog.data.z[halo]),
                              self.R_times * transform.arcmin2rad(
                                  self.catalog.data.R_th_200c[halo]),
                              inclusive=self.inclusive,
                              )

        def gen_pixel_ang(self, halo_list="All"):
            if halo_list is "All":
                halo_list = range(self.catalog.size)

            assert hasattr(halo_list, '__iter__')

            for index in self.gen_pixel_index(halo_list):
                yield hp.pix2ang(self.nside, index)

        def gen_pixel_vec(self, halo_list="All"):
            """
            generate the unit vectors pointing to the disc pixels

            Returns
            -------
            None
            """
            if halo_list is "All":
                halo_list = range(self.catalog.size)

            assert hasattr(halo_list, '__iter__')

            for index in self.gen_pixel_index(halo_list):
                yield np.asarray(hp.pix2vec(self.nside, index)).T

        def gen_cent2pix_rad(self, halo_list="All"):
            if halo_list is "All":
                halo_list = range(self.catalog.size)

            assert hasattr(halo_list, '__iter__')

            for (pixel_ang, center_ang) in zip(self.gen_pixel_ang(halo_list),
                                               self.gen_center_ang(halo_list)):
                yield hp.rotator.angdist(pixel_ang, center_ang)

        def gen_cent2pix_mpc(self, halo_list="All"):
            if halo_list is "All":
                halo_list = range(self.catalog.size)

            assert hasattr(halo_list, '__iter__')

            for (halo, pix2cent_rad) in zip(halo_list, self.gen_cent2pix_rad(halo_list)):
                yield self.center_D_a[halo] * pix2cent_rad

        def gen_cent2pix_hat(self, halo_list="All"):
            if halo_list is "All":
                halo_list = range(self.catalog.size)

            for (pix_vec, cent_vec) in zip(self.gen_pixel_vec(halo_list),
                                           self.gen_center_vec(halo_list)):
                yield Canvas._normalize_vec(pix_vec - cent_vec)

    def generate_discs(self):
        """instantiate the discs attribute using the Disc class
        Useful when the disc generators are exhausted and need to be reset"""

        self.discs = self.Disc(self.catalog, self.nside, self.R_times, self.inclusive)

    @staticmethod
    def _normalize_vec(vec, axis=-1):
        """normalize_vec the input vector along the given axis"""

        norm = np.linalg.norm(vec, axis=axis)

        return np.true_divide(vec, np.expand_dims(norm, axis=axis))

    def clean(self):
        """
        Clean the canvas and set all pixels to zero

        Returns
        -------
        None
        """

        self.pixels = np.zeros(self.npix)

    def get_Cl(self):
        """find the power spectrum of the map (.pixels)"""

        self._Cl = hp.anafast(self.pixels,
                              lmax=self.lmax)

        self._Cl_is_outdated = False

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

        if graticule:
            hp.graticule()

        if s is None:
            s = np.log(self.catalog.data.M_200c)

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

        [set_to_1(disc) for disc in self.discs.gen_pixel_index()]

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

    # ------------------------
    #  saving/loading  methods
    # ------------------------

    def save_map_to_file(self,
                         prefix=None,
                         suffix=None,
                         filename=None):
        """save the healpy map to file

        Parameters
        ----------
        prefix: str
            prefix string to be added to the beginning of the default file name

        suffix: str
            suffix string to be added to the end of default file name

        filename: str
            custom file name; overrides the prefix and suffix and default file name
        """

        if prefix:
            if str(prefix)[-1] != "_":
                prefix = "".join([prefix, "_"])
        if suffix:
            if str(suffix)[0] != "_":
                suffix = "".join(["_", suffix])

        if filename is None:
            #TODO: complete this
            filename = f"{str(prefix or '')}" \
                       f"{self.template_name}" \
                       f"_NSIDE={self.nside}" \
                       f"{str(suffix or '')}" \
                       f".fits"

        hp.write_map(filename,
                     self.pixels)


    def save_Cl_to_file(self,
                         prefix=None,
                         suffix=None,
                         filename=None):
        """save the map power spectrum to file

        Parameters
        ----------
        prefix: str
            prefix string to be added to the beginning of the default file name

        suffix: str
            suffix string to be added to the end of default file name

        filename: str
            custom file name; overrides the prefix and suffix and default file name
        """
        if prefix:
            if str(prefix)[-1] != "_":
                prefix = "".join([prefix, "_"])
        if suffix:
            if str(suffix)[0] != "_":
                suffix = "".join(["_", suffix])

        if filename is None:
            #TODO: complete this
            filename = f"{str(prefix or '')}" \
                       f"{self.template_name}" \
                       f"_NSIDE={self.nside}" \
                       f"{str(suffix or '')}"

        print(filename)
        np.savez(filename,
                 ell=self.ell,
                 Cl=self.Cl,
                 Dl=self.Dl)

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

        #
        canvas.template_name = self.template.__name__

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
            r_pix2cent = canvas.discs.gen_cent2pix_mpc
        elif distance_units.lower() in ["radians", "rad", "rads"]:
            r_pix2cent = canvas.discs.gen_cent2pix_rad
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

        #TODO: unify this with the other two conditions
        elif 'r_hat' in self.template_args_list:
            r_hat = canvas.discs.gen_cent2pix_hat

            [np.add.at(canvas.pixels,
                       pixel_index,
                       self.template(r,
                                     r_hat,
                                     **spray_df.loc[halo]))
            for halo, r, r_hat, pixel_index in zip(range(canvas.catalog.size),
                                             r_pix2cent(),
                                             r_hat(),
                                             canvas.discs.gen_pixel_index())]

        else:
            #FIXME: list comprehension
            [np.add.at(canvas.pixels,
                       pixel_index,
                       self.template(r,
                                     **spray_df.loc[halo]))
            for halo, r, pixel_index in zip(range(canvas.catalog.size),
                                             r_pix2cent(),
                                             canvas.discs.gen_pixel_index())]

        print("Your artwork is fininshed. Check it out with Canvas.show_map()")

        # acticate the canvas.pixels setter
        canvas.pixels = canvas.pixels


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



