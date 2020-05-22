from astropaint import Catalog, Canvas, Painter
from astropaint.profiles import NFW, Battaglia16
from astropaint.lib.utils import timeit
from astropy.cosmology import Planck18_arXiv_v2 as cosmo
import matplotlib.pyplot as plt

h = cosmo.h
nside = 8192
n_halos = None  # Set to an integer to sample the catalog (e.g. 5000)

lon_range = [0, 2]  # degrees
lat_range = [0, 2]  # degrees

mass_cut = 5E12  # M_sun

with timeit("painting the websky catalog"):

	catalog = Catalog(data="websky_2x2")

	# adjust the masses
	catalog.data.M_200c *= h

	# apply mass and fov cuts
	catalog.cut_M_200c(mass_min=mass_cut)
	catalog.cut_lon_lat(lon_range=lon_range, lat_range=lat_range)

	if n_halos is not None:
		catalog.data = catalog.data.sample(n=n_halos)
	print(f"Minimum mass of halos in the catalog: {catalog.data.M_200c.min():.2e}")
	print(f"Number of halos to paint: {catalog.size}")

	canvas = Canvas(catalog,
					nside=nside,
					R_times=5)

	## to use the NFW model use:
	painter = Painter(template=NFW.kSZ_T)

	## to use the Battaglia16 model use:
	# painter = Painter(template=Battaglia16.kSZ_T)

	painter.spray(canvas,
				  #cache=True, ## uncomment to cache the profiles to disc; makes the painting slower
				   )

	canvas.save_map_to_file(prefix="websky_NFW", suffix="2x2")

	canvas.show_map("cartesian",
					lonra=lon_range, latra=lat_range,
					title=f"{painter.template_name}",
					min=-1E-5,
					max=1E-5
					)
	plt.show()
