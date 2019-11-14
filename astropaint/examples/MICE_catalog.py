from astropaint import Catalog, Canvas, Painter
from astropaint import profile
from astropaint.lib.mytimer import timeit
nside = 256 
with timeit("painting the MICE catalog"):
	catalog = Catalog(data="MICE")

	# slice the catalog
	#catalog.data = catalog.data[(catalog.data.lon<5) & (catalog.data.lat<5)]
	catalog.data = catalog.data.sample(frac=0.0001)
	catalog.move_to_box_center()
	canvas = Canvas(catalog,
			nside=nside,
			R_times=2)

	painter = Painter(template=profile.kSZ_T_NFW)

	painter.spray(canvas)

	canvas.save_map_to_file(prefix="MICE_lite_", suffix="test")
	canvas.save_Cl_to_file(prefix="MICE_lite", suffix="_test")

