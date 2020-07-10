"""
This is a cool webapp for painting the sky using AstroPaint
https://github.com/syasini/AstroPaint
"""

__author__ = ["Siavash Yasini", "Shobeir K. S. Mazinani"]
__email__ = ["siavash.yasini@gmail.com"]

import streamlit as st
from matplotlib import cm
import matplotlib.pyplot as plt
import healpy as hp
from astropaint.profiles import art_gallery
from astropaint import Catalog, Canvas, Painter

nside = 256  # map resolution
fig_dpi = 150
template_dict = {"drops": art_gallery.drops,
                 "bacteria": art_gallery.bacteria,
                 "twilight": art_gallery.twilight}

st.markdown(
"""
# AstroPaint Art Gallery 
Created by [Siavash Yasini](https://github.com/syasini) & [Shobeir K. S. Mazinani](
https://github.com/S-KSM)


Check out [AstroPaint](https://github.com/syasini/AstroPaint) on GitHub and 
press that star button to show you support! ️🤩👉⭐ 
"""
    )

st.sidebar.markdown(
    """
    ## How to paint like a pro 
1. Randomly pick a configuration 🎨🖌
4. Press Paint! 😎
    """
    )


# ======================
# Catalog configurations
# ======================

st.sidebar.markdown(
    """
    ## Catalog
    """
    )

# choose the number of halos
n_tot = st.sidebar.number_input("number of halos",
                      value=500,
                      min_value=0,
                      max_value=10000)

# select the shell radius (distance from center)
shell_radius = st.sidebar.number_input("shell radius [mega parsecs]",
                      value=40,
                      min_value=10,
                      max_value=500)

# instantiate the catalog object with the given parameters
catalog = Catalog()
catalog.generate_random_shell(n_tot=n_tot, shell_radius=shell_radius)

# cutout angular sizes if prompted by the user
if st.sidebar.checkbox(label="limit angular sizes"):

    R_ang_min_max = st.sidebar.slider("angular radius [arcmins]",
                      value=[0, 120],
                      min_value=0,
                      max_value=240)
    catalog.cut_R_ang_200c(*R_ang_min_max)

if catalog.size == 0:
    st.sidebar.warning("No halos found in this range")


# =====================
# Canvas configurations
# =====================

st.sidebar.markdown(
    """
    ## Canvas
    """
    )

# determine how far from the center we want to paint
R_times = st.sidebar.slider("paint up to R x ",
                      value=4,
                      min_value=1,
                      max_value=10)

# instantiate the canvas
canvas = Canvas(catalog, nside=nside, R_times=R_times)


show_graticule = st.sidebar.checkbox("show canvas graticule")

if st.sidebar.button("Clear Canvas"):
    canvas.clean()


# ======================
# Painter configurations
# ======================
st.sidebar.markdown(
    """
    ## Painter
    """
    )

# select the template
template = st.sidebar.selectbox(label="radial profile template",
                                options=["drops",
                                         "bacteria",
                                         "twilight",
                                         "custom"])


if template == "custom":
    st.sidebar.warning("coming soon! select another profile...")
    template = template_dict["drops"]
else:
     template = template_dict[template]

painter = Painter(template)

# select the colormap
cmap = st.sidebar.selectbox(label="color map",
                            options=["Blues",
                                     "twilight",
                                     "Greys_r",
                                     "RdBu_r",
                                     "custom"])

if cmap == "custom":
    custom_cmap = st.sidebar.text_input("matplotlib cmap", "Oranges")
    canvas.cmap = cm.get_cmap(custom_cmap)
else:
    canvas.cmap = cm.get_cmap(cmap)


#
if st.sidebar.button("Paint"):
    with st.spinner(f"Painting {catalog.size} halos..."):
        painter.spray(canvas)
        st.sidebar.balloons()
    st.success("Done!")


plt.figure(dpi=fig_dpi)

cmap_dict = {"Blues": cm.Blues,
             "RdBu_r": cm.RdBu_r}



canvas.show_map("mollview", fig=1, title="", cbar=False)
#hp.graticule(dmer=360,dpar=360,alpha=0)


if show_graticule:
    hp.graticule()
st.pyplot()

#
# lonra = st.slider("longitude [degree]",
#           value=[-20, 20],
#           step=5,
#           min_value=-180,
#           max_value=180)
#
#
# latra = st.slider("latitude [degree]",
#           value=[-20, 20],
#           step=5,
#           min_value=-90,
#           max_value=90)
#
# lonra = [*lonra]
# latra = [*latra]

canvas.show_map("cart", lonra=[-40, 40], latra=[-40, 40], fig=1, title="", cbar=False)
if show_graticule:
    hp.graticule()
st.pyplot()

st.sidebar.markdown(
    """
    ---
    
    ### Examples to try
    
    
    **drops**
    - *number of halos*: 300
    - *shell radius*: 20
    - *R x* : 6 
    - *colormap*: Blues
    
    **twilight**
    - *number of halos*: 5000
    - *shell radius*: 350
    - *R x* : 10
    - *colormap*: guess!
    
    **bacteria**
    - *number of halos*: 1000
    - *shell radius*: 40
    - *R x* : 3 
    - *colormap*: Greys_r
    
    
    """
    )


