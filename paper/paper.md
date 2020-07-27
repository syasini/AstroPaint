---
title: 'AstroPaint: A Python Package for Painting the Sky'
tags:
  - python
  - astrophysics
  - simulation
  - visualization
  - extragalactic foregrounds
authors:
  - name: Siavash Yasini^[corresponding author]
    orcid: 0000-0003-1978-6325
    affiliation: 1 
  - name: Marcelo Alvarez 
    affiliation: "2, 3"
  - name: Emmanuel Schaan 
    orcid: 0000-0002-4619-8927
    affiliation: "2, 3"
  - name: Karime Maamari
    affiliation: "1, 5"
  - name: Shobeir K. S. Mazinani
    affiliation: 4
  - name: Nareg Mirzatuny
    affiliation: 1
  - name: Elena Pierpaoli
    affiliation: 1
affiliations:
 - name: University of Southern California 
   index: 1
 - name: University of California, Berkeley 
   index: 2
 - name: Lawrence Berkeley National Laboratory
   index: 3
 - name: Aetna Inc.
   index: 4
 - name: Argonne National Lab 
   index: 5
date:  31 July 2020
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Overview 

`AstroPaint` is a python package for simulating and visualizing
 mock maps of astrophysical signals. The
  code takes in a halo catalog and the radial profile of an astrophysical
   signal from the user, then combines the two to paint a whole-sky mock
    map of the observable at
    high resolution (see the workflow section for details
    ). `AstroPaint` also contains a
     variety of methods that
     facilitate analysis routines such as data filtering, map manipulation
     , and cutout stacking. The package has an Object-Oriented design and
      runs in parallel, making it both easy to use and scale for
       production of high resolution
       maps with large underlying catalogs.  
          

![Map of the Birkinshaw-Gull effect painted with AstroPaint on top of the
 WebSky catalog \label{fig:BG}](../images/BG_websky_cover.png)

# Statement of need 

Studying the large scale structure of the universe heavily relies on
 observations of astrophysical signals at various frequencies. Examples of such
  studies include detection or characterization of objects such as galaxies
  , clusters, or voids
   through either gravitaional lensing, electromagnetic scattering
   , absorption or emission events in the optical, radio, or x-ray
    frequency bands. Such studies typically require simulated high resolution
     maps of various astrophysical effects to emulate both the signal and
      noise (foregrounds) components. For example, in a study that aims
       to evaluate the detection significance of the Birkinshaw-Gull (BG)
       effect [@Birkinshaw:1983; @Yasini:2018] to measure transverse
        velocities of halos
        using the Simons
        Observatory [@SO:2019] or CMB-S4 [@CMB-S4:2019], one needs a mock
         map of the BG effect
         \autoref{fig:BG} 
        as well as maps of potential contaminants such as kinetic and
         thermal Sunyaev-Zeldovich effects (kSZ and tSZ) [@Sunyaev:1970] for the
         same set of objects. 
     
Creating realistic maps of astrophysical effects through
 hydrodynamical simulations [@Dolag:2015] can be prohibitive for large
  numbers of objects
 . An alternative strategy to creating mock observations of extended objects
  such as galaxies and galaxy cluster halos would be to simulate the
   positions of these objects (either semi-analytically or through N-body
    simulations [@Stein:2020; @Stein:2018; @Sehgal:2010]) and then synthetically
     paint the
     desired signal at the
     location of the halos. `AstroPaint` is developed to accomplish this
      latter step.  
 
# Package Structure and Workflow 


`AstroPaint` consists of three main objects that interact with each other
\: `Catalog`, `Canvas`, and `Painter`. 


`Catalog` contains the locations, velocities, and masses of the objects. 
`Canvas` contains the map of the astrophysical signal in HEALPix format
 [@Healpy:2019]. 
`Painter` contains the template for the radial profile of the signal to be
 painetd on the `Canvas` in circular discs centered at the location of the
  halos in the
  `Catalog`.   

 These objects are sequentially passed into each other according to the
  following workflow: 

```python
from astropaint import Catalog, Canvas, Painter

catalog = Catalog(data=input_data)
canvas = Canvas(catalog, nside)
painter = Painter(template=radial_profile)

painter.spray(canvas)
```

The output map array can be accessed via `canvas.pixels` or directly
 visualized using `canvas.show_map()`. Here `input_data` is the dataframe that
  hold the locations, velocities, and
 masses of the halos. `nside` is a parameter in `healpy` that determines the
  total number of pixels (`npix = 12 * nside ** 2)` and
   consequently the resolution of the map . Finally
  , `radial_profile` is a one-dimensional function that determines the shape
   of the profile. A mind map visualization of the package structure can be
    found in [here](https://www.mindmeister.com/1417665103/astropaint
    -astropaint-py?fullscreen=1).   


# Acknowledgements

SY is grateful to the BCCP group at UC Berkeley for their hospitality during
 Summer 2019 where this project was inaugurated.

# References