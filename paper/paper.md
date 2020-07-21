---
title: 'AstroPaint: A Python Package for Painting the Sky'
tags:
  - python
  - astrophysics
  - simulation
  - extragalactic foregrounds
authors:
  - name: Siavash Yasini^[corresponding author]
    orcid: 0000-0003-1978-6325
    affiliation: 1 
  - name: Marcelo Alvarez 
    affiliation: "2, 3"
  - name: Emmanuel Schaan 
    affiliation: "2, 3"
affiliations:
 - name: University of Southern California 
   index: 1
 - name: University of California, Berkeley 
   index: 2
 - name: Lawrence Berkeley National Labratory
   index: 3
date:  31 July 2020
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Definitions 

halo catalog:
profile: 


# Summary

`AstroPaint` is a python package for simulating and visualizing
 mock maps of astrophysical signals. The
  code takes in a halo catalog and the radial profile of an astrophysical
   signal from user and then combines the two to paint a whole-sky mock map of
    the observable at
    high resolution \autoref{sec:workflow}. `AstroPaint` also contains a
     variety of methods that
     facilitate analysis routines such as data filtering, map manipulation
     , and cutout stacking. The package has an Object-Oriented structure and
      runs in parallel, making it scalable for use with large catalogs and
       high resolutions.  
          
     

![BG](../images/BG_websky_cover.png)

# Statement of need 
whole-sky high resolution maps of astrophysical signals at various frequencies
 extragalactc foregrounds
 
correlations 
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

# Package Structure and Workflow \label{sec:workflow}

`AstroPaint` consists of three main objects that interact with each other
: `Catalog`, `Canvas`, and `Painter`. 
*`Catalog`* contains the locations, velocities, and masses of the objects. 
`Canvas` contains the map of the astrophysical signal in healpix format. 
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

Here `input_data` is the dataframe that holds the locations, velocities, and
 masses of the halos. `nside` is a parameter in `healpy` that determines the
  resolution of the map (`npix = 12 * nside ** 2)`. And finally
  , `radial_profile` is a one-dimensional function that determines the shape
   of the profile.   



# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Fenced code blocks are rendered with syntax highlighting:
```python
for n in range(10):
    yield f(n)
```	

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References