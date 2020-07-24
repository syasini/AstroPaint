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
    orcid: 0000-0002-4619-8927
    affiliation: "2, 3"
  - name: Karime Maamari
    affiliation: 1
  - name: Shobeir K. S. Mazinani
    affiliation: "4"
affiliations:
 - name: University of Southern California 
   index: 1
 - name: University of California, Berkeley 
   index: 2
 - name: Lawrence Berkeley National Laboratory
   index: 3
 - name: Aetna Inc.
   index: 4
date:  31 July 2020
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

`AstroPaint` is a python package for simulating and visualizing
 mock maps of astrophysical signals. The
  code takes in a halo catalog and the radial profile of an astrophysical
   signal from the user and then combines the two to paint a whole-sky mock
    map of
    the observable at
    high resolution (see the workflow \autoref{sec:workflow} for details
    ). `AstroPaint` also contains a
     variety of methods that
     facilitate analysis routines such as data filtering, map manipulation
     , and cutout stacking. The package has an Object-Oriented structure and
      runs in parallel, making it scalable for production of high resolution
       maps with large underlying catalogs.  
          
     

![Map of the Birkinshaw-Gull effect painted with astropaint on top of the
 WebSky catalog](../images/BG_websky_cover.png)

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
       effect to measure transverse velocities of halos using the Simons
        Observatory or CMB-S4, one needs a mock map of BG 
        as well as maps of potential contaminants such as kSZ or tSZ for the
         same set of objects.  
     
Creating realistic maps of astrophysical effects through
 hydrodynamical simulations can be prohibitive for large numbers of objects
 . An alternative strategy to creating mock observations of extended objects
  such as galaxies and galaxy cluster halos would be to simulate the
   positions of
   these objects (either semi-analytically or through N-body
    simulations)
  and then synthetically paint the desired signal at the location of the
   halos. `AstroPaint` is developed to accomplish this latter step.  
 

 

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