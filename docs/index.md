# Welcome to teglon's documentation

![Map](static/map.png)

Teglon is a pixel based gravitational wave search and analysis pipeline. It optimizes 
EM search strategies for gravitational wave sources given an input galaxy catalog. Teglon can 
also be used to calculate transient lightcurve model detections efficiencies at the pixel level as 
a function of a set of observations. Teglon is instrument agnostics and can model any telescope 
footprint. Teglon is easily extensible to contain and query any 
[HEALPix](https://healpix.jpl.nasa.gov/) based data (e.g., dustmaps, GW localizations and CMB data).

Essential Teglon functions

* Load LIGO HEALPix maps
* Extract observations
* Plot observation plans and LIGO localization probably maps 
* Ingest community EM observations
* Calculate detection efficiency of a Library transient lightcurves given a map and a set of observations
