Repository containing scripts related to _Reich et al., (2017, ???)_ publication
======================================================================
**Article:**  
Marvin Reich<sup>1</sup>, Michal Mikolaj<sup>1</sup>, Theresa Blume<sup>1</sup>, Andreas Güntner<sup>1,4</sup>: **_TITLE_**, published in: [???](http://link.com), **2017**

**Affiliation**  
<sup>1</sup>Helmholtz Centre Potsdam GFZ German Research Centre for Geosciences, Section 5.4 Hydrology, 14473 Potsdam, Germany  
<sup>2</sup>University of Potsdam, Institute of Earth and Environmental Science, 14476 Potsdam, Germany
> _Please cite this article when using here provided scripts_

## Description

This is a small R-package which aims at supporting the presented gravity data reduction method,
presented in the above mentioned paper.
With this package you can easily walk through the necessary steps in order to reduce your local
gravity observation data for the mass changes occurring below the SG building.
The procedure within the package is identical with the one discussed and presented in the publication.
Consequentially, all the needed parameters are included in in-package data.

For bug fixes, comments or further development please contact: mreich@posteo.de.

## Installation

1. Start R
2. Install package via devtools: 
`devtools::install_github("marcianito/UmbrellaEffect")`

4. load package: 
`library(UmbrellaEffect)`

## Dependencies

### Computationally
* r-base version 3.3.1
* following R-packages: devtools, dplyr, ggplot2, gstat, ptinpoly, raster, reshape2, viridis, xts, zoo
* system libraries for rgdal and devtools

in debian install using: 
`apt-get update && apt-get install libgdal-dev libproj-dev libssl-dev`

### Data-wise
It is necessary to have modeled or measured soil moisture data from directly next to the SG building.
This data should consists at least of one vertical soil moisture profile with sufficient vertical point information
to cover the site specific stratification of soil layers.

## Processing

1. Start R
2. Load reduction_example.r script
3. Modify according to description below
4. Run script and use output time series for reduction

## Reduction procedure (computations)
#### For more details, please look at the vignette or the corresponding help-files (within R).

All changes should be done in a new file following (or a copy of) the reduction_example.r file.

(Step 2 is only explanatory for what the script does; no modifications necessary.)

1. Setup: 
	* _Directory and configs (input / output, file extentions, enable plotting)_
	* _Gravimeter coordinates (x, y, z + height of sensor)_
	* _Model domain (x and y extensions)_
	* _Discretization and vertical model extent_
	* _Foundation of building (walls, baseplate, gravimeter pillar)_
	* _Gravimeter location within building (Center, Corner, Side, Trans or Wettzell)_
	* _Dominant hydrological site condition (Possible scenarios: Soil types, Ksat anisotrophy, climate region)_
	* _Set correct file to load for DEM input_
	* _Input file settings (data dimensions, order of data columns from input files)_
	* _Input file names (DEM, soil moisture data, observed gravity signal)_
2. Calculations: 
	* _Construct surface grid_
	* _Create gravity component grid_
	* _Correct gravity component grid for observatory building foundation stucture_
	* _Plot 2d segment of created 3d grid_
	* _Extrapolate soil moisture from input data to 3d grid_
	* _Calculate gravity response for area outside of observatory building_
	* _Calculate mean soil moisture at every timestep_
	* _Set reduction factor based on hydrological scenario and gravimeter location_
	* _Set reduction factor based footprint of observatory building_
	* _Convert gravity response outside of building to gravity response below building_
	* _Reduce observed gravity signal with this gravity response below buidling_
	* _Plot all time series_
5. Run the entire script and look at output files

