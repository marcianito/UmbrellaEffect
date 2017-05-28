Repository containing scripts related to _Reich et al., (2017, GWR)_ publication
======================================================================
**Article:**  
Marvin Reich<sup>1</sup>, Michal Mikolaj<sup>1</sup>, Theresa Blume<sup>1</sup>, Andreas Güntner<sup>1,4</sup>: **_TITLE_**, published in: [GWR](http://link.com), **2017**

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
2. Load library "devtools"
3. Install package "UmbrellaEffect"
4. load library "UmbrellaEffect"

## Dependencies

### Computationally
* r-base version 3.3.1
* following R-packages: ....name packages...

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
(Steps 2 to X are only explanatory for what the script does; no modifications necessary.)

1. Setup: WHAT IS DOES, WHAT YOU NEED
	* _Set output directory (and file extention !? .csv / .rData !?)_
	* _Set SG coordinates_
	* _Set correct file to load for DEM input_
	* _Choose discretization_
	* _Construct foundation of building_
	* _Set SG location variables (center, corner, etc.)_
	* _Choose dominant hydrological site condition_
	* __
2. Calculate gravity response time series of observation data: WHAT IS DOES 
	* __
3. Calculate corresponding gravity response below the building : WHAT IS DOES 
	* __
X. : WHAT IS DOES 
	* __
5. Run the entire script and look at output files
