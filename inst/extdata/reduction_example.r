#!/usr/bin/Rscript

########################################################
### example script: reduction due to umbrella effect ###
########################################################

####################
## Instructions
####################
# Please modify or (better) copy this file to start a new reduction routine
# All lines in the SETUP section have to be filled, according to SG setup
# On completion, the script can be sourced
# or independently run within an R console
# Once it is finished, output files are stored in the set folders
##
# Optional: If a gravity data observation time series was supplied,
# data will be reduced by output results (directly by this script)
##
# For questions, comments or bugs,
# please visit https://github.com/marcianito/UmbrellaEffect
# or write to mreich@posteo.de
####################

## developing
# library(HyGra)
library(devtools)
library(roxygen2)
setwd("/home/mreich/Dokumente/written/ResearchRepos/")
load_all("UmbrellaEffect")
# create docu
setwd("/home/mreich/Dokumente/written/ResearchRepos/UmbrellaEffect")
document()
# install package locally
# not working !? -> use load_all() above
# install()

####################
## load libraries
message("Loading necessary libraries..")
library(zoo); Sys.setenv(TZ = "GMT")
library(xts)
library(dplyr)
library(raster)
# library(UmbrellaEffect)
# library(reshape2)
library(ggplot2)
library(viridis)
library(gstat)
message("done.")
####################

# # library(dyplrExtras)
# library(grid)
# library(gridExtra)
# library(scales)
# # library(tables)
# library(RColorBrewer)
# # library(xtable)

#########################################
## SETUP
####################
message("Initializing setup..checking input data..")

## Output and input settings
# Directory
# path should be absolute
# (if not, it will be relative to the packages library path)
dir_input = "~/temp/UM/"
dir_output = "~/temp/UM/"
# Output file type
# set to "csv", if output should also be saved as .csv (besides .rData)
output_type = ""
# Plotting option: should a plot of all time series be shown (and saved) in the end?
plot_data = TRUE

## Gravimeter location
# in [m]
SG_x = 3
SG_y = 4
SG_Z = 0
SG_SensorHeight = 1.5 

## Model domain
# in [m]
# local grid or UTM, depending on the coordinates of the SG !
Building_x = c(0, 6) # min, max
Building_y = c(0, 6) # min, max
grid3d_depth = c(0, 2) # min, max

## Model discretization
# in [m]
grid3d_discr = data.frame(x = 1, y = 1, z = 1)

## Parameters for foundation of building
# these include baseplate, walls, SG pillar(s)
# please use same units as in DEM and model domain
Building_walls_x = .5 # extension
Building_walls_y = .2 # extension 
Building_walls_z = 2 # extension
Building_baseplate_z = c(0, .5) # min, max
Building_SGpillar_x = c(2, 4) # min, max
Building_SGpillar_y = c(3, 5) # min, max
Building_SGpillar_z = c(0, 1) # min, max

## SG position
# options are: Center, Corner, Side, Trans, Wettzell
SG_position = "Center"

## Hydrological site condition
# options are: standard, Ksat anisotropy, Climate [dry, normal, wet], ...
Hydro_condition = "Soil type sandy loam"

## Input files
## general settings
# in case using .csv data, the special information has to supplied, in which columns the spatial information is stored
# the settings below are valid for 1d data files
# in the vector, the order is: x, y, z
spatial_col = c(NA, NA, 2)
# in all cases, a column has to be specified, containing the observation data
# columns of observation data
data_col = 3
# if the .csv has special characters for separating columns
# or decimal places, etc.
# the have to be EXPLICITLY specified in the read_data-function
# using sep = "??"
# using dec = "??"
# for further usage see ?read.csv

## DEM input file
# file name including its path
# should be absolute
# if left empty, a flat topographie will be assumed
DEM_input_file = ""
DEM_input_file = "WE_UP_TO_300m_05m.asc"

## Soil moisture data time series (observed or modelled)
soilMoisture_input_file = "SMdata_TS_1d.rData"

## Observed gravity data time series
# this is optional and can be left empty if no automatized reduction is desired
gravityObservations_input_file = "SG030_TS_1month.tsf"

message("done.")
## end SETUP
#########################################

#########################################
## CALCULATIONS
####################
## nothing has to be changed from here on !!
message("Starting with calculation routine..")

# set working directory
dir_wd = system.file("data", package="UmbrellaEffect")
setwd(dir_wd)

#########################################
## Gravimeter location
#########################################
SG_z = SG_Z + SG_SensorHeight
SGloc = data.frame(x=SG_x, y=SG_y, z=SG_z)

#########################################
## Generate cropped DEM and surface grid
#########################################
message("Generate cropped DEM and surface grid..")

surface_grid = surface_grid(
            DEM = DEM_input_file,
            grid_domain_x = Building_x,
            grid_domain_y = Building_y,
            input_dir = dir_input,
            output_dir = dir_output
            # , sep = "a", etc.
)

message("done.")
#########################################
# Generate 3d gravity component grid 
#########################################
message("Generate 3d gravity component grid..")

## create 3d grid of gravity components
gravity_component_grid3d = gravity_comp_grid(
            surface = surface_grid,
            SG_coordinates = SGloc,
            grid_discretization = grid3d_discr,
            grid_depth = grid3d_depth
)

message("done.")
#########################################
# Generate foundation of SG building (baseplate, walls & SG pillar)
#########################################
message("Generate foundation of SG building (baseplate, walls & SG pillar)..")

SGbuilding_foundation = building_foundation(
            Bdwall_ext_x = Building_walls_x,
            Bdwall_ext_y = Building_walls_y,
            Bdwall_ext_z = Building_walls_z,
            Bdbase_x = Building_x,
            Bdbase_y = Building_y,
            Bdbase_z = Building_baseplate_z,
            Pillar_x = Building_SGpillar_x,
            Pillar_y = Building_SGpillar_y,
            Pillar_z = Building_SGpillar_z,
            grid_discretizaion = grid3d_discr
)

message("done.")
#########################################
## Correct gravity component grid for SG building foundation
#########################################
message("Correct gravity component grid for SG building foundation..")

gravity_component_grid3d = correct_SGbuilding_foundation(
            gravity_gcomp = gravity_component_grid3d,
            building_foundation = SGbuilding_foundation
)

message("done.")
#########################################
## Extrapolate soil moisture time series data (observed or modelled) to gravity grid domain
#########################################
message("Extrapolate soil moisture time series data (observed or modelled) to gravity grid domain..")

SMgrid3d_outside = SoilMoisture_grid3d(
            grid_domain = gravity_component_grid3d,
            soilMoisture_input = soilMoisture_input_file,
            grid_discretization = grid3d_discr,
            grid_depth = grid3d_depth,
            input_dir = dir_input
            # , sep = "a", etc..
)

message("done.")
#########################################
## Calculate gravity response (from outside of building)
#########################################
message("Calculate gravity response (from outside of building)..")

gravity_response_outside_building = calculate_gravity_response(
            gcomp_grid = gravity_component_grid3d,
            mass_input = SMgrid3d_outside
)

message("done.")
#########################################
## Calculate mean soil moisture within model domain for each timestep
#########################################
message("Calculate mean soil moisture within model domain for each timestep..")

SoilMoisture_mean_ts = dplyr::group_by(SMgrid3d_outside, datetime) %>%
                           dplyr::summarize(value = mean(value, na.rm=T))

message("done.")
#########################################
## Convert gravity response (outside) to gravity response below SG building
#########################################
message("Convert gravity response (outside) to gravity response below SG building..")

# reduction factor corresponding to chosen dominant hydrological scenario and SG position
reduction_factor_hydSen_SGloc = reduction_hydScen_SGloc(
            Scenario == Hydro_condition,
            SGlocation == SG_position,
            MeanSoilMoisture = SoilMoisture_mean_ts)

# calculate SG building size
buildingSize = BdSize(
           Bd_x = Building_x,
           Bd_y = Building_y
)

# reduction factor corresponding to size (area) of the SG building
reduction_factor_BdSize = reduction_BdSize(
            SG_BdSize = buildingSize
)

## convert gravity response from next to building
gravity_response_below_building = convert_gravity_response(
            gravity_input = gravity_response_outside_building,
            factor_hydScen_SGloc = reduction_factor_hydSen_SGloc,
            factor_BdSize = reduction_factor_BdSize
            )

message("done.")
#########################################
## Save gravity response of mass variations, occuring below SG building
#########################################
message("Save gravity response of mass variations, occuring below SG building..")

if(output_type == "csv"){
    save(gravity_response_below_building, file=paste0(dir_output, "gravity_response_below_building.rData"))
    write.table(...)
}else{
    save(gravity_response_below_building, file=paste0(dir_output, "gravity_response_below_building.rData"))
}

message("done.")
#########################################
## Correct gravity observation data (if supplied)
#########################################

if(gravityObservations_input_file == ""){
  message("No reduction of observed gravity data desired.")
}else{
  message("Reducing gravity observation data..")
  
  gravity_data_reduced = reduce_gravity(
            gravity_obs = gravityObservations_input_file,
            gravity_below = gravity_response_below_building,
            input_dir = dir_input
  )
  
  if(output_type == "csv"){
      save(gravity_data_reduced, file=paste0(dir_output, "gravity_data_reduced.rData"))
      write.table(...)
  }else{
      save(gravity_data_reduced, file=paste0(dir_output, "gravity_data_reduced.rData"))
  }
  
  message("done.")
}

#########################################
## Plot all time series
#########################################

if(plot_data){
  message("Plotting time series and saving plot to output directory..")
  if(gravityObservations_input_file == ""){
    plot_ts_data(
            gravity_outside = gravity_response_outside_building,
            gravity_below = gravity_response_below_building,
            input_dir = dir_input,
            output_dir = dir_output
  )
  }else{
    plot_ts_data(
            gravity_obs = gravityObservations_input_file,
            gravity_outside = gravity_response_outside_building,
            gravity_below = gravity_response_below_building,
            gravity_reduced = gravity_data_reduced,
            input_dir = dir_input,
            output_dir = dir_output
  )
  }

  message("done.")
}else{
  message("No plotting desired.")
}

## end CALCULATIONS
#########################################

message("ALL calculations have finished.")
message("Please have a look at the output file, located at: ")
message(dir_output)

message("If gravity observation data was supplied, the data has been recuded automatically by the UmbrellaEffect results,
        and stored as well in the output directory")


##########
### prepare input TS data

# sm = read_data(soilMoisture_input_file, dir_input)
# sm_ts = unique(sm$datetime)
# gg = read_data(gravityObservations_input_file, dir_input)
# gg_ts = unique(gg$datatime)
# sm = sm %>% dplyr::filter(datetime != 0)
# change_dates = data.frame(
#                           datetime = sm_ts[2:745],
#                           datatime = gg_ts)
# sm_mod = left_join(change_dates, sm)
# SoilMoisture_input_1d = data.frame(
#                                    datetime = sm_mod$datatime,
#                                    z = sm_mod$z,
#                                    data = sm_mod$value)
# write.table(SoilMoisture_input_1d, file="SMdata_TS_1d.csv", row.names=F)
# save(SoilMoisture_input_1d, file="SMdata_TS_1d.rdata")
# getwd()

