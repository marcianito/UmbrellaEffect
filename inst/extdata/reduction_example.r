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
# library(devtools)
# setwd("/home/mreich/Dokumente/written/ResearchRepos/")
# load_all("/home/mreich/Dokumente/written/ResearchRepos/UmbrellaEffect")
# create docu
# library(roxygen2)
# setwd("/home/mreich/Dokumente/written/ResearchRepos/UmbrellaEffect")
# devtools::document()
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
library(HyGra)
library(UmbrellaEffect)
library(reshape2)
library(ggplot2)
library(viridis)
library(gstat)
library(ptinpoly)
#
library(plot3D)
message("done.")
####################

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
# relativ, local coordinate sytem
SG_x = 3
SG_y = 3
SG_Z = 0
SG_SensorHeight = 0
# UTM coordinate system
SG_x = 4564041.87 
SG_y = 5445662.88 
SG_Z = 606.471
SG_SensorHeight = 1.05 

## Model domain
# in [m]
# local grid or UTM, depending on the coordinates of the SG !
Building_x = c(0, 6) # min, max
Building_y = c(0, 6) # min, max
# grid3d_depth = c(-3, 0) # min, max
# UTM
Building_x = c(SG_x - 3, SG_x + 3) # min, max
Building_y = c(SG_y - 3, SG_y + 3) # min, max
# grid3d_depth = c(SG_Z, SG_Z - 3) # min, max

## Model discretization
# in [m]
grid3d_discr = data.frame(x = .5, y = .5, z = .5)
grid3d_depth = c(-5, 0) # min, max

## Parameters for foundation of building
# these include baseplate, walls, SG pillar(s)
# please use same units as in DEM and model domain
Building_walls_x = .6 # extension
Building_walls_y = .6 # extension 
Building_walls_z = 1.5 # extension
# local grid
Building_baseplate_z = c(-.5, 0) # min, max
Building_SGpillar_x = c(2, 4) # min, max
Building_SGpillar_y = c(2, 4) # min, max
Building_SGpillar_z = c(-1, 0) # min, max
# UTM
Building_baseplate_z = c(SG_Z - .5, SG_Z) # min, max
Building_SGpillar_x = c(SG_x - 1, SG_x + 1) # min, max
Building_SGpillar_y = c(SG_y - 1, SG_y + 1) # min, max
Building_SGpillar_z = c(SG_Z - 1, SG_Z) # min, max

## SG position
# options are: Center, Corner, Side, Trans, Wettzell
SG_position = "Center"

## Hydrological site condition
# options are: 
# "Soil type sandy loam" , "Soil type silty loam", "Soil type clay loam" ; for soil type scenarios
# "Anisotropy sandy loam", "", "" ; for anisotropy of hydraulic conductivity scenarios
# "Climate AI=0.35", "Climate AI=0.575", "Climate AI=1.0", "Climate AI=2.0" ; for climate scenarios
Hydro_condition = "Soil type sandy loam"

## Vertical extent of reduction
# indicates the vertical extent until which a reduction factor should be calculated and applied
# this is valid for both Hydroglical site condition & SG position reduction as also building dimension reduction factors
# units in [m]
# maximum value: 5
VerticalExt_reduction = 5

## Input files
## general settings
# in case using .csv data, the special information has to supplied, in which columns the spatial information is stored
# the settings below are valid for 1d data files
# in the vector, the order is: x, y, z
spatial_col = c(NA, NA, 2)
# in all cases, a column has to be specified, containing the observation data
# columns of observation data (in .csv and .rData files)
data_col = 3
# columns of observation data (in .tsf files)
data_tsf = 13 
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
# ## test datasets
# # soil moisture
# load(paste0(dir_input, soilMoisture_input_file))
# SoilMoisture_input
# #gravity
# Gravity_input = read_data(
#                     data_in = gravityObservations_input_file,
#                     data_dir = dir_input,
#                     # spat_col = c(NA, NA, 2),
#                     # dat_col = 3,
#                     dat_tsf = data_tsf)
# dev.new()
# plot(Gravity_input)
# ##

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
# Generate 3d gravity component grid 
#########################################
message("Generate 3d gravity component grid..")
## with NEW FUNCTION "create_gravityGrid()"
load_all("/home/mreich/Dokumente/written/ResearchRepos/gravitySynth")
gravity_component_grid3d = create_gravityGrid(
            DEM_input_file = DEM_input_file,
            dir_input_DEM = "",
            SG_coordinates = SGloc,
            grid_discretization = grid3d_discr,
            grid_depth = grid3d_depth,
            range_coords_x = Building_x,
            range_coords_y = Building_y
)

message("done.")
#########################################
## Correct gravity component grid for SG building foundation
#########################################
message("Correct gravity component grid for SG building foundation..")

gravity_component_grid3d = correct_SGbuilding_foundation(
            gravity_gcomp = gravity_component_grid3d,
            Bdwall_ext_x = Building_walls_x,
            Bdwall_ext_y = Building_walls_y,
            Bdwall_ext_z = Building_walls_z,
            Bdbase_x = Building_x,
            Bdbase_y = Building_y,
            Bdbase_z = Building_baseplate_z,
            Pillar_x = Building_SGpillar_x,
            Pillar_y = Building_SGpillar_y,
            Pillar_z = Building_SGpillar_z,
            grid_discretization = grid3d_discr
)

if(plot_data){
  message("Plotting transect of gravity component grid and saving plot to output directory..")
  plot_gcomp_grid(
                  grid_input = gravity_component_grid3d,
                  yloc = SG_y,
                  output_dir = dir_output,
                  grid_discretization = grid3d_discr
)
}

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

## TESTING WITH SOME
# soil moisture and gravity data
SMplot = dplyr::filter(SMgrid3d_outside, Depth == 0) %>%
    dplyr::group_by(x,y) %>%
    ggplot(aes(x = datetime, y = value, colour = paste0(x,"_",y))) + geom_line()

unique(SMgrid3d_outside$Depth)
SMplot = SMgrid3d_outside %>%
    dplyr::group_by(datetime, Depth) %>%
    dplyr::summarize(value = mean(value)) %>%
    ggplot(aes(x = datetime, y = value, colour = Depth)) + geom_line() +
    facet_grid(Depth ~ .)

gcomps_depth = gravity_component_grid3d %>%
    dplyr::group_by(Depth) %>%
    dplyr::summarize(gcomp = sum(gcomp))

gravity_per_depth = SMgrid3d_outside %>%
    dplyr::group_by(datetime, Depth) %>%
    dplyr::summarize(value = mean(value)) %>%
    dplyr::inner_join(gcomps_depth) %>%
    dplyr::mutate(gravity = gcomp * value)

ggplot(gravity_per_depth, aes(x = datetime, y = gravity, colour = Depth)) + geom_line() +
    facet_grid(Depth ~ .)
    
gravity_over_time = gravity_per_depth %>%
  dplyr::group_by(datetime) %>%
  dplyr::summarize(value = sum(gravity))
dev.new()
## end TESTING WITH SOME

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
            setScenario = Hydro_condition,
            setSGlocation = SG_position,
            setVertLimit = VerticalExt_reduction,
            MeanSoilMoisture = SoilMoisture_mean_ts
)

# calculate SG building size
buildingSize = BdSize(
           Bd_x = Building_x,
           Bd_y = Building_y
)

# reduction factor corresponding to size (area) of the SG building
reduction_factor_BdSize = reduction_BdSize(
            SG_BdSize = buildingSize,
            VertLimit = VerticalExt_reduction
)

## convert gravity response from next to building
gravity_response_below_building = convert_gravity_response(
            gravity_input = gravity_response_outside_building,
            factor_hydScen_SGloc = reduction_factor_hydSen_SGloc,
            factor_BdSize = reduction_factor_BdSize
)

g_resp_below_ma = rollapply(gravity_response_below_building$value, 24*7, mean)
dev.new()
plot(g_resp_below_ma)

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
            input_dir = dir_input,
            dat_tsf = data_tsf
  )
  
  if(output_type == "csv"){
      save(gravity_data_reduced, file=paste0(dir_output, "gravity_data_reduced.rData"))
      write.table(gravity_data_reduced, file=paste0(dir_output, "gravity_data_reduced.csv"), row.names = F)
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
            output_dir = dir_output,
            dat_tsf = data_tsf
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
and stored as well in the output directory.")

