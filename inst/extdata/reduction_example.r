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
library(UmbrellaEffect)
message("done.")
####################

library(reshape2)
library(ggplot2)
library(viridis)
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
dir_input = ""
dir_output = "XX/YY/"
# Output file type
# set to "csv", if output should also be saved as .csv (besides .rData)
output_type = ""
# Plotting option: should a plot of all time series be shown (and saved) in the end?
plot_data = TRUE

## Gravimeter location
# in [m]
SG_x = 4564082.00
SG_y = 5445669.70
SG_Z = 609.755
SG_SensorHeight = 1.05 

## DEM input file
# file name including its path
# should be absolute
# if left empty, a flat topographie will be assumed
DEM_input_file = "WE_UP_TO_300m_05m.asc"

## Model domain
# in [m]
# local grid or UTM, depending on the coordinates of the SG !
Building_x = c(0,10) # min, max
Building_y = c(0,10) # min, max
grid3d_depth = c(0,3) # min, max

## Model discretization
# in [m]
grid3d_discr = data.frame(x=0.1,y=0.1,z=0.1)

## Parameters for foundation of building
# these include baseplate, walls, SG pillar(s)
Building_walls_x = c() # min, max
Building_walls_y = c() # min, max
Building_walls_z = c() # min, max
Building_baseplate_z = c() # min, max
Building_SGpillar_x = c() # min, max
Building_SGpillar_y = c() # min, max
Building_SGpillar_z = c() # min, max

## SG position
# options are: Center, Corner, Side, Trans, Wettzell
SG_position = "Center"

## Hydrological site condition
# options are: standard, Ksat anisotropy, Climate [dry, normal, wet], ...
Hydro_condition = "Soil texture sandy loam"

## Soil moisture data time series (observed or modelled)
soilMoisture_input_file = ""

## Observed gravity data time series
# this is optional and can be left empty if no automatized reduction is desired
gravityObservations_input_file = ""

message("done.")
## end SETUP
#########################################

#########################################
## CALCULATIONS
####################
## nothing has to be changed from here on !!
message("Starting with calculation routine..")

# set working directory
dir_wd = paste0(system.file(package="UmbrellaEffect"), "/data/")
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
## read DEM from ascii
dempath = ""
filename = "WE_UP_TO_300m_05m.asc" #r=300m, dxdy=0,5m
dem_raster = raster(paste(dempath,filename, sep=""))
# 15m x 15m = 225m²
grid_domain_distances = data.frame(x=c(-7.5,-7.5,7.5,7.5) , y=c(-7.5,7.5,7.5,7.-5))
grid_domain = data.frame(x = grid_domain_distances$x + igrav_x, y = grid_domain_distances$y + igrav_y)
extend_grid_domain = extent(grid_domain)
dem_grid_domain = crop(dem_raster, extend_grid_domain)
writeRaster(dem_grid_domain, "dem_grid", format="ascii", NAflag=-9999, overwrite=T)
# reload grid for verification
filename = "dem_grid.asc" 
read_dem(dempath,filename)
## plot DEM 
plot_dem(dem, dem.info, locations)
# convert to data.frame
surface_grid = convert_demtodf(dem, dem.info)
# plot grid
ggplot(surface_grid, aes(x=x,y=y)) + geom_tile(aes(fill=z)) + geom_point(data=igravLoc, aes(x=x,y=y))
# save


surface_grid = surface_grid(
            DEM = DEM_input_file,
            grid_domain_x = Building_x,
            grid_domain_y = Building_y,
            input_dir = dir_input
)

message("done.")
#########################################
# Generate 3d gravity component grid 
#########################################
message("Generate 3d gravity component grid..")
# load surface grid:
load(file="surface_grid.rdata")
# define 3d grid discretizations
# in [m]
grid3d_discr = data.frame(x=0.1,y=0.1,z=0.1)
grid3d_depth = c(0,3)
# set edge type
# options: first, last, both
# this defines how gravity components are computed at the upper and lower grid boundaries
grid3d_edges = "both"

# convert surface grid to 3d grid
grid3d = demgrid_to_gcompgrid_Edges(surface_grid, grid3d_discr, grid3d_depth, T, surface_grid)

## create 3d grid of gravity components
gravity_component_grid3d = gcomp_raw_Edges(
gravity_component_grid3d = gravity_comp_grid(
            surface = surface_grid,
            # grid_domain = grid3d,
            SG_coordinates = SGloc,
            grid_discretizaion = grid3d_discr,
            grid_depth = grid3d_depth
            #             grid3d_edges
)

message("done.")
#########################################
# Generate foundation of SG building (baseplate, walls & SG pillar)
#########################################
message("Generate foundation of SG building (baseplate, walls & SG pillar)..")

# building walls
x.walls = seq(5,16,by=.1)
yx.walls = c(seq(5,5.5,by=.1),seq(15.5,16,by=.1))
y.walls = seq(5,13, by=.1)
xy.walls = c(seq(5,5.5,by=.1),seq(12.5,13,by=.1))
z.walls = seq(0,1.5, by=.1)
walls.x = expand.grid(xrel=x.walls, yrel=xy.walls, zgrid=z.walls)
walls.y = expand.grid(xrel=yx.walls, yrel=y.walls, zgrid=z.walls)
# building baseplat
x.base = seq(5.6,15.4, by=.1)
y.base = seq(5.6,12.4, by=.1)
z.base = seq(0,.25, by=.1)
base = expand.grid(xrel=x.base, yrel=y.base, zgrid=z.base)
# SG pillar (in center of building)
## ! should add relative coordinates of SG real
## !
x.SGpillar = seq(10,11.1, by=.1)
y.SGpillar = seq(8.5,9.6, by=.1)
z.SGpillar = seq(0,2, by=.1)
SGpillar = expand.grid(xrel=x.SGpillar, yrel=y.SGpillar, zgrid=z.SGpillar)

# combine all parts
SGhouse_grid = cbind(rbind(walls.x, walls.y, base, SGpillar), house= TRUE)

# add UTM coordiinates
SGhouse_grid$x = SGhouse_grid$xrel + min(gcomp_grid_igrav$x)
SGhouse_grid$y = SGhouse_grid$yrel + min(gcomp_grid_igrav$y)

# save building correction grid
save(SGhouse_grid, file="SGbuilding_basement_walls_SGpillar_coords.rdata")

## visual check
ggplot(SGhouse_grid, aes(x=xrel, y=yrel)) + geom_point(aes(fill=house))



SGbuilding_foundation = building_foundation(
            Bdwall_x = Building_walls_x,
            Bdwall_y = Building_walls_y,
            Bdwall_z = Building_walls_z,
            Bdbase = Building_baseplate_z,
            Pillar_x = Building_SGpillar_x,
            Pillar_y = Building_SGpillar_y,
            Pillar_z = Building_SGpillar_z
)


message("done.")
#########################################
## Correct gravity component grid for SG building foundation
#########################################
message("Correct gravity component grid for SG building foundation..")

## removing SGPILLAR, WALLS, BASEMENT of building
gravity_component_grid3d = left_join(gravity_component_grid3d, select(SGbuilding_foundation,-xrel,-yrel,-Depth), by=c("x","y","zgrid")) %>%
			mutate(gcomp = ifelse(is.na(house) == T, gcomp, 0)) %>%
			select(-house)

# remove duplicated entries
dups = which(duplicated(gravity_component_grid3d[,1:3]) == T)
gravity_component_grid3d = gravity_component_grid3d[-dups,]

gravity_component_grid3d = correct_SGbuilding_foundation(
            gravity_gcomp = gravity_component_grid3d,
            building_foundation = SGbuilding_foundation
)

message("done.")
#########################################
## Extrapolate soil moisture time series data (observed or modelled) to gravity grid domain
#########################################
message("Extrapolate soil moisture time series data (observed or modelled) to gravity grid domain..")

# round data to make use join is performed corretely !
SMmod_NObd$Depth = round(SMmod_NObd$Depth, 1)
colnames(SMmod_NObd)[2] = "zgrid"
# create SM 3d grid with data from NEXT to building
SMgrid3d_outside = dplyr::left_join(dplyr::select(SMgrid3d, -value), SMmod_NObd)# %>%

SMgrid3d_outside = SoilMoisture_grid3d(
            grid_domain = gravity_component_grid3d,
            soilMoisture_input = soilMoisture_input_file,
            input_dir = dir_input
)

message("done.")
#########################################
## Calculate gravity response (from outside of building)
#########################################
message("Calculate gravity response (from outside of building)..")

# gravity_response_outside_building = gsignal_grids_3d(gravity_component_grid3d, SMgrid3d_outside, F)
gravity_response_outside_building = calculate_gravity_response(
            gravity_comp_grid = gravity_component_grid3d,
            mass_input = SMgrid3d_outside,
            ? = F
)

message("done.")
#########################################
## Calculate mean soil moisture within model domain for each timestep
#########################################
message("Calculate mean soil moisture within model domain for each timestep..")

SoilMoisture_mean_ts = dplyr::group_by(SMgrid3d_outside, Timestep) %>%
                           dplyr::summarize(water = mean(value, na.rm=T))

message("done.")
#########################################
## Convert gravity response (outside) to gravity response below SG building
#########################################
message("Convert gravity response (outside) to gravity response below SG building..")

## for hydrological scenario and SG location
# load reduction parameters
load(file="reduction_parameters.rData")
reduction_parameter_hydScen_SGloc = dplyr::filter(reduction_parameters,
                                                  Scenario == Hydro_condition,
                                                  SGlocation == SG_position)

reduction_factor_hydSen_SGloc = reduction_hydScen_SGloc(
            parameterset = reduction_parameter_hydScen_SGloc,
            MeanSoilMoisture = SoilMoisture_mean_ts)

## for hydrological scenario and SG location
# load reduction parameters
load(file="reduction_parameters_building.rData")

# calculate SG building size
buildingSize = BdSize(
            Bd_x = Building_x,
            Bd_y = Building_y
)

reduction_parameter_BdSize = dplyr::filter(reduction_parameters,
            BuildingDimension == buildingSize)

reduction_factor_BdSize = reduction_BdSize(
            parameterset = reduction_parameter_BdSize,
            MeanSoilMoisture = SoilMoisture_mean_ts)

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
            gravity_data = gravityObservations_input_file,
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

  plot_ts_data(
            gravity_obs = gravityObservations_input_file,
            gravity_outside = gravity_response_outside_building,
            gravity_below = gravity_response_below_building,
            gravity_reduced = gravity_data_reduced,
            input_dir = dir_input,
            output_dir = dir_output
  )

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



#' @title test
#'
#' @description test
#'
#' @param test
#' @param test
#' @param test
#' 
#' @return test
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

