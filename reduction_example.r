########################################################
### example script: reduction due to umbrella effect ###
########################################################

## please modify or (better) copy this file to start a new reduction routine


# load package
# library(HyGra)
library(devtools)
load_all("/home/mreich/HyGra")

## load libraries
library(zoo); Sys.setenv(TZ = "GMT")
library(xts)
library(raster)
library(reshape2)
library(ggplot2)
library(dplyr)
# library(dyplrExtras)
library(grid)
library(gridExtra)
library(scales)
# library(tables)
library(RColorBrewer)
# library(xtable)
library(viridis)
# library(R.matlab)
# library(h5)
# library(RcppOctave)





#########################################
## SETUP
####################

## Working directores
dir_input = "./data/"
dir_output = "XX/YY/"
setwd(dir_data)

## Gravimeter location
SG_x = 4564082.00
SG_y = 5445669.70
SG_Z = 609.755
SG_SensorHeight = 1.05
igrav_z = igrav_dem + igrav_sensor
igravLoc = data.frame(x=igrav_x, y=igrav_y, z=igrav_z)
locations = data.frame(x=igrav_x, y=igrav_y, name="iGrav")

## DEM input file
# complete path (or absolute or relative to input directory)
DEM_input_file = "WE_UP_TO_300m_05m.asc"

## Mdel domain
# in [m]
# local grid or UTM, depending on the coordinates of the SG !
Building_x = c() # min, max
Building_y = c() # min, max
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
# options are: center, corner, side, trans1
SG_position = "center"

## Hydrological site condition
# options are: standard, Ksat anisotropy, Climate [dry, normal, wet], ...
Hydro_condition = "standard"

## end SETUP
#########################################

#########################################
## CALCULATIONS
####################



#########################################
## Gravimeter location
#########################################
SG_x = 4564082.00
SG_y = 5445669.70
SG_Z = 609.755
SG_SensorHeight = 1.05
igrav_z = igrav_dem + igrav_sensor
igravLoc = data.frame(x=igrav_x, y=igrav_y, z=igrav_z)
locations = data.frame(x=igrav_x, y=igrav_y, name="iGrav")

#########################################
## read DEM data & plot DEM 
#########################################
## read DEM from ascii
dempath = ""
filename = "WE_UP_TO_300m_05m.asc" #r=300m, dxdy=0,5m
dem_raster = raster(paste(dempath,filename, sep=""))

#########################################
## generate grids for gravity component calculations
#########################################
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
save(surface_grid, file="surface_grid.rdata")

#########################################
# generate 3d grid with individual discretization
#########################################
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
gcomp_grid_igrav = gcomp_raw_Edges(grid3d, igravLoc, grid3d_discr, grid3d_edges)
save(gcomp_grid_igrav, file="gcomp_grid_igrav.rdata")

#########################################
# generate coordinates of SG building (including baseplate, walls & SG pillar)
#########################################

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


##!
HERE NOW RE-WRITTEN SCRIPTS FOR "REAL" CALCULATIONS
INCLUDE (and generate) CORRECTION PARAMETER DATA !!
##!


## end CALCULATIONS
#########################################

print("ALL calculations have finished.")
print("Please have a look at the output file, located at: ")
print(dir_output)

print("If direct reduction of a observed gravity data time series is desired,
      please run script "XX.r" and provide the necessary gravity input data")

