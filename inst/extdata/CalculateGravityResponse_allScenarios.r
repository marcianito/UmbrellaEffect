############################################################################
### Calculate gravity signals of umbrella space (data below and outside), on GFZ cluster server ###
############################################################################

#devtools + load package
library(devtools)
load_all("/home/hydro/mreich/R/hygraFuncs")
#library(vimcom)

#libraries
library(zoo); Sys.setenv(TZ = "GMT")
library(xts)
#library(SDMTools)
library(reshape2)
library(ggplot2)
library(dplyr)
#library(dyplrExtras)
library(grid); library(gridExtra)
library(scales)

# overview
# this script does the following:
# calculate gravity response for all grids and sum them to one gmod / timestep

#intro
print("Starting gravity response calculation..")

#start counting for timing of whole script
script_start = proc.time()
#defining directories for storing and loading

#########################
## settings
#########################
# set directories
hydrodir = "~/Scenarios_year2011/"
gravitydir = "~/gravity/"
SG_locations = c("real", "center", "corner", "trans1", "trans2")
#########################

# set all folders and then run the script in a loop to cover all folers with same calculations
folders=c(
    "Bodenart_ClayLoam_shortEdge/",
    "Bodenart_SiltyLoam_shortEdge/",
    "Bodenart_SandyLoam_shortEdge/",
    "Anisotropy_SandyLoam_shortEdge/",
    "Anisotropy_SiltyLoam_shortEdge/",
    "Anisotropy_ClayLoam_shortEdge/",
    "Climate_AI035_shortEdge/",
    "Climate_AI0575_shortEdge/",
    "Climate_AI10_shortEdge/",
    "Climate_AI20_shortEdge/"
	# "RealisticModel_gneiss_layerlimit_sy_shortEdge"
	  )

for(fol in folders){ #start for loop for reading and processing various hydrus scenarios

#########################
## calculate modeled gravity response for SGbuilding domain for one scenarios
#########################
# 1) data SM modeled from below building
# 2) data SM modeled from next to building

# set folder for data from NEXT to building (NObd-scenarios)
fol_NObd = paste0("NObd_", sub("_shortEdge","", fol))
scen = sub("_shortEdge/","", fol)

#########################
## load gravity grids
#########################

# SGbuilding model domain
load(file=paste(gravitydir,"grids/gcomp_SGbuilding_domain_withoutUmbrella.rdata",sep=""))
# hydrus model domain
load(file=paste(gravitydir,"grids/gcomp_hydrus_domain_withoutUmbrella.rdata",sep=""))

#########################
## load soil moisture data
#########################

# below building (3d)
load(file=paste(hydrodir,fol,"/output/SMgrid3d.rdata",sep=""))
print("SMgrid3d loaded")
print(str(SMgrid3d))

# next to building (1d)
load(file=paste(hydrodir,fol_NObd,"/output/SMdata_TS_1Dvertical_interpolated_allLayers_r50m_r300m.rdata",sep=""))
print(str(SMvertical_r50m_r300m))

#########################
## prepare SM data
#########################

# # add UTM coordiinates
# SMgrid3d$yDEM = SMgrid3d$y + min(gcomp_hydrus_domain$y)
# SMgrid3d$xDEM = SMgrid3d$x + min(gcomp_hydrus_domain$x)

# decrease size of SMgrid3d (limit to year 2011)
# SMgrid3d = dplyr::filter(SMgrid3d, Timestep >= 17520)
# free memory !?
# gc()

# filter for values only beneath SG builidng
hydrusx_min = 5
hydrusx_max = 16
hydrusy_min = 6
hydrusy_max = 14
SMgrid3d_bd = dplyr::select(SMgrid3d, x, y, Depth, Timestep, value) %>%
           dplyr::filter(y <= hydrusy_max & y >= hydrusy_min) %>%
		   dplyr::filter(x <= hydrusx_max & x >= hydrusx_min)
print("reduced to building size")
# free memory !?
rm(SMgrid3d)
gc()
SMgrid3d = dplyr::mutate(SMgrid3d_bd, xDEM = x + min(gcomp_hydrus_domain$x)) %>%
           dplyr::mutate(yDEM = y + min(gcomp_hydrus_domain$y)) %>%
           dplyr::select(xDEM, yDEM, Depth, Timestep, value)
print("added UTM coords")
# free memory !?
rm(SMgrid3d_bd)
gc()
# round data to make use join is performed corretely !
SMgrid3d$Depth = round(SMgrid3d$Depth, 1)
colnames(SMgrid3d)[1:3] = c("x","y","zgrid")
# colnames(SMgrid3d)[1:2] = c("x","y")

# # limit SM data to space below building
# SMgrid3d = left_join(dplyr::select(gcomp_SGbuilding_domain, x, y), SMgrid3d)
print(str(SMgrid3d))

SMmod_NObd = dplyr::filter(SMvertical_r50m_r300m, Timestep >= 17520)
# rename value column in order to merge with SMgrid3d
# colnames(SMmod_NObd)[3] = "value_NObd"
# round data to make use join is performed corretely !
SMmod_NObd$Depth = round(SMmod_NObd$Depth, 1)
colnames(SMmod_NObd)[2] = "zgrid"
print(str(SMmod_NObd))
# create SM 3d grid with data from NEXT to building
SMgrid3d_outside = dplyr::left_join(dplyr::select(SMgrid3d, -value), SMmod_NObd)# %>%
                   # dplyr::mutate(value = value_NObd) %>%
                   # dplyr::select(-value_NObd)
print(str(SMgrid3d_outside))

#########################
## start for loop for different SG locations within building
#########################

## exclude building basement, walls and SG pillar
for(SGloc in SG_locations){ 
switch(SGloc,
  real = { # SG pillar at real pillar location inside of building
         load(file=paste(gravitydir,"grids/SGbuilding_basement_walls_SGreal_coords.rdata",sep=""))
         load(file=paste(gravitydir,"grids/gcomp_SGbuilding_domain_withoutUmbrella_SGreal_rel.rdata",sep=""))
         gcomp_SGbuilding_domain_SGlocation = gcomp_SGbuilding_domain_SGreal_rel
         },
  center = { # SG pillar at real pillar location inside of building
         load(file=paste(gravitydir,"grids/SGbuilding_basement_walls_SGcenter_coords.rdata",sep=""))
         load(file=paste(gravitydir,"grids/gcomp_SGbuilding_domain_withoutUmbrella_SGcenter_rel.rdata",sep=""))
         gcomp_SGbuilding_domain_SGlocation = gcomp_SGbuilding_domain_SGcenter_rel
         },
  corner = { # SG pillar at real pillar location inside of building
         load(file=paste(gravitydir,"grids/SGbuilding_basement_walls_SGcorner_coords.rdata",sep=""))
         load(file=paste(gravitydir,"grids/gcomp_SGbuilding_domain_withoutUmbrella_SGcorner_rel.rdata",sep=""))
         gcomp_SGbuilding_domain_SGlocation = gcomp_SGbuilding_domain_SGcorner_rel
         },
  trans1 = { # SG pillar at real pillar location inside of building
         load(file=paste(gravitydir,"grids/SGbuilding_basement_walls_SGtrans1_coords.rdata",sep=""))
         load(file=paste(gravitydir,"grids/gcomp_SGbuilding_domain_withoutUmbrella_SGtrans1_rel.rdata",sep=""))
         gcomp_SGbuilding_domain_SGlocation = gcomp_SGbuilding_domain_SGtrans1_rel
         },
  trans2 = { # SG pillar at real pillar location inside of building
         load(file=paste(gravitydir,"grids/SGbuilding_basement_walls_SGtrans2_coords.rdata",sep=""))
         load(file=paste(gravitydir,"grids/gcomp_SGbuilding_domain_withoutUmbrella_SGtrans2_rel.rdata",sep=""))
         gcomp_SGbuilding_domain_SGlocation = gcomp_SGbuilding_domain_SGtrans2_rel
         }
) # end switch

#########################
# adjust gcomp grid to area below SGbuilding
#########################

# add UTM coordiinates to gcomp grid 
gcomp_SGbuilding_domain_SGlocation$y = gcomp_SGbuilding_domain_SGlocation$y + min(gcomp_SGbuilding_domain$y)
gcomp_SGbuilding_domain_SGlocation$x = gcomp_SGbuilding_domain_SGlocation$x + min(gcomp_SGbuilding_domain$x)

## removing SGPILLAR, WALLS, BASEMENT of building
gcomp_SGbuilding_domain_umbrella = left_join(gcomp_SGbuilding_domain_SGlocation, select(SGhouse_grid,-xrel,-yrel,-Depth), by=c("x","y","zgrid")) %>%
			mutate(gcomp = ifelse(is.na(house) == T, gcomp, 0)) %>%
			select(-house)

# remove duplicated entries
dups = which(duplicated(gcomp_SGbuilding_domain_umbrella[,1:3]) == T)
gcomp_SGbuilding_domain_umbrella = gcomp_SGbuilding_domain_umbrella[-dups,]

# show menory sizes of objects
#print("AFTER gmod r50m calculations")
#sort( sapply(mget(ls()),object.size) )
#print("system memory status:")
#print(system("free -m"))

#print(str(gcomp_SGbuilding_domain_umbrella))
#print(str(SMgrid3d))

##################################
## calculate gravity: use data from below building
##################################

print("calucaling gmod: SGbuilding domain, data below...")
gmod_SGbuilding_domain_umbrella = gsignal_grids_3d(gcomp_SGbuilding_domain_umbrella, SMgrid3d, F)
# save
save(gmod_SGbuilding_domain_umbrella, file=paste(gravitydir,"output/gmod_SGbuilding_domain_umbrella_",SGloc, "_", scen, ".rdata",sep=""))

##################################
## calculate gravity: use data from next to building
##################################

print("calucaling gmod: SGbuilding, data outside domain...")
gmod_SGbuilding_domain_umbrella_dataOutside = gsignal_grids_3d(gcomp_SGbuilding_domain_umbrella, SMgrid3d_outside, F)
# save
save(gmod_SGbuilding_domain_umbrella_dataOutside, file=paste(gravitydir,"output/gmod_SGbuilding_domain_umbrella_",SGloc, "_", scen, "_dataOutside.rdata",sep=""))

## clean up memory
rm(gcomp_SGbuilding_domain_umbrella); gc()

## end for loop for SG locations
}

## clean up memory
rm(SMgrid3d, SMgrid3d_outside); gc()

## end for folder
print(paste("finished folder: ",fol,sep=""))

## end for loop for different folders
} 

print("FINISHED")
#end couting runtime
script_end = proc.time()
script_end - script_start #total runtime in seconds

