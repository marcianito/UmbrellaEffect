#' @title test
#'
#' @description test
#'
#' @param test test
#' 
#' @return test
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing


#' @title Calculate gravity component, WITHOUT umbrella effect
#'
#' @description Calculates the gravity components (100percent values) for a given DEM using a nested approach
#'
#' @param g_grid grid for which to calculate gravity components. (generally derived from DEM of region in combination with hydrological modeling mesh/grid).
#' @param senloc coordinates of gravity sensor location. data.frame with x,y,z columns.
#' @param g_discr discretization of g_grid in differences in x,y,and z-direction. they have to be uniform steps!
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples example MISSING
#' @export

fill_gcompgrid <- function(g_grid, senloc, g_discr){
#constants
gama <- 6.673e-11 #m³/(kg*s²)
rho <- 1000 #kg/m³
#w <- 1e8 #gravity units µGal
w <- 1e9 #gravity units nm/s²
#radius criteria
#r2exac=2^2 #bei 2.5 x 2.5 cellsizes ~=   50 m
#r2macm=9^2 #bei 2.5 x 2.5 cellsizes ~= 1000 m
r_inner = 50 #[m]
r_outer = 1000 #[m]
#prepare SG-data-file
# gravity component is just a new column on g_grid
# containing the g_comp, conditionally to some causes, calculated using one of the 3 methods
# !! important
# due to maintaining original grid, gravity component cells (volumini) of the first row (surface) and last row (lower boundary),
# have to be adjusted to only HALF SIZE (value), due to cut off
# advantage of this approach:
# maintaining org. grid facilitates trasnformation of (hydrological) model output
# because no grid_transformation is needed
# if NOT, use layer "would be lost"
# due to midpoint of g component cell (volumina)grid
# here the calculations for each line (grid pointgrid)

# set edge to upper AND lower model domain boundary
# this is due to how to use hydrological model output
edge = "both"

layern_max = max(g_grid$layer)
# rowwise
g_grid = g_grid %>%
	rowwise() %>%
	mutate(gcomp = select_gMethod_Edges(x,y,z,senloc,g_discr, edge, layer, layern_max,r_inner,r_outer,gama,rho,w)) %>%
	ungroup()

return(g_grid)
}

#' @title choose method for computing gravity component for a cell of a grid
#'
#' @description After calculating the gravity components (gcomp()), this function can be used to get real values using soil moisture data (theta). The output of the gravity signal is possible in different dimensions. 
#'
#' @param gcompfiles vector containing the names (complete paths or relative) of the gravity component files to be used for calculations. These files are generated using gcomp().
#' @param theta data.frame with column structure $time, $timestep, $theta (value), $layer
#' @param output Defines the output to be a singel value, layer or grid (default is value)
#' ...
#' @details Usually one needs relative gravity signals, therefore the inputed theta timeseries should be acutally delta theta vaules per timestep. The
#' other option is calculating the differences in values per timestep of the output of this function.
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing 
#' @export
#' 
select_gMethod_Edges = function(xloc, yloc, zloc, gloc, gdiscr, edge, layer, layermax, r_inner, r_outer,gama,rho,w){
        #distances
        rad = sqrt((xloc-gloc$x)^2+(yloc-gloc$y)^2+(zloc-gloc$z)^2) #radial distance to SG
        r2=rad^2
        dr2=gdiscr$x^2+gdiscr$y^2+gdiscr$z^2 #radial "size" of DEM / coordinate-data system
        f2=r2/dr2 #abstand zelle-SG / diagonale aktueller berechnungs-quader
        # different methods after the distance from mass to SG
        #if (f2<=r2exac){ #very close to SG
        if (rad <= r_inner){ #very close to SG
		xl=xloc-(0.5*gdiscr$x);xr=xloc+(0.5*gdiscr$x)
		yl=yloc-(0.5*gdiscr$y);yr=yloc+(0.5*gdiscr$y)
		if(layer==1 | layer == layermax){
			if((edge=="first" & layer==1) | (edge=="both" & layer==1)){ # first z-layer
			 Zint =zloc 
			 Zend = zloc-(0.5*gdiscr$z)
			}
			if((edge =="last" & layer==layermax) | (edge=="both" & layer==layermax)){ # last z-layer
			 Zint = zloc+(0.5*gdiscr$z)
			 Zend = zloc
			}
			else{
			 Zint = zloc+(0.5*gdiscr$z)
			 Zend = zloc-(0.5*gdiscr$z)
			}
		}
		else{ # all other z-layers
		Zint = zloc+(0.5*gdiscr$z)
		Zend = zloc-(0.5*gdiscr$z)
		}
           gcomp_cell=forsberg_raw(gama,w,xl,xr,yl,yr,Zint,Zend,gloc$x,gloc$y,gloc$z,rho) #unit depends on w
        }
         #if(f2>r2macm){ #very far from SG
         if(rad >= r_outer){ #very far from SG
		if(layer==1 | layer == layermax){
			if((edge=="first" & layer==1) | (edge=="both" & layer==1)){ # first z-layer
		  	 zdiscr = 0.5*gdiscr$z
		 	 zloc_mid = zloc+0.5*zdiscr
			}
			if((edge =="last" & layer==layermax) | (edge=="both" & layer==layermax)){ # last z-layer
			 zdiscr = 0.5*gdiscr$z
		  	 zloc_mid = zloc-0.5*zdiscr
			}
			else{
		 	zloc_mid = zloc
		 	zdiscr = gdiscr$z
			}
		}
		else{ # all other z-layers
		 zloc_mid = zloc
		 zdiscr = gdiscr$z
		}
           gcomp_cell=pointmass(gama,zloc_mid,gloc$z,gdiscr$x,gdiscr$y,zdiscr,rad,w,rho) #unit depends on w
        }
        #if(f2>r2exac & f2<r2macm){ #in the "middlle"
        if(rad > r_inner & rad < r_outer){ #in the "middlle"
		if(layer==1 | layer == layermax){
			if((edge=="first" & layer==1) | (edge=="both" & layer==1)){ # first z-layer
		  	 zdiscr = 0.5*gdiscr$z
		 	 zloc_mid = zloc+0.5*zdiscr
			}
			if((edge =="last" & layer==layermax) | (edge=="both" & layer==layermax)){ # last z-layer
			 zdiscr = 0.5*gdiscr$z
		  	 zloc_mid = zloc-0.5*zdiscr
			}
			else{
		 	zloc_mid = zloc
		 	zdiscr = gdiscr$z
			}
		}
		else{ # all other z-layers
		 zloc_mid = zloc
		 zdiscr = gdiscr$z
		}
           gcomp_cell=macmillan_raw(gama,xloc,yloc,zloc_mid,gloc$x,gloc$y,gloc$z,gdiscr$x,gdiscr$y,zdiscr,rad,w,rho) #unit depends on w
        }
# output one value: gravity component for corresponding input cell
return (gcomp_cell)
}

#' @title Interpolate grid(df) to 3D grid
#'
#' @description interpolate DEM to necessary grid of corresponding layer
#'
#' @param grid_cords coordinates of the original hydrus mesh.
#' @param data_input data to be interpolated. in fomat of the hydrus mesh.
#' @param grid_discr vector of discretization of new grid in (x,y,z).
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

surface_to_grid3d = function(
            surface_grid,
            grid_discr,
            depth_split
){
    # library(gstat)
    if(!is.null(surface_grid)){
    	grid.z = seq(min(depth_split), max(depth_split), by=grid_discr$z)
    	grid.xyz = data.frame(x=rep(surface_grid$x,length(grid.z)),
                               y=rep(surface_grid$y,length(grid.z)),
                               z=rep(grid.z, each=length(surface_grid[,1])))
    }else{
    	grid.x = seq(min(Building_x), max(Building_x), by=grid_discr$x)
    	grid.y = seq(min(Building_y), max(Building_y), by=grid_discr$y)
    	grid.z = seq(min(depth_split), max(depth_split), by=grid_discr$z)
    	grid.xyz = expand.grid(x=grid.x, y=grid.y, z=grid.z)
    	surface_grid = expand.grid(x=grid.x, y=grid.y)
    }
    
    # interpolate dem data to new grid in 3d
    # idw.gstat = gstat(formula = z ~ 1, locations = ~ x + y, data = demgrid_in, nmax = 4, set = list(idp = 2))
    idw.gstat = gstat(formula = z ~ 1, locations = ~ x + y, data = surface_grid, nmax = 4, set = list(idp = 2))
    surface = predict(idw.gstat, surface_grid)
    grid3d = cbind(grid.xyz[,1:2],
                   z=surface[,3] - grid.xyz$z #,
                   # zgrid=grid.xyz$z, 
                   # layer=rep(seq(1,length(grid.z), by=1),each=length(surface_grid$x))
                   )
    return(grid3d)
}


#' @title Convert DEM to data.frame
#'
#' @description interpolate DEM to necessary grid of corresponding layer
#'
#' @param grid_cords coordinates of the original hydrus mesh.
#' @param data_input data to be interpolated. in fomat of the hydrus mesh.
#' @param grid_discr vector of discretization of new grid in (x,y,z).
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

convert_demtodf = function(dem, dem.info){
# convert (melt) dem
dem_x = seq(dem.info[3,1],length.out=length(dem[1,]),by=dem.info[5,1])
dem_y = seq(dem.info[4,1],length.out=length(dem[,1]),by=dem.info[5,1])
colnames(dem) = dem_x
rownames(dem) = rev(dem_y)
dem.melt =  melt(as.matrix(dem))#, id.vars=
colnames(dem.melt) = c("y","x","z")
return(dem.melt)
} #end function


#' @title Read Digital Elevation Models 
#'
#' @description Read DEM's and output information files for further usage. 
#'
#' @param dempath Path of input DEM-file. 
#' @param filename DEM-file in .acs-format (see details)
#' 
#' @details So far only it is only possible to read DEM-file in .asc arquitecture format.
#' @details Outputs files: data.DEM, info.DEM (as data.frames).
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing
#' 
read_dem <- function(dempath,filename){
info.DEM <- read.table(file=paste(dempath,filename,sep=""),nrows=6, row.names=1, header=F, sep="", stringsAsFactors=F, colClasses=c("character","numeric"))
nodata <- info.DEM[6,1]
data.DEM <- read.table(file=paste(dempath,filename,sep=""), sep="", dec=".", skip=6, stringsAsFactors=F, na.strings=nodata, header=F)
dem.info <<- info.DEM; dem<<- data.DEM
sucess = "output: dem.info and dem"
return(sucess)
} ### end read DEM ###


#' @title Calulate gravity signals 3D
#'
#' @description After calculating the gravity components (gcomp()), this function can be used to get real values using soil moisture data (theta). The output of the gravity signal is possible in different dimensions. 
#'
#' @param gcompfiles vector containing the names (complete paths or relative) of the gravity component files to be used for calculations. These files are generated using gcomp().
#' @param theta data.frame with column structure $time, $timestep, $theta (value), $layer
#' @param output Defines the output to be a singel value, layer or grid (default is value)
#' ...
#' @details Usually one needs relative gravity signals, therefore the inputed theta timeseries should be acutally delta theta vaules per timestep. The
#' other option is calculating the differences in values per timestep of the output of this function.
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing 
#' @export
#' 
# gsignal_grids_3d <- function(gcompfile,SMdif_mod, outputtype="signal"){
gsignal_grids_3d <- function(gcompfile,SMdif_mod, realdates=T){

#SMdif_mod$x = SMdif_mod$x + min(gcompfile$x)
#SMdif_mod$y = SMdif_mod$y + min(gcompfile$y)

##sort both data after same varianles
#gcomp_sort = arrange(gcompfile, x,y,zgrid)
#rm(gcompfile); gc()
#SMdif_sort = arrange(SMdif_mod, xDEM,yDEM,Depth)
#rm(SMdif_mod); gc()
#SMdif_comp = data.frame(Timestep = SMdif_sort$Timestep, SMobs=SMdif_sort$value, gcomp=gcomp_sort$gcomp)
#SMdif_comp$gsignal = SMdif_comp$SMobs * SMdif_comp$gcomp

#gsignals = left_join(SMdif_mod, gcompfile,by=c("xDEM" = "x", "yDEM" = "y", "Depth" = "zgrid" )) %>%
#gsignals = left_join(gcompfile, SMdif_mod, by=c("x" = "xDEM", "y" = "yDEM", "zgrid" = "Depth")) %>%
if(realdates){
gsignals = dplyr::left_join(gcompfile, SMdif_mod) %>%
		dplyr::mutate(gsignal = gcomp * value) %>%
		dplyr::group_by(datetime) %>%
		dplyr::summarize(gmod = sum(gsignal, na.rm=T))
}else{
gsignals = dplyr::left_join(gcompfile, SMdif_mod) %>%
		dplyr::mutate(gsignal = gcomp * value) %>%
		dplyr::group_by(Timestep) %>%
		dplyr::summarize(gmod = sum(gsignal, na.rm=T))
}
# select output type
#switch(outputtype,
       #signal = {
		#gsignals = group_by(SMdif_comp, Timestep) %>%
		#summarize(gmod = sum(gsignal, na.rm=T))
       #},
       #grid = {
		#gsignals = SMdif_comp
       #}
       #) # end switch

return(gsignals)
}

##############################
## some random


# grid2d_cords = read_2dgrid(fol_path)
# save(grid2d_cords, file=paste(homedir,fol,"output/Hydrus2d_grid.rdata",sep=""))

#########################
#########################

# ## vertical SM data for grids outside of SGbuilding (only 1D [vertical] information necessary)
# # radius 50 m & 300 m
# grid_r50m_r300m = c(0.1,0.1); depth_r50m_r300m = c(0,5) #until and including z=0.1
# nmax=1
# SMgrid_r50m_r300m = nodestogrid_2d(grid2d_cords, SMdata_hydrus2d, grid_r50m_r300m, depth_r50m_r300m,nmax)
# #get mean over all x locations
# SMvertical_r50m_r300m = group_by(SMgrid_r50m_r300m, Timestep,Depth) %>%
#             filter(Timestep > 0) %>%
#             summarize(value = mean(value, na.rm=T))
# #colnames(SMvertical_r50m_r300m)[3] = "value"
# save(SMvertical_r50m_r300m, file=paste(homedir,fol,"output/SMdata_TS_1Dvertical_interpolated_allLayers_r50m_r300m.rdata",sep=""))
# 
# # radius 5000 m
# # average 1d grid with higher resolution to 1 m vertical resolution
# transformGrid = data.frame(Depth = unique(SMvertical_r50m_r300m$Depth), layer = c(rep(1,11),rep(2,10),rep(3,10),rep(4,10),rep(5,10)))
# SMvertical_r5000m = inner_join(SMvertical_r50m_r300m, transformGrid, by="Depth") %>%
#         group_by(Timestep, layer) %>%
#         summarize(value=mean(value, na.rm=T))
# colnames(SMvertical_r5000m)[2] = "Depth"



#' @title test
#'
#' @description test
#'
#' @param test test
#' 
#' @return test
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

forsberg_raw <- function(gama,w,xl,xr,yl,yr,Zint,Zend,xs,ys,zs,rho){

  #distance mass to SG
  x=0
  x[1]=xl-xs
  x[2]=xr-xs
  y=0
  y[1]=yl-ys
  y[2]=yr-ys
  z=0
  z[1]=Zint-zs
  z[2]=Zend-zs

  sum=0
  for (i in 1:2){
    for (ii in 1:2){
      for (iii in 1:2){
        rf=sqrt(x[i]^2+y[ii]^2+z[iii]^2)
        sum=sum+(-1)^(i+ii+iii)*(x[i]*log(y[ii]+rf)+y[ii]*log(x[i]+rf)-z[iii]*atan(x[i]*y[ii]/z[iii]/rf))
      } #end iii
    } #end ii
  } #end i
  d_forsberg=-w*gama*rho*sum #in µGal, if w = 1e8 
  return(d_forsberg)

}

#' @title test
#'
#' @description test
#'
#' @param test test
#' 
#' @return test
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

macmillan_raw <- function(gama,xp,yp,zp,xs,ys,zs,dx,dy,dz,rad,w,rho){
  
  # calculations for distances
  alfa=2*dx^2-dy^2-dz^2
  beta=-dx^2+2*dy^2-dz^2
  ome=-dx^2-dy^2+2*dz^2
  abg=alfa*(xp-xs)^2+beta*(yp-ys)^2+ome*(zp-zs)^2
  # 3 different macmillan terms
  tm1=-((zp-zs)/rad^3)
  tm2=-5/24*(zp-zs)*abg/rad^7
  tm3=ome/12*(zp-zs)/rad^5 ##12!?! benjamin hatte hier mal eine 24 stehen...warum??
  # multiply together for final result for this layer, spacial extent R&C and for one step in time
  d_macmillan=w*gama*rho*dx*dy*dz*(tm1+tm2+tm3) #in µGal, if w = 1e8 #NEGATIVE SIGN REMOVED!!
  return(d_macmillan)
  
}

#' @title test
#'
#' @description test
#'
#' @param test test
#' 
#' @return test
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

pointmass <- function(gama,zp,zs,dx,dy,dz,rad,w,rho){
  
  d_pointmass=-w*gama*rho*dx*dy*dz*(zp-zs)/rad^3 #in µGal, if w = 1e8 
  return(d_pointmass)
  
}


#' @title Convert zoo-object to data frame
#'
#' @description Converting zoo-object to data frame with proper indexing of time column and column names
#'
#' @param value must be a zoo-object
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
# ' @examples examples MISSING
#' 
zootodf <- function(value) {
    df <- as.data.frame(value)
    df$time <- index(value) #create a Date column
    rownames(df) <- NULL #so row names not filled with dates
    df <- df[,c(ncol(df), 1:(ncol(df)-1))] #reorder columns so Date first
    return(df)
}
