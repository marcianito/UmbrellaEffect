#' @title Gravity component grid
#'
#' @description Generates a grid of gravity components.
#'
#' @param surface this
#' @param SG_coordinates this
#' @param grid_discretization this
#' @param grid_depth this
#' @param range_coords_x,y this
#' 
#' @return Returned is a data.frame with 4 columns.
#' 3 for coordinates in space (x,y,z) and one with
#' the gravity component corresponding to this gridpoint.
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

gravity_comp_grid = function(
            surface,
            SG_coordinates,
            grid_discretization,
            grid_depth,
            range_coords_x,
            range_coords_y
){
    # surface = surface_grid
    # SG_coordinates = SGloc
    # grid_discretization = grid3d_discr
    # grid_depth = grid3d_depth
    # range_coords_x = sprinklingArea_x
    # range_coords_y = sprinklingArea_y

    # generate 3d grid
    grid3d = surface_to_grid3d(
            surface_grid = surface,
            grid_discr = grid_discretization,
            depth_split = grid_depth,
            Bd_x = range_coords_x,
            Bd_y = range_coords_y
    )
    # generate gravity component grid
    gcomp_grid = fill_gcompgrid(
            g_grid = grid3d,
            senloc = SG_coordinates,
            g_discr = grid_discretization
    )
    # # take out column "layer", which was needed only for internal calculations
    # gcomp_grid = dplyr::select(gcomp_grid, -layer)

    # round all x,y,z to same decimal places
    # if not, joining might not be complete!
    round_x = decimalplaces(grid_discretization$x)
    round_y = decimalplaces(grid_discretization$y)
    round_z = decimalplaces(grid_discretization$z)
    gcomp_grid$x = round(gcomp_grid$x, round_x)
    gcomp_grid$y = round(gcomp_grid$y, round_y)
    gcomp_grid$z = round(gcomp_grid$z, round_z)
    gcomp_grid$Depth = round(gcomp_grid$Depth, round_z)

    return(gcomp_grid)
}

#' @title Calculate gravity component, WITHOUT umbrella effect
#'
#' @description Calculates the gravity components (100percent values) for a given DEM using a nested approach
#'
#' @param g_grid grid for which to calculate gravity components. (generally derived from DEM of region in combination with hydrological modeling mesh/grid).
#' @param senloc coordinates of gravity sensor location. data.frame with x,y,z columns.
#' @param g_discr discretization of g_grid in differences in x,y,and z-direction. they have to be uniform steps!
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples example MISSING
#' @export

fill_gcompgrid <- function(g_grid, senloc, g_discr){
#constants
gamma <- 6.673e-11 #m³/(kg*s²)
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
	mutate(gcomp = select_gMethod_Edges(x,y,z,senloc,g_discr, edge, layer, layern_max,r_inner,r_outer,gamma,rho,w)) %>%
	ungroup()

return(g_grid)
}

#' @title Choose method for computing gravity component for a cell of a grid
#'
#' @description test
#'
#' @param xloc  t
#' @param yloc  t
#' @param zloc  t
#' @param gloc  t
#' @param gdiscr  t
#' @param edge t
#' @param layer  t
#' @param layermax  t
#' @param r_inner  t
#' @param r_outer  t
#' @param gamma t
#' @param rho t
#' @param w t
#' ...
#' @details test
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing 
#' @export
#' 
select_gMethod_Edges = function(xloc, yloc, zloc, gloc, gdiscr, edge, layer, layermax, r_inner, r_outer, gamma, rho, w){
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
           gcomp_cell=forsberg_raw(gamma,w,xl,xr,yl,yr,Zint,Zend,gloc$x,gloc$y,gloc$z,rho) #unit depends on w
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
           gcomp_cell=pointmass(gamma,zloc_mid,gloc$z,gdiscr$x,gdiscr$y,zdiscr,rad,w,rho) #unit depends on w
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
           gcomp_cell=macmillan_raw(gamma,xloc,yloc,zloc_mid,gloc$x,gloc$y,gloc$z,gdiscr$x,gdiscr$y,zdiscr,rad,w,rho) #unit depends on w
        }
# output one value: gravity component for corresponding input cell
return (gcomp_cell)
}
