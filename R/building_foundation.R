#' @title Constructs a grid of the SG buidling foundation
#'
#' @description Gridpoints for all points in space which belong to the
#' SG building. This includes walls, the baseplate and one SG pillar.
#' This foundation grid can later on be used to correct the gravity component grid,
#' in order to set mass changes equal zero at the foundation.
#'
#' @param Bdwall_ext_x,y,z Numeric, giving the extent of the wall in x,y,z directions.
#' @param Bdbase_x,y,z Vector, giving min and max x,y,z values of the coordinates from the baseplate of the building.
#' @param Pillar_x,y,z Vector, giving min and max x,y,z values of the coordinates from the SG pillar within the building.
#' @param grid_discretization data.frame with columns (x,y,z), indicating the discretization in the corresponding direction.
#' 
#' @return A data frame including x,y and z coordinates of walls, baseplate and SG pillar of the obervatory building.
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

building_foundation = function(
            Bdwall_ext_x,
            Bdwall_ext_y,
            Bdwall_ext_z,
            Bdbase_x,
            Bdbase_y,
            Bdbase_z,
            Pillar_x,
            Pillar_y,
            Pillar_z,
            grid_discretizaion
){
    # building walls
    x.walls = seq(min(Bdbase_x),max(Bdbase_x),by=grid_discretizaion$x)
    yx.walls = c(seq(min(Bdbase_x),(min(Bdbase_x) + Bdwall_ext_y),by=grid_discretizaion$y),seq((max(Bdbase_x) - Bdwall_ext_y),max(Bdbase_x),by=grid_discretizaion$y))
    y.walls = seq(min(Bdbase_y),max(Bdbase_y), by=grid_discretizaion$y)
    xy.walls = c(seq(min(Bdbase_y),(min(Bdbase_y) + Bdwall_ext_x),by=grid_discretizaion$x),seq((max(Bdbase_y) - Bdwall_ext_x),max(Bdbase_y),by=grid_discretizaion$x))
    z.walls = seq(min(Bdbase_z),(min(Bdbase_z) + Bdwall_ext_z), by=grid_discretizaion$z)
    walls.x = expand.grid(x=x.walls, y=xy.walls, z=z.walls)
    walls.y = expand.grid(x=yx.walls, y=y.walls, z=z.walls)
    # building baseplat
    x.base = seq(min(Bdbase_x),max(Bdbase_x), by=grid_discretizaion$x)
    y.base = seq(min(Bdbase_y),max(Bdbase_y), by=grid_discretizaion$y)
    z.base = seq(min(Bdbase_z),max(Bdbase_z), by=grid_discretizaion$z)
    base = expand.grid(x=x.base, y=y.base, z=z.base)
    # SG pillar (in center of building)
    ## ! should add relative coordinates of SG real
    ## !
    x.SGpillar = seq(min(Pillar_x),max(Pillar_x), by=grid_discretizaion$x)
    y.SGpillar = seq(min(Pillar_y),max(Pillar_y), by=grid_discretizaion$y)
    z.SGpillar = seq(min(Pillar_z),max(Pillar_z), by=grid_discretizaion$z)
    SGpillar = expand.grid(x=x.SGpillar, y=y.SGpillar, z=z.SGpillar)
    
    # combine all parts
    BdFoundation = rbind(walls.x, walls.y, base, SGpillar)
    
    return(BdFoundation)
}

