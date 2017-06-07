#' @title Interpolate grid(df) to 3D grid
#'
#' @description interpolate DEM to necessary grid of corresponding layer
#'
#' @param grid_cords coordinates of the original hydrus mesh.
#' @param data_input data to be interpolated. in fomat of the hydrus mesh.
#' @param grid_discr vector of discretization of new grid in (x,y,z).
#' @param Bd_x,y Vector, containing min and max coordinates of buildings x and y coordinates.
#' They only need to be supplied, if no surface grid is given as input.
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

surface_to_grid3d = function(
            surface_grid,
            grid_discr,
            depth_split,
            Bd_x = Building_x,
            Bd_y = Building_y
){
    grid_discr = grid3d_discr
    depth_split = grid3d_depth

    # decide if a surface grid is supplied
    if(!is.null(surface_grid)){
    	grid.z = rev(seq(min(depth_split), max(depth_split), by=grid_discr$z))
    	grid.xyz = data.frame(x=rep(surface_grid$x,length(grid.z)),
                               y=rep(surface_grid$y,length(grid.z)),
                               z=rep(grid.z, each=length(surface_grid[,1])))
    }else{
    	grid.x = seq(min(Building_x), max(Building_x), by=grid_discr$x)
    	grid.y = seq(min(Building_y), max(Building_y), by=grid_discr$y)
    	grid.z = rev(seq(min(depth_split), max(depth_split), by=grid_discr$z))
    	grid.xyz = expand.grid(x=grid.x, y=grid.y, z=grid.z)
        surface_grid = expand.grid(x=grid.x, y=grid.y, z = max(depth_split))
   
    }
    # interpolate surface data to new grid in 3d
    idw.gstat = gstat(formula = z ~ 1, locations = ~ x + y, data = surface_grid, nmax = 4, set = list(idp = 2))
    surface = predict(idw.gstat, grid.xyz)
    grid3d = cbind(grid.xyz[,1:2],
                   z=surface[,3] + grid.xyz$z,
                   Depth=grid.xyz$z, 
                   layer=rep(seq(1,length(grid.z), by=1),each=(length(unique(grid.xyz$x)) * length(unique(grid.xyz$y))))
                   )

    return(grid3d)
}
