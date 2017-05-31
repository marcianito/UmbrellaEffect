#' @title Gravity component grid
#'
#' @description Generates a grid of gravity components.
#'
#' @param surface this
#' @param SG_coordinates this
#' @param grid_discretization this
#' @param grid_depth this
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
            # grid_domain = grid3d,
            SG_coordinates,
            grid_discretizaion,
            grid_depth
){
    # generate 3d grid
    grid3d = surface_to_grid3d(
            surface_grid = surface,
            grid_discr = grid_discretizaion,
            depth_split = grid_depth) #,
            # T, surface)
    # generate gravity component grid
    gcomp_grid = fill_gcompgrid(
            g_grid = grid3d,
            senloc = SG_coordinates,
            g_discr = grid_discretizaion
    )

    return(gcomp_grid)
}

