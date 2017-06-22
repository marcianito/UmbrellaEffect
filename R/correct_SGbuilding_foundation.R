#' @title Correct the gravity component grid for the SG building foundation
#'
#' @description In order to set all mass variations to zero where a firm and fix building can be found,
#' this function modifies the gravity component grid in the way that all gravity components laying physically
#' within the building foundations, are set to zero.
#'
#' @param gravity_gcomp Data.frame, 4 columns: spatial information (x,y,z) and corresponding gravity component (value).
#' @param building_foundation Data.frame, 3 columns: spatial information (x,y,z) about the stucture of the foundation of the SG building.
#' 
#' @return Returns a data.frame, consisting of 4 columns: 3 dimensions in space (x,y,Z) and value, holding the gravity component data.
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

correct_SGbuilding_foundation = function(
            gravity_gcomp3d,
            Bdwall_ext_x,
            Bdwall_ext_y,
            Bdwall_ext_z,
            Bdbase_x,
            Bdbase_y,
            Bdbase_z,
            Pillar_x,
            Pillar_y,
            Pillar_z,
            grid_discretization
){
    # gravity_gcomp3d = gravity_component_grid3d
    # Bdwall_ext_x = Building_walls_x
    # Bdwall_ext_y = Building_walls_y
    # Bdwall_ext_z = Building_walls_z
    # Bdbase_x = Building_x
    # Bdbase_y = Building_y
    # Bdbase_z = Building_baseplate_z
    # Pillar_x = Building_SGpillar_x
    # Pillar_y = Building_SGpillar_y
    # Pillar_z = Building_SGpillar_z
    # grid_discretization = grid3d_discr

    # construct faces
    faces = as.matrix(data.frame(
            x = c(1, 7, 1, 6, 1, 4, 2, 8, 3, 8, 5, 8),
            y = c(5, 3, 2, 5, 3, 2, 4, 6, 7, 4, 6, 7),
            z = c(3, 5, 5, 2, 2, 3, 6, 4, 4, 7, 7, 6)
    ))

    # round z values in case of UTM coordinates
    round_x = decimalplaces(grid_discretization$y)
    round_y = decimalplaces(grid_discretization$x)
    round_z = decimalplaces(grid_discretization$z)

    # construct vertices
    # for SG pillar
    vert_SGpillar = construct_vertices(
                        x_cords = Pillar_x,
                        y_cords = Pillar_y,
                        z_cords = Pillar_z
    )
    # for building baseplate
    vert_Bdbase = construct_vertices(
                        x_cords = Bdbase_x,
                        y_cords = Bdbase_y,
                        z_cords = c(min(Bdbase_z), max(round(Bdbase_z, round_z)) + grid_discretization$z)
    )
    # for building walls 
    vert_Bdwallxy1 = construct_vertices(
                        x_cords = Bdbase_x,
                        y_cords = c(min(Bdbase_y), min(Bdbase_y) + Bdwall_ext_y),
                        z_cords = c(min(Bdbase_z) - Bdwall_ext_z,min(Bdbase_z))
    )
    vert_Bdwallxy2 = construct_vertices(
                        x_cords = Bdbase_x,
                        y_cords = c(max(Bdbase_y), max(Bdbase_y) - Bdwall_ext_y),
                        z_cords = c(min(Bdbase_z) - Bdwall_ext_z,min(Bdbase_z))
    )
    vert_Bdwallyx1 = construct_vertices(
                        x_cords = c(min(Bdbase_x), min(Bdbase_x) + Bdwall_ext_x),
                        y_cords = Bdbase_y,
                        z_cords = c(min(Bdbase_z) - Bdwall_ext_z,min(Bdbase_z))
    )
    vert_Bdwallyx2 = construct_vertices(
                        x_cords = c(max(Bdbase_x), max(Bdbase_x) - Bdwall_ext_x),
                        y_cords = Bdbase_y,
                        z_cords = c(min(Bdbase_z) - Bdwall_ext_z,min(Bdbase_z))
    )

    # try out to play with ENLARGING area of vertices
    # not working propery !
    # for building walls 
    # vert_Bdwallxy1 = construct_vertices(
    #                     x_cords = Bdbase_x + c(- .9 * grid_discretization$x, .9 * grid_discretization$x),
    #                     y_cords = c(min(Bdbase_y), min(Bdbase_y) + Bdwall_ext_y) + c(- .9 * grid_discretization$y, .9 * grid_discretization$y),
    #                     z_cords = c(min(Bdbase_z) - Bdwall_ext_z,min(Bdbase_z)) + c(- .9 * grid_discretization$z, .9 * grid_discretization$z)
    # )
    # vert_Bdwallxy2 = construct_vertices(
    #                     x_cords = Bdbase_x+ c(- .9 * grid_discretization$x, .9 * grid_discretization$x),
    #                     y_cords = c(max(Bdbase_y), max(Bdbase_y) - Bdwall_ext_y)+ c(- .9 * grid_discretization$y, .9 * grid_discretization$y),
    #                     z_cords = c(min(Bdbase_z) - Bdwall_ext_z,min(Bdbase_z))+ c(- .9 * grid_discretization$z, .9 * grid_discretization$z)
    # )
    # vert_Bdwallyx1 = construct_vertices(
    #                     x_cords = c(min(Bdbase_x), min(Bdbase_x) + Bdwall_ext_x)+ c(- .9 * grid_discretization$x, .9 * grid_discretization$x),
    #                     y_cords = Bdbase_y+ c(- .9 * grid_discretization$y, .9 * grid_discretization$y),
    #                     z_cords = c(min(Bdbase_z) - Bdwall_ext_z,min(Bdbase_z))+ c(- .9 * grid_discretization$z, .9 * grid_discretization$z)
    # )
    # vert_Bdwallyx2 = construct_vertices(
    #                     x_cords = c(max(Bdbase_x), max(Bdbase_x) - Bdwall_ext_x)+ c(- .9 * grid_discretization$x, .9 * grid_discretization$x),
    #                     y_cords = Bdbase_y+ c(- .9 * grid_discretization$y, .9 * grid_discretization$y),
    #                     z_cords = c(min(Bdbase_z) - Bdwall_ext_z,min(Bdbase_z))+ c(- .9 * grid_discretization$z, .9 * grid_discretization$z)
    # )
    
    # check and remove if points of gravity component grid
    # lay in one of the above constructed rectangulars
    # = SG pillar, building baseplate and walls
    gcomp_grid3d = gravity_gcomp3d %>%
        dplyr::mutate(SGpillar = pip3d(vert_SGpillar, faces, as.matrix(gravity_gcomp3d[,1:3]))) %>%
        dplyr::mutate(Bd_baseplate = pip3d(vert_Bdbase, faces, as.matrix(gravity_gcomp3d[,1:3]))) %>%
        dplyr::mutate(Bd_wallsxy1 = pip3d(vert_Bdwallxy1, faces, as.matrix(gravity_gcomp3d[,1:3]))) %>%
        dplyr::mutate(Bd_wallsxy2 = pip3d(vert_Bdwallxy2, faces, as.matrix(gravity_gcomp3d[,1:3]))) %>%
        dplyr::mutate(Bd_wallsyx1 = pip3d(vert_Bdwallyx1, faces, as.matrix(gravity_gcomp3d[,1:3]))) %>%
        dplyr::mutate(Bd_wallsyx2 = pip3d(vert_Bdwallyx2, faces, as.matrix(gravity_gcomp3d[,1:3]))) %>%
        dplyr::mutate(gcomp = ifelse(SGpillar >= 0 | Bd_baseplate >= 0 | Bd_wallsxy1 >= 0 | Bd_wallsxy2 >= 0 | Bd_wallsyx1 >= 0 | Bd_wallsyx2 >= 0, 0, gcomp)) %>%
        dplyr::select(-SGpillar, -Bd_baseplate, -Bd_wallsxy1, -Bd_wallsxy2, -Bd_wallsyx1, -Bd_wallsyx2)

    return(gcomp_grid3d)
}

