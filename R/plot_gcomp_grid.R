#' @title Plot transect of gravity component grid
#'
#' @description Plots a transect of a supplied input grid, meaning that a 3d grid is reduced to a 2d
#' transect on basis of a given axis (coordinate).
#'
#' @param grid_input Data.frame, containing spatial coordinates in the columns (x, y, z). 
#' @param yloc Numeric, coordinate of the y-axis, along which the transect will be constructed.
#' @param output_dir a
#' 
#' @return Returns a ggplot figure. The plot is also saved to the specified output directory.
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

plot_gcomp_grid = function(
            grid_input,
            yloc,
            output_dir,
            grid_discretization,
            ...
){
    # grid_input = as.data.frame(gravity_component_grid3d)
    # yloc = round(SG_y, 1)
    # finding closest y-coordinate to supplied yloc
    y_coords = unique(grid_input$y)
    yloc_plot = y_coords[which.min(abs(y_coords - yloc))]
    # set x and z-differences for plotting correct size of tiles
    x_dif = unique(grid_input$x)[3] - unique(grid_input$x)[2]
    # z_dif = unique(gravity_component_grid3d$z)[3] - unique(gravity_component_grid3d$z)[4]
    z_dif = grid_discretization$z

    # reduce grid to one 2d transect along yloc
    grid_transect = grid_input %>%
        dplyr::filter(y == yloc_plot) %>%
        # dplyr::mutate(y_rnd = round(y, 1)) %>%
        # dplyr::mutate(y_exclude = y_rnd - yloc) %>%
        # dplyr::filter(abs(y_exclude) < .9) %>%
        dplyr::mutate(gcomp = ifelse(gcomp == 0, NA, gcomp)) %>%
        ggplot(aes(x = x, y = z)) + 
        geom_tile(aes(fill = gcomp), width = x_dif, height = 2 * z_dif) + 
        # geom_raster(aes(fill = gcomp)) + 
        scale_fill_gradient(low = viridis(10)[3], high = viridis(10)[8], na.value = "gray") + 
        # geom_point(aes(colour = gcomp)) +
        # scale_color_gradient(low = viridis(10)[3], high = viridis(10)[8], na.value = "gray") + 
        ylab("Depth") + xlab("Grid x-axis") + 
        labs(fill = "Gravity component")

    # save plot
    png(filename = paste0(output_dir, "Gravity_component_grid_transect.png"),
                      width = 600,
                      height = 400,
                      res = 150)
    print(grid_transect)
    dev.off()

    return(grid_transect)
}

