#' @title Soil moisture grid (3d)
#'
#' @description Generates a grid in 3d, containing time series
#' at each point of soil moisture. 
#'
#' @param grid_domain Data.frame, constructed with minimum 3 columns (x,y,z) which describe spatial dimensions of 3d grid.
#' @param soilMoisture_input Name of file, containing the soil moisture time series.
#' This file has to contain minimal 2 columns (datetime and z) but can optionally also contain (x and/or y), if a higher dimensional soil moisture dataset is available.
#' @param grid_discretization data.frame with columns (x,y,z), indicating the discretization in the corresponding direction.
#' @param input_dir Directory where
#' 
#' @return Returns a data.frame with 5 columns.
#' Information about time (datetime), spatial information (x,y,z) and value which is the soil moisture data.
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

SoilMoisture_grid3d = function(
            grid_domain,
            soilMoisture_input,
            grid_discretization,
            grid_depth,
            input_dir,
            ... # should be able to pass stuff like sep="", etc.
){
    # soilMoisture_input = soilMoisture_input_file
    # input_dir = dir_input

    # read input data
    SMdata_in = read_data(soilMoisture_input, input_dir)

    # create results data.frame
    SMdata = data.frame(
                        datetime = rep(unique(SMdata_in$datetime), each = length(unique(grid_domain$Depth))),
                        layer = rep(unique(grid_domain$layer), length(unique(SMdata_in$datetime))),
                        value = NA
    )
    # run for every timestep
    for(ts in unique(SMdata_in$datetime)){
    # subset SM data
    SM_sub = dplyr::filter(SMdata_in, datetime == ts) %>%
             dplyr::select(-datetime)
    # interpolate data to new resolution
    SM_int_tmp = approx(SM_sub, xout = unique(grid_domain$Depth))
    # SM_int_tmp = approx(SM_sub, xout = grid_new$Depth)
    # combine data
    SMdata$value[which(SMdata$datetime == ts)] = SM_int_tmp$y
    }

    # 2) combine x,y information of gcomp_grid with this inter/extrapolated data

    # create SM 3d grid with data from NEXT to building
    # based on spatial coordinates from gravity component grid
    SMgrid3d = dplyr::left_join(SMdata, dplyr::select(grid_domain, -gcomp)) %>%
               dplyr::select(-layer)

    # round all x,y,z to same decimal places
    # if not, joining might not be complete!
    round_x = decimalplaces(grid_discretization$x)
    round_y = decimalplaces(grid_discretization$y)
    round_z = decimalplaces(grid_discretization$z)
    SMgrid3d$x = round(SMgrid3d$x, round_x)
    SMgrid3d$y = round(SMgrid3d$y, round_y)
    SMgrid3d$z = round(SMgrid3d$z, round_z)

    return(SMgrid3d)
}

