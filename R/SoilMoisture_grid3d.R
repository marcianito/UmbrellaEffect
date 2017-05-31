#' @title Soil moisture grid (3d)
#'
#' @description Generates a grid in 3d, containing time series
#' at each point of soil moisture. 
#'
#' @param grid_domain Data.frame, constructed with minimum 3 columns (x,y,z) which describe spatial dimensions of 3d grid.
#' @param soilMoisture_input Name of file, containing the soil moisture time series.
#' This file has to contain minimal 2 columns (datetime and z) but can optionally also contain (x and/or y), if a higher dimensional soil moisture dataset is available.
#' @param input_dir Directory where
#' 
#' @return Returns a data.frame with 5 columns.
#' Information about time (datetime), spatial information (x,y,z) and value which is the soil moisture data.
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

SoilMoisture_grid3d = function(
            grid_domain = gravity_component_grid3d,
            soilMoisture_input = soilMoisture_input_file,
            input_dir = dir_input,
            ... # should be able to pass stuff like sep="", etc.
){
    # include various routines to read
    # different file types:
    # .rData, .csv, .tsf, 

    SMdata = 

    ## probably the best to average all available 1 d
    # data per depth layer and then
    # distribute (join) to existing 3d grid !?
    # how to process when someone has 3d data !?
    # how if someone has 2d model data?

    # round data to make use join is performed corretely !
    SMmod_NObd$Depth = round(SMmod_NObd$Depth, 1)
    # create SM 3d grid with data from NEXT to building
    SMgrid3d_outside = dplyr::left_join(dplyr::select(grid_domain, -value), SMdata)# %>%

    return(SMgrid3d_outside)
}

