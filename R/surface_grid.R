#' @title Surface grid
#'
#' @description Generates a 2d grid of the surface of topography at the observatory building.
#' This is done on basis of a DEM or assuming a flat surface.
#'
#' @param DEM Character string, containing the name of the DEM-file to use.
#' If left empty, a flat topography with value 0 (zero) will be assumed.
#' @param grid_domain_x test
#' @param input_dir Character string, specifying the path of the directory where the DEM-file is located.
#' @param output_dir Character string, specifying the path of the directory where output should be stored.
#' 
#' @return Returns a data.frame, which holds topographical information, including x,y coordinates,
#' from the area of the SG building.
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

surface_grid = function(
            DEM,
            grid_domain_x,
            grid_domain_y,
            input_dir,
            output_dir,
            ...
){
    # grid_domain_x = Building_x
    # grid_domain_y = Building_y
    # DEM = DEM_input_file
    # input_dir = dir_input

    ## check if DEM is non-empty, if not return NULL
    if(DEM == ""){
        surface_grid = NULL
    }else{
        ## implement reading-in routines for various DEM-file extensions !!
        ## !
        ## what to do when DEM is not supplied (== NULL)
        ## !
        
        # read DEM as raster
        dem_raster = raster(paste0(input_dir,DEM))
        
        # create local grid domain (extension),
        # based on building coordinates
        grid_domain = data.frame(x=c(min(grid_domain_x) ,min(grid_domain_x), max(grid_domain_x), max(grid_domain_x)),
                                 y=c(min(grid_domain_y), max(grid_domain_y), max(grid_domain_y), min(grid_domain_y))
                                 )
        extent_grid_domain = extent(grid_domain)
        # limit DEM to extent of SG building
        # dem_grid_domain = crop(dem_raster, extent_grid_domain)
        # DEM had to be snapped, in order to get the correct number of elements (its extent)
        # as a non-DEM surface grid on the basis of Building_x,y
        dem_grid_domain = crop(dem_raster, extent_grid_domain, snap= "out")
        writeRaster(dem_grid_domain, filename = paste0(output_dir, "dem_grid"), format="ascii", NAflag=-9999, overwrite=T)
        
        # reload grid for verification
        # in order to create the surface grid
        ## ! here has to be found a more straight forward way,
        # going directly from dem_grid_domain to surface_grid
        # WITHOUT writing a physical file !!
        dem_grid_file = "dem_grid.asc" 
        read_dem(output_dir, dem_grid_file)
        # # plot DEM 
        # plot_dem(dem, dem.info, locations)
        # convert DEM grid to data.frame
        surface_grid = convert_demtodf(dem, dem.info)
    }

    return(surface_grid)
}

