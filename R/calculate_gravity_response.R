#' @title Calulate gravity response for 3d grids
#'
#' @description Calculate integrative gravity response, resulting from the mass changes due to changes in soil moisture contents. 
#'
#' @param gcomp_grid Data.frame, containing the gravity component grid (x,y,z, gcomp).
#' @param mass_input Data.frame, containing columns structure $datetime, $value, $Depth.
#' Important: soil moisture values have to be in \%, not integer values. Otherwise the order of the resulting gravity response is by a factor 100.
#' 
#' @details test
#'
#' @return Returns a data.frame containing columns of time information and the corresponding integrative gravity response.
#'
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing 

calculate_gravity_response = function(
            gcomp_grid,
            mass_input
){
    # gcomp_grid = gravity_component_grid3d
    # mass_input = SMgrid3d_outside
    # exclude z-columns of soil moisture dataset for smoother joining
    mass_input = dplyr::select(mass_input, -z)
    # join datasets (by x,y,z) and multiply
    # gravity component with soil moisture content of each cell
    gsignals = dplyr::left_join(gcomp_grid, mass_input, by=c("x","y","Depth")) %>%
		dplyr::mutate(gsignal = gcomp * value) %>%
		dplyr::group_by(datetime) %>%
		dplyr::summarize(value = sum(gsignal, na.rm=T))
    # subtract first value
    # this is done to have a correct reference of cumsum
    # for comparison with other gravity (cumsum) time series
    gsignals$value = gsignals$value - gsignals$value[1]
    # return data
    return(gsignals)
}
