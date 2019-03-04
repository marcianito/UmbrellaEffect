#' @title Convert gravity data
#'
#' @description Converts gravity data from outside of the building to the gravity response
#' below the building. This is done multiplying with factors 1 & 2, both
#' obtained as a consequence of dominant hydrological scenario, SG position and SG building size.
#'
#' @param gravity_input Data.frame, time series containing time information (datetime) and values of the gravity response from outside of the building.
#' @param factor_hydScen_SGloc Data.frame, time series having two columns: time information (datetime) and
#' the temporal dynamic factor for reduction due to hydrological scenario and SG position.
#' @param factor_BdSize Numeric, one value containing the reduction factor due to the influence of the size of the SG building.
#' 
#' @return Returns a data.frame (time series), which has information about time (as the gravity_input time series) and
#' as value the converted gravity response for below the SG building. This converted gravity time series can now be 
#' used to reduce observed gravity data from this very SG building.
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing
#' @export

convert_gravity_response = function(
            gravity_input,
            factor_hydScen_SGloc,
            factor_BdSize
){
    gravity_converted = gravity_input %>%
        # dplyr::left_join(factor_hydScen_SGloc) %>%
        dplyr::inner_join(factor_hydScen_SGloc) %>%
        # new
        dplyr::filter(is.na(value) == F) %>%
        dplyr::filter(is.na(fac1) == F) %>%
        #
        dplyr::mutate(fac2 = factor_BdSize) %>%
        dplyr::mutate(value = value * fac1 * fac2) %>%
        dplyr::select(-fac1, -fac2)

    return(gravity_converted)
}

