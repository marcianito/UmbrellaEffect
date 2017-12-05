#' @title Receive reduction factor based on hydrological conditon and position of SG
#'
#' @description A dynamic correction ratio is created, based on mean soil moisture contents next to the observatory building and the chosen conditons
#' of a dominant hydrological scenario as well as the position of the SG within its building.
#'
#' @param Scenario Character vector, defining the dominant hydrological condition found at or around the observatory building.
#' For possible options please read the instructions.
#' @param SGlocation Character string, describing the position of the SG within its building.
#' For possible options please read the instructions.
#' @param VertLimit Numeric value of the vertical extent to be considered for the reduction. Values are in meters and have to be positive.
#' @param MeanSoilMoisture Data.frame, time series of the mean soil moisture content in the complete space outside of the observatory building.
#' It needs 2 columns: one for time information (datetime) and one with the soil moisture data (value).
#' 
#' @return Returns a time series of numerical values, which represents the dynamic correction factors based on a dominant hydrological condition
#' and the in-house position of the SG.
#' This can later be used to adjust the modeled gravity response from outside of the building.
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

reduction_hydScen_SGloc = function(
            Scenario,
            SGlocation,
            VertLimit,
            MeanSoilMoisture
){
    # load reduction parameters
    # old file: for complete 5 m vertical model extent
    # load(file="reduction_parameters_HydScen_SGloc.rData")
    # new file: dynamic selection of vertical reduction depth possible
    load(file="reduction_parameters_HydScen_SGloc_limitedDepth.rData")
    # filter parameterset for chosen settings
    redParam_hydScen_SGloc = dplyr::filter(reduction_parameters,
                                                  Scenario == Scenario,
                                                  SGlocation == SGlocation,
                                                  verticalLimit == VertLimit)
    factor_ts = MeanSoilMoisture %>%
                dplyr::mutate(fac1 = redParam_hydScen_SGloc$Intercept + redParam_hydScen_SGloc$Slope * value) %>%
                dplyr::select(-value)

    return(factor_ts)
}

