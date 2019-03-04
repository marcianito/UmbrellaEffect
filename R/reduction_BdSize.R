#' @title Receive reduction factor based on the dimension of the SG building
#'
#' @description test
#'
#' @param MeanSoilMoisture Data.frame, time series of the mean soil moisture content in the complete space outside of the observatory building.
#' It needs 2 columns: one for time information (datetime) and one with the soil moisture data (value).
#' @param VertLimit Numeric value of the vertical extent to be considered for the reduction. Values are in meters and have to be positive.
#' 
#' @return Returns a numerical values, which represents the correction factor based on the dimension of the SG building.
#' This can later be used to adjust the modeled gravity response from outside of the building.
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing
#' @export

reduction_BdSize = function(
            SG_BdSize,
            VertLimit
){
    # load reduction parameters
    # old file: for complete 5 m vertical model extent
    # load(file="reduction_parameters_BuildingSize.rData")
    # new file: dynamic selection of vertical reduction depth possible
    load(file="reduction_parameters_BuildingSize_limitedDepth.rData")

    Reg_params = reduction_parameters_building %>%
        dplyr::filter(verticalLimit == VertLimit)
    
    factor_num = Reg_params$A + Reg_params$B * exp(1 / SG_BdSize)

    return(factor_num)
}

