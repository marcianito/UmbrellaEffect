#' @title Receive reduction factor based on the dimension of the SG building
#'
#' @description test
#'
#' @param MeanSoilMoisture Data.frame, time series of the mean soil moisture content in the complete space outside of the observatory building.
#' It needs 2 columns: one for time information (datetime) and one with the soil moisture data (value).
#' 
#' @return Returns a numerical values, which represents the correction factor based on the dimension of the SG building.
#' This can later be used to adjust the modeled gravity response from outside of the building.
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

reduction_BdSize = function(
            SG_BdSize
){
    # load reduction parameters
    load(file="reduction_parameters_BuildingSize.rData")
    load(file="/home/mreich/Dokumente/written/ResearchRepos/UmbrellaEffect/data/reduction_parameters_BuildingSize.rData")

    factor_num = reduction_parameters_building$Intercept + reduction_parameters_building$Slope * exp(1 / SG_BdSize)

    return(factor_num)
}

