#' @title Calculate the size of the SG building
#'
#' @description Calculates the area of the SG building
#'
#' @param Bd_x Vector containing min and max x coordinate of the building.
#' @param Bd_y Vector containing min and max y coordinate of the building.
#' 
#' @return A numeric value is returned, providing the size of the
#'         area of the building. Units are as input values are supplied.
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

BdSize = function(
            Bd_x,
            Bd_y
){
    buildingSize = (max(Bd_x) - min(Bd_x)) * (max(Bd_y) - min(Bd_y))
    return(buildingSize)
}

