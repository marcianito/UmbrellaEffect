#' @title Correct the gravity component grid for the SG building foundation
#'
#' @description In order to set all mass variations to zero where a firm and fix building can be found,
#' this function modifies the gravity component grid in the way that all gravity components laying physically
#' within the building foundations, are set to zero.
#'
#' @param gravity_gcomp Data.frame, 4 columns: spatial information (x,y,z) and corresponding gravity component (value).
#' @param building_foundation Data.frame, 3 columns: spatial information (x,y,z) about the stucture of the foundation of the SG building.
#' 
#' @return Returns a data.frame, consisting of 4 columns: 3 dimensions in space (x,y,Z) and value, holding the gravity component data.
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

correct_SGbuilding_foundation = function(
            gravity_gcomp,
            building_foundation
){
    ## removing SGPILLAR, WALLS, BASEMENT of building
    # add column for use of filtering
    building_foundation_house = cbind(building_foundation, house = T)
    gcomp_grid3d = left_join(gcomp_grid3d, building_foundation_house) %>%
    			mutate(value = ifelse(is.na(house) == T, value, 0)) %>%
    			select(-house)
    
    # remove duplicated entries
    dups = which(duplicated(gcomp_grid3d[,1:3]) == T)
    gcomp_grid3d = gcomp_grid3d[-dups,]

    return(gcomp_grid3d)
}

