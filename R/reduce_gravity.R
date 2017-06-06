#' @title Reduce observed gravity data
#'
#' @description Reduces (corrects) observed gravity data by mass variations occuring and originating below SG observatory buildings.
#'
#' @param gravity_obs Vector, containing the filename of the observed gravity data, which is to be reduced.
#' @param gravity_below Data.frame, containing time (datatime) and gravity response data (value).
#' @param input_dir Vector, containing the directory of the gravity_data-file.
#' 
#' @return Returns a data.frame, consisting of a time series (2 columns, time info and data)
#' of the reduced gravity response. This means that observed gravity data was corrected (subtracted) by
#' the response due to mass variation from below observatory buildings.
#' This way improving (or cleaning) the observed gravtiy signal from near-field local hydrology due to the umbrella effect.
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

reduce_gravity = function(
            gravity_obs,
            gravity_below,
            input_dir,
            ...
){
    # load gravity observations
    gravity_obs_data = read_data(gravity_obs, input_dir)

    # reduce observed data by gravity reponse from below building
    gravity_reduced = dplyr::inner_join(gravity_obs_data, gravity_below) %>%
        dplyr::mutate(value = obs - value) %>%
        dplyr::select(-obs)

    return(gravity_reduced)
}
