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
            dat_tsf = 7,
            ...
){
    # load gravity observations
    gravity_obs_data = read_data(gravity_obs, input_dir, dat_tsf = dat_tsf)
    # subtract mean value 
    gravity_obs_data$value = gravity_obs_data$value - mean(gravity_obs_data$value, na.rm = T)
    # change column name to enable join
    gravity_obs_data = gravity_obs_data %>%
        dplyr::mutate(obs = value) %>%
        dplyr::select(-value)
    # subtract first value (first timestep)
    # in order to create relative mass change response
    # gravity_below$value = gravity_below$value - gravity_below$value[1]
    # subtract mean value 
    gravity_below$value = gravity_below$value - mean(gravity_below$value, na.rm = T)

    # reduce observed data by gravity reponse from below building
    gravity_reduced = dplyr::inner_join(gravity_obs_data, gravity_below) %>%
        dplyr::mutate(value = obs - value) %>%
        dplyr::select(-obs)

    return(gravity_reduced)
}
