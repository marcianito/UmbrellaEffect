#' @title Plot multiply time series with one same basis
#'
#' @description Plot all gravity time series used and calulated throught this reduction procedure in one joined plot.
#'
#' @param gravity_obs Character string with the filename of the gravity observation data to be used.
#' @param gravity_outside Data.frame, containing gravity response corresponding to mass variations from outside of the buildling.
#' It has 2 columns (datatime and value).
#' @param gravity_below Data.frame, containing gravity response corresponding to mass variations from below of the buildling.
#' It has 2 columns (datatime and value).
#' @param gravity_reduced Data.frame, containing the reduced gravity data (gravity observation data corrected for mass variations from below the SG building).
#' It has 2 columns (datatime and value).
#' @param input_dir Character string, containing the path to input file "gravity_obs".
#' @param output_dir Character string, containing the path to the directory where to save the plot.
#' 
#' @return Returns NULL as object. It creates a plot of all input time series.
#' This plot is shown on screen as well as saved to the specified output folder.
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

plot_ts_data = function(
            gravity_obs = gravityObservations_input_file,
            gravity_outside = gravity_response_outside_building,
            gravity_below = gravity_response_below_building,
            gravity_reduced = gravity_data_reduced,
            input_dir = dir_input,
            output_dir = dir_output
){
    ## read in routine for gravity obs data
    # same as in function "reduce_gravity"
    # this time resulting in columns $datetime, $value

    # combine datasets
    gravity_data_all = rbind(
            cbind(gravity_obs, Source = "Observed"),
            cbind(gravity_outside, Source = "Outside"),
            cbind(gravity_below, Source = "Below"),
            cbind(gravity_reduced, Source = "Reduced")
    )

    # plot
    gravity_ts.gg = ggplot(data = gravity_data_all,
                           aes(x = datetime, y = value, colour = Source)) +
                    geom_line() + 
                    ylab("Gravity response [nm/s²]") + 
                    xlab("Time") + 
                    colour_scale_manual(values = viridis(5)[4]) + 
                    theme(
                          legend.position = "right"
                          )

    # save plot
    png(file = paste0(output_dir, "Gravity_timeseries.png"),
                      width = 600,
                      height = 400,
                      res = 150)
    print(gravity_ts.gg)
    dev.off()
    # print on screen
    print(gravity_ts.gg)
    # return NULL
    return(NULL)
}

