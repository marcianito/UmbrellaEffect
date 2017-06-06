#' @title Read input data files
#'
#' @description Reads data files for usage of input data.
#' Various file formats are supported.
#'
#' @param data_in Vector, containing the filename of the data to read.
#' @param input_dir Vector, containing the directory of the file.
#' @param spat_col Vector, defines the columns of the statial coordinates in the stucture: vector(x, y, z).
#' If not all dimensions are supplied, the corresponding entry has to be NA.
#' @param dat_col Numeric, containing the number of the column which holds the measurment data.
#' @param ... additional parameters for reading .rData lines. (sep, dec, etc.). 
#' 
#' @return Returns a data.frame, consisting of a time series (time info and data).
#' 
#' @details If no columns are specified for spatial coordinates or data,
#' the simplest 1d case is assumed, providing time, z and data.
#' So far, input data is supported in the following file formats:
#' .rData, .csv, .tsf.
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

read_data = function(
                    data_in,
                    data_dir,
                    spat_col = c(NA, NA, 2),
                    dat_col = 3,
                    ...
){
    # .rData, .csv, .tsf, 
    # find input file type (extension)

    file_type = strsplit(data_in, ".", fixed = TRUE)
    f_type = file_type[[1]][length(file_type[[1]])]
    switch(f_type,
           rData = {
              read_rData(data_in, data_dir, spat_col, dat_col, ...)
           },
           csv = {
              read_csv(data_in, data_dir, spat_col, dat_col, ...)
           },
           tsf = {
              read_tsf(data_in, data_dir, dat_col, ...)
           }
    )

}

#' @title Read .rData files
#'
#' @description Reads .rData files and stores it as a data.frame.
#'
#' @param data_in Vector, containing the filename of the data to read.
#' @param input_dir Vector, containing the directory of the file.
#' @param ... additional parameters for reading .rData lines. (sep, dec, etc.). 
#' See documentation of read.table() for a detailed list of options.
#' 
#' @return Returns a data.frame, consisting of a time series (2 columns, time info and data).
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

read_rData = function(
                    data_in,
                    data_dir,
                    spat_col = c(NA, NA, 2),
                    dat_col = 3,
                    ...
){
    # read .rData file
    # data_in = "SMdata_TS_1d.rData"
    # data_dir = "~/temp/UM/"
    file_in = load(file = paste0(data_dir, data_in))
    data_read = get(file_in)

    # correct naming of columns
    colnames(data_read)[1] = "datetime"
    colnames(data_read)[spat_col[3]] = "z"
    colnames(data_read)[dat_col] = "value"
    # if x and y coordinates are provided
    if(!is.na(spat_col[1])) colnames(data_read)[spat_col[1]] = "x"
    if(!is.na(spat_col[2])) colnames(data_read)[spat_col[2]] = "y"

    return(data_read)
}

#' @title Read .csv files
#'
#' @description Reads .csv files and stores it as a data.frame.
#'
#' @param data_in Vector, containing the filename of the data to read.
#' @param input_dir Vector, containing the directory of the file.
#' @param ... additional parameters for reading .csv lines. (sep, dec, etc.). 
#' See documentation of read.csv() for a detailed list of options.
#' 
#' @return Returns a data.frame, consisting of a time series (2 columns, time info and data).
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

read_csv = function(
                    data_in,
                    data_dir,
                    spat_col = c(NA, NA, 2),
                    dat_col = 3,
                    ...
){
    # read .csv file
    data_read = read.csv(file = paste0(data_dir, data_in), header = T)

    # correct naming of columns
    colnames(data_read)[1] = "datetime"
    colnames(data_read)[spat_col[3]] = "z"
    colnames(data_read)[dat_col] = "value"
    # if x and y coordinates are provided
    if(!is.na(spat_col[1])) colnames(data_read)[spat_col[1]] = "x"
    if(!is.na(spat_col[2])) colnames(data_read)[spat_col[2]] = "y"

    return(data_read)
}

#' @title Read .tsf files
#'
#' @description Reads .tsf files and stores it as a data.frame.
#'
#' @param data_in Vector, containing the filename of the data to read.
#' Please note that the .tsf file should end with regular data and NOT a special statement.
#' @param input_dir Vector, containing the directory of the file.
#' 
#' @return Returns a data.frame, consisting of a time series (2 columns, time info and data).
#' 
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

read_tsf = function(
                    data_in,
                    data_dir,
                    dat_col = 7,
                    ...
){
    # determine number of rows to skip
    # this is defined by hitting "[DATA]" in the .tsf file
    data_lines = readLines(con = paste0(data_dir, data_in), n = 100)
    line_start = which(data_lines %in% "[DATA]")
    
    # read in .tsf file
    # and concatenate date information
    data_input = read.table(file = paste0(data_dir, data_in), skip=line_start, header=F, sep="", dec=".", na.strings=9999.999)
    data_input$datetime = strptime(paste(data_input[,3],"/",data_input[,2],"/",data_input[,1],"-", data_input[,4],":", data_input[,5], sep=""),  "%d/%m/%Y-%H:%M")
    data_input$datetime = as.POSIXct(data_input$datetime)
    # select column of gravity time series data
    data_read = data.frame(
                            datatime = data_input$datetime,
                            obs = data_input[,dat_col]
    )
    # transform 9999.999 values into NA, necessary!?
    # lysidata_raw[which(lysidata_raw[,7:10] == 9999.999,arr.ind=T)] = NA
    
    return(data_read)
}

