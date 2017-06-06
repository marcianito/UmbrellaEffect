#' @title Read Digital Elevation Models 
#'
#' @description Read DEM's and output information files for further usage. 
#'
#' @param dempath Path of input DEM-file. 
#' @param filename DEM-file in .acs-format (see details)
#' 
#' @details So far only it is only possible to read DEM-file in .asc arquitecture format.
#' @details Outputs files: data.DEM, info.DEM (as data.frames).
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing
 
read_dem = function(dempath, filename){
    info.DEM = read.table(file=paste(dempath,filename,sep=""),nrows=6, row.names=1, header=F, sep="", stringsAsFactors=F, colClasses=c("character","numeric"))
    nodata = info.DEM[6,1]
    data.DEM = read.table(file=paste(dempath,filename,sep=""), sep="", dec=".", skip=6, stringsAsFactors=F, na.strings=nodata, header=F)
    dem.info <<- info.DEM; dem<<- data.DEM
    sucess = "output: dem.info and dem"
    return(sucess)
}
