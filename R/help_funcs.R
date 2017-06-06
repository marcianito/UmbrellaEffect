#' @title Convert DEM to data.frame
#'
#' @description interpolate DEM to necessary grid of corresponding layer
#'
#' @param grid_cords coordinates of the original hydrus mesh.
#' @param data_input data to be interpolated. in fomat of the hydrus mesh.
#' @param grid_discr vector of discretization of new grid in (x,y,z).
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

convert_demtodf = function(dem, dem.info){
# convert (melt) dem
dem_x = seq(dem.info[3,1],length.out=length(dem[1,]),by=dem.info[5,1])
dem_y = seq(dem.info[4,1],length.out=length(dem[,1]),by=dem.info[5,1])
colnames(dem) = dem_x
rownames(dem) = rev(dem_y)
dem.melt =  melt(as.matrix(dem))#, id.vars=
colnames(dem.melt) = c("y","x","z")
return(dem.melt)
}

#' @title Convert zoo-object to data frame
#'
#' @description Converting zoo-object to data frame with proper indexing of time column and column names
#'
#' @param value must be a zoo-object
#' @references Marvin Reich (2017), mreich@@posteo.de
# ' @examples examples MISSING
#' 
zootodf <- function(value) {
    df <- as.data.frame(value)
    df$time <- index(value) #create a Date column
    rownames(df) <- NULL #so row names not filled with dates
    df <- df[,c(ncol(df), 1:(ncol(df)-1))] #reorder columns so Date first
    return(df)
}

#' @title Calculate number of decimal places in a number
#'
#' @description 
#'
#' @param x Numeric, the number to check decimal places.
#' @details missing
#' @references Marvin Reich (2017), mreich@@posteo.de
#' @examples missing

decimalplaces <- function(x) {
    if ((x %% 1) != 0) {
        nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
    } else {
        return(0)
    }
}

##############################
## some random


# grid2d_cords = read_2dgrid(fol_path)
# save(grid2d_cords, file=paste(homedir,fol,"output/Hydrus2d_grid.rdata",sep=""))

#########################
#########################

# ## vertical SM data for grids outside of SGbuilding (only 1D [vertical] information necessary)
# # radius 50 m & 300 m
# grid_r50m_r300m = c(0.1,0.1); depth_r50m_r300m = c(0,5) #until and including z=0.1
# nmax=1
# SMgrid_r50m_r300m = nodestogrid_2d(grid2d_cords, SMdata_hydrus2d, grid_r50m_r300m, depth_r50m_r300m,nmax)
# #get mean over all x locations
# SMvertical_r50m_r300m = group_by(SMgrid_r50m_r300m, Timestep,Depth) %>%
#             filter(Timestep > 0) %>%
#             summarize(value = mean(value, na.rm=T))
# #colnames(SMvertical_r50m_r300m)[3] = "value"
# save(SMvertical_r50m_r300m, file=paste(homedir,fol,"output/SMdata_TS_1Dvertical_interpolated_allLayers_r50m_r300m.rdata",sep=""))
# 
# # radius 5000 m
# # average 1d grid with higher resolution to 1 m vertical resolution
# transformGrid = data.frame(Depth = unique(SMvertical_r50m_r300m$Depth), layer = c(rep(1,11),rep(2,10),rep(3,10),rep(4,10),rep(5,10)))
# SMvertical_r5000m = inner_join(SMvertical_r50m_r300m, transformGrid, by="Depth") %>%
#         group_by(Timestep, layer) %>%
#         summarize(value=mean(value, na.rm=T))
# colnames(SMvertical_r5000m)[2] = "Depth"

