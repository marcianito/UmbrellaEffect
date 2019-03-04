######################################################################
### functions concerning hydrological software "hydrus" (1D,2D,3D) ###
######################################################################

#' @title Data preparation for HYDRUS
#'
#' @description dem TS 
#'
#' @param inputpath DEM-data (realtive or absolute paths)
#' @param inputfile data.frame containing information about the DEM (row wise): columns, rows, starting x, starting y, lengths of one cellsize
#' @param dz vector containing the thickness of each subsurface layer
#' @param limit.area vector of coordinates declaring the area of interest within the inputfile (DEM). Structure: X(left), Y(top), X(right), Y(bottom).
#' 
#' @details Both DEM-data and its information file can be created using readDEM(). 
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing 

hydrus_dem <- function(inputpath,inputfile,dz,limit.area){
  ##loading until package status..
      library(reshape2)
  #adjust DEM for hydrus input (local coordinates, vertical layer adjustments, etc...)
      readDEM(inputpath, inputfile)
      discre(dem.info[1,],dem.info[2,],dem.info[3,],dem.info[4,],dem.info[5,])
      #change DEM cols & rows to coordinates
      colnames(dem) = xp
      rownames(dem) = yp
      YXZ=melt(as.matrix(dem))
      XYZ = cbind(X=YXZ[,2],Y=YXZ[,1],Z=YXZ[,3])
      #XYZ.hydrus = cbind(X=(YXZ[,2]-min(YXZ[,2])),Y=(YXZ[,1]-min(YXZ[,1])),Z=YXZ[,3]) #real XYZ coordinates
      #XYZlayers.hydrus = cbind(X=(YXZ[,2]-min(YXZ[,2])),Y=(YXZ[,1]-min(YXZ[,1])),Z0=(min(YXZ[,3])-dzmax),Z1=(min(YXZ[,3])-dzmax), Zmax=YXZ[,3]) #modified coordinates for hydrus input
      outside.area=which((XYZ[,1]<limit.area[1] & XYZ[,2]>limit.area[2]) | (XYZ[,1]>limit.area[3] & XYZ[,2]<limit.area[4]))#, arr.ind=T)
      outside.area=which((XYZ[,1]<limit.area[1] | XYZ[,2]>limit.area[2]) | (XYZ[,1]>limit.area[3] | XYZ[,2]<limit.area[4]))#, arr.ind=T)
      XYZ.filter = XYZ[-outside.area,]
      #outside.X=which(XYZ[,1]<limit.area[1] | XYZ[,1]>limit.area[3])
      #outside.Y=which(XYZ[,2]>limit.area[2] | XYZ[,2]<limit.area[4])
      dzmax=sum(dz)
      #XYZlayers.hydrus = cbind(X=(XYZ.filter[,1]-min(XYZ[,1])),Y=(XYZ.filter[,2]-min(XYZ[,2])),Zbase=(min(XYZ.filter[,3])-dzmax),ZGW=(min(XYZ.filter[,3])-dzmax)+dz[4],Zsap=(min(XYZ.filter[,3])-dzmax)+sum(dz[3:4]),Zsoil=(min(XYZ.filter[,3])-dzmax)+dz[2:4],Zmax=XYZ.filter[,3]) #modified coordinates for hydrus input
      XYZlayers.hydrus = cbind(X=(XYZ.filter[,1]-min(XYZ[,1])),Y=(XYZ.filter[,2]-min(XYZ[,2])),Zbase=(min(XYZ.filter[,3])-dzmax))
      for(i in length(dz):1){
       XYZlayers.hydrus = cbind(XYZlayers.hydrus, (XYZ.filter[,3]-sum(dz[1:i])))
        #,ZGWlow=XYZ.filter[,3]-sum(dz[1:4]),ZGW=XYZ.filter[,3]-sum(dz[1:3]),Zsap=XYZ.filter[,3]-sum(dz[1:2]),Zsoil=XYZ.filter[,3]-dz[1],Zmax=XYZ.filter[,3]) #modified coordinates for hydrus input 
      }
      XYZlayers.hydrus = cbind(XYZlayers.hydrus, Zmax=XYZ.filter[,3])
      return(XYZlayers.hydrus)      
}

#' @title Data preparation for HYDRUS
#'
#' @description Create evapotranspiration time series 
#'
#' @param date.in starting date (POSIXct)
#' @param date.out ending date (POSIXct)
#' @param freq frequency of data. Possible values are "day", "hour", "min", "sec".
#' @param aszoo output either in data.frame (aszoo=F) or zoo (aszoo=T). Default value is FALSE.
#' 
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de

hydrus_eto <- function(date.in, date.out, freq="hour", aszoo=F){
      ## später mit readwettzell daten laden/kombinieren !!!
      #ETO from lysimeter data
      library(zoo); Sys.setenv(TZ = "GMT")
      load(file="/home/mreich/server/hygra/DataWettzell/Lysimeter/data_complete/eto_lysimeter_agg.rdata")
      #load(file="/home/mreich/server/hygra/DataWettzell/Lysimeter/data_complete/eto_lysimeter_raw.rdata")
      eto = eto.agg[which(index(eto.agg)==date.in):which(index(eto.agg)==date.out)]
      eto.agg = aggregate(eto, function(x) as.POSIXct(trunc(x, freq),ts="GMT"),sum, na.rm=T) #autocorrection of non-rounded timestampts
      if(aszoo==F){
      return(zootodf(eto.agg))
      }
      else {return(eto.agg)}
}

#' @title Data preparation for HYDRUS
#'
#' @description Create lysimeter drainage time series
#'
#' @param date.in starting date (POSIXct)
#' @param date.out ending date (POSIXct)
#' @param freq frequency of data. Possible values are "day", "hour", "min", "sec".
#' @param area bottom surface area of the lysimeter. In m².
#' @param aszoo output either in data.frame (aszoo=F) or zoo (aszoo=T). Default value is FALSE.
#'
#' @details The outout will be given in meters drained water column per timestep (differences).
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de

hydrus_lysidrain_ben <- function(date.in, date.out, freq="hour", area=1,winlength=10, aszoo=F){
      library(zoo); Sys.setenv(TZ = "GMT")
# load(file="/home/mreich/server/hygra/DataWettzell/Lysimeter/data_complete/weight_lysimeter_raw.rdata")
      load(file="/home/mreich/server/hygra/DataWettzell/Lysimeter/data_complete/weight_lysimeter20082011_raw.rdata")
      #adjust script to drainage from lysimeter
      lysimeter = weight_lysimeter_20082011[which(index(weight_lysimeter_20082011)==date.in):which(index(weight_lysimeter_20082011)==date.out)]
      # drain_kg = lysimeter$SiwaWaage #weight in[kg] of drained water
      drain_kg = lysimeter$SiwaWaageCUM #weight in[kg] of drained water
      #change to data.frame
      drain_kg_df = zootodf(drain_kg)
      #####
      #set parameters
      # lag = laglength/length(drain_kg_df$value)
      # robustness = 5
      win=rep(1/(winlength+1),(winlength+1))
      ## remove NA-values
      drain_kg_df$value[which(is.na(drain_kg_df$value)==T)] = 0
      ## smooth using rlowess
      # drain_kg_df$siwa_smooth = lowess(drain_kg_df$value, f=lag, iter=robustness)
      # drain_kg_df$siwa_smooth = lowess(drain_kg_df$time, drain_kg_df$value, f=lag, iter=robustness)$y
      # drain_kg_df$siwa_mv = rollmean(drain_kg_df$value, f=lag, iter=robustness)
      drain_kg_df$siwa_mv = stats::filter(coredata(drain_kg_df$value),win,sides=2)
      ## interpolate missing time steps
      # drain_kg_df$siwa = na.approx(drain_kg_df$siwa_mv, na.rm=T)
      ######
      #shift data one timestep
      drain_kg_df$siwaSHIFT = c(drain_kg_df$siwa_mv[-1],NA) #create new column with values shifted one timestep
      #calculate difference (t+1) - t
      drain_kg_df$siwaDIF = -1*(drain_kg_df$siwaSHIFT - drain_kg_df$siwa_mv)
      #eliminate NA-values
      # drain_kg_df$SiwaWaageDIF[which(is.na(drain_kg_df$SiwaWaageDIF)==T)] = 0
      #correct the times of pumping out the water tank
			# drain_kg_df[which(drain_kg_df$SiwaWaageDIF < -0.5)] = 0
      #change timeseries back to zoo
      drain_dif_kg = zoo(drain_kg_df$siwaDIF, order.by=drain_kg_df$time)
      #convert to [m] water column
      density_water = 999.9720 #[kg/m³]
      lysi_drain = (drain_dif_kg/lysi_area)*(1/density_water)
      # drain.agg = lysi_drain
      drain.agg = aggregate(lysi_drain, function(x) as.POSIXct(trunc(x, freq),ts="GMT"),sum, na.rm=T) #autocorrection of non-rounded timestampts
      if(aszoo==F){
      return(zootodf(drain.agg))
      }
      else {return(drain.agg)}
}

#' @title Data preparation for HYDRUS
#'
#' @description Create lysimeter drainage time series
#'
#' @param date.in starting date (POSIXct)
#' @param date.out ending date (POSIXct)
#' @param freq frequency of data. Possible values are "day", "hour", "min", "sec".
#' @param area bottom surface area of the lysimeter. In m².
#' @param aszoo output either in data.frame (aszoo=F) or zoo (aszoo=T). Default value is FALSE.
#'
#' @details The outout will be given in meters drained water column per timestep (differences).
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de

hydrus_lysidrain <- function(date.in, date.out, freq="hour", area=1, aszoo=F){
      library(zoo); Sys.setenv(TZ = "GMT")
      load(file="/home/mreich/server/hygra/DataWettzell/Lysimeter/data_complete/weight_lysimeter_raw.rdata")
      #adjust script to drainage from lysimeter
      lysimeter = weight.lysimeter[which(index(weight.lysimeter)==date.in):which(index(weight.lysimeter)==date.out)]
      drain_kg = lysimeter$SiwaWaage #weight in[kg] of drained water
      #aggregate to hourly data before processing
      # drain_kg = aggregate(drain_kg, function(x) as.POSIXct(trunc(x, freq),ts="GMT"),mean, na.rm=T) #autocorrection of non-rounded timestampts
      #change to data.frame
      drain_kg_df = zootodf(drain_kg)
      #shift data one timestep
      drain_kg_df$SiwaWaageSHIFT = c(drain_kg_df$value[-1],NA) #create new column with values shifted one timestep
      #calculate difference (t+1) - t
      drain_kg_df$SiwaWaageDIF = -1*(drain_kg_df$SiwaWaageSHIFT - drain_kg_df$value)
      #eliminate tank pumping times 
      drain_kg_df$SiwaWaageDIF[which(drain_kg_df$SiwaWaageDIF > 0.22)] = 0 #0.5; opt=0.22
      #eliminate unrealistic tank increases (too much drainage?) 
      # check value for reasonability !!!
      drain_kg_df$SiwaWaageDIF[which(drain_kg_df$SiwaWaageDIF < -0.22)] = 0 #-0.1; opt=-0.22
      #eliminate NA-values
      drain_kg_df$SiwaWaageDIF[which(is.na(drain_kg_df$SiwaWaageDIF)==T)] = 0
      #correct the times of pumping out the water tank
			# drain_kg_df[which(drain_kg_df$SiwaWaageDIF < -0.5)] = 0
      #change timeseries back to zoo
      drain_dif_kg = zoo(drain_kg_df$SiwaWaageDIF, order.by=drain_kg_df$time)
      #convert to [m] water column
      density_water = 999.9720 #[kg/m³]
      lysi_drain = (drain_dif_kg/lysi_area)*(1/density_water)
      # drain.agg = lysi_drain
      drain.agg = aggregate(lysi_drain, function(x) as.POSIXct(trunc(x, freq),ts="GMT"),sum, na.rm=T) #autocorrection of non-rounded timestampts
      if(aszoo==F){
      return(zootodf(drain.agg))
      }
      else {return(drain.agg)}
}

#' @title Data preparation for HYDRUS
#'
#' @description Create precipitation time series
#'
#' @param date.in starting date (POSIXct)
#' @param date.out ending date (POSIXct)
#' @param freq frequency of data. Possible values are "day", "hour", "min", "sec".
#' @param aszoo output either in data.frame (aszoo=F) or zoo (aszoo=T). Default value is FALSE.
#'
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de

hydrus_precip <- function(date.in, date.out, freq="hour", aszoo=F){
      library(zoo); Sys.setenv(TZ = "GMT")
      load(file="/home/mreich/server/hygra/DataWettzell/Climate/30min_wettzell/clima30min_raw.rdata")
      clima.all = clima.raw[which(index(clima.raw)==date.in):which(index(clima.raw)==date.out)]
      precip = (clima.all$Prec_Sen1)
      precip.agg = aggregate(precip, function(x) as.POSIXct(trunc(x, freq),ts="GMT"),sum, na.rm=T) #autocorrection of non-rounded timestampts
      if(aszoo==F){
      return(zootodf(precip.agg))
      }
      else {return(precip.agg)}
}

#' @title Data preparation for HYDRUS
#'
#' @description groundwater TS 
#'
#' @param date.in starting date (POSIXct)
#' @param date.out ending date (POSIXct)
#' @param name name of groundwater well to use. Possible so far "BK1", "BK2", "BK3", "BK14".
#' @param depthLB depths of lower model boundary. All groundwater level values will be subtracted from this value.
#' @param freq frequency of data. Possible values are "day", "hour", "min", "sec".
#' @param aszoo output either in data.frame (aszoo=F) or zoo (aszoo=T). Default value is FALSE.
#' @param aprox if the time series has some NA-values, they can be approximized. Default value is TRUE.
#' 
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de

hydrus_gw <- function(date.in, date.out, name, depthLB, freq="hour",aszoo=F){
      library(zoo); Sys.setenv(TZ = "GMT") #;library(xts)
      gw.name = load(file=paste("/home/mreich/server/hygra/DataWettzell/Groundwater/GW_Pegel_complete/",name,"_corrected.rdata",sep=""))
      gw.agg = aggregate(get(gw.name), function(x) as.POSIXct(trunc(x, freq),ts="GMT"),mean, na.rm=T) #autocorrection of non-rounded timestampts
      #gw.agg = period.apply(get(gw.name), endpoints(get(gw.name), freq), mean, na.rm=T)
      date.in.freq = as.POSIXct(trunc(date.in, freq),ts="GMT")
      date.out.freq = as.POSIXct(trunc(date.out, freq),ts="GMT")
      #gw = gw.agg[which(index(gw.agg)==date.in):which(index(gw.agg)==date.out)]
      gw = gw.agg[which(index(gw.agg)==date.in.freq):which(index(gw.agg)==date.out.freq)]
      gw.corrected = depthLB - gw #recalculate so output TS is in GW-table in reference to lower model boundary (LB)
      if(aszoo==F){
      return(zootodf(gw.corrected))
      }
      else {return(gw.corrected)}
}

#' @title Data preparation for HYDRUS
#'
#' @description soil moisture TS 
#'
#' @param date.in starting date (POSIXct)
#' @param date.out ending date (POSIXct)
#' @param name name of groundwater well to use. Possible so far "BK1", "BK2", "BK3", "BK14".
#' @param depthLB depths of lower model boundary. All groundwater level values will be subtracted from this value.
#' @param freq frequency of data. Possible values are "day", "hour", "min", "sec".
#' @param aszoo output either in data.frame (aszoo=F) or zoo (aszoo=T). Default value is FALSE.
#' @param aprox if the time series has some NA-values, they can be approximized. Default value is TRUE.
#' 
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de

hydrus_sm <- function(date.in, date.out, name, datacol, freq="hour",aszoo=F){
      library(zoo); Sys.setenv(TZ = "GMT") #;library(xts)
      sm.name = load(file=paste("/home/mreich/server/hygra/DataWettzell/SoilMoisture/Cluster_Data/Data_filtered/",name,"_filtered_6hourmean.rdata",sep=""))
      sm.agg = aggregate(get(sm.name)[,datacol], function(x) as.POSIXct(trunc(x, freq),ts="GMT"),mean, na.rm=T) #autocorrection of non-rounded timestampts
      date.in.freq = as.POSIXct(trunc(date.in, freq),ts="GMT")
      date.out.freq = as.POSIXct(trunc(date.out, freq),ts="GMT")
      #gw = gw.agg[which(index(gw.agg)==date.in):which(index(gw.agg)==date.out)]
      sm = sm.agg[which(index(sm.agg)==date.in.freq):which(index(sm.agg)==date.out.freq)]
      if(aszoo==F){
      return(zootodf(sm))
      }
      else {return(sm)}
}

#' @title Read Hydrus 1D output data
#'
#' @description test
#'
#' @param folder path to hydrus project folder
#' @param vertical_nodes mumber of nodes used in the vertical model discretisation
#' @param timesteps number of modeled timesteps
#' @param t_printout number of model printout times
#' @param plotting indicate if standard parameters should be plotted (TRUE); default is FALSE
#' @param timestamps vector of timestamps of the modeled timeseries. if not provided output time will be in counts of timesteps
#' 
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing

read_hydrus1d <- function(folder, timesteps, vertical_nodes, t_printout, plotting=F, timestamps, tprint_show=10){
  # set working directory to hydrus project folder
 setwd(folder) 
 ##read output files tlevel
 num_lines = length(readLines("T_Level.out"))
 tlevel = read.table(file="T_Level.out", header=F, skip=9, sep="",nrows=(num_lines-10), dec=".")
 colnames(tlevel) = c("Time","rTop","rRoot","vTop","vRoot","vBot","sum(rTop)","sum(rRoot)","sum(vTop)","sum(vRoot)","sum(vBot)","hTop","hRoot","hBot","RunOff","sum(RunOff)","Volume","sum(Infil)","sum(Evap)","TLevel","Cum(WTrans)","SnowLayer")
 #remove duplicated data lines (interpolation artefacts probably)
 tlevel$Timecut=round(tlevel$Time,0)
 tlevel = tlevel[!duplicated(tlevel$Timecut),]

system("sed -n '/^[[:space:]]*[T]/p' Nod_Inf.out > bash_time.out")
system("sed -n '/^[[:space:]]*[0-9]/p' Nod_Inf.out > bash_data.out")
system("sed -n '/^[[:space:]]*[0-9]/p' Profile.out > bash_profile.out")

printout = read.table(file="bash_time.out", dec=".")
time_print = printout[,2]
data_all = read.table(file="bash_data.out", dec=".")
profile_data = read.table(file="bash_profile.out", dec=".")
grid_cords = profile_data[,2]

th_h_profiles = data.frame(Time=rep(time_print, each=vertical_nodes),
		      #x=rep(grid_cords$x, each=length(time_print)),
		      #y=rep(grid_cords$y, each=length(time_print)),
		      Depth=rep(grid_cords, length(time_print)),
		      Moisture=data_all[,4],Head=data_all[,3],
		      K=data_all[,5], Flux=data_all[,7])
## add real timestamp information
 #th_h_profiles = mutate(th_h_profiles, datetime=as.POSIXct(Time, origin=timestamps[1], tz="GMT"))
## store variable "outside" of function for further usage
 th_h_profiles <<- th_h_profiles
 print("Output stored in th_h_profiles")

 #for(i in 2:(t_printout+1)){
 #if(i == (t_printout+1)) break
 #balanceTIME = read.table(file="Balance.out", header=F, skip=((i-1)*15+20), sep="", nrows=1, dec=".")
 #balanceVOL = read.table(file="Balance.out", header=F, skip=((i-1)*15+30), sep="", nrows=1, dec=".")
 #balancePER = read.table(file="Balance.out", header=F, skip=((i-1)*15+31), sep="", nrows=1, dec=".")
 #balanceDATA = cbind(balanceTIME[,3], balanceVOL[,3], balancePER[,3])
 #balance = rbind(balance, balanceDATA)
 #}
 
 #colnames(nodal_info) = c("Time","Node","Depth","Head","Moisture","K","C","Flux","Sink","Kappa","v/KsTop","Temp")
 #colnames(balance) = c("Time","WaterBalanceVolume","WaterBalancePercent")
 tlevel = mutate(tlevel, datetime=timestamps)
 tlevel <<- tlevel
 #nodal_info <<- nodal_info
 #balance <<- balance
 #print("Output stored in variables tlevel, nodal_info and balance")

 #plotting
 if(plotting){
 tlevel.plot = melt(tlevel, id.vars=c("datetime","Time"), measure.vars=colnames(tlevel)[2:22], variable.name = "Parameters")
 #nodal.plot = melt(nodal_info, id.vars=c("Time", "Depth"), measure.vars=colnames(nodal_info)[-1:-3], variable.name = "Parameters")
 nodal.plot = melt(th_h_profiles, id.vars=c("Time", "Depth"), measure.vars=c("Moisture","Head","K","Flux"), variable.name = "Parameters")

 printtimes = unique(nodal.plot$Time) 
 printdates = c(datetime[1],datetime[printtimes])
 prints = round(seq(1,t_printout, length.out=tprint_show),0)
 printtimes.gg = printtimes[prints]
 printdates.gg = printdates[prints]
 dates_times = data.frame(Time=printtimes, Dates=format(printdates, "%Y-%m-%d"))
 #dates_times = data.frame(Time=printtimes, Dates=printdates)

 #xcolors = colorRampPalette(c("blue","red", "orange","darkgreen"))(t_printout+1)
 xcolors = colorRampPalette(c("blue","red", "orange","darkgreen"))(length(printdates.gg))
 #  xcolors = colorRampPalette(c("blue","green"))(t_printout+1)

 #  balance.plot = melt(balance, id.vars="Time", measure.vars=colnames(balance)[-1])
 #subsetting to choose standard parameters
 choose.params.tlevel = c("vTop", "vRoot", "vBot", "hTop", "hRoot", "hBot")
 units = data.frame(type=c("Flux","Head"), Unit=c("[m/h]","[m]"))
 cols = data.frame(Parameters = choose.params.tlevel, Boundary = c("Top","Root","Bottom"))
 type = data.frame(Parameters = c("vTop","vBot","hTop","hBot"), type=c("Flux","Flux","Head","Head"))
 # tlevel.plot.filter = 	filter(tlevel.plot, grepl("vTop|vBot|hTop|hBot",Parameters)) %>%
 #                         filter(Parameters != "sum(vTop)" & Parameters != "sum(vRoot)" & Parameters != "sum(vBot)") %>%
 #                         mutate(type="Flux") %>%
 #                         mutate_if(Parameters == "hTop" | Parameters == "hRoot" | Parameters == "hBot", type=="Head") %>%
 #                         inner_join(units,by="type") %>%
 #                        mutate(typeUnits = paste(type,Unit)) %>%
 #                        inner_join(cols, by="Parameters")
 # browser()			
 tlevel.plot.filter = filter(tlevel.plot, grepl("vTop|vBot|hTop|hBot",Parameters)) %>%
 			inner_join(type,by="Parameters") %>%
 			inner_join(units,by="type") %>%
			mutate(typeUnits = paste(type,Unit)) %>%
			inner_join(cols, by="Parameters")

			# tlevel.plot.filter$ParamUnits = factor(tlevel.plot.filter$ParamUnits,levels=c("Head [m]","Moisture [%]","K [m/h]","Flux [m/h]"))
 tlevel.plot.filter$Boundary = factor(tlevel.plot.filter$Boundary,levels=c("Top","Root","Bottom"))
 #  tlevel.gg = ggplot(tlevel.plot.filter, aes(x=datetime, y=value, colour=Boundary)) + geom_line() + xlab("") + ylab ("") + facet_grid(typeUnits ~ ., scale="free_y") + theme(axis.text.x=element_text(colour=xcolors, face="bold"),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + geom_vline(xintercept = as.numeric(printdates), colour="grey") + scale_x_datetime(breaks=printdates, labels=date_format("%d-%m-%Y")) 
 #  tlevel.gg = ggplot(tlevel.plot.filter, aes(x=datetime, y=value)) + geom_line(aes(linetype=Boundary)) + xlab("") + ylab ("") + facet_grid(typeUnits ~ ., scale="free_y") + theme(axis.text.x=element_text(colour=xcolors, face="bold"),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + geom_vline(xintercept = as.numeric(printdates),colour="grey") + scale_x_datetime(breaks=printdates, labels=date_format("%d-%m-%Y")) 
 tlevel.gg = ggplot(tlevel.plot.filter, aes(x=datetime, y=value)) + geom_line(aes(linetype=Boundary)) + xlab("") + ylab ("") + 
	 facet_grid(typeUnits ~ ., scale="free_y") + 
	 theme(axis.text.x=element_text(face="bold"),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + 
	 #geom_vline(xintercept = as.numeric(printdates),colour=xcolors, size=2, alpha=0.5) + 
	 #scale_x_datetime(breaks=printdates, labels=date_format("%d-%m-%Y")) 
	 geom_vline(xintercept = as.numeric(printdates.gg),colour=rep(xcolors,2), size=2, alpha=0.5) + 
	 scale_x_datetime(breaks=printdates.gg, labels=date_format("%d-%m-%Y")) 
 
 choose.params.nodal = c("Head", "Moisture", "K", "Flux")
 units = data.frame(Parameters=choose.params.nodal, Unit=c("[m]", "[%]", "[m/h]", "[m/h]"))
 #  nodal.plot.filter = 	filter(nodal.plot, Parameters == choose.params.nodal) %>%
 nodal.plot.filter = 	filter(nodal.plot, grepl("Head|Moisture|K|Flux",Parameters)) %>%
			transform(Hours=as.factor(Time)) %>%
			transform(Days=as.factor(Time/24)) %>%
 			inner_join(units,by="Parameters") %>%
			mutate(ParamUnits = paste(Parameters,Unit)) %>%
 			inner_join(dates_times,by="Time") 
			
 nodal.plot.filter$ParamUnits = factor(nodal.plot.filter$ParamUnits,levels=c("Head [m]","Moisture [%]","K [m/h]","Flux [m/h]"))
 # force value to be numeric, in case a character unwanted changed class of value
 nodal.plot.filter$value = as.numeric(nodal.plot.filter$value)
 #  nodal.gg = ggplot(nodal.plot.filter, aes(x=value, y=Depth, colour=Dates, group=Dates)) + geom_path() + xlab("") + ylab ("Depth [m]") + facet_grid(.~ParamUnits, scale="free_x") + theme(legend.position = "bottom",legend.text=element_text(size=10),legend.title=element_blank())
 nodal.gg = ggplot(filter(nodal.plot.filter, Time == printtimes.gg), aes(x=value, y=Depth, colour=Dates, group=Dates)) + geom_path() + xlab("") + ylab ("Depth [m]") + 
	 	facet_grid(.~ParamUnits, scale="free_x") + theme(legend.position = "none") + 
		ylim(max(grid_cords),0) +
		scale_color_manual(values = xcolors)# + 
		#scale_x_continuous(breaks=pretty_breaks(n=6))

 #merge both plots together
 results_plot.gg = grid.arrange(tlevel.gg,nodal.gg)
 return(results_plot.gg)
}
 return(print("Output generated without a plot"))
}

#' @title Read Hydrus 2D output data (rectangular mesh)
#'
#' @description test
#'
#' @param folder foldername of project to read
#' @param vertical_nodes number of nodes used in the vertical model discretisation
#' @param horizontal_nodes number of nodes used in the horizontal model discretisation
#' @param timesteps number of modeled timesteps
#' @param t_printout number of model printout times
#' @param plotting indicate if standard parameters should be plotted (TRUE); default is FALSE
#' @param timestamps vector of timestamps of the modeled timeseries. if not provided output time will be in counts of timesteps
#' @param px number of horizontal node, where the vertical profile should be analyzed
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing

read_hydrus2d_rect <- function(folder, timesteps, vertical_nodes, horizontal_nodes, t_printout, plotting=F, timestamps, px){
	# library(dplyr)
	# library(dplyrExtras)

 path_hydrus = "/home/mreich/server/marvin_reich/hydrus2d/" 
 setwd(paste(path_hydrus, folder, "/", sep="")) 
 #read hydrus 2d output files
 cord_depthIN = read.table("MESHTRIA.TXT", header=F, skip=1 ,sep="", nrows=vertical_nodes*horizontal_nodes, dec=".") #read z coordinates of grid
 cord_depth = unique(cord_depthIN[,3]);
 cord_depth = cord_depth - max(cord_depth)
 h_num_lines = length(readLines("h_Mean.out"))
 v_num_lines = length(readLines("v_Mean.out"))
 h_Mean = read.table(file="h_Mean.out", header=F, skip=6, sep="",nrows=(h_num_lines-7), dec=".")
 v_Mean = read.table(file="v_Mean.out", header=F, skip=13, sep="",nrows=(v_num_lines-14), dec=".")
 colnames(h_Mean) = c("Time","hAtm","hRoot","hKode3","hKode1","hSeep","hKode5","hKode6","hKode7","hKode8","hKode9")
 colnames(v_Mean) = c("Time","rAtm","rRoot","vAtm","vRoot","vKode3","vKode1","vSeep","vDrain","vBottom","vKode7","vKode8","vKode9","RunOff","Evapor","Infiltr","SnowLayer")
 #remove duplicated data lines (interpolation artefacts probably)
 h_Mean$Timecut=round(h_Mean$Time,0)
 v_Mean$Timecut=round(v_Mean$Time,0)
 h_Mean = h_Mean[!duplicated(h_Mean$Timecut),]
 v_Mean = v_Mean[!duplicated(v_Mean$Timecut),]
 #joint datasets h_Mean and v_Mean
 h_v_Mean = inner_join(h_Mean, v_Mean, by="Time")
 h_v_Mean = mutate(h_v_Mean, datetime=timestamps) #add real dates along with timesteps

 #generate sequence of vertical profiles to read out; each vertical profile needs their own sequence!!
 vert1 = seq(px,by=horizontal_nodes,length.out=vertical_nodes)
 
 th_h_profiles = data.frame()
 for(i in 0:t_printout){
 TH_infoTIME = read.table(file="TH.TXT", header=F, skip=1 + i*(ceiling(horizontal_nodes*vertical_nodes/10) + 3), sep="", nrows=1, dec=".")[,3]
 TH_infoIN = scan(file="TH.TXT", skip=3+ i*(ceiling(horizontal_nodes*vertical_nodes/10) + 3), sep="", nmax=vertical_nodes*horizontal_nodes, dec=".")
 H_infoTIME = read.table(file="H.TXT", header=F, skip=1+ i*(ceiling(horizontal_nodes*vertical_nodes/10) + 3), sep="", nrows=1, dec=".")[,3]
 H_infoIN = scan(file="H.TXT", skip=3+ i*(ceiling(horizontal_nodes*vertical_nodes/10) + 3), sep="", nmax=vertical_nodes*horizontal_nodes, dec=".")
 #fetch times and values of each corresponding vertical profile! 
 TH_profile1 = TH_infoIN[vert1]
 H_profile1 = H_infoIN[vert1]
 profiles = data.frame(Time=H_infoTIME, x= px, z= (1:vertical_nodes), Depth= cord_depth, Head= H_profile1, Moisture=TH_profile1)
 th_h_profiles = rbind(th_h_profiles,profiles)
 }

 balanceTIME = read.table(file="Balance.out", header=F, skip=21, sep="", nrows=1, dec=".")
 balanceVOL = read.table(file="Balance.out", header=F, skip=29, sep="", nrows=1, dec=".")
 balancePER = read.table(file="Balance.out", header=F, skip=30, sep="", nrows=1, dec=".")
 balance = cbind(balanceTIME[,3], balanceVOL[,3], balancePER[,3])

 for(i in 2:(t_printout+1)){
 if(i == (t_printout+1)) break
 balanceTIME = read.table(file="Balance.out", header=F, skip=((i-1)*13+21), sep="", nrows=1, dec=".")
 balanceVOL = read.table(file="Balance.out", header=F, skip=((i-1)*13+29), sep="", nrows=1, dec=".")
 balancePER = read.table(file="Balance.out", header=F, skip=((i-1)*13+30), sep="", nrows=1, dec=".")
 balanceDATA = cbind(balanceTIME[,3], balanceVOL[,3], balancePER[,3])
 balance = rbind(balance, balanceDATA)
 }
 colnames(balance) = c("Time","WaterBalanceVolume","WaterBalancePercent")
 h_Mean <<- h_Mean
 v_Mean <<- v_Mean
 th_h_profiles <<- th_h_profiles 
 balance <<- balance

 #plotting
 if(plotting){
# library(ggplot2)
# library(reshape2)
# library(dplyrExtras)
# library(grid)
# library(gridExtra)
# library(scales)
 h_v_Mean.plot = melt(h_v_Mean, id.vars=c("datetime","Time"), measure.vars=colnames(h_v_Mean)[2:22], variable.name = "Parameters")
 nodal.plot = melt(th_h_profiles, id.vars=c("Time", "Depth"), measure.vars=colnames(th_h_profiles)[-1:-4], variable.name = "Parameters")

 printtimes = unique(th_h_profiles$Time)
 #  printtimes = unique(balance[,1]) 
 printdates = c(datetime[1],datetime[printtimes])
 #  printdates = c(datetime[printtimes+1], datetime[timesteps])
 dates_times = data.frame(Time=printtimes, Dates=format(printdates, "%Y-%m-%d:%H"))

 xcolors = colorRampPalette(c("blue","red", "orange","darkgreen"))(t_printout+1)
 #  xcolors = colorRampPalette(c("blue","green"))(t_printout+1)

 #  balance.plot = melt(balance, id.vars="Time", measure.vars=colnames(balance)[-1])
 #subsetting to choose standard parameters
 choose.params.h_v_Mean = c("vAtm", "vRoot", "vKode3", "hAtm", "hRoot", "hKode3")
 units = data.frame(type=c("Flux","Head"), Unit=c("[m²/h]","[m]"))
 type = data.frame(Parameters = c("vAtm","vKode3","hAtm","hKode3"), type=c("Flux","Flux","Head","Head"))
 cols = data.frame(Parameters = choose.params.h_v_Mean, Boundary = c("Top","Root","Bottom"))
 #  h_v_Mean.plot.filter = 	filter(h_v_Mean.plot, Parameters == choose.params.h_v_Mean) %>%
 # browser()
 # h_v_Mean.plot.filter = filter(h_v_Mean.plot, grepl("vAtm|vKode3|hAtm|hKode3",Parameters)) %>%
 #                         filter(Parameters != "sum(vAtm)" & Parameters != "sum(vRoot)" & Parameters != "sum(vKode3)") %>%
 #                         mutate(type="Flux") %>%
 #                         mutate_if(Parameters == "hAtm" | Parameters == "hRoot" | Parameters == "hKode3", type=="Head") %>%
 #                         inner_join(units,by="type") %>%
 #                        mutate(typeUnits = paste(type,Unit)) %>%
 #                        inner_join(cols, by="Parameters")
			
 h_v_Mean.plot.filter = filter(h_v_Mean.plot, grepl("vAtm|vKode3|hAtm|hKode3",Parameters)) %>%
 			inner_join(type,by="Parameters") %>%
 			inner_join(units,by="type") %>%
			mutate(typeUnits = paste(type,Unit)) %>%
			inner_join(cols, by="Parameters")

			# h_v_Mean.plot.filter$ParamUnits = factor(h_v_Mean.plot.filter$ParamUnits,levels=c("Head [m]","Moisture [%]","K [m/h]","Flux [m/h]"))
 h_v_Mean.plot.filter$Boundary = factor(h_v_Mean.plot.filter$Boundary,levels=c("Top","Root","Bottom"))

 #  h_v_Mean.gg = ggplot(h_v_Mean.plot.filter, aes(x=datetime, y=value, colour=Boundary)) + geom_line() + xlab("") + ylab ("") + facet_grid(typeUnits ~ ., scale="free_y") + theme(axis.text.x=element_text(colour=xcolors, face="bold"),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + geom_vline(xintercept = as.numeric(printdates), colour="grey") + scale_x_datetime(breaks=printdates, labels=date_format("%d-%m-%Y")) 
 #  h_v_Mean.gg = ggplot(h_v_Mean.plot.filter, aes(x=datetime, y=value)) + geom_line(aes(linetype=Boundary)) + xlab("") + ylab ("") + facet_grid(typeUnits ~ ., scale="free_y") + theme(axis.text.x=element_text(colour=xcolors, face="bold"),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + geom_vline(xintercept = as.numeric(printdates),colour="grey") + scale_x_datetime(breaks=printdates, labels=date_format("%d-%m-%Y")) 
 h_v_Mean.gg = ggplot(h_v_Mean.plot.filter, aes(x=datetime, y=value)) + geom_line(aes(linetype=Boundary)) + xlab("") + ylab ("") + facet_grid(typeUnits ~ ., scale="free_y") + theme(axis.text.x=element_text(face="bold"),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + geom_vline(xintercept = as.numeric(printdates),colour=xcolors, size=2, alpha=0.5) + scale_x_datetime(breaks=printdates, labels=date_format("%d-%m-%Y")) 
 
 choose.params.nodal = c("Head", "Moisture")
 units = data.frame(Parameters=choose.params.nodal, Unit=c("[m]", "[%]"))
 #  nodal.plot.filter = 	filter(nodal.plot, Parameters == choose.params.nodal) %>%
 nodal.plot.filter = 	filter(nodal.plot, grepl("Head|Moisture",Parameters)) %>%
			transform(Hours=as.factor(Time)) %>%
			transform(Days=as.factor(Time/24)) %>%
 			inner_join(units,by="Parameters") %>%
			mutate(ParamUnits = paste(Parameters,Unit)) %>%
 			inner_join(dates_times,by="Time") 
			
 nodal.plot.filter$ParamUnits = factor(nodal.plot.filter$ParamUnits,levels=c("Head [m]","Moisture [%]"))

 #  nodal.gg = ggplot(nodal.plot.filter, aes(x=value, y=Depth, colour=Dates, group=Dates)) + geom_path() + xlab("") + ylab ("Depth [m]") + facet_grid(.~ParamUnits, scale="free_x") + theme(legend.position = "bottom",legend.text=element_text(size=10),legend.title=element_blank())
 nodal.gg = ggplot(nodal.plot.filter, aes(x=value, y=Depth, colour=Dates, group=Dates)) + geom_path() + xlab("") + ylab ("Depth [m]") + facet_grid(.~ParamUnits, scale="free_x") + theme(legend.position = "none") + scale_color_manual(values = xcolors)
 #merge both plots together
 grid.arrange(h_v_Mean.gg,nodal.gg,main=paste(folder,": vertical profile at x =", px, sep=" "))
}
 return(print("Output stored in variables h_Mean, v_Mean, th_h_profiles and balance"))
}

#' @title Read Hydrus 2D output data
#'
#' @description test
#'
#' @param folder foldername of project to read
#' @param timesteps number of modeled timesteps
#' @param t_printout string of dates of time series where informations are printed
#' @param plotting indicate if standard parameters should be plotted (TRUE); default is FALSE
#' @param timestamps vector of timestamps of the modeled timeseries. if not provided output time will be in counts of timesteps
#' @param px number of horizontal node, where the vertical profile should be analyzed
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing

# read_hydrus2d_irregular <- function(folder, timesteps, vertical_nodes, horizontal_nodes, t_printout, plotting=F, timestamps, px){
read_hydrus2d <- function(folder, timesteps, factorforprintout, plotting=F, timestamps, profile_loc, vertical=T){
	# library(dplyr)
	# library(dplyrExtras)

 #path_hydrus = "/home/mreich/server/marvin_reich/hydrus2d/experiments/" 
 #path_hydrus = "/home/mreich/server/hydro72/hydrus2d/experiments/" 
 path_hydrus = "/home/mreich/server/sec54c139/Documents/hydrus2d/experiments/" 
 setwd(paste(path_hydrus, folder, "/", sep="")) 
 #setwd(paste("./", folder, "/", sep="")) 
 #read hydrus 2d output files
 nodes_meta = read.table("MESHTRIA.TXT", header=F, sep="", nrows=1, dec=".") #read grid meta data
 nodes_max = nodes_meta[1,2]
 cord_depthIN = read.table("MESHTRIA.TXT", header=F, skip=1 ,sep="", nrows=nodes_max, dec=".") #read z coordinates of grid
 #order_x = cord_depthIN[,2]
 #order_z = cord_depthIN[,3]
 grid_cords = data.frame(x=cord_depthIN[,2] ,z=cord_depthIN[,3])
 # cord_depth = unique(cord_depthIN[,3]);
 # cord_depth = cord_depth - max(cord_depth)
 h_num_lines = length(readLines("h_Mean.out"))
 v_num_lines = length(readLines("v_Mean.out"))
 h_Mean = read.table(file="h_Mean.out", header=F, skip=6, sep="",nrows=(h_num_lines-7), dec=".")
 v_Mean = read.table(file="v_Mean.out", header=F, skip=13, sep="",nrows=(v_num_lines-14), dec=".")
 colnames(h_Mean) = c("Time","hAtm","hRoot","hKode3","hKode1","hSeep","hKode5","hKode6","hKode7","hKode8","hKode9")
 colnames(v_Mean) = c("Time","rAtm","rRoot","vAtm","vRoot","vKode3","vKode1","vSeep","vDrain","vBottom","vKode7","vKode8","vKode9","RunOff","Evapor","Infiltr","SnowLayer")
 #remove duplicated data lines (interpolation artefacts probably)
 h_Mean$Timecut=round(h_Mean$Time,0)
 v_Mean$Timecut=round(v_Mean$Time,0)
 h_Mean = h_Mean[!duplicated(h_Mean$Timecut),]
 v_Mean = v_Mean[!duplicated(v_Mean$Timecut),]
 #joint datasets h_Mean and v_Mean
 h_v_Mean = inner_join(h_Mean[-1,], v_Mean[-1,], by="Time")
 h_v_Mean = mutate(h_v_Mean, datetime=timestamps) #add real dates along with timesteps

 #generate sequence of vertical profiles to read out; each vertical profile needs their own sequence!!
 #vertiklales profil bei x=2 & x=10m
 #vert1 = seq(profile_loc,by=horizontal_nodes,length.out=vertical_nodes)
 
system("sed -n '/^[[:space:]]*[T]/p' TH.TXT > bash_time.out")
system("sed -n '/^[[:space:]]*[0-9]/p' TH.TXT > bash_dataTH.out")
system("sed -n '/^[[:space:]]*-*[0-9]/p' H.TXT > bash_dataH.out")

printout = read.table(file="bash_time.out", dec=".")
data.th = scan(file="bash_dataTH.out", dec=".")
data.h = scan(file="bash_dataH.out", dec=".")
time_print = printout[,3]
th_h_profiles = data.frame(Time=rep(time_print, each=nodes_max),
		      x=rep(grid_cords$x, each=length(time_print)),
		      #y=rep(grid_cords$y, each=length(time_print)),
		      Depth=rep(grid_cords$z, each=length(time_print)),
		      Moisture=data.th,Head=data.h)
##interpolate hydrus mesh output to regular grid
#library(fields)
library(gstat)
#generate regular-spaced grid
grid.x <- seq(min(grid_cords$x), max(grid_cords$x), by=1)
#grid.y <- seq(min(grid_cords$y), max(grid_cords$y), by=1)
grid.z <- seq(min(grid_cords$z), max(grid_cords$z), by=0.05)
#grid.xyz <- expand.grid(x=grid.x, y=grid.y, Depth=grid.z)
grid.xz <- expand.grid(x=grid.x, Depth=grid.z)

#interpolate and "stack" for each timestep
nodal.plot=data.frame()
for(i in time_print){
th_h_data = filter(th_h_profiles, Time==i) #filter for one timestep
#theta
idw.gstat <- gstat(formula = Moisture ~ 1, locations = ~ x + Depth, data = th_h_data, nmax = 15, set = list(idp = 2))
theta_interpolated <- predict(idw.gstat, grid.xz)
data_interpolated = cbind(Time = i, theta_interpolated[,-4])
colnames(data_interpolated)[4] = "Moisture"
theta_t = melt(data_interpolated, id.vars=c("Time","Depth","x"), measure.vars=c("Moisture"), variable.name = "Parameters")
#head
idw.gstat <- gstat(formula = Head ~ 1, locations = ~ x + Depth, data = th_h_data, nmax = 15, set = list(idp = 2))
head_interpolated <- predict(idw.gstat, grid.xz)
data_interpolated = cbind(Time = i, head_interpolated[,-4])
colnames(data_interpolated)[4] = "Head"
head_t = melt(data_interpolated, id.vars=c("Time","Depth","x"), measure.vars=c("Head"), variable.name = "Parameters")
#join data together
nodal.plot = rbind(nodal.plot, theta_t, head_t)
}

 balanceTIME = read.table(file="Balance.out", header=F, skip=21, sep="", nrows=1, dec=".")
 balanceVOL = read.table(file="Balance.out", header=F, skip=29, sep="", nrows=1, dec=".")
 balancePER = read.table(file="Balance.out", header=F, skip=30, sep="", nrows=1, dec=".")
 balance = cbind(balanceTIME[,3], balanceVOL[,3], balancePER[,3])
for(i in 2:length(time_print)){
 if(i == length(time_print)) break
 balanceTIME = read.table(file="Balance.out", header=F, skip=((i-1)*13+21), sep="", nrows=1, dec=".")
 balanceVOL = read.table(file="Balance.out", header=F, skip=((i-1)*13+29), sep="", nrows=1, dec=".")
 balancePER = read.table(file="Balance.out", header=F, skip=((i-1)*13+30), sep="", nrows=1, dec=".")
 balanceDATA = cbind(balanceTIME[,3], balanceVOL[,3], balancePER[,3])
 balance = rbind(balance, balanceDATA)
 }
 colnames(balance) = c("Time","WaterBalanceVolume","WaterBalancePercent")
 h_Mean <<- h_Mean
 v_Mean <<- v_Mean
 th_h_profiles <<- th_h_profiles 
 balance <<- balance

 #plotting
 if(plotting){
# library(ggplot2)
# library(reshape2)
# library(dplyrExtras)
# library(grid)
# library(gridExtra)
# library(scales)
 h_v_Mean.plot = melt(h_v_Mean, id.vars=c("datetime","Time"), measure.vars=colnames(h_v_Mean)[2:22], variable.name = "Parameters")
 #nodal.plot = melt(th_h_profiles, id.vars=c("Time", "Depth","x"), measure.vars=c("Moisture","Head"), variable.name = "Parameters")

 #printtimes = unique(th_h_profiles$Time)
 #printdates = c(datetime[1],datetime[printtimes])
 printdates = c(datetime[1],datetime[time_print])
 #  printdates = c(datetime[printtimes+1], datetime[timesteps])
 #dates cut to "days", better for displaying in x-axis
 #but problematic if 2 dates are selected within one day,then use second option
 #dates_times = data.frame(Time=time_print, Dates=as.POSIXct(format(printdates, "%Y-%m-%d:%H:%M")))
 dates_times = data.frame(Time=time_print, Dates=printdates)
 t_printout = time_print[seq(1, length(time_print), factorforprintout)]
 dates_times_print = filter(dates_times, Time%in%t_printout)

 xcolors = colorRampPalette(c("blue","red", "orange","darkgreen"))(length(t_printout))

 #  balance.plot = melt(balance, id.vars="Time", measure.vars=colnames(balance)[-1])
 #subsetting to choose standard parameters
 choose.params.h_v_Mean = c("vAtm", "vRoot", "vKode3", "hAtm", "hRoot", "hKode3")
 units = data.frame(type=c("Flux","Head"), Unit=c("[m²/h]","[m]"))
 type = data.frame(Parameters = c("vAtm","vKode3","hAtm","hKode3"), type=c("Flux","Flux","Head","Head"))
 cols = data.frame(Parameters = choose.params.h_v_Mean, Boundary = c("Top","Root","Bottom"))
 #  h_v_Mean.plot.filter = 	filter(h_v_Mean.plot, Parameters == choose.params.h_v_Mean) %>%
 # browser()
 # h_v_Mean.plot.filter = filter(h_v_Mean.plot, grepl("vAtm|vKode3|hAtm|hKode3",Parameters)) %>%
 #                         filter(Parameters != "sum(vAtm)" & Parameters != "sum(vRoot)" & Parameters != "sum(vKode3)") %>%
 #                         mutate(type="Flux") %>%
 #                         mutate_if(Parameters == "hAtm" | Parameters == "hRoot" | Parameters == "hKode3", type=="Head") %>%
 #                         inner_join(units,by="type") %>%
 #                        mutate(typeUnits = paste(type,Unit)) %>%
 #                        inner_join(cols, by="Parameters")
			
 h_v_Mean.plot.filter = filter(h_v_Mean.plot, grepl("vAtm|vKode3|hAtm|hKode3",Parameters)) %>%
 			inner_join(type,by="Parameters") %>%
 			inner_join(units,by="type") %>%
			mutate(typeUnits = paste(type,Unit)) %>%
			inner_join(cols, by="Parameters")

			# h_v_Mean.plot.filter$ParamUnits = factor(h_v_Mean.plot.filter$ParamUnits,levels=c("Head [m]","Moisture [%]","K [m/h]","Flux [m/h]"))
 h_v_Mean.plot.filter$Boundary = factor(h_v_Mean.plot.filter$Boundary,levels=c("Top","Root","Bottom"))

 #  h_v_Mean.gg = ggplot(h_v_Mean.plot.filter, aes(x=datetime, y=value, colour=Boundary)) + geom_line() + xlab("") + ylab ("") + facet_grid(typeUnits ~ ., scale="free_y") + theme(axis.text.x=element_text(colour=xcolors, face="bold"),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + geom_vline(xintercept = as.numeric(printdates), colour="grey") + scale_x_datetime(breaks=printdates, labels=date_format("%d-%m-%Y")) 
 #  h_v_Mean.gg = ggplot(h_v_Mean.plot.filter, aes(x=datetime, y=value)) + geom_line(aes(linetype=Boundary)) + xlab("") + ylab ("") + facet_grid(typeUnits ~ ., scale="free_y") + theme(axis.text.x=element_text(colour=xcolors, face="bold"),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + geom_vline(xintercept = as.numeric(printdates),colour="grey") + scale_x_datetime(breaks=printdates, labels=date_format("%d-%m-%Y")) 
 #h_v_Mean.gg = ggplot(h_v_Mean.plot.filter, aes(x=datetime, y=value)) + geom_line(aes(linetype=Boundary)) + xlab("") + ylab ("") + facet_grid(typeUnits ~ ., scale="free_y") + theme(axis.text.x=element_text(face="bold"),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + geom_vline(xintercept = as.numeric(printdates),colour=xcolors, size=2, alpha=0.5) + scale_x_datetime(breaks=printdates, labels=date_format("%d-%m-%Y")) 
 h_v_Mean.gg = ggplot(h_v_Mean.plot.filter, aes(x=datetime, y=value)) + geom_line(aes(linetype=Boundary)) + xlab("") + ylab ("") + facet_grid(typeUnits ~ ., scale="free_y") + theme(axis.text.x=element_text(face="bold"),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + 
	 geom_vline(xintercept = as.numeric(dates_times_print$Dates),color=rep(xcolors,2), size=2, alpha=0.5) +
	 scale_x_datetime(breaks=dates_times_print$Dates, labels=date_format("%d-%m-%Y")) 
 
 choose.params.nodal = c("Head", "Moisture")
 units = data.frame(Parameters=choose.params.nodal, Unit=c("[m]", "[%]"))
 #  nodal.plot.filter = 	filter(nodal.plot, Parameters == choose.params.nodal) %>%
 if(vertical==T){
 nodal.plot.filter = 	filter(nodal.plot, grepl("Head|Moisture",Parameters)) %>%
			filter(x==profile_loc) %>%
			filter(Time%in%t_printout) %>%
			transform(Hours=as.factor(Time)) %>%
			transform(Days=as.factor(Time/24)) %>%
 			inner_join(units,by="Parameters") %>%
			mutate(ParamUnits = paste(Parameters,Unit)) %>%
 			inner_join(dates_times,by="Time") %>%
			arrange(Time, Parameters, Depth, x)
 nodal.plot.filter$ParamUnits = factor(nodal.plot.filter$ParamUnits,levels=c("Head [m]","Moisture [%]"))
 nodal.plot.filter$Dates = as.factor(nodal.plot.filter$Dates)
 profile_type="vertical"
 nodal.gg = ggplot(nodal.plot.filter, aes(x=value, y=Depth, colour=Dates, group=Dates)) + geom_path() + xlab("") + ylab ("Depth [m]") + facet_grid(.~ParamUnits, scale="free_x") + theme(legend.position = "none") + scale_color_manual(values = xcolors)
 }else{
 nodal.plot.filter = 	filter(nodal.plot, grepl("Head|Moisture",Parameters)) %>%
			filter(Depth==profile_loc) %>%
			filter(Time%in%t_printout) %>%
			transform(Hours=as.factor(Time)) %>%
			transform(Days=as.factor(Time/24)) %>%
 			inner_join(units,by="Parameters") %>%
			mutate(ParamUnits = paste(Parameters,Unit)) %>%
 			inner_join(dates_times,by="Time") %>%
			arrange(Time, Parameters, Depth, x)
 nodal.plot.filter$ParamUnits = factor(nodal.plot.filter$ParamUnits,levels=c("Head [m]","Moisture [%]"))
 nodal.plot.filter$Dates = as.factor(nodal.plot.filter$Dates)
 profile_type="horizontal"
 nodal.gg = ggplot(nodal.plot.filter, aes(x=x, y=value, colour=Dates, group=Dates)) + geom_path() + xlab("Cross-section [m]") + ylab ("") + facet_grid(ParamUnits ~ ., scale="free_y") + theme(legend.position = "none") + scale_color_manual(values = xcolors)
 }
 #  nodal.gg = ggplot(nodal.plot.filter, aes(x=value, y=Depth, colour=Dates, group=Dates)) + geom_path() + xlab("") + ylab ("Depth [m]") + facet_grid(.~ParamUnits, scale="free_x") + theme(legend.position = "bottom",legend.text=element_text(size=10),legend.title=element_blank())
 #nodal.gg = ggplot(nodal.plot.filter, aes(x=value, y=Depth, colour=Dates, group=Dates)) + geom_path() + xlab("") + ylab ("Depth [m]") + facet_grid(.~ParamUnits, scale="free_x") + theme(legend.position = "none") + scale_color_manual(values = xcolors)
 #merge both plots together
 grid.arrange(h_v_Mean.gg,nodal.gg,top=textGrob(paste(folder,":", profile_type," profile at: ", profile_loc,"[m]", sep="")))
}
 return(print("Output stored in variables h_Mean, v_Mean, th_h_profiles and balance"))
}

#' @title Read Hydrus 3D output data
#'
#' @description test
#'
#' @param folder foldername of project to read
#' @param timesteps number of modeled timesteps
#' @param t_printout number of model printout times
#' @param plotting indicate if standard parameters should be plotted (TRUE); default is FALSE
#' @param timestamps vector of timestamps of the modeled timeseries. if not provided output time will be in counts of timesteps
#' @param px number of horizontal node, where the vertical profile should be analyzed
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing

read_hydrus3d <- function(folder, timesteps, t_printout, plotting=F, timestamps, px){
	# library(dplyr)
	# library(dplyrExtras)

 #path_hydrus = "/home/mreich/Dokumente/Wettzell/hydrologicalmodelling/hydrus3d_out/" 
 path_hydrus = "/home/mreich/server/sec54c139/Documents/hydrus3d/hydrus3D_experiments/" 
 setwd(paste(path_hydrus, folder, "/", sep="")) 
 #read hydrus 2d output files
 nodes_meta = read.table("MESHTRIA.TXT", header=F, sep="",skip=5, nrows=1, dec=".") #read grid meta data
 nodes_max = nodes_meta[1,1]
 grid_cords = read.table("MESHTRIA.TXT", header=T, skip=6 ,sep="", nrows=nodes_max, dec=".") #read z coordinates of grid

 h_num_lines = length(readLines("h_Mean.out"))
 v_num_lines = length(readLines("v_Mean.out"))
 h_Mean = read.table(file="h_Mean.out", header=F, skip=6, sep="",nrows=(h_num_lines-7), dec=".")
 v_Mean = read.table(file="v_Mean.out", header=F, skip=13, sep="",nrows=(v_num_lines-14), dec=".")
 colnames(h_Mean) = c("Time","hAtm","hRoot","hKode3","hKode1","hSeep","hKode5","hKode6","hKode7","hKode8","hKode9")
 colnames(v_Mean) = c("Time","rAtm","rRoot","vAtm","vRoot","vKode3","vKode1","vSeep","vDrain","vBottom","vKode7","vKode8","vKode9","RunOff","Evapor","Infiltr","SnowLayer")
 #remove duplicated data lines (interpolation artefacts probably)
 h_Mean$Timecut=round(h_Mean$Time,0)
 v_Mean$Timecut=round(v_Mean$Time,0)
 h_Mean = h_Mean[!duplicated(h_Mean$Timecut),]
 v_Mean = v_Mean[!duplicated(v_Mean$Timecut),]
 #joint datasets h_Mean and v_Mean
 h_v_Mean = inner_join(h_Mean, v_Mean, by="Time")
 h_v_Mean = mutate(h_v_Mean, datetime=timestamps) #add real dates along with timesteps

 #generate sequence of vertical profiles to read out; each vertical profile needs their own sequence!!
 #vertiklales profil bei x=2 & x=10m
 #vert1 = seq(px,by=horizontal_nodes,length.out=vertical_nodes)
 
system("sed -n '/^[[:space:]]*[T]/p' TH.TXT > bash_time.out")
system("sed -n '/^[[:space:]]*[0-9]/p' TH.TXT > bash_dataTH.out")
system("sed -n '/^[[:space:]]*-*[0-9]/p' H.TXT > bash_dataH.out")

printout = read.table(file="bash_time.out", dec=".")
data.th = scan(file="bash_dataTH.out", dec=".")
data.h = scan(file="bash_dataH.out", dec=".")
time_print = printout[,3]
th_h_profiles = data.frame(Time=rep(time_print, each=nodes_max),
		      x=rep(grid_cords$x, each=length(time_print)),
		      y=rep(grid_cords$y, each=length(time_print)),
		      Depth=rep(grid_cords$z, each=length(time_print)),
		      Moisture=data.th,Head=data.h)

##interpolate hydrus mesh output to regular grid
#library(fields)
library(gstat)

#generate regular-spaced grid
grid.x <- seq(min(grid_cords$x), max(grid_cords$x), by=1)
grid.y <- seq(min(grid_cords$y), max(grid_cords$y), by=1)
grid.z <- seq(min(grid_cords$z), max(grid_cords$z), by=0.05)
grid.xyz <- expand.grid(x=grid.x, y=grid.y, Depth=grid.z)

#interpolate and "stack" for each timestep
nodal.plot=data.frame()
for(i in time_print){
th_h_data = filter(th_h_profiles, Time==i) #filter for one timestep
#theta
idw.gstat <- gstat(formula = Moisture ~ 1, locations = ~ x + y + Depth, data = th_h_data, nmax = 15, set = list(idp = 2))
theta_interpolated <- predict(idw.gstat, grid.xyz)
data_interpolated = cbind(Time = i, theta_interpolated[,-5])
colnames(data_interpolated)[5] = "Moisture"
theta_t = melt(data_interpolated, id.vars=c("Time","Depth","x","y"), measure.vars=c("Moisture"), variable.name = "Parameters")
#head
idw.gstat <- gstat(formula = Head ~ 1, locations = ~ x + y + Depth, data = th_h_data, nmax = 15, set = list(idp = 2))
head_interpolated <- predict(idw.gstat, grid.xyz)
data_interpolated = cbind(Time = i, head_interpolated[,-5])
colnames(data_interpolated)[5] = "Head"
head_t = melt(data_interpolated, id.vars=c("Time","Depth","x","y"), measure.vars=c("Head"), variable.name = "Parameters")
#join data together
nodal.plot = rbind(nodal.plot, theta_t, head_t)
}

#th_h_profiles = data.frame()
 #for(i in 0:t_printout){
 #TH_infoTIME = read.table(file="TH.TXT", header=F, skip=1 + i*(ceiling(nodes_max/10) + 3), sep="", nrows=1, dec=".")[,3]
 #TH_infoIN = scan(file="TH.TXT", skip=3+ i*(ceiling(nodes_max/10) + 3), sep="", nmax=nodes_max, dec=".")
 #H_infoTIME = read.table(file="H.TXT", header=F, skip=1+ i*(ceiling(nodes_max/10) + 3), sep="", nrows=1, dec=".")[,3]
 #H_infoIN = scan(file="H.TXT", skip=3+ i*(ceiling(nodes_max/10) + 3), sep="", nmax=nodes_max, dec=".")
 #profiles = data.frame(x=order_x ,Depth=(order_z-max(order_z)), Moisture=TH_infoIN, Head=H_infoIN,Time=TH_infoTIME)
 #th_h_profiles = rbind(th_h_profiles,profiles)
 #}

 balanceTIME = read.table(file="Balance.out", header=F, skip=21, sep="", nrows=1, dec=".")
 balanceVOL = read.table(file="Balance.out", header=F, skip=29, sep="", nrows=1, dec=".")
 balancePER = read.table(file="Balance.out", header=F, skip=30, sep="", nrows=1, dec=".")
 balance = cbind(balanceTIME[,3], balanceVOL[,3], balancePER[,3])
for(i in 2:(t_printout+1)){
 if(i == (t_printout+1)) break
 balanceTIME = read.table(file="Balance.out", header=F, skip=((i-1)*13+21), sep="", nrows=1, dec=".")
 balanceVOL = read.table(file="Balance.out", header=F, skip=((i-1)*13+29), sep="", nrows=1, dec=".")
 balancePER = read.table(file="Balance.out", header=F, skip=((i-1)*13+30), sep="", nrows=1, dec=".")
 balanceDATA = cbind(balanceTIME[,3], balanceVOL[,3], balancePER[,3])
 balance = rbind(balance, balanceDATA)
 }
 colnames(balance) = c("Time","WaterBalanceVolume","WaterBalancePercent")
 h_Mean <<- h_Mean
 v_Mean <<- v_Mean
 th_h_profiles <<- th_h_profiles 
 balance <<- balance

 #plotting
 if(plotting){
# library(ggplot2)
# library(reshape2)
# library(dplyrExtras)
# library(grid)
# library(gridExtra)
# library(scales)
 h_v_Mean.plot = melt(h_v_Mean, id.vars=c("datetime","Time"), measure.vars=colnames(h_v_Mean)[2:22], variable.name = "Parameters")
 #nodal.plot = melt(th_h_profiles, id.vars=c("Time", "Depth","x","y"), measure.vars=c("Moisture","Head"), variable.name = "Parameters")

 printtimes = unique(th_h_profiles$Time)
 #printdates = c(datetime[1],datetime[printtimes])
 printdates = c(timestamps[1],timestamps[printtimes])
 dates_times = data.frame(Time=printtimes, Dates=format(printdates, "%Y-%m-%d:%H"))

 xcolors = colorRampPalette(c("blue","red", "orange","darkgreen"))(t_printout+1)
 #  xcolors = colorRampPalette(c("blue","green"))(t_printout+1)

 #  balance.plot = melt(balance, id.vars="Time", measure.vars=colnames(balance)[-1])
 #subsetting to choose standard parameters
 choose.params.h_v_Mean = c("vAtm", "vRoot", "vKode3", "hAtm", "hRoot", "hKode3")
 units = data.frame(type=c("Flux","Head"), Unit=c("[m²/h]","[m]"))
 type = data.frame(Parameters = c("vAtm","vKode3","hAtm","hKode3"), type=c("Flux","Flux","Head","Head"))
 cols = data.frame(Parameters = choose.params.h_v_Mean, Boundary = c("Top","Root","Bottom"))
 #  h_v_Mean.plot.filter = 	filter(h_v_Mean.plot, Parameters == choose.params.h_v_Mean) %>%
 # browser()
 # h_v_Mean.plot.filter = filter(h_v_Mean.plot, grepl("vAtm|vKode3|hAtm|hKode3",Parameters)) %>%
 #                         filter(Parameters != "sum(vAtm)" & Parameters != "sum(vRoot)" & Parameters != "sum(vKode3)") %>%
 #                         mutate(type="Flux") %>%
 #                         mutate_if(Parameters == "hAtm" | Parameters == "hRoot" | Parameters == "hKode3", type=="Head") %>%
 #                         inner_join(units,by="type") %>%
 #                        mutate(typeUnits = paste(type,Unit)) %>%
 #                        inner_join(cols, by="Parameters")
			
 h_v_Mean.plot.filter = filter(h_v_Mean.plot, grepl("vAtm|vKode3|hAtm|hKode3",Parameters)) %>%
 			inner_join(type,by="Parameters") %>%
 			inner_join(units,by="type") %>%
			mutate(typeUnits = paste(type,Unit)) %>%
			inner_join(cols, by="Parameters")

			# h_v_Mean.plot.filter$ParamUnits = factor(h_v_Mean.plot.filter$ParamUnits,levels=c("Head [m]","Moisture [%]","K [m/h]","Flux [m/h]"))
 h_v_Mean.plot.filter$Boundary = factor(h_v_Mean.plot.filter$Boundary,levels=c("Top","Root","Bottom"))

 #  h_v_Mean.gg = ggplot(h_v_Mean.plot.filter, aes(x=datetime, y=value, colour=Boundary)) + geom_line() + xlab("") + ylab ("") + facet_grid(typeUnits ~ ., scale="free_y") + theme(axis.text.x=element_text(colour=xcolors, face="bold"),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + geom_vline(xintercept = as.numeric(printdates), colour="grey") + scale_x_datetime(breaks=printdates, labels=date_format("%d-%m-%Y")) 
 #  h_v_Mean.gg = ggplot(h_v_Mean.plot.filter, aes(x=datetime, y=value)) + geom_line(aes(linetype=Boundary)) + xlab("") + ylab ("") + facet_grid(typeUnits ~ ., scale="free_y") + theme(axis.text.x=element_text(colour=xcolors, face="bold"),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + geom_vline(xintercept = as.numeric(printdates),colour="grey") + scale_x_datetime(breaks=printdates, labels=date_format("%d-%m-%Y")) 
 h_v_Mean.gg = ggplot(h_v_Mean.plot.filter, aes(x=datetime, y=value)) + geom_line(aes(linetype=Boundary)) + xlab("") + ylab ("") + facet_grid(typeUnits ~ ., scale="free_y") + theme(axis.text.x=element_text(face="bold"),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + geom_vline(xintercept = as.numeric(printdates),colour=xcolors, size=2, alpha=0.5) + scale_x_datetime(breaks=printdates, labels=date_format("%d-%m-%Y")) 
 
 choose.params.nodal = c("Head", "Moisture")
 units = data.frame(Parameters=choose.params.nodal, Unit=c("[m]", "[%]"))
 units$Parameters = factor(units$Parameters,levels=c("Moisture","Head"))

 #  nodal.plot.filter = 	filter(nodal.plot, Parameters == choose.params.nodal) %>%
 nodal.plot.filter = 	filter(nodal.plot, grepl("Head|Moisture",Parameters)) %>%
			filter(x==5) %>%
			filter(y==5) %>%
			#filter(x==px[1]) %>%
			#filter(y==px[1]) %>%
			transform(Hours=as.factor(Time)) %>%
			transform(Days=as.factor(Time/24)) %>%
		        inner_join(units,by="Parameters") %>%
			mutate(ParamUnits = paste(Parameters,Unit)) %>%
 			inner_join(dates_times,by="Time") %>%
			arrange(Time, Parameters, Depth, x)
			
 nodal.plot.filter$ParamUnits = factor(nodal.plot.filter$ParamUnits,levels=c("Head [m]","Moisture [%]"))

 #  nodal.gg = ggplot(nodal.plot.filter, aes(x=value, y=Depth, colour=Dates, group=Dates)) + geom_path() + xlab("") + ylab ("Depth [m]") + facet_grid(.~ParamUnits, scale="free_x") + theme(legend.position = "bottom",legend.text=element_text(size=10),legend.title=element_blank())
 nodal.gg = ggplot(nodal.plot.filter, aes(x=value, y=Depth, colour=Dates, group=Dates)) + geom_path() + xlab("") + ylab ("Depth [m]") + facet_grid(.~ParamUnits, scale="free_x") + theme(legend.position = "none") + scale_color_manual(values = xcolors)
 #merge both plots together
 grid.arrange(h_v_Mean.gg,nodal.gg,main=paste(folder,": vertical profile at x =", px, sep=" "))
}
 return(print("Output stored in variables h_Mean, v_Mean, th_h_profiles and balance"))
}

#' @title Select time series of one modeled node
#'
#' @description test
#'
#' @param nodal_info_in input datasets. should be a data.frame
#' @param loc_hor horizontal coordinate (in model coordintaes)
#' @param loc_vert vertical coordinate (in model coordinates)
#' @param sensorname name of sensor (for plotting and structuring)
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing

obsNode_hydrus = function(nodal_info_in, loc_hor, loc_vert,sensorname){
nodal.filter = 	filter(nodal_info_in, grepl("Moisture",Parameters)) %>%
			filter(x==loc_hor) %>%
			filter(Depth==loc_vert) %>%
			mutate(sensor = sensorname) %>%
			mutate(type = "modeled") %>%
 			inner_join(dates_times,by="Time") #%>%
			#arrange(Dates, sensor, type, value)
		node_out = data.frame(datetime = nodal.filter$Dates, sensor = nodal.filter$sensor, value = nodal.filter$value, type = nodal.filter$type)
		return(node_out)
}

#' @title Read Hydrus 2D observation node data 
#'
#' @description test
#'
#' @param folder foldername of project to read
#' @param realTime logical. Output in UTC (TRUE; default) or model time steps (FALSE).
#' @param startdate starting date for output in UTC. Format is POSIXct.
#' @param plotting do you want node time series plotted (default is FALSE).
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing

# read_hydrus2d_irregular <- function(folder, timesteps, vertical_nodes, horizontal_nodes, t_printout, plotting=F, timestamps, px){
read_obsNode2d <- function(folder,realTime=T,startdate,plotting=F){
 library(stringr)
 #path_hydrus = "/home/mreich/server/hydro72/hydrus2d/experiments/" 
 #path_hydrus = "/home/mreich/server/sec54c139/Documents/hydrus2d/experiments/" 
 path_hydrus = "/home/mreich/Dokumente/Wettzell/hydrologicalmodelling/hydrus2_out/"
 setwd(paste(path_hydrus, folder, "/", sep="")) 
 #read mesh
 nodes_meta = read.table("MESHTRIA.TXT", header=F, sep="", nrows=1, dec=".") #read grid meta data
 nodes_max = nodes_meta[1,2]
 nodes_cords = read.table("MESHTRIA.TXT", header=F, skip=1 ,sep="", nrows=nodes_max, dec=".") #read z coordinates of grid

 #read obsveration node output file
 num_lines = length(readLines("ObsNod.out"))
 obsNodeData = read.table(file="ObsNod.out", header=T, skip=5, sep="",nrows=(num_lines-7), dec=".")
 nodeNames.in = read.table(file="ObsNod.out", header=F, skip=3, sep="",nrows=1)
 nodeNames = vector()
 for(i in 1:(length(nodeNames.in)/2)){
 nodeNames[i] = as.numeric(str_extract(nodeNames.in[,i*2], "[0-9]+"))
 }
 #select columns
 theta_data = select(obsNodeData, contains("theta"))
 #generate zoo-TS
 if(realTime){
 obsNodeTheta = zoo(theta_data, order.by=as.POSIXct(obsNodeData$time*3600, format="%H", origin=startdate))
 }
 else{obsNodeTheta = zoo(theta_data, order.by=obsNodeData$time)}
 colnames(obsNodeTheta) = nodeNames
 #offer possibility to get output as pivot-table with node coordinates!
 obsNode.melt = melt(zootodf(obsNodeTheta), id="time", variable.name="node", value.name="theta")
 #extract node cordinates from cords
 ## dont know if this line works, but its something like this..!
 #obsNode_cords = match(nodeNames, nodes_cords)
 #obsNode.melt = inner_join(obsNode_cords, by=)
 #...
 if(plotting){
	ggplot(obsNode.melt, aes(x=time, y=theta, colour=node)) + geom_line() + facet_grid(node~.)
 }
 return(obsNodeTheta)
} # end function read obsNode data

#' @title Read Hydrus 3D observation node data 
#'
#' @description test
#'
#' @param folder foldername of project to read
#' @param realTime logical. Output in UTC (TRUE; default) or model time steps (FALSE).
#' @param startdate starting date for output in UTC. Format is POSIXct.
#' @param plotting do you want node time series plotted (default is FALSE).
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing

read_obsNode3d <- function(folder,realTime=T,startdate,plotting=F){
 #library(stringr)
 #path_hydrus = "/home/mreich/server/sec54c139/Documents/hydrus3d/hydrus3D_experiments/" 
 path_hydrus = "/home/mreich/Dokumente/Wettzell/hydrologicalmodelling/hydrus3d_out/" 
 setwd(paste(path_hydrus, folder, "/", sep="")) 
 #read mesh
 nodes_meta = read.table("MESHTRIA.TXT", header=F, sep="", nrows=1, dec=".", skip=5) #read grid meta data
 nodes_max = nodes_meta[1,1]
 nodes_cords = read.table("MESHTRIA.TXT", header=F, skip=7 ,sep="", nrows=nodes_max, dec=".") #read z coordinates of grid

 #read obsveration node output file
 num_lines = length(readLines("ObsNod.out"))
 obsNodeData = read.table(file="ObsNod.out", header=T, skip=5, sep="",nrows=(num_lines-7), dec=".")
 nodeNames.in = read.table(file="ObsNod.out", header=F, skip=3, sep="",nrows=1)
 #this 2d routine doesnt work because observation nodes in hydrus 3d are named WRONG (only ***)
 #so the alternative are dummy_names:
 #in the order of appearance, numbered trough
 #nodeNames = vector()
 #for(i in 1:(length(nodeNames.in)/2)){
 #nodeNames[i] = as.numeric(str_extract(nodeNames.in[,i*2], "[0-9]+"))
 #}
 nodeNames = seq(1,length(nodeNames.in))
 #select columns
 theta_data = select(obsNodeData, contains("theta"))
 #generate zoo-TS
 if(realTime){
 obsNodeTheta = zoo(theta_data, order.by=as.POSIXct(obsNodeData$time*3600, format="%H", origin=startdate))
 }
 else{obsNodeTheta = zoo(theta_data, order.by=obsNodeData$time)}
 colnames(obsNodeTheta) = nodeNames
 #offer possibility to get output as pivot-table with node coordinates!
 obsNode.melt = melt(zootodf(obsNodeTheta), id="time", variable.name="node", value.name="theta")
 #extract node cordinates from cords
 ## dont know if this line works, but its something like this..!
 #obsNode_cords = match(nodeNames, nodes_cords)
 #obsNode.melt = inner_join(obsNode_cords, by=)
 #...
 if(plotting){
	ggplot(obsNode.melt, aes(x=time, y=theta, colour=node)) + geom_line() + facet_grid(node~.)
 }
 return(obsNodeTheta)
} # end function read obsNode data

#' @title Read Hydrus 1D observation node data 
#'
#' @description test
#'
#' @param folder path to hydrus project folder
#' @param realTime logical. Output in UTC (TRUE; default) or model time steps (FALSE).
#' @param startdate starting date for output in UTC. Format is POSIXct.
#' @param plotting do you want node time series plotted (default is FALSE).
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing

# read_hydrus2d_irregular <- function(folder, timesteps, vertical_nodes, horizontal_nodes, t_printout, plotting=F, timestamps, px){
read_obsNode1d <- function(folder,startdate,plotting=F){
  # set working directory to hydrus project folder
  setwd(folder) 
  #read mesh
  #nodes_meta = read.table("MESHTRIA.TXT", header=F, sep="", nrows=1, dec=".") #read grid meta data
  #nodes_max = nodes_meta[1,2]
  #nodes_cords = read.table("MESHTRIA.TXT", header=F, skip=1 ,sep="", nrows=nodes_max, dec=".") #read z coordinates of grid

  #read obsveration node output file
  num_lines = length(readLines("Obs_Node.out"))
  obsNodeData = read.table(file="Obs_Node.out", header=T, skip=10, sep="",nrows=(num_lines-12), dec=".")
  nodeNames.in = read.table(file="Obs_Node.out", header=F, skip=8, sep=")",nrows=1)
  nodeNames = vector()
  #for(i in 1:(length(nodeNames.in)/2)){
  #nodeNames[i] = as.numeric(str_extract(nodeNames.in[,i*2], "[0-9]+"))
  #}
  for(i in 1:(length(nodeNames.in))){
    nodeNames[i] = as.numeric(str_extract(nodeNames.in[,i], "[0-9]+"))
  }
  nodeNames = nodeNames[-length(nodeNames)]
  #select columns
  theta_data = select(obsNodeData, contains("theta"))
  #generate zoo-TS
  obsNodeTheta = zoo(theta_data, order.by=as.POSIXct(obsNodeData$time*3600, format="%H", origin=startdate))
  ## ??
  colnames(obsNodeTheta) = nodeNames
  ## !!
  #offer possibility to get output as pivot-table with node coordinates!
  obsNode.melt = melt(zootodf(obsNodeTheta), id="time", variable.name="node", value.name="theta")
  #extract node cordinates from cords
  ## dont know if this line works, but its something like this..!
  #obsNode_cords = match(nodeNames, nodes_cords)
  #obsNode.melt = inner_join(obsNode_cords, by=)
  #...
  if(plotting){
     print(
         ggplot(obsNode.melt, aes(x=time, y=theta, colour=node)) + geom_line() + facet_grid(node~.)
         )
  }
  return(obsNodeTheta)
} # end function read obsNode data

#' @title Plot modeled nodes and observed soil moisture sensors
#'
#' @description plot time series of modeled observation nodes and measured soil moisture. This works with both hydrus 2D and 3D models.
#'
#' @param filename filename of time series of the observation nodes (see details).
#' @param normdata logical. compare absolute (FALSE) or normalized values (TRUE; default).
#' @param dim2d TRUE for 2d, FALSE for 3D input
#' @details time series of observation nodes can be extracted and saved using function read_obsNode2d() or read_obsNode3d()
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

obsNodePlot = function(filename, normdata=T, dim2d=T, time_start, time_end, filtersen=F, hydrusnodes){
#load SM; filtered at 6 hourly intervalls
#change for cluster version to maintain this file locally
#load(file="/home/mreich/server/hygra/DataWettzell/SoilMoisture/Cluster_Data/Data_filtered/SGnew_filtered_6hourmean.rdata")
#beneathBuilding = SGnew.filter[,11:18] #get sensor beneath SG building
#besidesBuilding = SGnew.filter[,19:25] #get closest sensors outside SG building
load(file="/home/mreich/server/hygra/DataWettzell/SoilMoisture/Cluster_Data/Data_2ndPhaseFiltering/SGnew_filtered_1hourLowPass.rdata")
beneathBuilding = SGnew.filter_1hlowPass[,1:8] #get sensor beneath SG building
besidesBuilding = SGnew.filter_1hlowPass[,9:15] #get closest sensors outside SG building
colnames(beneathBuilding) = c("c shallow","b shallow","d shallow","c middle high","c middle low","b deep","c deep","d deep")
colnames(besidesBuilding) = c("a shallow","a03","a04","a06","a10","a14","a deep")

SMsensors = merge.zoo(beneathBuilding,besidesBuilding, all=T, fill=NA)
if(normdata){#normalize theta data
SMsensors.norm = zoo(apply(coredata(SMsensors),2,normalize_mean), order.by=index(SMsensors))
SMsensors.melt = cbind(melt(zootodf(SMsensors.norm), id="time",variable.name="sensor"),type="observed")
}else{
SMsensors.melt = cbind(melt(zootodf(SMsensors), id="time",variable.name="sensor"),type="observed")
}

#load precipitation TS
load("/home/mreich/server/hygra/DataWettzell/Climate/30min_wettzell/clima30min_raw.rdata")
precip = clima.raw$Prec_Sen1
precip.df = zootodf(precip); colnames(precip.df)[2] = "Precipitation"
precip.melt = cbind(melt(precip.df, id="time", variable.name="sensor"), type="observed")
#adjust time series length to SM observations
precip.melt = filter(precip.melt, time > as.POSIXct("2010-03-10 16:00:00"))

##!! change file path for cluster version; exclude 2d option
#read hydrus observation nodes
if(dim2d) folpath = "/home/mreich/Dokumente/Wettzell/hydrologicalmodelling/hydrus2_out/"
else folpath = "/home/mreich/Dokumente/Wettzell/hydrologicalmodelling/hydrus3d_out/"
nodes=load(file=paste(folpath, filename, sep=""))
if(normdata){#normalize theta data
nodes.norm = zoo(apply(coredata(get(nodes)),2,normalize_mean), order.by=index(get(nodes)))
nodes.melt = cbind(melt(zootodf(nodes.norm), id="time", variable.name="node"), type="modeled")
}else{
nodes.melt = cbind(melt(zootodf(get(nodes)), id="time", variable.name="node"), type="modeled")
}
#IMPORTANT!!
#sensor names have to be in the ORDER of node number names from hydrus
#this is somewhat caotic due to internal hydrus numbering
#this is for 2D only!!
#if(dim2d) sensors_nodes = data.frame(sensor = c("b shallow","b deep","c shallow","c middle high","c middle low","d shallow","d deep","a02","a03","a04","a06","a10","a18","a14","c deep"), node = colnames(get(nodes)))
#above was OLD hydrus 2d model obs Nodes !!
#old structure
#if(dim2d) sensors_nodes = data.frame(sensor = c("a deep","c deep","d deep","a02","a04","a03","b shallow","b deep","a14","c middle low","d shallow","c shallow","c middle high","a shallow","a10"), node = colnames(get(nodes)))
#if(dim2d) sensors_nodes = data.frame(sensor = c("b shallow","b deep","c shallow","c middle high","c middle low","c deep","d shallow","d deep","a02","a03","a04","a shallow","a10","a deep","a14"), node = colnames(get(nodes)))
if(dim2d) sensors_nodes = data.frame(sensor = hydrusnodes, node = colnames(get(nodes)))
#3d
else sensors_nodes = data.frame(sensor = c("a04","a shallow","a03","a06","a10","a14","b shallow","b deep","c shallow","c middle high","c middle low","d shallow","a deep","c deep","d deep"), node = colnames(get(nodes)))

nodes.melt = inner_join(nodes.melt, sensors_nodes, by="node") %>%
		select(-node) %>%
		select(time, sensor, value, type)
#without precipitation
#SMdata = rbind(SMsensors.melt, nodes.melt)
#with precipitation
SMdata = rbind(SMsensors.melt, nodes.melt, precip.melt) %>%
	filter(time > time_start) %>%
	filter(time < time_end)
if(filtersen == T){
SMdata = filter(SMdata, sensor == "Precipitation" | sensor == "a shallow" | sensor == "b shallow" | sensor == "d shallow" | sensor == "a deep" | sensor == "d deep")
}
SMdata.gg = ggplot(SMdata, aes(x=time, y=value, colour=type)) + geom_line() + ylab("Soil moisture [%VWC]") + xlab("") +
		facet_grid(sensor~., scale="free_y")
#save plot
#png(file="/home/mreich/Dokumente/Wettzell/hydrologicalmodelling/compare/NORMALIZED_SMobs-mod_ts49025_2d.png", width=1800, height=1200, res=100)
return(SMdata.gg)
#dev.off()
}


#' @title Read Hydrus 2D theta / head data
#'
#' @description test
#'
#' @param folder foldername of project to read
#' @param timestamps vector of timestamps of the modeled timeseries.
#' @param datatype character. possible values are: "moisture", "phead".
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

read_hydrus2d_data = function(folder, timestamps, datatype="moisture"){

 #path_hydrus = "/home/mreich/Dokumente/Wettzell/hydrologicalmodelling/hydrus2_out/" 
 #path_hydrus = "/home/mreich/server/sec54c139/Documents/hydrus3d/hydrus3D_expRuns/" 
 #setwd(paste(path_hydrus, folder, "/", sep="")) 
 #setwd(paste("./", folder, "/", sep="")) 
 setwd(paste("./", folder, sep="")) 
 #read hydrus 2d output files
 nodes_meta = read.table("MESHTRIA.TXT", header=F, sep="", nrows=1, dec=".") #read grid meta data
 nodes_max = nodes_meta[1,2]
 grid_cords = read.table("MESHTRIA.TXT", header=F, skip=1 ,sep="", nrows=nodes_max, dec=".") #read z coordinates of grid
 colnames(grid_cords) = c("n","x","z")

#read output data via bash due to speed! 
system("sed -n '/^[[:space:]]*[T]/p' TH.TXT > bash_time.out")
printout = read.table(file="bash_time.out", dec=".")
time_print = printout[,3]
timestamps_select = c(timestamps[1]-3600,timestamps[time_print])
#timestamps_select = timestamps[time_print]
#select output data: soil moisture / pressure head
switch(datatype,
       moisture = {
	system("sed -n '/^[[:space:]]*[0-9]/p' TH.TXT > bash_dataTH.out")
	data.th = scan(file="bash_dataTH.out", dec=".")
	data_nodes = data.frame(Timestep=rep(time_print, each=nodes_max),
		      Time = rep(timestamps_select, each=nodes_max),
		      x=rep(grid_cords$x, length(time_print)),
		      Depth=rep(grid_cords$z, length(time_print)),
		      dataraw=data.th) #soil moisture
       },
       phead = {
	system("sed -n '/^[[:space:]]*-*[0-9]/p' H.TXT > bash_dataH.out")
	data.h = scan(file="bash_dataH.out", dec=".")
	data_nodes = data.frame(Timestep=rep(time_print, each=nodes_max),
		      Time = rep(timestamps_select, each=nodes_max),
		      x=rep(grid_cords$x, length(time_print)),
		      Depth=rep(grid_cords$z, length(time_print)),
		      dataraw=data.h) #soil moisture
       }) #end switch
return(data_nodes)
} #end function

#' @title Read Hydrus 3D theta / head data
#'
#' @description test
#'
#' @param folder foldername of project to read
#' @param timestamps vector of timestamps of the modeled timeseries.
#' @param datatype character. possible values are: "moisture", "phead".
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

read_hydrus3d_data = function(folder, timestamps, datatype="moisture"){

 #path_hydrus = "/home/mreich/Dokumente/Wettzell/hydrologicalmodelling/hydrus3d_out/" 
 path_hydrus = "/home/mreich/server/sec54c139/Documents/hydrus3d/hydrus3D_expRuns/" 
 #path_hydrus = "/home/mreich/server/sec54c139/Documents/hydrus3d/hydrus3D_experiments/" 
 setwd(paste(path_hydrus, folder, "/", sep="")) 
 #read hydrus 3d output files
 nodes_meta = read.table("MESHTRIA.TXT", header=F, sep="",skip=5, nrows=1, dec=".") #read grid meta data
 nodes_max = nodes_meta[1,1]
 grid_cords = read.table("MESHTRIA.TXT", header=T, skip=6 ,sep="", nrows=nodes_max, dec=".") #read z coordinates of grid

#read output data via bash due to speed! 
system("sed -n '/^[[:space:]]*[T]/p' TH.TXT > bash_time.out")
#system("sed -n '/^[[:space:]]*[0-9]/p' TH.TXT > bash_dataTH.out")
#system("sed -n '/^[[:space:]]*-*[0-9]/p' H.TXT > bash_dataH.out")

printout = read.table(file="bash_time.out", dec=".")
#data.th = scan(file="bash_dataTH.out", dec=".")
#data.h = scan(file="bash_dataH.out", dec=".")
time_print = printout[,3]
timestamps_select = c(timestamps[1]-3600,timestamps[time_print])
#select output data: soil moisture / pressure head
switch(datatype,
       moisture = {
	system("sed -n '/^[[:space:]]*[0-9]/p' TH.TXT > bash_dataTH.out")
	data.th = scan(file="bash_dataTH.out", dec=".")
	data_nodes = data.frame(Timestep=rep(time_print, each=nodes_max),
		      Time = rep(timestamps_select, each=nodes_max),
		      x=rep(grid_cords$x, length(time_print)),
		      y=rep(grid_cords$y, length(time_print)),
		      Depth=rep(grid_cords$z, length(time_print)),
		      dataraw=data.th) #soil moisture
       },
       phead = {
	system("sed -n '/^[[:space:]]*-*[0-9]/p' H.TXT > bash_dataH.out")
	data.h = scan(file="bash_dataH.out", dec=".")
	data_nodes = data.frame(Timestep=rep(time_print, each=nodes_max),
		      Time = rep(timestamps_select, each=nodes_max),
		      x=rep(grid_cords$x, length(time_print)),
		      y=rep(grid_cords$y, length(time_print)),
		      Depth=rep(grid_cords$z, length(time_print)),
		      dataraw=data.h) #soil moisture
       }) #end switch
return(data_nodes)
} #end function


#' @title Read Hydrus 3D theta / head data version FOR BIG DATA
#'
#' @description test
#'
#' @param filenames character vector of filenames for model output (ordered): time, soil moisture, pressure head.
#' @param timestamps vector of timestamps of the modeled timeseries.
#' @param datatype character. possible values are: "moisture", "phead".
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

read_hydrus3d_dataBD = function(filenames, timestamps, datatype="moisture"){

 #read hydrus 3d output files
 nodes_meta = read.table("MESHTRIA.TXT", header=F, sep="",skip=5, nrows=1, dec=".") #read grid meta data
 nodes_max = nodes_meta[1,1]
 grid_cords = read.table("MESHTRIA.TXT", header=T, skip=6 ,sep="", nrows=nodes_max, dec=".") #read z coordinates of grid

printout = read.table(file=filenames[1], dec=".")
time_print = printout[,3]
timestamps_select = timestamps[time_print]
#timestamps_select = c(timestamps[1]-3600,timestamps[time_print])
#select output data: soil moisture / pressure head
switch(datatype,
       moisture = {
	data.th = scan(file=filenames[2], dec=".")
	data_nodes = data.frame(Timestep=rep(time_print, each=nodes_max),
		      Time = rep(timestamps_select, each=nodes_max),
		      x=rep(grid_cords$x, length(time_print)),
		      y=rep(grid_cords$y, length(time_print)),
		      Depth=rep(grid_cords$z, length(time_print)),
		      dataraw=data.th) #soil moisture
       },
       phead = {
	data.h = scan(file=filenames[3], dec=".")
	data_nodes = data.frame(Timestep=rep(time_print, each=nodes_max),
		      Time = rep(timestamps_select, each=nodes_max),
		      x=rep(grid_cords$x, length(time_print)),
		      y=rep(grid_cords$y, length(time_print)),
		      Depth=rep(grid_cords$z, length(time_print)),
		      dataraw=data.h) #soil moisture
       }) #end switch
return(data_nodes)
} #end function

#' @title Read Hydrus 3D grid coordinates 
#'
#' @description test
#'
#' @param NA no parameter has to be inputed
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

read_3dgrid = function(){

 #read hydrus 3d output files
 nodes_meta = read.table("MESHTRIA.TXT", header=F, sep="",skip=5, nrows=1, dec=".") #read grid meta data
 nodes_max = nodes_meta[1,1]
 grid_cords = read.table("MESHTRIA.TXT", header=T, skip=6 ,sep="", nrows=nodes_max, dec=".") #read z coordinates of grid

 return(grid_cords)
} #end function

#' @title Read Hydrus 2D grid coordinates 
#'
#' @description test
#'
#' @param NA no parameter has to be inputed
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

read_2dgrid = function(path){
 setwd(path)
 #read hydrus 2d output files
 nodes_meta = read.table("MESHTRIA.TXT", header=F, sep="", nrows=1, dec=".") #read grid meta data
 nodes_max = nodes_meta[1,2]
 grid_cords = read.table("MESHTRIA.TXT", header=F, skip=1 ,sep="", nrows=nodes_max, dec=".") #read z coordinates of grid
 colnames(grid_cords)=c("n","x","z")

 return(grid_cords)
} #end function

#' @title Interpolate 3D nodes to regular grid
#'
#' @description interpolate hydrus mesh output to regular grid
#'
#' @param grid_cords coordinates of the original hydrus mesh.
#' @param data_input data to be interpolated. in fomat of the hydrus mesh.
#' @param grid_discr vector of discretization of new grid in (x,y,z).
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

nodestogrid = function(grid_cords, data_input, grid_discr, depth_split){
library(gstat)

#generate regular-spaced grid
grid.x <- seq(min(grid_cords$x), max(grid_cords$x), by=grid_discr[1])
grid.y <- seq(min(grid_cords$y), max(grid_cords$y), by=grid_discr[2])
grid.z <- seq(min(depth_split), max(depth_split), by=grid_discr[3])
#grid.z <- seq(min(grid_cords$z), max(grid_cords$z), by=grid_discr[3])
grid.xyz <- expand.grid(x=grid.x, y=grid.y, Depth=grid.z)

#filter input data for the corrsponding horizon of the model (in z-direction)
#this is done for faster computing:
#splitting up dataset allows different discretizations of each horizon
#data_input = 

#interpolate and "stack" for each timestep
data_grid=data.frame()
for(i in unique(data_input$Timestep)){
data_in = filter(data_input, Timestep==i) #filter for one timestep
#interpolate data to new grid
idw.gstat = gstat(formula = dataraw ~ 1, locations = ~ x + y + Depth, data = data_in, nmax = 10, set = list(idp = 2))
#idw.gstat = gstat(formula = dataraw ~ 1, locations = ~ x + y + Depth, data = data_in, nmax = 10, nmin=3, maxdist=1.1, set = list(idp = 2))
data_convert <- predict(idw.gstat, grid.xyz)
data_interpolated = cbind(Timestep = i, data_convert[,-5])
#colnames(data_interpolated)[5] = "data"
#data_interpolated.melt = melt(data_interpolated, id.vars=c("Timestep","Depth","x","y"), measure.vars=c("data"), variable.name = "Parameters")
data_interpolated.melt = melt(data_interpolated, id.vars=c("Timestep","Depth","x","y"))
#join data together in the loop
data_grid = rbind(data_grid, data_interpolated.melt)
}
return(data_grid)
} #end function

#' @title Interpolate 2D nodes to regular, already FIXED grid
#'
#' @description interpolate hydrus mesh output to regular grid
#'
#' @param grid_cords coordinates of the original hydrus mesh.
#' @param data_input data to be interpolated. in fomat of the hydrus mesh.
#' @param grid_discr vector of discretization of new grid in (x,y,z).
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

nodestoFIXgrid_2d = function(fixgrid, data_input,nintpmax){
library(gstat)

#interpolate and "stack" for each timestep
#data_grid=data.frame()
#interpolate data to new grid
idw.gstat = gstat(formula = value ~ 1, locations = ~ x + Depth, data = as.data.frame(data_input), nmax = nintpmax, set = list(idp = 2))
#idw.gstat = gstat(formula = dataraw ~ 1, locations = ~ x + y + Depth, data = data_in, nmax = 10, nmin=3, maxdist=1.1, set = list(idp = 2))
data_convert <- predict(idw.gstat, fixgrid)
data_interpolated = data_convert[,-4]
#colnames(data_interpolated)[5] = "data"
#data_interpolated.melt = melt(data_interpolated, id.vars=c("Timestep","Depth","x","y"), measure.vars=c("data"), variable.name = "Parameters")
data_interpolated.melt = melt(data_interpolated, id.vars=c("Depth","x"))
#join data together in the loop
data_grid = data_interpolated.melt

return(data_grid)
} #end function

#' @title Interpolate 1D nodes to regular grid
#'
#' @description interpolate hydrus mesh output to regular grid
#'
#' @param grid_cords coordinates of the original hydrus mesh.
#' @param data_input data to be interpolated. in fomat of the hydrus mesh.
#' @param grid_discr vector of discretization of new grid in (x,y,z).
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

nodestogrid_1d = function(data_input, grid_discr, depth_split){
library(gstat)
library(sp)

#generate regular-spaced grid
grid.z <- seq(min(depth_split), max(depth_split), by=grid_discr)
#grid.xz <- expand.grid(x=1, Depth=grid.z)

#filter input data for the corrsponding horizon of the model (in z-direction)
#this is done for faster computing:
#splitting up dataset allows different discretizations of each horizon

#interpolate and "stack" for each timestep
data_grid=data.frame()
j=1
for(i in unique(data_input$time)){
data_in = filter(data_input, time==i) #filter for one timestep
	#mutate(x=1)
#interpolate data to new grid
data_interpolated = data.frame(Depth = grid.z,
			       # now manually add up values and interpolate if necessary
			      value = c(data_in$value[1],data_in$value[1],data_in$value[1],data_in$value[2],data_in$value[3], # 0.0 - 0.4
					mean(c(data_in$value[3],data_in$value[4]),na.rm=T), # 0.5
					#data_in$value[4], # 0.6
					seq(data_in$value[4],data_in$value[5], length.out=5), # 0.7 - 0.9
					#data_in$value[5], # 1.0
					seq(data_in$value[5],data_in$value[6], length.out=4), # 1.1 -1.3
					#data_in$value[6], # 1.4
					seq(data_in$value[6],data_in$value[7], length.out=4), # 1.5 -1.7
					#data_in$value[7], # 1.8
					rep(data_in$value[7],32)) # 1.9 - 5.0
			      )
data_interpolated$Timestep = j

#idw.gstat = gstat(formula = value ~ 1, locations = ~ x + Depth, data = data_in, nmax = nintpmax, set = list(idp = 2))
#data_convert <- predict(idw.gstat, grid.xz)
#data_interpolated = cbind(Timestep = j, data_convert[,-4])
data_interpolated.melt = melt(data_interpolated, id.vars=c("Timestep","Depth"))
data_grid = rbind(data_grid, data_interpolated.melt)
j = j+1
}
#return(grid_interpolated)
return(data_grid)
} #end function

#' @title Interpolate 2D nodes to regular grid
#'
#' @description interpolate hydrus mesh output to regular grid
#'
#' @param grid_cords coordinates of the original hydrus mesh.
#' @param data_input data to be interpolated. in fomat of the hydrus mesh.
#' @param grid_discr vector of discretization of new grid in (x,y,z).
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

nodestogrid_2d = function(grid_cords, data_input, grid_discr, depth_split,nintpmax){
library(gstat)

#generate regular-spaced grid
grid.x <- seq(min(grid_cords$x), max(grid_cords$x), by=grid_discr[1])
grid.z <- seq(min(depth_split), max(depth_split), by=grid_discr[2])
grid.xz <- expand.grid(x=grid.x, Depth=grid.z)

#filter input data for the corrsponding horizon of the model (in z-direction)
#this is done for faster computing:
#splitting up dataset allows different discretizations of each horizon

#interpolate and "stack" for each timestep
data_grid=data.frame()
for(i in unique(data_input$Timestep)){
data_in = filter(data_input, Timestep==i) #filter for one timestep
#interpolate data to new grid
idw.gstat = gstat(formula = dataraw ~ 1, locations = ~ x + Depth, data = data_in, nmax = nintpmax, set = list(idp = 2))
data_convert <- predict(idw.gstat, grid.xz)
data_interpolated = cbind(Timestep = i, data_convert[,-4])
data_interpolated.melt = melt(data_interpolated, id.vars=c("Timestep","Depth","x"))
data_grid = rbind(data_grid, data_interpolated.melt)
}
#return(grid_interpolated)
return(data_grid)
} #end function

#' @title Interpolate 2D nodes to regular grid, using kriging
#'
#' @description interpolate hydrus mesh output to regular grid
#'
#' @param grid_cords coordinates of the original hydrus mesh.
#' @param data_input data to be interpolated. in fomat of the hydrus mesh.
#' @param grid_discr vector of discretization of new grid in (x,y,z).
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

nodestogrid_2d_krige = function(grid_cords, data_input, grid_discr, depth_split,tstart){
library(gstat)
library(spacetime)
library(sp)
library(automap)

# estimate kriging exemplarily for SMrange data (only one point in time)
# afterwards, this should be done for the complete (year 2011) timeseries of SM_obs_2011
# !!

#generate regular-spaced grid
grid.x <- seq(min(grid_cords), max(grid_cords), by=grid_discr[1])
grid.z <- seq(min(depth_split), max(depth_split), by=grid_discr[2])
grid.xz <- expand.grid(x=grid.x, z=grid.z)

griddata_hydrus_domain = cbind(grid.xz, dummyt=2, dummyb=4)
gridded(griddata_hydrus_domain) = ~ x + z
#filter input data for the corrsponding horizon of the model (in z-direction)
#this is done for faster computing:
#splitting up dataset allows different discretizations of each horizon

#interpolate and "stack" for each timestep
data_grid=data.frame()
j = tstart
for(i in unique(data_input$time)){
data_in = filter(data_input, time==i) #filter for one timestep
#data preparation
colnames(data_in)[4]="z"
coordinates(data_in)= ~x+z
#hydrus domain grid
vario = autofitVariogram(value ~ 1, data_in,
	      model = c("Sph", "Exp", "Gau", "Ste", "Mat"),
              #model = c("Mat"),
	      verbose=T)
#interpolate data to new grid
data_pred = krige(value ~ 1, data_in, griddata_hydrus_domain, vario$var_model)
# construct data.frame
data_interpolated = data.frame(Timestep = j, Depth = coordinates(data_pred)[,2] , x = coordinates(data_pred)[,1], variable = "var1.pred", value = data_pred$var1.pred)
data_grid = rbind(data_grid, data_interpolated)
j = j + 1
}
#return(grid_interpolated)
return(data_grid)
} #end function

#' @title PARALLEL Interpolate 2D nodes to regular grid
#'
#' @description interpolate hydrus mesh output to regular grid
#'
#' @param grid_cords coordinates of the original hydrus mesh.
#' @param data_input data to be interpolated. in fomat of the hydrus mesh.
#' @param grid_discr vector of discretization of new grid in (x,y,z).
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

nodestogrid_2dPP = function(grid_cords, data_input, grid_discr, depth_split,nintpmax){
library(gstat)
library(foreach)
library(doParallel, quiet=T)

#generate regular-spaced grid
grid.x <- seq(min(grid_cords$x), max(grid_cords$x), by=grid_discr[1])
grid.z <- seq(min(depth_split), max(depth_split), by=grid_discr[2])
grid.xz <- expand.grid(x=grid.x, Depth=grid.z)

#filter input data for the corrsponding horizon of the model (in z-direction)
#this is done for faster computing:
#splitting up dataset allows different discretizations of each horizon
#data_input = 

#settings for parallel computing
cores = detectCores() #detect cores
cluster = makeCluster(cores) #create cluster
registerDoParallel(cluster) #register cluster

#interpolate and "stack" for each timestep
grid_interpolated = foreach(i=unique(data_input$Timestep),.combine=rbind,.packages=c('dplyr','gstat','reshape2')) %dopar% {
data_in = filter(data_input, Timestep==i) #filter for one timestep
#interpolate data to new grid
idw.gstat = gstat(formula = dataraw ~ 1, locations = ~ x + Depth, data = data_in, nmax = nintpmax, set = list(idp = 2))
data_convert <- predict(idw.gstat, grid.xz)
data_interpolated = cbind(Timestep = i, data_convert[,-4])
data_interpolated.melt = melt(data_interpolated, id.vars=c("Timestep","Depth","x"))
}

stopCluster(cluster)
return(grid_interpolated)
} #end function

#' @title PARALLEL Interpolate 3D nodes to regular grid
#'
#' @description interpolate hydrus mesh output to regular grid
#'
#' @param grid_cords coordinates of the original hydrus mesh.
#' @param data_input data to be interpolated. in fomat of the hydrus mesh.
#' @param grid_discr vector of discretization of new grid in (x,y,z).
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

nodestogridPP = function(grid_cords, data_input, grid_discr, depth_split){
library(gstat)
library(foreach)
library(doParallel, quiet=T)

#generate regular-spaced grid
grid.x <- seq(min(grid_cords$x), max(grid_cords$x), by=grid_discr[1])
grid.y <- seq(min(grid_cords$y), max(grid_cords$y), by=grid_discr[2])
grid.z <- seq(min(depth_split), max(depth_split), by=grid_discr[3])
#grid.z <- seq(min(grid_cords$z), max(grid_cords$z), by=grid_discr[3])
grid.xyz <- expand.grid(x=grid.x, y=grid.y, Depth=grid.z)

#filter input data for the corrsponding horizon of the model (in z-direction)
#this is done for faster computing:
#splitting up dataset allows different discretizations of each horizon
#data_input = 

#settings for parallel computing
cores = detectCores() #detect cores
cluster = makeCluster(cores) #create cluster
registerDoParallel(cluster) #register cluster

#interpolate and "stack" for each timestep
grid_interpolated = foreach(i=unique(data_input$Timestep),.combine=rbind,.packages=c('dplyr','gstat','reshape2')) %dopar% {
data_in = filter(data_input, Timestep==i) #filter for one timestep
#interpolate data to new grid
idw.gstat = gstat(formula = dataraw ~ 1, locations = ~ x + y + Depth, data = data_in, nmax = 10, set = list(idp = 2))
data_convert <- predict(idw.gstat, grid.xyz)
data_interpolated = cbind(Timestep = i, data_convert[,-5])
data_interpolated.melt = melt(data_interpolated, id.vars=c("Timestep","Depth","x","y"))
}
stopCluster(cluster)
return(grid_interpolated)
} #end function

#' @title Calculate delta values of 2 grids at different time steps
#' 
#' @description 2D
#'
#' @param data1,data2,data3 input data. has to be in the format XX.
#' @param timestep vector of timesteps to be considered for calculations.
#' @param valmax,valmin values where to cut off resulting differces, leaving a data-area of special consideration
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

grid2d_difs = function(data1,data2,data3=NA,timestep, valmin, valmax){
#active for all 3 compartments
#grid_difs = function(data1,data2,data3,timestep){
	if(is.na(data3)){data3 = data.frame(Timestep = 0)}
	ind = which(data1$Timestep == timestep)
	TSnext = data1$Timestep[(last(ind)+1)]
	#ind2 = which(data2$Timestep == timestep)
	#TSnext2 = data2$Timestep[(last(ind2)+1)]
	#TS2 = rbind(filter(data1, Timestep == TSnext),filter(data2, Timestep == TSnext))
	TS2 = rbind(filter(data1, Timestep == TSnext),filter(data2, Timestep == TSnext),filter(data3, Timestep == TSnext))
	#TS1 = rbind(filter(data1, Timestep == (timestep)),filter(data2, Timestep == (timestep)))
	TS1 = rbind(filter(data1, Timestep == (timestep)),filter(data2, Timestep == (timestep)),filter(data3, Timestep == (timestep)))
	grid_dif = cbind(TS1,DifValue=(TS2$value - TS1$value))
	grid_dif$dTheta = ifelse(grid_dif$DifValue >= valmax, "increase > ThetaMax", ifelse(grid_dif$DifValue <= valmin, "decrease < ThetaMin", "umbrella"))
return(grid_dif)
}


#' @title Calculate delta values of 2 grids at different time steps
#' 
#' @description 3D
#'
#' @param data1,data2,data3 input data. has to be in the format XX.
#' @param timestep vector of timesteps to be considered for calculations.
#' @param valmax,valmin values where to cut off resulting differces, leaving a data-area of special consideration
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

grid3d_difs = function(data1,data2,data3=NA,timestep, valmin, valmax){
#active for all 3 compartments
#grid_difs = function(data1,data2,data3,timestep){
	if(is.na(data3)){data3 = data.frame(Timestep = 0)}
	ind = which(data1$Timestep == timestep)
	TSnext = data1$Timestep[(last(ind)+1)]
	#TS2 = rbind(filter(data1, Timestep == TSnext),filter(data2, Timestep == TSnext))
	TS2 = rbind(filter(data1, Timestep == TSnext),filter(data2, Timestep == TSnext),filter(data3, Timestep == TSnext))
	#TS1 = rbind(filter(data1, Timestep == (timestep)),filter(data2, Timestep == (timestep)))
	TS1 = rbind(filter(data1, Timestep == (timestep)),filter(data2, Timestep == (timestep)),filter(data3, Timestep == (timestep)))
	#grid_dif = cbind(TS1[,-6],DifValue=(TS2$value - TS1$value))
	grid_dif = cbind(TS1,DifValue=(TS2$value - TS1$value))
	grid_dif$dTheta = ifelse(grid_dif$DifValue >= valmax, "increase > ThetaMax", ifelse(grid_dif$DifValue <= valmin, "decrease < ThetaMin", "umbrella"))
return(grid_dif)
}

#' @title Calculate delta values of 2 grids at different time steps
#' 
#' @description 3D
#'
#' @param data1,data2,data3 input data. has to be in the format XX.
#' @param timestep vector of timesteps to be considered for calculations.
#' @param valmax,valmin values where to cut off resulting differces, leaving a data-area of special consideration
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

#grid3d_difs_profile = function(data1,data2,timestep,yloc, valmin, valmax){
#active for all 3 compartments
grid3d_difs_profile = function(data1,data2,data3=NA,timestep,yloc, valmin, valmax){
	if(is.na(data3)){data3 = data.frame(Timestep = 0, y=0)}
	ind = which(data1$Timestep == timestep)
	TSnext = data1$Timestep[(last(ind)+1)]
	TS2 = rbind(filter(filter(data1, Timestep == TSnext),y==yloc),filter(filter(data2, Timestep == TSnext),y==yloc),filter(filter(data3, Timestep == TSnext),y==yloc))
	TS1 = rbind(filter(filter(data1, Timestep == (timestep)),y==yloc),filter(filter(data2, Timestep == (timestep)),y==yloc),filter(filter(data3, Timestep == (timestep)),y==yloc))
	#grid_dif = cbind(TS1[,-6],DifValue=(TS2$value - TS1$value))
	grid_dif = cbind(TS1,DifValue=(TS2$value - TS1$value))
	grid_dif$dTheta = ifelse(grid_dif$DifValue >= valmax, "increase > ThetaMax", ifelse(grid_dif$DifValue <= valmin, "decrease < ThetaMin", "umbrella"))
return(grid_dif)
}


#' @title Calculate delta values of 2 grids at different time steps
#' 
#' @description 3D: SM differences each timestep
#'
#' @param data1,data2,data3 input data. has to be in the format XX.
#' @param timestep vector of timesteps to be considered for calculations.
#' @param valmax,valmin values where to cut off resulting differces, leaving a data-area of special consideration
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

grid3d_difs_allTS = function(data1,data2=NA,data3=NA){
#active for all 3 compartments
	if(is.na(data2)){data2 = data.frame(Timestep = 0)}
	if(is.na(data3)){data3 = data.frame(Timestep = 0)}
	for(i in unique(data1$Timestep)){
		ind = which(data1$Timestep == i)
		TSnext = data1$Timestep[(last(ind)+1)]
		#TS2 = rbind(filter(data1, Timestep == TSnext),filter(data2, Timestep == TSnext))
		TS2 = rbind(filter(data1, Timestep == TSnext),filter(data2, Timestep == TSnext),filter(data3, Timestep == TSnext))
		#TS1 = rbind(filter(data1, Timestep == (timestep)),filter(data2, Timestep == (timestep)))
		TS1 = rbind(filter(data1, Timestep == i),filter(data2, Timestep == i),filter(data3, Timestep == i))
		#grid_dif = cbind(TS1[,-6],DifValue=(TS2$value - TS1$value))
		grid_dif_TS = cbind(TS2,DifValue=(TS2$value - TS1$value))
		if(i == unique(data1$Timestep)[1]){
			grid_dif = grid_dif_TS
			next
		}
		grid_dif = rbind(grid_dif, grid_dif_TS)
	}
return(grid_dif)
}

#' @title make a sequence (movie) from model results
#' 
#' @description or differences between timesteps or absolute values each timestep
#'
#' @param data1,data2,data3 input data. has to be in the format XX.
#' @param timestep vector of timesteps to be considered for calculations.
#' @param valmax,valmin values where to cut off resulting differces, leaving a data-area of special consideration
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

ModOut_animate = function(input1,input2,counter,type,outputname,valmax,valmin){
#later also add 3rd input!! (gneiss)
#ModOut_animate = function(input1,input2,input3,counter,type,outputname){
msg = "Starting to create Images.."
print(msg)
saveGIF({
	switch(type,
		difference = {
		for(i in counter){
		if(i == max(counter)) break
		print(
		ggplot(data=grid3d_difs(input1,input2,i,valmin,valmax), aes(x=x, y=Depth)) + geom_tile(aes(fill=dTheta, width=.5, height=.5)) + ylim(5,0) + xlab("Horizontal [m]") + ylab("Depth [m]") + ggtitle(paste("Timestep: ",i,sep="")) +
		scale_fill_manual(values=c("increase > ThetaMax"="#151B54","decrease < ThetaMin"="#8C001A","umbrella"="white"))
		) #end print
		}},
		real = {
		for(i in counter){
		if(i == max(counter)) break
		print(
		ggplot(data=rbind(filter(input1, Timestep == i),filter(input2, Timestep == i)), aes(x=x, y=Depth)) + geom_tile(aes(fill=value, width=.5, height=.5)) + ylim(5,0) + xlab("Horizontal [m]") + ylab("Depth [m]")+ ggtitle(paste("Timestep: ",i,sep="")) + 
		scale_fill_gradientn(colours=c("#8C001A","#437C17","#151B54"), breaks=c(0.01,0.1,0.2,0.25,0.3,0.4,0.5)) 
		) #end print
		}}
	)
	},
	movie.name=paste("./",outputname,".gif",sep=""), clean=T, interval=0.05, ani.width=600, ani.height=600)
	#movie.name="SMgrid_test.gif", clean=T, interval=0.05, ani.width=600, ani.height=600)

#return(print(paste("finished: ",outputname,".gif"," generated in working directory",sep="")))
fin="finished!"
return(fin)
}


#' @title plot 3D grid
#' 
#' @description or differences between timesteps or absolute values each timestep
#'
#' @param data1,data2,data3 input data. has to be in the format XX.
#' @param timestep vector of timesteps to be considered for calculations.
#' @param valmax,valmin values where to cut off resulting differces, leaving a data-area of special consideration
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

grid3d_plot = function(input,type){
#later also add 3rd input!! (gneiss)
switch(type,
		difference = {
		grid_plot = ggplot(data=input, aes(x=x, y=Depth)) + geom_tile(aes(fill=dTheta, width=.5, height=.5)) + ylim(5,0) + xlab("Horizontal [m]") + ylab("Depth [m]") + ggtitle(paste("Timestep: ",input$Timestep," - Y= ",unique(input$y),sep="")) +
		scale_fill_manual(values=c("increase > ThetaMax"="#151B54","decrease < ThetaMin"="#8C001A","umbrella"="white"))
		},
		real = {
		grid_plot = ggplot(data=input, aes(x=x, y=Depth)) + geom_tile(aes(fill=value, width=.5, height=.5)) + ylim(5,0) + xlab("Horizontal [m]") + ylab("Depth [m]")+ ggtitle(paste("Timestep: ",input$Timestep," - Y= ",unique(input$y),sep="")) + 
		scale_fill_gradientn(colours=c("#8C001A","#437C17","#151B54"), breaks=c(0.01,0.1,0.2,0.25,0.3,0.4,0.5)) 
		},
		real_mean = {
		grid_plot = ggplot(data=input, aes(x=x, y=Depth)) + geom_tile(aes(fill=val_mean, width=.5, height=.5)) + ylim(5,0) + xlab("Horizontal [m]") + ylab("Depth [m]")+ ggtitle(paste("Timestep: mean"," - Y= ",unique(input$y),sep="")) + 
		scale_fill_gradientn(colours=c("#8C001A","#437C17","#151B54"), breaks=c(0.01,0.1,0.2,0.25,0.3,0.4,0.5)) 
		},
		SD = {
		grid_plot = ggplot(data=input, aes(x=x, y=Depth)) + geom_tile(aes(fill=value, width=.5, height=.5)) + ylim(5,0) + xlab("Horizontal [m]") + ylab("Depth [m]")+ ggtitle(paste("Timestep: sd"," - Y= ",unique(input$y),sep="")) + 
		scale_fill_gradientn(colours=c("#8C001A","#437C17","#151B54"))#, breaks=c(0.01,0.1,0.2,0.25,0.3,0.4,0.5)) 
		}
	)
return(grid_plot)
}


#' @title estimate limiting points of umbrella space / curve for FIXED y-profiles
#' 
#' @description scanning area from beneath, looking for depth with smallest difference
#'
#' @param data1,data2,data3 input data. has to be in the format XX.
#' @param timestep vector of timesteps to be considered for calculations.
#' @param valmax,valmin values where to cut off resulting differces, leaving a data-area of special consideration
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

umbrella_points_fixedY = function(data1,data2,data3=NA,timestep,yloc,valmin,valmax){
	#if(is.na(data3)){data3 = data.frame(Timestep = 0, y=0)}
	umbrellapoints = data.frame(x=NA,Depth=NA,y=NA)
	#setup grid
	SMgrid=grid3d_difs(data1,data2,data3,timestep,valmin,valmax)
	#choose the right counter
	#all x values wont work because layer bvcv does not have the same discretization
	#therefore in every other column, depth are only chosen from layer ah
	#this limits to a depth of ah = 0.1; and is not what we want
	#counter = unique(SMgrid_test$x)
	#instead we chose x discretization of layer bvcv:
	counter = unique(SMgrid$x)[c(TRUE, FALSE)]
	#set up counter for y-profiles
	ycount = unique(SMgrid$y)[c(TRUE, FALSE)]
	# get umbrellapoints for FIXED Y
	j=1
	ycol =2 
	for(i in counter){
		depth.filter = 	filter(SMgrid, x == i) %>%
				filter(y == ycol) %>%
				filter(DifValue > valmin & DifValue < valmax)
		umbrellapoints[j,2] = max(depth.filter$Depth)
		range(depth.filter$Depth)
		umbrellapoints[j,1] = i
		umbrellapoints[j,3] = ycol
		j=j+1
	}
	
return(umbrellapoints)
}

#' @title estimate limiting points of umbrella space / curve for ALL y-profiles
#' 
#' @description scanning area from beneath, looking for depth with smallest difference
#'
#' @param data1,data2,data3 input data. has to be in the format XX.
#' @param timestep vector of timesteps to be considered for calculations.
#' @param valmax,valmin values where to cut off resulting differces, leaving a data-area of special consideration
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

umbrella_points = function(data1,data2,data3=NA,timestep,valmin,valmax){
	#if(is.na(data3)){data3 = data.frame(Timestep = 0, y=0)}
	umbrellapoints = data.frame(x=NA,Depth=NA,y=NA)
	#setup grid
	SMgrid=grid3d_difs(data1,data2,data3,timestep,valmin,valmax)
	counter = unique(SMgrid$x)[c(TRUE, FALSE)]
	ycount = unique(SMgrid$y)[c(TRUE, FALSE)]
	j=1
	for (k in ycount){
	SMgrid=grid3d_difs_profile(SMgrid_ah,SMgrid_bvcv,timestep=time_i,yloc=k,valmin=valmin,valmax=valmax)
	for(i in counter){
		depth.filter = 	filter(SMgrid, x == i) %>%
				filter(DifValue > valmin & DifValue < valmax)
		umbrellapoints[j,2] = max(depth.filter$Depth)
		range(depth.filter$Depth)
		umbrellapoints[j,1] = i
		umbrellapoints[j,3] = k
		j=j+1
	}
	#plot(umbrellapoints, ylim=c(5,0))
	}
	
return(umbrellapoints)
}


#' @title estimate limiting points of umbrella space / curve for ALL y-profiles with MEAN soil moisture values
#' 
#' @description scanning area from beneath, looking for depth with smallest difference
#'
#' @param data1,data2,data3 input data. has to be in the format XX.
#' @param timestep vector of timesteps to be considered for calculations.
#' @param valmax,valmin values where to cut off resulting differces, leaving a data-area of special consideration
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

umbrella_points_mean = function(data1,data2,data3=NA,valmin,valmax){
	#if(is.na(data3)){data3 = data.frame(Timestep = 0, y=0)}
	umbrellapoints = data.frame(x=NA,Depth=NA,y=NA)
	#setup grid
	#SMgrid=grid3d_difs(data1,data2,data3,timestep,valmin,valmax)
	SMgrid_difs = rbind(data1,data2)
	counter = unique(data1$x)[c(TRUE, FALSE)]
	ycount = unique(data2$y)[c(TRUE, FALSE)]
	j=1
	for (k in ycount){
	#SMgrid=grid3d_difs_profile(SMgrid_ah,SMgrid_bvcv,timestep=time_i,yloc=k,valmin=limlow,valmax=limup)
	for(i in counter){
		depth.filter = 	filter(SMgrid_difs, x == i) %>%
				filter(y == k) %>%
				filter(DifValue > valmin & DifValue < valmax)
		umbrellapoints[j,2] = max(depth.filter$Depth)
		umbrellapoints[j,1] = i
		umbrellapoints[j,3] = k
		j=j+1
	}
	#plot(umbrellapoints, ylim=c(5,0))
	}
	
return(umbrellapoints)
}

#' @title estimate limiting points of umbrella space / curve for ALL y-profiles with SD of soil moisture point time series
#' 
#' @description scanning area from beneath, looking for depth with smallest difference
#'
#' @param data1,data2,data3 input data. has to be in the format XX.
#' @param timestep vector of timesteps to be considered for calculations.
#' @param valmax,valmin values where to cut off resulting differces, leaving a data-area of special consideration
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

umbrella_points_sd = function(data1,data2=NA,data3=NA,valmax,counters){
	#if(is.na(data3)){data3 = data.frame(Timestep = 0, y=0)}
	umbrellapoints = data.frame(x=NA,Depth=NA,y=NA)
	#setup grid
	#SMgrid=grid3d_difs(data1,data2,data3,timestep,valmin,valmax)
	#SMgrid_sd = rbind(data1,data2,data3)
	SMgrid_sd = data1
	#counter = unique(data1$x)[c(TRUE, FALSE)]
	#ycount = unique(data2$y)[c(TRUE, FALSE)]
	j=1
	for (k in counters$y){
	#SMgrid=grid3d_difs_profile(SMgrid_ah,SMgrid_bvcv,timestep=time_i,yloc=k,valmin=limlow,valmax=limup)
	for(i in counters$x){
		depth.filter = 	filter(SMgrid_sd, x == i) %>%
				filter(y == k) %>%
				filter(value < valmax)
		if(length(depth.filter$Depth) < 1) umbrellapoints[j,2] = NA
		else umbrellapoints[j,2] = max(depth.filter$Depth)
		umbrellapoints[j,1] = i
		umbrellapoints[j,3] = k
		j=j+1
	}
	#plot(umbrellapoints, ylim=c(5,0))
	}
	
return(umbrellapoints)
}


#' @title estimate limiting points of umbrella space / curve for ALL y-profiles with SD of soil moisture point time series
#' 
#' @description scanning area from beneath, looking for depth with smallest difference
#'
#' @param data1,data2,data3 input data. has to be in the format XX.
#' @param timestep vector of timesteps to be considered for calculations.
#' @param valmax,valmin values where to cut off resulting differces, leaving a data-area of special consideration
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

umbrella_points_cluster = function(data1,clusternum, counters){
	#if(is.na(data3)){data3 = data.frame(Timestep = 0, y=0)}
	umbrellapoints = data.frame(x=NA,Depth=NA,y=NA)
	#setup grid
	#SMgrid=grid3d_difs(data1,data2,data3,timestep,valmin,valmax)
	SMgrid = data1
	#counter = unique(data1$x)[c(TRUE, FALSE)]
	#ycount = unique(data1$y)[c(TRUE, FALSE)]
	j=1
	for (k in counters$y){
	#SMgrid=grid3d_difs_profile(SMgrid_ah,SMgrid_bvcv,timestep=time_i,yloc=k,valmin=limlow,valmax=limup)
	for(i in counters$x){
		depth.filter = 	filter(SMgrid, x == i) %>%
				filter(y == k) %>%
				filter(cluster == clusternum)
		umbrellapoints[j,2] = max(depth.filter$Depth)
		umbrellapoints[j,1] = i
		umbrellapoints[j,3] = k
		j=j+1
	}
	#plot(umbrellapoints, ylim=c(5,0))
	}
	
return(umbrellapoints)
}

#' @title estimate limiting points of umbrella space / curve for ALL y-profiles with SD of soil moisture point time series
#' 
#' @description scanning area from beneath, looking for depth with smallest difference
#'
#' @param data1,data2,data3 input data. has to be in the format XX.
#' @param timestep vector of timesteps to be considered for calculations.
#' @param valmax,valmin values where to cut off resulting differces, leaving a data-area of special consideration
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

umbrella_points2D_cluster = function(data1,clusternum,counter){
	#if(is.na(data3)){data3 = data.frame(Timestep = 0, y=0)}
	umbrellapoints = data.frame(x=NA,Depth=NA)
	#setup grid
	#SMgrid=grid3d_difs(data1,data2,data3,timestep,valmin,valmax)
	SMgrid = data1
	#counter = unique(data1$x)[c(TRUE, FALSE)]
	#ycount = unique(data1$y)[c(TRUE, FALSE)]
	j=1
	#for (k in ycount){
	#SMgrid=grid3d_difs_profile(SMgrid_ah,SMgrid_bvcv,timestep=time_i,yloc=k,valmin=limlow,valmax=limup)
	for(i in counter){
		depth.filter = 	filter(SMgrid, x == i) %>%
				#filter(y == k) %>%
				filter(cluster == clusternum)
		if(length(depth.filter$Depth) < 1) umbrellapoints[j,2] = NA
		else umbrellapoints[j,2] = max(depth.filter$Depth)
		umbrellapoints[j,1] = i
		#umbrellapoints[j,3] = k
		j=j+1
	}
	#plot(umbrellapoints, ylim=c(5,0))
	#}
	
return(umbrellapoints)
}

#' @title estimate function for umbrella space
#' 
#' @description using lm()
#'
#' @param data1,data2,data3 input data. has to be in the format XX.
#' @param timestep vector of timesteps to be considered for calculations.
#' @param valmax,valmin values where to cut off resulting differces, leaving a data-area of special consideration
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

predict_umbrellaspace = function(inputdata, displ_x, displ_y, polyorder, fitgeometry,plotting=T){

#displ_x = 9.5
#displ_y = 10
x = inputdata[,1] - displ_x
y = inputdata[,2] - displ_y
z = inputdata[,3]
#generate model
switch(fitgeometry,
       paraboloid = {model = lm(z ~ poly(x,polyorder,raw=T) + poly(y,polyorder,raw=T))}
       #pyramid = {model = lm(z ~ poly(x,polyorder,raw=T) + poly(y,polyorder,raw=T))},
       #cone = {model = lm(z ~ poly(x,polyorder,raw=T) + poly(y,polyorder,raw=T))}
       )
#generate output grid
x.pred = unique(x)
y.pred = unique(y)
xy = expand.grid(x = x.pred, y = y.pred)
#use model to predict depth
z.vals = predict(model, newdata = xy)
z.pred = matrix(-1*z.vals,
		nrow = length(x.pred), ncol = length(y.pred))
#points !?
#fitpoints = predict(model_parapoloid)
#predicted parameters
#parameter_estimated = coef(model)
parameter_estimated = summary(model)
if(plotting){
	library(plot3D)
	library(plot3Drgl)
	#static 3d plot
	scatter3D(x,y,-1*z,theta=20,phi=20,ticktype="detailed",surf = list(x = x.pred, y = y.pred, z = z.pred,facets=NA),cex.axis = 2, cex.lab = 2.5,bty = "b2", colkey = FALSE)
	#interactive 3d plot
	#plotrgl()
}
#save plot
#png(file="/home/mreich/Dokumente/Wettzell/hydrologicalmodelling/hydrus3d_out/Bodenart_SandyLoam_Example/umbrellacurve3D_TS26200_paraboloid.png", width=1000,height=1000,res=150)
#scatter3D(x,y,-1*z,theta=20,phi=20,ticktype="detailed",surf = list(x = x.pred, y = y.pred, z = z.pred,facets=NA),bty="g")
#dev.off()
#UmbrellaSpacePlot <<- UmbrellaSpacePlot
return(parameter_estimated)
}

#' @title estimate function for umbrella space in 2D
#' 
#' @description using lm()
#'
#' @param data1,data2,data3 input data. has to be in the format XX.
#' @param timestep vector of timesteps to be considered for calculations.
#' @param valmax,valmin values where to cut off resulting differces, leaving a data-area of special consideration
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

predict_umbrellaspace_2d = function(inputdata, displ_x, polyorder, fitgeometry,plotting=T){

#displ_x = 9.5
#displ_y = 10
x = inputdata[,1] - displ_x
z = inputdata[,2]
#generate model
switch(fitgeometry,
       parabel = {model = lm(z ~ poly(x,polyorder,raw=T))}
       #pyramid = {model = lm(z ~ poly(x,polyorder,raw=T) + poly(y,polyorder,raw=T))},
       #cone = {model = lm(z ~ poly(x,polyorder,raw=T) + poly(y,polyorder,raw=T))}
       )
#generate output grid
x.pred = unique(x)
#xy = expand.grid(x = x.pred, y = y.pred)
#use model to predict depth
z.vals = predict(model, newdata = x.pred)
z.pred = matrix(-1*z.vals,
		nrow = length(x.pred), ncol = 1)
#points !?
#fitpoints = predict(model_parapoloid)
#predicted parameters
#parameter_estimated = coef(model)
parameter_estimated = summary(model)
if(plotting){
	library(plot3D)
	library(plot3Drgl)
	#static 3d plot
	plot(x,z)
	points(x.pred,z.pred, col="blue")
	#scatter3D(x,y,-1*z,theta=20,phi=20,ticktype="detailed",surf = list(x = x.pred, z = z.pred,facets=NA),cex.axis = 2, cex.lab = 2.5,bty = "b2", colkey = FALSE)
	#interactive 3d plot
	#plotrgl()
}
#save plot
#png(file="/home/mreich/Dokumente/Wettzell/hydrologicalmodelling/hydrus3d_out/Bodenart_SandyLoam_Example/umbrellacurve3D_TS26200_paraboloid.png", width=1000,height=1000,res=150)
#scatter3D(x,y,-1*z,theta=20,phi=20,ticktype="detailed",surf = list(x = x.pred, y = y.pred, z = z.pred,facets=NA),bty="g")
#dev.off()
#UmbrellaSpacePlot <<- UmbrellaSpacePlot
return(parameter_estimated)
}

#' @title aggregate data from hourly to daily (weekly, monthly) averages
#' 
#' @description using sorting and a moving window average
#'
#' @param data1,data2,data3 input data. has to be in the format XX.
#' @param timestep vector of timesteps to be considered for calculations.
#' @param valmax,valmin values where to cut off resulting differces, leaving a data-area of special consideration
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

grid_aggregate = function(input, win_width, every_n){
	data_sort = arrange(input, x,y,Depth)
	grid.agg = data.frame(Timestep = data_sort$Timestep[seq(1,length(data_sort[,1]),24)],
				Depth = data_sort$Depth[seq(1,length(data_sort[,1]),24)],
				x = data_sort$x[seq(1,length(data_sort[,1]),24)],
				y = data_sort$y[seq(1,length(data_sort[,1]),24)],
				value=rollapply(data_sort$value, width=win_width,by=every_n,FUN=mean,align="left"))
return(grid.agg)
}

#' @title compute grouped standard deviation, range, max from values averaged in y-direction (of 3D grid)
#' 
#' @description using sorting and a moving window average
#'
#' @param data1,data2,data3 input data. has to be in the format XX.
#' @param timestep vector of timesteps to be considered for calculations.
#' @param valmax,valmin values where to cut off resulting differces, leaving a data-area of special consideration
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

grid_stats = function(input, method, dimension){
	# chose method for statistics
	switch(method,
	       MEAN = {
		grid.stats = group_by(input, x,y,Depth) %>%
			summarize(value=mean(value, na.rm=T))
	       },
	       SD = {
		grid.stats = group_by(input, x,y,Depth) %>%
			summarize(value=sd(value, na.rm=T))
	       },
	       RANGE = {
		grid.stats = group_by(input, x,y,Depth) %>%
			summarize(value=abs(range(value)[2]-range(value)[1]))
	       },
	       MAX = {
		grid.stats = group_by(input, x,y,Depth) %>%
			summarize(value=max(value, na.rm=T))
	       }
	       )

	# if dim = 2, all values will be agregated (mean) over y
	if(dimension==2){
	grid.stats = group_by(grid.stats, x,Depth) %>%
	       		summarize(value=mean(value))
	}
return(grid.stats)
}


#' @title compute grouped standard deviation, range, max, mean from values of 2D grid
#' 
#' @description using sorting and a moving window average
#'
#' @param data1,data2,data3 input data. has to be in the format XX.
#' @param timestep vector of timesteps to be considered for calculations.
#' @param valmax,valmin values where to cut off resulting differces, leaving a data-area of special consideration
#' @details missing
#' @references Marvin Reich (2016), mreich@@gfz-potsdam.de
#' @examples missing

grid_stats_2d = function(input, method){
	# chose method for statistics
	switch(method,
	       MEAN = {
		grid.stats = group_by(input, x,Depth) %>%
			summarize(value=mean(value, na.rm=T))
	       },
	       SD = {
		grid.stats = group_by(input, x,Depth) %>%
			summarize(value=sd(value, na.rm=T))
	       },
	       RANGE = {
		grid.stats = group_by(input, x,Depth) %>%
			summarize(value=abs(range(value, na.rm=T)[2]-range(value, na.rm=T)[1]))
	       },
	       MAX = {
		grid.stats = group_by(input, x,Depth) %>%
			summarize(value=max(value, na.rm=T))
	       }
	       )
return(grid.stats)
}

#' @title Plot timeseries input / output : Hydrus 1D
#'
#' @description test
#'
#' @param folder foldername of project to read
#' @param vertical_nodes mumber of nodes used in the vertical model discretisation
#' @param timesteps number of modeled timesteps
#' @param t_printout number of model printout times
#' @param plotting indicate if standard parameters should be plotted (TRUE); default is FALSE
#' @param datetime vector of timestamps of the modeled timeseries. if not provided output time will be in counts of timesteps
#' 
#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples missing

plot_hydrus1d_InOutput <- function(input_ts, tlevel_out, obsNode, timesteps, vertical_nodes, t_printout, tprint_show=10){
  
 datetime = timeseries$Datetime	

 tlevel.plot = melt(tlevel_out, id.vars=c("datetime","Time"), measure.vars=colnames(tlevel)[2:22], variable.name = "Parameters")
 #nodal.plot = melt(nodal_out, id.vars=c("Time", "Depth"), measure.vars=colnames(nodal_info)[-1:-3], variable.name = "Parameters")
 obsNode.melt = melt(zootodf(obsNode), id="time", variable.name="node", value.name="theta")
 
 #printtimes = unique(nodal.plot$Time) 
 #printdates = c(datetime[1],datetime[printtimes])
 #prints = round(seq(1,t_printout, length.out=tprint_show),0)
 #printtimes.gg = printtimes[prints]
 #printdates.gg = printdates[prints]
 #dates_times = data.frame(Time=printtimes, Dates=format(printdates, "%Y-%m-%d"))
 ##dates_times = data.frame(Time=printtimes, Dates=printdates)

 ##xcolors = colorRampPalette(c("blue","red", "orange","darkgreen"))(t_printout+1)
 #xcolors = colorRampPalette(c("blue","red", "orange","darkgreen"))(length(printdates.gg))
 ##  xcolors = colorRampPalette(c("blue","green"))(t_printout+1)

 ##  balance.plot = melt(balance, id.vars="Time", measure.vars=colnames(balance)[-1])
 ##subsetting to choose standard parameters
 #choose.params.tlevel = c("vTop", "vRoot", "vBot", "hTop", "hRoot", "hBot")
 #units = data.frame(type=c("Flux","Head"), Unit=c("[m/h]","[m]"))
 #cols = data.frame(Parameters = choose.params.tlevel, Boundary = c("Top","Root","Bottom"))
 #type = data.frame(Parameters = c("vTop","vBot","hTop","hBot"), type=c("Flux","Flux","Head","Head"))
 #tlevel.plot.filter = filter(tlevel.plot, grepl("vTop|vBot|hTop|hBot",Parameters)) %>%
                         #inner_join(type,by="Parameters") %>%
                         #inner_join(units,by="type") %>%
			#mutate(typeUnits = paste(type,Unit)) %>%
			#inner_join(cols, by="Parameters") %>%
			#select(datetime, Boundary, value, Unit, typeUnits)

			## tlevel.plot.filter$ParamUnits = factor(tlevel.plot.filter$ParamUnits,levels=c("Head [m]","Moisture [%]","K [m/h]","Flux [m/h]"))
 #tlevel.plot.filter$Boundary = factor(tlevel.plot.filter$Boundary,levels=c("Top","Root","Bottom"))

 ##tlevel.gg = ggplot(tlevel.plot.filter, aes(x=datetime, y=value)) + geom_line(aes(linetype=Boundary)) + 
 #tlevel.gg = ggplot(level_data, aes(x=datetime, y=value)) + geom_line(aes(colour=Boundary)) + 
          #xlab("") + ylab ("") + 
	 #facet_grid(typeUnits ~ ., scale="free_y") + 
	 #theme(axis.text.x=element_text(face="bold"),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) 

 ## tlevel data (hydrus)
 type = data.frame(Parameters = c("vTop","vBot","hTop","hBot"), Unit=c("[m/h]","[m/h]","[m]","[m]"), WaterType=c("Flux","Flux","Head","Head"))
 tlevel.plot.filter = filter(tlevel.plot, grepl("vTop|vBot|hTop|hBot",Parameters)) %>%
	 		select(-Time) %>%
 			inner_join(type,by="Parameters") %>%
			mutate(ParamUnits = paste(Parameters,Unit))
 ## hydrus input time series
 timeseries.melt = melt(timeseries, id.vars="Datetime")
 #units_tsinput = data.frame(variable=c("Precipitation","ET0","Groundwater"), Unit=c("[mm/h]","[mm/h]","[m above lower model boundary"))
 units_tsinput = data.frame(variable=c("Precipitation","ET0","Groundwater"), Unit=c("[mm/h]","[mm/h]","[m above LmodelB"), WaterType=c("Flux","Flux","Head"))
 timeseries.filter = inner_join(timeseries.melt, units_tsinput,by="variable") %>%
			mutate(ParamUnits = paste(variable,Unit))
 colnames(timeseries.filter)[1:2] = c("datetime","Parameters")

 ## observation nodes (hydrus)
 colnames(obsNode.melt) = c("datetime","node","value")
 obsNode.filter = mutate(obsNode.melt, Unit = "[VWC %]") %>%
	 		mutate(WaterType = "State") %>%
			mutate(Parameters = paste("obs node",node)) %>%
			mutate(ParamUnits = paste("obs node",Unit)) %>%
			select(-node)

 # combine datasets
 level_data = rbind(tlevel.plot.filter, timeseries.filter, obsNode.filter)
 # further classification for better plotting
 #subtypes = data.frame(Parameters=c("Precipitation","ET0","Groundwater", "vTop","vBot","hTop","hBot"),
		       #subclass=c("Atmo","Atmo2","GW","Atmoflux","GWflux","Atmo2","GW"))
 #level_data = inner_join(level_data, subtypes, by="Parameters")

 data.gg = ggplot(level_data, aes(x=datetime, y=value)) + geom_line(aes(colour=Parameters)) + 
 	 xlab("") + ylab ("") + 
	 facet_grid(ParamUnits ~ ., scale="free_y") + 
	 #facet_grid(WaterType ~ ., scale="free_y") + 
	 #facet_grid(subclass ~ ., scale="free_y") + 
	 theme(axis.text.x=element_text(face="bold"),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) 
 
 #choose.params.nodal = c("Head", "Moisture", "K", "Flux")
 #units = data.frame(Parameters=choose.params.nodal, Unit=c("[m]", "[%]", "[m/h]", "[m/h]"))
 ##  nodal.plot.filter = 	filter(nodal.plot, Parameters == choose.params.nodal) %>%
 #nodal.plot.filter = 	filter(nodal.plot, grepl("Head|Moisture|K|Flux",Parameters)) %>%
			#transform(Hours=as.factor(Time)) %>%
			#transform(Days=as.factor(Time/24)) %>%
                         #inner_join(units,by="Parameters") %>%
			#mutate(ParamUnits = paste(Parameters,Unit)) %>%
                         #inner_join(dates_times,by="Time") 
			
 #nodal.plot.filter$ParamUnits = factor(nodal.plot.filter$ParamUnits,levels=c("Head [m]","Moisture [%]","K [m/h]","Flux [m/h]"))

 ##  nodal.gg = ggplot(nodal.plot.filter, aes(x=value, y=Depth, colour=Dates, group=Dates)) + geom_path() + xlab("") + ylab ("Depth [m]") + facet_grid(.~ParamUnits, scale="free_x") + theme(legend.position = "bottom",legend.text=element_text(size=10),legend.title=element_blank())
 #nodal.gg = ggplot(filter(nodal.plot.filter, Time == printtimes.gg), aes(x=value, y=Depth, colour=Dates, group=Dates)) + geom_path() + xlab("") + ylab ("Depth [m]") + 
		 #facet_grid(.~ParamUnits, scale="free_x") + theme(legend.position = "none") + 
		#scale_color_manual(values = xcolors)
 ##merge both plots together
 #plot.gg = grid.arrange(tlevel.gg,nodal.gg)

 #return(plot.gg)
 return(data.gg)
}

