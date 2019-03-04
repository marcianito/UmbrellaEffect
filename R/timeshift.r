#' @title Time lag estimation via cross-correlation of two timeseries
#'
#' @description This function is based on standard time series decomposition and afterwards cross-correlation estimation between the two input datasets.
#' @description The main aim is to find dominant time shifts between two time series (e.g. soil moisture, climate data, precipitation stations, etc).
#'
#' @param data1,data2 Timeseries of type .zoo.
#' @param nmax Number of correlation values output.
#' @param norm logical. Use normalized data (default) or not.
#' @param stl logical. Use decomposed data for cross-correlation (default) or not.
#' @param plotting logical. Output cross-correlation plots. Default is F, own output plots are provided. Should be kept FALSE.
#' @param swin Time window used for seasonality estimation. Only used if stl=T.
#' @param twin Time window used for trend estimation. Only used if stl=T.
#' @param ... Further parameters passed to internal functions.

#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples outside.cor = ccf.zoo(besidesBuilding$mux43_04,besidesBuilding$mux43_08,5,T,T,T,24,17521)
#' 

timeshift = function(data1,data2,nmax,norm=T,stl=F,plotting=T,swin,twin,sensor1="data1", sensor2="data2",...){
# library(zoo);
  Sys.setenv(TZ = "GMT")
# library(xts)
# library(gridExtra)
# if(plotting==TRUE){
# library(reshape2)
# library(ggplot2)}
  #get both datasets to use the same timestamp-frequency
  data1_sameF = aggregate(data1, function(x) as.POSIXct(trunc(x, "hour"), tz = "GMT"), mean, na.rm=T)
  data2_sameF = aggregate(data2, function(x) as.POSIXct(trunc(x, "hour"), tz = "GMT"), mean, na.rm=T)
  #merging all data together can result in problems when na.approx doesn't find 2 non-NA values for interpolation
  #therefor merge has to be used with option all=F
  #data.merge = merge(data1_sameF, data2_sameF, all=T, fill= NA)
  data.merge = merge(data1_sameF, data2_sameF, all=F)
  #merge the datasets, depending on the analysis time
  #this is necessary to avoid na.approx-problems with too many NA's
  #and on the other hand don't loose any data if stl=F
  #if(stl==T) data.merge=merge.zoo(na.approx(data1),na.approx(data2), all=T, fill=NA) 
  #else data.merge=merge.zoo(data1,data2, all=F) 
  #aggregate to hourly time
  #   data1.agg = period.apply(data1, endpoints(data1,"hours"), mean, na.rm=T)
  #   data2.agg = period.apply(data2, endpoints(data2,"hours"), mean, na.rm=T)
  #split in yearly intervalls
  ts.years = endpoints(data.merge, "years") 
  data1.year=list()
  data2.year=list()
  for(i in 1:(length(ts.years)-1)){
	data1.year[[i]]=data.merge[ts.years[i]:ts.years[i+1],1]
	data2.year[[i]]=data.merge[ts.years[i]:ts.years[i+1],2]
  }
  if(norm==T){ #use normalize data as input
  data1.year.norm = lapply(data1.year, function(x) (x - mean(x, na.rm=T))/sd(x,na.rm=T))
  data2.year.norm = lapply(data2.year, function(x) (x - mean(x, na.rm=T))/sd(x,na.rm=T))
  }
  else{ #use original data as input
  data1.year.norm = data1.year
  data2.year.norm = data2.year
  }
  #cross-correlation for each year
  correlation=data.frame(matrix(ncol=(length(ts.years)-1)*2,nrow=nmax)) 
  res.plots = list(); data.plots= list(); data1.stl = list(); data2.stl = list()
  j=1;k=2; cornames=NA
  for(i in 1:(length(correlation)/2)){
  	p1.flag=1;p2.flag=1 #setting standard plotting parameter
    	timedif.data = as.numeric(difftime(index(data1.year.norm[[i]][2]),index(data1.year.norm[[i]][1]), units="hours"))
	freq = 24/timedif.data #frequency in hours
	year=unique(format(index(data1.year.norm[[i]]), "%Y")) #get year from current TS
	if(stl==TRUE){ #decompose TS
	#print(i) #only for debugging
	check1=try(stl(ts(data1.year.norm[[i]],frequency=freq),s.window=swin,t.window=twin),TRUE)
	if(class(check1)=="try-error"){
		data1.na= na.approx(data1.year.norm[[i]])
		#if(i==2){browser()}
		check1.2=try(stl(ts(data1.na,frequency=freq),s.window=swin,t.window=twin),TRUE)
		if(class(check1.2)=="try-error"){
		correlation[,j] = NA; cornames[j]= paste("cor_",max(year),sep="") ;j=j+2
		correlation[,k] = NA; cornames[k]=paste("lagHOUR_",max(year),sep="") ;k=k+2
		print(paste(deparse(substitute(data1)),"Year",max(year), ": calculation not possible. Too much NA data.", sep=" "))
	       	next #skip year if stl is not possible because of missing data
		}
		data1.stl[[i]]= stl(ts(data1.na,frequency=freq),s.window=swin,t.window=twin)
		data1.in= data1.stl[[i]]$time.series[,3] #residuals of stl()
		print(paste(deparse(substitute(data1)),"Year",max(year), ": some NA data was approximized.", sep=" "))
		p1.flag=NA #set flag for choosing what data to plot
	}
	else{
	data1.stl[[i]]= stl(ts(data1.year.norm[[i]],frequency=freq),s.window=swin,t.window=twin)
	data1.in= data1.stl[[i]]$time.series[,3] #residuals of stl()
	} #end check-run data1
	check2=try(stl(ts(data2.year.norm[[i]],frequency=freq),s.window=swin,t.window=twin),TRUE)
	if(class(check2)=="try-error"){ #skip year if stl is not possible because of missing data
		data2.na= na.approx(data2.year.norm[[i]])
		check2.2=try(stl(ts(data2.na,frequency=freq),s.window=swin,t.window=twin),TRUE)
		if(class(check2.2)=="try-error"){
		correlation[,j] = NA; cornames[j]= paste("cor_",max(year),sep="") ;j=j+2
		correlation[,k] = NA; cornames[k]=paste("lagHOUR_",max(year),sep="") ;k=k+2
		print(paste(deparse(substitute(data2)),"Year",max(year), ": calculation not possible. Too much NA data.", sep=" "))
	       	next #skip year if stl is not possible because of missing data
		}
		data2.stl[[i]]= stl(ts(data2.na,frequency=freq),s.window=swin,t.window=twin)
		data2.in= data2.stl[[i]]$time.series[,3] #residuals of stl()
		print(paste(deparse(substitute(data2)),"Year",max(year), ": some NA data was approximized.", sep=" "))
		p2.flag=NA #set flag for choosing what data to plot
	}
	else{
	data2.stl[[i]]= stl(ts(data2.year.norm[[i]],frequency=freq),s.window=swin,t.window=twin)
	data2.in= data2.stl[[i]]$time.series[,3] #residuals of stl()
	} #end check-run data2
	}
	else{ #use NON-decomposed TS
	#detrend both time series, using linear fit
	#method uses least-squares fit by an regression analysis
	#datasets needs to be approximized
	#if not, lm() returns shorter datasets (if NAs are within the data)
	#this would lead to wrong detrending when resting data(org)-trend!
	data1.approx = na.approx(coredata(data1.year.norm[[i]]))
	data1.model = lm(data1.approx ~ I(1:length(data1.approx)))
	data1.trend = predict(data1.model)
	#summary(model)
	#removing the trend from "raw"-data
	data1.in = data1.approx - data1.trend
	#plot(coredata(data1.year.norm[[i]]))
	#abline(data1.model)
	#lines(data1.detrended, col="blue")
	data2.approx = na.approx(coredata(data2.year.norm[[i]]))
	data2.model = lm(data2.approx ~ I(1:length(data2.approx)))
	data2.trend = predict(data2.model)
	#summary(model)
	#removing the trend from "raw"-data
	data2.in = data2.approx - data2.trend
	#plot(coredata(data2.year.norm[[i]]))
	#abline(data2.model)
	#lines(data2.detrended, col="blue")
	#old version without detrending if no "stl" is used
	#data1.in = coredata(data1.year.norm[[i]])
	#data2.in = coredata(data2.year.norm[[i]])
	}
	# cross-correlation
	#browser()
	res = ccf(data1.in,data2.in,na.action=na.pass,plot=F,main=max(year),...)
	# browser()
    	acf.maxs=NA; cornmax=NA
    	for(ii in 1:nmax){ #finding nth max values
		 if(is.na(sort(res$acf, TRUE)[ii])==T){ cornmax[ii]=NA; acf.maxs[ii]=NA; next}
		 cornmax[ii] = sort(res$acf, TRUE)[ii]
     		 acf.maxs[ii] = which(res$acf == cornmax[ii])
   	}
    	lags.dom = res$lag[acf.maxs]
	#convert output from frequency to time units[hours]
	if(stl==TRUE){lags.real = lags.dom*freq*timedif.data}
	#if(stl==TRUE){lags.real = lags.dom/(timedif.data)} #old
	#else{lags.real = lags.dom*frequency(data1.year.norm[[i]])/(timedif.data)} #old
	else{lags.real = lags.dom*timedif.data}
	correlation[,j] = cornmax; cornames[j]= paste("cor_",max(year),sep="") ;j=j+2
	correlation[,k] = lags.real; cornames[k]=paste("lagHOUR_",max(year),sep="") ;k=k+2
	#plotting
	if(plotting==TRUE){
		if(is.na(p1.flag)==T){ #plot approximized data1
			dataplot1.ts = data.frame(Time=index(data1.na),Timeseries = coredata(data1.na), Residuals = data1.in) #use already processed (na.aproxx) data to show in as timeseries
			# dataplot1.ts = data.frame(Time=index(data1.year[[i]]),Timeseries = coredata(data1.year[[i]]), Residuals = data1.in)
		}
		else{ #plot original normalized data1
			dataplot1.ts = data.frame(Time=index(data1.year[[i]]),Timeseries = coredata(data1.year[[i]]), Residuals = data1.in)
		}
		if(is.na(p2.flag)==T){ #plot approximized data2
			dataplot2.ts = data.frame(Time=index(data2.na),Timeseries = coredata(data2.na), Residuals = data2.in) #use already processed (na.aproxx) data to show in as timeseries
			# dataplot2.ts = data.frame(Time=index(data2.year[[i]]),Timeseries = coredata(data2.year[[i]]), Residuals = data2.in)
		}
		else{ #plot original normalized data2
			dataplot2.ts = data.frame(Time=index(data2.year[[i]]),Timeseries = coredata(data2.year[[i]]), Residuals = data2.in)
		}
		# browser()
		# dataplot1.ts = merge(data1.year[[i]],data1.in, all=T, fill=NA); colnames(dataplot1.ts)=c("Time","Timeseries","Residuals")
		# dataplot2.ts = merge(data2.year[[i]],data2.in, all=T, fill=NA); colnames(dataplot2.ts)=c("Time","Timeseries","Residuals")
		# dataplot1.ts = data.frame(Time=index(data1.year[[i]]),Timeseries = coredata(data1.year[[i]]), Residuals = data1.in)
		# dataplot2.ts = data.frame(Time=index(data2.year[[i]]),Timeseries = coredata(data2.year[[i]]), Residuals = data2.in)
	data1.ts.resh = melt(dataplot1.ts, id="Time")
	data2.ts.resh = melt(dataplot2.ts, id="Time")
	# datasets = rbind(cbind(data1.ts.resh, Sensor=factor(deparse(substitute(data1)))),cbind(data2.ts.resh, Sensor=factor(deparse(substitute(data2)))))
	datasets = rbind(cbind(data1.ts.resh, Sensor=factor(sensor1)),cbind(data2.ts.resh, Sensor=factor(sensor2)))
	TS.plot= ggplot(datasets, aes(x=Time, y=value, colour=Sensor)) + geom_line() + facet_grid(variable ~ ., scale="free_y") + xlab("") + ylab("Signal") # + labs(title=paste(deparse(substitute(data1)),"&",deparse(substitute(data2)), ":", max(year), sep=" "))
	dataccf = data.frame(Lag=res$lag[,,1], Correlation=res$acf[,,1])
	ccf.plot= ggplot(dataccf, aes(x=Lag, y=Correlation)) + geom_bar(stat = "identity")
	#combining plots
	# grid.arrange(TS.plot, ccf.plot,main=textGrob(paste(deparse(substitute(data1)),"&",deparse(substitute(data2)), ":", max(year), sep=" "), gp=gpar(cex=1.5)))
	#grid.arrange(TS.plot, ccf.plot,main=textGrob(paste(sensor1,"&", sensor2, ":", max(year), sep=" "), gp=gpar(cex=1.5)))
	grid.arrange(TS.plot, ccf.plot,top=textGrob(paste(sensor1,"&", sensor2, ":", max(year), sep=" "), gp=gpar(cex=1.5)))
	data.plots[[i]] = recordPlot() #store plot in variable
	}
   }
    colnames(correlation) = cornames #set col names as specific years 
    data.plots <<- data.plots
    data1.stl <<- data1.stl
    data2.stl <<- data2.stl
    ccf.results <<-res 
    return(correlation)
}
### end function timeshift


#' @title Event-based calulation of correlation coefficients, lagtimes and signal to noise ratios
#'
#' @description So far hard-coded for SM data of wettzell SGnew TS
#'
#' @param event vector of start and end date of event.
#' @param event_name character, name of event (mainly used for the plot).
#' @param mv_win length of window for moving average to smooth data in order to calculate signal to noise ratios.
#' @param plots logical, should the results be plotted? (default is TRUE).
#' @param nmax number of outputed results per correlation.
#'
#' @details Output is a data.frame with correlation coeficients, lagtime and signal to soise ratios for each sensor.
#' @references Marvin Reich (2015), mreich@@gfz-potsdam.de
#' @examples
#' 

ccf_events = function(event, event_name, mv_win, nmax, stl_flag=F,plots=T,...){
#load data
#soil moisture
load(file="/home/mreich/server/hygra/DataWettzell/SoilMoisture/Cluster_Data/Data_filtered/SGnew_filtered_6hourmean.rdata")
beneathBuilding = SGnew.filter[,11:18] #get sensor beneath SG building
besidesBuilding = SGnew.filter[,19:25] #get closest sensors outside SG building
#GW data
load(file="/home/mreich/server/hygra/DataWettzell/Groundwater/GW_Pegel_complete/BK14_corrected.rdata") #GW 14
BK14.agg = aggregate(BK14, function(x) as.POSIXct(trunc(x, "hour"), tz = "GMT"), mean, na.rm=T)
load(file="/home/mreich/server/hygra/DataWettzell/Groundwater/GW_Pegel_complete/BK3_corrected.rdata") #GW 14
BK3.agg = aggregate(BK3, function(x) as.POSIXct(trunc(x, "hour"), tz = "GMT"), mean, na.rm=T)
gw = merge(BK14.agg,BK3.agg=(BK3.agg-2.01),all=T,fill=NA)
BK14_na = which(is.na(gw$BK14)==T)
gw$BK14[BK14_na] = gw$BK3[BK14_na]
BK14 = gw$BK14 * -1
#precipitation
load("/home/mreich/server/hygra/DataWettzell/Climate/30min_wettzell/clima30min_raw.rdata")
precip = clima.raw$Prec_Sen1
#trim ts to event period
a06_event = besidesBuilding$mux43_05[which(index(besidesBuilding$mux43_05)==event[1]):which(index(besidesBuilding$mux43_05)==event[2])]
a10_event = besidesBuilding$mux43_06[which(index(besidesBuilding$mux43_06)==event[1]):which(index(besidesBuilding$mux43_06)==event[2])]
a14_event = besidesBuilding$mux43_07[which(index(besidesBuilding$mux43_07)==event[1]):which(index(besidesBuilding$mux43_07)==event[2])]
a18_event = besidesBuilding$mux43_08[which(index(besidesBuilding$mux43_08)==event[1]):which(index(besidesBuilding$mux43_08)==event[2])]
b06_event = beneathBuilding$mux44_02[which(index(beneathBuilding$mux44_02)==event[1]):which(index(beneathBuilding$mux44_02)==event[2])]
b16_event = beneathBuilding$mux44_06[which(index(beneathBuilding$mux44_06)==event[1]):which(index(beneathBuilding$mux44_06)==event[2])]
c06_event = beneathBuilding$mux44_01[which(index(beneathBuilding$mux44_01)==event[1]):which(index(beneathBuilding$mux44_01)==event[2])]
c08_event = beneathBuilding$mux44_04[which(index(beneathBuilding$mux44_04)==event[1]):which(index(beneathBuilding$mux44_04)==event[2])]
c14_event = beneathBuilding$mux44_05[which(index(beneathBuilding$mux44_05)==event[1]):which(index(beneathBuilding$mux44_05)==event[2])]
c17_event = beneathBuilding$mux44_07[which(index(beneathBuilding$mux44_07)==event[1]):which(index(beneathBuilding$mux44_07)==event[2])]
d07_event = beneathBuilding$mux44_03[which(index(beneathBuilding$mux44_03)==event[1]):which(index(beneathBuilding$mux44_03)==event[2])]
d18_event = beneathBuilding$mux44_08[which(index(beneathBuilding$mux44_08)==event[1]):which(index(beneathBuilding$mux44_08)==event[2])]
gw_event = BK14[which(index(BK14)==event[1]):which(index(BK14)==event[2])]
precip_event = precip[which(index(precip)==event[1]):which(index(precip)==event[2])]
#merging all dataset for a visual check
ts_event = merge(a06_event,a10_event,a14_event,a18_event,b06_event,b16_event,c06_event,c08_event,c14_event,c17_event,d07_event,d18_event,gw_event,precip_event,all=T,fill=NA)
ts_event$gw_event = na.approx(ts_event$gw_event)
ts_event$precip_event = na.approx(ts_event$precip_event)
colnames(ts_event) = c("a shallow","a middle high","a middle low","a deep","b shallow","b deep","c shallow","c middle high","c middle low","c deep","d shallow","d deep","gw","precip")
if(plots==T) plot(ts_event, xlab="", main=event_name)
#SNR: ratios and filtering
plotting=F
#ratios
b06_event.snr = round(snr_ratio(b06_event,mv_win,T,plotting),2)
c06_event.snr = round(snr_ratio(c06_event,mv_win,T,plotting),2)
d07_event.snr = round(snr_ratio(d07_event,mv_win,T,plotting),2)
c08_event.snr = round(snr_ratio(c08_event,mv_win,T,plotting),2)
c14_event.snr = round(snr_ratio(c14_event,mv_win,T,plotting),2)
b16_event.snr = round(snr_ratio(b16_event,mv_win,T,plotting),2)
c17_event.snr = round(snr_ratio(c17_event,mv_win,T,plotting),2)
d18_event.snr = round(snr_ratio(d18_event,mv_win,T,plotting),2)
#testwise also for profile a
a06_event.snr = round(snr_ratio(a06_event,mv_win,T,plotting),2)
a10_event.snr = round(snr_ratio(a10_event,mv_win,T,plotting),2)
a14_event.snr = round(snr_ratio(a14_event,mv_win,T,plotting),2)
a18_event.snr = round(snr_ratio(a18_event,mv_win,T,plotting),2)
#filtering TS (use mv as signal instead of origial "raw" data
#also filter profile a: apply same filtering on ALL ts included in analysis!!
b06_event.filter = snr_filter(b06_event,mv_win,T,plotting)
c06_event.filter = snr_filter(c06_event,mv_win,T,plotting)
d07_event.filter = snr_filter(d07_event,mv_win,T,plotting)
c08_event.filter = snr_filter(c08_event,mv_win,T,plotting)
c14_event.filter = snr_filter(c14_event,mv_win,T,plotting)
b16_event.filter = snr_filter(b16_event,mv_win,T,plotting)
c17_event.filter = snr_filter(c17_event,mv_win,T,plotting)
d18_event.filter = snr_filter(d18_event,mv_win,T,plotting)
a06_event.filter = snr_filter(a06_event,mv_win,T,plotting)
a10_event.filter = snr_filter(a10_event,mv_win,T,plotting)
a14_event.filter = snr_filter(a14_event,mv_win,T,plotting)
a18_event.filter = snr_filter(a18_event,mv_win,T,plotting)
#Correlation anaylsis
#stl_flag=F;
plotting=F
swin=17521 ;twin=17521
#shallow
a06_b06_event.cor = timeshift(a06_event.filter,b06_event.filter,nmax,T,stl_flag,plotting,swin,twin,"","")
a06_c06_event.cor = timeshift(a06_event.filter,c06_event.filter,nmax,T,stl_flag,plotting,swin,twin,"","")
a06_d07_event.cor = timeshift(a06_event.filter,d07_event.filter,nmax,T,stl_flag,plotting,swin,twin,"","")
#in between sensors (only 2 from profile c)
a10_c08_event.cor = timeshift(a10_event.filter,c08_event.filter,nmax,T,stl_flag,plotting,swin,twin,"","")
a14_c14_event.cor = timeshift(a14_event.filter,c14_event.filter,nmax,T,stl_flag,plotting,swin,twin,"","")
#deep
a18_b16_event.cor = timeshift(a18_event.filter,b16_event.filter,nmax,T,stl_flag,plotting,swin,twin,"","")
a18_c17_event.cor = timeshift(a18_event.filter,c17_event.filter,nmax,T,stl_flag,plotting,swin,twin,"","")
a18_d18_event.cor = timeshift(a18_event.filter,d18_event.filter,nmax,T,stl_flag,plotting,swin,twin,"","")
#correlation results
results_event = rbind(cbind(Profile="shallow",a06_b06_event.cor,sensor="b", depth=0.6, event=event_name, snr=b06_event.snr),
		cbind(Profile="shallow",a06_c06_event.cor,sensor="c", depth=0.6, event=event_name, snr=c06_event.snr),
		cbind(Profile="shallow",a06_d07_event.cor,sensor="d", depth=0.7, event=event_name, snr=d07_event.snr),
		cbind(Profile="middle high",a10_c08_event.cor,sensor="c", depth=0.9, event=event_name, snr=c06_event.snr),
		cbind(Profile="middle low",a14_c14_event.cor,sensor="c", depth=1.4, event=event_name, snr=c14_event.snr),
		cbind(Profile="deep",a18_b16_event.cor,sensor="b", depth=1.6, event=event_name, snr=b16_event.snr),
		cbind(Profile="deep",a18_c17_event.cor,sensor="c", depth=1.7, event=event_name, snr=c17_event.snr),
		cbind(Profile="deep",a18_d18_event.cor,sensor="d", depth=1.8, event=event_name, snr=d18_event.snr))
colnames(results_event)[2:3] = c("correlation","lagtime")
#print snr ratios for profile a
print("SNR ratios for profile a:")
print(paste("Event:",event_name,sep=" "))
print(paste("a06:",a06_event.snr, sep=" "))
print(paste("a10:",a10_event.snr, sep=" "))
print(paste("a14:",a14_event.snr, sep=" "))
print(paste("a18:",a18_event.snr, sep=" "))
#return table of ccf results, together with other meta information and snr-ratios
return(results_event)
} #end of function

#' @title Multi-event based ccf approach, including lagtimes and SNR
#'
#' @description Input events are filtered by SNR and optionally correlation and lagtime thresholds. The filtered datasets are then cross-correlated.
#' @description The resulting correlation coefficents, lagtimes and statistics are plotted and returned in form of a table (data frame).
#' @description So far hard-coded for SM data of wettzell SGnew TS
#'
#' @param events data.frame of events; each row limits one event, event[,1] = start dates, event[,2] = end dates, event[,3] = event name.
#' @param mv_win length of window for moving average to smooth data in order to calculate signal to noise ratios.
#' @param limit_cor threshold for correlation coefficient values; all values below limit_cor will be discarted.
#' @param limit_lag threshold for lagtime values; all values above limit_lag will be discarted.
#' @param correct_cor,correct_lag logical. should these quality filters be applied (TRUE, default) or only be included into statistics (done always).
#' @param plotting logical, should the results be plotted? (default is TRUE).
#' @param nmax number of outputs of correlation.
#'
#' @details Output is a data.frame with mean values for correlation coeficients, lagtime and quality of signal to soise ratio.
#' @references Marvin Reich (2015), mreich@@gfz-potsdam.de
#' @examples
#' 

ccf_multievent = function(events, mv_win, nmax, limit_cor, correct_cor=T, limit_lag, correct_lag=T, plotting=T){

	numberofevents = length(events[,1])
	results_all = data.frame()
	#creating results for each event
	for(i in 1:numberofevents){
		event = c(events[i,1:2])
		name=as.character(events[i,3])
		#res = ccf_events(event, name, mv_win, nmax, stl_flag, F)
		res = cclag_events(event, name, mv_win, nmax, F)
		#filter results for bad snr AND count occurances of bad snr
		res$snrtest = 1; res$cortest = 1; res$lagtest = 1
		#filter bad SNR
		snr_out = which(res$snr < 1)
		#filter bad correlation
		cor_out = which(res$correlation < limit_cor)
		#filter bad lagtime
		lag_out = which(res$lagtime > limit_lag)
		#apply snr filter; not optional
		res$correlation[snr_out] = NA
		res$lagtime[snr_out] = NA
		#apply cor and lag filter; optional
		if(correct_cor == T){
		res$correlation[cor_out] = NA
		res$lagtime[cor_out] = NA
		}
		if(correct_lag == T){
		res$correlation[lag_out] = NA
		res$lagtime[lag_out] = NA
		}
		#for statistics
		res$snrtest[snr_out] = 0
		res$cortest[cor_out] = 0 
		res$lagtest[lag_out] = 0
		#combine actual with all results
		results_all = rbind(results_all, res)
		sensor_metadata = data.frame(Profile=res$Profile,sensor=res$sensor,depth=res$depth)
		sensor_metadata = mutate(sensor_metadata, sensorprofile = paste(sensor,Profile,sep=" "))
	}
	#statisical merging
	#adding unique site-identification (profile + sensor)
	results_ordered = mutate(results_all, sensorprofile = paste(sensor,Profile,sep=" ")) %>%
	 	  group_by(sensorprofile)
	#calculate means and check stats of bad data quality (=snr > limit_snr)
	cor_mean = summarise(results_ordered, cor_mean=mean(correlation, na.rm=T),cor_min=min(correlation, na.rm=T),cor_max=max(correlation, na.rm=T))
	lag_mean = summarise(results_ordered, lag_mean=mean(lagtime, na.rm=T),lag_min=min(lagtime, na.rm=T),lag_max=max(lagtime, na.rm=T))
	snrpass = summarise(results_ordered, SNRpass = sum(snrtest, na.rm=T)/nmax)
	corpass = summarise(results_ordered, CORpass = sum(cortest, na.rm=T)/nmax)
	lagpass = summarise(results_ordered, LAGpass = sum(lagtest, na.rm=T)/nmax)
	filtered_abs = summarise(results_ordered, filtered =  length(which(is.na(correlation)==T))/nmax)
	used_abs = summarise(results_ordered, dataused =  numberofevents - length(which(is.na(correlation)==T))/nmax)
	#putting together with sensor meta data
	results_mean = cbind(cor_mean[,1:2], lag_mean[,2], snrpass[,2], corpass[,2], lagpass[,2], filtered_abs[,2])
	valuerange = cbind(cor_mean[,-2], lag_mean[,3:4])
	results_mean = inner_join(results_mean, sensor_metadata, by="sensorprofile") %>%
		       mutate(SNR = SNRpass / numberofevents) %>%
		       mutate(COR = CORpass / numberofevents) %>%
		       mutate(LAG = LAGpass / numberofevents) %>%
	       	       mutate(EventsUsed = numberofevents - filtered)
	results_table = data.frame(Sensor = results_mean$sensorprofile, Depth = results_mean$depth, Correlation = results_mean$cor_mean, Lagtime = results_mean$lag_mean, SNRpass = results_mean$SNR, CORpass = results_mean$COR,LAGpass = results_mean$LAG, EventsUsed = results_mean$EventsUsed)
	if(plotting ==T){
	#prepare for plotting
	results_mean_plot = melt(results_mean, id=c("Profile","sensor","depth","sensorprofile"))
	#actual plotting
	# creating a transparent theme for all plots
	theme_trans_X = theme(#legend.position = "bottom",
	 legend.position="none",
	 panel.background = element_rect(fill="white"),
	 #panel.background = element_blank(),
	 #panel.border = element_rect(fill="white"),
	 #panel.border = element_blank(),
	 #panel.grid.major = element_line(colour = "black", linetype = "dotted"),
	 panel.grid.minor = element_blank(),
	 panel.grid.major.x = element_blank(),
	 panel.grid.major.y = element_blank(),
	 panel.grid.minor.y = element_blank(),
	 #axis.line = element_line(size = .5, linetype = "solid", colour = "black"),
	 #axis.line.x = element_blank(),
	 axis.title.y=element_text(margin=margin(0,20,0,0)),
	 axis.ticks.x = element_line(size=2),
 	 axis.text=element_text(size=15),
         axis.title=element_text(size=26,face="bold"),
	 strip.text=element_text(size=16,face="bold"),
	 plot.background = element_rect(fill = "transparent",colour = NA)
	 #axis.line = element_line(colour = "black")
	 )
	#correlation
	#results_cor = mutate(results_mean_plot, sensor = factor(sensor, levels=c("d","c","b"))) %>%
	results_cor = mutate(results_mean_plot, sensor = factor(paste("Profile ",sensor,sep=""), levels=c("Profile d","Profile c","Profile b"))) %>%
		      filter(variable == "cor_mean") %>%
		      inner_join(valuerange, by="sensorprofile") %>%
		      inner_join(used_abs, by="sensorprofile") #%>%
		      #group_by(sensor, depth)
	correlation.gg = ggplot(data=results_cor, aes(x=value, y=(depth*-1))) + ylab("Depth [m]") + xlab("Correlation coefficient [%]") + geom_point(size=4.5) + theme_bw() + ggtitle("") +
	geom_errorbarh(aes(xmax=cor_max,xmin=cor_min)) +
	geom_text(aes(label=dataused), vjust=1.6, hjust=1.3) +
	scale_y_continuous(breaks = c(-0.6,-1.0,-1.5,-1.8)) +
	facet_grid(.~sensor) + 
	theme_trans_X + 
	scale_colour_gradientn("",colours=c("black","blue","red"))#, breaks=breaks_val)
	#lagtime
	results_lag = mutate(results_mean_plot, sensor = factor(paste("Profile ",sensor,sep=""), levels=c("Profile d","Profile c","Profile b"))) %>%
		      filter(variable == "lag_mean") %>%
		      inner_join(valuerange, by="sensorprofile") %>%
		      inner_join(used_abs, by="sensorprofile") #%>%
	lagtime.gg = ggplot(data=results_lag, aes(x=value, y=(depth*-1))) + ylab("Depth [m]") + xlab("Lagtime [h]") + geom_point(size=4.5) + theme_bw() + ggtitle("") +
	geom_errorbarh(aes(xmax=lag_max,xmin=lag_min)) +
	geom_text(aes(label=dataused), vjust=1.7, hjust=1.3) +
	scale_y_continuous(breaks = c(-0.6,-1.0,-1.5,-1.8)) +
	facet_grid(.~sensor) + 
	theme_trans_X + 
	scale_colour_gradientn("",colours=c("black","blue","red"))#, breaks=breaks_val)
	#snrtest failure
	results_snr = mutate(results_mean_plot, sensor = factor(paste("Profile ",sensor,sep=""), levels=c("Profile d","Profile c","Profile b"))) %>%
		      filter(variable == "SNR" | variable == "COR" | variable == "LAG")
	snrpass.gg = ggplot(data=results_snr, aes(x=value, y=(depth*-1),colour=variable)) + ylab("Depth [m]") + xlab("Quality criteria passed [%]") + geom_point(size=4.5, position = position_jitter(w=0,h=0.07)) + theme_bw() + ggtitle("") +
	scale_y_continuous(breaks = c(-0.6,-1.0,-1.5,-1.8)) +
	facet_grid(.~sensor) + 
	guides(colour = guide_legend(title="", keywidth=1.5,label.position="bottom")) + 
	theme_trans_X + 
	theme(legend.position="bottom",legend.background = element_rect(fill = "transparent", colour = "transparent"))
	#all plots together
	grid.arrange(correlation.gg,lagtime.gg,snrpass.gg)
	}
	#
	#return table of mean result
	#return(results_table)
	#slightly differently formated data.frame, used for plotting
	#valuerange <<- valuerange
	#used_abs <<- used_abs
	return(results_mean)
}

#' @title Estimation of correlation coefficients and time lags of two timeseries (via cross-correlation )
#'
#' @description THIS IS A NEW VERSION..ojo!! has to be decided and merged afterwards!!!!
#' @description or replace old timeshift()-function !?
#' @description This function is based on standard time series decomposition and afterwards cross-correlation estimation between the two input datasets.
#' @description The main aim is to find dominant time shifts between two time series (e.g. soil moisture, climate data, precipitation stations, etc).
#'
#' @param data1,data2 Timeseries of type .zoo.
#' @param nmax Number of correlation values output.
#' @param norm logical. Use normalized data (default) or not.
#' @param mv_win window width used for signal to noise detection. value is numerical and represents data-timesteps.
#' @param ccfplot logical. Plots for each time series interval the cross-correlation plots. Default is F, own output plots are provided. Should be kept FALSE.
#' @param resplot logical. Plots an overall plot with results of correlation coefficients, lagtimes and signal to noise relationships of the complete datasets.
#' @param sensor1,sensor2 names used for the input datasets for plotting and in the output table.
#' @param ... Further parameters passed to internal functions.

#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples outside.cor = ccf.zoo(besidesBuilding$mux43_04,besidesBuilding$mux43_08,5,T,T,T,24,17521)
#' 

cclag = function(data1,data2,nmax,normdata=T,mv_win=96,ccfplot=F,resplot=F,sensor1="data1", sensor2="data2",...){
  #stucture of function:
	#* merge dataset
	#* calculate snr ratio
	#* continue with snr_signal as time series
	#* normalize data
	#* detrend data
	#* cross-correlate datasets
  Sys.setenv(TZ = "GMT")
  #get both datasets to use the same timestamp-frequency
  data1_sameF = aggregate(data1, function(x) as.POSIXct(trunc(x, "hour"), tz = "GMT"), mean, na.rm=T)
  data2_sameF = aggregate(data2, function(x) as.POSIXct(trunc(x, "hour"), tz = "GMT"), mean, na.rm=T)
  #merging all data together can result in problems when na.approx doesn't find 2 non-NA values for interpolation
  #therefor merge has to be used with option all=F
  #data.merge = merge(data1_sameF, data2_sameF, all=T, fill= NA)
  data.merge = merge(data1_sameF, data2_sameF, all=F)
  #merge the datasets, depending on the analysis time
  #this is necessary to avoid na.approx-problems with too many NA's
  #and on the other hand don't loose any data if stl=F
  #aggregate to hourly time
  #   data1.agg = period.apply(data1, endpoints(data1,"hours"), mean, na.rm=T)
  #   data2.agg = period.apply(data2, endpoints(data2,"hours"), mean, na.rm=T)
  #split in yearly intervalls
  ts.years = endpoints(data.merge, "years")
  data1.year=list()
  data2.year=list()
  for(i in 1:(length(ts.years)-1)){
	data1.year[[i]]=data.merge[(ts.years[i]+1):ts.years[i+1],1]
	data2.year[[i]]=data.merge[(ts.years[i]+1):ts.years[i+1],2]
  }
  #run this loop for all time intervalls (standard = years)
  data_plots = list() # prepare list for possible single time series interval plots
  for(i in 1:(length(ts.years)-1)){
	#use NON-decomposed TS; stl = decomposition was excluded in this newer version!
	#SNR = soil to noise relationship (not ratio; different criteria)
	SNRratio_data1 = snr_ratio(data1.year[[i]],mv_win,T,F) #rmNA=T, plotting=F 
	SNRratio_data2 = snr_ratio(data2.year[[i]],mv_win,T,F) #rmNA=T, plotting=F
	data1_SNRsignal = snr_filter(data1.year[[i]],mv_win,T,F) #rmNA=T, plotting=F
	data2_SNRsignal = snr_filter(data2.year[[i]],mv_win,T,F) #rmNA=T, plotting=F
	#normalize data
  	if(normdata==T){ #use normalize data as input
          #data1.year.norm = lapply(data1.year, function(x) (x - mean(x, na.rm=T))/sd(x,na.rm=T))
          #data2.year.norm = lapply(data2.year, function(x) (x - mean(x, na.rm=T))/sd(x,na.rm=T))
  		data1_year_norm = normalize(data1.year[[i]])
  		data2_year_norm = normalize(data2.year[[i]])
  	}else{ #dont normalize data
  		data1_year_norm = data1_SNRsignal
  		data2_year_norm = data2_SNRsignal
  	}
	#detrend both time series, using linear fit
	#method uses least-squares fit by an regression analysis
	#datasets needs to be approximized
	#if not, lm() returns shorter datasets (if NAs are within the data)
	#this would lead to wrong detrending when resting data(org)-trend!
	#data1.approx = na.approx(coredata(data1_year_norm))
	data1.approx = na.approx(data1_year_norm)
	data1.model = lm(data1.approx ~ I(1:length(data1.approx)))
	data1.trend = predict(data1.model)
	#summary(model)
	#removing the trend from "raw"-data
	data1.in = data1.approx - data1.trend
	#abline(data1.model)
	#lines(data1.detrended, col="blue")
	#data2.approx = na.approx(coredata(data2_year_norm))
	data2.approx = na.approx(data2_year_norm)
	data2.model = lm(data2.approx ~ I(1:length(data2.approx)))
	data2.trend = predict(data2.model)
	#summary(model)
	#removing the trend from "raw"-data
	data2.in = data2.approx - data2.trend
	#abline(data2.model)
	#lines(data2.detrended, col="blue")
	#time frequency of dataset
    	timedif.data = as.numeric(difftime(index(data1.year[[i]][2]),index(data1.year[[i]][1]), units="hours"))
	freq = 24/timedif.data #frequency in hours
	year=unique(format(index(data1.year[[i]]), "%Y")) #get year from current TS
	# cross-correlation
	#browser()
	#input CANNOT be zoo-object!! no error messge but wrong results!
	#OLD: res = ccf(data1.in,data2.in,na.action=na.pass,plot=T,main=max(year))
	res = ccf(coredata(data1.in),coredata(data2.in),na.action=na.pass,plot=F,main=max(year))
	# browser()
    	acf.maxs=NA; cornmax=NA
    	for(ii in 1:nmax){ #finding nth max values
		 if(is.na(sort(res$acf, TRUE)[ii])==T){ cornmax[ii]=NA; acf.maxs[ii]=NA; next}
		 cornmax[ii] = sort(res$acf, TRUE)[ii]
     		 acf.maxs[ii] = which(res$acf == cornmax[ii])
   	}
    	lags.dom = res$lag[acf.maxs]
	#convert output from frequency to time units[hours]
	lags.real = lags.dom*timedif.data
	#construct result output table
	if(i == 1){ results = data.frame(year=rep(year,nmax),SNRdata1=rep(SNRratio_data1,nmax),SNRdata2=rep(SNRratio_data2,nmax),correlation=cornmax,lagtime=lags.real)
	}else{ #in all years / intervalls but the first time
  	res_year = data.frame(year=rep(year,nmax),SNRdata1=rep(SNRratio_data1,nmax),SNRdata2=rep(SNRratio_data2,nmax),correlation=cornmax,lagtime=lags.real) #define data structure to show results
	results = rbind(results, res_year)
	}
	#plotting
	if(ccfplot==TRUE){
		plotdata1_year = merge.zoo(original=data1.year[[i]],residual=data1.in, all=T, fill=NA)
		plotdata2_year = merge.zoo(original=data2.year[[i]],residual=data2.in, all=T, fill=NA)
	data1.ts.resh = melt(zootodf(plotdata1_year), id="time")
	data2.ts.resh = melt(zootodf(plotdata2_year), id="time")
	datasets = rbind(cbind(data1.ts.resh, Sensor=factor(sensor1)),cbind(data2.ts.resh, Sensor=factor(sensor2)))
	TS.plot= ggplot(datasets, aes(x=time, y=value, colour=Sensor)) + geom_line() + facet_grid(variable ~ ., scale="free_y") + xlab("") + ylab("Signal") 
	dataccf = data.frame(Lag=res$lag[,,1], Correlation=res$acf[,,1])
	ccf.plot= ggplot(dataccf, aes(x=Lag, y=Correlation)) + geom_bar(stat = "identity")
	#combining plots
	grid.arrange(TS.plot, ccf.plot,top=textGrob(paste(sensor1,"&", sensor2, ":", max(year), sep=" "), gp=gpar(cex=1.5)))
	data_plots[[i]] = recordPlot() #store plot in variable
	}
   }
   if(resplot == TRUE){
	   #here goes ggplot of result plot
   }
    ccf_plots <<- data_plots
    return(results)
}
### end function cclag

#' @title Estimation of correlation coefficients and time lags of two timeseries (via cross-correlation ) - NO SNR included
#'
#' @description THIS IS A NEW VERSION..ojo!! has to be decided and merged afterwards!!!!
#' @description or replace old timeshift()-function !?
#' @description This function is based on standard time series decomposition and afterwards cross-correlation estimation between the two input datasets.
#' @description The main aim is to find dominant time shifts between two time series (e.g. soil moisture, climate data, precipitation stations, etc).
#'
#' @param data1,data2 Timeseries of type .zoo.
#' @param nmax Number of correlation values output.
#' @param norm logical. Use normalized data (default) or not.
#' @param mv_win window width used for signal to noise detection. value is numerical and represents data-timesteps.
#' @param ccfplot logical. Plots for each time series interval the cross-correlation plots. Default is F, own output plots are provided. Should be kept FALSE.
#' @param resplot logical. Plots an overall plot with results of correlation coefficients, lagtimes and signal to noise relationships of the complete datasets.
#' @param sensor1,sensor2 names used for the input datasets for plotting and in the output table.
#' @param ... Further parameters passed to internal functions.

#' @details missing
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples outside.cor = ccf.zoo(besidesBuilding$mux43_04,besidesBuilding$mux43_08,5,T,T,T,24,17521)
#' 

cclag_nosnr = function(data1,data2,data1.raw,data2.raw,nmax,normdata=T,mv_win=96,ccfplot=F,resplot=F,sensor1="data1", sensor2="data2",...){
  #stucture of function:
	#* merge dataset
	#* calculate snr ratio
	#* continue with snr_signal as time series
	#* normalize data
	#* detrend data
	#* cross-correlate datasets
  Sys.setenv(TZ = "GMT")
  #get both datasets to use the same timestamp-frequency
  data1_sameF = aggregate(data1, function(x) as.POSIXct(trunc(x, "hour"), tz = "GMT"), mean, na.rm=T)
  data2_sameF = aggregate(data2, function(x) as.POSIXct(trunc(x, "hour"), tz = "GMT"), mean, na.rm=T)
  data1_sameF.raw = aggregate(data1.raw, function(x) as.POSIXct(trunc(x, "hour"), tz = "GMT"), mean, na.rm=T)
  data2_sameF.raw = aggregate(data2.raw, function(x) as.POSIXct(trunc(x, "hour"), tz = "GMT"), mean, na.rm=T)
  #merging all data together can result in problems when na.approx doesn't find 2 non-NA values for interpolation
  #therefor merge has to be used with option all=F
  #data.merge = merge(data1_sameF, data2_sameF, all=T, fill= NA)
  data.merge = merge(data1_sameF, data2_sameF, all=F)
  data.merge.raw = merge(data1_sameF.raw, data2_sameF.raw, all=F)
  #merge the datasets, depending on the analysis time
  #this is necessary to avoid na.approx-problems with too many NA's
  #and on the other hand don't loose any data if stl=F
  #aggregate to hourly time
  #   data1.agg = period.apply(data1, endpoints(data1,"hours"), mean, na.rm=T)
  #   data2.agg = period.apply(data2, endpoints(data2,"hours"), mean, na.rm=T)
  #split in yearly intervalls
  ts.years = endpoints(data.merge, "years")
  data1.year=list()
  data2.year=list()
  data1.year.raw=list()
  data2.year.raw=list()
  for(i in 1:(length(ts.years)-1)){
	data1.year[[i]]=data.merge[(ts.years[i]+1):ts.years[i+1],1]
	data2.year[[i]]=data.merge[(ts.years[i]+1):ts.years[i+1],2]
	data1.year.raw[[i]]=data.merge.raw[(ts.years[i]+1):ts.years[i+1],1]
	data2.year.raw[[i]]=data.merge.raw[(ts.years[i]+1):ts.years[i+1],2]
  }
  #run this loop for all time intervalls (standard = years)
  data_plots = list() # prepare list for possible single time series interval plots
  for(i in 1:(length(ts.years)-1)){
	#use NON-decomposed TS; stl = decomposition was excluded in this newer version!

    SNRratio_data1 = snr_ratio_2TS(data1.year.raw[[i]],data1.year[[i]],T,F) #rmNA=T, plotting=F 
    SNRratio_data2 = snr_ratio_2TS(data2.year.raw[[i]],data2.year[[i]],T,F) #rmNA=T, plotting=F

	#data1_signal_sd = sd(data1.year[[i]], na.rm=T)
	#data1_noise_sd = sd(data1.year.raw[[i]], na.rm=T)
	#SNRratio_data1 = data1_noise_sd/data1_signal_sd
	#data2_signal_sd = sd(data2.year[[i]], na.rm=T)
	#data2_noise_sd = sd(data2.year.raw[[i]], na.rm=T)
	#SNRratio_data2 = data2_noise_sd/data2_signal_sd

    # data1_signal_mean = mean(data1.year[[i]], na.rm=T)
    # data1_noise_sd = sd(data1.year.raw[[i]], na.rm=T)
    # SNRratio_data1 = data1_signal_mean/data1_noise_sd
    # data2_signal_mean = mean(data2.year[[i]], na.rm=T)
    # data2_noise_sd = sd(data2.year.raw[[i]], na.rm=T)
    # SNRratio_data2 = data2_signal_mean/data2_noise_sd

	#normalize data
  	if(normdata==T){ #use normalize data as input
          #data1.year.norm = lapply(data1.year, function(x) (x - mean(x, na.rm=T))/sd(x,na.rm=T))
          #data2.year.norm = lapply(data2.year, function(x) (x - mean(x, na.rm=T))/sd(x,na.rm=T))
  		data1_year_norm = normalize(data1.year[[i]])
  		data2_year_norm = normalize(data2.year[[i]])
  	}else{ #dont normalize data
  		data1_year_norm = data1.year[[i]]
  		data2_year_norm = data2.year[[i]]
  	}
	#detrend both time series, using linear fit
	#method uses least-squares fit by an regression analysis
	#datasets needs to be approximized
	#if not, lm() returns shorter datasets (if NAs are within the data)
	#this would lead to wrong detrending when resting data(org)-trend!
	#data1.approx = na.approx(coredata(data1_year_norm))
	data1.approx = na.approx(data1_year_norm)
	data1.model = lm(data1.approx ~ I(1:length(data1.approx)))
	data1.trend = predict(data1.model)
	#summary(model)
	#removing the trend from "raw"-data
	data1.in = data1.approx - data1.trend
	#abline(data1.model)
	#lines(data1.detrended, col="blue")
	#data2.approx = na.approx(coredata(data2_year_norm))
	data2.approx = na.approx(data2_year_norm)
	data2.model = lm(data2.approx ~ I(1:length(data2.approx)))
	data2.trend = predict(data2.model)
	#summary(model)
	#removing the trend from "raw"-data
	data2.in = data2.approx - data2.trend
	#abline(data2.model)
	#lines(data2.detrended, col="blue")
	#time frequency of dataset
    	timedif.data = as.numeric(difftime(index(data1.year[[i]][2]),index(data1.year[[i]][1]), units="hours"))
	freq = 24/timedif.data #frequency in hours
	year=unique(format(index(data1.year[[i]]), "%Y")) #get year from current TS
	# cross-correlation
	#browser()
	#input CANNOT be zoo-object!! no error messge but wrong results!
	#OLD: res = ccf(data1.in,data2.in,na.action=na.pass,plot=T,main=max(year))
	res = ccf(coredata(data1.in),coredata(data2.in),na.action=na.pass,plot=F,main=max(year))
	# browser()
    	acf.maxs=NA; cornmax=NA
    	for(ii in 1:nmax){ #finding nth max values
		 if(is.na(sort(res$acf, TRUE)[ii])==T){ cornmax[ii]=NA; acf.maxs[ii]=NA; next}
		 cornmax[ii] = sort(res$acf, TRUE)[ii]
     		 acf.maxs[ii] = which(res$acf == cornmax[ii])
   	}
    	lags.dom = res$lag[acf.maxs]
	#convert output from frequency to time units[hours]
	lags.real = lags.dom*timedif.data
	#construct result output table
	if(i == 1){ results = data.frame(year=rep(year,nmax),SNRdata1=rep(SNRratio_data1,nmax),SNRdata2=rep(SNRratio_data2,nmax),correlation=cornmax,lagtime=lags.real)
	}else{ #in all years / intervalls but the first time
  	res_year = data.frame(year=rep(year,nmax),SNRdata1=rep(SNRratio_data1,nmax),SNRdata2=rep(SNRratio_data2,nmax),correlation=cornmax,lagtime=lags.real) #define data structure to show results
	results = rbind(results, res_year)
	}
	#plotting
	if(ccfplot==TRUE){
		plotdata1_year = merge.zoo(original=data1.year[[i]],residual=data1.in, all=T, fill=NA)
		plotdata2_year = merge.zoo(original=data2.year[[i]],residual=data2.in, all=T, fill=NA)
	data1.ts.resh = melt(zootodf(plotdata1_year), id="time")
	data2.ts.resh = melt(zootodf(plotdata2_year), id="time")
	datasets = rbind(cbind(data1.ts.resh, Sensor=factor(sensor1)),cbind(data2.ts.resh, Sensor=factor(sensor2)))
	TS.plot= ggplot(datasets, aes(x=time, y=value, colour=Sensor)) + geom_line() + facet_grid(variable ~ ., scale="free_y") + xlab("") + ylab("Signal") 
	dataccf = data.frame(Lag=res$lag[,,1], Correlation=res$acf[,,1])
	ccf.plot= ggplot(dataccf, aes(x=Lag, y=Correlation)) + geom_bar(stat = "identity")
	#combining plots
	grid.arrange(TS.plot, ccf.plot,top=textGrob(paste(sensor1,"&", sensor2, ":", max(year), sep=" "), gp=gpar(cex=1.5)))
	data_plots[[i]] = recordPlot() #store plot in variable
	}
   }
   if(resplot == TRUE){
	   #here goes ggplot of result plot
   }
    ccf_plots <<- data_plots
    return(results)
}
### end function cclag

#' @title Event-based calulation of correlation coefficients, lagtimes and signal to noise ratios
#'
#' @description NEW VERSION of ccf_events !! which uses cclag() instead of timeshift()
#' @description So far hard-coded for SM data of wettzell SGnew TS
#'
#' @param event vector of start and end date of event.
#' @param event_name character, name of event (mainly used for the plot).
#' @param mv_win length of window for moving average to smooth data in order to calculate signal to noise ratios.
#' @param plots logical, should the results be plotted? (default is TRUE).
#' @param nmax number of outputed results per correlation.
#'
#' @details Output is a data.frame with correlation coeficients, lagtime and signal to soise ratios for each sensor.
#' @references Marvin Reich (2015), mreich@@gfz-potsdam.de
#' @examples
#' 

cclag_events = function(event, event_name, mv_win, nmax,plots=T,...){
#load data
#soil moisture
load(file="/home/mreich/server/hygra/DataWettzell/SoilMoisture/Cluster_Data/Data_filtered/SGnew_filtered_6hourmean.rdata")
beneathBuilding = SGnew.filter[,11:18] #get sensor beneath SG building
besidesBuilding = SGnew.filter[,19:25] #get closest sensors outside SG building
#GW data
load(file="/home/mreich/server/hygra/DataWettzell/Groundwater/GW_Pegel_complete/BK14_corrected.rdata") #GW 14
BK14.agg = aggregate(BK14, function(x) as.POSIXct(trunc(x, "hour"), tz = "GMT"), mean, na.rm=T)
load(file="/home/mreich/server/hygra/DataWettzell/Groundwater/GW_Pegel_complete/BK3_corrected.rdata") #GW 14
BK3.agg = aggregate(BK3, function(x) as.POSIXct(trunc(x, "hour"), tz = "GMT"), mean, na.rm=T)
gw = merge(BK14.agg,BK3.agg=(BK3.agg-2.01),all=T,fill=NA)
BK14_na = which(is.na(gw$BK14)==T)
gw$BK14[BK14_na] = gw$BK3[BK14_na]
BK14 = gw$BK14 * -1
#precipitation
load("/home/mreich/server/hygra/DataWettzell/Climate/30min_wettzell/clima30min_raw.rdata")
precip = clima.raw$Prec_Sen1
#trim ts to event period
a06_event = besidesBuilding$mux43_05[which(index(besidesBuilding$mux43_05)==event[1]):which(index(besidesBuilding$mux43_05)==event[2])]
a10_event = besidesBuilding$mux43_06[which(index(besidesBuilding$mux43_06)==event[1]):which(index(besidesBuilding$mux43_06)==event[2])]
a14_event = besidesBuilding$mux43_07[which(index(besidesBuilding$mux43_07)==event[1]):which(index(besidesBuilding$mux43_07)==event[2])]
a18_event = besidesBuilding$mux43_08[which(index(besidesBuilding$mux43_08)==event[1]):which(index(besidesBuilding$mux43_08)==event[2])]
b06_event = beneathBuilding$mux44_02[which(index(beneathBuilding$mux44_02)==event[1]):which(index(beneathBuilding$mux44_02)==event[2])]
b16_event = beneathBuilding$mux44_06[which(index(beneathBuilding$mux44_06)==event[1]):which(index(beneathBuilding$mux44_06)==event[2])]
c06_event = beneathBuilding$mux44_01[which(index(beneathBuilding$mux44_01)==event[1]):which(index(beneathBuilding$mux44_01)==event[2])]
c08_event = beneathBuilding$mux44_04[which(index(beneathBuilding$mux44_04)==event[1]):which(index(beneathBuilding$mux44_04)==event[2])]
c14_event = beneathBuilding$mux44_05[which(index(beneathBuilding$mux44_05)==event[1]):which(index(beneathBuilding$mux44_05)==event[2])]
c17_event = beneathBuilding$mux44_07[which(index(beneathBuilding$mux44_07)==event[1]):which(index(beneathBuilding$mux44_07)==event[2])]
d07_event = beneathBuilding$mux44_03[which(index(beneathBuilding$mux44_03)==event[1]):which(index(beneathBuilding$mux44_03)==event[2])]
d18_event = beneathBuilding$mux44_08[which(index(beneathBuilding$mux44_08)==event[1]):which(index(beneathBuilding$mux44_08)==event[2])]
gw_event = BK14[which(index(BK14)==event[1]):which(index(BK14)==event[2])]
precip_event = precip[which(index(precip)==event[1]):which(index(precip)==event[2])]
#merging all dataset for a visual check
ts_event = merge(a06_event,a10_event,a14_event,a18_event,b06_event,b16_event,c06_event,c08_event,c14_event,c17_event,d07_event,d18_event,gw_event,precip_event,all=T,fill=NA)
ts_event$gw_event = na.approx(ts_event$gw_event)
ts_event$precip_event = na.approx(ts_event$precip_event)
colnames(ts_event) = c("a shallow","a middle high","a middle low","a deep","b shallow","b deep","c shallow","c middle high","c middle low","c deep","d shallow","d deep","gw","precip")
if(plots==T) plot(ts_event, xlab="", main=event_name)
#SNR: ratios and filtering
#Correlation anaylsis
plotccf = F
plotresults = F
normData = T
#shallow
a06_b06_event.cor = cclag(a06_event,b06_event,nmax,normData,mv_win,plotccf,plotresults,"a shallow","b shallow")
a06_c06_event.cor = cclag(a06_event,c06_event,nmax,normData,mv_win,plotccf,plotresults,"a shallow","c shallow")
a06_d07_event.cor = cclag(a06_event,d07_event,nmax,normData,mv_win,plotccf,plotresults,"a shallow","d shallow")
#in between sensors 
a10_c08_event.cor = cclag(a10_event,c08_event,nmax,normData,mv_win,plotccf,plotresults,"a middle high","c middle high")
a14_c14_event.cor = cclag(a14_event,c14_event,nmax,normData,mv_win,plotccf,plotresults,"a middle low","c middle low")
#deep
a18_b16_event.cor = cclag(a18_event,b16_event,nmax,normData,mv_win,plotccf,plotresults,"a deep","b deep")
a18_c17_event.cor = cclag(a18_event,c17_event,nmax,normData,mv_win,plotccf,plotresults,"a deep","c deep")
a18_d18_event.cor = cclag(a18_event,d18_event,nmax,normData,mv_win,plotccf,plotresults,"a deep","d deep")
#correlation results
results_event = rbind(cbind(Profile="shallow",a06_b06_event.cor,sensor="b", depth=0.6, event=event_name),
		cbind(Profile="shallow",a06_c06_event.cor,sensor="c", depth=0.6, event=event_name),
		cbind(Profile="shallow",a06_d07_event.cor,sensor="d", depth=0.7, event=event_name),
		cbind(Profile="middle high",a10_c08_event.cor,sensor="c", depth=0.9, event=event_name),
		cbind(Profile="middle low",a14_c14_event.cor,sensor="c", depth=1.4, event=event_name),
		cbind(Profile="deep",a18_b16_event.cor,sensor="b", depth=1.6, event=event_name),
		cbind(Profile="deep",a18_c17_event.cor,sensor="c", depth=1.7, event=event_name),
		cbind(Profile="deep",a18_d18_event.cor,sensor="d", depth=1.8, event=event_name))

a06_event.snr = results_event$SNRdata1[which(results_event$Profile=="shallow")[1]]
a10_event.snr = results_event$SNRdata1[which(results_event$Profile=="middle high")[1]]
a14_event.snr = results_event$SNRdata1[which(results_event$Profile=="middle low")[1]]
a18_event.snr = results_event$SNRdata1[which(results_event$Profile=="deep")[1]]
#re-structure "results_event" to match other functions/code reading-structures
results = cbind(results_event[,1],results_event[,5:9],results_event[,4])
colnames(results) = c("Profile","correlation","lagtime","sensor","depth","event","snr")
#print snr ratios for profile a
print("SNR ratios for profile a:")
print(paste("Event:",event_name,sep=" "))
print(paste("a06:",a06_event.snr, sep=" "))
print(paste("a10:",a10_event.snr, sep=" "))
print(paste("a14:",a14_event.snr, sep=" "))
print(paste("a18:",a18_event.snr, sep=" "))
#return table of ccf results, together with other meta information and snr-ratios
return(results)
} #end of function



