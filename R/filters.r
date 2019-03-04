
#' @title Filtering raw soil moisture data
#'
#' @description Filtering raw soil moisture data with adjustable filter criterias.
#'
#' @param file.input input dataset (zoo)
#' @param val.max/min cuts all data greater/smaller than this value
#' @param val.incr/decr maximum increase/decrease in absolute m³/m³ within timesteps lag.max
#' @param lag.max number of timesteps to consider for maximal increase / decrease
#' @param per.max/min value in percent for maximum / minimum deviation between actual value and a mean value of the time series
#' @param lag.mean period for moving mean calculation
#' @param k splits datasets at every k-th element of lag.mean
#' @details The output will be as zoo-object.
#' @details Default values for val.max=1 and val.min=0.
#' @details Possible values for lag.mean are "hours", "days", "weeks", "months", "quarters" or "years".
#' @references Marvin Reich (2014), mreich@@gfz-potsdam.de
#' @examples example MISSING
#' @export
 
filter_soildata <- function(file.input,val.max=1,val.min=0,val.incr,val.decr,lag.max, per.max, per.min, lag.mean,k){
# source("/home/mreich/R-gravity-package/package09/hygra_test/R/convertzootodataframe.r")
# require(xts)
file.filter = zootodf(file.input)
#filtering
high=which(file.filter[,-1] > val.max, arr.ind=T); high[,2] = high[,2]+1 #unrealistic values upper boundary
low=which(file.filter[,-1] < val.min, arr.ind=T);low[,2] = low[,2]+1 #unrealistic values lower boundary
file.filter[high] = NA
file.filter[low] = NA
num.high=length(high[,1])
num.low=length(low[,1])

# lag.max=10
# num.increase=0; num.decrease=0
#for(comlag in 1:lag.max){
#decrease=which(sapply(file.filter[,-1],diff,comlag) < -val.decr, arr.ind=T); decrease[,1] = decrease[,1]+comlag; decrease[,2] = decrease[,2]+1
#increase=which(sapply(file.filter[,-1],diff,comlag) > val.incr, arr.ind=T); increase[,1] = increase[,1]+comlag; increase[,2] = increase[,2]+1
#file.filter[decrease] = NA
#file.filter[increase] = NA
#num.increase=num.increase + length(increase[,1])
#num.decrease=num.decrease + length(decrease[,1])
#}  

breaks=endpoints(file.filter$time,lag.mean,k) #create breaks in the chosen intervall for splitting up data-analysis (mean calculation)
num.rollingmean = 0 #initialize counting number of corrections
for(i in 2:dim(file.filter)[2]){ #mux (column)-wise
  for(j in 1:(length(breaks)-1)){ #week-wise
    mean.actual = mean(file.filter[breaks[j]:breaks[j+1],i], na.rm=T)
    filter.out = which(file.filter[breaks[j]:breaks[j+1],i]/mean.actual > per.max | file.filter[breaks[j]:breaks[j+1],i]/mean.actual < per.min, arr.ind=T)
    file.filter[(filter.out+breaks[j]-1),i] = NA
    num.rollingmean = num.rollingmean + length(filter.out)
  }
}

# stats = data.frame(number=c(num.high,num.low,num.increase,num.decrease,num.rollingmean)); rownames(stats) = c("greater max", "lower min", "increasing too fast", "decreasing too fast", "deviation from moving mean")
stats = data.frame(number=c(num.high,num.low,num.rollingmean)); rownames(stats) = c("greater max", "lower min", "deviation from moving mean")
print(stats)

file.out = read.zoo(file.filter,format="%Y-%m-%d %H:%M:%S",tz="GMT",index.column=1) #convert back to .zoo object
return(file.out)
}

### end filtering soil data ###

#' @title Signal to noise ratio
#'
#' @description Calculates signal to noise ratio (SNR) based on moving average window signal correction.
#'
#' @param data.in input time series dataset (zoo)
#' @param mv_win width of moving average window, in observation time steps (no explicit time unit)
#' @param rmNA logical. should NA values be remove within the forming of mean values? default is TRUE
#' @param plotting logical. if true (not default) datasets are plotted
#' @details The output is a numeric number, representative for the complete input time period.
#' @references Marvin Reich (2015), mreich@@gfz-potsdam.de
#' @examples example MISSING
#' @export
#' 

snr_ratio = function(data.in, mv_win, rmNA = T, plotting=F){
	data_mv = rollapply(data.in,mv_win, mean, na.rm=rmNA) #mv of 24 hours
	dataset = merge(data.in,data_mv,all=F)
	colnames(dataset) = c("raw","signal")
	dataset$noise = dataset$raw - dataset$signal
	signal_amp = abs(range(dataset$signal, na.rm=T)[1] - range(dataset$signal, na.rm=T)[2])
	noise_sd = sd(dataset$noise, na.rm=T)
	#noise_amp = abs(range(dataset$noise, na.rm=T)[1] - range(dataset$noise, na.rm=T)[2])
	#data_snr = signal_amp/noise_amp
	snr_ratio = signal_amp / (2*noise_sd)
	if(snr_ratio > 1) snr_pass = 1 #snr test passed, signal seems reasonable
	else snr_pass = 0 #snr test failed, signal is too small
	if(plotting==T) plot(dataset,xlab = "", main=paste(deparse(substitute(data.in)),": SNR = ",format(snr_ratio,digits=2)," (mavg window: ",mv_win/4," hours)",sep=""))
	return(snr_ratio)
}

#' @title Signal to noise ratio (inputing directly 2 time series (raw and signal)
#'
#' @description Calculates signal to noise ratio (SNR) based on moving average window signal correction.
#'
#' @param data.in input time series dataset (zoo)
#' @param mv_win width of moving average window, in observation time steps (no explicit time unit)
#' @param rmNA logical. should NA values be remove within the forming of mean values? default is TRUE
#' @param plotting logical. if true (not default) datasets are plotted
#' @details The output is a numeric number, representative for the complete input time period.
#' @references Marvin Reich (2015), mreich@@gfz-potsdam.de
#' @examples example MISSING
#' @export
#' 

snr_ratio_2TS = function(data.in.raw, data.in.signal, rmNA = T, plotting=F){
	dataset = merge(data.in.raw,data.in.signal,all=F)
	colnames(dataset) = c("raw","signal")
	dataset$noise = dataset$raw - dataset$signal
	signal_amp = abs(range(dataset$signal, na.rm=T)[1] - range(dataset$signal, na.rm=T)[2])
	signal_sd = sd(dataset$signal, na.rm=T)
	noise_sd = sd(dataset$noise, na.rm=T)
	#noise_amp = abs(range(dataset$noise, na.rm=T)[1] - range(dataset$noise, na.rm=T)[2])
	#data_snr = signal_amp/noise_amp
	snr_ratio = signal_amp / (2*noise_sd)
	#snr_ratio = noise_sd/signal_sd

	if(snr_ratio > 1) snr_pass = 1 #snr test passed, signal seems reasonable
	else snr_pass = 0 #snr test failed, signal is too small
	if(plotting==T) plot(dataset,xlab = "", main=paste(deparse(substitute(data.in.raw)),": SNR = ",format(snr_ratio,digits=2),sep=""))
	return(snr_ratio)
}

#' @title Signal to noise filtering 
#'
#' @description Filters input data based on signal to noise ratios (SNR).
#'
#' @param data.in input time series dataset (zoo)
#' @param mv_win width of moving average window, in observation time steps (no explicit time unit)
#' @param rmNA logical. should NA values be remove within the forming of mean values? default is TRUE
#' @param plotting logical. if true (not default) datasets are plotted
#' @details The output is a time series: signal =  the input time series (raw data) - noise.
#' @references Marvin Reich (2015), mreich@@gfz-potsdam.de
#' @examples example MISSING
#' @export
#' 

snr_filter = function(data.in, mv_win, rmNA = T, plotting=F){
	data_mv = rollapply(data.in,mv_win, mean, na.rm=rmNA) #mv of 24 hours
	dataset = merge(data.in,data_mv,all=F)
	colnames(dataset) = c("raw","signal")
	dataset$noise = dataset$raw - dataset$signal
	signal_amp = abs(range(dataset$signal, na.rm=T)[1] - range(dataset$signal, na.rm=T)[2])
	noise_sd = sd(dataset$noise, na.rm=T)
	snr_ratio = signal_amp / (2*noise_sd)
	#if(snr_ratio > 1) snr_signal = dataset$signal #snr test passed, signal seems reasonable
	#else snr_signal = data.in #snr test failed, signal is too small
	snr_signal = dataset$signal
	if(plotting==T) plot(dataset,xlab = "", main=paste(deparse(substitute(data.in)),": SNR = ",format(snr_ratio,digits=2)," (mavg window: ",mv_win/4," hours)",sep=""))
	return(snr_signal)
}
