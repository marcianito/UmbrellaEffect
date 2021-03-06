% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filters.r
\name{snr_ratio_2TS}
\alias{snr_ratio_2TS}
\title{Signal to noise ratio (inputing directly 2 time series (raw and signal)}
\usage{
snr_ratio_2TS(data.in.raw, data.in.signal, rmNA = T, plotting = F)
}
\arguments{
\item{rmNA}{logical. should NA values be remove within the forming of mean values? default is TRUE}

\item{plotting}{logical. if true (not default) datasets are plotted}

\item{data.in}{input time series dataset (zoo)}

\item{mv_win}{width of moving average window, in observation time steps (no explicit time unit)}
}
\description{
Calculates signal to noise ratio (SNR) based on moving average window signal correction.
}
\details{
The output is a numeric number, representative for the complete input time period.
}
\examples{
example MISSING

}
\references{
Marvin Reich (2015), mreich@gfz-potsdam.de
}
