% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hydrus.r
\name{grid_aggregate}
\alias{grid_aggregate}
\title{aggregate data from hourly to daily (weekly, monthly) averages}
\usage{
grid_aggregate(input, win_width, every_n)
}
\arguments{
\item{data1, data2, data3}{input data. has to be in the format XX.}

\item{timestep}{vector of timesteps to be considered for calculations.}

\item{valmax, valmin}{values where to cut off resulting differces, leaving a data-area of special consideration}
}
\description{
using sorting and a moving window average
}
\details{
missing
}
\examples{
missing
}
\references{
Marvin Reich (2016), mreich@gfz-potsdam.de
}
