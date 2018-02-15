####################
## update man pages and NAMESPACE file
####################
library(devtools)
# load_all("/home/mreich/Dokumente/written/ResearchRepos/UmbrellaEffect")
library(roxygen2)
setwd("/home/mreich/Dokumente/written/ResearchRepos/UmbrellaEffect")
devtools::document()
