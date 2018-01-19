library(proftools)
# source("http://bioconductor.org/biocLite.R")
# biocLite("Rgraphviz")



Rprof(filename="profile_outputdata.R")
## some code to be profiled
source("EnSRF_data.R")
Rprof()


# summarize the results
summaryRprof("profile_data.R")
plotProfileCallGraph(readProfileData("profile_data.R"),
                     score = "total")

