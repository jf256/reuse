# copy manually to climstor after generation because script has no permission to do so

setwd('/scratch3/joerg/projects/reuse/git')
wd <- getwd()
setwd('/mnt/climstor/giub/EKF400/assimil_data/v2/docum_v2_2019')
#setwd('/Users/joerg/unibe/projects/reuse/data/docum_v2_2019')
year_min <- 1600
year_max <- 2019
n_stations <- 12

x.data <- array(dim=c((year_max-year_min+1)*12,n_stations))
x.lon <- array(dim=c(n_stations))
x.lat <- array(dim=c(n_stations))
x.name <- array(data="temp2",dim=c(n_stations))
x.height <- array(dim=c(n_stations))
x.time <- seq(year_min + 1/24, by=1/12, length=nrow(x.data))
i_staz <- 0


### READ
files <- read.table("list_temp")
for (filename in t(files)) {
  i_staz <- i_staz+1
#  header <- readLines(filename,n=3)
  header <- read.table(filename,skip=34,nrow=1,sep='\t',comment.char="|")
  lon <- read.table(filename,skip=23,nrow=1,sep='\t',comment.char="|")[2]
  lat <- read.table(filename,skip=22,nrow=1,sep='\t',comment.char="|")[2]
#  x.lon[i_staz] <- as.numeric(substr(header[2],9,12))
#  x.lat[i_staz] <- as.numeric(substr(header[3],9,12))
  x.lon[i_staz] <- as.numeric(lon) #as.numeric(header[1,2])
  x.lat[i_staz] <- as.numeric(lat) #as.numeric(header[2,2])
#  x.height[i_staz] <- as.integer(substr(header[3],49,52))
  data <- read.table(filename,na.strings="-999.9",skip=35)
#  # normalize
#  data <- cbind(data[,1],scale(data[,2:13],center=T,scale=T))
  i_data <- (data[1,1]-year_min)*12+1
  for (i in 1:dim(data)[1]) {
    for (j in 2:13) {
      if (i_data>0) x.data[i_data,i_staz] <- data[i,j]
      i_data <- i_data+1
    }
  }
}

print(dim(x.data))
t <- list(data=scale(x.data), lon=x.lon, lat=x.lat, names=x.name, height=x.height, time=x.time)
setwd(wd)
save(t,file="t_docu_monthly_angie2019.Rdata")
#write.table(cbind(x.lon,x.lat),file="list_coord",sep="\t",col.names=FALSE,row.names=FALSE)
#system("./plot_stations.gmt")
