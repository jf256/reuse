## make the same structure as t35.R for PAGES_2018_01

# 1. Load the files 
pages = read.csv2('/scratch3/veronika/reuse/assimil_data/pages/Data_export_Veronika/metaMatrix_v1_6_1_Oct-Sep_merged_Ant2k.csv')
pages_data =read.table('/scratch3/veronika/reuse/assimil_data/pages/Data_export_Veronika/dataMatrix_v1_6_1_Oct-Sep_merged_Ant2k.txt', header = T, sep=";")

# 2. Combine the colnames of pages_data with pages -> they have the same order luckily
combined_pages = data.frame(colnames(pages_data)[2:length(colnames(pages_data))],pages)

# 3. Make lonlat
lon = as.numeric(as.character(combined_pages$longitude))
lat = as.numeric(as.character(combined_pages$latitude))
lonlat = matrix(rbind(lon,lat),2,length(lon))

# 4. Plot the location of the proxies
library(maps)
# world
map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-90, 90), mar=c(0,0,0,0))
points(lon,lat, col="red", pch=16)

# 5. Elevation
elevation = matrix(-999.99, nrow = nrow(pages_data[1]), ncol=length(pages_data[2:dim(pages_data)[2]]))

# 6. Archivetype
combined_pages$archiveType = as.character(combined_pages$archiveType)

# 7. Year
# year = matrix(as.numeric(as.character(pages_data[1])), nrow(pages_data[1]),1)
year = as.matrix(pages_data[1], nrow(pages_data[1]),1)

# 8. Create a list and save it
pages_proxies = list(chronologies=as.matrix(pages_data[2:dim(pages_data)[2]]), year=year, lonlat=lonlat, archivetype=combined_pages$archiveType, elevation=elevation)
save(pages_proxies, file ="/scratch3/veronika/reuse/pages_proxies.RData")