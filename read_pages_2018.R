# JHF 201905 added part for seasonal coral extraction
## make the same structure as t35.R for PAGES_2018_01

# 1. Load the files 
pages_dir = '/mnt/climstor/giub/EKF400/assimil_data/proxies/PAGES/PAGES_DB_extraction_201801/'
pages = read.csv2(paste0(pages_dir,'metaMatrix_v1_6_1_Oct-Sep_merged_Ant2k.csv'))
pages_data = read.table(paste0(pages_dir,'dataMatrix_v1_6_1_Oct-Sep_merged_Ant2k.txt'), header=T, sep=";")
# read subannual version for coral data assimilation (JF 201905)
pages_oct_mar = read.csv2(paste0(pages_dir,'metaMatrix_v1_6_1_Oct-Mar_subannual.csv'))
pages_data_oct_mar = read.table(paste0(pages_dir,'dataMatrix_v1_6_1_Oct-Mar_subannual.txt'), header=T, sep=";")
pages_apr_sep = read.csv2(paste0(pages_dir,'metaMatrix_v1_6_1_Apr-Sep_subannual.csv'))
pages_data_apr_sep = read.table(paste0(pages_dir,'dataMatrix_v1_6_1_Apr-Sep_subannual.txt'), header=T, sep=";")

# 2. Combine the colnames of pages_data with pages -> they have the same order luckily
combined_pages = data.frame(colnames(pages_data)[2:length(colnames(pages_data))],pages)
# same for seasonal data (JF 201905)
combined_pages_oct_mar = data.frame(colnames(pages_data_oct_mar)[2:length(colnames(pages_data_oct_mar))],pages_oct_mar)
combined_pages_apr_sep = data.frame(colnames(pages_data_apr_sep)[2:length(colnames(pages_data_apr_sep))],pages_apr_sep)

# 2.1 FDR screened series
fdr <- read.table('/mnt/climstor/giub/EKF400/assimil_data/proxies/PAGES/PAGES_FDR-screened/metadata_2.0.0_PAGES-crit-regional+FDR.txt',sep='\t')
fdr_id <- fdr[6,-1]
fdr_pos <- (as.character(combined_pages[,'paleoData.TSid']) %in% as.character(t(fdr_id)))

# 3. Make lonlat
lon = as.numeric(as.character(combined_pages$longitude))
lat = as.numeric(as.character(combined_pages$latitude))
lon_oct_mar = as.numeric(as.character(combined_pages_oct_mar$longitude))
lat_oct_mar = as.numeric(as.character(combined_pages_oct_mar$latitude))
lon_apr_sep = as.numeric(as.character(combined_pages_apr_sep$longitude))
lat_apr_sep = as.numeric(as.character(combined_pages_apr_sep$latitude))
lonlat = matrix(rbind(lon,lat),2,length(lon))
lonfdr = as.numeric(as.character(combined_pages$longitude[fdr_pos]))
latfdr = as.numeric(as.character(combined_pages$latitude[fdr_pos]))
lonlatfdr = matrix(rbind(lonfdr,latfdr),2,length(lonfdr))
lonlat_oct_mar = matrix(rbind(lon_oct_mar,lat_oct_mar),2,length(lon_oct_mar))
lonlat_apr_sep = matrix(rbind(lon_apr_sep,lat_apr_sep),2,length(lon_apr_sep))

# 4. Plot the location of the proxies
library(maps)
# world
map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-90,90), mar=c(0,0,0,0))
points(lon,lat, col="red", pch=16)
points(lonfdr,latfdr,col='blue',pch=4)
points(lon_oct_mar,lat_oct_mar, col="cyan", pch=1)
points(lon_apr_sep,lat_apr_sep, col="orange", pch=2)
# result: only coral, speleothem and instrumental data in subannual files.

# 5. Elevation
elevation = matrix(-999.99, nrow = nrow(pages_data[1]), ncol=length(pages_data[2:dim(pages_data)[2]]))
elevation_oct_mar = matrix(-999.99, nrow = nrow(pages_data_oct_mar[1]), ncol=length(pages_data_oct_mar[2:dim(pages_data_oct_mar)[2]]))
elevation_apr_sep = matrix(-999.99, nrow = nrow(pages_data_apr_sep[1]), ncol=length(pages_data_apr_sep[2:dim(pages_data_apr_sep)[2]]))

# 6. Archivetype
combined_pages$archiveType = as.character(combined_pages$archiveType)
combined_pages_oct_mar$archiveType = as.character(combined_pages_oct_mar$archiveType)
combined_pages_apr_sep$archiveType = as.character(combined_pages_apr_sep$archiveType)

# 7. Year
# year = matrix(as.numeric(as.character(pages_data[1])), nrow(pages_data[1]),1)
year = as.matrix(pages_data[1], nrow(pages_data[1]),1)
year_oct_mar = as.matrix(pages_data_oct_mar[1], nrow(pages_data_oct_mar[1]),1)
year_apr_sep = as.matrix(pages_data_apr_sep[1], nrow(pages_data_apr_sep[1]),1)

# 8. Create a list and save it
pages_proxies = list(chronologies=as.matrix(pages_data[2:dim(pages_data)[2]]), year=year, lonlat=lonlat, archivetype=combined_pages$archiveType, elevation=elevation)
save(pages_proxies, file="/scratch3/joerg/projects/reuse/data/pages/pages_proxies.RData")

pages_proxies_oct_mar = list(chronologies=as.matrix(pages_data_oct_mar[2:dim(pages_data_oct_mar)[2]]), 
                             year=year_oct_mar, lonlat=lonlat_oct_mar, archivetype=combined_pages_oct_mar$archiveType, 
                             elevation=elevation_oct_mar)
save(pages_proxies_oct_mar, file="/scratch3/joerg/projects/reuse/data/pages/pages_proxies_oct_mar.RData")

pages_proxies_apr_sep = list(chronologies=as.matrix(pages_data_apr_sep[2:dim(pages_data_apr_sep)[2]]), 
                             year=year_apr_sep, lonlat=lonlat_apr_sep, archivetype=combined_pages_apr_sep$archiveType, 
                             elevation=elevation_apr_sep)
save(pages_proxies_apr_sep, file="/scratch3/joerg/projects/reuse/data/pages/pages_proxies_apr_sep.RData")
