#---------------------------------------------------------------------
# Whittaker (part02)
#---------------------------------------------------------------------
# IN:  2 CSV files (from Whittaker part01): ts.csv & pr.csv
#      Size:64800
#
#
# OUT: intermediate csv file (whittaker.csv)
#      Size:64800
#---------------------------------------------------------------------

library(plotbiomes)

#---------------------------------------------------------------------
# reading temperature and precipitation intermediate files
#---------------------------------------------------------------------
tsin <- read.table("ts.csv",header=FALSE)     #AT
prin <- read.table("pr.csv",header=FALSE)     #APP

varsize=64800

sites<-data.frame(MAT=tsin,MAP=prin)  #creating dataframe
colnames(sites) <- c('MAT', 'MAP')    # MAP=precipitation MAT=temperature

df_masked <- sites[!is.na(sites$MAT), ]


mm<-Whittaker_biomes_poly # get original polygon biomes from plotbiomes


#---------------------------------------------------------------------
#masking and preparing data
#---------------------------------------------------------------------
df_masked[df_masked == -9999] <- NA
df_masked$MAP<-df_masked$MAP/10 # biomeplot use cm not mm (precipitation)


#---------------------------------------------------------------------
# Computing Whittaker's biomes
#---------------------------------------------------------------------

library(rgeos)

ff<-data.frame(df_masked$MAT,df_masked$MAP)
ff2<-na.omit(ff)
dat <- SpatialPoints(ff2) #convert to spatialpoints
distances2 <- gDistance(dat, mm, byid = TRUE)   #computing euclidean distance

col_min_index <- apply(distances2, 2, which.min)  #get idx of biomes
col_min_index2<-as.vector(col_min_index)          # convert to vector

not_na_rows <- which(!is.na(ff$df_masked.MAT))
ff[not_na_rows, "whittaker"] <- col_min_index2[1:length(not_na_rows)]

ff <- replace(ff, is.na(ff), -9999)   #fillvalues -9999


#---------------------------------------------------------------------
# save to intermediate csv file to be read by whittaker_plot.ncl-> (part03)
#---------------------------------------------------------------------
write.csv(ff, "whittaker.csv", row.names=FALSE)
