
# Open libraries
library(sp)
library(dismo)
library(raster)
library(GISTools)
library(rgdal)
library(maptools)
library(rgeos)
library(rgl)
library(rasterVis)
library( rgdal )                                    
library( adehabitatHR )                             
library( maptools )                                
library( sp )
library( move )
   
# Define workspace and folder sturcture
work_wd <- "D:/OneDrive/_BioMoveGISKurs/_BioMove_GISIntro_2019/Geodata"
      # relative to wd
raster_wd <- paste(work_wd,"/","rasterfiles",sep='')
shapes_wd <- paste(work_wd,"/","shapefiles",sep='')
output_wd <- paste(work_wd,"/","output",sep='')

setwd(shapes_wd)                                      

#load data
Tiere <- readOGR("DieTiere.shp", layer="DieTiere")  # animals
Landn <- readOGR("Landnutzung.shp", layer="Landnutzung",stringsAsFactors=FALSE)   #land use
River <- readOGR("linearegew.shp", layer="linearegew",stringsAsFactors=FALSE)   #waterways
Kettl <- read.table('KettleHoles.txt', header=TRUE, sep='\t') # sep='\t' tab delimited

head(Kettl); class(Kettl)  # data.frame

# convert txt-file to SpatialPointsDataFrame
coordinates(Kettl) <- ~MyX + MyY
head(Kettl); class(Kettl) # ok 

                                    
    
# 0 - ########################################################################## 
## - here no sf adjustment as all neded commands are availabe in script part 1 ## sf update 2019
### Map embellishment - a colour code for column NS1
mycol <- c('red','antiquewhite','chartreuse4','cornsilk4','blue')

plot(Landn[Landn$NS1 == 1,],col= mycol[1],border='transparent')
plot(Landn[Landn$NS1 == 2,],col= mycol[2],border='transparent',add=T)
plot(Landn[Landn$NS1 == 3,],col= mycol[3],border='transparent',add=T)
plot(Landn[Landn$NS1 == 4,],col= mycol[4],border='transparent',add=T)
plot(Landn[Landn$NS1 == 5,],col= mycol[5],border='transparent',add=T)

plot(Tiere,pch=Tiere[[1]],add=T) ###oops-where are they??
identicalCRS(Landn,Tiere)        # FALSE...

plot(River, col= 'blue', lwd=0.1,add=T)
identicalCRS(Landn,River)

plot(Kettl, col= 'lightblue',cex=2, pch=15,add=T)

# show projection
proj4string(Tiere)   
proj4string(Landn)  # = crs(Landn) = projection(Landn)
crs(River)
crs(Kettl)  # aesj - not even defined

# Tiere and River have the same - what about KettleHoles?
head(coordinates(River)[[1]])
head(coordinates(Landn))

head(coordinates(Kettl))  # seems to be the same as River!

# thus, we transform the Landn and define the projection for Kettl :
# Transform data
crs_t <- proj4string(River) #define the CRS
Landn <- spTransform(Landn,CRS(crs_t)) # project the shapefile, !Landn overwritten!
crs(Kettl) <-  crs_t # define the projection

# replot everything
plot(Landn[Landn$NS1 == 1,],col= mycol[1],border='transparent')
plot(Landn[Landn$NS1 == 2,],col= mycol[2],border='transparent',add=T)
plot(Landn[Landn$NS1 == 3,],col= mycol[3],border='transparent',add=T)
plot(Landn[Landn$NS1 == 4,],col= mycol[4],border='transparent',add=T)
plot(Landn[Landn$NS1 == 5,],col= mycol[5],border='transparent',add=T)
plot(Tiere,pch=Tiere[[1]],add=T) ### here we go!
plot(River, col= 'blue', lwd=0.1,add=T)
plot(Kettl, col= 'cyan',cex=2, pch=19,add=T)



# 1 - ##########################################################################
####################  Biodiversity Index #######################################

# create a buffer around kettleholes

?gBuffer

myBuf <- gBuffer(spgeom=Kettl, byid=TRUE, width=1000, quadsegs=50)
plot(myBuf, col='cyan', add=T)
plot(Kettl, col= 'pink', cex=2, pch=19, add=T)

head(myBuf) 

myBufLandn <- gIntersection(spgeom1=Landn, spgeom2=myBuf, byid=TRUE, drop_not_poly=TRUE) 
# gIntersection does not keep  the data.frame !   What a sh....
# there is also crop, intersect, gIntersects... None keeps both dataframes or 
# separates overlapping buffers...

df <- data.frame(id = getSpPPolygonsIDSlots(myBufLandn))
row.names(df) <- getSpPPolygonsIDSlots(myBufLandn)
spdf <- SpatialPolygonsDataFrame(myBufLandn, data =df)  ######## this is my new buffered and clipped polygon
plot(spdf,col='lightblue')
spdf@data

### the following is retrieving step by step the corresponding data frames of Landn and Kettl
myids <- getSpPPolygonsIDSlots(myBufLandn) # length 22 polygons
myids2=as.numeric(as.character(unlist(strsplit(myids, ' ')))) #now 44
seq_uneven <- seq(1,43,by=2)
seq_even   <- seq(2,44,by=2)
id_Landn <- myids2[c(seq_uneven)]
id_Kettl <- myids2[c(seq_even)]
mydf <- data.frame(id_Kettl,id_Landn)
spdf@data <- mydf
area_km2 <- round(gArea(spdf, byid=TRUE)/1e6,digits=4)
spdf@data$area_km2 = area_km2

a <- match(spdf@data$id_Landn,as.numeric(row.names(Landn@data))) 
newdat <- Landn@data[a,]
spdf@data <- data.frame(spdf@data,newdat) 


### calculate the sum of the area of each landuse NS1 per kettlehole
for (i in 1 : 3)
{
  myKettl <- subset(spdf@data,spdf@data$id_Kettl == i)
  res <- tapply(X=myKettl$area_km2, INDEX=myKettl$NS1,FUN=sum)
  print(i);print(res); print(sum(res))
} 

# I am giving up here...
# forget about rivers, forget about clipping in R .....
# if you need to do it in R, the install SAGA-GIS, the RSAGA package and use the
# rsaga.geoprocessor interface!

# compare with results  done in ArcGIS
writeOGR(obj=spdf, dsn=output_wd, layer='R_BufInter_Landn',
         driver ='ESRI Shapefile',overwrite=TRUE)

library(foreign)
agf <- read.dbf(paste(output_wd,'/BufferIntersect_Landcov.dbf',sep=''))  #agf =your ArcGisFile
head(agf)

for (i in 1 : 3)
{
  myKettl <- subset(agf,agf$No == i) # adapt 'No' ?
  res <- tapply(X=myKettl$new_area2/1e6, INDEX=myKettl$NS1,FUN=sum) #adapt 'new_area2'
  print(i);print(res); print(sum(res))
} 

# A = pi r^2, so the sum of the area per Buffer should be  exactly 3.14 km² !

# so, both ways are fine, but R is really complicated here!!! 

################################################################################


# 1.B - ##########################################################################
####################  Biodiversity Index (SF version) #######################################  ## sf update 2019

# create a buffer around kettleholes

Kettl.sf <- as(Kettl, "sf")
myBuf.sf <- st_buffer(x = Kettl.sf, dist = 1000, nQuadSegs = 50)

plot(myBuf.sf, col='cyan', add=T)
plot(Kettl, col= 'pink', cex=2, pch=19, add=T)

head(myBuf.sf) 
Landn.sf <- as(Landn, "sf")
myBufLandn.sf <- st_intersection(Landn.sf, myBuf.sf) 

## st_intersection does keep the data! 

### calculate the sum of the area of each landuse NS1 per kettlehole using sf
class(myBufLandn.sf$NS1)
class(myBufLandn.sf$NS1)

## compute the area of each polygon
myBufLandn.sf$area <- st_area(myBufLandn.sf)

## summarize the area of each LU class within each buffer
library(dplyr)
my.LU.result <- myBufLandn.sf %>% ## data 
  group_by(No, NS1) %>% ## group the data by the columns 'No' and "NS1" 
  summarise( area = sum(area),
             SPHEROID = first(SPHEROID),
             PERIMETER = mean(PERIMETER)) ## summarize all columns names using the function named


## plot the result as control
## plot all buffer
plot(my.LU.result[, "area"]) 
## plot each buffer separate 
plot(my.LU.result[my.LU.result$No == 1, "area"]) 
plot(my.LU.result[my.LU.result$No == 2, "area"])
plot(my.LU.result[my.LU.result$No == 3, "area"])


################################################################################






# 2 ############################################################################
### Animal Home Ranges    ######################################################

# How many locations
table(Tiere@data$ID)


# MCP
cp <- mcp(Tiere[,1], percent=95) # Minimum Convex Polygon (95%) per TierID # set unout for units!
plot(cp)                         
plot(Tiere, add=TRUE, col= Tiere$ID)            

# nice diagram! 
hrs <- mcp.area(Tiere[,1], percent=seq(50, 100, by = 5),unout='km2') 
hrs                                                      

writePolyShape(cp, paste(output_wd,"/My_MCP95",sep=''))     # another save command
# take care! Has no prj!

# Kernel density
kud <- kernelUD(Tiere[,1], h="href")               # Kernel berechnen mit h="href"
       image(kud)                                        # den Kernel anzeigen
       kud$'1'@h                                          # Parameter 'h' für ID = 1 anzeigen
kernel.area(kud,unout='km2')
     
    
kud1 <- kernelUD(Tiere[,1], h="LSCV")               # Kernel berechnen mit h="LSCV"
image(kud1)
kernel.area(kud1,unout='km2') ## LSCV

# Umwandeln von density in contour
homerange1 <- getverticeshr(kud1,percent = 90)      # 90% Home Range für kudl
gArea(homerange1,byid=T)/1e6
homerange2 <- getverticeshr(kud,percent = 90)
gArea(homerange2,byid=T)/1e6                        # same size?

plot(Landn,border='transparent')      
plot(homerange1, border='red',col='transparent',add=T)                          # Home Range anzeigen
plot(homerange2, border='blue',col='transparent',add=T)                          # Home Range anzeigen

# density plot as contour
vud <- getvolumeUD(kud1)                             

     image(kud1[[1]])
       title("Output of kernelUD")
     xyz <- as.image.SpatialGridDataFrame(kud1[[1]])   
       contour(xyz, add=TRUE)                           

     image(vud[[1]])
     title("Output of getvolumeUD")
     xyzv <- as.image.SpatialGridDataFrame(vud[[1]])
     contour(xyzv, add=TRUE)


writePolyShape(homerange1, paste(output_wd,"/kernel_hr1",sep=''))   
writePolyShape(homerange2, paste(output_wd,"/kernel_hr2",sep=''))  

#LoCoH nd Brownian Bridge in move-package


#### density as 3D surface plot ------------------------------------------------
install.packages("devtools")
devtools::install_github("ropensci/plotly")
library(plotly)
kd <- with(MASS::geyser, MASS::kde2d(duration, waiting, n = 50))
persp(kd)
with(kd, plot_ly(x = x, y = y, z = z, type = "surface"))
str(kd)
#List of 3
# $ x: num [1:50] 0.833 0.928 1.022 1.116 1.21 ...
# $ y: num [1:50] 43 44.3 45.7 47 48.3 ...
# $ z: num [1:50, 1:50] 9.07e-13 1.81e-12 3.43e-12 6.11e-12 1.03e-11 ...

xy <- coordinates(kud[[1]])
z <- kud[[1]]@data$ud  

df_kud <- data.frame(x=xy[,1],y=xy[,2],z=z)
persp(x=df_kud$x, y=df_kud$y, z= df_kud$z) 
plot3d(x=xy[,1],y=xy[,2],z=z)

df_kud_pt <- df_kud

kud[[1]]@grid@cells.dim
r <- raster(ncols=60,nrows=58)
coordinates(df_kud_pt) <- ~x + y 
r_kud <- rasterize(df_kud_pt,r,'z',fun=max)


myz <- matrix(z ,nrow= kud[[1]]@grid@cells.dim[[2]], 
                    ncol= kud[[1]]@grid@cells.dim[[1]],byrow=F)
image(myz)

kd.list <- list(x=xy[,1],y=xy[,2],z=myz)
with(kd.list, plot_ly(x = x, y = y, z = z, type = "surface"))

plot_ly(z = myz, type = "surface")






# 3 - ##########################################################################
### Preference of LandUse categories ########################################### (SF version not yet added;  ## sf update 2019 )
### Spatial Join
spatjoin <- over(Tiere,Landn)
Tiere_Landn <- cbind(Tiere@data,spatjoin)  
class(Tiere_Landn)

### Analyse: Landnutzungskategorie per location of animal 
myres <- table(Tiere_Landn$NS1,Tiere_Landn$ID)
myres
#     1   2   3
#  2  48  25  95
#  3 147  64   9

barplot(myres,beside=T)

### create random points in home range 
# My_MCP95.shp has no CRS   
# shapefile MCP einlesen
MCP95 <- readOGR(dsn=output_wd,layer="My_MCP95",stringsAsFactors=FALSE) 
crs_tmp <- proj4string(Tiere)
proj4string(MCP95) <- crs_tmp     
summary(MCP95) # oder:   MCP_t@data

plot(Landn,border='transparent')      
plot(MCP95, add=T, border = 'grey')

#locations per animal 
table(Tiere[[1]])   
  #   1   2   3 
  # 195  89 104 

# Gleiche Anzahl an Zufallspunkten ins MCP des ersten Tieres:
rp1 <- spsample(MCP95[1,],n=195,type='random') #erste Zeile [1,] ist ID = 1
  points(rp1,col='red',pch=3)
  class(rp1)
  
ID_r <- as.data.frame(cbind(rep(1, length(rp1)))) #data.frame kreieren, um aus SpatialPoints object
                                # ein SpatialPointsDataFrame zu machen (Anhängen der attribute table)
  
rp1_final <- SpatialPointsDataFrame(coords=coordinates(rp1),data=ID_r,proj4string = CRS(crs_t))

### erneut Spatial Join
spatjoin_rp <- over(rp1_final,Landn)
rp_Landn <- cbind(rp1_final@data,spatjoin_rp) #### !!!!!!!!!!!! 
#class(rp_Landn)

### Analyse: Landnutzungskategorie pro Peilungspunkt und Tier
myres_rp <- table(rp_Landn$NS1)
myres_rp
#  2   3 
# 136  59 

#  Tier 1 im Vergleich:  
#  2  48  
#  3 147

# 3  - #########################################################################
# Was wir wollen:

                  #LANDNUTZUNGSTYP
# TYP             #Feld(NS1=2)      #Wald(NS1=3)
#Angebot(random)     136              59
#Tier (observed)      48             147

chi_matrix <- matrix(data=c(136,59,48,147),ncol=2)
chisq.test(chi_matrix)

# 'von Hand' ausrechnen:
aa <- chi_matrix[1,1]
bb <- chi_matrix[1,2]
cc <- chi_matrix[2,1]
dd <- chi_matrix[2,2]

Na <- aa+bb
Nb <- cc+dd
Ns <- aa+cc
Nf <- bb+dd
N <- Na+Nb

mychi <- (N*(abs((aa*dd) - (bb*cc)) - (N/2))^2) / (Ns*Nf*Na*Nb)
mychi

################################################################################



# 4 - ##########################################################################
#### AnimalMovement

library(move)
getwd()  # zum Überprüfen
Fuchs <- read.table("Fuxi.txt", header=TRUE, sep=',')            
head(Fuchs)
class(Fuchs) # Achtung! not spatial, not a move Object!

coordinates(Fuchs) <-  ~lon + lat
crs_kurs <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")    # geographical,datum WGS84
proj4string(Fuchs) <- crs_kurs     # Definiere Projektion der Daten
summary(Fuchs)

plot(Fuchs, add=T)  # warum funktioniert es nicht????? -> Antwort siehe oben!
Fuchs_t <- spTransform(Fuchs,crs(Tiere))
plot(Fuchs_t,pch=18,col='grey',cex=1,add=T)
zoom(Landn)

# mal abspeichern
writeOGR(obj=Fuchs_t,dsn=output_wd, layer='Fuxi_aus_R',driver ='ESRI Shapefile',overwrite=TRUE)

# create move object   - ad a datetime stamp DT

MoveFuchs <- read.table("Fuxi.txt", header=TRUE, sep=',')      
dates <- as.character(MoveFuchs$Date)
times <- as.character(MoveFuchs$Time)
DateTime <- paste(dates, times)
DT <- strptime(DateTime, "%m/%d/%Y %H:%M:%S")
MoveFuchs <- cbind(MoveFuchs, DT)
head(MoveFuchs)

#  Datum WGS84      Lat      Lon       Date     Time                  DT
#1    TP     D 52.91443 13.43366 03/18/2015 09:14:29 2015-03-18 09:14:29
#2    TP     D 52.91448 13.43363 03/18/2015 09:15:45 2015-03-18 09:15:45
#3    TP     D 52.91448 13.43363 03/18/2015 09:16:47 2015-03-18 09:16:47
#4    TP     D 52.91448 13.43363 03/18/2015 09:17:17 2015-03-18 09:17:17
#5    TP     D 52.91449 13.43363 03/18/2015 09:18:02 2015-03-18 09:18:02
#6    TP     D 52.91449 13.43363 03/18/2015 09:18:39 2015-03-18 09:18:39
 

my_fox <- move(x=MoveFuchs$lon, y= MoveFuchs$lat, 
time=as.POSIXct(MoveFuchs$DT,format="%Y-%m-%d %H:%M:%S"),
proj= CRS("+proj=longlat"),data=MoveFuchs,animal="FuchsX"  )  #now a Move Object!

show(my_fox)
plot(my_fox, type="o", col=3, lwd=2, pch=20, xlab="location_long",
 ylab="location_lat")


#  Steplength and turning angles
install.packages("circular") 
library(circular)

turnang <- angle(my_fox)
rose.diag(turnang, bins=24,shrink=0.6,xlim=c(-2,2),ylim=c(-2,2),axes=T,prop=2,units='degrees')

steplength <- distance(my_fox)
hist(steplength)

zespeed <- speed(my_fox)
hist(zespeed)



####################    THE END   ##############################################