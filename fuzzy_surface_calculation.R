################################################
### prerequisites
################################################

#check for existance of necessary packages and if they are not present in R instalation
#then install them
list.of.packages <- c("gstat", "qualityTools", "rgdal", "sp", 
                      "reshape", "ggplot2", "maptools", "MASS")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
#load all the packages
lapply(list.of.packages,function(x){library(x,character.only=TRUE)}) 
#remove the variables related to packages
remove(list.of.packages, new.packages)

#load additional functions
source("_functions.R", encoding = "UTF-8")

#to load all the for the analysis please uncoment the line below
#load(file = "all_data_precalculated.RData")



################################################
### data loading and processing
################################################

#load the dataset of stations with latitude, longitude and PM10 concetrations for October 2013
stations <- read.table("data/pm10-october-2013.txt", header=TRUE, quote="\"", stringsAsFactors=FALSE)

#remove stations without measurement and duplicate stations
stations <- stations[!is.na(stations[,3]),]
stations <- stations[!duplicated(stations[,1:2]),]
#change format to data.frame for further processing
stations <- as.data.frame(stations)

#change data type to spatial data with correct projection
coordinates(stations) <- ~lon + lat
proj4string(stations) <- CRS('+proj=longlat +datum=WGS84 +no_defs')

#simple visualization of the orginal dataset
bubble(stations, 'pm10', fill=F, do.sqrt=F, col = 'red', main="October data", maxsize = 2)

#normality check
qqPlot(stations$pm10)
shapiro.test(stations$pm10)
#normality check for log data
qqPlot(log(stations$pm10))
shapiro.test(log(stations$pm10))
#results show that is more convinient to use logarithm of pm10 concentrations





################################################
### empirical variogram calculation
### theoretical fuzzy variogram estimation
################################################

#calculation of emprical variogram, printing values and plotting the variogram
pm10.vgm <- variogram(log(pm10) ~ 1, stations)
print(pm10.vgm)
plot(pm10.vgm)

#estimation of fuzzy variogram - three components
#low estimate has highest value of range, and lowest values of nugget and sill
pm10.vgm.low <- vgm(psill = 0.07, model="Sph",range= 130, nugget=0.02)
plot(pm10.vgm, pm10.vgm.low, ylim = c(0,0.15))
#modal estimate has modal value of range, nugget and sill (can be also fitted automaticaly)
pm10.vgm.mid <- vgm(psill = 0.075, model="Sph",range= 100, nugget=0.025)
plot(pm10.vgm, pm10.vgm.mid, ylim = c(0,0.15))
#high estimate has lowest value of range, and highest values of nugget and sill
pm10.vgm.high <- vgm(psill = 0.085, model="Sph",range= 50, nugget=0.03)
plot(pm10.vgm, pm10.vgm.high, ylim = c(0,0.15))

#plot the variogram in one image as shown in the paper
plotVariograms(pm10.vgm,pm10.vgm.low,pm10.vgm.mid,pm10.vgm.high)





############################################################
#### preparation of kriging grid
############################################################

#load SHP file with Czech rep. borders
#the shapefile can be obtained as open data from https://osm.wno-edv-service.de/boundaries/
state <- readOGR("data/border.shp", "border")
state <- spTransform(state, CRS('+proj=longlat +datum=WGS84 +no_defs'))
state@data <- state@data[1]

#grid for kriging predictions (in the shape of Czech rep.)
bb <- bbox(state)

x <- seq(bb[1,1],bb[1,2],by = 0.0125)
y <- seq(bb[2,1],bb[2,2],by = 0.0125)
xy <- expand.grid(x,y)

grid <- SpatialPoints(xy)
grid <- SpatialPointsDataFrame(as.data.frame(grid), data = as.data.frame(rep(1,nrow(as.data.frame(grid)))))
proj4string(grid) <- CRS('+proj=longlat +datum=WGS84 +no_defs')

#the selection of prediction points in Czech rep. may take a while
ov <- over(grid, state)
grid@data$border <- ov$OBJECTID
grid <- na.exclude(as.data.frame(grid))

coordinates(grid) <- ~Var1 + Var2
gridded(grid) <- TRUE
proj4string(grid) <- CRS('+proj=longlat +datum=WGS84 +no_defs')

#remove the variables that are no longer needed
remove(xy,x,y,bb,state,ov)





################################################
### kriging
### fuzzy kriging
################################################

#calculate the modal value of kriging (logarithmic values)
pm10.kriged <- krige(log(pm10) ~ 1, stations, grid, model = pm10.vgm.mid)
#exponentiate the results back to get true predictions
pm10.kriged@data$pred.exp <- exp(pm10.kriged$var1.pred)

#remove data from the kriging that are not needed
pm10.kriged$var1.pred <- NULL
pm10.kriged$var1.var <- NULL

#plot predictions and locations of stations
spplot(pm10.kriged["pred.exp"], sp.layout=list("sp.points", stations, pch="+", cex = 2, col = "black"))

#obtain limits of sill, range and nugget for fuzzy predictions
sills = c(pm10.vgm.low$psill[2], pm10.vgm.high$psill[2])
ranges = c(pm10.vgm.low$range[2], pm10.vgm.high$range[2])
nuggets = c(pm10.vgm.low$psill[1], pm10.vgm.high$psill[1])

#estimation of minimum and maximum values of fuzzy kriging
#optimization scheme suggested by Loquin and Dubois (2010)
pm10.MinMax = calculateMinMaxPoints(sills, ranges, nuggets, c("Sph"), pm10.vgm.mid, 
                                    log(pm10) ~ 1, stations, grid, logResults=TRUE)

#create sequences of sill, range and nugget values of simulated annealing
sillsSeq = seq(pm10.vgm.low$psill[2], pm10.vgm.high$psill[2], length.out=15)
rangesSeq = seq(pm10.vgm.low$range[2], pm10.vgm.high$range[2], length.out=32)
nuggetsSeq = seq(pm10.vgm.low$psill[1], pm10.vgm.high$psill[1], length.out=10)

#estimation of minimum and maximum values of fuzzy kriging by simulated annealing 
#takes a very long time (4800 calculations of kriging) !!!
pm10.MinMaxSimulatedAnnealing = calculateMinMaxPoints(sillsSeq, rangesSeq, nuggetsSeq, c("Sph"),
                                                      pm10.vgm.mid, log(pm10) ~ 1, 
                                                      stations, grid, logResults=TRUE)

#visualization of differences amongst estimates by optimisation scheme and simulated annealing
hist(pm10.MinMax$min.data - pm10.MinMaxSimulatedAnnealing$min.data) 
hist(pm10.MinMax$max.data - pm10.MinMaxSimulatedAnnealing$max.data)


#attach values of prediction limits to result of kriging to spatial class
pm10.kriged[['modal']] <- pm10.kriged@data$pred.exp
pm10.kriged[['maximal']] <- pm10.MinMax$max.data
pm10.kriged[["minimal"]] <- pm10.MinMax$min.data
pm10.kriged[['maximalSA']] <- pm10.MinMaxSimulatedAnnealing$max.data
pm10.kriged[["minimalSA"]] <- pm10.MinMaxSimulatedAnnealing$min.data



#calculate the width of estimates for optimisation scheme and simulated annealing
pm10.kriged[["differenceSA"]] <- pm10.kriged[['maximalSA']] - pm10.kriged[["minimalSA"]]
pm10.kriged[["difference"]] <- pm10.kriged[['maximal']] - pm10.kriged[["minimal"]]

#visualize the range of estimates for optimisation scheme and simulated annealing
spplot(pm10.kriged, c("difference","differenceSA"),
       names.attr= c("difference according to optimisation scheme","difference according to simulated annealing"),
       colorkey=list(space="bottom"), layout=c(2,1), pretty=TRUE)
#based on this result it is concluded that simulated annealing provides more complete estimates

#plot limits of fuzzy surface estimated by simulated annealing
spplot(pm10.kriged, c("minimalSA","maximalSA"), names.attr= c("minimal value","maximal value"), colorkey=list(space="bottom"),
       layout=c(2,1), pretty=TRUE)





################################################
### determination of exceedance of threshold
################################################

#threshold value is determined based on literature review
threshold = 30

#determination of possibility and necessity of excedance
pm10.kriged[['possSA']] <- possibilityExceedance(pm10.kriged$modal, pm10.kriged$maximalSA, threshold)
pm10.kriged[['necSA']] <- necessityExceedance(pm10.kriged$minimalSA, pm10.kriged$modal, threshold)

#visualization of results
spplot(pm10.kriged, c("possSA","necSA"), names.attr= c("possibility of exceedance","necessity of exceedance"), colorkey=list(space="bottom"),
                layout=c(2,1), pretty=TRUE)