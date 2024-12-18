library(tidyverse)
library(data.table)
library(leaflet)
library(sf)

bikeRentals <- fread("bikeshare/daily_rent_detail.csv")
stationList <- read_csv("bikeshare/station_list.csv")
usageFreq <- read_csv("bikeshare/usage_frequency.csv")
weatherFrame <- read_csv("bikeshare/weather.csv")
 
table(bikeRentals$rideable_type) # Classic, Docked, and electric bike
bikeRentals$classic <- factor(bikeRentals$rideable_type, 
                              levels = c("classic_bike", "electric_bike"))
## Do electir bikes have longer trips than classic bikes
bikeRentals <- bikeRentals %>% mutate(tripDuration = difftime(ended_at, started_at, units = "mins"))
wilcox.test(as.numeric(tripDuration) ~ classic, data = bikeRentals)
hist(bikeRentals[classic == "classic_bike", tripDuration] %>% as.numeric, breaks = 50)

## Exponential Model analysis
sumClassic <- sum(bikeRentals[classic == "classic_bike", tripDuration]) %>% as.numeric
lengthClassic <- length(bikeRentals[classic == "classic_bike", tripDuration])
sumElectric <- sum(bikeRentals[classic == "electric_bike", tripDuration]) %>% as.numeric
lengthElectric <- length(bikeRentals[classic == "electric_bike", tripDuration])
priorAlpha <- 1
priorBeta <- 1/2

classicLambdas <- rgamma(1e+5, lengthClassic + priorAlpha, sumClassic + priorBeta)
electricLambdas <- rgamma(1e+5, lengthElectric + priorAlpha, sumElectric + priorBeta)
classicSamps <- rexp(1e+5, classicLambdas)
electricSamps <- rexp(1e+5, electricLambdas)

## Testing Fit

testReps <- 1e+3
postYs <- matrix(rexp(1e+4 * testReps, rep(classicLambdas, each = testReps)), 
                 nrow = testReps)
tenthStat <- apply(postYs, 2, quantile, 0.1)
medianStat <- apply(postYs, 2, median)
hist(tenthStat) ; abline(v = quantile(bikeRentals[classic == "classic_bike", tripDuration], 0.1) %>% as.numeric())
hist(medianStat) ; abline(v = median(bikeRentals[classic == "classic_bike", tripDuration]) %>% as.numeric())





## Creating Map

## Change to only look at the summer of 2023
## Memorial Day to Labor Day
### Adding Coordinates to stationList
stationIndices <- match(stationList$station_id, as.numeric(bikeRentals$start_station_id))
stationCoords <- bikeRentals[stationIndices[-915:-914], .(start_lat, start_lng)] %>% drop_na %>% as.matrix
stationList <- data.frame(stationID = as.numeric(stationList$station_id[-915:-913]), stationName = stationList$station_name[-915:-913], stationCoords)
dailyFrequency <- bikeRentals[, .(pickup = .N), by = .(year(started_at), yday(started_at), start_station_id)]
## Adding Coords to usageFreq
coordMatch <- function() {
  
  latVec <- numeric(nrow(dailyFrequency))
  lonVec <- numeric(nrow(dailyFrequency))
  for (i in 1:nrow(dailyFrequency)) {
    
    nameInd <- which(dailyFrequency$start_station_id[i] == stationList$stationID)[1]
    latVec[i] <- stationList$start_lat[nameInd]
    lonVec[i] <- stationList$start_lng[nameInd]
    
  }
  dailyFrequency$latCoord <- latVec
  dailyFrequency$lonCoord <- lonVec
  return(dailyFrequency)
  
}

coordFrame <- coordMatch()
specInd <- sample(1:nrow(coordFrame), 1e+2)
geoR::variog(coords = coordFrame[specInd, 5:6], data = coordFrame$pickup[specInd],
             uvec = seq(0, 0.25 * max(iDist(coordFrame[specInd, 5:6])), length = 50))

data.frame(stationCoords, stationName = stationList$station_name[-915:-913]) %>% leaflet %>% 
  addTiles() %>% addMarkers(lng = ~start_lng, lat = ~start_lat, 
                            label = ~stationName, size = 0.1)

## Testing generic EXP
splitData <- tapply(bikeRentals$tripDuration, bikeRentals$classic, c)
nVec <- map_dbl(splitData, length)
splitMat <- do.call(cbind, splitData)
splitMat <- splitMat[1:1000, ]
nVec <- rep(1000, 2)
expTest <- jags.model("testingExp.txt", data = list(y = splitMat, 
                                                    n = nVec))
expSamps <- coda.samples(expTest, c("lambda", "newY"), 1e+4)

expJags <- function(formula, data, coordMat2) {
  
  # Gather response variable and design matrix
  responseVariable <- model.frame(formula, data = data)[[1]]
  regressionMatrix <- model.matrix(formula, data = data) 
  
  # Information Vars
  responseSize <- length(responseVariable)
  predictorSize <- ncol(regressionMatrix)
  
  # Compute Distances
  dMat <- spBayes::iDist(coordMat2)
  
  # JAGS model
  spatialJags <- jags.model(file = "~/Documents/Fall 2023/Spatial Stuffs/spatialJags.txt", 
                            data = list(D = dMat,
                                        n = responseSize,
                                        p = predictorSize,
                                        y = responseVariable,
                                        X = regressionMatrix,
                                        iMat = diag(nrow = nrow(dMat)) # Identity matrix of same dims as coordMat
                            ), 
                            n.chains = 1)
  
  return(spatialJags)
  
}

spaitalNB <- function(formula, data, coordMat2) {
  
  # Gather response variable and design matrix
  responseVariable <- model.frame(formula, data = data)[[1]]
  regressionMatrix <- model.matrix(formula, data = data) 
  
  # Information Vars
  responseSize <- length(responseVariable)
  predictorSize <- ncol(regressionMatrix)
  
  # Compute Distances
  dMat <- spBayes::iDist(coordMat2)
  
  # JAGS model
  spatialJags <- jags.model(file = "nbJags.txt", 
                            data = list(D = dMat,
                                        n = responseSize,
                                        p = predictorSize,
                                        y = responseVariable,
                                        X = regressionMatrix,
                                        iMat = diag(nrow = nrow(dMat)),
                                        zVec = rep(0, times = responseSize)# Identity matrix of same dims as coordMat
                            ), 
                            n.chains = 1)
  
  return(spatialJags)
  
}

coordFrame23 <- coordFrame[year == 2023]

## Test Duration Vec
summerBikes <- bikeRentals[year(started_at) == 2023 & yday(started_at) >= 149 & yday(started_at) <= 247][order(started_at)]
summerBikes$tripDuration <- difftime(summerBikes$ended_at, summerBikes$started_at, units = "min") %>% as.numeric()
summerBikes <- summerBikes[tripDuration >= 0]

weatherFrame <- as.data.table(weatherFrame)
summerWeather <- weatherFrame[year(datetime) == 2023 & yday(datetime) >= 149 & yday(datetime) <= 247]
weatherVars <- c(2, 3, 4, 5, 6, 8:14, 17:22, 25, 27:29)
summerWeather <- summerWeather[, ..weatherVars]
summerWeather$sunrise <- hour(summerWeather$sunrise)
summerWeather$sunset <- hour(summerWeather$sunset)

weatherGen <- function() {
  
  weatherMat <- matrix(0, nrow = nrow(summerBikes), 
                       ncol = length(weatherVars[-1]))
  for (i in 1:nrow(summerBikes)) {
    
    daySelect <- which(summerBikes$started_at[i] == yday(summerWeather$datetime))
    weatherMat[i, ] <- as.numeric(summerWeather[daySelect, -1])
    
  }
  names(weatherMat) <- names(summerWeather)[-1]
  return(weatherMat)
  
}
newVars <- weatherGen()

which(yday(summerBikes$started_at[1]) == yday(summerWeather$datetime))


plot(density(summerBikes$tripDuration), xlim = c(0, 300))
