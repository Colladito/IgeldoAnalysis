rm(list = ls())
options(java.parameters = "-Xmx8000m")

# Ordenador Grande
workPath <- "D:/Victor/Master/TFM/R/2024_VictorCollado_UC3/"
setwd(paste0(workPath))
########################### DATA LOAD ####################################
library(MASS)
library(abind)
library(loadeR)
library(downscaleR)
library(loadeR.2nc)
library(climate4R.indices)
library(RCurl)
library(visualizeR)
library(ggplot2)
library(POT)
library(ismev)
# source("AuxiliarFunctions.R")

if (!file.exists(paste0(workPath, "prIgeldo.Rdata"))) {
  pr.Igeldo <- read.csv(paste0(workPath,"2023_TFM_Ana_Roberto/CorreoTFM/MATLAB/igeldo.csv"), header = TRUE, sep = ",")
  aux.dates <- strsplit(pr.Igeldo$X, " ")
  pr.day.dates <- lapply(1:length(pr.Igeldo$X), function(i){
    l <- aux.dates[[i]][1]
  })
  pr.day.dates <- unlist(pr.day.dates)
  pr.min.dates <- pr.day.dates
  pr.day.dates <- unique(pr.day.dates)
  pr.day.data <- lapply(1:length(pr.day.dates), function(i){
    ind <- which(pr.min.dates == pr.day.dates[i])
    l <- sum(pr.Igeldo$prec[ind], na.rm = TRUE)
  })
  pr.day.data <- unlist(pr.day.data)
  min <- list(dates = pr.Igeldo$X, data = pr.Igeldo$prec)
  day <- list(dates = pr.day.dates, data = pr.day.data)
  aux.dates <- strsplit(pr.Igeldo$X, " ")
  pr.hour.dates <- lapply(1:length(pr.Igeldo$X), function(i){
    a <- strsplit(aux.dates[[i]][2], ":")
    l <- paste0(aux.dates[[i]][1]," ",a[[1]][1],":00:00 GMT")
  })
  pr.hour.dates <- unlist(pr.hour.dates)
  pr.min.dates <- pr.hour.dates
  pr.hour.dates <- unique(pr.hour.dates)
  pr.hour.data <- lapply(1:length(pr.hour.dates), function(i){
    ind <- which(pr.min.dates == pr.hour.dates[i])
    l <- sum(pr.Igeldo$prec[ind], na.rm = TRUE)
  })
  pr.hour.data <- unlist(pr.hour.data)
  hour <- list(dates = pr.hour.dates, data = pr.hour.data)
  pr.Igeldo <- list(min = min, hour = hour, day = day)
  save(pr.Igeldo, file = paste0(workPath, "2023_TFM_Ana_Roberto/prIgeldo.Rdata"))
}else{
  load(paste0(workPath, "prIgeldo.Rdata"))
}
## Esto es con los datos descargados pero no afecta al script:
# era5dataset <- paste0(workPath,"ERA5/ECMWF_ERA5_SFC_clustering.ncml")
# di.era5 <- dataInventory(era5dataset)
if (!file.exists(paste0(workPath, "cinIgeldo.Rdata"))) {
  cin <- lapply(1940:2016, function(i){
    l <- lapply(1:12, function(m){
      if (m<10){
        era5dataset <- paste0(workPath,"ERA5/SFC_CIN.128/ECMWF_ERA5_", i,"0", m, "_SFC_CIN.128.nc")
      }else{
        era5dataset <- paste0(workPath,"ERA5/SFC_CIN.128/ECMWF_ERA5_", i, m, "_SFC_CIN.128.nc")
      }
      ll <- loadGridData(era5dataset, var = "cin", latLim = c(42.5, 43.5), lonLim = c(-2.5, 1.5), years = i, season = m)
    })
    l <- bindGrid(l, dimension = "time")
    l <- redim(l, drop = TRUE)
  })
  cin <- bindGrid(cin, dimension = "time")
  cin <- redim(cin, drop = TRUE)
  cin$Data[which(is.infinite(cin$Data))] <- NA
  save(cin, file = paste0(workPath, "2023_TFM_Ana_Roberto/cinIgeldo.Rdata"))
}else{
  load(paste0(workPath, "cinIgeldo.Rdata"))
}



if (!file.exists(paste0(workPath, "capeIgeldo.Rdata"))) {
  cape <- lapply(1940:2016, function(i){
    l <- lapply(1:12, function(m){
      if (m<10){
        era5dataset <- paste0(workPath,"ERA5/SFC_CAPE.128/ECMWF_ERA5_", i,"0", m, "_SFC_CAPE.128.nc")
      }else{
        era5dataset <- paste0(workPath,"ERA5/SFC_CAPE.128/ECMWF_ERA5_", i, m, "_SFC_CAPE.128.nc")
      }
      ll <- loadGridData(era5dataset, var = "cape", latLim = c(42.5, 43.5), lonLim = c(-2.5, 1.5), years = i, season = m)
    })
    l <- bindGrid(l, dimension = "time")
    l <- redim(l, drop = TRUE)
  })
  cape <- bindGrid(cape, dimension = "time")
  cape <- redim(cape, drop = TRUE)
  cape$Data[which(is.infinite(cape$Data))] <- NA
  save(cape, file = paste0(workPath, "2023_TFM_Ana_Roberto/capeIgeldo.Rdata"))
}else{
  load(paste0(workPath, "capeIgeldo.Rdata"))
}


if (!file.exists(paste0(workPath, "kIgeldo.Rdata"))) {
  K <- lapply(1940:2016, function(i){
    era5dataset <- paste0(workPath,"ERA5/SFC_K.128/ECMWF_ERA5_", i,"_SFC_K.128.nc")
    ll <- loadGridData(era5dataset, var = "kx", latLim = c(42.5, 43.5), lonLim = c(-2.5, 1.5), years = i)
  })
  K <- bindGrid(K, dimension = "time")
  K <- redim(K, drop = TRUE)
  K$Data[which(is.infinite(K$Data))] <- NA
  save(K, file = paste0(workPath, "2023_TFM_Ana_Roberto/kIgeldo.Rdata"))
}else{
  load(paste0(workPath, "kIgeldo.Rdata"))
}

if (!file.exists(paste0(workPath, "slpIgeldo.Rdata"))) {
  slp <- lapply(1940:2016, function(i){
    l <- lapply(1:12, function(m){
      if (m<10){
        era5dataset <- paste0(workPath,"2023_TFM_Ana_Roberto/ERA5/SFC_151.128/ECMWF_ERA5_", i,"0", m, "_SFC_151.128.nc")
      }else{
        era5dataset <- paste0(workPath,"2023_TFM_Ana_Roberto/ERA5/SFC_151.128/ECMWF_ERA5_", i, m, "_SFC_151.128.nc")
      }
      ll <- loadGridData(era5dataset, var = "msl", latLim = c(42.5, 43.5), lonLim = c(-2.5, 1.5), years = i, season = m)
    })
    l <- bindGrid(l, dimension = "time")
    l <- redim(l, drop = TRUE)
  })
  
  slp <- bindGrid(slp, dimension = "time")
  slp <- redim(slp, drop = TRUE)
  slp$Data[which(is.infinite(slp$Data))] <- NA
  save(slp, file = paste0(workPath, "2023_TFM_Ana_Roberto/slpIgeldo.Rdata"))
}else{
  load(paste0(workPath, "slpIgeldo.Rdata"))
}

if (!file.exists(paste0(workPath, "ttIgeldo.Rdata"))) {
  TT <- lapply(1940:2016, function(i){
    era5dataset <- paste0(workPath,"ERA5/SFC_TT.128/ECMWF_ERA5_", i,"_SFC_TT.128.nc")
    ll <- loadGridData(era5dataset, var = "totalx", latLim = c(42.5, 43.5), lonLim = c(-2.5, 1.5), years = i)
  })
  TT <- bindGrid(TT, dimension = "time")
  TT <- redim(TT, drop = TRUE)
  TT$Data[which(is.infinite(TT$Data))] <- NA
  save(TT, file = paste0(workPath, "2023_TFM_Ana_Roberto/ttIgeldo.Rdata"))
}else{
  load(paste0(workPath, "ttIgeldo.Rdata"))
}


## Homogenizamos las fechas:
fechas <- intersect(cin$Dates$start,cape$Dates$start)
I1 <- which(is.element(cin$Dates$start, fechas))
I2 <- which(is.element(cape$Dates$start, fechas))
cin$Data[] <- cin$Data[I1,,]
cin$Dates$start <- cin$Dates$start[I1]
cin$Dates$end <- cin$Dates$end[I1]

cape$Data <- cape$Data[I2,,]
attr(cape$Data, "dimensions") <- c("time","lat","lon")
cape$Dates$start <- cape$Dates$start[I2]
cape$Dates$end <- cape$Dates$end[I2]

Jk <- which(is.element(K$Dates$start, fechas))
Jt <- which(is.element(TT$Dates$start, fechas))
js <- which(is.element(slp$Dates$start, fechas))
K$Data[] <- K$Data[Jk,,]
attr(K$Data, "dimensions") <- c("time","lat","lon")
K$Dates$start <- K$Dates$start[Jk]
K$Dates$end <- K$Dates$end[Jk]

TT$Data <- TT$Data[Jt,,]
attr(TT$Data, "dimensions") <- c("time","lat","lon")
TT$Dates$start <- TT$Dates$start[Jt]
TT$Dates$end <- TT$Dates$end[Jt]

slp$Data <- slp$Data[Jt,,]
attr(slp$Data, "dimensions") <- c("time","lat","lon")
slp$Dates$start <- slp$Dates$start[Jt]
slp$Dates$end <- slp$Dates$end[Jt]


## Make the data daily

cin.day <- aggregateGrid(grid = cin, aggr.d = list(FUN = "mean" , na.rm = TRUE))
cin.day$Dates <- strftime(as.POSIXct(cin.day$Dates$start), format="%Y-%m-%d")

cape.day <- aggregateGrid(grid = cape, aggr.d = list(FUN = "mean" , na.rm = TRUE))
cape.day$Dates <- strftime(as.POSIXct(cape.day$Dates$start), format="%Y-%m-%d")

k.day <- aggregateGrid(grid = K, aggr.d = list(FUN = "mean" , na.rm = TRUE))
k.day$Dates <- strftime(as.POSIXct(k.day$Dates$start), format="%Y-%m-%d")

slp.day <- aggregateGrid(grid = slp, aggr.d = list(FUN = "mean" , na.rm = TRUE))
slp.day$Dates <- strftime(as.POSIXct(slp.day$Dates$start), format="%Y-%m-%d")

TT.day <- aggregateGrid(grid = TT, aggr.d = list(FUN = "mean" , na.rm = TRUE))
TT.day$Dates <- strftime(as.POSIXct(TT.day$Dates$start), format="%Y-%m-%d")

### Seleccion de los indices de precipitacion
fechas <- intersect(cin.day$Dates, cape.day$Dates)

### Indice x horas
## Indice año completo
timelapse <- "All-Year"
indDates <- which(is.element(pr.Igeldo$day$dates,fechas))
I2 <- which(is.element(fechas, pr.Igeldo$day$dates))

### Indices entre Octubre y Marzo
# timelapse <- "Aut-Wint"
# IJanMarch <- which(lubridate::month(as.POSIXct(pr.Igeldo$day$dates, tz = "GMT")) <= 3)
# IOctDec <- which(lubridate::month(as.POSIXct(pr.Igeldo$day$dates, tz = "GMT")) >= 10)
# IOctMarch <- union(IJanMarch, IOctDec)
# I2 <- which(is.element(fechas, pr.Igeldo$day$dates[IOctMarch]))
# indDates <- which(is.element(pr.Igeldo$day$dates, fechas[I2]))

### Indices entre Abril y Septiembre
# timelapse <- "Spr-Summ"
# IApr <- which(lubridate::month(as.POSIXct(pr.Igeldo$day$dates, tz = "GMT")) >= 4)
# ISep <- which(lubridate::month(as.POSIXct(pr.Igeldo$day$dates, tz = "GMT")) <= 9)
# IAprSept <- intersect(IApr, ISep)
# I2 <- which(is.element(fechas, pr.Igeldo$day$dates[IAprSept]))
# indDates <- which(is.element(pr.Igeldo$day$dates, fechas[I2]))



igeldo <- c(-2.009944 ,43.322886)
iCoord <- which.min(abs(cin$xyCoords$x - igeldo[1]))
jCoord <- which.min(abs(cin$xyCoords$y - igeldo[2]))

### Fijamos el threshold
threshold <- 5

#### Using POT method 503304
prec <- data.frame(obs = pr.Igeldo$day$data[indDates])
prec$time <- as.POSIXct(pr.Igeldo$day$dates[indDates], format = "%Y-%m-%d")
# Dividimos entre 86400 y 365.25 para que cada año corresponda con la unidad 
prec$time <- as.numeric(prec$time)/86400/365.25 
POTclust <- POT::clust(prec, threshold, tim.cond = 1/365.25, clust.max = T)
indPOT <- POTclust[,3]

## Añadimos delay si hace falta
delayslp <- 0
delay <- 0

nanIndex <- intersect(indPOT, which(is.na(cin.day$Data[I2, jCoord,iCoord])) + delay)
## nanIndex <- indPOT + delay
nanCAPE <- nanIndex - delay


#### Normal threshold sin POT
# indG90 <- which(pr.Igeldo$hour$data[indDates] >= threshold)
# nanIndex <- intersect(indG90, which(is.na(cin$Data[I2, jCoord,iCoord])))

#====================================================#

################ Clustering and Plots ################ 

#====================================================#

#### Function to pass a grid into a data.frame
# grid: data in grid which we want to pass to data fram
# timeind: time indices to set the name of rows
# delay: delay of de grid to the precipitation
gtdf <- function(grid, timeind, delay = 0){
  longitud <- dim(grid$Data)[2]
  latitud <- dim(grid$Data)[3]
  # Añadimos las fechas como indices del dataframe
  grid_df <- data.frame(row.names = grid$Dates$start[timeind + delay])
  
  for (long in 1:longitud) {
    for (lat in 1:latitud) {
      label <- paste("C", long, "x", lat, sep = "")
      grid_df[label] <- grid$Data[timeind, long, lat]
    }
  }
  return(grid_df)
}

gtdf2 <- function(grid, grid1, grid2, grid3, timeind, delay = 0, delayslp = 0){
  longitud <- dim(grid$Data)[2]
  latitud <- dim(grid$Data)[3]
  # Añadimos las fechas como indices del dataframe
  grid_df <- data.frame(row.names = grid$Dates[timeind + delay])
  # grid <- NULL
  # grid1 <- NULL
  # grid2 <- NULL
  # grid3 <- NULL
  for (long in c(-1,0,1)) {
    for (lat in c(-1,0,1)) {
      # label <- paste("CAPE", iCoord + long, "x", jCoord + lat, sep = "")
      # grid_df[label] <- grid$Data[timeind, jCoord+lat, iCoord+long]
      label <- paste("K", iCoord + long, "x", jCoord + lat, sep = "")
      grid_df[label] <- grid1$Data[timeind, jCoord+lat, iCoord+long]
      label <- paste("TT", iCoord + long, "x", jCoord + lat, sep = "")
      grid_df[label] <- grid2$Data[timeind, jCoord+lat, iCoord+long]
    }
  }
  # label <- "CAPE"
  # grid_df[label] <- scale(grid$Data[timeind, jCoord, iCoord], center = TRUE, scale = FALSE)
  # label <- "K"
  # grid_df[label] <- scale(grid1$Data[timeind, jCoord, iCoord], center = TRUE, scale = FALSE)
  # label <- "TT"
  # grid_df[label] <- scale(grid2$Data[timeind, jCoord, iCoord], center = TRUE, scale = FALSE)
  
  for (long in 1:longitud) {
    for (lat in 1:latitud) {
      label <- paste("CAPE", long, "x", lat, sep = "")
      grid_df[label] <- grid$Data[timeind, long, lat]
      label <- paste("SLP", long, "x", lat, sep = "")
      grid_df[label] <- grid3$Data[timeind-delayslp, long, lat]
    }
  }
  return(grid_df)
}

# # Transformamos el grid de CAPE a dataframe
# # CAPE with Delay
# cape.df <- gtdf(cape, I2[nanCAPE], delay)
# # K with Delay
# cape.df <- gtdf(K, I2[nanCAPE], delay)
# TT with Delay
# cape.df <- gtdf(TT, I2[nanCAPE], delay)

# Usando todo el grid
cape.df <- gtdf2(cape.day, k.day, TT.day, slp.day, I2[nanCAPE], delay = delay, delayslp = delayslp)

# Dataframe completo con el CAPE, la precipitacion y el tiempo asociado a ella.
compl.df <- cape.df
compl.df$prec <- pr.Igeldo$day$data[indDates[nanIndex]]
# compl.df$time <- prec$time[indDates[nanIndex]]
compl.df$time <- nanIndex
# Para el tiempo considero los indices ya que cada unidad corresponde a una hora

#====================================================#

##################### Clustering ##################### 

#====================================================#

# Creamos el vector de clusters
km <- kmeans(cape.df, centers = 2, iter.max = 10000)
clust <- data.frame(row.names = row.names(cape.df),
                    cluster = km$cluster)
summary(as.factor(clust$cluster))


#====================================================#

################## Aproximacion GPD ################## 

#====================================================#

# Usamos una vector de thresholds donde analizamos los parametros de forma para 
# cada valor del threshold

# Secuencia de tresholds que consideramos
secthres <- seq(threshold, 40, 0.1)

# Creamos data frames
xis <- data.frame(Threshold = secthres,
                  row.names = paste0("Threshold", secthres))
alphas <- data.frame(Threshold = secthres,
                     row.names = paste0("Threshold", secthres))

for (k in 1:2) {
  label <- paste0("Cluster", k)
  xis[label] <- rep(NA, length(secthres))
  alphas[label] <- rep(NA, length(secthres))
  label <- paste0("LowCluster", k)
  xis[label] <- rep(NA, length(secthres))
  alphas[label] <- rep(NA, length(secthres))
  label <- paste0("UpperCluster", k)
  xis[label] <- rep(NA, length(secthres))
  alphas[label] <- rep(NA, length(secthres))
}

### Aproximamos la GPD no estacionaria
# Si no trabajamos con el periodo de Octubre-Marzo habria que sustituir npy por 24*365.25
for (thres in seq_along(secthres)){
  for (k in 1:2) {
    # Non-Stationary
    out <- gpd.fit(compl.df$prec[which(clust$cluster == k)], threshold = secthres[thres], npy = 182, 
                   ydat = as.matrix(compl.df$time[which(clust$cluster == k)]), 
                   sigl = 1, siglink = identity, show = F)
    xis[thres,k^2+1] <- out$mle[3]
    alphas[thres,k^2+1] <- out$mle[2]
    xis[thres,k^2+2] <- out$mle[3] - 1.96*out$se[3]
    xis[thres,k^2+3] <- out$mle[3] + 1.96*out$se[3]
    alphas[thres,k^2+2] <- out$mle[2] - 1.96*out$se[2]
    alphas[thres,k^2+3] <- out$mle[2] + 1.96*out$se[2]
  }
}


### Plot del parámetro de forma
nonstatshapeplt <- ggplot(xis, aes(x = Threshold)) +
  geom_line(aes(y = Cluster1, color = "Cluster 1"), linewidth = 0.75) +
  geom_line(aes(y = Cluster2, color = "Cluster 2"), linewidth = 0.75) +
  # Confidence interval
  geom_ribbon(aes(ymin = LowCluster1, ymax = UpperCluster1, fill = "Cluster 1"), alpha = 0.2)+
  geom_ribbon(aes(ymin = LowCluster2, ymax = UpperCluster2, fill = "Cluster 2"), alpha = 0.2)+
  # geom_line(aes(y = LowCluster1, color = "CI 1"), linewidth = 0.75)+
  # geom_line(aes(y = UpperCluster1, color = "CI 1"), linewidth = 0.75)+
  # geom_line(aes(y = LowCluster2, color = "CI 2"), linewidth = 0.75)+
  # geom_line(aes(y = UpperCluster2, color = "CI 2"), linewidth = 0.75)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_color_manual(values = c("Cluster 1" = "#517CD3", 
                                "Cluster 2" = "#C95F9E")) +
  scale_fill_manual(values = c("Cluster 1" = "#003399", 
                               "Cluster 2" = "#8B0A50")) +
  theme_minimal() +
  labs(color = "MLE", y = "Shape Parameter", x = "Threshold", fill = "Conf. Interval"
       # title = "Evolution of Shape Parameter Along Threshold"
  )+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(color = guide_legend(order = 1), fill = guide_legend(order = 2))
print(nonstatshapeplt)

label <- paste0("pr_ALL_non-stat_shape_delay", delay, "_threshold.png")
ggsave(label, plot = nonstatshapeplt, 
       path = paste0("D:/Victor/Master/TFM/Graficos/Final Plots/Daily/", timelapse, "/SLP Delay", delayslp, "/NonStat_shape"),
       dpi = 300, bg = "white", width = 7, height = 5, units = "in")

### Plot del parámetro no estacionario
nonstatparamplt <- ggplot(alphas, aes(x = Threshold)) +
  geom_line(aes(y = Cluster1, color = "Cluster 1"), linewidth = 0.75) +
  geom_line(aes(y = Cluster2, color = "Cluster 2"), linewidth = 0.75) +
  # Confidence interval
  geom_ribbon(aes(ymin = LowCluster1, ymax = UpperCluster1, fill = "Cluster 1"), alpha = 0.2)+
  geom_ribbon(aes(ymin = LowCluster2, ymax = UpperCluster2, fill = "Cluster 2"), alpha = 0.2)+
  # geom_line(aes(y = LowCluster1, color = "CI 1"), linewidth = 1, alpha = 1)+
  # geom_line(aes(y = UpperCluster1, color = "CI 1"), linewidth = 1, alpha = 1)+
  # geom_line(aes(y = LowCluster2, color = "CI 2"), linewidth = 1, alpha = 1)+
  # geom_line(aes(y = UpperCluster2, color = "CI 2"), linewidth = 1, alpha = 1)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_color_manual(values = c("Cluster 1" = "#517CD3", 
                                "Cluster 2" = "#C95F9E")) +
  scale_fill_manual(values = c("Cluster 1" = "#003399", 
                               "Cluster 2" = "#8B0A50")) +
  theme_minimal() +
  labs(color = "MLE", y = "Non-Stationary Parameter", x = "Threshold", fill = "Conf. Interval"
       # title = "Evolution of Non-Stationary Parameter Along Threshold"
  )+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(color = guide_legend(order = 1), fill = guide_legend(order = 2))
print(nonstatparamplt)

label <- paste0("pr_ALL_non-stat_param_delay", delay, "_threshold.png")
ggsave(label, plot = nonstatparamplt, 
       path = paste0("D:/Victor/Master/TFM/Graficos/Final Plots/Daily/", timelapse, "/SLP Delay", delayslp, "/NonStat_param"), 
       dpi = 300, bg = "white", width = 7, height = 5, units = "in")

###### Stationary Model

xis <- data.frame(Threshold = secthres,
                  row.names = paste0("Threshold", secthres))

for (k in 1:2) {
  label <- paste0("Cluster", k)
  xis[label] <- rep(NA, length(secthres))
  label <- paste0("LowCluster", k)
  xis[label] <- rep(NA, length(secthres))
  label <- paste0("UpperCluster", k)
  xis[label] <- rep(NA, length(secthres))
}


### Aproximamos la GPD estacionaria
# Si no trabajamos con el periodo de Octubre-Marzo habria que sustituir npy por 24*365.25
for (thres in seq_along(secthres)){
  for (k in 1:2) {
    # Non-Stationary
    out <- gpd.fit(compl.df$prec[which(clust$cluster == k)], threshold = secthres[thres], npy = 365.25, show = F)
    xis[thres,k^2+1] <- out$mle[2]
    xis[thres,k^2+2] <- out$mle[2] - 1.96*out$se[2]
    xis[thres,k^2+3] <- out$mle[2] + 1.96*out$se[2]
  }
}


### Plot del parametro de forma
statshapeplt <- ggplot(xis, aes(x = Threshold)) +
  geom_line(aes(y = Cluster1, color = "Cluster 1"), linewidth = 0.75) +
  geom_line(aes(y = Cluster2, color = "Cluster 2"), linewidth = 0.75) +
  # Confidence interval
  geom_ribbon(aes(ymin = LowCluster1, ymax = UpperCluster1, fill = "Cluster 1"), alpha = 0.2)+
  geom_ribbon(aes(ymin = LowCluster2, ymax = UpperCluster2, fill = "Cluster 2"), alpha = 0.2)+
  # geom_line(aes(y = LowCluster1, color = "CI 1"), linewidth = 0.75)+
  # geom_line(aes(y = UpperCluster1, color = "CI 1"), linewidth = 0.75)+
  # geom_line(aes(y = LowCluster2, color = "CI 2"), linewidth = 0.75)+
  # geom_line(aes(y = UpperCluster2, color = "CI 2"), linewidth = 0.75)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_color_manual(values = c("Cluster 1" = "#517CD3", 
                                "Cluster 2" = "#C95F9E")) +
  scale_fill_manual(values = c("Cluster 1" = "#003399", 
                               "Cluster 2" = "#8B0A50")) +
  theme_minimal() +
  labs(color = "MLE", y = "Shape Parameter", x = "Threshold", fill = "Conf. Interval"
       # title = "Evolution of Shape Parameter Along Threshold"
  )+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(color = guide_legend(order = 1), fill = guide_legend(order = 2))
print(statshapeplt)

label <- paste0("pr_ALL_stat_shape_delay", delay, "_threshold.png")
ggsave(label, plot = statshapeplt, 
       path = paste0("D:/Victor/Master/TFM/Graficos/Final Plots/Daily/", timelapse, "/SLP Delay", delayslp, "/Stat_shape"), 
       dpi = 300, bg = "white", width = 7, height = 5, units = "in")



#### Boxplots de los peaks en cada cluster

compl.df$cluster <- as.factor(clust$cluster)

uppind <- which(compl.df$prec >= threshold)
box_all_year <- ggplot(compl.df[uppind,], aes(y = prec, x = cluster, group = cluster, fill = cluster))+
  geom_boxplot()+
  scale_fill_manual(values = c("#517CD3", "#C95F9E"))+
  scale_y_continuous(trans='log2')+
  theme_minimal()+
  labs(fill = "Cluster", y = "Precipitation Over Threshold (mm)", 
       # title = "Boxplots of POT in Each Cluster", 
       x = "Cluster", color = "Station")+
  theme(plot.title = element_text(hjust = 0.5, size = 15))
print(box_all_year)

label <- paste0("pr_ALL_stat_shape_delay", delay, "_", timelapse, "_boxplot_logscale.png")
ggsave(label, plot = box_all_year, 
       path = paste0("D:/Victor/Master/TFM/Graficos/Final Plots/Daily/", timelapse, "/SLP Delay", delayslp), 
       dpi = 300, bg = "white", width = 7, height = 5, units = "in")


# Diferencia cluster entre estaciones

statind <- union(which(lubridate::month(as.POSIXct(row.names(compl.df), tz = "GMT")) <= 3),
             which(lubridate::month(as.POSIXct(row.names(compl.df), tz = "GMT")) >= 10))
compl.df$station <- "Spr-Summ"
compl.df$station[statind] <- "Aut-Wint"
compl.df$station <- as.factor(compl.df$station)


uppind <- which(compl.df$prec >= threshold)
box_stations <- ggplot(compl.df[uppind,], aes(y = prec, x = cluster, group = interaction(station, cluster), fill = interaction(station, cluster)))+
  geom_boxplot()+
  scale_fill_manual(values = c("#517CD3", "#003399", "#C95F9E", "#8B0A50"))+
  scale_y_continuous(trans='log2')+
  theme_minimal()+
  labs(fill = "Station/Cluster", y = "Precipitation Over Threshold (mm)",
       # title = "Boxplots of POT in Each Cluster",
       x = "Cluster")+
  theme(plot.title = element_text(hjust = 0.5, size = 15))
print(box_stations)

label <- paste0("pr_ALL_stat_shape_delay", delay, "_stations_boxplot_logscale.png")
ggsave(label, plot = box_stations,
       path = paste0("D:/Victor/Master/TFM/Graficos/Final Plots/Daily/", timelapse, "/SLP Delay", delayslp),
       dpi = 300, bg = "white", width = 7, height = 5, units = "in")


##### Fitted GPD-CDF of each clusters

#====================================================================#

################## Aproximacion GPD on each cluster ################## 

#====================================================================#

# Auxiliar functions
aggregatedGPD_CDF <- function(x, th1, sc1, sh1, th2, sc2, sh2) {
  pgpd(x, loc = th1, scale = sc1, shape = sh1) * pgpd(x, loc = th2, scale = sc2, shape = sh2)
}

# Invierte la CDF de las dos GPD, i.e., da los quantiles dadas unas probabilades ordenadas
aggregatedGPD_INV <- function(probs, th1, sc1, sh1, th2, sc2, sh2) {
  # Inicializamos el vector con los cuantiles
  x <- numeric(length(probs))
  
  for (i in 1:length(x)) {
    # Tomamos el punto medio del quantil para cada probabilidad entre las dos GEV 
    ini <- (qgpd(probs[i], loc = th1, scale = sc1, shape = sh1) + qgpd(probs[i], loc = th2, scale = sc2, shape = sh2)) / 2
    
    # Obtenemos la raiz (el 0) que iguale la probabilidad de las GEV juntas 
    # con la probabilidad definida que tenemos
    result <- uniroot(function(y){
      probs[i] - aggregatedGPD_CDF(y, th1, sc1, sh1, th2, sc2, sh2)},
      interval = c(ini,10^10))
    
    x[i] <- result$root
  }
  
  return(x)
}





# We use the threshold as 5 for each cluster because the cluster 2 is negative.
compl.gpd <- gpd.fit(compl.df$prec, threshold = 5, npy = 365.25, show = T)
gpd1 <- gpd.fit(compl.df$prec[which(clust$cluster == 1)], threshold = 5, npy = 365.25, show = T)
gpd2 <- gpd.fit(compl.df$prec[which(clust$cluster == 2)], threshold = 5, npy = 365.25, show = T)

gpd.diag(compl.gpd)
gpd.diag(gpd1)
gpd.diag(gpd2)



df1 <- data.frame(prec = compl.df$prec[which(clust$cluster == 1)], 
                  row.names = row.names(compl.df[which(clust$cluster == 1),]))
df2 <- data.frame(prec = compl.df$prec[which(clust$cluster == 2)], 
                  row.names = row.names(compl.df[which(clust$cluster == 2),]))

gpd.cdf<- function(param, thres, data){
  1 - (1 + (param[2] * (data - thres))/param[1])^(-1/param[2])
}

gpd.q<- function(param, thres, probs){
  thres + (param[1] * (probs^( - param[2]) - 1))/param[2]
}



######### PP-plot ######### 
ppgpd <-  function(param, thres, data){
  # Create the values
  x_values <- (1:length(data))/length(data)
  y_values <- gpd.cdf(param, thres, sort(data))
  df <- data.frame(Empirical = x_values, Model = y_values)
  
  ggplot(df, aes(x = Empirical, y = Model)) +
    geom_point(color = "#517CD3") +
    geom_abline(slope = 1, intercept = 0, color = "#C95F9E", linewidth = 0.75) +
    labs(x = "Empirical", y = "Model", title = "PP-Plot")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5))
}

ppgpd(gpd1$mle, 5, compl.df$prec[which(compl.df$cluster == 1)])

label <- "ppplot_clust1.png"
ggsave(label, 
       path = paste0("D:/Victor/Master/TFM/Graficos/Final Plots/Daily/", timelapse, "/SLP Delay", delayslp, "/Diagnostics"), 
       dpi = 300, bg = "white", width = 7, height = 5, units = "in")


ppgpd(gpd2$mle, 5, compl.df$prec[which(compl.df$cluster == 2)])
label <- "ppplot_clust2.png"
ggsave(label, 
       path = paste0("D:/Victor/Master/TFM/Graficos/Final Plots/Daily/", timelapse, "/SLP Delay", delayslp, "/Diagnostics"), 
       dpi = 300, bg = "white", width = 7, height = 5, units = "in")


######### QQ-plot ######### 
qqgpd <-  function(param, thres, data){
  # Create the values
  x_values <- gpd.q(param, thres, 1 - (1:length(data)/(length(data) + 1)))
  y_values <- sort(data)
  df <- data.frame(Empirical = x_values, Theorical = y_values)
  
  ggplot(df, aes(x = Empirical, y = Theorical)) +
    geom_point(color = "#517CD3") +
    geom_abline(slope = 1, intercept = 0, color = "#C95F9E", size = 0.75) +
    labs(x = "Empirical Quantiles", y = "Theorical Quantiles", title = "QQ-Plot")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5))
}

qqgpd(gpd1$mle, 5, compl.df$prec[which(compl.df$cluster == 1)])
label <- "qqplot_clust1.png"
ggsave(label, 
       path = paste0("D:/Victor/Master/TFM/Graficos/Final Plots/Daily/", timelapse, "/SLP Delay", delayslp, "/Diagnostics"), 
       dpi = 300, bg = "white", width = 7, height = 5, units = "in")

qqgpd(gpd2$mle, 5, compl.df$prec[which(compl.df$cluster == 2)])
label <- "qqplot_clust2.png"
ggsave(label, 
       path = paste0("D:/Victor/Master/TFM/Graficos/Final Plots/Daily/", timelapse, "/SLP Delay", delayslp, "/Diagnostics"), 
       dpi = 300, bg = "white", width = 7, height = 5, units = "in")



######### CDF ######### 

Prob2 <- seq(1, 10000, length = 1000) / (10000 + 1)
cdf1 <- qgpd(Prob2, loc = 5, scale = gpd1$mle[1], shape = gpd1$mle[2])
cdf2 <- qgpd(Prob2, loc = 5, scale = gpd2$mle[1], shape = gpd2$mle[2])
datacdf <- data.frame(prob = Prob2,
                      compl.gpd = aggregatedGPD_INV(Prob2, 5, gpd1$mle[1], gpd1$mle[2], 5, gpd2$mle[1], gpd2$mle[2]),
                      clust1 = cdf1,
                      clust2 = cdf2)
com.cdf <- data.frame(prec = sort(compl.df$prec), 
                      cdf = aggregatedGPD_CDF(sort(compl.df$prec), 5, gpd1$mle[1], gpd1$mle[2], 5, gpd2$mle[1], gpd2$mle[2]))

dfclust1 <- data.frame(prec = sort(compl.df$prec[which(compl.df$cluster == 1)]),
                       cdf = pgpd(sort(compl.df$prec[which(compl.df$cluster == 1)]), 
                                  loc = 5, scale = gpd1$mle[1], shape = gpd1$mle[2]))
dfclust2 <- data.frame(prec = sort(compl.df$prec[which(compl.df$cluster == 2)]),
                       cdf = pgpd(sort(compl.df$prec[which(compl.df$cluster == 2)]), 
                                  loc = 5, scale = gpd2$mle[1], shape = gpd2$mle[2]))

ggplot(datacdf, aes(y = prob)) +
  # Ploteamos las distribuciones aproximadas 
  geom_line(aes(x = compl.gpd, color = "Aggregated"), size = 0.75)+
  geom_line(aes(x = clust1, color = "Cluster 1"), size = 0.75)+
  geom_line(aes(x = clust2, color = "Cluster 2"), size = 0.75)+
  geom_point(data = com.cdf, aes(x = prec, y = cdf), alpha = 0.3, shape = 21, color = "#F08F46")+
  geom_point(data = dfclust1, aes(x = prec, y = cdf), alpha = 0.3, shape = 21, color = "#517CD3", fill = "#517CD3")+
  geom_point(data = dfclust2, aes(x = prec, y = cdf), alpha = 0.3, shape = 21, color = "#C95F9E", fill = "#C95F9E")+
  scale_color_manual(values = c("Cluster 1" = "#517CD3", "Cluster 2" = "#C95F9E", "Aggregated" = "#F08F46")) +
  labs(x = "Precipitation (mm)", y = "Probability", color = "")+
  theme_minimal()

label <- "cdf.png"
ggsave(label, 
       path = paste0("D:/Victor/Master/TFM/Graficos/Final Plots/Daily/", timelapse, "/SLP Delay", delayslp, "/Diagnostics"), 
       dpi = 300, bg = "white", width = 7, height = 5, units = "in")






######### Return periods ######### 

Prob2 <- seq(1, 1000000, length = 10000) / (1000000 + 1)
T2 <- 1 / (1 - Prob2)


library(scales)

Rfit1 <- qgpd(Prob2, loc = 5, scale = gpd1$mle[1], shape = gpd1$mle[2])
Rfit2 <- qgpd(Prob2, loc = 5, scale = gpd2$mle[1], shape = gpd2$mle[2])
Rfit3 <- aggregatedGPD_INV(Prob2, 5, gpd1$mle[1], gpd1$mle[2], 5, gpd2$mle[1], gpd2$mle[2])

data <- data.frame(years = T2,
                   Clust1 = Rfit1,
                   Clust2 = Rfit2,
                   RealCombWei = Rfit3)

Rsort <- sort(compl.df$prec[which(compl.df$cluster == 1)])
t1 <- 1/(1-pgpd(Rsort, loc = 5, scale = gpd1$mle[1], shape = gpd1$mle[2]))

Rsort2 <- sort(compl.df$prec[which(compl.df$cluster == 2)])
t2 <- 1/(1-pgpd(Rsort2, loc = 5, scale = gpd2$mle[1], shape = gpd2$mle[2]))

Rsort3 <- sort(compl.df$prec)
t3 <- 1/(1-pgpd(Rsort3, loc = 5, scale = compl.gpd$mle[1], shape = compl.gpd$mle[2]))

# Data frames
data1 <- data.frame(years = t1,
                    R1 = Rsort)
data2 <- data.frame(years = t2,
                    R2 = Rsort2)
asint2 <- 5-gpd2$mle[1]/gpd2$mle[2]

data3 <- data.frame(years = t3,
                    R3 = Rsort3)



ggplot(data, aes(x = years)) +
  # Ploteamos las distribuciones aproximadas 
  geom_line(aes(y = Clust1, color = "Cluster 1"), size = 0.75) +
  geom_line(aes(y = Clust2, color = "Cluster 2"), size = 0.75) +
  geom_line(aes(y = RealCombWei, color = "Aggregated"), size = 0.75) +

  # Ploteamos la simulacion con 100 valores para ver los puntos 
  geom_point(data = data1, aes(x = years, y = R1), color = "#517CD3",
             alpha = 0.5, size = 2)+
  geom_point(data = data2, aes(x = years, y = R2), color = "#C95F9E",
             alpha = 0.5, size = 2)+
  # geom_point(data = data3, aes(x = years, y = R3), color = "brown2",
  #            size = 2, shape = 21, alpha = 0.5)+
  # Asintota para el Cluster 2
  geom_hline(yintercept = asint2, color = "#C95F9E")+
  theme_minimal() +
  annotate(geom="text", x=10, y=asint2+6, label=paste("Asymptotic limit:", round(asint2)),
           color="#C95F9E", size=5)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = c("Cluster 1" = "#517CD3", "Cluster 2" = "#C95F9E",
                                "Aggregated" = "#F08F46")) +
  scale_linetype_manual(values = c("Real" = "solid", "Fitted" = "dashed")) +
  labs(color = "", linetype = "Type") +
  labs(x = "Years", y = "Precipitation Over Threshold (mm)")

label <- "return_period.png"
ggsave(label, 
       path = paste0("D:/Victor/Master/TFM/Graficos/Final Plots/Daily/", timelapse, "/SLP Delay", delayslp, "/Diagnostics"), 
       dpi = 300, bg = "white", width = 7, height = 5, units = "in")


#======================================================================#

################## PRUEBAS ADICIONALES CON MAS CLUSTERS ################

#======================================================================#


### Make another cluster in the frechet type tail

indClust1 <- which(compl.df$cluster == 1)
km2 <- kmeans(cape.df[indClust1,], centers = 2)
summary(as.factor(km2$cluster))
indClust1_1 <- which(km2$cluster == 1)
indClust1_2 <- which(km2$cluster == 2)


# Fit stationary model
secthres <- seq(threshold, 35, 0.1)
xis2 <- data.frame(Threshold = secthres,
                  row.names = paste0("Threshold", secthres))

for (k in 1:2) {
  label <- paste0("Cluster", k)
  xis2[label] <- rep(NA, length(secthres))
  label <- paste0("LowCluster", k)
  xis2[label] <- rep(NA, length(secthres))
  label <- paste0("UpperCluster", k)
  xis2[label] <- rep(NA, length(secthres))
}


### Aproximamos la GPD estacionaria
for (thres in seq_along(secthres)){
  for (k in 1:2) {
    # Non-Stationary
    out <- gpd.fit(compl.df$prec[indClust1[which(km2$cluster == k)]], threshold = secthres[thres], npy = 365.25, show = F)
    xis2[thres,k^2+1] <- out$mle[2]
    xis2[thres,k^2+2] <- out$mle[2] - 1.96*out$se[2]
    xis2[thres,k^2+3] <- out$mle[2] + 1.96*out$se[2]
  }
}


### Plot del parametro de forma
statshapeplt2 <- ggplot(xis2, aes(x = Threshold)) +
  geom_line(aes(y = Cluster1, color = "Cluster 1"), linewidth = 0.75) +
  geom_line(aes(y = Cluster2, color = "Cluster 2"), linewidth = 0.75) +
  # Confidence interval
  geom_ribbon(aes(ymin = LowCluster1, ymax = UpperCluster1, fill = "Cluster 1"), alpha = 0.2)+
  geom_ribbon(aes(ymin = LowCluster2, ymax = UpperCluster2, fill = "Cluster 2"), alpha = 0.2)+
  # geom_line(aes(y = LowCluster1, color = "CI 1"), linewidth = 0.75)+
  # geom_line(aes(y = UpperCluster1, color = "CI 1"), linewidth = 0.75)+
  # geom_line(aes(y = LowCluster2, color = "CI 2"), linewidth = 0.75)+
  # geom_line(aes(y = UpperCluster2, color = "CI 2"), linewidth = 0.75)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_color_manual(values = c("Cluster 1" = "#517CD3", 
                                "Cluster 2" = "#C95F9E")) +
  scale_fill_manual(values = c("Cluster 1" = "#003399", 
                               "Cluster 2" = "#8B0A50")) +
  theme_minimal() +
  labs(color = "MLE", y = "Shape Parameter", x = "Threshold", fill = "Conf. Interval"
       # title = "Evolution of Shape Parameter Along Threshold"
  )+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(color = guide_legend(order = 1), fill = guide_legend(order = 2))
print(statshapeplt2)

label <- paste0("pr_ALL_stat_shape_delay", delay, "_threshold_2.png")
ggsave(label, plot = statshapeplt2, 
       path = paste0("D:/Victor/Master/TFM/Graficos/Final Plots/Daily/", timelapse, "/SLP Delay", delayslp, "/Stat_shape"), 
       dpi = 300, bg = "white", width = 7, height = 5, units = "in")


# Cluster 1 negativo para el siguiente threshold
th2 <- 35.5

indClust1 <- which(compl.df$cluster == 1)
km3 <- kmeans(cape.df[indClust1[indClust1_2],], centers = 2)
summary(as.factor(km3$cluster))
indClust1_2_1 <- which(km3$cluster == 1)
indClust1_2_2 <- which(km3$cluster == 2)


# Fit stationary model
xis3 <- data.frame(Threshold = secthres,
                   row.names = paste0("Threshold", secthres))

for (k in 1:2) {
  label <- paste0("Cluster", k)
  xis3[label] <- rep(NA, length(secthres))
  label <- paste0("LowCluster", k)
  xis3[label] <- rep(NA, length(secthres))
  label <- paste0("UpperCluster", k)
  xis3[label] <- rep(NA, length(secthres))
}


### Aproximamos la GPD estacionaria
for (thres in seq_along(secthres)){
  for (k in 1:2) {
    # Non-Stationary
    out <- gpd.fit(compl.df$prec[indClust1[indClust1_2[which(km3$cluster == k)]]], 
                   threshold = secthres[thres], npy = 365.25, show = F)
    xis3[thres,k^2+1] <- out$mle[2]
    xis3[thres,k^2+2] <- out$mle[2] - 1.96*out$se[2]
    xis3[thres,k^2+3] <- out$mle[2] + 1.96*out$se[2]
  }
}


### Plot del parametro de forma
statshapeplt3 <- ggplot(xis3, aes(x = Threshold)) +
  geom_line(aes(y = Cluster1, color = "Cluster 1"), linewidth = 0.75) +
  geom_line(aes(y = Cluster2, color = "Cluster 2"), linewidth = 0.75) +
  # Confidence interval
  geom_ribbon(aes(ymin = LowCluster1, ymax = UpperCluster1, fill = "Cluster 1"), alpha = 0.2)+
  geom_ribbon(aes(ymin = LowCluster2, ymax = UpperCluster2, fill = "Cluster 2"), alpha = 0.2)+
  # geom_line(aes(y = LowCluster1, color = "CI 1"), linewidth = 0.75)+
  # geom_line(aes(y = UpperCluster1, color = "CI 1"), linewidth = 0.75)+
  # geom_line(aes(y = LowCluster2, color = "CI 2"), linewidth = 0.75)+
  # geom_line(aes(y = UpperCluster2, color = "CI 2"), linewidth = 0.75)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_color_manual(values = c("Cluster 1" = "#517CD3", 
                                "Cluster 2" = "#C95F9E")) +
  scale_fill_manual(values = c("Cluster 1" = "#003399", 
                               "Cluster 2" = "#8B0A50")) +
  theme_minimal() +
  labs(color = "MLE", y = "Shape Parameter", x = "Threshold", fill = "Conf. Interval"
       # title = "Evolution of Shape Parameter Along Threshold"
  )+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(color = guide_legend(order = 1), fill = guide_legend(order = 2))
print(statshapeplt3)

