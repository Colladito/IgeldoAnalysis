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

# cin.day <- aggregateGrid(grid = cin, aggr.d = list(FUN = "mean" , na.rm = TRUE))

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


##### Spatial plot for Sea Level Pressure & CAPE

# library(sp)
# test <- data.frame(x = -2.009944, y = 43.322886)
# coordinates(test) <- ~ x + y
# pts <- list("sp.points", test, pch=19, cex = 1, col="black")
# slpclim <- climatology(slp)
# # slpclim$Data <- slpclim$Data[1,,]
# 
# # Plot con la correlacion de Spearman en el grid cercano a Igeldo
# spatialPlot(slpclim, backdrop.theme = "coastline",
#             color.theme = "jet.colors", main = "", 
#             sp.layout = list(pts), xlim =c(-3,1), ylim = c(42,44))
# 
# 
# capeclim <- climatology(cape)
# 
# # Plot con la correlacion de Spearman en el grid cercano a Igeldo
# spatialPlot(capeclim, backdrop.theme = "coastline",
#             color.theme = "jet.colors", main = "", 
#             sp.layout = list(pts), xlim =c(-3,1), ylim = c(42,44))



### Seleccion de los indices de precipitacion

### Indice x horas
## Indice año completo
timelapse <- "All-Year"
indDates <- which(is.element(pr.Igeldo$hour$dates,fechas))
I2 <- which(is.element(fechas, pr.Igeldo$hour$dates))

### Indices entre Octubre y Marzo
# timelapse <- "Aut-Wint"
# IJanMarch <- which(lubridate::month(as.POSIXct(pr.Igeldo$hour$dates, tz = "GMT")) <= 3)
# IOctDec <- which(lubridate::month(as.POSIXct(pr.Igeldo$hour$dates, tz = "GMT")) >= 10)
# IOctMarch <- union(IJanMarch, IOctDec)
# I2 <- which(is.element(fechas, pr.Igeldo$hour$dates[IOctMarch]))
# indDates <- which(is.element(pr.Igeldo$hour$dates, fechas[I2]))

### Indices entre Abril y Septiembre
# timelapse <- "Spr-Summ"
# IApr <- which(lubridate::month(as.POSIXct(pr.Igeldo$hour$dates, tz = "GMT")) >= 4)
# ISep <- which(lubridate::month(as.POSIXct(pr.Igeldo$hour$dates, tz = "GMT")) <= 9)
# IAprSept <- intersect(IApr, ISep)
# I2 <- which(is.element(fechas, pr.Igeldo$hour$dates[IAprSept]))
# indDates <- which(is.element(pr.Igeldo$hour$dates, fechas[I2]))



# cape.day <- aggregateGrid(grid = cape, aggr.d = list(FUN = "mean" , na.rm = TRUE))

igeldo <- c(-2.009944 ,43.322886)
iCoord <- which.min(abs(cin$xyCoords$x - igeldo[1]))
jCoord <- which.min(abs(cin$xyCoords$y - igeldo[2]))

### Fijamos el threshold
threshold <- 5

#### Using POT method 503304
prec <- data.frame(obs = pr.Igeldo$hour$data[indDates])
prec$time <- as.POSIXct(pr.Igeldo$hour$dates[indDates], format = "%Y-%m-%d %H:%M:%S", tz = "GMT")
# Dividimos entre 86400 y 365.25 para que cada año corresponda con la unidad 
prec$time <- as.numeric(prec$time)/86400/365.25 
POTclust <- POT::clust(prec, threshold, tim.cond = 1/365.25, clust.max = T)
indPOT <- POTclust[,3]

## Añadimos delay si hace falta
delayslp <- 24
delay <- 2
# for (delayslp in c(0,1,2,3,4,5,12,24)){
#   for (delay in 0:5) {
#   

nanIndex <- intersect(indPOT, which(is.na(cin$Data[I2, jCoord,iCoord])) + delay)
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
  grid_df <- data.frame(row.names = grid$Dates$start[timeind + delay])
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
cape.df <- gtdf2(cape, K, TT, slp, I2[nanCAPE], delay = delay, delayslp = delayslp)

# Dataframe completo con el CAPE, la precipitacion y el tiempo asociado a ella.
compl.df <- cape.df
compl.df$prec <- pr.Igeldo$hour$data[indDates[nanIndex]]
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
secthres <- seq(threshold, 10, 0.01)

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
    out <- gpd.fit(compl.df$prec[which(clust$cluster == k)], threshold = secthres[thres], npy = 24*365.25, 
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
       path = paste0("D:/Victor/Master/TFM/Graficos/Final Plots/ComplCAPE/", timelapse, "/SLP Delay", delayslp, "/NonStat_shape"),
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
# print(nonstatparamplt)

label <- paste0("pr_ALL_non-stat_param_delay", delay, "_threshold.png")
ggsave(label, plot = nonstatparamplt, 
       path = paste0("D:/Victor/Master/TFM/Graficos/Final Plots/ComplCAPE/", timelapse, "/SLP Delay", delayslp, "/NonStat_param"), 
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
    out <- gpd.fit(compl.df$prec[which(clust$cluster == k)], threshold = secthres[thres], npy = 24*181.25, show = F)
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
# print(statshapeplt)

label <- paste0("pr_ALL_stat_shape_delay", delay, "_threshold.png")
ggsave(label, plot = statshapeplt, 
       path = paste0("D:/Victor/Master/TFM/Graficos/Final Plots/ComplCAPE/", timelapse, "/SLP Delay", delayslp, "/Stat_shape"), 
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
       path = paste0("D:/Victor/Master/TFM/Graficos/Final Plots/ComplCAPE/", timelapse, "/SLP Delay", delayslp), 
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
       path = paste0("D:/Victor/Master/TFM/Graficos/Final Plots/", timelapse, "/SLP Delay", delayslp),
       dpi = 300, bg = "white", width = 7, height = 5, units = "in")

# }
# }
