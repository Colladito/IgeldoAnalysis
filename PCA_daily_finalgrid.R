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
# timelapse <- "All-Year"
# indDates <- which(is.element(pr.Igeldo$day$dates,fechas))
# I2 <- which(is.element(fechas, pr.Igeldo$day$dates))

### Indices entre Octubre y Marzo
timelapse <- "Aut-Wint"
IJanMarch <- which(lubridate::month(as.POSIXct(pr.Igeldo$day$dates, tz = "GMT")) <= 3)
IOctDec <- which(lubridate::month(as.POSIXct(pr.Igeldo$day$dates, tz = "GMT")) >= 10)
IOctMarch <- union(IJanMarch, IOctDec)
I2 <- which(is.element(fechas, pr.Igeldo$day$dates[IOctMarch]))
indDates <- which(is.element(pr.Igeldo$day$dates, fechas[I2]))

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

################ Principal Components ################ 

#====================================================#

library(irlba)
pc <- prcomp_irlba(cape.df, n = 10)
summary(pc)

pc1 <- as.matrix(cape.df)%*%pc$rotation[,1]
pc2 <- as.matrix(cape.df)%*%pc$rotation[,2]
pc3 <- as.matrix(cape.df)%*%pc$rotation[,3]
pc4 <- as.matrix(cape.df)%*%pc$rotation[,4]
pc5 <- as.matrix(cape.df)%*%pc$rotation[,5]
pc6 <- as.matrix(cape.df)%*%pc$rotation[,6]
pc7 <- as.matrix(cape.df)%*%pc$rotation[,7]
pc8 <- as.matrix(cape.df)%*%pc$rotation[,8]
pc9 <- as.matrix(cape.df)%*%pc$rotation[,9]
pc10 <- as.matrix(cape.df)%*%pc$rotation[,10]



pca_df <- data.frame(pc1 = pc1, pc2 = pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10)
pca_df$prec <- compl.df$prec

ggplot(pca_df, aes(x = '', y = pc4))+
  geom_boxplot(fill ='#517CD3')+
  geom_jitter(data = subset(pca_df, prec >= 5), aes(size = prec), color = '#C95F9E', alpha = 0.5)+
  theme_minimal()+
  labs(x = '', y = 'PC4', size = 'Prec. (mm)')
lab <- paste0('PC',4,'_daily.png')
ggsave(lab,
       path = "D:/Victor/Master/TFM/Graficos/Final Plots/DailyPCA/", 
       dpi = 300, bg = "white", width = 7, height = 5, units = "in")



################### PRUEBAS ###################

#====================================================#

##################### Clustering  #################### 

#====================================================#

num_pcs <- 1
pc_ind <- 1:num_pcs

km <- kmeans(pca_df[,pc_ind], centers = 2, iter.max = 10000)
clust <- data.frame(row.names = row.names(pca_df),
                    cluster = km$cluster)
pca_df$cluster <- as.factor(clust$cluster)

#====================================================#

################## Aproximacion GPD ################## 

#====================================================#

# Usamos un vector de thresholds donde analizamos los parametros de forma para 
# cada valor del threshold

# Secuencia de tresholds que consideramos
secthres <- seq(threshold, 40, 0.1)

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
    out <- gpd.fit(pca_df$prec[which(clust$cluster == k)], threshold = secthres[thres], npy = 365.25, show = F)
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

# label <- paste0("pr_ALL_stat_shape_delay", delay, "_threshold.png")
# ggsave(label, plot = statshapeplt, 
#        path = paste0("D:/Victor/Master/TFM/Graficos/Final Plots/Daily/", timelapse, "/SLP Delay", delayslp, "/Stat_shape"), 
#        dpi = 300, bg = "white", width = 7, height = 5, units = "in")



#### Boxplots de los peaks en cada cluster


uppind <- which(pca_df$prec >= threshold)
box_all_year <- ggplot(pca_df[uppind,], aes(y = prec, x = cluster, group = cluster, fill = cluster))+
  geom_boxplot()+
  scale_fill_manual(values = c("#517CD3", "#C95F9E"))+
  theme_minimal()+
  labs(fill = "Cluster", y = "Precipitation Over Threshold (mm)", 
       # title = "Boxplots of POT in Each Cluster", 
       x = "Cluster", color = "Station")+
  theme(plot.title = element_text(hjust = 0.5, size = 15))
print(box_all_year)

# label <- paste0("pr_ALL_stat_shape_delay", delay, "_", timelapse, "_boxplot.png")
# ggsave(label, plot = box_all_year, 
#        path = paste0("D:/Victor/Master/TFM/Graficos/Final Plots/Daily/", timelapse, "/SLP Delay", delayslp), 
#        dpi = 300, bg = "white", width = 7, height = 5, units = "in")



### Diferencia cluster entre estaciones

statind <- union(which(lubridate::month(as.POSIXct(row.names(pca_df), tz = "GMT")) <= 3),
                 which(lubridate::month(as.POSIXct(row.names(pca_df), tz = "GMT")) >= 10))
pca_df$station <- "Spr-Summ"
pca_df$station[statind] <- "Aut-Wint"
pca_df$station <- as.factor(pca_df$station)


uppind <- which(pca_df$prec >= threshold)
box_stations <- ggplot(pca_df[uppind,], aes(y = prec, x = cluster, group = interaction(station, cluster), fill = interaction(station, cluster)))+
  geom_boxplot()+
  scale_fill_manual(values = c("#517CD3", "#003399", "#C95F9E", "#8B0A50"))+
  theme_minimal()+
  labs(fill = "Station/Cluster", y = "Precipitation Over Threshold (mm)",
       # title = "Boxplots of POT in Each Cluster",
       x = "Cluster")+
  theme(plot.title = element_text(hjust = 0.5, size = 15))
print(box_stations)

# label <- paste0("pr_ALL_stat_shape_delay", delay, "_stations_boxplot.png")
# ggsave(label, plot = box_stations,
#        path = paste0("D:/Victor/Master/TFM/Graficos/Final Plots/Daily/", timelapse, "/SLP Delay", delayslp),
#        dpi = 300, bg = "white", width = 7, height = 5, units = "in")



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



### NON STATIONARY GPD
# Si no trabajamos con el periodo de Octubre-Marzo habria que sustituir npy por 24*365.25
for (thres in seq_along(secthres)){
  for (k in 1:2) {
    # Non-Stationary
    out <- gpd.fit(pca_df$prec[which(clust$cluster == k)], threshold = secthres[thres], npy = 182, 
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

# label <- paste0("pr_ALL_non-stat_shape_delay", delay, "_threshold.png")
# ggsave(label, plot = nonstatshapeplt, 
#        path = paste0("D:/Victor/Master/TFM/Graficos/Final Plots/Daily/", timelapse, "/SLP Delay", delayslp, "/NonStat_shape"),
#        dpi = 300, bg = "white", width = 7, height = 5, units = "in")

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

# label <- paste0("pr_ALL_non-stat_param_delay", delay, "_threshold.png")
# ggsave(label, plot = nonstatparamplt, 
#        path = paste0("D:/Victor/Master/TFM/Graficos/Final Plots/Daily/", timelapse, "/SLP Delay", delayslp, "/NonStat_param"), 
#        dpi = 300, bg = "white", width = 7, height = 5, units = "in")



#====================================================#

################## Dividing Cluster ################## 

#====================================================#


### Make another cluster in the frechet type tail

indClust1 <- which(pca_df$cluster == 2)
km2 <- kmeans(pca_df[indClust1,pc_ind], centers = 2)
summary(as.factor(km2$cluster))
indClust1_1 <- which(km2$cluster == 1)
indClust1_2 <- which(km2$cluster == 2)
# indClust1_3 <- which(km2$cluster == 3)


# Fit stationary model
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
    out <- gpd.fit(pca_df$prec[indClust1[which(km2$cluster == k)]], threshold = secthres[thres], npy = 365.25, show = F)
    xis2[thres,(3*(k-1))+2] <- out$mle[2]
    xis2[thres,(3*(k-1))+3] <- out$mle[2] - 1.96*out$se[2]
    xis2[thres,(3*(k-1))+4] <- out$mle[2] + 1.96*out$se[2]
  }
}


### Plot del parametro de forma
statshapeplt2 <- ggplot(xis2, aes(x = Threshold)) +
  geom_line(aes(y = Cluster1, color = "Cluster 1"), linewidth = 0.75) +
  geom_line(aes(y = Cluster2, color = "Cluster 2"), linewidth = 0.75) +
  # geom_line(aes(y = Cluster3, color = "Cluster 3"), linewidth = 0.75) +
  # Confidence interval
  geom_ribbon(aes(ymin = LowCluster1, ymax = UpperCluster1, fill = "Cluster 1"), alpha = 0.2)+
  geom_ribbon(aes(ymin = LowCluster2, ymax = UpperCluster2, fill = "Cluster 2"), alpha = 0.2)+
  # geom_ribbon(aes(ymin = LowCluster3, ymax = UpperCluster3, fill = "Cluster 3"), alpha = 0.2)+
  # geom_line(aes(y = LowCluster1, color = "CI 1"), linewidth = 0.75)+
  # geom_line(aes(y = UpperCluster1, color = "CI 1"), linewidth = 0.75)+
  # geom_line(aes(y = LowCluster2, color = "CI 2"), linewidth = 0.75)+
  # geom_line(aes(y = UpperCluster2, color = "CI 2"), linewidth = 0.75)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_color_manual(values = c("Cluster 1" = "#517CD3", 
                                "Cluster 2" = "#C95F9E",
                                "Cluster 3" = "#F08F46")) +
  scale_fill_manual(values = c("Cluster 1" = "#003399", 
                               "Cluster 2" = "#8B0A50",
                               "Cluster 3" = "#D0660E")) +
  theme_minimal() +
  labs(color = "MLE", y = "Shape Parameter", x = "Threshold", fill = "Conf. Interval"
       # title = "Evolution of Shape Parameter Along Threshold"
  )+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(color = guide_legend(order = 1), fill = guide_legend(order = 2))
print(statshapeplt2)

