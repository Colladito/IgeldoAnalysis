options(java.parameters = "-Xmx10g")
library(evd)
library(POT)
library(dplyr)
library(tidyr)
library(loadeR)
library(transformeR)
library(visualizeR)
library(ggplot2)
library(MASS)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(fpc)

################## Igeldo ################## 
# Latitud: 43º 18' 00" N (43.303) [8]
# Longitud: 02º 01' 59" O (-2.045) [33]
# En las coordenadas de la cape, esta ubicado en [33,8]

igeldo <- read.csv("Datos/igeldo.csv")
igeldo$prec <- igeldo$prec %>%
  replace_na(0)

igeldo$X <- as.POSIXlt(igeldo$X, format = "%Y-%m-%d %H:%M:%S")

# Extract year, month, day, hour, minute, and second
igeldo$Year <- format(igeldo$X, "%Y")
igeldo$Month <- format(igeldo$X, "%m")
igeldo$Day <- format(igeldo$X, "%d")
igeldo$Hour <- format(igeldo$X, "%H")
igeldo$Minute <- format(igeldo$X, "%M")

#Transform the data from 10-minutes to 1-hour
igeldo_hour <- igeldo %>%
  group_by(Year, Month, Day, Hour) %>%
  summarise(prec = max(prec))
igeldo_hour$X <- as.POSIXlt(paste(igeldo_hour$Year, 
                                  igeldo_hour$Month, 
                                  igeldo_hour$Day, 
                                  igeldo_hour$Hour, sep = "-"), 
                            format = "%Y-%m-%d-%H")
igeldo_hour$X <- paste(igeldo_hour$X, "GMT", sep = " ")

# Extract Igeldo data for 1940-2016(ultimo año)
igeldo_hour1940 <- igeldo_hour %>%
  filter(Year >= 1940) %>%
  filter(Year <= 2015)


# Load CAPE in Igeldo
ncmlFile <- "Datos/ECMWF_ERA5_SFC_clustering.ncml"
# dataInventory(ncmlFile)
cape_igeldo <- loadGridData(dataset = ncmlFile, var = "cape", years = 1940:2015, 
             latLim = 43.25, lonLim = -2)

# Threshold for POTs
# threshold <- quantile(igeldo_hour1940$prec, 0.999)
# threshold <- sort(igeldo_hour1940$prec, decreasing = TRUE)[770]
threshold <- 4


# Indices for POTs
ind_pot <- which(igeldo_hour1940$prec >= threshold)
# Select DELTA T
delta_t <- 4
ind_pot_delta <- ind_pot - delta_t
# Time indices
time_ind_pot <- igeldo_hour1940$X[ind_pot_delta]

# Indices for CAPE
ind_cape_pot <- numeric(length = length(time_ind_pot))
for (i in 1:length(time_ind_pot)) {
  ind <- which(cape_igeldo$Dates[[1]] == time_ind_pot[i])
  ind_cape_pot[i] <- ind
}

## BOXCOX transformations
# Precipitation
b_prec <- boxcox(lm(igeldo_hour1940$prec[ind_pot] ~ 1), plotit = FALSE)
lambda_prec <- b_prec$x[which.max(b_prec$y)]
# lambda_prec
prec_boxcox <- (igeldo_hour1940$prec[ind_pot] ^ lambda_prec - 1) / lambda_prec

# CAPE
cape_igel_pos <- cape_igeldo$Data[ind_cape_pot] + 0.00000001
b_cape <- boxcox(lm(cape_igel_pos ~ 1), plotit = FALSE)
lambda_cape <- b_cape$x[which.max(b_cape$y)]
# lambda_cape
cape_boxcox <- (cape_igel_pos ^ lambda_cape - 1) / lambda_cape



df <- data.frame(cape = cape_igeldo$Data[ind_cape_pot],
                 cape_boxcox = cape_boxcox, 
                 prec = igeldo_hour1940$prec[ind_pot],
                 prec_boxcox = prec_boxcox)

# Mirar para aplicar alguna transformacion a las variables (log, sqrt,...)
corr <- cor(df$prec, df$cape_boxcox)
tit <- paste("Scatterplot betwenn POTs in Prec and CAPE with delay = ", delta_t , " cor = ", corr, sep = "")
ggplot(df, aes(x = cape, y = prec))+
  geom_point(color = "deepskyblue3")+
  theme_minimal()+
  labs(title = tit)+
  xlab("CAPE")+
  ylab("Precipitation")

## La mayor correlacion se alcanza con delta_t = 2, es decir 2 horas de delay.



### ScatterPlot con Clusters
resultados <- regmix(df$prec, df$cape_boxcox, nclust = 3)

# Visualizar los resultados del clustering
df$g=resultados$g
df$g=as.factor(df$g)

ggplot(df)+
  aes(x=cape,y=prec,color=g,shape=g)+
  geom_point()+
  labs(title = "ScatterPlot, delay 4")+
  theme_minimal()





# Bucle para ver la correlacion
df2 <- data.frame(Delta = 0:20)
ind_pot <- which(igeldo_hour1940$prec >= threshold)
for (delta in 0:20) {
  # Select DELTA T
  delta_t <- delta
  ind_pot_delta <- ind_pot - delta_t
  # Time indices
  time_ind_pot <- igeldo_hour1940$X[ind_pot_delta]
  
  # Indices for CAPE
  ind_cape_pot <- numeric(length = length(time_ind_pot))
  for (i in 1:length(time_ind_pot)) {
    ind <- which(cape_igeldo$Dates[[1]] == time_ind_pot[i])
    ind_cape_pot[i] <- ind
  }
  
  cape_igel_pos <- cape_igeldo$Data[ind_cape_pot] + 0.00000001
  b_cape <- boxcox(lm(cape_igel_pos ~ 1), plotit = FALSE)
  lambda_cape <- b_cape$x[which.max(b_cape$y)]
  lambda_cape
  
  if (lambda_cape <= 0.001 && lambda_cape >= -0.001) {
    cape_boxcox <- log(cape_igel_pos)
  } else {
    cape_boxcox <- (cape_igel_pos ^ lambda_cape - 1) / lambda_cape
  }
  
  df <- data.frame(cape = cape_boxcox, 
                   prec = igeldo_hour1940$prec[ind_pot])
  
  corr <- cor(igeldo_hour1940$prec[ind_pot], cape_igeldo$Data[ind_cape_pot])
  corr1 <- cor(df$prec, df$cape)
  corr2 <- cor(prec_boxcox ,df$cape)
  
  df2$corr[delta+1] <- corr
  df2$corr_boxcox[delta+1] <- corr1
  df2$corr_boxcox2[delta+1] <- corr2
}

write.csv(df2, file = "Correlacion_Igeldo.csv", row.names = FALSE)


#### NO FUNCIONA BIEN
######## Matrix Loop ######
# Bucle para todas las posiciones del grid

longitud <- seq(-6, 0, by = 0.25)
latitud <- seq(46, 41.5, by = -0.25)

corrmatrix <- array(dim = c(17, length(longitud), length(latitud)))
corrmatrix_boxcox1 <- array(dim = c(17, length(longitud), length(latitud)))
corrmatrix_boxcox2 <- array(dim = c(17, length(longitud), length(latitud)))

ncmlFile <- "Datos/ECMWF_ERA5_SFC_clustering.ncml"

ind_pot <- which(igeldo_hour1940$prec >= threshold)

b_prec <- boxcox(lm(igeldo_hour1940$prec[ind_pot] ~ 1), plotit = FALSE)
lambda_prec <- b_prec$x[which.max(b_prec$y)]
#lambda_prec
prec_boxcox <- (igeldo_hour1940$prec[ind_pot] ^ lambda_prec - 1) / lambda_prec

for (long in 1:length(longitud)) {
  for (lat in 1:length(latitud)) {
    cape_igeldo <- loadGridData(dataset = ncmlFile, var = "cape", years = 1940:2015, 
                                latLim = latitud[lat], lonLim = longitud[long])
    
    for (delta in 0:16) {
      # Select DELTA T
      delta_t <- delta
      ind_pot_delta <- ind_pot - delta_t
      # Time indices
      time_ind_pot <- igeldo_hour1940$X[ind_pot_delta]
      
      # Indices for CAPE
      ind_cape_pot <- numeric(length = length(time_ind_pot))
      for (i in 1:length(time_ind_pot)) {
        ind <- which(cape_igeldo$Dates[[1]] == time_ind_pot[i])
        ind_cape_pot[i] <- ind
      }
      
      # Transformacion boxcox
      cape_igel_pos <- cape_igeldo$Data[ind_cape_pot] + 0.00000001
      b_cape <- boxcox(lm(cape_igel_pos ~ 1), plotit = FALSE)
      lambda_cape <- b_cape$x[which.max(b_cape$y)]
      lambda_cape
      
      if (lambda_cape <= 0.001 && lambda_cape >= -0.001) {
        cape_boxcox <- log(cape_igel_pos)
      } else {
        cape_boxcox <- (cape_igel_pos ^ lambda_cape - 1) / lambda_cape
      }
      
      
      # Calculo de correlaciones
      
      # Correlacion sin transformaciones
      corr <- cor(igeldo_hour1940$prec[ind_pot], cape_igeldo$Data[ind_cape_pot])
      # Correlacion con transformacion en el CAPE
      corr1 <- cor(igeldo_hour1940$prec[ind_pot], cape_boxcox)
      # Correlacion con transformacion en PREC y CAPE
      corr2 <- cor(prec_boxcox, cape_boxcox)
      
      corrmatrix[delta+1, lat, long] <- corr
      corrmatrix_boxcox1[delta+1, lat, long] <- corr1
      corrmatrix_boxcox2[delta+1, lat, long] <- corr2
    }
  }
}



############# Bucle 2 ###################
longitud <- seq(-10, 0, by = 0.25)
latitud <- seq(50, 41.5, by = -0.25)


data <- expand.grid(lat = latitud, long = longitud)
# Select DELTA T
delta_t <- 6
j <- 1

ind_pot_delta <- ind_pot - delta_t
# Time indices
time_ind_pot <- igeldo_hour1940$X[ind_pot_delta]

# Indices for CAPE
ind_cape_pot <- numeric(length = length(time_ind_pot))
for (i in 1:length(time_ind_pot)) {
  ind <- which(cape_igeldo$Dates[[1]] == time_ind_pot[i])
  ind_cape_pot[i] <- ind
}

for (long in 1:length(longitud)) {
  for (lat in 1:length(latitud)) {
    cape_igeldo <- loadGridData(dataset = ncmlFile, var = "cape", years = 1940:2015, 
                                latLim = latitud[lat], lonLim = longitud[long])
    
    # Transformacion boxcox
    cape_igel_pos <- cape_igeldo$Data[ind_cape_pot] + 0.00000001
    b_cape <- boxcox(lm(cape_igel_pos ~ 1), plotit = FALSE)
    lambda_cape <- b_cape$x[which.max(b_cape$y)]
    lambda_cape
    
    if (lambda_cape <= 0.001 && lambda_cape >= -0.001) {
      cape_boxcox <- log(cape_igel_pos)
    } else {
      cape_boxcox <- (cape_igel_pos ^ lambda_cape - 1) / lambda_cape
    }
    
    
    # Calculo de correlaciones y añadimos al dataframe

    data$corr[j] <- cor(igeldo_hour1940$prec[ind_pot], cape_igeldo$Data[ind_cape_pot])
    data$corr_box1[j] <- cor(igeldo_hour1940$prec[ind_pot], cape_boxcox)
    data$corr_box2[j] <- cor(prec_boxcox, cape_boxcox)

    print(j)
    
    j <- j+1
    
  }
}

## Para 475 datos del grid
# Start: 13:53:56
# Finish: 16:22:01


## Para el grid completo (7:15h)
# Start: 17:40:02
# Finish: 00:55:48

write.csv(data, file = "Graficos/Corr&Delay/Datos/corr_grid_delay6.csv", row.names = FALSE)

# Leemos los datos de la correlacion ya generados
# data <- read.csv("Graficos/Corr&Delay/Datos/corr_grid_delay0.csv")

datos_sf <- st_as_sf(data, coords = c("long", "lat"), crs = 4326)

# Cargar el mapa del mundo
mundo <- ne_countries(scale = "large", continent = "Europe", returnclass = "sf")

# Crear el mapa
ggplot(data = mundo) +
  geom_sf() +
  geom_sf(data = datos_sf, aes(color = corr), size = 3) +
  geom_point(data = data.frame(lat = -2.045, long = 43.303), aes(x = lat, y = long), alpha = 0.4)+
  coord_sf(xlim = c(-12, 5), ylim = c(41, 51)) +
  scale_color_viridis_c(limits=c(-0.2,0.2)) +
  theme_minimal() +
  labs(title = "Delay 0")

ggplot(data = mundo) +
  geom_sf() +
  geom_sf(data = datos_sf, aes(color = corr_box1), size = 3) +
  geom_point(data = data.frame(lat = -2.045, long = 43.303), aes(x = lat, y = long), alpha = 0.4)+
  coord_sf(xlim = c(-12, 5), ylim = c(41, 51)) +
  scale_color_viridis_c(limits=c(-0.2,0.2)) +
  theme_minimal() +
  labs(title = "Delay 0")

ggplot(data = mundo) +
  geom_sf() +
  geom_sf(data = datos_sf, aes(color = corr_box2), size = 3) +
  geom_point(data = data.frame(lat = -2.045, long = 43.303), aes(x = lat, y = long), alpha = 0.4)+
  coord_sf(xlim = c(-12, 5), ylim = c(41, 51)) +
  scale_color_viridis_c(limits=c(-0.2,0.25)) +
  theme_minimal() +
  labs(title = "Delay 0")

  

