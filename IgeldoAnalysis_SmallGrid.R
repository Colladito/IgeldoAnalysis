options(java.parameters = "-Xmx8000m")
library(dplyr)
library(tidyr)
library(loadeR)
library(transformeR)
library(visualizeR)
library(ggplot2)


# load("Workspace/SmallGridPCA.RData")

################## Igeldo ################## 

# Latitud: 43º 18' 00" N (43.303) [8]
# Longitud: 02º 01' 59" O (-2.045) [33]
# En las coordenadas de la cape, esta ubicado en [33,8]

igeldo <- read.csv("D:/Victor/Master/TFM/Datos/igeldo.csv")
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


################## CAPE ################## 

# Cargamos CAPE en un grid cercano a Igeldo
ncmlFile <- "D:/Victor/Master/TFM/Datos/ECMWF_ERA5_SFC_clustering.ncml"
# dataInventory(ncmlFile)
cape_igeldo <- loadGridData(dataset = ncmlFile, var = "cape", years = 1940:2015,
                            latLim = c(42, 44), lonLim = c(-3, -1))



################## Analisis ################## 

# Threshold for POTs
# threshold <- quantile(igeldo_hour1940$prec, 0.999)
threshold <- sort(igeldo_hour1940$prec, decreasing = TRUE)[200]
# threshold <- 4


# Indices for POTs
ind_pot <- which(igeldo_hour1940$prec >= threshold)
# Select DELTA T
delta_t <- 0
ind_pot_delta <- ind_pot - delta_t
# Time indices
time_ind_pot <- igeldo_hour1940$X[ind_pot_delta]

# Indices for CAPE
ind_cape_pot <- numeric(length = length(time_ind_pot))
for (i in 1:length(time_ind_pot)) {
  ind <- which(cape_igeldo$Dates[[1]] == time_ind_pot[i])
  ind_cape_pot[i] <- ind
}


# Creamos el objeto climatologico para despues realizar el grafico de correlacion
dataClim <- climatology(cape_igeldo)
for (i in 1:9) {
  for (j in 1:9) {
    dataClim$Data[1,i,j] <- cor(cape_igeldo$Data[ind_cape_pot,i,j], igeldo_hour1940$prec[ind_pot],
                                method = "spearman")
  }
}


## Añadimos el punto de Igeldo para tenerlo localizado en el mapa
test <- data.frame(x = -2.045, y = 43.303)
coordinates(test) <- ~ x + y
pts <- list("sp.points", test, pch=19, cex = 1, col="black")

# Plot con la correlacion de Spearman en el grid cercano a Igeldo
spatialPlot(dataClim, backdrop.theme = "coastline",
            color.theme = "jet.colors", main = "Delay 6",
            scales = list(draw = TRUE),
            at = seq(-.2, .3, by = .025),
            sp.layout = list(pts))
    

################## PCA ################## 

pca <- prinComp(cape_igeldo, keep.orig = FALSE)


# Explained variance 
vexp <- attributes(pca$cape[[1]])$explained_variance
df_vexp <- data.frame(PC = 1:length(vexp), vexp = vexp)
ggplot(df_vexp, aes(PC, 1-vexp))+
  geom_point()+
  theme_minimal()+
  ylab("Fraction of unexplained variability")

ggplot(df_vexp, aes(PC, vexp))+
  geom_point()+
  theme_minimal()+
  ylab("Fraction of explained variability")+
  geom_hline(yintercept = 0.99, color = "darkturquoise")

# Para ver cuantas PC explican el 80% de la varianza
which(vexp <= 0.9)

df_pca <- data.frame(PC1 = pca$cape[[1]]$PCs[,1],
                     PC2 = pca$cape[[1]]$PCs[,2],
                     PC3 = pca$cape[[1]]$PCs[,3],
                     PC4 = pca$cape[[1]]$PCs[,4],
                     PC5 = pca$cape[[1]]$PCs[,5],
                     PC6 = pca$cape[[1]]$PCs[,6],
                     PC7 = pca$cape[[1]]$PCs[,7],
                     PC8 = pca$cape[[1]]$PCs[,8],
                     PC9 = pca$cape[[1]]$PCs[,9],
                     PC10 = pca$cape[[1]]$PCs[,10],
                     PC11 = pca$cape[[1]]$PCs[,11],
                     PC12 = pca$cape[[1]]$PCs[,12])

# Añadimos otra variable que indique si es POT
df_pca$pot <- ifelse(row.names(df_pca) %in% ind_cape_pot, TRUE, FALSE)
df_pca$Time <- attributes(pca)$dates_start
# Añadimos el tamaño de la precipitacion en ese POT (solo hacerlo con delta = 0)
df_pca$prec <- ifelse(df_pca$Time %in% igeldo_hour1940$X, igeldo_hour1940$prec, 1000)




# Boxplot + Jitter con el tamaño de la precipitacion (solo hacerlo con delta = 0)
ggplot(data = df_pca, aes(x = "", y = PC1))+
  geom_boxplot()+
  geom_jitter(data = subset(df_pca, pot == TRUE), aes(color = pot, size = prec)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "deepskyblue3")) +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+ 
  theme(legend.position="none")    

# Boxplot + Jitter sin tamaño de la precipitacion por punto
ggplot(data = df_pca, aes(x = "", y = PC9))+
  geom_boxplot()+
  geom_jitter(data = subset(df_pca, pot == TRUE), aes(color = pot)) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "deepskyblue3")) +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+ 
  theme(legend.position="none")        
## AUMENTAR EL THRESHOLD PARA TENER MENOS PUNTOS QUE SINO NO SE VE NADA Y A VER SI SE VE MEJOR ALGO

# Densidad de cada PC y de los POTs para compararlas
ggplot(data = df_pca, aes(x = PC9)) +
  stat_density(geom = "line") +
  geom_point(data = subset(df_pca, pot == TRUE), aes(y = 0), color = "deepskyblue3") +
  stat_density(data = subset(df_pca, pot == TRUE), geom = "line", color = "deepskyblue3")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")




#### Scatter plots para los POTs y las PCs

# Creamos el dataframe con la precipitacion y las PCs
df_scat <- data.frame(prec = igeldo_hour1940$prec[ind_pot],
                      PC1 = pca$cape[[1]]$PCs[ind_cape_pot,1],
                      PC2 = pca$cape[[1]]$PCs[ind_cape_pot,2],
                      PC3 = pca$cape[[1]]$PCs[ind_cape_pot,3],
                      PC4 = pca$cape[[1]]$PCs[ind_cape_pot,4],
                      PC5 = pca$cape[[1]]$PCs[ind_cape_pot,5],
                      PC6 = pca$cape[[1]]$PCs[ind_cape_pot,6],
                      PC7 = pca$cape[[1]]$PCs[ind_cape_pot,7],
                      PC8 = pca$cape[[1]]$PCs[ind_cape_pot,8],
                      PC9 = pca$cape[[1]]$PCs[ind_cape_pot,9],
                      PC10 = pca$cape[[1]]$PCs[ind_cape_pot,10],
                      PC11 = pca$cape[[1]]$PCs[ind_cape_pot,11],
                      PC12 = pca$cape[[1]]$PCs[ind_cape_pot,12])

# Creamos el Scatterplot
ggplot(df_scat, aes(x = PC6, y = prec))+
  geom_point(color = "deepskyblue3")+
  theme_minimal()+
  labs(title = "ScatterPlot PCs")+
  ylab("Precipitation")
