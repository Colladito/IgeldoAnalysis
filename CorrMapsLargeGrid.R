data <- read.csv("D:/Victor/Master/TFM/Graficos/Corr&Delay/Datos/corr_grid_delay4.csv")

library(sp)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)

datos_sf <- st_as_sf(df_cape1, coords = c("longitud", "latitud"), crs = 4326)

# Cargar el mapa del mundo
mundo <- ne_countries(scale = "large", continent = "Europe", returnclass = "sf")

custom_colors <- c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")

# Create the plot
ggplot(data = mundo) +
  # Base layer with mundo data
  geom_sf() +
  # Layer with datos_sf data, colored by 'corr'
  geom_sf(data = datos_sf, aes(color = TT), size = 18.8,shape = 15) +
  # Adding a specific point with transparency and custom color
  geom_point(data = data.frame(lat = -2.045, long = 43.303), aes(x = lat, y = long), size = 3,alpha = 0.8, color = "black", shape = 17) +
  # Setting coordinate limits
  coord_sf(xlim = c(-3.25, 0.75), ylim = c(42.1, 43.9), expand = F) +
  # Applying a custom color gradient
  scale_color_gradientn(colors = custom_colors) +
  # Applying minimal theme
  theme_minimal() +
  # Adding title and labels
  labs(
    title = "",
    x = "Longitude",
    y = "Latitude",
    color = "CAPE"
  ) +
  # Customizing the theme
  theme(plot.title = element_text(hjust = 0.5)) 


ggsave("GridMeanTT.png",
       path = "D:/Victor/Master/TFM/Graficos/Final Plots/", 
       dpi = 300, bg = "white", width = 8, height = 5, units = "in")




# Day with precipitation "1985-10-25 13:00:00 GMT" of 40.4mm
# Indice I2[nanCAPE[448]]
# Sino poner la media y ya esta mas facil
long <- length(cape$xyCoords$x)
lat <- length(cape$xyCoords$y)

# Create an empty data frame with columns for longitude, latitude, and the values from matrices
df_cape1 <- data.frame(longitud = rep(cape$xyCoords$x, each = lat),
                       latitud = rep(cape$xyCoords$y, long),
                       cape = numeric(lat * long),
                       slp = numeric(lat * long),
                       K = numeric(lat * long),
                       TT = numeric(lat * long))

# Fill in the values from the matrices into the data frame
for (lo in 1:long) {
  for (la in 1:lat) {
    ind <- which(df_cape1$longitud == cape$xyCoords$x[lo] & df_cape1$latitud == cape$xyCoords$y[la])
    df_cape1$cape[ind] <- cape$Data[I2[nanCAPE[448]], la, lo]
    df_cape1$slp[ind] <- slp$Data[I2[nanCAPE[448]], la, lo]
    df_cape1$K[ind] <- K$Data[I2[nanCAPE[448]], la, lo]
    df_cape1$TT[ind] <- TT$Data[I2[nanCAPE[448]], la, lo]
  }
}

long <- length(K$xyCoords$x)
lat <- length(K$xyCoords$y)
kcoordx <- K$xyCoords$x[iCoord+c(-1,0,1)]
kcoordy <- K$xyCoords$y[jCoord+c(-1,0,1)]

# Create an empty data frame with columns for longitude, latitude, and the values from matrices
df_cape1 <- data.frame(longitud = rep(kcoordx, each = 3),
                       latitud = rep(kcoordy, 3),
                       K = numeric(3 * 3),
                       TT = numeric(3 * 3))

# Fill in the values from the matrices into the data frame
for (lo in 1:3) {
  for (la in 1:3) {
    ind <- which(df_cape1$longitud == kcoordx[lo] & df_cape1$latitud == kcoordy[la])
    df_cape1$K[ind] <- mean(K$Data[, la, lo], na.rm = T)
    df_cape1$TT[ind] <- mean(TT$Data[, la, lo], na.rm = T)
  }
