#===============================================================#

################### AUXILIAR FUNCTIONS FOR TFM ##################

#===============================================================#


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



#===============================================================#

####################### NON-STATIONARY GPD ###################### 

#===============================================================#

# Auxiliar functions related to non-stationary GPD


#### Cumulative Distribution Function for non-stationary GPD 
gpd.cdf <- function(data, time, sigma, alpha, xi){
  f <- 1 - (1 + (xi*data)/(sigma+alpha*time))^(-1/xi)
  return(f)
}

#### Quantile function for non-stationary GPD
gpd.quant <- function(prob, time, sigma, alpha, xi){
  f <- ((sigma + alpha*time) / xi)*((1-prob)^(-xi)-1)
  return(f)
}



#### Function to compute the mle estimators for a non stationary gpd
# data: data which has to contain 2 columns, "time" and "prec"
# params = (sigma_0, alpha, chi)
# alpha: significance level for confidence interval
# ci: if TRUE make the confidence interval
# stat: function for the scale parameter, default identity 
#       with exp we make the scale parameter always positive
gpd.mle <- function(data, threshold, params, alpha = 0.05, ci = F, stat = identity){
  time <- data$time
  prec <- data$prec
  
  z <- list()
  
  # Negative likelihood function
  op.loglike <- function(params) {
    # Extract columns from the data frame
    sc <- stat(params[1] + params[2] * time)
    xi <- params[3]
    y <- 1 + xi*(prec - threshold)/sc
    
    if (min(sc) <= 0) {
      l <- 10^6  # if sc <0 unfeasible
    }
    else {
      if (min(y) <= 0) { # if y<0 unfeasible
        l <- 10^6
      } else {
        # log-likelihood calculation based on parameters and data
        l <- sum(log(sc)) + (1 / xi + 1) * sum(log(y))
      }
    }
    
    # return negative log-likelihood 
    return(l)
  }
  
  result <- optim(par = params, fn = op.loglike, hessian = T)
  
  z$convergence <- result$convergence
  z$mle <- result$par
  z$hessian <- result$hessian
  z$nllh <- result$value
  
  # Confidence interval
  if(ci){
    # Fisher Info (Cov matrix)
    fish_inf <- solve(result$hessian)
    se <- sqrt(diag(fish_inf))
    upper <- result$par + qnorm(alpha/2, lower.tail = F)*se
    lower <- result$par - qnorm(alpha/2, lower.tail = F)*se
    
    z$ci <- list(sig0 = c(lower[1], upper[1]),
                 alpha = c(lower[2], upper[2]),
                 xi = c(lower[3], upper[3]))
  }
  
  
  return(z)
}


#===============================================================#

################ FUNCTION TO COMPUTE EVERYTHING  ################ 

#===============================================================#

# This function compute the mle estimators for the GPD for all the
# clusters. Also, we can change the number of clusters, threshold
# and delay of the CAPE variable with the precipitation. 

# If we add more variables, the function should change.

# gridvars: list containing all the grid variables (CAPE, CIN) 
#          the first element containing the CAPE and second element CIN
# precip: data containing the precipitation
# threshold: threshold where to make the POTs 
# delay: between the gridvars and precipitation default 0
# k: number of clusters to compute using kmeans default 2
# params: initial parameters for computing the MLE by default:
#          0.1  scale parameter, sigma_0
#          0    alpha (time)
#          0.1) shape parameter, xi
CAPEclustPrec <- function(cape, cin, precip, threshold = 0.8, delay = 0, clusters, k = 2,
                          params = c(0.1,0,0.1), ci = F){
  
  
  pr.Igeldo <- precip
  
  # Homogeneizamos las fechas del grid
  fechas <- intersect(cin$Dates$start,cape$Dates$start)
  I1 <- which(is.element(cin$Dates$start, fechas))
  I2 <- which(is.element(cape$Dates$start, fechas))
  
  # Cogemos las fechas para cada variable
  cin$Data[] <- cin$Data[I1,,]
  cin$Dates$start <- cin$Dates$start[I1]
  cin$Dates$end <- cin$Dates$end[I1]
  
  cape$Data <- cape$Data[I2,,]
  attr(cape$Data, "dimensions") <- c("time","lat","lon")
  cape$Dates$start <- cape$Dates$start[I2]
  cape$Dates$end <- cape$Dates$end[I2]
  
  
  # Indice x horas
  indDates <- which(is.element(pr.Igeldo$hour$dates,fechas))
  I2 <- which(is.element(fechas, pr.Igeldo$hour$dates))
  # cape.day <- aggregateGrid(grid = cape, aggr.d = list(FUN = "mean" , na.rm = TRUE))
  
  igeldo <- c(-2.009944 ,43.322886)
  iCoord <- which.min(abs(cin$xyCoords$x - igeldo[1]))
  jCoord <- which.min(abs(cin$xyCoords$y - igeldo[2]))
  
  
  # Seleccionamos los POTS
  prec <- data.frame(obs = pr.Igeldo$hour$data[indDates])
  
  min_date <- as.POSIXct(min(pr.Igeldo$hour$dates[indDates]), format = "%Y-%m-%d %H:%M:%S", tz = "GMT")
  prec$time <- as.numeric(difftime(as.POSIXct(pr.Igeldo$hour$dates[indDates], 
                                              format = "%Y-%m-%d %H:%M:%S", tz = "GMT"), 
                                   min_date, units = "days") / 365.25)
  
  # Hacemos el declustering para nuestros datos
  POTclust <- POT::clust(prec, threshold, tim.cond = 1/365.25, clust.max = T)
  indPOT <- POTclust[,3]
  
  
  # Añadimos el delay a los indices
  nanIndex <- intersect(indPOT, which(is.na(cin$Data[I2, jCoord,iCoord])) + delay)
  nanCAPE <- nanIndex - delay
  
  
  # Creamos el dataframe con el cape (Funcion creada arriba)
  cape.df <- gtdf(cape, I2[nanCAPE], delay)
  
  # Creamos df auxiliar para hacer usar el clustering posteriormente
  compl.df <- cape.df
  compl.df$prec <- pr.Igeldo$hour$data[indDates[nanIndex]]
  compl.df$time <- prec$time[indDates[nanIndex]]
  
  ### Clustering
  # Realizamos el clustering
  # km <- kmeans(cape.df, centers = k, iter.max = 100000)
  # # DF auxiliar incluyendo los clusters
  # compl.df1 <- compl.df
  # compl.df1$cluster <- km$cluster
  
  ### AÑadimos los clusters
  compl.df1 <- compl.df
  # indcapeclust <- intersect(row.names(cape.df), row.names(clusters))
  compl.df1$cluster <- clusters$cluster[I2[nanCAPE]]
  
  compl.df1$datetime_index <- as.POSIXct(rownames(compl.df1), format = "%Y-%m-%d %H:%M:%S", tz = "GMT")
  min_date <- as.POSIXct(min(compl.df1$datetime_index))
  compl.df1$time <- as.numeric(difftime(compl.df1$datetime_index, min_date, units = "days") / 365.25)
  
  z <- list()
  z$mle <- data.frame(sigma = rep(NA, k),
                      alpha = rep(NA, k),
                      xi = rep(NA, k),
                      row.names = paste0("Cluster", 1:k))
  
  if (ci) {
    z$ci <- list()
    z$ci$sigma <- data.frame(lower = rep(NA,k),
                             upper = rep(NA,k),
                             row.names = paste0("Cluster", 1:k))
    z$ci$alpha <- data.frame(lower = rep(NA,k),
                             upper = rep(NA,k),
                             row.names = paste0("Cluster", 1:k))
    z$ci$xi <- data.frame(lower = rep(NA,k),
                          upper = rep(NA,k),
                          row.names = paste0("Cluster", 1:k))
  }
  
  for (i in 1:k) {
    df <- data.frame(time = compl.df1$time[which(compl.df1$cluster == i)], 
                      prec = compl.df1$prec[which(compl.df1$cluster == i)], 
                      row.names = rownames(compl.df1[which(compl.df1$cluster == i),]))
 
    mle1 <- gpd.mle(df, threshold, params, ci = ci)
    z$mle[i,] <- mle1$mle
    if (ci) {
      z$ci$sigma[i,] <- mle1$ci$sig0
      z$ci$alpha[i,] <- mle1$ci$alpha
      z$ci$xi[i,] <- mle1$ci$xi
    }
  }
  
  
  return(z)
} 





#===============================================================#

################ AUTOMATIC CLUSTERING FOR MRLP ################ 

#===============================================================#



################# KMEANS + MRLP ################# 

# Function which computes the kmeans algorithm for the data you introduce and 
# plot the MRLP for the different clusters associated to the precipitation
kmeans.mrlp <- function(data, prec, k, plot = TRUE){
  km <- kmeans(data, centers = k)
  
  data.prec <- data
  data.prec$prec <- prec
  data.prec$cluster <- km$cluster
  
  if(plot){
    for (i in 1:k) {
      label <- paste("MRLP of Cluster", i)
      POT::mrlplot(data.prec$prec[which(data.prec$cluster == i)],
                   main = label)
    }
  }
  return(data.prec)
}




################# KMEDOIDS + MRLP ################# 
library(cluster)
kmed.mrlp <- function(data, prec, k, plot = TRUE, metric = "euclidean"){
  kmed <- pam(data, k = k, metric = metric)
  
  data.prec <- data
  data.prec$prec <- prec
  data.prec$cluster <- kmed$cluster
  
  if(plot){
    for (i in 1:k) {
      label <- paste("MRLP of Cluster", i)
      POT::mrlplot(data.prec$prec[which(data.prec$cluster == i)],
                   main = label)
    }
  }
  return(data.prec)
}



################# SOM + MRLP ################# 

library(kohonen)
som.mrlp <- function(data, prec, k, som.plot.prec = TRUE, som.plot = FALSE, mrlp.plot = FALSE,
                     xdimgrid = 4, ydimgrid = 4, topology = "hexagonal"){
  
  # Make the som object
  som <- som(as.matrix(data), 
             grid = somgrid(xdim = xdimgrid, ydim = xdimgrid, topology))
  
  # Prepare the data
  data.prec <- data
  data.prec$prec <- prec
  data.prec$cluster <- som$unit.classif
  
  
  # SOM plot of the cape in the nearest point to igeldo
  if(som.plot){
    ngrid <- xdimgrid*ydimgrid
    centroids <- numeric(ngrid)
    for (i in 1:ngrid){
      centroids[i] <- mean(data.prec[which(data.prec$clust == i), 26])
    }
    plot(som, type = "property", property = centroids, 
         main = paste("Cape Nearliest Point to Igeldo", colnames(getCodes(som))[26]), 
         shape = "straight", palette.name = heat.colors)
  }
  
  # SOM plot for precipitation
  if(som.plot.prec){
    ngrid <- xdimgrid*ydimgrid
    prec.centroids <- numeric(ngrid)
    for (i in 1:ngrid){
      prec.centroids[i] <- mean(data.prec$prec[which(data.prec$clust == i)])
    }
    plot(som, type = "property", property = prec.centroids, 
         main = "Precipitation in Igeldo", 
         shape = "straight", palette.name = heat.colors)
    
    
  }
  
  
  # MRLP plots for the k greatest average precipitation clusters 
  if(mrlp.plot){
    ngrid <- xdimgrid*ydimgrid
    prec.centroids <- numeric(ngrid)
    for (i in 1:ngrid){
      prec.centroids[i] <- mean(data.prec$prec[which(data.prec$clust == i)])
    }
    
    k.great.clust <- order(prec.centroids, decreasing = TRUE)[1:k]
    
    for (i in k.great.clust) {
      label <- paste("MRLP of Cluster", i)
      POT::mrlplot(data.prec$prec[which(data.prec$cluster == i)],
                   main = label)
    }
  }
  return(list(data = data.prec, kclust = k.great.clust))
}


################# MAX-DISS + MRLP ################# 



### Auxiliar functions
# Función que calcula la distancia normalizada
# M es la matriz 
distancia_normalizada <- function(M, D, escalar) {
  # Get the dimensions of matrix M
  dim_M <- dim(M)
  
  # Initialize a matrix 'dif' with zeros
  dif <- matrix(0, nrow = dim_M[1], ncol = dim_M[2])
  
  # Calculate differences for specified columns
  for (i in escalar) {
    dif[, i] <- D[, i] - M[, i]
  }
  
  # Rest of your function implementation...
  # You can add more code here based on the logic you want to implement
  
  # For example, if you want to calculate the normalized distance, you can do something like:
  Dist <- sqrt(rowSums(dif^2))
  
  # Return the result
  return(Dist)
}



## Funcion que calcula los centros de los clusters de la max diss

MaxDiss <- function(semilla, num, train_n, escalar) {
  # Check if there is more than one seed
  if (length(semilla) > 1) {
    print(paste("You have more than one semilla", length(semilla)))
    semilla <- semilla[1]
  }
  
  # Initialization of centroids
  subset <- train_n[semilla,]
  
  # Remove the selected centroid from the training data
  train_n2 <- train_n[-semilla,]
  
  # Repeat the process until the desired number of centroids is reached
  ncenters <- 1
  qerr <- 0
  
  while (ncenters < num) {
    dim_train_n2 <- dim(train_n2)
    m <- rep(1, dim_train_n2[1])
    
    dim_subset <- dim(subset)
    if (dim_subset[1] == 1) {
      xx2 <- subset[m, ]
      dultima <- distancia_normalizada(train_n2, xx2, escalar)
    } else {
      xx <- subset[dim_subset[1], ]
      xx2 <- xx[m, ]
      danterior <- distancia_normalizada(train_n2, xx2, escalar)
      dultima <- min(danterior, dultima)
    }
    
    bmu <- which.max(dultima)
    qerr <- max(dultima)
    
    if (!is.na(qerr)) {
      subset <- rbind(subset, train_n2[bmu, , drop = FALSE])
      train_n2 <- train_n2[-bmu, , drop = FALSE]
    }
    
    ncenters <- dim(subset)[1]
  }
  
  return(list(subset = subset, ncenters = ncenters))
}


###### Complete maxdiss algorith with MRLP

maxdiss.mrlp <- function(data, prec, k, plot = TRUE){
  data.prec <- data
  data.prec$prec <- prec
  indmax <- which.max(prec)
  
  sample <- MaxDiss(indmax, k, data, 1:dim(data)[2])
  
  dist <- data.frame(row.names = row.names(data))
  for (clust in 1:k) {
    label <- paste0("dist",clust)
    dist[label]<- sapply(1:dim(data.prec)[1], function(y){sqrt(sum((data[y,]-sample$subset[clust,])^2))})
  }
  clust <- apply(dist, 1, which.min)
  
  data.prec$cluster <- clust
  
  if(plot){
    for (i in 1:k) {
      label <- paste("MRLP of Cluster", i)
      POT::mrlplot(data.prec$prec[which(data.prec$cluster == i)],
                   main = label)
    }
  }
  return(data.prec)
}
