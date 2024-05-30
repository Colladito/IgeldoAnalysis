
#==============================================#
################## SIMULATION ##################
#==============================================#


# Set the different parameters for each distribution
mu <- 15
sig <- 2
xi <- -0.1

mu2 <- 5
sig2 <- 6
xi2 <- -0.05

mu-sig/xi

# Number of samples
ny <- 100

# Use fExtremes to use basic functions of GEV distribution
library(fExtremes)

# Generate the samples of size ny
R <- rgev(ny, xi = xi, mu = mu, beta = sig)
R2 <- rgev(ny, xi = xi2, mu = mu2, beta = sig2)

# Create probability vector (length of the samples)
prob <- seq(1:ny)/(ny+1)

# Time stamp
t <-  1/(1-prob)

# Sort the data
Rsort <- sort(R)
Rsort2 <- sort(R2)
R3 <- pmax(Rsort, Rsort2)
Rsort3 <- sort(R3)

# Create probability vector
Prob2 <- seq(1, 100000, length = 1000) / (100000 + 1)
T2 <- 1 / (1 - Prob2)


# Defining auxiliarly functions
# Aggregated cdf (independent events)
aggregatedGEV_CDF <- function(x, xi, sig, mu, xi2, sig2, mu2) {
  pgev(x, xi, mu, sig) * pgev(x, xi2, mu2, sig2)
}


# Invert the CDF of both GEV to obtain the quantiles
aggregatedGEV_INV <- function(probs, xi, sig, mu, xi2, sig2, mu2) {
  # Initialize the vector of quantiles
  x <- numeric(length(probs))
  
  for (i in 1:length(x)) {
    # Set middle point of both quantiles
    ini <- (qgev(probs[i], xi, mu, sig) + qgev(probs[i], xi2, mu2, sig2)) / 2
    
    # Find the root of the probability and aggregated dist
    result <- uniroot(function(y){
      probs[i] - aggregatedGEV_CDF(y, xi, sig, mu, xi2, sig2, mu2)},
      interval = c(ini,10^10))

    x[i] <- result$root
  }
  
  return(x)
}



# Compute the quantiles of the fixed probabilities for all GEV and aggregated
Rfit <- qgev(Prob2, xi, mu, sig)
Rfit2 <- qgev(Prob2, xi2, mu2, sig2)
Rfit3 <- aggregatedGEV_INV(Prob2, xi, sig, mu, xi2, sig2, mu2)


# Library ismev to fit the data
library(ismev)

# Fit GEV distribution for R3 max of two Weibull
fit_R3 <- gev.fit(R3)
parmhat <- fit_R3$mle
# Since xi > 0 Frechet, without asymptotic limit
Rfit4 <- qgev(Prob2, mu = parmhat[1], beta = parmhat[2], xi = parmhat[3])


# Fit GEV distribution for R 
fit_R <- gev.fit(R)
parmhat5 <- fit_R$mle
# Como xi < 0 Weibull 
Rfit5 <- qgev(Prob2, mu = parmhat5[1], beta = parmhat5[2], xi = parmhat5[3])
# Asymptotic limit of Fitted Weibull 1
asint1 <- parmhat5[1]+parmhat5[2]/abs(parmhat5[3])

# Fit GEV distribution for R2
fit_R2 <- gev.fit(R2)
parmhat6 <- fit_R2$mle
# xi < 0 Weibull
Rfit6 <- qgev(Prob2, mu = parmhat6[1], beta = parmhat6[2], xi =  parmhat6[3])
# Asymptotic limit of Fitted Weibull 2
asint2 <- parmhat6[1]+parmhat6[2]/abs(parmhat6[3])




######## Plots ######## 
# Auxiliarly Data Frames
data <- data.frame(years = T2,
                   RealWei1 = Rfit,
                   RealWei2 = Rfit2,
                   RealCombWei = Rfit3,
                   AdjFrech = Rfit4,
                   AdjWei1 = Rfit5,
                   AdjWei2 = Rfit6)
data2 <- data.frame(years = t,
                    R1 = Rsort,
                    R2 = Rsort2,
                    MaxR1R2 = Rsort3)

library(ggplot2)
library(scales)
library(latex2exp)

azul <- "#517CD3"
morado <- "#C95F9E"
naranja <- "#F08F46"

ggplot(data, aes(x = years)) +
  # Plot real distributions
  geom_line(aes(y = RealWei1, color = "Weibull 1", linetype = "Real"), size = 0.75) +
  geom_line(aes(y = RealWei2, color = "Weibull 2", linetype = "Real"), size = 0.75) +
  geom_line(aes(y = RealCombWei, color = "Aggregated", linetype = "Real"), size = 0.75) +
  # Plot fitted distributions
  geom_line(aes(y = AdjFrech, color = "Frechet", linetype = "Fitted"), size = 0.75) +
  geom_line(aes(y = AdjWei1, color = "Weibull 1", linetype = "Fitted"), size = 0.75) +
  geom_line(aes(y = AdjWei2, color = "Weibull 2", linetype = "Fitted"), size = 0.75) +
  # Plot samples 
  geom_point(data = data2, aes(x = years, y = R1), color = azul, 
             alpha = 0.5, size = 1)+
  geom_point(data = data2, aes(x = years, y = R2), color = morado, 
             alpha = 0.5, size = 1)+
  geom_point(data = data2, aes(x = years, y = MaxR1R2), color = naranja,
             size = 2, shape = 21)+
  # Asymptotic limit
  geom_hline(yintercept = asint1, color = azul)+
  theme_bw() +
  annotate(geom="text", x=3*10^4, y=65, label=paste("Asymptotic limit:", round(asint2)),
           color=morado, size=5)+
  annotate(geom="text", x=4, y=asint1-7, label=paste("Asymptotic limit:", round(asint1)),
           color=azul, size=5)+
  annotate(geom="text", x=2*10^4, y=125, label=TeX(paste("Fitted tail $\\xi$:", round(parmhat[3], 4))),
           color="black", size=5)+
  #ylim(0, 60) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = c("Weibull 1" = azul, "Weibull 2" = morado,
                                "Aggregated" = naranja, "Frechet" = "black",
                                "Fitted Weibull 1" = azul, "Fitted Weibull 2" = morado)) +
  scale_linetype_manual(values = c("Real" = "solid", "Fitted" = "dashed")) +
  labs(color = "Distribution", linetype = "Type") +
  labs(title = "", x = "Return Years", y = "Annual maximum precipitation (mm)")+
  theme(plot.title = element_text(hjust = 0.5))


# Uncoment to save the before plot  
# ggsave("ExamplePMP.png", dpi = 300, path="Graficos/Final Plots/",
#        bg = "white", width = 9, height = 5, units = "in")
