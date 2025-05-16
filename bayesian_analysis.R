#--- Bayesian Inference ----

setwd("D:/Github/Thesis_job")

lista_de_paquetes <- c("readxl", "ggplot2", "readr", "tidyverse") # Reemplaza con tus paquetes

for (paquete in lista_de_paquetes) {
  if(!require(paquete, character.only = TRUE)){
    install.packages(paquete)
    require(paquete, character.only = TRUE)
  }
}

data1 <- read_excel("Acceptance Sampling MIL-STD 105E for Quality Control.xlsx",
                    sheet = "tab")
data2 <- read_excel("Jurnal Rekavasi.xlsx", sheet = "tabla")

#----1. Metropolis-Hastings ----

x0 <- c(0.5); burn_in <- 5000; n_samples <- 1000; thinning <- 1

samples = NULL
x_currente <- array(data = NA, dim = c(1, length(x0)), dimnames = NULL)

total_iterations = burn_in + n_samples * thinning 
accepted_samples <- NULL

i <- 0

while(i < total_iterations){
  # prior propuesta
  x_proposed <- proposal_sampler(x_current)
  # Compute denominator of the acceptance ratio
  denominator <- target_pdf(x_current) * proposal_pdf(x_proposed, x_current)
  
  # Check if the denominator is nonzero
  if (denominator > 0){
    acceptance_ratio = (target_pdf(x_proposed) * proposal_pdf(x_current, x_proposed)) / denominator
  acceptance_ratio = min(1, acceptance_ratio)
  } else{
  # Accept or reject the proposal
  if (runif(0, 1) < acceptance_ratio){
    x_current = x_proposed  # Accept proposal
  accepted_samples <- accepted_samples + 1  # Increment acceptance counter
  
  # Store sample only after the burn-in period and at every 'thinning' interval
  } else{
  if (i >= burn_in && (i - burn_in) %% thinning == 0){
    samples.append(x_current)
  }
  }
  # Only increment i if the calculation was valid
  i <- i + 1  # Advance iteration
}  
  # Compute acceptance rate
  acceptance_rate = accepted_samples / total_iterations
  
  return(samples, acceptance_rate)
}

#----2. Instrtumental Distribution ----

# Proposal: Product of two log-noraml distributions
proposal_sampler <- function(y){
  mu1 <- log(y[1]); sigma1 <- 0.1 # Parameters of the Log-Normal for Alpha
  mu2 <- log(y[2]); sigma2 <- 0.1 # Parameters of the Log-Normal for Beta
  alpha_proposed <- rlnorm(mean=mu1, sigma=sigma1)
  beta_proposed <- rlnorm(mean=mu2, sigma=sigma2)
  return(alpha_proposed, beta_proposed)
}

# Probability density function of the proposal distribution (Product of Log-Normal distributions)
# Given vector y, conditioned on vector z.

proposal_pdf <- function(y, z){
  mu1 <- log(z[1]); sigma1 <- 0.1 # Parameters for Alpha
  mu2 <- log(z[2]); sigma2 <- 0.1 # Parameters for Beta
  pdf_prop <- dlnorm(mean=y[1], sigma=sigma1, exp(mu1))*dlnorm(mean=y[2], sigma=sigma2, exp(mu2))
  return(pdf_prop) 
}


# Datos
xn <- data1 %>% select(Defect) %>% data.matrix() %>% array(c(nrow(.), 1, 1, ncol(.)))
# class(xn)
# Estimación de los parámetros usando momentos: Se utilizará como semilla del método MCMC.
alpha0 = mean(xn)*(mean(xn)*(1-mean(xn))/var(xn)-1)
beta0 = (1-mean(xn))*(mean(xn)*(1-mean(xn))/var(xn)-1)

#---- Function for Monitoring MCMC Convergence ----

monitor_convergence <- function(x, y){
  # Plot the trace of the samples
  #.plot1 <- # Asumiendo que 'x' e 'y' son vectores de datos ya definidos
    
    par(mfrow = c(2, 1), mar = c(4, 4, 2, 1) + 0.1) # Configura el panel de gráficos 2x1 y los márgenes
  
  plot(x, type = "l", col = "blue", alpha = 0.7,
       main = "Trace plot of the chain for Alpha",
       xlab = "", ylab = "") # Grafica la primera serie
  
  plot(y, type = "l", col = "red", alpha = 0.7,
       main = "Trace plot of the chain for Beta",
       xlab = "Iterations", ylab = "") # Grafica la segunda serie con la etiqueta del eje x

  # Asumiendo que 'x' e 'y' son vectores de datos ya definidos
  
  plot_alpha_hist <- ggplot(data.frame(value = x), aes(x = value)) +
    geom_histogram(bins = 50, fill = "blue", alpha = 0.7, aes(y = ..density..)) +
    geom_density(color = "darkblue", linewidth = 1) +
    labs(title = "Histogram of Alpha", x = "Value", y = "Density") +
    theme_minimal()
  
  plot_beta_hist <- ggplot(data.frame(value = y), aes(x = value)) +
    geom_histogram(bins = 50, fill = "red", alpha = 0.7, aes(y = ..density..)) +
    geom_density(color = "darkred", linewidth = 1) +
    labs(title = "Histogram of Beta", x = "Value", y = "Density") +
    theme_minimal()
  
  grid.arrange(plot_alpha_hist, plot_beta_hist, ncol = 2)
  
  # Asumiendo que 'x' e 'y' son vectores de datos ya definidos
  
  par(mfrow = c(1, 2), mar = c(4, 4, 2, 1) + 0.1) # Configura el panel de gráficos 1x2 y los márgenes
  
  acf(x, lag.max = 50, main = "Autocorrelation of Alpha", xlab = "Lag", ylab = "ACF")
  acf(y, lag.max = 50, main = "Autocorrelation of Beta", xlab = "Lag", ylab = "ACF")
  
}

#---- 4. Posterior distribution using Gammas prior Product

post_prior_gammas <- function(x, alpha, beta, a, b, c, d){
  
}