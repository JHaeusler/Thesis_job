#--- Bayesian Inference ----

setwd("D:/Github/Thesis_job")

lista_de_paquetes <- c("readxl", "ggplot2", "readr", "tidyr") # Reemplaza con tus paquetes

for (paquete in lista_de_paquetes) {
  if(!require(paquete, character.only = TRUE)){
    install.packages(paquete)
    require(paquete, character.only = TRUE)
  }
}

data1 <- read_excel("Acceptance Sampling MIL-STD 105E for Quality Control.xlsx",
                    sheet = "tab")
data2 <- read_excel("Jurnal Rekavasi.xlsx", sheet = "tabla")

#---- Metropolis-Hastings ----

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

#---- Instrtumental Distribution ----

# Proposal: Product of two log-noraml distributions
proposal_sampler <- function(y){
  mu1 <- log(y[1]); sigma1 <- 0.1 # Parameters of the Log-Normal for Alpha
  mu2 <- log(y[2]); sigma1 <- 0.1 # Parameters of the Log-Normal for Beta
  alpha_proposed <- rlognormal(mean=mu1, sigma=sigma1)
  beta_proposed <- rlognormal(mean=mu2, sigma=sigma2)
}


