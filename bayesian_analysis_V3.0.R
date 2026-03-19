# ============================================================================== #
# SCRIPT 1: ESTIMACIÓN BAYESIANA DE LA CALIDAD (OBJETIVO ESPECÍFICO 1)
# Metodología: Torres, Tovar y Bran (2025) vs Priors Clásicas
# Autor: Juan Sebastián Haeusler
# ============================================================================== #

# --- 1. Carga de Paquetes ----
paquetes <- c("readxl", "readr", "dplyr", "stats", "boot", "coda", "matrixStats", "AcceptanceSampling")

for (p in paquetes) {
  if(!require(p, character.only = TRUE)){
    install.packages(p)
    require(p, character.only = TRUE)
  }
}

# ============================================================================== #
# --- 2. Carga de Datos, Simulación y Visualización Empírica (R Base) ---
# ============================================================================== #
# setwd("~/Thesis_job") 

# Fuente 1: Defiatri & Damayanti (2023)
data1_excel <- read_excel("Acceptance Sampling MIL-STD 105E for Quality Control.xlsx", sheet = "tab")
vec_data1 <- data1_excel %>% select(Defect) %>% pull()

# Fuente 2: Isnanto et al. (2019)
data2_excel <- read_excel("Jurnal Rekavasi.xlsx", sheet = "tabla")
vec_data2 <- data2_excel %>% select(p) %>% pull()

# Fuente 3: Simulación Propia
set.seed(060722)
N <- 1000; AQL <- 0.05; LTPD <- 0.10; alpha_des <- 0.01; beta_des <- 0.05
plan_kier <- find.plan(PRP = c(AQL, 1 - alpha_des), CRP = c(LTPD, beta_des), N = N, type = "hypergeom")
n_sample <- plan_kier$n; c_critico <- plan_kier$c

lotes <- 35; sum_n <- numeric(lotes); decis <- numeric(lotes)
for (k in 1:lotes) {
  L <- ifelse(runif(N) <= AQL, 1, 0)
  sample_n <- sample(L, n_sample)
  sum_n[k] <- sum(sample_n)
  decis[k] <- ifelse(sum_n[k] > c_critico, 1, 0)
}
resultados <- data.frame(Lote = 1:lotes, Defectuosos = sum_n, Rechazado = decis) %>%
  mutate(p_est = Defectuosos / n_sample)
vec_data3 <- resultados %>% pull(p_est)

# --- GRÁFICO 1: DENSIDADES EMPÍRICAS COMPARADAS (R BASE) ----
# Sanitizamos ligeramente para el cálculo de densidad general
# d1 <- density(na.omit(vec_data1))
# d2 <- density(na.omit(vec_data2))
# d3 <- density(na.omit(vec_data3))

# Calcular límites para que todo encaje en el lienzo
# xlim_val <- range(c(d1$x, d2$x, d3$x))
# ylim_val <- c(0, max(c(d1$y, d2$y, d3$y)) * 1.1)
# 
# plot(NULL, xlim = xlim_val, ylim = ylim_val, bty = "l",
#      main = "Densidades Empíricas de la Proporción de Defectos (p)",
#      xlab = "Proporción de Defectuosos", ylab = "Densidad")
# 
# # Dibujar polígonos con transparencia (alpha)
# polygon(d1, col = rgb(0.1, 0.5, 0.8, alpha = 0.4), border = "dodgerblue", lwd = 2)
# polygon(d2, col = rgb(0.8, 0.2, 0.2, alpha = 0.4), border = "firebrick", lwd = 2)
# polygon(d3, col = rgb(0.2, 0.8, 0.2, alpha = 0.4), border = "forestgreen", lwd = 2)
# 
# legend("topright", legend = c("Data 1 (Defiatri)", "Data 2 (Isnanto)", "Data 3 (Simulación)"),
#        fill = c(rgb(0.1, 0.5, 0.8, alpha = 0.4), rgb(0.8, 0.2, 0.2, alpha = 0.4), rgb(0.2, 0.8, 0.2, alpha = 0.4)),
#        border = c("dodgerblue", "firebrick", "forestgreen"), bty = "n", cex = 0.9)
# 
# # ----------------------------------------------------------------------------
# SELECTOR DINÁMICO DE DATOS (1, 2 o 3) 
FUENTE_ACTIVA <- 1
# ----------------------------------------------------------------------------

if(FUENTE_ACTIVA == 1) { data_empirica <- vec_data1; nombre_fuente <- "Data 1: Defiatri & Damayanti"
} else if(FUENTE_ACTIVA == 2) { data_empirica <- vec_data2; nombre_fuente <- "Data 2: Isnanto et al."
} else { data_empirica <- vec_data3; nombre_fuente <- "Data 3: Simulación Propia" }

# Sanitización de Datos (Corrección log(0))
data_empirica <- as.numeric(na.omit(data_empirica))
data_empirica <- pmax(pmin(data_empirica, 1 - 1e-5), 1e-5)
message(paste("\n>>> FUENTE DE DATOS CARGADA:", nombre_fuente, "<<<\n"))

# ============================================================================== #
# --- 3. Funciones Base MCMC ---
# ============================================================================== #
get_initial_values <- function(data) {
  mean_xn <- mean(data); var_xn <- var(data)
  a_0 <- mean_xn * (mean_xn * (1 - mean_xn) / var_xn - 1)
  b_0 <- (1 - mean_xn) * (mean_xn * (1 - mean_xn) / var_xn - 1)
  if (a_0 <= 0 || b_0 <= 0) { a_0 <- 0.001; b_0 <- 0.001 }
  return(c(a_0, b_0))
}

proposal_sampler <- function(y) {
  sigma1 <- 0.1; sigma2 <- 0.1 
  mu1 <- if (y[1] > 0) log(y[1]) else log(0.0001)
  mu2 <- if (y[2] > 0) log(y[2]) else log(0.0001)
  return(c(rlnorm(1, mu1, sigma1), rlnorm(1, mu2, sigma2)))
}

proposal_log_pdf <- function(y, z) {
  sigma1 <- 0.1; sigma2 <- 0.1
  mu1 <- if (z[1] > 0) log(z[1]) else log(0.0001)
  mu2 <- if (z[2] > 0) log(z[2]) else log(0.0001)
  return(dlnorm(y[1], mu1, sigma1, log = TRUE) + dlnorm(y[2], mu2, sigma2, log = TRUE))
}

# ============================================================================== #
# --- 4. Definición de las 3 Funciones Objetivo (Priors) ---
# ============================================================================== #
target_log_pdf_gamma <- function(current_params, x_data, prior_hyperparams) {
  a <- current_params[1]; b <- current_params[2]
  if (a <= 0 || b <= 0) return(-Inf)
  n <- length(x_data)
  log_likelihood <- n * (lgamma(a + b) - lgamma(a) - lgamma(b)) + 
    (a - 1) * sum(log(x_data)) + (b - 1) * sum(log(1 - x_data))
  log_prior <- dgamma(a, prior_hyperparams[1], prior_hyperparams[2], log=TRUE) + 
    dgamma(b, prior_hyperparams[3], prior_hyperparams[4], log=TRUE)
  return(log_likelihood + log_prior)
}

target_log_pdf_jeffrey <- function(current_params, x_data, prior_hyperparams = NULL) {
  a <- current_params[1]; b <- current_params[2]
  if (a <= 0 || b <= 0) return(-Inf)
  psi_a <- psigamma(a, 1); psi_b <- psigamma(b, 1); psi_ab <- psigamma(a + b, 1)
  term <- psi_a * psi_b - (psi_a + psi_b) * psi_ab
  if (term <= 0) return(-Inf)
  n <- length(x_data)
  log_likelihood <- n * (lgamma(a + b) - lgamma(a) - lgamma(b)) + 
    (a - 1) * sum(log(x_data)) + (b - 1) * sum(log(1 - x_data))
  return(log_likelihood + 0.5 * log(term))
}

target_log_pdf_bivariada <- function(current_params, x_data, prior_hyperparams) {
  a <- current_params[1]; b <- current_params[2]
  if (a <= 0 || b <= 0) return(-Inf)
  n <- length(x_data)
  log_likelihood <- n * (lgamma(a + b) - lgamma(a) - lgamma(b)) + 
    (a - 1) * sum(log(x_data)) + (b - 1) * sum(log(1 - x_data))
  log_prior <- (prior_hyperparams[1] - 1) * log(a) + (prior_hyperparams[2] - 1) * log(b) - 
    (prior_hyperparams[1] + prior_hyperparams[2]) * log(a + b) - 
    (prior_hyperparams[3] + prior_hyperparams[4]) * log(a + b + 1)
  return(log_likelihood + log_prior)
}

# ============================================================================== #
# --- 5. Motor MCMC y Diagnósticos ---
# ============================================================================== #
metropolis_hastings <- function(target_log_pdf, proposal_sampler, proposal_log_pdf,
                                x0, n_samples, burn_in, thinning, x_data, prior_hyperparams=NULL) {
  samples <- matrix(NA, nrow = n_samples, ncol = 2)
  x_current <- x0
  accepted <- 0
  for (i in 1:(burn_in + n_samples * thinning)) {
    x_prop <- proposal_sampler(x_current)
    log_r <- target_log_pdf(x_prop, x_data, prior_hyperparams) + proposal_log_pdf(x_current, x_prop) -
      (target_log_pdf(x_current, x_data, prior_hyperparams) + proposal_log_pdf(x_prop, x_current))
    if (is.nan(log_r)) log_r <- -Inf
    if (log(runif(1)) < min(0, log_r)) { x_current <- x_prop; accepted <- accepted + 1 }
    if (i > burn_in && (i - burn_in) %% thinning == 0) samples[(i - burn_in) / thinning, ] <- x_current
  }
  return(list(samples = samples, acceptance_rate = accepted / (burn_in + n_samples * thinning)))
}

Measure_Diagnostic <- function(chain_a, chain_b) {
  mcmc_a <- as.mcmc(chain_a); mcmc_b <- as.mcmc(chain_b)
  ess_a <- effectiveSize(mcmc_a); ess_b <- effectiveSize(mcmc_b)
  c_int_a <- quantile(chain_a, c(0.025, 0.975)); c_int_b <- quantile(chain_b, c(0.025, 0.975))
  cov_ab <- cov(chain_a, chain_b); cor_ab <- cor(chain_a, chain_b)
  
  res <- c(mean(chain_a), var(chain_a), ess_a, c_int_a[1], c_int_a[2],
           mean(chain_b), var(chain_b), ess_b, c_int_b[1], c_int_b[2], cov_ab, cor_ab)
  names <- c("Mean_a", "Var_a", "ESS_a", "CI_a_Lower", "CI_a_Upper",
             "Mean_b", "Var_b", "ESS_b", "CI_b_Lower", "CI_b_Upper", "Covariance_ab", "Correlation_ab")
  return(data.frame(Measure = names, Value = round(res, 4)))
}

# monitor_convergence <- function(chain_a, chain_b, prior_name, dataset_name) {
#   par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
#   plot(chain_a, type="l", col="blue", main=paste("Traza a -", prior_name), ylab="a")
#   plot(chain_b, type="l", col="red", main=paste("Traza b -", prior_name), ylab="b")
#   acf(chain_a, main="ACF a")
#   acf(chain_b, main="ACF b")
#   mtext(dataset_name, outer = TRUE, cex = 1.2, font = 2)
#   par(mfrow = c(1, 1))
# }
# ============================================================================== #
# --- FUNCIÓN DE DIAGNÓSTICO GRÁFICO (Panel 4 Pantallas - Torres et al. 2025) ---
# ============================================================================== #

monitor_convergence <- function(chain_a, chain_b, prior_name, dataset_name) {
  
  # Sub-función para generar el panel de 4 cuadrantes para un parámetro específico
  plot_4_panel <- function(chain, param_name) {
    # Configurar el panel 2x2 con márgenes ajustados
    old_par <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
    
    # 1. Arriba Izquierda: Histograma y Densidad
    hist(chain, prob = TRUE, col = rgb(0.8, 0.8, 0.8, 0.5), border = "gray",
         main = "Histogram and density", xlab = param_name, ylab = "Density")
    lines(density(chain), col = "dodgerblue", lwd = 2)
    
    # 2. Arriba Derecha: Traza de la muestra
    plot(chain, type = "l", col = "gray40",
         main = paste("Trace of the random sample of size", length(chain)),
         xlab = "Iteration", ylab = param_name)
    
    # 3. Abajo Izquierda: Control de Convergencia (Media Ergódica)
    ergodic_mean <- cumsum(chain) / seq_along(chain)
    plot(ergodic_mean, type = "l", col = "firebrick", lwd = 2,
         main = "Convergence control using averaging",
         xlab = "Iteration", ylab = "Averaging")
    # Línea asintótica ideal (el valor final donde converge)
    abline(h = ergodic_mean[length(ergodic_mean)], col = "black", lty = 2, lwd = 1.5)
    
    # 4. Abajo Derecha: Función de Autocorrelación (ACF)
    acf(chain, main = "Autocorrelation function", ylab = "ACF", xlab = "Lag")
    
    # Título principal de la ventana (para saber a qué escenario pertenece)
    mtext(paste("Diagnósticos MCMC Parámetro '", param_name, "' | Prior:", prior_name, "-", dataset_name),
          outer = TRUE, cex = 1.1, font = 2)
    
    par(old_par) # Restaurar parámetros originales
  }
  
  # Al llamar a la función, RStudio generará dos "ventanas" (plots) consecutivas:
  # Puede navegar entre ellas usando las flechas de "Plots" en RStudio.
  
  # Generar Panel 4x4 para el parámetro 'a'
  plot_4_panel(chain_a, "a")
  
  # Generar Panel 4x4 para el parámetro 'b'
  plot_4_panel(chain_b, "b")
}
# ============================================================================== #
# --- 6. Ejecución del Análisis de Sensibilidad ---
# ============================================================================== #
theta_0 <- get_initial_values(data_empirica)
n_s <- 5000; b_in <- 1000; thin <- 2

message("Ejecutando MCMC: Prior Gamma...")
res_gamma <- metropolis_hastings(target_log_pdf_gamma, proposal_sampler, proposal_log_pdf, 
                                 theta_0, n_s, b_in, thin, data_empirica, c(0.1, 0.1, 0.1, 0.1))
diag_gamma <- Measure_Diagnostic(res_gamma$samples[,1], res_gamma$samples[,2])

message("Ejecutando MCMC: Prior Jeffreys...")
res_jeffrey <- metropolis_hastings(target_log_pdf_jeffrey, proposal_sampler, proposal_log_pdf, 
                                   theta_0, n_s, b_in, thin, data_empirica, NULL)
diag_jeffrey <- Measure_Diagnostic(res_jeffrey$samples[,1], res_jeffrey$samples[,2])

message("Ejecutando MCMC: Prior Bivariada (Torres et al.)...")
res_biv <- metropolis_hastings(target_log_pdf_bivariada, proposal_sampler, proposal_log_pdf, 
                               theta_0, n_s, b_in, thin, data_empirica, c(2, 2, 2, 2))
diag_biv <- Measure_Diagnostic(res_biv$samples[,1], res_biv$samples[,2])

# ============================================================================== #
# --- 7. Consolidación y Gráficos Posteriores (R BASE) ---
# ============================================================================== #
message("\n=== TABLA COMPARATIVA DE DIAGNÓSTICOS ===")
tabla_final <- diag_gamma %>% rename(Gamma = Value) %>%
  mutate(Jeffreys = diag_jeffrey$Value, Bivariada = diag_biv$Value)
print(tabla_final)

monitor_convergence(res_biv$samples[,1], res_biv$samples[,2], "Bivariada Torres", nombre_fuente)

# --- GRÁFICO 2: DENSIDADES POSTERIORES COMPARADAS (a y b) EN R BASE ---
par(mfrow = c(1, 2), oma = c(0, 0, 2, 0), mar = c(4, 4, 3, 1))

# Densidades para el parámetro 'a'
da_g <- density(res_gamma$samples[,1])
da_j <- density(res_jeffrey$samples[,1])
da_b <- density(res_biv$samples[,1])

plot(NULL, xlim = range(c(da_g$x, da_j$x, da_b$x)), ylim = c(0, max(c(da_g$y, da_j$y, da_b$y)) * 1.1),
     main = "Densidad Posterior 'a'", xlab = "Valor de a", ylab = "Densidad", bty = "l")
polygon(da_g, col = rgb(0.5, 0.5, 0.5, alpha = 0.4), border = "gray30", lwd = 2)
polygon(da_j, col = rgb(0.8, 0.2, 0.2, alpha = 0.4), border = "firebrick", lwd = 2)
polygon(da_b, col = rgb(0.1, 0.5, 0.8, alpha = 0.4), border = "dodgerblue", lwd = 2)

legend("topright", legend = c("Gamma", "Jeffreys", "Bivariada"),
       fill = c(rgb(0.5, 0.5, 0.5, alpha = 0.4), rgb(0.8, 0.2, 0.2, alpha = 0.4), rgb(0.1, 0.5, 0.8, alpha = 0.4)),
       border = c("gray30", "firebrick", "dodgerblue"), bty = "n", cex = 0.8)

# Densidades para el parámetro 'b'
db_g <- density(res_gamma$samples[,2])
db_j <- density(res_jeffrey$samples[,2])
db_b <- density(res_biv$samples[,2])

plot(NULL, xlim = range(c(db_g$x, db_j$x, db_b$x)), ylim = c(0, max(c(db_g$y, db_j$y, db_b$y)) * 1.1),
     main = "Densidad Posterior 'b'", xlab = "Valor de b", ylab = "Densidad", bty = "l")
polygon(db_g, col = rgb(0.5, 0.5, 0.5, alpha = 0.4), border = "gray30", lwd = 2)
polygon(db_j, col = rgb(0.8, 0.2, 0.2, alpha = 0.4), border = "firebrick", lwd = 2)
polygon(db_b, col = rgb(0.1, 0.5, 0.8, alpha = 0.4), border = "dodgerblue", lwd = 2)

mtext(paste("Comparación de Priors -", nombre_fuente), outer = TRUE, cex = 1.2, font = 2)
par(mfrow = c(1, 1)) # Restaurar panel

# --- EXTRACCIÓN DE PARÁMETROS PARA LA FASE II ---
a_est <- mean(res_biv$samples[, 1])
b_est <- mean(res_biv$samples[, 2])
message(sprintf("\n--- PARÁMETROS DEFINITIVOS FASE II: a = %.4f | b = %.4f ---", a_est, b_est))