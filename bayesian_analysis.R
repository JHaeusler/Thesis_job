#==============================================================================#
#                   Script Consolidado de Inferencia Bayesiana MCMC          #
#==============================================================================#

# --- 0. Carga de Paquetes Necesarios ---
# Instala los paquetes si no los tienes:
# install.packages(c("readxl", "ggplot2", "readr", "tidyverse", "dplyr",
#                    "stats", "boot", "coda", "matrixStats", "gridExtra"))

# 1. Carga e instalación de paquetes ----
lista_de_paquetes <- c("readxl", "ggplot2", "readr", "tidyverse", "dplyr",
                       "stats", "boot", "coda", "matrixStats", "gridExtra")
# Simplificado a los esenciales para el mapa de calor

for (paquete in lista_de_paquetes) {
  if(!require(paquete, character.only = TRUE)){
    install.packages(paquete)
    require(paquete, character.only = TRUE)
  }
}

#==============================================================================#
# --- 1. Definición de Datos y Cálculo de Parámetros Iniciales (Semilla) ---
#==============================================================================#
setwd("Thesis_job")
# Configura el directorio de trabajo (asegúrate de que sea la ruta correcta)
# setwd("Thesis_job")

# Carga de datos
# Los datos 'Defect' se asumen como proporciones o valores de una distribución Beta.
# Asegúrate de que tu archivo "Acceptance Sampling MIL-STD 105E for Quality Control.xlsx"
# esté en el directorio de trabajo o proporciona la ruta completa.
data1_excel <- read_excel("Acceptance Sampling MIL-STD 105E for Quality Control.xlsx", sheet = "tab")
data2_excel <- read_excel("Jurnal Rekavasi.xlsx", sheet = "tabla") # No usado en este script

# Extrae la columna 'Defect' como un vector numérico.
xn <- data1_excel %>% dplyr::select(Defect) %>% pull()

# Calcula los parámetros iniciales (alpha0, beta0) usando el método de los momentos
# para la distribución Beta. Estos servirán como el punto de partida (semilla)
# para la cadena de Markov de MCMC.
mean_xn <- mean(xn)
var_xn <- var(xn)

alpha0 <- mean_xn * (mean_xn * (1 - mean_xn) / var_xn - 1)
beta0 <- (1 - mean_xn) * (mean_xn * (1 - mean_xn) / var_xn - 1)

# Valida que los parámetros iniciales sean positivos para evitar errores de logaritmo.
# Si son no positivos, se ajustan a un valor pequeño y se emite una advertencia.
if (alpha0 <= 0 || beta0 <= 0) {
  warning("Los valores iniciales calculados (alpha0 o beta0) son no positivos. Ajustando a un valor mínimo.")
  alpha0 <- max(alpha0, 0.001)
  beta0 <- max(beta0, 0.001)
}

# Define el punto inicial de la cadena MCMC (x0)
x0 <- c(alpha0, beta0)

message(paste("Valores iniciales (alpha0, beta0): (", round(alpha0, 4), ", ", round(beta0, 4), ")"))

#==============================================================================#
# --- 2. Definición de Funciones de Propuesta (Instrumental Distribution) ---
#==============================================================================#

# proposal_sampler: Genera una propuesta para el siguiente estado de la cadena.
# Utiliza un producto de dos distribuciones Log-Normal.
proposal_sampler <- function(y) {
  sigma1 <- 0.1 # Desviación estándar del logaritmo para Alpha
  sigma2 <- 0.1 # Desviación estándar del logaritmo para Beta
  
  # Calcula la media del logaritmo basada en el valor actual (y).
  # Se incluye un pequeño valor para manejar casos donde y[1] o y[2] podrían ser no positivos.
  mu1 <- if (y[1] > 0) log(y[1]) else log(0.0001)
  mu2 <- if (y[2] > 0) log(y[2]) else log(0.0001)
  
  # Genera una nueva propuesta de Alpha y Beta usando la distribución Log-Normal.
  alpha_proposed <- rlnorm(1, meanlog = mu1, sdlog = sigma1)
  beta_proposed <- rlnorm(1, meanlog = mu2, sdlog = sigma2)
  
  # Retorna la propuesta como un vector.
  return(c(alpha_proposed, beta_proposed))
}

# proposal_log_pdf: Calcula el logaritmo de la densidad de probabilidad de la propuesta.
# Crucial para la estabilidad numérica en el algoritmo Metropolis-Hastings.
proposal_log_pdf <- function(y, z) {
  sigma1 <- 0.1
  sigma2 <- 0.1
  
  # Calcula la media del logaritmo basada en el estado actual (z).
  mu1 <- if (z[1] > 0) log(z[1]) else log(0.0001)
  mu2 <- if (z[2] > 0) log(z[2]) else log(0.0001)
  
  # Suma de las log-densidades de las dos Log-Normal (log del producto).
  log_pdf_prop <- dlnorm(y[1], meanlog = mu1, sdlog = sigma1, log = TRUE) +
    dlnorm(y[2], meanlog = mu2, sdlog = sigma2, log = TRUE)
  
  return(log_pdf_prop)
}

#==============================================================================#
# --- 3. Definición de Funciones de Densidad Posterior (Target) ---
#==============================================================================#

# target_log_pdf: Log-densidad posterior con priors Gamma para Alpha y Beta.
# Esta es la función objetivo que el MCMC intentará muestrear.
target_log_pdf <- function(current_params, x_data, prior_hyperparams) {
  alpha <- current_params[1]
  beta <- current_params[2]
  a <- prior_hyperparams[1] # shape para Alpha
  b <- prior_hyperparams[2] # rate para Alpha
  c <- prior_hyperparams[3] # shape para Beta
  d <- prior_hyperparams[4] # rate para Beta
  
  # Penalización si los parámetros no son válidos (Alpha o Beta <= 0).
  if (alpha <= 0 || beta <= 0) {
    return(-Inf)
  }
  
  n <- length(x_data)
  
  # Log-verosimilitud para datos que siguen una distribución Beta.
  # 'lgamma' es el logaritmo de la función Gamma.
  log_beta_ratio <- n * (lgamma(alpha + beta) - lgamma(alpha) - lgamma(beta))
  log_likelihood_data_terms <- (alpha - 1) * sum(log(x_data)) + (beta - 1) * sum(log(1 - x_data))
  log_likelihood <- log_beta_ratio + log_likelihood_data_terms
  
  # Log-priors Gamma para Alpha y Beta.
  log_prior_alpha <- dgamma(alpha, shape = a, rate = b, log = TRUE)
  log_prior_beta <- dgamma(beta, shape = c, rate = d, log = TRUE)
  log_prior <- log_prior_alpha + log_prior_beta
  
  # Log-densidad posterior total (verosimilitud + priors).
  log_posterior_density <- log_likelihood + log_prior
  return(log_posterior_density)
}

# target_log_pdf_jeffrey: Log-densidad posterior con Prior de Jeffreys.
# Utiliza la prior de Jeffreys para los parámetros de la distribución Beta.
target_log_pdf_jeffrey <- function(current_params, x_data, prior_hyperparams = NULL) {
  alpha <- current_params[1]
  beta <- current_params[2]
  
  # Penalización si los parámetros no son válidos.
  if (alpha <= 0 || beta <= 0) {
    return(-Inf)
  }
  
  # Cálculo de la Prior de Jeffreys usando la función Polygamma.
  # psigamma(x, deriv=1) es la primera derivada de la función digamma (trigamma).
  psi_alpha <- psigamma(alpha, deriv = 1)
  psi_beta <- psigamma(beta, deriv = 1)
  psi_alpha_beta <- psigamma(alpha + beta, deriv = 1)
  prior_term_inside_sqrt <- psi_alpha * psi_beta - (psi_alpha + psi_beta) * psi_alpha_beta
  
  # Penaliza si el término dentro de la raíz cuadrada es no positivo.
  if (prior_term_inside_sqrt <= 0) {
    return(-Inf)
  }
  log_prior_jeffrey <- 0.5 * log(prior_term_inside_sqrt)
  
  # Log-verosimilitud (misma que en target_log_pdf).
  n <- length(x_data)
  log_beta_ratio <- n * (lgamma(alpha + beta) - lgamma(alpha) - lgamma(beta))
  log_likelihood_data_terms <- (alpha - 1) * sum(log(x_data)) + (beta - 1) * sum(log(1 - x_data))
  log_likelihood <- log_beta_ratio + log_likelihood_data_terms
  
  # Log-densidad posterior total (verosimilitud + prior de Jeffreys).
  log_posterior_density <- log_likelihood + log_prior_jeffrey
  return(log_posterior_density)
}

#==============================================================================#
# --- 4. Algoritmo Metropolis-Hastings ---
#==============================================================================#

# Función principal para ejecutar la cadena de Markov.
metropolis_hastings <- function(target_log_pdf, proposal_sampler, proposal_log_pdf,
                                x0, n_samples, burn_in=1000, thinning=1,
                                x_data, prior_hyperparams = NULL) { # prior_hyperparams puede ser NULL para Jeffreys
  
  samples <- matrix(NA, nrow = n_samples, ncol = length(x0))
  x_current <- x0
  total_iterations <- burn_in + n_samples * thinning
  accepted_samples <- 0
  i <- 0
  
  while (i < total_iterations) {
    x_proposed <- proposal_sampler(x_current)
    
    log_target_current <- target_log_pdf(x_current, x_data, prior_hyperparams)
    log_target_proposed <- target_log_pdf(x_proposed, x_data, prior_hyperparams)
    
    log_proposal_proposed_to_current <- proposal_log_pdf(x_current, x_proposed)
    log_proposal_current_to_proposed <- proposal_log_pdf(x_proposed, x_current)
    
    log_acceptance_ratio <- log_target_proposed + log_proposal_proposed_to_current -
      (log_target_current + log_proposal_current_to_proposed)
    
    if (is.nan(log_acceptance_ratio) || is.infinite(log_acceptance_ratio)) {
      log_acceptance_ratio <- -Inf
    }
    
    if (log(runif(1)) < min(0, log_acceptance_ratio)) {
      x_current <- x_proposed
      accepted_samples <- accepted_samples + 1
    }
    
    if (i >= burn_in && (i - burn_in) %% thinning == 0) {
      sample_index <- (i - burn_in) / thinning + 1
      if (sample_index <= n_samples) {
        samples[sample_index, ] <- x_current
      }
    }
    i <- i + 1
  }
  
  acceptance_rate <- accepted_samples / total_iterations
  return(list(samples = samples, acceptance_rate = acceptance_rate))
}

#==============================================================================#
# --- 5. Funciones de Diagnóstico y Métricas ---
#==============================================================================#

# Measure_Diagnostic: Calcula varias métricas de diagnóstico para las cadenas MCMC.
Measure_Diagnostic <- function(data1, data2, burnin=0, thin=1, digits=5, cred_level=0.95, batch_size=100) {
  
  # Sub-función para el error estándar de la varianza por medias de bloques.
  batch_means_stderr_var_R <- function(chain, batch_size) {
    n <- length(chain)
    num_batches <- floor(n / batch_size)
    if (num_batches < 2) {
      warning("Se requieren al menos 2 bloques para estimar el error estándar de la varianza.")
      return(NA)
    }
    batch_vars <- unlist(lapply(0:(num_batches-1), function(i) {
      block_start <- i * batch_size + 1
      block_end <- (i + 1) * batch_size
      var(chain[block_start:block_end])
    }))
    return(sd(batch_vars) / sqrt(num_batches))
  }
  
  # Sub-función para el intervalo de credibilidad de la media.
  Cred_Interval_R <- function(x, level=0.95) {
    lower <- quantile(x, (1 - level)/2)
    upper <- quantile(x, 1 - (1 - level)/2)
    return(c(lower, upper))
  }
  
  # Sub-función para el intervalo de credibilidad de la varianza (bootstrap).
  Cred_Interval_V_R <- function(x, level=0.95) {
    n_bootstrap <- 10000
    variance_stat <- function(data, indices) {
      var(data[indices])
    }
    boot_obj <- boot(data = x, statistic = variance_stat, R = n_bootstrap)
    ci_results <- boot.ci(boot_obj, type = "perc", conf = level)$percent[4:5]
    return(round(ci_results, digits))
  }
  
  # Sub-función para el error estándar de la covarianza por medias de bloques.
  batch_means_stderr_cov_R <- function(chain, batch_size=100) {
    n <- nrow(chain)
    num_batches <- floor(n / batch_size)
    if (num_batches < 2) {
      warning("Se requieren al menos 2 bloques para estimar el error estándar de la covarianza.")
      return(NA)
    }
    batch_covs <- unlist(lapply(0:(num_batches-1), function(i) {
      block_start <- i * batch_size + 1
      block_end <- (i + 1) * batch_size
      cov(chain[block_start:block_end, 1], chain[block_start:block_end, 2])
    }))
    return(sd(batch_covs) / sqrt(num_batches))
  }
  
  # Sub-función para el intervalo de credibilidad de la covarianza (bootstrap).
  Cred_Interval_Cov_R <- function(x, level=0.95) {
    n_bootstrap <- 10000
    covariance_stat <- function(data, indices) {
      cov_matrix <- cov(data[indices, ])
      return(cov_matrix[1, 2])
    }
    boot_obj <- boot(data = x, statistic = covariance_stat, R = n_bootstrap, sim = "ordinary")
    ci_results <- boot.ci(boot_obj, type = "perc", conf = level)$percent[4:5]
    return(round(ci_results, digits))
  }
  
  # Aplica burn-in y thinning a los datos de entrada (si no se hizo ya)
  start_index <- burnin + 1
  seq_indices <- seq(from = start_index, to = length(data1), by = thin)
  new_data1 <- data1[seq_indices]
  new_data2 <- data2[seq_indices]
  
  mcmc_data1 <- as.mcmc(new_data1)
  mcmc_data2 <- as.mcmc(new_data2)
  
  ess_x1 <- effectiveSize(mcmc_data1)
  ess_x2 <- effectiveSize(mcmc_data2)
  
  CIM1 <- Cred_Interval_R(new_data1, cred_level)
  CIM2 <- Cred_Interval_R(new_data2, cred_level)
  CIV1 <- Cred_Interval_V_R(new_data1, cred_level)
  CIV2 <- Cred_Interval_V_R(new_data2, cred_level)
  
  combined_data <- cbind(new_data1, new_data2)
  CIcov <- Cred_Interval_Cov_R(combined_data, cred_level)
  cov_val <- cov(new_data1, new_data2)
  
  new_names <- c("Mean_X1", "Var_X1", "ESS_X1", "STDERR_Mean_X1", "STDERR_Var_X1",
                 "CI_Mean1_Lower", "CI_Mean1_Upper", "CI_Var1_Lower", "CI_Var1_Upper",
                 "Mean_X2", "Var_X2", "ESS_X2", "STDERR_Mean_X2", "STDERR_Var_X2",
                 "CI_Mean2_Lower", "CI_Mean2_Upper", "CI_Var2_Lower", "CI_Var2_Upper",
                 "Cov", "STDERR_Cov", "CI_Cov_Lower", "CI_Cov_Upper", "Length")
  
  results_vector <- c(
    mean(new_data1), var(new_data1), ess_x1, sd(new_data1) / sqrt(ess_x1),
    batch_means_stderr_var_R(new_data1, batch_size), CIM1[1], CIM1[2], CIV1[1], CIV1[2],
    mean(new_data2), var(new_data2), ess_x2, sd(new_data2) / sqrt(ess_x2),
    batch_means_stderr_var_R(new_data2, batch_size), CIM2[1], CIM2[2], CIV2[1], CIV2[2],
    cov_val, batch_means_stderr_cov_R(combined_data, batch_size), CIcov[1], CIcov[2],
    length(new_data1)
  )
  
  Numerical_results_vertical <- data.frame(
    Measure = new_names,
    Value = round(results_vector, digits)
  )
  
  return(list("Numerical" = Numerical_results_vertical))
}

#==============================================================================#
# --- 6. Monitoreo de Convergencia (Visualizaciones) ---
#==============================================================================#

# monitor_convergence: Genera gráficos para diagnosticar la convergencia de la cadena.
monitor_convergence <- function(x, y) {
  # Gráficos de Traza (Trace Plots)
  plot_trace_alpha <- ggplot(data.frame(Iteration = seq_along(x), Alpha = x), aes(x = Iteration, y = Alpha)) +
    geom_line(color = "blue", alpha = 0.7) +
    labs(title = "Gráfico de Traza de la Cadena para Alpha", x = "Iteraciones", y = "Valor de Alpha") +
    theme_minimal()
  plot_trace_beta <- ggplot(data.frame(Iteration = seq_along(y), Beta = y), aes(x = Iteration, y = Beta)) +
    geom_line(color = "red", alpha = 0.7) +
    labs(title = "Gráfico de Traza de la Cadena para Beta", x = "Iteraciones", y = "Valor de Beta") +
    theme_minimal()
  grid.arrange(plot_trace_alpha, plot_trace_beta, ncol = 1)
  
  # Histogramas de las Distribuciones Marginales
  plot_hist_alpha <- ggplot(data.frame(Value = x), aes(x = Value)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "blue", alpha = 0.7) +
    geom_density(color = "darkblue", linewidth = 1) +
    labs(title = "Histograma de Alpha", x = "Valor", y = "Densidad") +
    theme_minimal()
  plot_hist_beta <- ggplot(data.frame(Value = y), aes(x = Value)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "red", alpha = 0.7) +
    geom_density(color = "darkred", linewidth = 1) +
    labs(title = "Histograma de Beta", x = "Valor", y = "Densidad") +
    theme_minimal()
  grid.arrange(plot_hist_alpha, plot_hist_beta, ncol = 2)
  
  # Diagrama de Dispersión (Scatter Plot)
  plot_scatter <- ggplot(data.frame(Alpha = x, Beta = y), aes(x = Alpha, y = Beta)) +
    geom_point(alpha = 0.3, color = "darkgreen") +
    labs(title = "Diagrama de Dispersión de las Muestras (Alpha, Beta)", x = "Alpha", y = "Beta") +
    theme_minimal() + theme(aspect.ratio = 1)
  print(plot_scatter)
  
  # Gráficos de Autocorrelación (ACF)
  old_par <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1) + 0.1)
  acf(x, lag.max = 50, main = "Autocorrelación de Alpha", xlab = "Retardo (Lag)", ylab = "ACF")
  acf(y, lag.max = 50, main = "Autocorrelación de Beta", xlab = "Retardo (Lag)", ylab = "ACF")
  par(old_par)
}

#==============================================================================#
# --- 7. Ejecución de la Simulación MCMC y Análisis con Prior Gamma ---
#==============================================================================#

message("\n--- Ejecutando MCMC con Priors Gamma ---")

# Hiperparámetros de los priors Gamma
prior_hyperparams_gammas <- c(a_shape = 0.1, b_rate = 0.1, c_shape = 0.1, d_rate = 0.1)

# Configuración de la simulación
n_samples_mcmc <- 20000
burn_in_mcmc <- 5000
thinning_mcmc <- 5

# Ejecuta el MCMC
results_mcmc_gammas <- metropolis_hastings(
  target_log_pdf = target_log_pdf, # Usa la función objetivo con prior Gamma
  proposal_sampler = proposal_sampler,
  proposal_log_pdf = proposal_log_pdf,
  x0 = x0,
  n_samples = n_samples_mcmc,
  burn_in = burn_in_mcmc,
  thinning = thinning_mcmc,
  x_data = xn,
  prior_hyperparams = prior_hyperparams_gammas
)

# Extrae las muestras y la tasa de aceptación
samples_gammas <- results_mcmc_gammas$samples
acceptance_rate_gammas <- results_mcmc_gammas$acceptance_rate

# Imprime la tasa de aceptación
message(paste("Tasa de Aceptación (Gammas):", round(acceptance_rate_gammas, 6)))
message(paste0("Media (alpha, beta) (Gammas): (", round(colMeans(samples_gammas)[1], 6), ", ", round(colMeans(samples_gammas)[2], 6), ")"))
message(paste0("Varianza (alpha, beta) (Gammas): (", round(apply(samples_gammas, 2, var)[1], 6), ", ", round(apply(samples_gammas, 2, var)[2], 6), ")"))

# Separa las muestras para el monitoreo y diagnóstico
alpha_samples_gammas <- samples_gammas[, 1]
beta_samples_gammas <- samples_gammas[, 2]

# Monitorea la convergencia visualmente
message("\nMostrando gráficos de convergencia para Priors Gamma...")
monitor_convergence(alpha_samples_gammas, beta_samples_gammas)

# Calcula y muestra las métricas de diagnóstico
results_gammas <- Measure_Diagnostic(
  data1 = alpha_samples_gammas,
  data2 = beta_samples_gammas,
  burnin = 0, thin = 1, digits = 5, cred_level = 0.95, batch_size = 100
)
message("\nResultados de Diagnóstico (Priors Gamma):")
print(results_gammas$Numerical)

#==============================================================================#
# --- 8. Ejecución de la Simulación MCMC y Análisis con Prior de Jeffreys ---
#==============================================================================#

message("\n--- Ejecutando MCMC con Prior de Jeffreys ---")

# Ejecuta el MCMC con la prior de Jeffreys
results_mcmc_jeffreys <- metropolis_hastings(
  target_log_pdf = target_log_pdf_jeffrey, # Usa la función objetivo con prior de Jeffreys
  proposal_sampler = proposal_sampler,
  proposal_log_pdf = proposal_log_pdf,
  x0 = x0, # Usa los mismos puntos iniciales
  n_samples = n_samples_mcmc,
  burn_in = burn_in_mcmc,
  thinning = thinning_mcmc,
  x_data = xn,
  prior_hyperparams = NULL # No se usan hiperparámetros para la prior de Jeffreys
)

# Extrae las muestras y la tasa de aceptación
samples_jeffreys <- results_mcmc_jeffreys$samples
acceptance_rate_jeffreys <- results_mcmc_jeffreys$acceptance_rate

# Imprime la tasa de aceptación
message(paste("Tasa de Aceptación (Jeffreys):", round(acceptance_rate_jeffreys, 6)))
message(paste0("Media (alpha, beta) (Jeffreys): (", round(colMeans(samples_jeffreys)[1], 6), ", ", round(colMeans(samples_jeffreys)[2], 6), ")"))
message(paste0("Varianza (alpha, beta) (Jeffreys): (", round(apply(samples_jeffreys, 2, var)[1], 6), ", ", round(apply(samples_jeffreys, 2, var)[2], 6), ")"))

# Separa las muestras para el monitoreo y diagnóstico
alpha_samples_jeffreys <- samples_jeffreys[, 1]
beta_samples_jeffreys <- samples_jeffreys[, 2]

# Monitorea la convergencia visualmente
message("\nMostrando gráficos de convergencia para Prior de Jeffreys...")
monitor_convergence(alpha_samples_jeffreys, beta_samples_jeffreys)

# Calcula y muestra las métricas de diagnóstico
results_jeffreys <- Measure_Diagnostic(
  data1 = alpha_samples_jeffreys,
  data2 = beta_samples_jeffreys,
  burnin = 0, thin = 1, digits = 5, cred_level = 0.95, batch_size = 100
)
message("\nResultados de Diagnóstico (Prior de Jeffreys):")
print(results_jeffreys$Numerical)

#==============================================================================#
# --- 9. Consolidación y Comparación de Resultados ---
#==============================================================================#

message("\n--- Consolidando Resultados para Comparación ---")

# Se crea un conjunto de resultados 'new' para demostrar la combinación,
# por ejemplo, podría ser otra ejecución o simplemente una copia de uno de los anteriores.
# Para este script consolidado, lo crearemos como una copia de los resultados de Jeffreys
# para que el bloque de consolidación funcione sin problemas.
# En un escenario real, 'results_new' provendría de una tercera ejecución MCMC.
alpha_samples_new <- alpha_samples_jeffreys # Placeholder
beta_samples_new <- beta_samples_jeffreys   # Placeholder
results_new <- Measure_Diagnostic(
  data1 = alpha_samples_new,
  data2 = beta_samples_new,
  burnin = 0, thin = 1, digits = 5, cred_level = 0.95, batch_size = 100
)

# Extraer el data frame 'Numerical' de cada objeto de resultados
results_gammas_df <- results_gammas$Numerical
results_jeffreys_df <- results_jeffreys$Numerical
results_new_df <- results_new$Numerical

# Asegurar que los DataFrames estén ordenados por 'Measure' para la combinación
results_gammas_df <- results_gammas_df %>% dplyr::arrange(Measure)
results_jeffreys_df <- results_jeffreys_df %>% dplyr::arrange(Measure)
results_new_df <- results_new_df %>% dplyr::arrange(Measure)

# Combinar las tablas en un único data frame comparativo
combined_results <- results_gammas_df %>%
  dplyr::rename(Gammas = Value) %>%
  dplyr::mutate(Jeffreys = results_jeffreys_df$Value,
                New = results_new_df$Value) # Usamos mutate para añadir columnas

message("\nTabla Comparativa de Resultados Diagnósticos:")
print(combined_results)

#==============================================================================#
#                           Fin del Script                                   #
#==============================================================================#
