# ==============================================================================
# --- 0. Configuración Inicial y Carga de Paquetes ---
# ==============================================================================

lista_de_paquetes <- c("readxl", "ggplot2", "readr", "tidyverse", "dplyr",
                       "stats", "boot", "coda", "matrixStats", "gridExtra", "pracma")

for (paquete in lista_de_paquetes) {
  if(!require(paquete, character.only = TRUE)){
    install.packages(paquete)
    require(paquete, character.only = TRUE)
  }
}

# La librería 'boot' y 'coda' ya están en la lista, pero las siguientes
# líneas aseguran que no haya errores si no se cargan explícitamente.
library(boot)
library(coda)

setwd("D:/Github/Thesis_job")

# ==============================================================================
# --- 1. Definición de Datos y Cálculo de Parámetros Iniciales (Semilla) ---
# ==============================================================================



### primera densidad

### método de Tovar

theta_exp <- sort(runif(2))




data <- read_excel("Acceptance Sampling MIL-STD 105E for Quality Control.xlsx")
xn <- data$Defect

# Evitar valores de 0 o 1, que causan problemas en el logaritmo.
epsilon <- .Machine$double.eps
xn[xn == 0] <- epsilon
xn[xn == 1] <- 1 - epsilon

# Calcula los parámetros iniciales (alpha0, beta0) usando el método de los momentos
mean_xn <- mean(xn)
var_xn <- var(xn)
denominator <- (mean_xn * (1 - mean_xn) / var_xn) - 1
if (is.infinite(denominator) || is.nan(denominator) || denominator <= 0) {
  warning("El cálculo de los parámetros iniciales para la distribución Beta es inestable. Ajustando a valores predeterminados.")
  alpha0 <- 1 # Valores predeterminados seguros
  beta0 <- 1
} else {
  alpha0 <- mean_xn * denominator
  beta0 <- (1 - mean_xn) * denominator
}

# Valida que los parámetros iniciales sean positivos.
if (alpha0 <= 0 || beta0 <= 0) {
  warning("Los valores iniciales calculados (alpha0 o beta0) son no positivos. Ajustando a un valor mínimo.")
  alpha0 <- max(alpha0, 0.001)
  beta0 <- max(beta0, 0.001)
}

# Define el punto inicial de la cadena MCMC (x0)
x0 <- c(alpha0, beta0)

# ==============================================================================
# --- 2. Definición de Funciones de Propuesta (Instrumental Distribution) ---
# ==============================================================================

# Utiliza un producto de dos distribuciones Log-Normal.
proposal_sampler <- function(y) {
  sigma1 <- 0.1 # Desviación estándar del logaritmo para Alpha
  sigma2 <- 0.1 # Desviación estándar del logaritmo para Beta
  
  mu1 <- if (y[1] > 0) log(y[1]) else log(0.0001)
  mu2 <- if (y[2] > 0) log(y[2]) else log(0.0001)
  
  # Genera una nueva propuesta de Alpha y Beta usando la distribución Log-Normal.
  alpha_proposed <- rlnorm(1, meanlog = mu1, sdlog = sigma1)
  beta_proposed <- rlnorm(1, meanlog = mu2, sdlog = sigma2)
  
  # Retorna la propuesta como un vector.
  return(c(alpha_proposed, beta_proposed))
}

# Versión corregida de la log-densidad de la función de propuesta.
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

# ==============================================================================
# --- 3. A prioris consideradas en el estudio ---
# ==============================================================================

# Log-densidad posterior con priors Gamma para Alpha y Beta.
post_prior_gammas <- function(x_data, current_params, prior_hyperparams) {
  alpha <- current_params[1]
  beta <- current_params[2]
  a <- prior_hyperparams[1] # shape para Alpha
  b <- prior_hyperparams[2] # rate para Alpha
  c <- prior_hyperparams[3] # shape para Beta
  d <- prior_hyperparams[4] # rate para Beta
  
  if (alpha <= 0 || beta <= 0) {
    return(-Inf)
  }
  
  n <- length(x_data)
  
  log_beta_ratio <- n * (lgamma(alpha + beta) - lgamma(alpha) - lgamma(beta))
  log_likelihood <- (alpha - 1) * sum(log(x_data)) + (beta - 1) * sum(log(1 - x_data))
  log_prior <- (a - 1) * log(alpha) + (c - 1) * log(beta) - b * alpha - d * beta
  log_post <- log_beta_ratio + log_likelihood + log_prior
  return(log_post)
}

# Versión corregida para la prior de Jeffreys.
post_prior_Jeffrey <- function(x_data, current_params, prior_hyperparams) {
  alpha <- current_params[1]
  beta <- current_params[2]
  
  if (alpha <= 0 || beta <= 0) {
    return(-Inf)
  }
  
  # Cálculo de la Prior de Jeffreys.
  psi_alpha <- psigamma(alpha, deriv = 1)
  psi_beta <- psigamma(beta, deriv = 1)
  psi_alpha_beta <- psigamma(alpha + beta, deriv = 1)
  
  # Se toma el logaritmo de la prior de Jeffreys.
  # Se usa pmax para evitar log(número negativo) si el cálculo de la prior resulta en un valor no positivo.
  prior_val <- psi_alpha * psi_beta - (psi_alpha + psi_beta) * psi_alpha_beta
  if (prior_val <= 0) {
    log_prior <- -Inf
  } else {
    log_prior <- 0.5 * log(prior_val)
  }
  
  n <- length(x_data)
  log_beta_ratio <- n * (lgamma(alpha + beta) - lgamma(alpha) - lgamma(beta))
  log_exp_term <- (alpha - 1) * sum(log(x_data)) + (beta - 1) * sum(log(1 - x_data))
  
  log_post <- log_beta_ratio + log_exp_term + log_prior
  return(log_post)
}

# La función post_new_prior se mantiene como estaba.
post_new_prior <- function(x_data, current_params, prior_hyperparams) {
  alpha <- current_params[1]
  beta <- current_params[2]
  
  if (alpha <= 0 || beta <= 0) {
    return(-Inf)
  }
  
  a <- prior_hyperparams[1]
  b <- prior_hyperparams[2]
  c <- prior_hyperparams[3]
  d <- prior_hyperparams[4]
  
  log_prior <- (a - 1) * log(alpha) + (b - 1) * log(beta) +
    (d - (a + b)) * log(alpha + beta) - (c + d) * log(alpha + beta + 1)
  
  n <- length(x_data)
  log_beta_ratio <- n * (lgamma(alpha + beta) - lgamma(alpha) - lgamma(beta))
  log_exp_term <- (alpha - 1) * sum(log(x_data)) + (beta - 1) * sum(log(1 - x_data))
  
  log_post <- log_beta_ratio + log_exp_term + log_prior
  return(log_post)
}

# ==============================================================================
# --- 4. Algoritmo Metropolis-Hastings Mejorado ---
# ==============================================================================

metropolis_hastings <- function(target_log_pdf, proposal_sampler, proposal_log_pdf,
                                x0, n_samples, burn_in = 1000, thinning = 1,
                                x_data, prior_hyperparams) {
  
  samples <- matrix(NA, nrow = n_samples, ncol = length(x0))
  x_current <- x0
  total_iterations <- burn_in + n_samples * thinning
  accepted_samples <- 0
  i <- 0
  
  message(paste0("Iniciando simulación MCMC con ", total_iterations, " iteraciones (", burn_in, " de burn-in, thinning de ", thinning, ")."))
  
  # Calcular la log-densidad inicial fuera del bucle para mayor eficiencia
  log_target_current <- target_log_pdf(x_data = x_data, current_params = x_current, prior_hyperparams = prior_hyperparams)
  
  while (i < total_iterations) {
    if (i > 0 && i %% (total_iterations / 10) == 0) {
      message(paste0("Progreso MCMC: ", round(i / total_iterations * 100), "% completado."))
    }
    
    x_proposed <- proposal_sampler(x_current)
    
    # Calcular las log-densidades
    log_target_proposed <- target_log_pdf(x_data = x_data, current_params = x_proposed, prior_hyperparams = prior_hyperparams)
    log_proposal_current_to_proposed <- proposal_log_pdf(y = x_proposed, z = x_current)
    log_proposal_proposed_to_current <- proposal_log_pdf(y = x_current, z = x_proposed)
    
    # Calcular la razón de aceptación en la escala logarítmica
    log_acceptance_ratio <- (log_target_proposed + log_proposal_proposed_to_current) -
      (log_target_current + log_proposal_current_to_proposed)
    
    # Aceptar o rechazar la propuesta
    if (is.nan(log_acceptance_ratio) || is.infinite(log_acceptance_ratio) || log(runif(1)) < log_acceptance_ratio) {
      x_current <- x_proposed
      log_target_current <- log_target_proposed
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
  message(paste0("Simulación MCMC finalizada. Tasa de aceptación: ", round(acceptance_rate, 4)))
  return(list(samples = samples, acceptance_rate = acceptance_rate))
}

# ==============================================================================
# --- 5. Funciones de Diagnóstico y Métricas ---
# ==============================================================================

# Medida y diagnóstico de las cadenas MCMC
Measure_Diagnostic <- function(data1, data2, burnin=0, thin=1, digits=5, cred_level=0.95, batch_size=100) {
  batch_means_stderr_var_R <- function(chain, batch_size) {
    n <- length(chain)
    num_batches <- floor(n / batch_size)
    if (num_batches < 2) {
      warning("Se requieren al menos 2 bloques para estimar el error estándar de la varianza. Retornando NA.")
      return(NA)
    }
    batch_vars <- unlist(lapply(0:(num_batches-1), function(i) {
      block_start <- i * batch_size + 1
      block_end <- (i + 1) * batch_size
      var(chain[block_start:block_end])
    }))
    return(sd(batch_vars) / sqrt(num_batches))
  }
  
  Cred_Interval_R <- function(x, level=0.95) {
    lower <- quantile(x, (1 - level)/2, na.rm = TRUE)
    upper <- quantile(x, 1 - (1 - level)/2, na.rm = TRUE)
    return(c(lower, upper))
  }
  
  Cred_Interval_V_R <- function(x, level=0.95) {
    x_clean <- na.omit(x)
    if (length(x_clean) < 2) {
      warning("Datos insuficientes para el bootstrap de la varianza después de NA. Retornando NA.")
      return(c(NA, NA))
    }
    n_bootstrap <- 1000
    variance_stat <- function(data, indices) {
      var(data[indices])
    }
    boot_obj <- tryCatch({
      boot(data = x_clean, statistic = variance_stat, R = n_bootstrap)
    }, error = function(e) {
      warning(paste("Error en bootstrap para varianza:", e$message))
      return(NULL)
    })
    if (is.null(boot_obj)) return(c(NA, NA))
    ci_results <- tryCatch({
      boot.ci(boot_obj, type = "perc", conf = level)$percent[4:5]
    }, error = function(e) {
      warning(paste("Error al calcular CI de bootstrap para varianza:", e$message))
      return(c(NA, NA))
    })
    return(round(ci_results, digits))
  }
  
  batch_means_stderr_cov_R <- function(chain, batch_size=100) {
    n <- nrow(chain)
    num_batches <- floor(n / batch_size)
    if (num_batches < 2) {
      warning("Se requieren al menos 2 bloques para estimar el error estándar de la covarianza. Retornando NA.")
      return(NA)
    }
    chain_clean <- na.omit(chain)
    if (nrow(chain_clean) < 2) {
      warning("Datos insuficientes para la covarianza de bloques después de NA. Retornando NA.")
      return(NA)
    }
    batch_covs <- unlist(lapply(0:(num_batches-1), function(i) {
      block_start <- i * batch_size + 1
      block_end <- (i + 1) * batch_size
      current_block <- na.omit(chain[block_start:block_end, ])
      if (nrow(current_block) < 2) return(NA)
      cov(current_block[, 1], current_block[, 2])
    }))
    batch_covs <- na.omit(batch_covs)
    if (length(batch_covs) < 2) {
      warning("Datos de covarianza de bloques insuficientes para sd. Retornando NA.")
      return(NA)
    }
    return(sd(batch_covs) / sqrt(length(batch_covs)))
  }
  
  Cred_Interval_Cov_R <- function(x, level=0.95) {
    x_clean <- na.omit(x)
    if (nrow(x_clean) < 2) {
      warning("Datos insuficientes para el bootstrap de la covarianza después de NA. Retornando NA.")
      return(c(NA, NA))
    }
    n_bootstrap <- 1000
    covariance_stat <- function(data, indices) {
      cov_matrix <- cov(data[indices, ])
      return(cov_matrix[1, 2])
    }
    boot_obj <- tryCatch({
      boot(data = x_clean, statistic = covariance_stat, R = n_bootstrap, sim = "ordinary")
    }, error = function(e) {
      warning(paste("Error en bootstrap para covarianza:", e$message))
      return(NULL)
    })
    if (is.null(boot_obj)) return(c(NA, NA))
    ci_results <- tryCatch({
      boot.ci(boot_obj, type = "perc", conf = level)$percent[4:5]
    }, error = function(e) {
      warning(paste("Error al calcular CI de bootstrap para covarianza:", e$message))
      return(c(NA, NA))
    })
    return(round(ci_results, digits))
  }
  
  if (!is.matrix(data1) && !is.data.frame(data1)) data1 <- as.matrix(data1)
  if (!is.matrix(data2) && !is.data.frame(data2)) data2 <- as.matrix(data2)
  
  combined_raw_data <- cbind(data1, data2)
  start_index <- burnin + 1
  if (start_index > nrow(combined_raw_data)) {
    warning("El burn-in es mayor que el número total de muestras. No hay datos para el análisis.")
    return(list("Numerical" = data.frame(Measure = character(), Value = numeric())))
  }
  seq_indices <- seq(from = start_index, to = nrow(combined_raw_data), by = thin)
  if (length(seq_indices) == 0) {
    warning("Después de burn-in y thinning, no quedan muestras. No hay datos para el análisis.")
    return(list("Numerical" = data.frame(Measure = character(), Value = numeric())))
  }
  
  new_data1 <- combined_raw_data[seq_indices, 1]
  new_data2 <- combined_raw_data[seq_indices, 2]
  new_data1_clean <- na.omit(new_data1)
  new_data2_clean <- na.omit(new_data2)
  combined_data_clean <- na.omit(cbind(new_data1, new_data2))
  
  if (length(new_data1_clean) < 2 || length(new_data2_clean) < 2) {
    warning("Datos insuficientes para calcular estadísticas después de eliminar NAs. Retornando NA para algunas métricas.")
    mean_x1 <- var_x1 <- ess_x1 <- stderr_mean_x1 <- stderr_var_x1 <- NA
    CIM1 <- c(NA, NA)
    CIV1 <- c(NA, NA)
    mean_x2 <- var_x2 <- ess_x2 <- stderr_mean_x2 <- stderr_var_x2 <- NA
    CIM2 <- c(NA, NA)
    CIV2 <- c(NA, NA)
    cov_val <- stderr_cov <- CIcov <- c(NA, NA)
    len_data <- length(new_data1)
  } else {
    mcmc_data1 <- as.mcmc(new_data1_clean)
    mcmc_data2 <- as.mcmc(new_data2_clean)
    
    mean_x1 <- mean(new_data1_clean)
    var_x1 <- var(new_data1_clean)
    ess_x1 <- effectiveSize(mcmc_data1)
    stderr_mean_x1 <- sd(new_data1_clean) / sqrt(ess_x1)
    stderr_var_x1 <- batch_means_stderr_var_R(new_data1_clean, batch_size)
    CIM1 <- Cred_Interval_R(new_data1_clean, cred_level)
    CIV1 <- Cred_Interval_V_R(new_data1_clean, cred_level)
    
    mean_x2 <- mean(new_data2_clean)
    var_x2 <- var(new_data2_clean)
    ess_x2 <- effectiveSize(mcmc_data2)
    stderr_mean_x2 <- sd(new_data2_clean) / sqrt(ess_x2)
    stderr_var_x2 <- batch_means_stderr_var_R(new_data2_clean, batch_size)
    CIM2 <- Cred_Interval_R(new_data2_clean, cred_level)
    CIV2 <- Cred_Interval_V_R(new_data2_clean, cred_level)
    
    cov_val <- cov(combined_data_clean[, 1], combined_data_clean[, 2])
    stderr_cov <- batch_means_stderr_cov_R(combined_data_clean, batch_size)
    CIcov <- Cred_Interval_Cov_R(combined_data_clean, cred_level)
    len_data <- length(new_data1)
  }
  
  new_names <- c("Mean_X1", "Var_X1", "ESS_X1", "STDERR_Mean_X1", "STDERR_Var_X1",
                 "CI_Mean1_Lower", "CI_Mean1_Upper", "CI_Var1_Lower", "CI_Var1_Upper",
                 "Mean_X2", "Var_X2", "ESS_X2", "STDERR_Mean_X2", "STDERR_Var_X2",
                 "CI_Mean2_Lower", "CI_Mean2_Upper", "CI_Var2_Lower", "CI_Var2_Upper",
                 "Cov", "STDERR_Cov", "CI_Cov_Lower", "CI_Cov_Upper", "Length")
  
  results_vector <- c(
    mean_x1, var_x1, ess_x1, stderr_mean_x1, stderr_var_x1, CIM1[1], CIM1[2], CIV1[1], CIV1[2],
    mean_x2, var_x2, ess_x2, stderr_mean_x2, stderr_var_x2, CIM2[1], CIM2[2], CIV2[1], CIV2[2],
    cov_val, stderr_cov, CIcov[1], CIcov[2], len_data
  )
  
  Numerical_results_vertical <- data.frame(
    Measure = new_names,
    Value = round(results_vector, digits)
  )
  
  return(list("Numerical" = Numerical_results_vertical))
}

# ==============================================================================
# --- 6. Monitoreo de Convergencia (Visualizaciones) ---
# ==============================================================================

# Genera gráficos para diagnosticar la convergencia de la cadena.
monitor_convergence <- function(x, y) {
  df_plot <- data.frame(Iteration = seq_along(x), Alpha = x, Beta = y) %>%
    na.omit()
  if (nrow(df_plot) == 0) {
    message("ADVERTENCIA: No hay datos válidos para generar gráficos de convergencia después de eliminar NAs.")
    return(invisible(NULL))
  }
  
  plot_trace_alpha <- ggplot(df_plot, aes(x = Iteration, y = Alpha)) +
    geom_line(color = "blue", alpha = 0.7) +
    labs(title = "Gráfico de Traza de la Cadena para Alpha", x = "Iteraciones", y = "Valor de Alpha") +
    theme_minimal()
  plot_trace_beta <- ggplot(df_plot, aes(x = Iteration, y = Beta)) +
    geom_line(color = "red", alpha = 0.7) +
    labs(title = "Gráfico de Traza de la Cadena para Beta", x = "Iteraciones", y = "Valor de Beta") +
    theme_minimal()
  gridExtra::grid.arrange(plot_trace_alpha, plot_trace_beta, ncol = 1)
  
  plot_hist_alpha <- ggplot(df_plot, aes(x = Alpha)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "blue", alpha = 0.7) +
    geom_density(color = "darkblue", linewidth = 1) +
    labs(title = "Histograma de Alpha", x = "Valor", y = "Densidad") +
    theme_minimal()
  plot_hist_beta <- ggplot(df_plot, aes(x = Beta)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "red", alpha = 0.7) +
    geom_density(color = "darkred", linewidth = 1) +
    labs(title = "Histograma de Beta", x = "Valor", y = "Densidad") +
    theme_minimal()
  gridExtra::grid.arrange(plot_hist_alpha, plot_hist_beta, ncol = 2)
  
  plot_scatter <- ggplot(df_plot, aes(x = Alpha, y = Beta)) +
    geom_point(alpha = 0.3, color = "darkgreen") +
    labs(title = "Diagrama de Dispersión de las Muestras (Alpha, Beta)", x = "Alpha", y = "Beta") +
    theme_minimal() + theme(aspect.ratio = 1)
  print(plot_scatter)
  
  old_par <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1) + 0.1)
  acf(df_plot$Alpha, lag.max = 50, main = "Autocorrelación de Alpha", xlab = "Retardo (Lag)", ylab = "ACF", na.action = na.pass)
  acf(df_plot$Beta, lag.max = 50, main = "Autocorrelación de Beta", xlab = "Retardo (Lag)", ylab = "ACF", na.action = na.pass)
  par(old_par)
}

# ==============================================================================
# --- 7. Ejecución de la Simulación MCMC y Análisis con Prior Gamma ---
# ==============================================================================

message("\n--- Ejecutando MCMC con Priors Gamma ---")
prior_hyperparams_gammas <- c(a_shape = 0.1, b_rate = 0.1, c_shape = 0.1, d_rate = 0.1)
n_samples_mcmc <- 20000
burn_in_mcmc <- 5000
thinning_mcmc <- 5

results_mcmc_gammas <- metropolis_hastings(
  target_log_pdf = post_prior_gammas,
  proposal_sampler = proposal_sampler,
  proposal_log_pdf = proposal_log_pdf,
  x0 = x0,
  n_samples = n_samples_mcmc,
  burn_in = burn_in_mcmc,
  thinning = thinning_mcmc,
  x_data = xn,
  prior_hyperparams = prior_hyperparams_gammas
)

samples_gammas <- results_mcmc_gammas$samples
acceptance_rate_gammas <- results_mcmc_gammas$acceptance_rate

if (!is.null(samples_gammas) && nrow(samples_gammas) > 0) {
  alpha_samples_gammas <- samples_gammas[, 1]
  beta_samples_gammas <- samples_gammas[, 2]
  message(paste("Tasa de Aceptación (Gammas):", round(acceptance_rate_gammas, 6)))
  message(paste0("Media (alpha, beta) (Gammas): (", round(colMeans(samples_gammas, na.rm = TRUE)[1], 6), ", ", round(colMeans(samples_gammas, na.rm = TRUE)[2], 6), ")"))
  message(paste0("Varianza (alpha, beta) (Gammas): (", round(apply(samples_gammas, 2, var, na.rm = TRUE)[1], 6), ", ", round(apply(samples_gammas, 2, var, na.rm = TRUE)[2], 6), ")"))
  message("\nMostrando gráficos de convergencia para Priors Gamma...")
  monitor_convergence(alpha_samples_gammas, beta_samples_gammas)
} else {
  message("ADVERTENCIA: No hay muestras válidas para el análisis (Gammas).")
}

# ==============================================================================
# --- 8. Ejecución de la Simulación MCMC y Análisis con Prior de Jeffreys ---
# ==============================================================================

message("\n--- Ejecutando MCMC con Prior de Jeffreys ---")

results_mcmc_jeffreys <- metropolis_hastings(
  target_log_pdf = post_prior_Jeffrey,
  proposal_sampler = proposal_sampler,
  proposal_log_pdf = proposal_log_pdf,
  x0 = x0,
  n_samples = n_samples_mcmc,
  burn_in = burn_in_mcmc,
  thinning = thinning_mcmc,
  x_data = xn,
  prior_hyperparams = NULL
)

samples_jeffreys <- results_mcmc_jeffreys$samples
acceptance_rate_jeffreys <- results_mcmc_jeffreys$acceptance_rate

if (!is.null(samples_jeffreys) && nrow(samples_jeffreys) > 0) {
  alpha_samples_jeffreys <- samples_jeffreys[, 1]
  beta_samples_jeffreys <- samples_jeffreys[, 2]
  message(paste("Tasa de Aceptación (Jeffreys):", round(acceptance_rate_jeffreys, 6)))
  message(paste0("Media (alpha, beta) (Jeffreys): (", round(colMeans(samples_jeffreys, na.rm = TRUE)[1], 6), ", ", round(colMeans(samples_jeffreys, na.rm = TRUE)[2], 6), ")"))
  message(paste0("Varianza (alpha, beta) (Jeffreys): (", round(apply(samples_jeffreys, 2, var, na.rm = TRUE)[1], 6), ", ", round(apply(samples_jeffreys, 2, var, na.rm = TRUE)[2], 6), ")"))
  message("\nMostrando gráficos de convergencia para Prior de Jeffreys...")
  monitor_convergence(alpha_samples_jeffreys, beta_samples_jeffreys)
} else {
  message("ADVERTENCIA: No hay muestras válidas para el análisis (Jeffreys).")
}

# ==============================================================================
# --- 8b. Ejecución de la Simulación MCMC y Análisis con New Prior (Llerzy & Tovar) ---
# ==============================================================================

message("\n--- Ejecutando MCMC con New Prior (Llerzy & Tovar) ---")

prior_hyperparams_new <- c(a = 1, b = 1, c = 1, d = 1)

results_mcmc_new <- metropolis_hastings(
  target_log_pdf = post_new_prior,
  proposal_sampler = proposal_sampler,
  proposal_log_pdf = proposal_log_pdf,
  x0 = x0,
  n_samples = n_samples_mcmc,
  burn_in = burn_in_mcmc,
  thinning = thinning_mcmc,
  x_data = xn,
  prior_hyperparams = prior_hyperparams_new
)

samples_new <- results_mcmc_new$samples
acceptance_rate_new <- results_mcmc_new$acceptance_rate

if (!is.null(samples_new) && nrow(samples_new) > 0) {
  alpha_samples_new <- samples_new[, 1]
  beta_samples_new <- samples_new[, 2]
  
  message(paste("Tasa de Aceptación (New Prior):", round(acceptance_rate_new, 6)))
  message(paste0("Media (alpha, beta) (New Prior): (", 
                 round(colMeans(samples_new, na.rm = TRUE)[1], 6), ", ", 
                 round(colMeans(samples_new, na.rm = TRUE)[2], 6), ")"))
  message(paste0("Varianza (alpha, beta) (New Prior): (", 
                 round(apply(samples_new, 2, var, na.rm = TRUE)[1], 6), ", ", 
                 round(apply(samples_new, 2, var, na.rm = TRUE)[2], 6), ")"))
  
  message("\nMostrando gráficos de convergencia para New Prior...")
  monitor_convergence(alpha_samples_new, beta_samples_new)
  
  results_new <- Measure_Diagnostic(
    data1 = alpha_samples_new,
    data2 = beta_samples_new,
    burnin = 0, thin = 1, digits = 5, cred_level = 0.95, batch_size = 100
  )
  message("\nResultados de Diagnóstico (New Prior):")
  print(results_new$Numerical)
} else {
  message("ADVERTENCIA: No hay muestras válidas para el análisis (New Prior).")
}

# ==============================================================================
# --- 9. Consolidación y Comparación de Resultados ---
# ==============================================================================

message("\n--- Consolidando Resultados para Comparación ---")

# Define las variables de resultados como NULL por defecto
results_gammas_df <- NULL
results_jeffreys_df <- NULL
results_new_df <- NULL

# Asigna los resultados solo si la variable existe y no es NULL
if (exists("results_gammas") && !is.null(results_gammas)) {
  results_gammas_df <- Measure_Diagnostic(
    data1 = samples_gammas[, 1],
    data2 = samples_gammas[, 2],
    burnin = 0, thin = 1, digits = 5, cred_level = 0.95, batch_size = 100
  )$Numerical
}

if (exists("results_mcmc_jeffreys") && !is.null(results_mcmc_jeffreys)) {
  results_jeffreys_df <- Measure_Diagnostic(
    data1 = samples_jeffreys[, 1],
    data2 = samples_jeffreys[, 2],
    burnin = 0, thin = 1, digits = 5, cred_level = 0.95, batch_size = 100
  )$Numerical
}

if (exists("results_mcmc_new") && !is.null(results_mcmc_new)) {
  results_new_df <- Measure_Diagnostic(
    data1 = samples_new[, 1],
    data2 = samples_new[, 2],
    burnin = 0, thin = 1, digits = 5, cred_level = 0.95, batch_size = 100
  )$Numerical
}

# Realiza la unión de los dataframes solo si todos existen
if (!is.null(results_gammas_df) && !is.null(results_jeffreys_df) && !is.null(results_new_df)) {
  combined_results <- results_gammas_df %>%
    dplyr::rename(Gammas = Value) %>%
    dplyr::mutate(Jeffreys = results_jeffreys_df$Value,
                  New = results_new_df$Value)
  
  message("\nTabla Comparativa de Resultados Diagnósticos:")
  print(combined_results)
} else {
  message("ADVERTENCIA: No se pueden comparar los resultados porque una o más simulaciones no produjeron muestras válidas.")
}

# ==============================================================================
# --- 10. Simulación de la Distribución Predictiva Posterior ---
# ==============================================================================
message("\n--- Realizando simulación de la Distribución Predictiva Posterior ---")

posterior_alpha_samples <- samples_jeffreys[, 1]
posterior_beta_samples <- samples_jeffreys[, 2]

posterior_alpha_samples <- na.omit(posterior_alpha_samples)
posterior_beta_samples <- na.omit(posterior_beta_samples)

if (length(posterior_alpha_samples) == 0 || length(posterior_beta_samples) == 0) {
  message("ADVERTENCIA: No hay muestras válidas de Alpha o Beta para la simulación predictiva posterior. Saltando la simulación.")
} else {
  n_simulations <- min(length(posterior_alpha_samples), length(posterior_beta_samples))
  
  if (length(posterior_alpha_samples) != length(posterior_beta_samples)) {
    message(paste0("ADVERTENCIA: Las longitudes de las muestras de Alpha y Beta difieren. Usando el mínimo (", n_simulations, ") para la simulación."))
    posterior_alpha_samples <- posterior_alpha_samples[1:n_simulations]
    posterior_beta_samples <- posterior_beta_samples[1:n_simulations]
  }
  
  simulated_proportions <- rbeta(n_simulations,
                                 shape1 = posterior_alpha_samples,
                                 shape2 = posterior_beta_samples)
  
  plot_simulated_proportions <- ggplot(data.frame(Proportion = simulated_proportions), aes(x = Proportion)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "purple", alpha = 0.7) +
    geom_density(color = "darkmagenta", linewidth = 1) +
    labs(title = "Distribución Predictiva Posterior de Proporciones Defectuosas",
         x = "Proporción Defectuosa Simulada", y = "Densidad") +
    theme_minimal() +
    xlim(0, 1)
  print(plot_simulated_proportions)
  message("Gráfico de la distribución predictiva posterior generado.")
  
  sample_size_new <- 100
  simulated_defects_in_sample <- rbinom(n_simulations, size = sample_size_new, prob = simulated_proportions)
}
