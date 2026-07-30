# ============================================================================== #
# SCRIPT 1: ESTIMACIÓN BAYESIANA DE LA CALIDAD (OBJETIVO 1B)
# Pura extracción matemática. Las gráficas se delegan a otro script.
# ============================================================================== #

get_initial_values <- function(data) {
  mean_xn <- mean(data); var_xn <- var(data)
  if (is.na(var_xn) || var_xn == 0) { return(c(0.1, 0.1)) }
  
  a_0 <- mean_xn * (mean_xn * (1 - mean_xn) / var_xn - 1)
  b_0 <- (1 - mean_xn) * (mean_xn * (1 - mean_xn) / var_xn - 1)
  if (is.na(a_0) || is.na(b_0) || a_0 <= 0 || b_0 <= 0) { a_0 <- 0.001; b_0 <- 0.001 }
  return(c(a_0, b_0))
}

proposal_sampler <- function(y, sigma) {
  mu1 <- if (y[1] > 0) log(y[1]) else log(0.0001)
  mu2 <- if (y[2] > 0) log(y[2]) else log(0.0001)
  return(c(rlnorm(1, mu1, sigma), rlnorm(1, mu2, sigma)))
}

proposal_log_pdf <- function(y, z, sigma) {
  mu1 <- if (z[1] > 0) log(z[1]) else log(0.0001)
  mu2 <- if (z[2] > 0) log(z[2]) else log(0.0001)
  return(dlnorm(y[1], mu1, sigma, log = TRUE) + dlnorm(y[2], mu2, sigma, log = TRUE))
}

target_log_pdf_gamma <- function(current_params, x_data, prior_hyperparams) {
  a <- current_params[1]; b <- current_params[2]
  if (a <= 0 || b <= 0 || a > 1e4 || b > 1e4) return(-Inf)
  log_likelihood <- length(x_data) * (lgamma(a + b) - lgamma(a) - lgamma(b)) + 
    (a - 1) * sum(log(x_data)) + (b - 1) * sum(log(1 - x_data))
  log_prior <- dgamma(a, prior_hyperparams[1], prior_hyperparams[2], log=TRUE) + 
    dgamma(b, prior_hyperparams[3], prior_hyperparams[4], log=TRUE)
  return(log_likelihood + log_prior)
}

target_log_pdf_jeffrey <- function(current_params, x_data, prior_hyperparams = NULL) {
  a <- current_params[1]; b <- current_params[2]
  if (a <= 0 || b <= 0 || a > 1e4 || b > 1e4) return(-Inf)
  term <- psigamma(a, 1) * psigamma(b, 1) - (psigamma(a, 1) + psigamma(b, 1)) * psigamma(a + b, 1)
  if (is.na(term) || term <= 0) return(-Inf) 
  log_likelihood <- length(x_data) * (lgamma(a + b) - lgamma(a) - lgamma(b)) + 
    (a - 1) * sum(log(x_data)) + (b - 1) * sum(log(1 - x_data))
  return(log_likelihood + 0.5 * log(term))
}

target_log_pdf_bivariada <- function(current_params, x_data, prior_hyperparams) {
  a <- current_params[1]; b <- current_params[2]
  if (a <= 0 || b <= 0 || a > 1e4 || b > 1e4) return(-Inf)
  log_likelihood <- length(x_data) * (lgamma(a + b) - lgamma(a) - lgamma(b)) + 
    (a - 1) * sum(log(x_data)) + (b - 1) * sum(log(1 - x_data))
  log_prior <- (prior_hyperparams[1] - 1) * log(a) + (prior_hyperparams[2] - 1) * log(b) - 
    (prior_hyperparams[1] + prior_hyperparams[2]) * log(a + b) - 
    (prior_hyperparams[3] + prior_hyperparams[4]) * log(a + b + 1)
  return(log_likelihood + log_prior)
}

metropolis_hastings <- function(target_log_pdf, x0, n_samples, burn_in,
                                thinning, x_data, prior_hyperparams=NULL, 
                                usar_adaptativo = FALSE) {
  
  total_iters <- burn_in + n_samples * thinning
  samples <- matrix(NA, nrow = n_samples, ncol = 2)
  
  # Trabajaremos en escala logarítmica para evitar restricciones > 0
  log_x_current <- log(x0)
  x_current <- x0
  
  # --- INICIO LÓGICA: ADAPTIVE METROPOLIS (SALTO CORRELACIONADO) ---
  cov_mat <- diag(c(0.1, 0.1))   # Matriz de covarianza inicial (independiente)
  eps_mat <- diag(c(1e-6, 1e-6)) # Para evitar que la matriz colapse a 0
  escala <- 1.0                  # Multiplicador del tamaño de salto global
  
  lote_adap <- 50
  aceptados_lote <- 0
  
  # Memoria para aprender la "forma" del pasillo durante el burn-in
  historia_log <- matrix(NA, nrow = burn_in, ncol = 2)
  # -----------------------------------------------------------------
  
  for (i in 1:total_iters) {
    
    # 1. Proponer nuevo valor (Caminata Aleatoria Correlacionada)
    # chol() extrae la raíz cuadrada de la matriz de covarianza
    L <- chol(escala * cov_mat + eps_mat)
    salto <- as.numeric(t(L) %*% rnorm(2)) # Salto multivariado
    
    log_x_prop <- log_x_current + salto
    x_prop <- exp(log_x_prop)
    
    # 2. Ratio de Aceptación (Incluye ajuste del Jacobiano logarítmico)
    # Al hacer sum(log_x) nos ahorramos calcular dlnorm y proposal_log_pdf. Es mucho más rápido.
    log_r <- target_log_pdf(x_prop, x_data, prior_hyperparams) + sum(log_x_prop) -
      (target_log_pdf(x_current, x_data, prior_hyperparams) + sum(log_x_current))
    
    if (is.nan(log_r) || is.na(log_r)) log_r <- -Inf
    
    # 3. Decisión
    if (log(runif(1)) < log_r) {
      x_current <- x_prop
      log_x_current <- log_x_prop
      if (i <= burn_in) aceptados_lote <- aceptados_lote + 1
    }
    
    # 4. Aprendizaje / Adaptación (Solo durante burn-in)
    if (i <= burn_in) {
      historia_log[i, ] <- log_x_current # Guardamos el paso para la historia
      
      if (usar_adaptativo && i %% lote_adap == 0) {
        tasa <- aceptados_lote / lote_adap
        
        # A) Calibramos la distancia del salto
        if (tasa < 0.23) {
          escala <- escala * 0.8
        } else if (tasa > 0.35) {
          escala <- escala * 1.2
        }
        
        # B) Aprendemos la dirección: actualizamos la matriz de covarianza empírica
        # Esperamos iteración 100 para tener suficientes puntos y evitar errores matemáticos
        if (i >= 100) {
          cov_mat <- cov(historia_log[1:i, ])
        }
        
        aceptados_lote <- 0
      } else if (!usar_adaptativo && i %% lote_adap == 0) {
        aceptados_lote <- 0
      }
    }
    
    # 5. Guardar la muestra oficial
    if (i > burn_in && (i - burn_in) %% thinning == 0) {
      samples[(i - burn_in) / thinning, ] <- x_current
    }
  }
  
  return(samples)
}
# ============================================================================== #
# --- 2. EJECUCIÓN PRINCIPAL: MOTOR MCMC ---
# ============================================================================== #

diccionario_emp_temp <- data.frame()
Cadenas_MCMC_List <- list() # Aquí guardaremos las cadenas para el script gráfico
#n_s <- 7500; b_in <- 1500; thin <- 3

correr_3_priors <- function(nombre_fuente, datos_crudos, usar_adaptativo = FALSE) {
  datos_limpios <- pmax(pmin(as.numeric(na.omit(datos_crudos)), 1 - 1e-5), 1e-5)
  theta_0 <- get_initial_values(datos_limpios)
  
  res_gamma <- metropolis_hastings(target_log_pdf_gamma, theta_0, n_s, b_in, thin, datos_limpios, c(0.1, 0.1, 0.1, 0.1), usar_adaptativo)
  res_jeffrey <- metropolis_hastings(target_log_pdf_jeffrey, theta_0, n_s, b_in, thin, datos_limpios, NULL, usar_adaptativo)
  res_biv <- metropolis_hastings(target_log_pdf_bivariada, theta_0, n_s, b_in, thin, datos_limpios, c(2, 2, 2, 2), usar_adaptativo)
  
  # Preparar la fila para el diccionario
  fila <- data.frame(
    Fuente = nombre_fuente,
    Prov_Gamma_a = mean(res_gamma[, 1]), Prov_Gamma_b = mean(res_gamma[, 2]),
    Prov_Jeffreys_a = mean(res_jeffrey[, 1]), Prov_Jeffreys_b = mean(res_jeffrey[, 2]),
    Prov_TTB_a = mean(res_biv[, 1]), Prov_TTB_b = mean(res_biv[, 2])
  )
  
  # Preparar las cadenas para el script gráfico
  cadenas <- list(Gamma = res_gamma, Jeffreys = res_jeffrey, TTB = res_biv)
  
  return(list(fila = fila, cadenas = cadenas))
}

# --- 2.1 PROVEEDOR 1 (Clásico / Empírico) ---
cat("\n>> Ejecutando MCMC: Proveedor 1...\n")
n_s <- 5000; b_in <- 2000; thin <- 2
res_p1 <- correr_3_priors("Prov1", vec_data1, usar_adaptativo = FALSE)
diccionario_emp_temp <- rbind(diccionario_emp_temp, res_p1$fila)
Cadenas_MCMC_List[["Prov1"]] <- res_p1$cadenas

# --- 2.2 PROVEEDOR 2 (Clásico / Empírico) ---
cat(">> Ejecutando MCMC: Proveedor 2...\n")
n_s <- 5000; b_in <- 2000; thin <- 2
res_p2 <- correr_3_priors("Prov2", vec_data2, usar_adaptativo = FALSE)
diccionario_emp_temp <- rbind(diccionario_emp_temp, res_p2$fila)
Cadenas_MCMC_List[["Prov2"]] <- res_p2$cadenas

# --- 2.3 PROVEEDOR 3 (ARWMH + Fuerza Bruta / Malla Extrema) ---
cat("\n>> Ejecutando MCMC: Proveedor 3 (Malla de Escenarios)...\n")
n_s <- 6000; b_in <- 3000; thin <- 10; USAR_ADAPTATIVO <- TRUE

# Asegúrese de que este ciclo FOR esté presente y activo:
for (k in 1:nrow(Esce_Global)) {
  nombre_esce <- paste0("Prov3_Esce_", k)
  
  # Imprimimos en consola para saber en qué escenario va (muy útil para monitorear)
  cat("   -> Procesando:", nombre_esce, "\n") 
  
  # Extraemos los 35 lotes simulados directamente desde la matriz maestra
  vec_data3 <- as.numeric(Matriz_Data3[k, ])
  
  res_p3 <- correr_3_priors(nombre_esce, vec_data3, usar_adaptativo = TRUE)
  diccionario_emp_temp <- rbind(diccionario_emp_temp, res_p3$fila)
  Cadenas_MCMC_List[[nombre_esce]] <- res_p3$cadenas
}

# Guardar en memoria global
diccionario_empirico <<- diccionario_emp_temp
cat("=== RESULTADO FASE I: 'diccionario_empirico' generado ===\n")