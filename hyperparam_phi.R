# ==============================================================================
# PROCEDIMIENTO COMPLETO: BOOTSTRAP + TOVAR CON INCERTIDUMBRE ARTIFICIAL
# ==============================================================================

vect_phi <- function(data, B = 10000 ){
  
  data_clean <- data[!is.na(data)]
  
  n <- length(data_clean)
  
  if (n < 2) return(c(NA, NA, NA, NA))
  
  boot_means <- numeric(B)
  boot_vars <- numeric(B)
  boot_W <- numeric(B) # Nuevo vector para la varianza escalada
  
  for (i in 1:B) {
    muestra_virtual <- sample(data_clean, size = n, replace = TRUE)
    
    m <- mean(muestra_virtual)
    v <- var(muestra_virtual)
    
    boot_means[i] <- m
    boot_vars[i] <- v
    
    # Calculamos W (Varianza Escalada) para cada muestra virtual
    boot_W[i] <- v / (m * (1 - m))
  }
  
  # ==============================================================================
  # PASO 1: HIPERPARÁMETROS a Y b (La Media) - Método de Momentos
  # ==============================================================================
  mu_media  <- mean(boot_means, na.rm = TRUE)
  var_media <- var(boot_means, na.rm = TRUE)
  
  b_est <- (1 - mu_media) * ((mu_media * (1 - mu_media) / var_media) - 1)
  a_est <- mu_media * b_est / (1 - mu_media)
  
  # ==============================================================================
  # PASO 2 DEFINITIVO: CUANTIL FIJO Y PENALIZACIÓN DE CHEBYSHEV
  # ==============================================================================
  # Usamos el vector boot_W (Varianza Escalada) que calculamos dentro del Bootstrap
  theta_0 <- mean(boot_W, na.rm = TRUE)
  
  # 1. LA CLAVE: Sacamos UN SOLO límite superior estándar (IC del 95%)
  # El percentil 97.5% deja 2.5% de probabilidad en la cola derecha
  theta_1_fijo <- as.numeric(quantile(boot_W, probs = 0.975, na.rm = TRUE))
  
  # 2. La distancia empírica al cuadrado se queda fija para todos los escenarios
  distancia_sq <- (theta_1_fijo - theta_0)^2
  
  # 3. Los escenarios actúan como el factor de penalización (alpha) de Tovar
  niveles_alpha <- 0.10 #c(0.05, 0.10, 0.20)
  
  for (alpha in niveles_alpha) {
    
    # Aplicamos la fórmula de Tovar: Varianza = alpha * distancia^2
    sigma_2_tovar <- alpha * distancia_sq
    
    # Calculamos c y d
    d_est <- (1 - theta_0) * ((theta_0 * (1 - theta_0) / sigma_2_tovar) - 1)
    c_est <- theta_0 * d_est / (1 - theta_0)
    
  }
  return(c(a_est, b_est, c_est, d_est))
}

phi_prov1 <- vect_phi(vec_data1)
phi_prov2 <- vect_phi(vec_data2)

phi_prov3 <- matrix(0, nrow = nrow(Matriz_Data3), ncol = 4)

for (phis in 1:nrow(Matriz_Data3)) {
  phi_prov3[phis, ] <- vect_phi(Matriz_Data3[phis, ])
}
