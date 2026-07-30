
# Probabilidad de Aceptación (PA) - Curva CO tipo A (Hypergeométrica)
Pa <- function(n, c, p, N){
  # Se utiliza phyper para la distribución hipergeométrica (muestreo sin reemplazo)
  Pa_val <- phyper(c, N * p, N * (1 - p), n)
  return(Pa_val)
}

# Integrando el Riesgo del Productor (RP): (1 - PA) * f(p)
# Error tipo I (Rechazar lote bueno): p en [0, AQL]
rap_p <- function(p, n, c, a_b, b_b, N_){
  f.p <- dbeta(p, a_b, b_b) # Densidad del historico de calidad
  rap_p_ <- (1 - Pa(n, c, p, N_)) * f.p
  return(rap_p_)
}

# Integrando el Riesgo del Consumidor (RC): PA * f(p)
# Error tipo II (Aceptar lote malo): p en [LTPD, 1]
rap_c <- function(p, n, c, a_b, b_b, N_){
  f.p <- dbeta(p, a_b, b_b) # Densidad a priori Beta
  rap_c_r <- Pa(n, c, p, N_) * f.p
  return(rap_c_r)
}

# calc_rap(N_, n, c, a_b=diccionario_base %>% select(a) %>% )
# Función principal para calcular el Riesgo Ponderado Integrado (RP y RC)
calc_rap <- function(N_, n, c, a_b, b_b, AQL, LTPD) {
  
  rap_p_val <- integral(f = function(p) rap_p(p, n, c, a_b, b_b, N_), xmin = 0, xmax = AQL, method = "Kron")
  rap_c_val <- integral(f = function(p) rap_c(p, n, c, a_b, b_b, N_), xmin = LTPD, xmax = 1, method = "Kron")
  rap_t_val <- rap_p_val + rap_c_val
  
  # CORRECCIÓN AQUÍ: Se retornan las variables correctas
  return(c(RAP_p_val = rap_p_val, RAP_c_val = rap_c_val, RAP_t_val = rap_t_val)) 
}