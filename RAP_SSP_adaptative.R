# ============================================================================== #
# SCRIPT 3: OPTIMIZACIÓN DEL PLAN DE MUESTREO (FASE II)
# Metodología: Búsqueda Secuencial Determinística (Riesgo Acumulado Ponderado)
# Autores Base: Fernández et al. (2020) | Adaptación: Juan Sebastián Haeusler
# ============================================================================== #

# --- 1. Carga de Librerías ---
paquetes <- c("pracma", "AcceptanceSampling", "sjstats")
for (p in paquetes) {
  if(!require(p, character.only = TRUE)) {
    install.packages(p)
    require(p, character.only = TRUE)
  }
}

# ============================================================================== #
# --- 2. Funciones de Riesgo Ponderado y Masa de Probabilidad ---
# ============================================================================== #

# Función de Método de Cook (¡DESCOMENTADA! Vital para los escenarios de expertos)
find_beta <- function(x1, p1, x2, p2) {
  loss_function <- function(params) {
    a <- params[1]; b <- params[2]
    if(a <= 0 || b <= 0) return(Inf)
    error <- (pbeta(x1, a, b) - p1)^2 + (pbeta(x2, a, b) - p2)^2
    return(error)
  }
  opt <- optim(par = c(1, 1), fn = loss_function, method = "Nelder-Mead")
  if(opt$convergence == 0) return(list(shape1 = opt$par[1], shape2 = opt$par[2]))
  else return(list(shape1 = NA, shape2 = NA))
}

# Probabilidad de Aceptación (PA) - Hipergeométrica
Pa <- function(n, c, p, N){
  Pa_val <- phyper(c, round(N * p), N - round(N * p), n)
  return(Pa_val)
}

# Riesgo del Productor Ponderado (RPP)
wr_p <- function(p, n, c, alpha_b, beta_b, N_){
  f.p <- dbeta(p, alpha_b, beta_b)
  prod_wr <- (1 - Pa(n, c, p, N_)) * f.p
  return(prod_wr)
}

# Riesgo del Consumidor Ponderado (RCP)
wr_c <- function(p, n, c, alpha_b, beta_b, N_){
  f.p <- dbeta(p, alpha_b, beta_b)
  cons_wr <- Pa(n, c, p, N_) * f.p
  return(cons_wr)
}

# Masa de Probabilidad (Zonas de Calidad)
calc_prob_mass <- function(alpha_b, beta_b, AQL, LTPD) {
  P_Good <- pbeta(AQL, shape1 = alpha_b, shape2 = beta_b)
  P_Bad <- 1 - pbeta(LTPD, shape1 = alpha_b, shape2 = beta_b)
  return(c(P_Good = P_Good, P_Bad = P_Bad))
}

# Integración Numérica con Pracma (Método Kronrod)
calc_wr <- function(N_, n, c, alpha_b, beta_b, AQL, LTPD) {
  wrp_val <- integral(f = function(p) wr_p(p, n, c, alpha_b, beta_b, N_),
                      xmin = 0, xmax = AQL, method = "Kron")
  wrc_val <- integral(f = function(p) wr_c(p, n, c, alpha_b, beta_b, N_),
                      xmin = LTPD, xmax = 1, method = "Kron")
  wrt_val <- wrp_val + wrc_val
  return(c(WRP_val = wrp_val, WRC_val = wrc_val, WRT_val = wrt_val))
}

# ============================================================================== #
# --- 3. Variables Globales, Escenarios y Perfiles ---
# ============================================================================== #

# N <- 1000                      
# alpha <- c(0.01, 0.05)         
# beta <- c(0.05, 0.10, 0.20)    
# AQL <- c(0.01, 0.02, 0.05)     
# LTPD <- c(0.08, 0.10, 0.15, 0.20)    
# 
# Esce <- expand.grid(N = N, alpha = alpha, beta = beta, AQL = AQL, LTPD = LTPD)
# Esce <- Esce[Esce$AQL < Esce$LTPD, ]
# rownames(Esce) <- NULL
# 
# # Perfiles Históricos (Cook)
p1 <- c(0.95, 0.75, 0.50, 0.25, 0.10, 0.05, 0.01)
p2 <- c(1e-4, 0.05, 0.10, 0.25, 0.40, 0.05, 0.01)
# 
# # Parámetros MCMC (Traídos del Objetivo 1 - Valores exactos del Script 1)
# # CORRECCIÓN: Se definen los 3 MCMC para que coincida con las 11 etiquetas
# a_mcmc_1 <- 5.6976;   b_mcmc_1 <- 102.6596  # MCMC 1: Gamma
# a_mcmc_2 <- 14.8335;  b_mcmc_2 <- 281.6133  # MCMC 2: Jeffreys
# a_mcmc_3 <- 10.1686;  b_mcmc_3 <- 190.0793  # MCMC 3: Bivariada (OTB)

lista_tablas_resultados <- list()

# ============================================================================== #
# --- 4. Motor de Búsqueda Determinística ---
# ============================================================================== #

cases <- 1:3 # Evaluando los primeros 3 escenarios para prueba

for (esce in cases) {
  
  cat(sprintf("\n--- Procesando Escenario %d de %d ---\n", esce, max(cases)))
  cat(sprintf("AQL: %.2f | LTPD: %.2f | Alpha: %.2f | Beta: %.2f\n", 
              Esce[esce, 4], Esce[esce, 5], Esce[esce, 2], Esce[esce, 3]))
  
  # 1. Plan Clásico (Naive)
  plan_clasic <- find.plan(PRP = c(Esce[esce, 4], 1 - Esce[esce, 2]),
                           CRP = c(Esce[esce, 5], Esce[esce, 3]),
                           N = Esce[esce, 1], type = "hypergeom")
  n_clasic <- plan_clasic$n
  c_clasic <- plan_clasic$c
  
  # 2. Configurar las Distribuciones A Priori
  alpha_beta_params <- data.frame(alpha_b = rep(NA, length(p1)), beta_b = rep(NA, length(p1)))
  for (i in 1:nrow(alpha_beta_params)) {
    Shape <- find_beta(x1 = Esce[esce, 4], p1 = p1[i], x2 = Esce[esce, 5], p2 = 1 - p2[i]) 
    alpha_beta_params[i, "alpha_b"] <- Shape$shape1
    alpha_beta_params[i, "beta_b"] <- Shape$shape2
  }
  
  # CORRECCIÓN: Añadir Prior Naive y las 3 Priors MCMC explícitamente (Total = 11 filas)
  alpha_beta_params <- rbind(alpha_beta_params, 
                             c(1, 1),               # Fila 8: Naive
                             c(a_mcmc_1, b_mcmc_1), # Fila 9: MCMC 1
                             c(a_mcmc_2, b_mcmc_2), # Fila 10: MCMC 2
                             c(a_mcmc_3, b_mcmc_3)) # Fila 11: MCMC 3
  
  # Etiquetas corregidas por usted
  esce_provee <- c("Excelente", "Bueno", "Regular", "Malo",
                   "Muy Malo", "Acum_ZI 5%", "Acum_ZI 1%",
                   "Naive", "MCMC 1 (Gamma)", "MCMC 2 (Jeffreys)", "MCMC 3 (Bivariada)")
  
  # 3. Inicializar Tabla de Resultados del Escenario
  resultados_riesgo <- data.frame(
    Escenario = 1:nrow(alpha_beta_params),
    Proveedor = esce_provee,
    n_clasic = n_clasic, c_clasic = c_clasic,
    P_Mass_Good_Beta = NA, P_Mass_Bad_Beta = NA,
    WRP_clasic = NA, WRC_clasic = NA, WRT_clasic = NA,
    n_opt = NA, c_opt = NA,
    WRP_opt = NA, WRC_opt = NA, WRT_opt = NA,
    Ahorro_n = NA, Porcentaje_Ahorro = NA
  )
  
  # 4. Evaluación de cada Proveedor
  for (j in 1:nrow(alpha_beta_params)) {
    
    alpha_b_val <- alpha_beta_params[j, "alpha_b"]
    beta_b_val <- alpha_beta_params[j, "beta_b"]
    
    # I. Masas de Probabilidad
    prob_mass_beta <- calc_prob_mass(alpha_b_val, beta_b_val, Esce[esce, 4], Esce[esce, 5])
    k_p_ <- as.numeric(prob_mass_beta["P_Good"])
    k_c_ <- as.numeric(prob_mass_beta["P_Bad"])
    
    resultados_riesgo$P_Mass_Good_Beta[j] <- k_p_
    resultados_riesgo$P_Mass_Bad_Beta[j] <- k_c_
    
    # II. Riesgos del Plan Clásico
    risks_clasic <- calc_wr(Esce[esce, 1], n_clasic, c_clasic, alpha_b_val, beta_b_val, Esce[esce, 4], Esce[esce, 5])
    resultados_riesgo$WRT_clasic[j] <- risks_clasic["WRT_val"]
    resultados_riesgo$WRP_clasic[j] <- risks_clasic["WRP_val"]
    resultados_riesgo$WRC_clasic[j] <- risks_clasic["WRC_val"]
    
    # III. Riesgo Tolerable Meta
    des_WRT <- Esce[esce, 2] * k_p_ + Esce[esce, 3] * k_c_
    
    # IV. BÚSQUEDA SECUENCIAL DETERMINÍSTICA
    n_opt_found <- NA
    c_opt_found <- NA
    
    for (n_ in 1:Esce[esce, 1]) {
      cumple <- FALSE
      for (c_ in 0:(n_ - 1)) {
        risks_opt <- calc_wr(Esce[esce, 1], n_, c_, alpha_b_val, beta_b_val, Esce[esce, 4], Esce[esce, 5])
        
        if (risks_opt["WRT_val"] <= des_WRT) {
          n_opt_found <- n_
          c_opt_found <- c_
          resultados_riesgo$WRT_opt[j] <- risks_opt["WRT_val"]
          resultados_riesgo$WRP_opt[j] <- risks_opt["WRP_val"]
          resultados_riesgo$WRC_opt[j] <- risks_opt["WRC_val"]
          cumple <- TRUE
          break
        }
      }
      if (cumple) break
    }
    
    # V. Guardar Optimización
    resultados_riesgo$n_opt[j] <- n_opt_found
    resultados_riesgo$c_opt[j] <- c_opt_found
    resultados_riesgo$Ahorro_n[j] <- n_clasic - n_opt_found
    resultados_riesgo$Porcentaje_Ahorro[j] <- round(((n_clasic - n_opt_found) / n_clasic) * 100, 2)
    
  }
  
  lista_tablas_resultados[[esce]] <- resultados_riesgo
  cat("Resultados Consolidados:\n")
  print(resultados_riesgo[, c("Proveedor", "n_clasic", "n_opt", "c_opt", "Ahorro_n", "Porcentaje_Ahorro")])
}

# # ============================================================================== #
# # SCRIPT 3: OPTIMIZACIÓN DEL PLAN DE MUESTREO (FASE II)
# # Metodología: Búsqueda Secuencial Determinística (Riesgo Acumulado Ponderado)
# # Autores Base: Fernández et al. (2020) | Adaptación: Juan Sebastián Haeusler
# # ============================================================================== #
# 
# # --- 1. Carga de Librerías ---
# paquetes <- c("pracma", "AcceptanceSampling", "sjstats")
# for (p in paquetes) {
#   if(!require(p, character.only = TRUE)) {
#     install.packages(p)
#     require(p, character.only = TRUE)
#   }
# }
# 
# # ============================================================================== #
# # --- 2. Funciones de Riesgo Ponderado y Masa de Probabilidad ---
# # ============================================================================== #
# 
# # Función de Método de Cook (Necesaria para los escenarios de expertos)
# # find_beta <- function(x1, p1, x2, p2) {
# #   loss_function <- function(params) {
# #     a <- params[1]; b <- params[2]
# #     if(a <= 0 || b <= 0) return(Inf)
# #     error <- (pbeta(x1, a, b) - p1)^2 + (pbeta(x2, a, b) - p2)^2
# #     return(error)
# #   }
# #   opt <- optim(par = c(1, 1), fn = loss_function,
# #                method = "Nelder-Mead")
# #   if(opt$convergence == 0) return(list(shape1 = opt$par[1], shape2 = opt$par[2]))
# #   else return(list(shape1 = NA, shape2 = NA))
# # }
# 
# # Probabilidad de Aceptación (PA) - Hipergeométrica
# Pa <- function(n, c, p, N){
#   Pa_val <- phyper(c, round(N * p), N - round(N * p), n)
#   return(Pa_val)
# }
# 
# # Riesgo del Productor Ponderado (RPP)
# wr_p <- function(p, n, c, alpha_b, beta_b, N_){
#   f.p <- dbeta(p, alpha_b, beta_b)
#   prod_wr <- (1 - Pa(n, c, p, N_)) * f.p
#   return(prod_wr)
# }
# 
# # Riesgo del Consumidor Ponderado (RCP)
# wr_c <- function(p, n, c, alpha_b, beta_b, N_){
#   f.p <- dbeta(p, alpha_b, beta_b)
#   cons_wr <- Pa(n, c, p, N_) * f.p
#   return(cons_wr)
# }
# 
# # Masa de Probabilidad (Zonas de Calidad)
# calc_prob_mass <- function(alpha_b, beta_b, AQL, LTPD) {
#   P_Good <- pbeta(AQL, shape1 = alpha_b, shape2 = beta_b)
#   P_Bad <- 1 - pbeta(LTPD, shape1 = alpha_b, shape2 = beta_b)
#   return(c(P_Good = P_Good, P_Bad = P_Bad))
# }
# 
# # Integración Numérica con Pracma (Método Kronrod)
# calc_wr <- function(N_, n, c, alpha_b, beta_b, AQL, LTPD) {
#   wrp_val <- integral(f = function(p) wr_p(p, n, c,
#                                            alpha_b, beta_b, N_),
#                       xmin = 0, xmax = AQL, method = "Kron")
#   wrc_val <- integral(f = function(p) wr_c(p, n, c,
#                                            alpha_b, beta_b, N_),
#                       xmin = LTPD, xmax = 1, method = "Kron")
#   wrt_val <- wrp_val + wrc_val
#   return(c(WRP_val = wrp_val, WRC_val = wrc_val,
#            WRT_val = wrt_val))
# }
# 
# # ============================================================================== #
# # --- 3. Variables Globales, Escenarios y Perfiles ---
# # ============================================================================== #
# 
# N <- 1000                      
# alpha <- c(0.01, 0.05)         
# beta <- c(0.05, 0.10, 0.20)    
# AQL <- c(0.01, 0.02, 0.05)     
# LTPD <- c(0.08, 0.10, 0.15, 0.20)    
# 
# Esce <- expand.grid(N = N, alpha = alpha, beta = beta,
#                     AQL = AQL, LTPD = LTPD)
# Esce <- Esce[Esce$AQL < Esce$LTPD, ]
# rownames(Esce) <- NULL
# 
# # Perfiles Históricos (Cook)
# p1 <- c(0.95, 0.75, 0.50, 0.25, 0.10, 0.05, 0.01) 
# p2 <- c(1e-4, 0.05, 0.10, 0.25, 0.40, 0.05, 0.01) 
# 
# # Parámetros MCMC (Traídos del Objetivo 1 - SCRIPT 1)
# a_est_mcmc_cook <- 10.1686
# b_est_mcmc_cook <- 190.0793
# 
# lista_tablas_resultados <- list()
# 
# # ============================================================================== #
# # --- 4. Motor de Búsqueda Determinística ---
# # ============================================================================== #
# 
# # NOTA: Para no saturar la consola en pruebas, evaluaremos los primeros 3 escenarios.
# # Cambie a `1:nrow(Esce)` para correr la tesis completa.
# cases <- 1:3 
# 
# for (esce in cases) {# esce <- 1 + esce
#   
#   cat(sprintf("\n--- Procesando Escenario %d de %d ---\n", esce, max(cases)))
#   cat(sprintf("AQL: %.2f | LTPD: %.2f | Alpha: %.2f | Beta: %.2f\n", 
#               Esce[esce, 4], Esce[esce, 5], Esce[esce, 2], Esce[esce, 3]))
#   
#   # 1. Plan Clásico (Naive)
#   plan_clasic <- find.plan(PRP = c(Esce[esce, 4], 1 - Esce[esce, 2]),
#                            CRP = c(Esce[esce, 5], Esce[esce, 3]),
#                            N = Esce[esce, 1], type = "hypergeom")
#   n_clasic <- plan_clasic$n
#   c_clasic <- plan_clasic$c
#   
#   # 2. Configurar las 8 Distribuciones A Priori
#   alpha_beta_params <- data.frame(alpha_b = rep(NA, length(p1)), beta_b = rep(NA, length(p1)))
#   for (i in 1:nrow(alpha_beta_params)) {
#     Shape <- find_beta(x1 = Esce[esce, 4], p1 = p1[i], x2 = Esce[esce, 5], p2 = 1 - p2[i]) 
#     alpha_beta_params[i, "alpha_b"] <- Shape$shape1
#     alpha_beta_params[i, "beta_b"] <- Shape$shape2
#   }
#   
#   # Añadir Prior Naive (Uniforme) y Prior OTB_est (MCMC Conjunta Bivariada)
#   alpha_beta_params <- rbind(alpha_beta_params, c(1, 1),
#                              c(a_est_mcmc_cook, b_est_mcmc_cook))
#   
#   esce_provee <- (c("Excelente", "Bueno", "Regular", "Malo",
#                     "Muy Malo", "Acum_ZI 5%", "Acum_ZI 1%",
#                     "Naive", paste("MCMC",1:3)))
#   
#   # 3. Inicializar Tabla de Resultados del Escenario
#   resultados_riesgo <- data.frame(
#     Escenario = 1:nrow(alpha_beta_params),
#     Proveedor = esce_provee,
#     n_clasic = n_clasic, c_clasic = c_clasic,
#     P_Mass_Good_Beta = NA, P_Mass_Bad_Beta = NA,
#     WRP_clasic = NA, WRC_clasic = NA, WRT_clasic = NA,
#     n_opt = NA, c_opt = NA,
#     WRP_opt = NA, WRC_opt = NA, WRT_opt = NA,
#     Ahorro_n = NA, Porcentaje_Ahorro = NA
#   )
#   
#   # 4. Evaluación de cada Proveedor
#   for (j in 1:nrow(alpha_beta_params)) {
#     
#     alpha_b_val <- alpha_beta_params[j, "alpha_b"]
#     beta_b_val <- alpha_beta_params[j, "beta_b"]
#     
#     # I. Masas de Probabilidad (El corazón de la propuesta de Fernández et al. 2020)
#     prob_mass_beta <- calc_prob_mass(alpha_b_val, beta_b_val, Esce[esce, 4], Esce[esce, 5])
#     k_p_ <- as.numeric(prob_mass_beta["P_Good"])
#     k_c_ <- as.numeric(prob_mass_beta["P_Bad"])
#     
#     resultados_riesgo$P_Mass_Good_Beta[j] <- k_p_
#     resultados_riesgo$P_Mass_Bad_Beta[j] <- k_c_
#     
#     # II. Riesgos del Plan Clásico evaluado bajo la lupa Bayesiana
#     risks_clasic <- calc_wr(Esce[esce, 1], n_clasic, c_clasic, alpha_b_val, beta_b_val, Esce[esce, 4], Esce[esce, 5])
#     resultados_riesgo$WRT_clasic[j] <- risks_clasic["WRT_val"]
#     resultados_riesgo$WRP_clasic[j] <- risks_clasic["WRP_val"]
#     resultados_riesgo$WRC_clasic[j] <- risks_clasic["WRC_val"]
#     
#     # III. Riesgo Tolerable Meta (Límite dinámico)
#     des_WRT <- Esce[esce, 2] * k_p_ + Esce[esce, 3] * k_c_
#     
#     # IV. BÚSQUEDA SECUENCIAL DETERMINÍSTICA
#     n_opt_found <- NA
#     c_opt_found <- NA
#     
#     for (n_ in 1:Esce[esce, 1]) {
#       cumple <- FALSE
#       for (c_ in 0:(n_ - 1)) {
#         risks_opt <- calc_wr(Esce[esce, 1], n_, c_, alpha_b_val, beta_b_val, Esce[esce, 4], Esce[esce, 5])
#         
#         # Criterio de parada: Encontrar el primer plan que cumpla el límite
#         if (risks_opt["WRT_val"] <= des_WRT) {
#           n_opt_found <- n_
#           c_opt_found <- c_
#           resultados_riesgo$WRT_opt[j] <- risks_opt["WRT_val"]
#           resultados_riesgo$WRP_opt[j] <- risks_opt["WRP_val"]
#           resultados_riesgo$WRC_opt[j] <- risks_opt["WRC_val"]
#           cumple <- TRUE
#           break
#         }
#       }
#       if (cumple) break
#     }
#     
#     # V. Guardar Optimización
#     resultados_riesgo$n_opt[j] <- n_opt_found
#     resultados_riesgo$c_opt[j] <- c_opt_found
#     resultados_riesgo$Ahorro_n[j] <- n_clasic - n_opt_found
#     resultados_riesgo$Porcentaje_Ahorro[j] <- round(((n_clasic - n_opt_found) / n_clasic) * 100, 2)
#     
#   }
#   
#   lista_tablas_resultados[[esce]] <- resultados_riesgo
#   cat("Resultados Consolidados:\n")
#   print(resultados_riesgo[, c("Proveedor", "n_clasic", "n_opt", "c_opt", "Ahorro_n", "Porcentaje_Ahorro")])
# }