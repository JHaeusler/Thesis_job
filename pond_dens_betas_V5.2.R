# Propuesta para diseño de planes de muestreo simple para atributos
# considerando el histórico de la proporción estimada de unidades
# defectuosas en la inspección de lotes

# Calidad del proveedor según el histórico
# Las prioris de calidad de los proveedores se definen por pares (PA en AQL, PA en LTPD)
# Excelente proveedor: PA(AQL)=0.95, PA(LTPD)=1-0.0001
# Buen proveedor: PA(AQL)=0.75, PA(LTPD)=1-0.05
# Regular: PA(AQL)=0.50, PA(LTPD)=1-0.10
# Malo: PA(AQL)=0.25, PA(LTPD)=1-0.25
# Muy malo: PA(AQL)=0.10, PA(LTPD)=1-0.40

# --- Carga de Librerías ---
# install.packages("pracma") # Si no lo tienes instalado
# install.packages("sjstats")

library(pracma)
library(sjstats)
library(AcceptanceSampling)
library(GA)

# --- Funciones de Riesgo Ponderado y Masa de Probabilidad ---

# Función para calcular la masa de probabilidad (densidad acumulada) en las regiones de riesgo
calc_prob_mass <- function(alpha_b, beta_b, AQL, LTPD) {
  P_Good <- pbeta(AQL, shape1 = alpha_b, shape2 = beta_b)
  P_Bad <- 1 - pbeta(LTPD, shape1 = alpha_b, shape2 = beta_b)
  return(c(P_Good = P_Good, P_Bad = P_Bad))
}

# Probabilidad de Aceptación (PA) - Curva CO tipo A (Hypergeométrica)
Pa <- function(n, c, p, N){
  # Se utiliza phyper para la distribución hipergeométrica (muestreo sin reemplazo)
  Pa_val <- phyper(c, N * p, N * (1 - p), n)
  return(Pa_val)
}

# Integrando el Riesgo del Productor (RP): (1 - PA) * f(p)
# Error tipo I (Rechazar lote bueno): p en [0, AQL]
wr_p <- function(p, n, c, alpha_b, beta_b, N_){
  f.p <- dbeta(p, alpha_b, beta_b) # Densidad del historico de calidad
  prod_wr <- (1 - Pa(n, c, p, N)) * f.p
  return(prod_wr)
}

# Integrando el Riesgo del Consumidor (RC): PA * f(p)
# Error tipo II (Aceptar lote malo): p en [LTPD, 1]
wr_c <- function(p, n, c, alpha_b, beta_b, N_){
  f.p <- dbeta(p, alpha_b, beta_b) # Densidad a priori Beta
  cons_wr <- Pa(n, c, p, N) * f.p
  return(cons_wr)
}

# Función para calcular la masa de probabilidad (densidad acumulada) en las regiones de riesgo
calc_prob_mass <- function(alpha_b, beta_b, AQL, LTPD) {
  P_Good <- pbeta(AQL, shape1 = alpha_b, shape2 = beta_b)
  P_Bad <- 1 - pbeta(LTPD, shape1 = alpha_b, shape2 = beta_b)
  return(c(P_Good = P_Good, P_Bad = P_Bad))
}

# Función principal para calcular el Riesgo Ponderado Integrado (RP y RC)
calc_wr <- function(N_, n, c, alpha_b, beta_b, AQL, LTPD, k_p, k_c) {
  
  # 1. Calcular RAT (RAT_val)
  wrp_val <- k_p*integral(f = function(p) wr_p(p, n, c, alpha_b, beta_b, N_), xmin = 0, xmax = AQL, method = "Kron")
  wrc_val <- k_c*integral(f = function(p) wr_c(p, n, c, alpha_b, beta_b, N_), xmin = LTPD, xmax = 1, method = "Kron")
  wrt_val <- wrp_val + wrc_val
  
  return(c(WRP_val = wrp_val, WRC_val = wrc_val, WRT_val = wrt_val))
}

decode <- function(string){
  string <- gray2binary(string)
  n <- binary2decimal(string[1:l1])
  c <- min(n, binary2decimal(string[(l1 + 1):(l1 + l2)]))
  
  return(c(n,c))
}

fitness <- function(string){
  par <- decode(string)
  n <- par[1]
  c <- par[2]
  Pa_p <- phyper(c, Esce[esce, 1]*Esce[esce, 4], Esce[esce, 1]*(1 - Esce[esce, 4]), n)
  Pa_c <- phyper(c, Esce[esce, 1]*Esce[esce, 5], Esce[esce, 1]*(1 - Esce[esce, 5]), n)
  Loss <- (Pa_p - (1 - Esce[esce, 2]))^2 + (Pa_c - Esce[esce, 3])^2
  -Loss
}

# --- Variables Globales y Parámetros del Problema ---

N <- 1000 #c(200, 1000)        # Tamaño del lote
alpha <- c(0.01, 0.05)   # Riesgo del productor (para el plan clásico)
beta <- c(0.05, 0.10, 0.20)    # Riesgo del consumidor (para el plan clásico)
AQL <- c(0.01, 0.02, 0.05)     # Nivel de Calidad Aceptable
LTPD <- c(0.08, 0.10, 0.15, 0.20)    # Tolerancia de Porcentaje Defectuoso de Lote

Esce <- expand.grid(N = N, alpha = alpha, beta = beta, AQL = AQL, LTPD = LTPD)
Esce <- Esce[Esce$AQL < Esce$LTPD,]

# Vectores de probabilidades para los 5 escenarios
# p1 = Probabilidad de Aceptación en el AQL
# p2 = Probabilidad de No Aceptación (o Rechazo) en el LTPD
p1 <- c(0.95, 0.75, 0.50, 0.25, 0.10) #, 0.1, 0.3)
p2 <- c(1e-4, 0.05, 0.10, 0.25, 0.40) #, 0.89, 0.69) # 1 - PA(LTPD) -> Beta (Riesgo del Consumidor)

# Inicializar las listas maestras
lista_tablas_resultados <- list()
lista_parametros_beta   <- list()

cases <- seq(1, dim(Esce)[1], 6)

for (esce in cases) { # esce <- 1 + esce
  
  # 1. Determinar el Plan Clásico (basado en alpha y beta fijos)
  plan_clasic <- find.plan(PRP = c(Esce[esce, 4], 1 - Esce[esce, 2]),
                           CRP = c(Esce[esce, 5], Esce[esce, 3]),
                           N = Esce[esce, 1], type = "hypergeom")
  
  n_clasic <- plan_clasic$n
  c_clasic <- plan_clasic$c
  
  # n_ran <- 2:Esce[esce, 1]
  # c_ran <- 0:(max(n_ran) - 1)
  # 
  # b1 <- decimal2binary(max(n_ran)); l1 <- length(b1)
  # b2 <- decimal2binary(max(c_ran)); l2 <- length(b2)
  # 
  # plan_genetico <- ga(type = "binary", nBits = l1 + l2,
  #                     fitness = fitness, popSize = 200,
  #                     maxiter = 200, run = 100, seed = 060722)
  # 
  # plan_ga <- decode(plan_genetico@solution)
  # 
  # n_ga <- plan_ga[1]
  # c_ga <- plan_ga[2]
  
  # 2. Obtener los pares de (alpha_b, beta_b) para la distribución Beta
  alpha_beta_params <- data.frame(alpha_b = rep(NA, length(p1)), beta_b = rep(NA, length(p1)))
  
  for (i in 1:dim(alpha_beta_params)[1]) { # i <- 1 + i
    
    # find_beta encuentra los parámetros de la distribución Beta
    Shape <- find_beta(x1 = Esce[esce, 4], p1 = p1[i], x2 = Esce[esce, 5], p2 = 1 - p2[i]) 
    alpha_beta_params[i, "alpha_b"] <- Shape$shape1
    alpha_beta_params[i, "beta_b"] <- Shape$shape2
  }
  alpha_beta_params <- rbind(alpha_beta_params, c(1, 1), c(4.775, 75.553))
  
  # 3. Inicializar la tabla de resultados
  resultados_riesgo <- data.frame(
    Escenario = 1:dim(alpha_beta_params)[1],
    Proveedor = c("Excelente", "Bueno", "Regular", "Malo", "Muy Malo", "Naive", "OTB_est"),# paste("otro", 1:(abs(5-length(p1))))),
    n_clasic = n_clasic,
    c_clasic = c_clasic,
    
    n_ga = NA, #n_ga,
    c_ga = NA, #c_ga,
    
    P_Mass_Good_Beta = NA,                 # P(p < AQL | Beta)
    P_Mass_Bad_Beta = NA,                  # P(p > LTPD | Beta)
    
    WRP_clasic = NA,
    WRC_clasic = NA,
    WRT_clasic = NA,
    # Plan Óptimo (Minimiza TWR bajo Información)
    n_opt = NA,
    c_opt = NA,
    WRP_opt = NA,
    WRC_opt = NA,
    WRT_opt = NA,
    
    # Ganancias (Plan Clásico Naive vs. Plan Óptimo Informado)
    WRP_Ganancia = NA, 
    WRC_Ganancia = NA,
    WRT_Ganancia = NA
  )
  
  # 4. Iterar sobre los escenarios y calcular los riesgos
  
  for (j in 1:dim(alpha_beta_params)[1]) { # j <-6 1 + j
    
    # Usar los valores calculados de alpha y beta específicos para el proveedor i
    alpha_b_val <- alpha_beta_params[j, "alpha_b"]
    beta_b_val <- alpha_beta_params[j, "beta_b"]
    
    # I. Calcular la masa de probabilidad para el Prior Beta específico (Informado)
    prob_mass_beta <- calc_prob_mass(alpha_b = alpha_b_val, beta_b = beta_b_val, Esce[esce, 4], Esce[esce, 5])
    resultados_riesgo[j, "P_Mass_Good_Beta"] <- prob_mass_beta["P_Good"]
    resultados_riesgo[j, "P_Mass_Bad_Beta"] <- prob_mass_beta["P_Bad"]
    
    k_p_ <- as.numeric(prob_mass_beta["P_Good"])
    k_c_ <- as.numeric((prob_mass_beta["P_Bad"]))
    
    # II. Calcular los riesgos para el PLAN CLÁSICO (Bajo el Prior Beta ESPECÍFICO)
    risks_clasic <- calc_wr(N_ = Esce[esce, 1], n = n_clasic, c = c_clasic, 
                            alpha_b = alpha_b_val, beta_b = beta_b_val, 
                            Esce[esce, 4], Esce[esce, 5], k_p = k_p_, k_c = k_c_)
    
    resultados_riesgo[j, "WRT_clasic"] <- risks_clasic["WRT_val"]
    resultados_riesgo[j, "WRP_clasic"] <- risks_clasic["WRP_val"]
    resultados_riesgo[j, "WRC_clasic"] <- risks_clasic["WRC_val"]
    
    
    # --- III. Búsqueda del Plan Óptimo (Minimizar TWR Puro) para el Escenario i ---
    
    dens_a <- dbeta(Esce[esce, 4], alpha_b_val, beta_b_val)
    dens_b <- dbeta(Esce[esce, 5], alpha_b_val, beta_b_val)
    
    # Ancho para el riesgo del productor: AQL - 0
    # Ancho para el riesgo del consumidor: 1 - LTPD
    des_WRT <- k_p_ * ((Esce[esce, 2]) * dens_a * (Esce[esce, 4]/8)) + 
      (k_c_ * Esce[esce, 3] * dens_b * ((1 - Esce[esce, 5])/8))
    
    n_opt_found <- NA
    c_opt_found <- NA
    cumple <- FALSE
    
    # Búsqueda exhaustiva
    for (n_ in 1:Esce[esce, 1]) { # n_ <- 200 + n_
      for (c_ in 0:(n_ - 1)) { # c_ <- 5 + 1 + c_
        
        risks_opt <- calc_wr(N_ = Esce[esce, 1], n = n_, c = c_, 
                             alpha_b = alpha_b_val, beta_b = beta_b_val, 
                             Esce[esce, 4], Esce[esce, 5], k_p = k_p_, k_c = k_c_)
        
        WRP_current <- risks_opt["WRP_val"]
        WRC_current <- risks_opt["WRC_val"]
        WRT_current <- risks_opt["WRT_val"]
        
        # Usando la lógica de tu código: buscar un plan 'mejor' que el plan Naive
        if (WRT_current <= des_WRT) {
          n_opt_found <- n_
          c_opt_found <- c_
          cumple <- TRUE
          break
        }
      }
      if(cumple) break
    } 
    
    
    # IV. Almacenar los resultados del plan óptimo (n_opt, c_opt)
    resultados_riesgo[j, "n_opt"] <- n_opt_found
    resultados_riesgo[j, "c_opt"] <- c_opt_found
    
    # Recalcular los riesgos para el análisis final
    if (!is.na(n_opt_found)) {
      risks_opt_final <- calc_wr(N_ = Esce[esce, 1], n = n_opt_found, c = c_opt_found, 
                                 alpha_b = alpha_b_val, beta_b = beta_b_val, 
                                 Esce[esce, 4], Esce[esce, 5], k_p = k_p_, k_c = k_c_) 
      
      resultados_riesgo[j, "WRT_opt"] <- risks_opt_final["WRT_val"]
      resultados_riesgo[j, "WRP_opt"] <- risks_opt_final["WRP_val"]
      resultados_riesgo[j, "WRC_opt"] <- risks_opt_final["WRC_val"]
      
      # V. ¡Cálculo de Ganancias Añadido!
      # Ganancia = Riesgo Naive - Riesgo Óptimo
      resultados_riesgo[j, "WRP_Ganancia"] <- resultados_riesgo[j, "WRP_clasic"] - resultados_riesgo[j, "WRP_opt"]
      resultados_riesgo[j, "WRC_Ganancia"] <- resultados_riesgo[j, "WRC_clasic"] - resultados_riesgo[j, "WRC_opt"]
      resultados_riesgo[j, "WRT_Ganancia"] <- resultados_riesgo[j, "WRT_clasic"] - resultados_riesgo[j, "WRT_opt"]
      
      
    } else {
      resultados_riesgo[j, "WRT_opt"] <- NA
      resultados_riesgo[j, "WRP_opt"] <- NA
      resultados_riesgo[j, "WRC_opt"] <- NA
      
      resultados_riesgo[j, "WRP_Ganancia"] <- NA
      resultados_riesgo[j, "WRC_Ganancia"] <- NA
      resultados_riesgo[j, "WRT_Ganancia"] <- NA
    }
  }
  
  
  # Guardar la tabla de resultados del escenario actual en la lista maestra
  lista_tablas_resultados[[esce]] <- resultados_riesgo
  
  # Guardar los parámetros de las densidades Beta calculadas para este escenario
  lista_parametros_beta[[esce]] <- alpha_beta_params
  
}