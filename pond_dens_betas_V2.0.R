# Script para comparar planes de muestreo de aceptación clásicos
# con planes óptimos basados en el Riesgo Ponderado (Enfoque Bayesiano).

# Las prioris de calidad de los proveedores se definen por pares (PA en AQL, PA en LTPD)
# Excelente proveedor: PA(AQL)=0.95, PA(LTPD)=1-0.0001
# Buen proveedor: PA(AQL)=0.75, PA(LTPD)=1-0.05
# Regular: PA(AQL)=0.50, PA(LTPD)=1-0.10
# Malo: PA(AQL)=0.25, PA(LTPD)=1-0.25
# Muy malo: PA(AQL)=0.10, PA(LTPD)=1-0.40

# --- Carga de Librerías ---
# install.packages("pracma") # Si no lo tienes instalado
library(pracma)
# install.packages("sjstats")
library(sjstats)
library(AcceptanceSampling)

# --- Variables Globales y Parámetros del Problema ---
N <- 200        # Tamaño del lote
alpha <- 0.05   # Riesgo del productor (para el plan clásico)
beta <- 0.10    # Riesgo del consumidor (para el plan clásico)
AQL <- 0.05     # Nivel de Calidad Aceptable
LTPD <- 0.10    # Tolerancia de Porcentaje Defectuoso de Lote


# LAMBDA <- (alpha/AQL) * (beta/LTPD)

# Vectores de probabilidades para los 5 escenarios
# p1 = Probabilidad de Aceptación en el AQL
# p2 = Probabilidad de No Aceptación (o Rechazo) en el LTPD
p1 <- c(0.95, 0.75, 0.50, 0.25, 0.10)
p2 <- c(1e-4, 0.05, 0.10, 0.25, 0.40) # 1 - PA(LTPD) -> Beta (Riesgo del Consumidor)

# --- Funciones de Riesgo Ponderado y Masa de Probabilidad ---

# Probabilidad de Aceptación (PA) - Curva CO tipo A (Hypergeométrica)
Pa <- function(n, c, p, N){
  # Se utiliza phyper para la distribución hipergeométrica (muestreo sin reemplazo)
  Pa_val <- phyper(c, N * p, N * (1 - p), n)
  return(Pa_val)
}

# Integrando el Riesgo del Productor (RP): (1 - PA) * f(p)
# Error tipo I (Rechazar lote bueno): p en [0, AQL]
wr_p <- function(p, n, c, alpha_b, beta_b, N){
  f.p <- dbeta(p, alpha_b, beta_b) # Densidad a priori Beta
  prod_wr <- (1 - Pa(n, c, p, N)) * f.p
  return(prod_wr)
}

# Integrando el Riesgo del Consumidor (RC): PA * f(p)
# Error tipo II (Aceptar lote malo): p en [LTPD, 1]
wr_c <- function(p, n, c, alpha_b, beta_b, N){
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
calc_wr <- function(n, c, alpha_b, beta_b, AQL, LTPD) {

  LAMBDA <- ((n / (1 + c)) * (1 / N)) * (1/(alpha * beta))
  
  # 1. Calcular RAT (RAT_val)
  rp_val <- integral(f = function(p) wr_p(p, n, c, alpha_b, beta_b, N), xmin = 0, xmax = AQL, method = "Kron")
  rc_val <- integral(f = function(p) wr_c(p, n, c, alpha_b, beta_b, N), xmin = LTPD, xmax = 1, method = "Kron")
  RAT_val <- rp_val + rc_val
  
  # 2. Calcular Penalty 1 (ZI)
  zi_val <- integral(f = function(p) wr_c(p, n, c, alpha_b, beta_b, N), xmin = AQL, xmax = LTPD, method = "Kron")
  Penalty1_val <- LAMBDA * zi_val
  
  # 3. Calcular Penalty 2 (Costo intrínseco)
  prob_masses <- calc_prob_mass(alpha_b = alpha_b, beta_b = beta_b, AQL, LTPD)
  P_Good_Bad_Sum <- prob_masses["P_Good"] + prob_masses["P_Bad"]
  Penalty2_val <- LAMBDA * P_Good_Bad_Sum
  
  # Calcular L(n,c) total
  L_val <- RAT_val + Penalty1_val + Penalty2_val
  
  return(c(RP_val = rp_val, RC_val = rc_val, L_val = L_val))
}

# 1. Determinar el Plan Clásico (basado en alpha y beta fijos)
plan_clasic <- find.plan(PRP = c(AQL, 1 - alpha),
                         CRP = c(LTPD, beta),
                         N = N, type = "hypergeom")

n_clasic <- plan_clasic$n
c_clasic <- plan_clasic$c

# --- CÁLCULO DE LÍNEA BASE BAJO IGNORANCIA (UNIFORME Beta(1, 1)) ---
# Calculamos los riesgos individuales RP y RC para el Plan Clásico bajo la densidad Uniforme
# ESTE ES EL RIESGO DE REFERENCIA: Plan Clásico (Naive)
risks_uniforme <- calc_wr(n = n_clasic, c = c_clasic, 
                          alpha_b = 1, beta_b = 1, # Uniforme U(0,1)
                          AQL, LTPD)
RP_U01_Clasic <- risks_uniforme["RP_val"]
RC_U01_Clasic <- risks_uniforme["RC_val"]
L_Clasic <- risks_uniforme["L_val.P_Good"]

# Y la masa de probabilidad (densidad acumulada) para el prior uniforme
prob_mass_uniforme <- calc_prob_mass(alpha_b = 1, beta_b = 1, AQL, LTPD)
P_Mass_Good_Naive <- prob_mass_uniforme["P_Good"]
P_Mass_Bad_Naive <- prob_mass_uniforme["P_Bad"]
# -------------------------------------------------------------------

# 2. Obtener los 5 pares de (alpha_b, beta_b) para la distribución Beta
alpha_beta_params <- data.frame(alpha_b = rep(NA, 5), beta_b = rep(NA, 5))

for (i in 1:5) {
  # find_beta encuentra los parámetros de la distribución Beta
  Shape <- find_beta(x1 = AQL, p1 = p1[i], x2 = LTPD, p2 = 1 - p2[i]) 
  alpha_beta_params[i, "alpha_b"] <- Shape$shape1
  alpha_beta_params[i, "beta_b"] <- Shape$shape2
}

# 3. Inicializar la tabla de resultados
resultados_riesgo <- data.frame(
  Escenario = 1:5,
  Proveedor = c("Excelente", "Bueno", "Regular", "Malo", "Muy Malo"),
  n_clasic = n_clasic,
  c_clasic = c_clasic,
  # Densidad Acumulada (Masa de Probabilidad) del Prior
  P_Mass_Good_Naive = P_Mass_Good_Naive, # P(p < AQL | Naive)
  P_Mass_Bad_Naive = P_Mass_Bad_Naive,   # P(p > LTPD | Naive)
  P_Mass_Good_Beta = NA,                 # P(p < AQL | Beta)
  P_Mass_Bad_Beta = NA,                  # P(p > LTPD | Beta)
  # Riesgos de la Línea Base (Plan Clásico, Densidad Uniforme)
  RP_U01_Clasic = RP_U01_Clasic, 
  RC_U01_Clasic = RC_U01_Clasic,
  L_Clasic = L_Clasic, # Mantenemos TWR para el cálculo intermedio
  # Riesgos bajo Información (Plan Clásico - Solo para referencia)
  RP_clasic = NA,
  RC_clasic = NA,
  L_clasic = NA, 
  # Plan Óptimo (Minimiza TWR bajo Información)
  n_opt = NA,
  c_opt = NA,
  RP_opt = NA, 
  RC_opt = NA,
  L_opt = NA, 
  # Ganancias (Plan Clásico Naive vs. Plan Óptimo Informado)
  RP_Ganancia = NA, 
  RC_Ganancia = NA,
  L_Ganancia = NA
)

# 4. Iterar sobre los 5 escenarios y calcular los riesgos
for (i in 1:5) {
  
  # Usar los valores calculados de alpha y beta específicos para el proveedor i
  alpha_b_val <- alpha_beta_params[i, "alpha_b"]
  beta_b_val <- alpha_beta_params[i, "beta_b"]
  
  # I. Calcular la masa de probabilidad para el Prior Beta específico (Informado)
  prob_mass_beta <- calc_prob_mass(alpha_b = alpha_b_val, beta_b = beta_b_val, AQL, LTPD)
  resultados_riesgo[i, "P_Mass_Good_Beta"] <- prob_mass_beta["P_Good"]
  resultados_riesgo[i, "P_Mass_Bad_Beta"] <- prob_mass_beta["P_Bad"]
  
  # II. Calcular los riesgos para el PLAN CLÁSICO (Bajo el Prior Beta ESPECÍFICO)
  risks_clasic <- calc_wr(n = n_clasic, c = c_clasic, 
                          alpha_b = alpha_b_val, beta_b = beta_b_val, 
                          AQL, LTPD)
  
  resultados_riesgo[i, "L_clasic"] <- risks_clasic["L_val.P_Good"]
  resultados_riesgo[i, "RP_clasic"] <- risks_clasic["RP_val"]
  resultados_riesgo[i, "RC_clasic"] <- risks_clasic["RC_val"]
  
  # --- III. Búsqueda del Plan Óptimo (Minimizar TWR Puro) para el Escenario i ---
  
  n_opt_found <- NA
  c_opt_found <- NA
  min_L <- Inf 
  
  # Búsqueda exhaustiva: Itera n de 1 hasta n_clasic (máximo tamaño de muestra del plan clásico)
  for (n_ in 1:n_clasic) { 
    for (c_ in 0:(n_ - 1)) {
      
      risks_opt <- calc_wr(n = n_, c = c_, 
                           alpha_b = alpha_b_val, beta_b = beta_b_val, 
                           AQL, LTPD)
      
      L_current <- risks_opt["L_val.P_Good"]
      
      # *** CRITERIO DE SELECCIÓN DEL PLAN ÓPTIMO (n_opt, c_opt) ***
      # Se utiliza el Riesgo Total Ponderado (TWR) para encontrar el único plan 
      # que ofrece el mejor balance entre el Riesgo del Productor (RP) y el 
      # Riesgo del Consumidor (RC). 
      # NOTA: Minimizar TWR es necesario para encontrar UN solo plan.
      # La comparación separada de RP y RC se hace en la tabla final (RP_Ganancia, RC_Ganancia).
      if (L_current < min_L) {
        min_L <- L_current
        n_opt_found <- n_
        c_opt_found <- c_
      }
    }
  }
  
  # IV. Almacenar los resultados del plan óptimo (n_opt, c_opt)
  resultados_riesgo[i, "n_opt"] <- n_opt_found
  resultados_riesgo[i, "c_opt"] <- c_opt_found
  
  # Recalcular los riesgos para el análisis final
  if (!is.na(n_opt_found)) {
    risks_opt_final <- calc_wr(n = n_opt_found, c = c_opt_found, 
                               alpha_b = alpha_b_val, beta_b = beta_b_val, 
                               AQL, LTPD) 
    
    resultados_riesgo[i, "L_opt"] <- risks_opt_final["L_val.P_Good"]
    resultados_riesgo[i, "RP_opt"] <- risks_opt_final["RP_val"]
    resultados_riesgo[i, "RC_opt"] <- risks_opt_final["RC_val"]
    
    # Cálculo de las ganancias (Plan Clásico Naive vs. Óptimo con Información)
    # Ganancia RP = Riesgo RP Naive - Riesgo RP Óptimo
    # Ganancia RC = Riesgo RC Naive - Riesgo RC Óptimo
    # resultados_riesgo[i, "RP_Ganancia"] <- resultados_riesgo[i, "RP_U01_Clasic"] - resultados_riesgo[i, "RP_opt"]
    # resultados_riesgo[i, "RC_Ganancia"] <- resultados_riesgo[i, "RC_U01_Clasic"] - resultados_riesgo[i, "RC_opt"]
    # resultados_riesgo[i, "L_Ganancia"] <- resultados_riesgo[i, "L_U01_Clasic"] - resultados_riesgo[i, "L_opt"]

  } else {
    resultados_riesgo[i, "L_opt"] <- NA
    resultados_riesgo[i, "RP_opt"] <- NA
    resultados_riesgo[i, "RC_opt"] <- NA
    resultados_riesgo[i, "RP_Ganancia"] <- NA
    resultados_riesgo[i, "RC_Ganancia"] <- NA
    resultados_riesgo[i, "L_Ganancia"] <- NA
  }
}

# 5. Mostrar la tabla de resultados final
options(digits = 6) 
print("Tabla de Comparación de Riesgos Individuales (RP y RC): Naive vs. Bayesiano Óptimo")
print("------------------------------------------------------------------------------------")

# SELECCIÓN Y REORDENAMIENTO: Se incluyen las columnas de masa de probabilidad a priori
resultados_final_comparado <- resultados_riesgo[, c(
  "Escenario", "Proveedor", 
  "n_clasic", "c_clasic", 
  "n_opt", "c_opt", 
  "P_Mass_Good_Naive", "P_Mass_Bad_Naive", 
  "P_Mass_Good_Beta", "P_Mass_Bad_Beta",
  "RP_U01_Clasic", "RP_opt", "RP_Ganancia",
  "RC_U01_Clasic", "RC_opt", "RC_Ganancia"
)]

print(resultados_final_comparado)
