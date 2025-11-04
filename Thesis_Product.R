
# Carga de librerías necesarias
library(AcceptanceSampling)
library(tidyverse)
library(dplyr)
library(tidyr)
library(gridExtra)

set.seed(123)

N <- 200  # Tamaño del lote
alpha <- 0.05 # Riesgo del productor
beta <- 0.10 # Riesgo del consumidor

AQL <- 0.05# Nivel de calidad aceptable
RQL <- 0.10# Nivel de calidad de rechazo
digs <- 5
delta_p <- 10^-digs # Tolerancia para el paso de p (1e-5)

# Espacio de calidades (p)
prop <- seq(0, 1, by = delta_p)

# Índices para las zonas de evaluación
ind_aql <- which(prop <= AQL) # Zona AQL (Aceptación)
ind_rql <- which(prop >= RQL) # Zona RQL (Rechazo)
ind_zi <- which(prop > AQL & prop < RQL) # Zona de Indiferencia (ZI)

# Parámetros de la distribución Hipergeométrica para cada p
m_hyp = N * prop
n_hyp = N * (1 - prop)

# Plan Clásico (Hipergeométrico)
plan_clasic <- find.plan(PRP = c(AQL, 1 - alpha),
                         CRP = c(RQL, beta),
                         N = N, type = "hypergeom")

n_clasic <- plan_clasic$n
c_clasic <- plan_clasic$c

# FACTOR DE ESCALA CRÍTICO (lambda)
SCALING_FACTOR_N <- (AQL / RQL) * (beta / N) 

# optimizar_loss_combinada <- function(min_p, max_p, escenario) {
  min_p <- 0
  max_p <- 0.03
  # a. Definir la Densidad a Priori (f(p)) para este escenario
  dens_eval_p <- dunif(prop, min = min_p, max = max_p)
  
  # b. Calcular la CO para el PLAN CLÁSICO
  CO_clasic <- phyper(q = c_clasic, m = m_hyp, n = n_hyp, k = n_clasic, lower.tail = TRUE)
  
  # c1. RIESGO DE ACEPTACIÓN ZI PONDERADO (Factor para Penalización 1)
  # Se usa dens_eval_p. Este es el componente que sí depende del prior.
  zi_acceptance_risk_clasic <- sum(CO_clasic[ind_zi] * dens_eval_p[ind_zi]) * delta_p
  
  # c2. ÁREA DE ACEPTACIÓN ZI NO PONDERADA (Factor para Penalización 2)
  # Este es el componente que NO depende del prior (propiedad intrínseca del plan).
  zi_acceptance_area_clasic <- sum(CO_clasic[ind_zi]) * delta_p

  # e. Riesgos del plan clásico (para referencia bajo esta distribución prior)
  rp_clasic <- sum((1 - CO_clasic[ind_aql]) * dens_eval_p[ind_aql]) * delta_p
  rc_clasic <- sum(CO_clasic[ind_rql] * dens_eval_p[ind_rql]) * delta_p
  
  # Penalización 1 (ZI-Severidad)
  penalizacion_clasic_zi_severidad <- SCALING_FACTOR_N * (n_clasic / (c_clasic + 1)) * zi_acceptance_risk_clasic
  
  # Penalización 2 (Costo Estructural Intrínseco) - ¡USANDO EL ÁREA NO PONDERADA!
  
  # densidad total en las zonas de riesgo
  total_density_p <- sum(dens_eval_p[ind_aql]*delta_p) + sum(dens_eval_p[ind_zi]*delta_p)
  total_density_c <- sum(dens_eval_p[ind_rql]*delta_p) + sum(dens_eval_p[ind_zi]*delta_p)
  
  
  penalizacion_clasic_costo <- SCALING_FACTOR_N*(n_clasic + c_clasic)*(sum(dens_eval_p[ind_aql]*delta_p) + sum(dens_eval_p[ind_rql]*delta_p))
  
  L_clasico <- (rp_clasic + rc_clasic) + penalizacion_clasic_zi_severidad + penalizacion_clasic_costo
  RAT_Clasico_Base <- rp_clasic + rc_clasic
  
  # INICIALIZACIÓN CRÍTICA
  L_minimo_global <- L_clasico
  
  n_optimo <- NA
  c_optimo <- NA
  RAT_optimo_final <- NA
  
  found_better_plan <- FALSE 
  
  # Búsqueda: Se detiene en la PRIMERA combinación (n, c) donde L_actual < L_clasico
  for (n_actual in 1:n_clasic) {
    for (c_actual in 0:(n_actual - 1)) {
      
      # 1. Calcular Riesgo Total Esperado (RAT)
      CO_actual <- phyper(q = c_actual, m = m_hyp, n = n_hyp, k = n_actual, lower.tail = TRUE)
      rp_actual <- sum((1 - CO_actual[ind_aql]) * dens_eval_p[ind_aql]) * delta_p
      rc_actual <- sum(CO_actual[ind_rql] * dens_eval_p[ind_rql]) * delta_p
      RAT_actual <- rp_actual + rc_actual
      
      # Riesgo ZI actual PONDERADO (para Penalización 1)
      zi_acceptance_risk_actual <- sum(CO_actual[ind_zi] * dens_eval_p[ind_zi]) * delta_p
      
      # Área ZI actual NO PONDERADA (para Penalización 2)
      zi_acceptance_area_actual <- sum(CO_actual[ind_zi]) * delta_p
      
      # 2. Calcular la Función de Pérdida Combinada (L)
      
      # Penalización 1: Severidad del Plan Actual * Riesgo ZI Actual * lambda
      penalizacion_actual_zi_severidad <- SCALING_FACTOR_N * (n_actual / (c_actual + 1)) * zi_acceptance_risk_actual
      
      # Penalización 2: Costo Estructural Intrínseco (USANDO EL ÁREA NO PONDERADA)
      penalizacion_actual_costo <- SCALING_FACTOR_N*(n_actual + c_actual) * (sum(dens_eval_p[ind_aql]*delta_p) + sum(dens_eval_p[ind_rql]*delta_p))
      
      L_actual <- RAT_actual + penalizacion_actual_zi_severidad + penalizacion_actual_costo
      
      # 3. Condición de Parada: Encontrar la primera mejora respecto al clásico
      if (L_actual < L_clasico) {
        L_minimo_global <- L_actual
        n_optimo <- n_actual
        c_optimo <- c_actual
        RAT_optimo_final <- RAT_actual
        found_better_plan <- TRUE
        break # Sale del bucle de c_actual (Búsqueda interna)
      }
    }
    
    if (found_better_plan) {
      break # Sale del bucle de n_actual (Búsqueda externa)
    }
  }
  
  # Devolver componentes desglosados (Penalizacion_Clasica solo incluye P1 para simplificar la tabla)
  return(data.frame(
    Escenario = escenario,
    n_Clasico = n_clasic,
    c_Clasico = c_clasic,
    L_Clasico = L_clasico,
    RAT_Clasico = RAT_Clasico_Base,
    Penalizacion_Clasica = penalizacion_clasic_zi_severidad,
    n_Optimo_Combinado = n_optimo,
    c_Optimo_Combinado = c_optimo,
    RAT_Optimo_Combinado = RAT_optimo_final,
    L_Optimo_Combinado = L_minimo_global
  ))
}

#-------------------------------------------------------------------------------
# --- 4. Ejecución para los 5 Escenarios de Densidad a Priori ---
#-------------------------------------------------------------------------------

message("\n-----------------------------------------------------")
message("--- 4. Ejecutando la Optimización para 5 Escenarios ---")
message("-----------------------------------------------------")

# Definición de los 5 escenarios (Distribuciones uniformes)
escenarios <- list(
  # E1: Calidad Excelente (Riesgo bajo)
  list(min_p = 0.00, max_p = 0.045, nombre = "E1 - Excelente"),
  # E2: Calidad Buena (Riesgo medio-bajo)
  list(min_p = 0.03, max_p = 0.07, nombre = "E2 - Bueno"),
  # E3: Calidad Regular (Cubre ZI, Riesgo medio-alto)
  list(min_p = 0.04, max_p = 0.12, nombre = "E3 - Regular"),
  # E4: Calidad Mala (Riesgo alto)
  list(min_p = 0.075, max_p = 0.13, nombre = "E4 - Malo"),
  # E5: Calidad Muy Mala (Riesgo muy alto, solo rechazo)
  list(min_p = 0.15, max_p = 0.20, nombre = "E5 - Muy Malo")
)

resultados_optimizacion <- do.call(rbind, lapply(escenarios, function(e) {
  optimizar_loss_combinada(e$min_p, e$max_p, e$nombre)
}))

# Limpieza y formateo de la tabla
resultados_optimizacion <- resultados_optimizacion %>%
  mutate(
    L_Clasico = round(L_Clasico, 6),
    RAT_Clasico = round(RAT_Clasico, 6),
    Penalizacion_Clasica = round(Penalizacion_Clasica, 6),
    RAT_Optimo_Combinado = round(RAT_Optimo_Combinado, 6),
    L_Optimo_Combinado = round(L_Optimo_Combinado, 6)
  ) %>%
  select(Escenario, n_Clasico, c_Clasico, L_Clasico, RAT_Clasico, Penalizacion_Clasica,
         n_Optimo_Combinado, c_Optimo_Combinado, L_Optimo_Combinado, RAT_Optimo_Combinado)

message("\n-----------------------------------------------------")
message("--- RESULTADOS FINALES DE LA OPTIMIZACIÓN L(n, c) ---")
message("-----------------------------------------------------")
print(resultados_optimizacion)
