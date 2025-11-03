# ==============================================================================
# SCRIPT DE OPTIMIZACIÓN BAYESIANA DE FASE ÚNICA CON FUNCIÓN DE PÉRDIDA COMBINADA
# OBJETIVO: Minimizar L = RAT + lambda * (n / (c + 1)) * P_ZI_Aceptacion_Clasica
# CALIBRACIÓN ACTUAL: Se utiliza la fórmula paramétrica sin el factor alpha. 
# ==============================================================================

# Instalación y carga de librerías necesarias
lista_de_paquetes <- c("AcceptanceSampling", "tidyverse", "dplyr", "tidyr")

for (paquete in lista_de_paquetes) {
  if(!require(paquete, character.only = TRUE, quietly = TRUE)){
    message(paste("Instalando paquete:", paquete))
    install.packages(paquete)
    require(paquete, character.only = TRUE, quietly = TRUE)
  }
}

# Establecer semilla para reproducibilidad
set.seed(123)

message("-----------------------------------------------------")
message("--- 1. Definición de Parámetros y Constantes ---")
message("-----------------------------------------------------")

N <- 200      # Tamaño del lote
alpha <- 0.05 # Riesgo del productor
beta <- 0.10  # Riesgo del consumidor

AQL <- 0.05   # Nivel de calidad aceptable
RQL <- 0.10   # Nivel de calidad de rechazo
digs <- 5
delta_p <- 10^-digs

# Espacio de calidades (p)
prop <- seq(0, 1, by = delta_p)

# Índices para las zonas de evaluación
ind_aql <- which(prop <= AQL)         # Zona AQL (Aceptación)
ind_rql <- which(prop >= RQL)         # Zona RQL (Rechazo)
ind_zi <- which(prop > AQL & prop < RQL) # Zona de Indiferencia (ZI)

# Parámetros de la distribución Hipergeométrica para cada p
m_hyp = N * prop
n_hyp = N * (1 - prop)

#-------------------------------------------------------------------------------
# --- 2. Cálculo del plan de muestreo clásico (Referencia) ---
#---------------------------------------------------------------
message("2. Calculando el plan de muestreo clásico...")

plan_clasic <- find.plan(PRP = c(AQL, 1 - alpha),
                         CRP = c(RQL, beta),
                         N = N, type = "hypergeom")

n_clasic <- plan_clasic$n 
c_clasic <- plan_clasic$c
message(paste("Plan Clásico de Referencia: n =", n_clasic, ", c =", c_clasic))

# FACTOR DE ESCALA CRÍTICO (lambda)
# Fórmula base paramétrica modificada: (AQL/RQL) * (beta / N)
lambda_base <- (AQL / RQL) * (beta / N)

# SCALING_FACTOR_N es igual al lambda_base (Sin alpha)
SCALING_FACTOR_N <- lambda_base

message(paste("FACTOR LAMBDA BASE (Sin alpha):", round(lambda_base, 8)))
message(paste("FACTOR LAMBDA ESCALADO (USADO):", round(SCALING_FACTOR_N, 8), "(≈ 2.5e-4)"))

#-------------------------------------------------------------------------------
# --- 3. Nota sobre el Espacio de Búsqueda ---
# La búsqueda se realiza mediante bucles anidados, no con una matriz precalculada.
# El espacio es n en [1, n_clasic] y c en [0, n - 1].
#-------------------------------------------------------------------------------
message("3. El espacio de búsqueda (n <= n_clasic) se explora mediante bucles anidados.")


#-------------------------------------------------------------------------------
# --- 4. Función de Optimización de Pérdida Combinada (L) ---
#-------------------------------------------------------------------------------

optimizar_loss_combinada <- function(min_p, max_p, escenario) {
  
  # a. Definir la Densidad a Priori (f(p)) para este escenario
  dens_eval_p <- dunif(prop, min = min_p, max = max_p)
  
  # b. Calcular la CO para el PLAN CLÁSICO 
  CO_clasic <- phyper(q = c_clasic, m = m_hyp, n = n_hyp, k = n_clasic, lower.tail = TRUE)
  
  # c. Calcular el RIESGO DE ACEPTACIÓN EN ZI DEL PLAN CLÁSICO (Factor de penalización)
  # Este es el factor fijo que activa la penalización del costo 'n/(c+1)'
  zi_acceptance_risk_clasic <- sum(CO_clasic[ind_zi] * dens_eval_p[ind_zi]) * delta_p
  
  message(paste("\nEscenario:", escenario))
  message(paste("  - Riesgo de Aceptación ZI Plan Clásico:", round(zi_acceptance_risk_clasic, 6)))
  message("  - Buscando el plan que minimiza L = RAT + lambda * (n/(c+1))*P_ZI_Aceptacion_Clasica...")
  
  L_minimo_global <- Inf
  n_optimo <- NA
  c_optimo <- NA
  RAT_optimo_final <- NA
  
  # Búsqueda exhaustiva en una fase con bucles anidados (n y c)
  # n va desde 1 hasta n_clasic
  for (n_actual in 1:n_clasic) {
    # c va desde 0 hasta n_actual - 1
    for (c_actual in 0:(n_actual - 1)) {
      
      # 1. Calcular Riesgo Total Esperado (RAT)
      CO_actual <- phyper(q = c_actual, m = m_hyp, n = n_hyp, k = n_actual, lower.tail = TRUE)
      rp_actual <- sum((1 - CO_actual[ind_aql]) * dens_eval_p[ind_aql]) * delta_p
      rc_actual <- sum(CO_actual[ind_rql] * dens_eval_p[ind_rql]) * delta_p
      RAT_actual <- rp_actual + rc_actual
      
      # 2. Calcular la Función de Pérdida Combinada (L)
      penalizacion <- SCALING_FACTOR_N * (n_actual / (c_actual + 1)) * zi_acceptance_risk_clasic
      L_actual <- RAT_actual + penalizacion
      
      # 3. Condición de Minimización ABSOLUTA
      if (L_actual < L_minimo_global) {
        L_minimo_global <- L_actual
        n_optimo <- n_actual
        c_optimo <- c_actual
        RAT_optimo_final <- RAT_actual
      }
    }
  }
  
  # e. Riesgos del plan clásico (para referencia bajo esta distribución prior)
  rp_clasic <- sum((1 - CO_clasic[ind_aql]) * dens_eval_p[ind_aql]) * delta_p
  rc_clasic <- sum(CO_clasic[ind_rql] * dens_eval_p[ind_rql]) * delta_p
  
  # ---------------------------------------------------------------------------
  # NUEVA SECCIÓN: Cálculo de la Pérdida Total (L) del Plan Clásico (Coherencia)
  # ---------------------------------------------------------------------------
  penalizacion_clasic <- SCALING_FACTOR_N * (n_clasic / (c_clasic + 1)) * zi_acceptance_risk_clasic
  L_clasico <- (rp_clasic + rc_clasic) + penalizacion_clasic
  
  # El campo 'RAT_Clasico' ahora contendrá el L_Clasico
  return(data.frame(
    Escenario = escenario,
    n_Clasico = n_clasic,
    c_Clasico = c_clasic,
    RAT_Clasico = L_clasico, # Ahora representa L_Clasico (Riesgo + Costo)
    n_Optimo_Combinado = n_optimo,
    c_Optimo_Combinado = c_optimo,
    RAT_Optimo_Combinado = RAT_optimo_final,
    L_Optimo_Combinado = L_minimo_global
  ))
}

#-------------------------------------------------------------------------------
# --- 5. Ejecución de los 5 Escenarios (Distribución Prior Específica) ---
#-------------------------------------------------------------------------------
escenarios_params <- list(
  list(min_p = 0.00, max_p = 0.03, nombre = "1. Proveedor Excelente (p en AQL)"),
  list(min_p = 0.00, max_p = 0.08, nombre = "2. Proveedor Bueno (cerca de ZI)"),
  list(min_p = 0.03, max_p = 0.13, nombre = "3. Proveedor Regular (cubre ZI)"),
  list(min_p = 0.08, max_p = 0.13, nombre = "4. Proveedor Malo (cubre RQL)"),
  list(min_p = 0.13, max_p = 0.20, nombre = "5. Proveedor Muy Malo (p en RQL)")
)

resultados <- data.frame()

for (esc in escenarios_params) {
  res <- optimizar_loss_combinada(esc$min_p, esc$max_p, esc$nombre)
  resultados <- bind_rows(resultados, res)
}

#-------------------------------------------------------------------------------
# --- 6. Impresión de Resultados ---
#-------------------------------------------------------------------------------
message("\n=======================================================================================")
message("  RESULTADOS FINALES: OPTIMIZACIÓN BAYESIANA CON PÉRDIDA COMBINADA (LAMBDA SIN ALPHA)")
message("=======================================================================================")
message(paste("Fórmula de Lambda (Sin alpha):", "(AQL/RQL) * (beta / N) =", round(lambda_base, 8)))
message("Factor de Refuerzo (K) aplicado: N/A")
message(paste("Factor de Escala (lambda) USADO:", round(SCALING_FACTOR_N, 8)))
message("Fórmula de Pérdida Minimizada (L): RAT + lambda * (n/(c+1)) * P_ZI_Aceptacion_Clasica")
message("---------------------------------------------------------------------------------------")
message("NOTA: El factor de costo es 2.5 veces MÁS FUERTE que la calibración óptima. Se espera un fuerte sesgo a muestras pequeñas.")
message("---------------------------------------------------------------------------------------")

print(resultados)
