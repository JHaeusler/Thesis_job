# Carga de la librería necesaria
library(AcceptanceSampling)

# ----------------------------------------------------------------------
# 1. Parámetros Globales (Fijos para todos los escenarios)
# ----------------------------------------------------------------------
N <- 200 # Tamaño del lote (finito)
AQL <- 0.05 # Nivel de Calidad Aceptable
RQL <- 0.10 # Nivel de Calidad Límite de Tolerancia

alpha <- 0.05 # Riesgo del Productor
beta <- 0.10 # Riesgo del Consumidor

# Parámetros de Discretización (para aproximación de la integral)
digs <- 5
delta_p <- 10^-digs # Incremento
prop <- seq(0, 1, delta_p) # Vector de posibles proporciones defectuosas (p)

# Vectores para la Curva CO (Hipergeométrica, asumimos lote finito)
m_hyp <- N * prop # Número de defectuosos en el lote si p es 'prop'
n_hyp <- N - m_hyp # Número de no defectuosos

# Índices para la Integración (RAT)
ind_aql <- which(prop <= AQL) # Región Aceptable [0, AQL]
ind_rql <- which(prop >= RQL) # Región Rechazable [RQL, 1]

# Densidad Naive (Uniforme no informativa)
dens_naive <- dunif(prop, min = 0, max = 1) 

# ----------------------------------------------------------------------
# 2. Función para Ejecutar el Análisis por Escenario
# ----------------------------------------------------------------------

#' Ejecuta el análisis de la búsqueda secuencial para un escenario de calidad dado.
#' El plan óptimo busca mejorar el RAT Clásico evaluado con la densidad NAIVE (Uniforme[0,1]).
#' La evaluación del plan óptimo se realiza con la densidad del escenario (dens_eval).
#' @param nombre_escenario Nombre del escenario (e.g., "Excelente").
#' @param p_min Límite inferior de la distribución Uniforme para f(p).
#' @param p_max Límite superior de la distribución Uniforme para f(p).
ejecutar_analisis_escenario <- function(nombre_escenario, p_min, p_max) {
  
  cat("\n========================================================\n")
  cat(sprintf("ESCENARIO: %s (f(p) ~ Uniforme[%.3f, %.3f])\n", nombre_escenario, p_min, p_max))
  cat("========================================================\n")
  
  # 2.1. Definición de la Densidad del Escenario (f(p)) - Densidad de Evaluación
  dens_eval <- dunif(prop, min = p_min, max = p_max)
  
  # 2.2. Cálculo del Plan Clásico (Referencia)
  plan_c <- find.plan(PRP = c(AQL, 1 - alpha),
                      CRP = c(RQL, beta),
                      N = N, type = "hypergeom")
  
  n_clasico <- plan_c$n
  c_clasico <- plan_c$c
  
  # Curva CO del Plan Clásico (P_a(p))
  CO_clasic <- phyper(q = c_clasico, m = m_hyp, k = n_clasico, n = n_hyp, lower.tail = TRUE)
  
  # 2.3. Cálculo de los RAT Clásicos (Baselines)
  
  # a) RAT Naive (RAT Clásico evaluado con la densidad Naive Uniforme[0,1])
  rp_naive <- sum((1 - CO_clasic[ind_aql]) * dens_naive[ind_aql] * delta_p)
  rc_naive <- sum(CO_clasic[ind_rql] * dens_naive[ind_rql] * delta_p)
  wr_c_naive <- rp_naive + rc_naive
  
  # b) RAT Escenario (RAT Clásico evaluado con la densidad Real del Escenario)
  rp_scenario <- sum((1 - CO_clasic[ind_aql]) * dens_eval[ind_aql] * delta_p)
  rc_scenario <- sum(CO_clasic[ind_rql] * dens_eval[ind_rql] * delta_p)
  wr_c_scenario <- rp_scenario + rc_scenario
  
  cat(sprintf("Plan Clásico: n=%d, c=%d\n", n_clasico, c_clasico))
  cat(sprintf("   - RAT Naive (Baseline de Búsqueda): %f\n", wr_c_naive))
  cat(sprintf("   - RAT Escenario (Riesgo Real): %f\n", wr_c_scenario))
  
  # 2.4. Búsqueda Secuencial Discreta
  
  # Definición del espacio de búsqueda (n <= n_clasico)
  n_busq <- 1:n_clasico
  c_busq <- 0:n_clasico
  space_busq <- expand.grid(n = n_busq, c = c_busq)
  
  # Filtramos la grilla para asegurar la validez: n > c (como solicitaste) y n > 0
  space_busq_filt <- space_busq[space_busq$n > space_busq$c & space_busq$n > 0, ]
  
  # Ordenamos (prioridad a la eficiencia: menor n, luego menor c)
  space_busq_filt <- space_busq_filt[order(space_busq_filt$n, space_busq_filt$c), ]
  
  # Inicialización de la búsqueda
  mejor_plan <- NULL
  min_risk <- wr_c_naive # <--- BASELINE DE BÚSQUEDA AHORA ES wr_c_naive
  found_better <- FALSE
  
  for (i_ in 1:nrow(space_busq_filt)) {
    n_b <- space_busq_filt[i_, "n"]
    c_b <- space_busq_filt[i_, "c"]
    
    # 1. Calcular la Curva CO para el plan (n_b, c_b)
    CO_b <- phyper(q = c_b, m = m_hyp, k = n_b, n = n_hyp, lower.tail = TRUE)
    
    # 2. Calcular el RAT (WR_B) utilizando la densidad evaluada (dens_eval)
    rp_b <- sum((1 - CO_b[ind_aql]) * dens_eval[ind_aql] * delta_p)
    rc_b <- sum(CO_b[ind_rql] * dens_eval[ind_rql] * delta_p)
    wr_b <- rp_b + rc_b
    
    # 3. Criterio de Parada: ¿El nuevo riesgo (evaluado con dens_eval) es mejor que el RAT Naive?
    if (wr_b < min_risk) {
      min_risk <- wr_b
      mejor_plan <- c(n = n_b, c = c_b)
      found_better <- TRUE # Detener al encontrar el primero mejorado
      break # Salir del loop for, ya que Kirmeier busca el "primer"
    }
  }
  
  # 2.5. Presentación de Resultados
  
  if (found_better) {
    cat("\n--- Plan Óptimo Encontrado (Búsqueda Secuencial) ---\n")
    cat(sprintf("Plan Secuencial: n=%d, c=%d\n", mejor_plan["n"], mejor_plan["c"]))
    cat(sprintf("RAT Secuencial (Evaluado con f(p)): %f\n", min_risk))
    cat(sprintf("Mejora de Riesgo vs Naive Baseline: %f\n", wr_c_naive - min_risk))
    
    if (mejor_plan["n"] < n_clasico) {
      cat(sprintf("¡Mejora en Eficiencia! El tamaño de muestra se redujo de %d a %d.\n", n_clasico, mejor_plan["n"]))
    } else if (mejor_plan["n"] == n_clasico) {
      cat("Mejora en riesgo sin reducción de tamaño de muestra.\n")
    }
  } else {
    cat("\nNo se encontró un plan que mejore el RAT Naive clásico en el espacio de búsqueda restringido.\n")
  }
}

# ----------------------------------------------------------------------
# 3. Ejecución de los 5 Casos de Proveedores (Densidades Uniformes)
# ----------------------------------------------------------------------

# i) PROVEEDOR EXCELENTE: Densidad concentrada antes de AQL (0.05)
ejecutar_analisis_escenario("I. EXCELENTE", p_min = 0.000, p_max = 0.020)

# ii) PROVEEDOR BUENO: Densidad concentrada alrededor de AQL (0.05)
ejecutar_analisis_escenario("II. BUENO", p_min = 0.030, p_max = 0.070)

# iii) PROVEEDOR REGULAR: La densidad contiene AQL (0.05) y RQL (0.10)
ejecutar_analisis_escenario("III. REGULAR", p_min = 0.030, p_max = 0.120)

# iv) PROVEEDOR MALO: Densidad concentrada alrededor de RQL (0.10)
ejecutar_analisis_escenario("IV. MALO", p_min = 0.080, p_max = 0.150)

# v) PROVEEDOR MUY MALO: Densidad por encima de RQL (0.10)
ejecutar_analisis_escenario("V. MUY MALO", p_min = 0.120, p_max = 0.200)
