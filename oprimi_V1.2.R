# ==============================================================================

# SCRIPT FINAL: Optimizacion basada en L = RAT + Penalización_ZI_Severidad + Penalización_Costo_Operacional

# - Penalización Costo: lambda*(n/(c+1))*(n/N). Fuerza a encontrar el menor 'n'.

# - Metodología: Búsqueda Secuencial (Min. n)

# - LÓGICA DE PARADA: Primera mejora respecto a L_clasico.

# ==============================================================================


# Carga de librerías necesarias

library(AcceptanceSampling)

library(tidyverse)

library(dplyr)

library(tidyr)

library(gridExtra)


# Establecer semilla para reproducibilidad

set.seed(123)


message("-----------------------------------------------------")

message("--- 1. Definición de Parámetros y Constantes ---")

message("-----------------------------------------------------")


N <- 200       # Tamaño del lote

alpha <- 0.05 # Riesgo del productor

beta <- 0.10  # Riesgo del consumidor


AQL <- 0.05   # Nivel de calidad aceptable

RQL <- 0.10   # Nivel de calidad de rechazo

digs <- 5

delta_p <- 10^-digs # Tolerancia para el paso de p (1e-5)


# Espacio de calidades (p)

prop <- seq(0, 1, by = delta_p)


# Índices para las zonas de evaluación

ind_aql <- which(prop <= AQL)          # Zona AQL (Aceptación)

ind_rql <- which(prop >= RQL)          # Zona RQL (Rechazo)

ind_zi <- which(prop > AQL & prop < RQL) # Zona de Indiferencia (ZI)


# Parámetros de la distribución Hipergeométrica para cada p

m_hyp = N * prop

n_hyp = N * (1 - prop)


#-------------------------------------------------------------------------------

# --- 2. Cálculo del plan de muestreo clásico (Referencia) ---

#---------------------------------------------------------------

message("2. Calculando el plan de muestreo clásico...")


# Plan Clásico (Hipergeométrico)

plan_clasic <- find.plan(PRP = c(AQL, 1 - alpha),
                         
                         CRP = c(RQL, beta),
                         
                         N = N, type = "hypergeom")


n_clasic <- plan_clasic$n 

c_clasic <- plan_clasic$c

message(paste("Plan Clásico de Referencia: n =", n_clasic, ", c =", c_clasic))


# FACTOR DE ESCALA CRÍTICO (lambda)

SCALING_FACTOR_N <- (AQL / RQL) * (alpha * beta / N)


message(paste("FACTOR LAMBDA ESCALADO (USADO):", round(SCALING_FACTOR_N, 8)))


#-------------------------------------------------------------------------------

# --- 3. Función de Optimización de Pérdida Combinada (L) [VERSIÓN CON COSTO] ---

#-------------------------------------------------------------------------------


optimizar_loss_combinada <- function(min_p, max_p, escenario) {
  
  
  
  # a. Definir la Densidad a Priori (f(p)) para este escenario
  
  dens_eval_p <- dunif(prop, min = min_p, max = max_p)
  
  
  
  # b. Calcular la CO para el PLAN CLÁSICO (usado para la penalización de anclaje) 
  
  CO_clasic <- phyper(q = c_clasic, m = m_hyp, n = n_hyp, k = n_clasic, lower.tail = TRUE)
  
  
  
  # c. Calcular el RIESGO DE ACEPTACIÓN EN ZI DEL PLAN CLÁSICO (para Penalización 1)
  
  zi_acceptance_risk_clasic <- sum(CO_clasic[ind_zi] * dens_eval_p[ind_zi] * delta_p)
  
  
  
  message(paste("\nEscenario:", escenario))
  
  message(paste("  - Riesgo de Aceptación ZI Plan Clásico:", round(zi_acceptance_risk_clasic, 6)))
  
  message("  - Buscando el plan que minimiza L = RAT + Penalizacion_ZI_Severidad + Penalizacion_Costo...")
  
  
  
  # e. Riesgos del plan clásico (para referencia bajo esta distribución prior)
  
  rp_clasic <- sum((1 - CO_clasic[ind_aql]) * dens_eval_p[ind_aql] * delta_p)
  
  rc_clasic <- sum(CO_clasic[ind_rql] * dens_eval_p[ind_rql] * delta_p)
  
  
  
  # ---------------------------------------------------------------------------
  
  # Cálculo de la Pérdida Total (L) del Plan Clásico (Para comparación)
  
  # L = RAT + Penalización 1 + Penalización 2 (Costo)
  
  # ---------------------------------------------------------------------------
  
  # Penalización 1 (ZI-Severidad)
  
  penalizacion_clasic_zi_severidad <- SCALING_FACTOR_N * (n_clasic / (c_clasic + 1)) * zi_acceptance_risk_clasic
  
  
  
  # Penalización 2 (Costo Operacional)
  
  # penalizacion_clasic_costo <- SCALING_FACTOR_N * (n_clasic / (c_clasic + 1)) * sum((CO_clasic*delta_p))

  
  L_clasico <- (rp_clasic + rc_clasic) + penalizacion_clasic_zi_severidad #+ penalizacion_clasic_costo
  
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
      
      rp_actual <- sum((1 - CO_actual[ind_aql]) * dens_eval_p[ind_aql] * delta_p)
      
      rc_actual <- sum(CO_actual[ind_rql] * dens_eval_p[ind_rql] * delta_p)
      
      RAT_actual <- rp_actual + rc_actual
      
      
      
      # Riesgo ZI actual (incluye f(p) y delta_p) para Penalización 1
      
      zi_acceptance_risk_actual <-  sum(CO_actual[ind_zi] * dens_eval_p[ind_zi] * delta_p)
      
      
      
      # 2. Calcular la Función de Pérdida Combinada (L)
      
      
      
      # Penalización 1: Severidad del Plan Actual * Riesgo ZI Actual * lambda
      
      penalizacion_actual_zi_severidad <- SCALING_FACTOR_N * (n_actual / (c_actual + 1)) * zi_acceptance_risk_actual
      
      
      
      # Penalización 2: Costo Operacional (Nuevo término)
      
      # penalizacion_actual_costo <- SCALING_FACTOR_N * (n_actual / (c_actual + 1)) * sum((CO_actual*delta_p))
      
      
      
      L_actual <- RAT_actual + penalizacion_actual_zi_severidad + SCALING_FACTOR_N * (n_actual / (c_actual + 1))*(sum((CO_clasic*delta_p)) - sum((CO_actual*delta_p))) #+ penalizacion_actual_costo
      
      
      
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
  
  
  
  # Nota: Para mantener la coherencia con el output anterior, Penalizacion_Clasica
  
  # solo incluye el término de Riesgo ZI (Penalización 1).
  
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

#---------------------------------------------------------------


message("\n-----------------------------------------------------")

message("--- 4. Ejecutando la Optimización para 5 Escenarios ---")

message("-----------------------------------------------------")


# Definición de los 5 escenarios (Distribuciones uniformes)

escenarios <- list(
  
  # E1: Calidad Excelente (Riesgo bajo)
  
  list(min_p = 0.00, max_p = 0.03, nombre = "E1 - Excelente"), 
  
  # E2: Calidad Buena (Riesgo medio-bajo)
  
  list(min_p = 0.03, max_p = 0.08, nombre = "E2 - Bueno"),       
  
  # E3: Calidad Regular (Cubre ZI, Riesgo medio-alto)
  
  list(min_p = 0.03, max_p = 0.13, nombre = "E3 - Regular"),   
  
  # E4: Calidad Mala (Riesgo alto)
  
  list(min_p = 0.08, max_p = 0.13, nombre = "E4 - Malo"),       
  
  # E5: Calidad Muy Mala (Riesgo muy alto, solo rechazo)
  
  list(min_p = 0.13, max_p = 0.20, nombre = "E5 - Muy Malo")    
  
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


# Preparación para la visualización (Se asume la Curva CO Tipo B para suavizado)

get_CO_curve_smooth <- function(n, c) {
  
  # Usamos pbinom para una curva suave (Tipo B)
  
  if (is.na(n)) return(rep(NA, length(prop)))
  
  return(pbinom(q = c, size = n, prob = prop, lower.tail = TRUE))
  
}


# Tabla de planes para la visualización (usa los resultados óptimos)

planes_visualizacion <- data.frame(
  
  Escenario = resultados_optimizacion$Escenario,
  
  n = resultados_optimizacion$n_Optimo_Combinado,
  
  c = resultados_optimizacion$c_Optimo_Combinado,
  
  Tipo = "Óptimo Secuencial (Min n)"
  
)

planes_visualizacion <- rbind(
  
  data.frame(Escenario = "Clásico (Referencia)", n = n_clasic, c = c_clasic, Tipo = "Clásico"), 
  
  planes_visualizacion
  
) %>%
  
  filter(!is.na(n)) # Eliminar filas donde n es NA (si ocurre un error de optimización)


# Generación de curvas CO y densidades

curvas_co_df <- planes_visualizacion %>%
  
  rowwise() %>%
  
  mutate(CO = list(get_CO_curve_smooth(n, c))) %>%
  
  unnest(CO) %>%
  
  mutate(Proporcion = rep(prop, times = nrow(planes_visualizacion)))


# 5. Densidades a priori para el gráfico inferior

densidades_df <- data.frame(
  
  Proporcion = prop,
  
  E1 = dunif(prop, min = 0.00, max = 0.03),
  
  E2 = dunif(prop, min = 0.03, max = 0.08),
  
  E3 = dunif(prop, min = 0.03, max = 0.13),
  
  E4 = dunif(prop, min = 0.08, max = 0.13),
  
  E5 = dunif(prop, min = 0.13, max = 0.20)
  
)

densidades_long <- densidades_df %>%
  
  pivot_longer(cols = starts_with("E"), names_to = "Escenario", values_to = "Densidad") %>%
  
  mutate(Escenario = factor(Escenario, levels = c("E1", "E2", "E3", "E4", "E5")))



message("\n-----------------------------------------------------")

message("--- 5. Generación de Gráficos ---")

message("-----------------------------------------------------")


# Gráfico 1: Curvas CO (Parte superior)

plot_co <- ggplot(curvas_co_df, aes(x = Proporcion, y = CO, color = Escenario, linetype = Tipo)) +
  
  geom_line(size = 1.2) +
  
  geom_vline(xintercept = AQL, linetype = "dashed", color = "gray50") +
  
  geom_vline(xintercept = RQL, linetype = "dashed", color = "gray50") +
  
  annotate("text", x = AQL, y = 1.05, label = "AQL", size = 3, color = "gray30") +
  
  annotate("text", x = RQL, y = 1.05, label = "RQL", size = 3, color = "gray30") +
  
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
  
  labs(title = "Curvas CO para Planes Secuenciales Óptimos vs. Clásico",
       
       y = "Probabilidad de Aceptación (Pa)",
       
       x = "") +
  
  theme_minimal() +
  
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5))


# Gráfico 2: Densidades a Priori (Parte inferior)

plot_dens <- ggplot(densidades_long, aes(x = Proporcion, y = Densidad, fill = Escenario)) +
  
  geom_area(position = "identity", alpha = 0.6) +
  
  geom_vline(xintercept = AQL, linetype = "dashed", color = "gray50") +
  
  geom_vline(xintercept = RQL, linetype = "dashed", color = "gray50") +
  
  labs(y = "Densidad a Priori f(p)",
       
       x = "Proporción de Defectuosos (p)") +
  
  theme_minimal() +
  
  theme(legend.position = "none") +
  
  scale_x_continuous(limits = c(0, 0.2), breaks = seq(0, 0.2, 0.02))


# Combinar gráficos

grid.arrange(plot_co, plot_dens, ncol = 1, heights = c(3, 1.5)) 