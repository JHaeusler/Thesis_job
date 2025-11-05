# ==============================================================================
# SCRIPT FINAL: Optimizacion basada en L = RAT + Penalización_ZI_Severidad + Penalización_Costo_Intrínseco
# - Penalización Costo Intrínseco: (n+c)/N * Area_ZI_Acc. Se basa en el área de Pa en ZI (NO ponderada por f(p)).
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

N <- 200 # Tamaño del lote
alpha <- 0.05 # Riesgo del productor
beta <- 0.10 # Riesgo del consumidor

AQL <- 0.05 # Nivel de calidad aceptable
RQL <- 0.10 # Nivel de calidad de rechazo
digs <- 5
delta_p <- 10^-digs # Tolerancia para el paso de p (1e-5)

# Espacio de calidades (p)
prop <- seq(0, 1, by = delta_p)

# Índices para las zonas de evaluación
ind_aql <- which(prop <= AQL) # Zona AQL (Riesgo del Productor)
ind_rql <- which(prop >= RQL) # Zona RQL (Riesgo del Consumidor)
ind_zi <- which(prop > AQL & prop < RQL) # Zona de Indiferencia (Riesgo ZI)

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

#-------------------------------------------------------------------------------
# --- 3. Búsqueda del Plan Óptimo por Dominación de Riesgo ---
#-------------------------------------------------------------------------------
message("\n3. Ejecutando la búsqueda de plan óptimo por dominación de riesgo...")

# Definición del Escenario (p ~ U(0.00, 0.03))
min_p <- 0.0 
max_p <- 0.03

# a. Definir la Densidad a Priori (f(p))
dens_eval_p <- dunif(prop, min = min_p, max = max_p)

# b. Calcular la CO para el PLAN CLÁSICO
CO_clasic <- phyper(q = c_clasic, m = m_hyp, n = n_hyp, k = n_clasic, lower.tail = TRUE)

# c. Riesgos del plan clásico (para referencia bajo esta distribución prior)
rp_clasic <- sum((1 - CO_clasic[ind_aql]) * dens_eval_p[ind_aql]) * delta_p
rc_clasic <- sum(CO_clasic[ind_rql] * dens_eval_p[ind_rql]) * delta_p
rzi_clasic <- sum((1 - CO_clasic[ind_zi]) * dens_eval_p[ind_zi]) * delta_p

message(paste("   - Riesgo Productor Clásico (RP):", round(rp_clasic, 8)))
message(paste("   - Riesgo Consumidor Clásico (RC):", round(rc_clasic, 8)))
message(paste("   - Riesgo ZI Clásico (RZI):", round(rzi_clasic, 8)))
message("   - Buscando el plan más pequeño (min n) que cumpla RP_act <= RP_clasic y RC_act <= RC_clasic y RZI_act <= RZI_clasic.")

n_optimo <- NA
c_optimo <- NA
found_better_plan <- FALSE

# Búsqueda: Se detiene en la PRIMERA combinación (n, c) que domina.
# Iteramos desde n=1 hasta n_clasic
for (n_actual in 1:n_clasic) {
  # Iteramos c desde 0 hasta n_actual - 1
  for (c_actual in 0:(n_actual - 1)) { 
    
    # 1. Calcular Curva CO actual
    CO_actual <- phyper(q = c_actual, m = m_hyp, n = n_hyp, k = n_actual, lower.tail = TRUE)
    
    # 2. Calcular Riesgos Esperados Actuales
    rp_actual <- sum((1 - CO_actual[ind_aql]) * dens_eval_p[ind_aql]) * delta_p
    rc_actual <- sum(CO_actual[ind_rql] * dens_eval_p[ind_rql]) * delta_p
    rzi_actual <- sum((1 - CO_actual[ind_zi]) * dens_eval_p[ind_zi]) * delta_p
    
    # 3. CONDICIÓN DE DOMINACIÓN
    if (rp_actual <= rp_clasic && rc_actual <= rc_clasic && rzi_actual <= rzi_clasic) {
      n_optimo <- n_actual
      c_optimo <- c_actual
      found_better_plan <- TRUE
      break # Sale del bucle de c_actual (Búsqueda interna)
    }
  }
  
  if (found_better_plan) {
    break # Sale del bucle de n_actual (Se encontró el n más pequeño)
  }
}

# 4. Resultados de la optimización
results <- data.frame(
  n_Clasico = n_clasic,
  c_Clasico = c_clasic,
  n_Optimo = n_optimo,
  c_Optimo = c_optimo
)

message(paste("\n   -> Plan Óptimo (Dominación de Riesgo): n =", n_optimo, ", c =", c_optimo))
message(paste("   -> Reducción en el tamaño de muestra:", round(((n_clasic - n_optimo) / n_clasic) * 100, 2), "%"))

#-------------------------------------------------------------------------------
# --- 4. Generación de Gráfica de Curvas CO Comparativas ---
#-------------------------------------------------------------------------------
message("\n4. Generando la gráfica de Curvas CO comparativas...")

# Calcular la CO del Plan Óptimo
CO_optimo <- phyper(q = c_optimo, m = m_hyp, n = n_hyp, k = n_optimo, lower.tail = TRUE)

# Preparar datos para ggplot
df_co <- data.frame(p = prop,
                    CO_Clasico = CO_clasic,
                    CO_Optimo = CO_optimo) %>%
  # Transformar a formato largo para ggplot
  pivot_longer(cols = starts_with("CO_"), names_to = "Plan", values_to = "Pa")

# Generar la gráfica
plot_co <- df_co %>%
  ggplot(aes(x = p, y = Pa, color = Plan)) +
  geom_line(linewidth = 1.2) +
  
  # Líneas de referencia AQL y RQL
  geom_vline(xintercept = AQL, linetype = "dashed", color = "darkgreen") +
  geom_vline(xintercept = RQL, linetype = "dashed", color = "red") +
  
  # Puntos de referencia (Riesgos Clásicos)
  geom_point(aes(x = AQL, y = 1 - alpha), color = "darkgreen", size = 3) +
  geom_point(aes(x = RQL, y = beta), color = "red", size = 3) +
  
  # Etiquetas de las referencias
  annotate("text", x = AQL, y = 0.0, label = paste("AQL=", AQL), color = "darkgreen", hjust = -0.1, vjust = 0, angle = 90, size=3) +
  annotate("text", x = RQL, y = 0.0, label = paste("RQL=", RQL), color = "red", hjust = -0.1, vjust = 0, angle = 90, size=3) +
  
  scale_color_manual(values = c("CO_Clasico" = "#0072B2", "CO_Optimo" = "#D55E00"),
                     labels = c(paste0("Clásico (n=", n_clasic, ", c=", c_clasic, ")"), 
                                paste0("Óptimo (n=", n_optimo, ", c=", c_optimo, ")"))) +
  
  labs(title = "Curvas de Funcionamiento (CO) del Plan Clásico vs. Plan Óptimo (Dominación de Riesgo)",
       subtitle = paste0("Distribución a Priori: Uniforme(", min_p, ", ", max_p, ")"),
       x = "Proporción de Defectuosos (p)",
       y = "Probabilidad de Aceptación (Pa)",
       color = "Plan de Muestreo") +
  
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

# Imprimir la gráfica
print(plot_co)