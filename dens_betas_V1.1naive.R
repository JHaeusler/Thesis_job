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

# Función principal para calcular el Riesgo Ponderado Integrado (RP y RC)
calc_wr <- function(n, c, alpha_b, beta_b, AQL, LTPD) {
  
  # Riesgo del Productor (RP) Integrado [0, AQL]
  rp_val <- integral(f = function(p) wr_p(p, n, c, alpha_b, beta_b, N), 
                     xmin = 0, xmax = AQL, 
                     method = "Kron") # Método de Gauss-Kronrod (preciso)
  
  # Riesgo del Consumidor (RC) Integrado [LTPD, 1]
  rc_val <- integral(f = function(p) wr_c(p, n, c, alpha_b, beta_b, N), 
                     xmin = LTPD, xmax = 1, 
                     method = "Kron")
  
  # TWR_val es el Riesgo Total Ponderado puro (RP + RC)
  return(c(RP_val = rp_val, RC_val = rc_val, TWR_val = rp_val + rc_val))
}

# Función para calcular la masa de probabilidad (densidad acumulada) en las regiones de riesgo
calc_prob_mass <- function(alpha_b, beta_b, AQL, LTPD) {
  # P_Good: Masa de probabilidad en la región del Riesgo del Productor [0, AQL]
  P_Good <- pbeta(AQL, shape1 = alpha_b, shape2 = beta_b)
  
  # P_Bad: Masa de probabilidad en la región del Riesgo del Consumidor [LTPD, 1]
  P_Bad <- 1 - pbeta(LTPD, shape1 = alpha_b, shape2 = beta_b)
  
  return(c(P_Good = P_Good, P_Bad = P_Bad))
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
TWR_U01_Clasic <- risks_uniforme["TWR_val"]

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
  TWR_U01_Clasic = TWR_U01_Clasic, # Mantenemos TWR para el cálculo intermedio
  # Riesgos bajo Información (Plan Clásico - Solo para referencia)
  RP_clasic = NA,
  RC_clasic = NA,
  TWR_clasic = NA, 
  # Plan Óptimo (Minimiza TWR bajo Información)
  n_opt = NA,
  c_opt = NA,
  RP_opt = NA, 
  RC_opt = NA,
  TWR_opt = NA, 
  # Ganancias (Plan Clásico Naive vs. Plan Óptimo Informado)
  RP_Ganancia = NA, 
  RC_Ganancia = NA,
  TWR_Ganancia = NA
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
  
  resultados_riesgo[i, "TWR_clasic"] <- risks_clasic["TWR_val"]
  resultados_riesgo[i, "RP_clasic"] <- risks_clasic["RP_val"]
  resultados_riesgo[i, "RC_clasic"] <- risks_clasic["RC_val"]
  
  # --- III. Búsqueda del Plan Óptimo (Minimizar TWR Puro) para el Escenario i ---
  
  n_opt_found <- NA
  c_opt_found <- NA
  min_twr <- Inf 
  
  # Búsqueda exhaustiva: Itera n de 1 hasta n_clasic (máximo tamaño de muestra del plan clásico)
  for (n_ in 1:n_clasic) { 
    for (c_ in 0:(n_ - 1)) {
      
      risks_opt <- calc_wr(n = n_, c = c_, 
                           alpha_b = alpha_b_val, beta_b = beta_b_val, 
                           AQL, LTPD)
      
      twr_current <- risks_opt["TWR_val"]
      
      # *** CRITERIO DE SELECCIÓN DEL PLAN ÓPTIMO (n_opt, c_opt) ***
      # Se utiliza el Riesgo Total Ponderado (TWR) para encontrar el único plan 
      # que ofrece el mejor balance entre el Riesgo del Productor (RP) y el 
      # Riesgo del Consumidor (RC). 
      # NOTA: Minimizar TWR es necesario para encontrar UN solo plan.
      # La comparación separada de RP y RC se hace en la tabla final (RP_Ganancia, RC_Ganancia).
      if (twr_current < min_twr) {
        min_twr <- twr_current
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
    
    resultados_riesgo[i, "TWR_opt"] <- risks_opt_final["TWR_val"]
    resultados_riesgo[i, "RP_opt"] <- risks_opt_final["RP_val"]
    resultados_riesgo[i, "RC_opt"] <- risks_opt_final["RC_val"]
    
    # Cálculo de las ganancias (Plan Clásico Naive vs. Óptimo con Información)
    # Ganancia RP = Riesgo RP Naive - Riesgo RP Óptimo
    # Ganancia RC = Riesgo RC Naive - Riesgo RC Óptimo
    resultados_riesgo[i, "RP_Ganancia"] <- resultados_riesgo[i, "RP_U01_Clasic"] - resultados_riesgo[i, "RP_opt"]
    resultados_riesgo[i, "RC_Ganancia"] <- resultados_riesgo[i, "RC_U01_Clasic"] - resultados_riesgo[i, "RC_opt"]
    resultados_riesgo[i, "TWR_Ganancia"] <- resultados_riesgo[i, "TWR_U01_Clasic"] - resultados_riesgo[i, "TWR_opt"]
    
  } else {
    resultados_riesgo[i, "TWR_opt"] <- NA
    resultados_riesgo[i, "RP_opt"] <- NA
    resultados_riesgo[i, "RC_opt"] <- NA
    resultados_riesgo[i, "RP_Ganancia"] <- NA
    resultados_riesgo[i, "RC_Ganancia"] <- NA
    resultados_riesgo[i, "TWR_Ganancia"] <- NA
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

# =========================================================================
# === SECCIÓN DE VISUALIZACIÓN  ===
# =========================================================================

proveedores_labels <- c("Excelente", "Bueno", "Regular", "Malo", "Muy Malo")
p_range <- seq(0, 0.3, length.out = 500) # Rango de p hasta 30%

# --- Definición del Layout para 2 Gráficos Grandes (1 CO, 1 Densidad) ---

# Matriz de layout: 2 filas y 1 columna
layout(matrix(c(1, 2), nrow = 2, ncol = 1), heights = c(1.5, 1)) # Más espacio para el gráfico CO
par(mar = c(4, 4, 2, 1), oma = c(0, 0, 0, 0)) # Márgenes para el gráfico superior

# -------------------------------------------------------------------------
# VISUALIZACIÓN 1 (PARTE SUPERIOR): Curvas CO Agrupadas (1 Gráfico)
# -------------------------------------------------------------------------

# Calcular la curva CO para el Plan Clásico (uniforme/naive) una sola vez
Pa_clasic_base <- Pa(n = n_clasic, c = c_clasic, p = p_range, N = N)

# Inicializar el plot con la curva Clásica
plot(p_range, Pa_clasic_base, type = "l", lty = 2, lwd = 2, col = "gray50",
     main = "Comparación de Curvas CO: Clásico (Naive) vs. Óptimo (Bayesiano) por Proveedor",
     xlab = "Proporción Defectuosa (p)",
     ylab = "Probabilidad de Aceptación (PA)",
     ylim = c(0, 1),
     cex.main = 1, cex.lab = 1)

# Añadir líneas de referencia AQL y LTPD
segments(AQL, 1 - alpha, AQL, -0.1, lty = 3, col = "black")
segments(LTPD, beta, LTPD, -0.1, lty = 3, col = "black")
text(AQL, 0.05, "AQL", cex = 0.8, pos = 4)
text(LTPD, 0.05, "LTPD", cex = 0.8, pos = 4)


legend_labels <- c()
legend_colors <- c()
legend_ltys <- c()
legend_lwds <- c()

# 1. Agregar el Plan Clásico a la leyenda
legend_labels <- c(legend_labels, paste0("Clásico Naive (", n_clasic, ", ", c_clasic, ")"))
legend_colors <- c(legend_colors, "gray50")
legend_ltys <- c(legend_ltys, 2)
legend_lwds <- c(legend_lwds, 2)

# 2. Definir paleta de colores para los 5 escenarios Óptimos
optimal_colors <- c("darkgreen", "blue", "orange", "red", "darkred")

# Loop para añadir las 5 curvas Óptimas
for (i in 1:5) {
  n_opt_val <- resultados_riesgo[i, "n_opt"]
  c_opt_val <- resultados_riesgo[i, "c_opt"]

  Pa_opt <- Pa(n = n_opt_val, c = c_opt_val, p = p_range, N = N)

  # Dibujar Curva CO Óptima
  lines(p_range, Pa_opt, lty = 1, lwd = 2, col = optimal_colors[i])

  # Agregar Plan Óptimo a la leyenda
  legend_labels <- c(legend_labels, paste0("Óptimo (", proveedores_labels[i], ": ", n_opt_val, ", ", c_opt_val, ")"))
  legend_colors <- c(legend_colors, optimal_colors[i])
  legend_ltys <- c(legend_ltys, 1)
  legend_lwds <- c(legend_lwds, 2)
}

# Mostrar Leyenda Completa
legend("topright",
       legend = legend_labels,
       col = legend_colors,
       lty = legend_ltys,
       lwd = legend_lwds,
       cex = 0.8,
       title = "Planes de Muestreo")


# -------------------------------------------------------------------------
# VISUALIZACIÓN 2 (PARTE INFERIOR): Densidades a Priori Agrupadas (1 Gráfico Grande)
# -------------------------------------------------------------------------

# Cambiar márgenes para el gráfico inferior
par(mar = c(5, 4, 3, 2) + 0.1)

# Calcular el límite Y máximo
max_y_limit <- max(dbeta(p_range, alpha_beta_params$alpha_b, alpha_beta_params$beta_b)) * 1.1

plot(p_range, dbeta(p_range, 1, 1), type = "l", col = "gray50", lty = 2, lwd = 2,
     main = "Distribuciones a Priori de la Calidad del Proveedor (p)",
     xlab = "Proporción Defectuosa (p)",
     ylab = "Densidad de Probabilidad f(p)",
     ylim = c(0, max_y_limit))

colores <- c("darkgreen", "blue", "orange", "red", "darkred")
for (i in 1:5) {
  lines(p_range, dbeta(p_range, alpha_beta_params[i, "alpha_b"], alpha_beta_params[i, "beta_b"]),
        col = colores[i], lwd = 2)
}

# Líneas de referencia para AQL y LTPD
abline(v = AQL, col = "black", lty = 3)
text(AQL, max_y_limit * 0.95,
     "AQL", cex = 0.8, pos = 4)
abline(v = LTPD, col = "black", lty = 3)
text(LTPD, max_y_limit * 0.95,
     "LTPD", cex = 0.8, pos = 4)

# legend("topright",
#        legend = c("Naive (Uniforme)", proveedores_labels),
#        col = c("gray50", colores),
#        lty = c(2, rep(1, 5)),
#        lwd = 2,
#        cex = 0.9,
#        title = "Prior (Beta)")

# Resetear el layout (Importante para la estabilidad del entorno)
layout(matrix(1), widths = 1, heights = 1)
par(mfrow = c(1, 1))

# Imprimir parámetros Beta finales para referencia
print("Parámetros (alpha_b, beta_b) de la Distribución Beta para cada Proveedor:")
print(alpha_beta_params)