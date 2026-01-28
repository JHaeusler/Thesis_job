<<<<<<< HEAD
# Script para comparar planes de muestreo de aceptación clásicos
# con planes óptimos basados en el Riesgo Ponderado.

=======
# Propuesta para diseño de planes de muestreo simple para atributos
# considerando el histórico de la proporción estimada de unidades
# defectuosas en la inspección de lotes

# Calidad del proveedor según el histórico
>>>>>>> 16a177554b09c394ca6f03b22b399b8fb6966908
# Las prioris de calidad de los proveedores se definen por pares (PA en AQL, PA en LTPD)
# Excelente proveedor: PA(AQL)=0.95, PA(LTPD)=1-0.0001
# Buen proveedor: PA(AQL)=0.75, PA(LTPD)=1-0.05
# Regular: PA(AQL)=0.50, PA(LTPD)=1-0.10
# Malo: PA(AQL)=0.25, PA(LTPD)=1-0.25
# Muy malo: PA(AQL)=0.10, PA(LTPD)=1-0.40

# --- Carga de Librerías ---
# install.packages("pracma") # Si no lo tienes instalado
<<<<<<< HEAD
library(pracma)
# install.packages("sjstats")
library(sjstats)
library(AcceptanceSampling)
library(GA)
# --- Variables Globales y Parámetros del Problema ---
N <- c(200, 1000)        # Tamaño del lote
alpha_des <- c(0.01, 0.03, 0.05)   # Riesgo del productor (para el plan clásico)
beta_des <- c(0.07, 0.10, 0.15, 0.20)    # Riesgo del consumidor (para el plan clásico)
AQL <- c(0.01, 0.03, 0.05)      # Nivel de Calidad Aceptable
LTPD <- c(0.05, 0.10, 0.15)   # Tolerancia de Porcentaje Defectuoso de Lote

escenarios <- expand.grid(N = N, AQL = AQL, LTPD = LTPD,
                          alpha_des = alpha_des, beta_des = beta_des)


# Vectores de probabilidades para los 5 escenarios
# p1 = Probabilidad de Aceptación en el AQL
# p2 = Probabilidad de No Aceptación (o Rechazo) en el LTPD
p1 <- c(0.95, 0.75, 0.50, 0.25, 0.10)
p2 <- c(1e-4, 0.05, 0.10, 0.25, 0.40) # 1 - PA(LTPD) -> Beta (Riesgo del Consumidor)
=======
# install.packages("sjstats")

library(pracma)
library(sjstats)
library(AcceptanceSampling)
library(GA)

# --- Funciones de Riesgo Ponderado y Masa de Probabilidad ---
>>>>>>> 16a177554b09c394ca6f03b22b399b8fb6966908

# Función para calcular la masa de probabilidad (densidad acumulada) en las regiones de riesgo
calc_prob_mass <- function(alpha_b, beta_b, AQL, LTPD) {
  P_Good <- pbeta(AQL, shape1 = alpha_b, shape2 = beta_b)
  P_Bad <- 1 - pbeta(LTPD, shape1 = alpha_b, shape2 = beta_b)
  return(c(P_Good = P_Good, P_Bad = P_Bad))
}

<<<<<<< HEAD
# Función principal para calcular el Riesgo Ponderado Integrado (RP y RC)
calc_wr <- function(n, c, alpha_b, beta_b, AQL, LTPD, k_p, k_c) {
  # 1. Calcular RAT (RAT_val)
  #rp_val <- integral(f = function(p) wr_p(p, n, c, alpha_b, beta_b, N), xmin = 0, xmax = AQL, method = "Kron")
  wrp_val <- k_p*(1-phyper(c, N * AQL, N * (1 - AQL), n))
  wrc_val <- k_c*phyper(c, N * LTPD, N * (1 - LTPD), n)
  WR_val <- wrp_val + wrc_val

  return(c(WRP_val = wrp_val, WRC_val = wrc_val, WR_val = WR_val))
}

# 1. Determinar el Plan Clásico (basado en alpha y beta fijos)
plan_clasic <- find.plan(PRP = c(AQL, 1 - alpha),
                         CRP = c(LTPD, beta),
                         N = N, type = "hypergeom")

n_clasic <- plan_clasic$n
c_clasic <- plan_clasic$c

=======
# Función principal para calcular el Riesgo Ponderado
calc_wr <- function(N_, n, c, alpha_b, beta_b, AQL_, LTPD_, k_p, k_c) {

  wrp_val <- k_p*(1-phyper(c, N_ * AQL_, N_ * (1 - AQL_), n))
  wrc_val <- k_c*phyper(c, N_ * LTPD_, N_ * (1 - LTPD_), n)
  wrt_val <- wrp_val + wrc_val

  return(c(WRP_val = wrp_val, WRC_val = wrc_val, WRT_val = wrt_val))
}

>>>>>>> 16a177554b09c394ca6f03b22b399b8fb6966908
decode <- function(string){
  string <- gray2binary(string)
  n <- binary2decimal(string[1:l1])
  c <- min(n, binary2decimal(string[(l1 + 1):(l1 + l2)]))
<<<<<<< HEAD

=======
  
>>>>>>> 16a177554b09c394ca6f03b22b399b8fb6966908
  return(c(n,c))
}

fitness <- function(string){
  par <- decode(string)
  n <- par[1]
  c <- par[2]
<<<<<<< HEAD
  Pa_p <- phyper(c, N*AQL, N*(1 - AQL), n)
  Pa_c <- phyper(c, N*LTPD, N*(1 - LTPD), n)
  Loss <- (Pa_p - (1 - alpha))^2 + (Pa_c - beta)^2
  -Loss
}

n_ran <- 1:N
c_ran <- 0:(max(n_ran) - 1)

b1 <- decimal2binary(max(n_ran)); l1 <- length(b1)
b2 <- decimal2binary(max(c_ran)); l2 <- length(b2)

plan_genetico <- ga(type = "binary", nBits = l1 + l2,
                    fitness = fitness, popSize = 200,
                    maxiter = 200, run = 100, seed = 060722)

plan_ga <- decode(plan_genetico@solution)

n_ga <- plan_ga[1]
c_ga <- plan_ga[2]

# Y la masa de probabilidad (densidad acumulada) para el prior uniforme
prob_mass_uniforme <- calc_prob_mass(alpha_b = 1, beta_b = 1, AQL, LTPD)
P_Mass_Good_Naive <- prob_mass_uniforme["P_Good"]
P_Mass_Bad_Naive <- prob_mass_uniforme["P_Bad"]

# -------------------------------------------------------------------

# 2. Obtener los 5 pares de (alpha_b, beta_b) para la distribución Beta
alpha_beta_params <- data.frame(alpha_b = rep(NA, length(p1)), beta_b = rep(NA, length(p1)))

for (i in 1:length(p1)) { # i <- 1 + i

  # find_beta encuentra los parámetros de la distribución Beta
  Shape <- find_beta(x1 = AQL, p1 = p1[i], x2 = LTPD, p2 = 1 - p2[i]) 
  alpha_beta_params[i, "alpha_b"] <- Shape$shape1
  alpha_beta_params[i, "beta_b"] <- Shape$shape2
}

alpha_beta_params <- rbind(alpha_beta_params, c(1, 1))

# 3. Inicializar la tabla de resultados
resultados_riesgo <- data.frame(
  Escenario = 1:dim(alpha_beta_params)[1],
  Proveedor = c("Excelente", "Bueno", "Regular", "Malo", "Muy Malo", "naive"),
  n_clasic = n_clasic,
  c_clasic = c_clasic,

  n_ga = n_ga,
  c_ga = c_ga,
  
  # Densidad Acumulada (Masa de Probabilidad) del Prior
  P_Mass_Good_Naive = P_Mass_Good_Naive, # P(p < AQL | Naive)
  P_Mass_Bad_Naive = P_Mass_Bad_Naive,   # P(p > LTPD | Naive)
  P_Mass_Good_Beta = NA,                 # P(p < AQL | Beta)
  P_Mass_Bad_Beta = NA,                  # P(p > LTPD | Beta)
  # Riesgos de la Línea Base (Plan Clásico, Densidad Uniforme)

  RP_clasic = NA,
  RC_clasic = NA,
  RAT_clasic = NA,
  # Plan Óptimo (Minimiza TWR bajo Información)
  n_opt = NA,
  c_opt = NA,
  RP_opt = NA,
  RC_opt = NA,
  RAT_opt = NA,

  # Ganancias (Plan Clásico Naive vs. Plan Óptimo Informado)
  RP_Ganancia = NA, 
  RC_Ganancia = NA,
  RAT_Ganancia = NA
)

# 4. Iterar sobre los 5 escenarios y calcular los riesgos

for (i in 1:dim(alpha_beta_params)[1]) { # i <- 1 + i
  
  # Usar los valores calculados de alpha y beta específicos para el proveedor i
  alpha_b_val <- alpha_beta_params[i, "alpha_b"]
  beta_b_val <- alpha_beta_params[i, "beta_b"]
  
  # I. Calcular la masa de probabilidad para el Prior Beta específico (Informado)
  prob_mass_beta <- calc_prob_mass(alpha_b = alpha_b_val, beta_b = beta_b_val, AQL, LTPD)
  resultados_riesgo[i, "P_Mass_Good_Beta"] <- prob_mass_beta["P_Good"]
  resultados_riesgo[i, "P_Mass_Bad_Beta"] <- prob_mass_beta["P_Bad"]
  
  k_p_ <- as.numeric(prob_mass_beta["P_Good"])
  k_c_ <- as.numeric((prob_mass_beta["P_Bad"]))
  
  # II. Calcular los riesgos para el PLAN CLÁSICO (Bajo el Prior Beta ESPECÍFICO)
  risks_clasic <- calc_wr(n = n_clasic, c = c_clasic, 
                          alpha_b = alpha_b_val, beta_b = beta_b_val, 
                          AQL, LTPD, k_p = k_p_, k_c = k_c_)
  
  resultados_riesgo[i, "WR_clasic"] <- risks_clasic["WR_val"]
  resultados_riesgo[i, "WRP_clasic"] <- risks_clasic["WRP_val"]
  resultados_riesgo[i, "WRC_clasic"] <- risks_clasic["WRC_val"]

  # --- III. Búsqueda del Plan Óptimo (Minimizar TWR Puro) para el Escenario i ---

  n_opt_found <- NA
  c_opt_found <- NA
  min_RAT <- k_p_ * alpha + k_c_ * beta
  cumple <- FALSE
  # Búsqueda exhaustiva: Itera n de 1 hasta n_clasic (máximo tamaño de muestra del plan clásico)
  for (n_ in 1:N) { # n_ <- 1 + n_
    for (c_ in 0:(n_ - 1)) { # c_ <- 0 + 1 + c_

      risks_opt <- calc_wr(n = n_, c = c_, 
                           alpha_b = alpha_b_val, beta_b = beta_b_val, 
                           AQL, LTPD, k_p = k_p_, k_c = k_c_)

      RAT_current <- risks_opt["WR_val"]

      # Usando la lógica de tu código: buscar un plan 'mejor' que el plan Naive
      if (RAT_current <= min_RAT) {
        n_opt_found <- n_
        c_opt_found <- c_
        cumple <- TRUE
        break
      }
    }
    if(cumple) break
  } 


# IV. Almacenar los resultados del plan óptimo (n_opt, c_opt)
resultados_riesgo[i, "n_opt"] <- n_opt_found
resultados_riesgo[i, "c_opt"] <- c_opt_found

# Recalcular los riesgos para el análisis final
if (!is.na(n_opt_found)) {
  risks_opt_final <- calc_wr(n = n_opt_found, c = c_opt_found, 
                             alpha_b = alpha_b_val, beta_b = beta_b_val, 
                             AQL, LTPD, k_p = k_p_, k_c = k_c_) 
  
  resultados_riesgo[i, "WR_opt"] <- risks_opt_final["WR_val"]
  resultados_riesgo[i, "WRP_opt"] <- risks_opt_final["WRP_val"]
  resultados_riesgo[i, "WRC_opt"] <- risks_opt_final["WRC_val"]
  
  # V. ¡Cálculo de Ganancias Añadido!
  # Ganancia = Riesgo Naive - Riesgo Óptimo

} else {
  resultados_riesgo[i, "WR_opt"] <- NA
  resultados_riesgo[i, "WRP_opt"] <- NA
  resultados_riesgo[i, "WRC_opt"] <- NA
  
  resultados_riesgo[i, "WRP_Ganancia"] <- NA
  resultados_riesgo[i, "WRC_Ganancia"] <- NA
  resultados_riesgo[i, "WR_Ganancia"] <- NA
  
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
  "RP_naive", "RP_opt", "RP_Ganancia",
  "RC_naive", "RC_opt", "RC_Ganancia"
)]

print(resultados_final_comparado)



# -------------------------------------------------------------------------
# VISUALIZACIÓN 1 (PARTE SUPERIOR): Curvas CO Agrupadas (1 Gráfico)
# -------------------------------------------------------------------------
delta_p <- 1e-5
p_range <- seq(0, 1, by = delta_p)
# Calcular la curva CO para el Plan Clásico (uniforme/naive) una sola vez
Pa_clasic_base <- Pa(n = n_clasic, c = c_clasic, p = p_range, N = N)
Pa_ga_base <- Pa(n = n_ga, c = c_ga, p = p_range, N = N)
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

plot(p_range, dbeta(p_range, alpha_beta_params$alpha_b[4], alpha_beta_params$beta_b[4] ))
abline(h =c(0,1), v =c(AQL, LTPD), lty = 2)


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
lines(p_range, Pa_ga_base, lty = 1, lwd = 2, col = optimal_colors[i])
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
=======
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
p1 <- c(0.95, 0.75, 0.50, 0.25, 0.10, 0.1, 0.3)
p2 <- c(1e-4, 0.05, 0.10, 0.25, 0.40, 0.89, 0.69) # 1 - PA(LTPD) -> Beta (Riesgo del Consumidor)

# Inicializar las listas maestras
lista_tablas_resultados <- list()
lista_parametros_beta   <- list()

for (esce in 1:dim(Esce)[1]) { # esce <- 1 + esce

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
  
  for (i in 1:length(p1)) { # i <- 1 + i
  
  # find_beta encuentra los parámetros de la distribución Beta
  Shape <- find_beta(x1 = Esce[esce, 4], p1 = p1[i], x2 = Esce[esce, 5], p2 = 1 - p2[i]) 
  alpha_beta_params[i, "alpha_b"] <- Shape$shape1
  alpha_beta_params[i, "beta_b"] <- Shape$shape2
  }
  
  # 3. Inicializar la tabla de resultados
  resultados_riesgo <- data.frame(
  Escenario = 1:length(p1),
  Proveedor = c(c("Excelente", "Bueno", "Regular", "Malo", "Muy Malo"), paste("otro", 1:(abs(5-length(p1))))),
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
  
  for (j in 1:length(p1)) { # j <- 1 + j
    
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
    
    n_opt_found <- NA
    c_opt_found <- NA
    des_WRT <- k_p_ * Esce[esce, 2] + k_c_ * Esce[esce, 3]
    cumple <- FALSE
    # Búsqueda exhaustiva
    for (n_ in 1:Esce[esce, 1]) { # n_ <- 1 + n_
      for (c_ in 0:(n_ - 1)) { # c_ <- 0 + 1 + c_
        
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
>>>>>>> 16a177554b09c394ca6f03b22b399b8fb6966908
