#### Simulación de inspección de lotes para atributos, comparando dos
#### distribuciones para la proporción de defectuosos: Uniforme y Beta

#### paquetes
# Lista de paquetes necesarios
lista_de_paquetes <- c("pracma", "GA", "AcceptanceSampling", "tidyverse")

# Itera sobre la lista e instala/carga los paquetes que falten
for (paquete in lista_de_paquetes) {
  if(!require(paquete, character.only = TRUE)){
    install.packages(paquete)
    require(paquete, character.only = TRUE)
  }
}

# semilla para reproducibilidad
set.seed(123)

# Parámetros de entrada de la simulación
N <- 200 # Tamaño del lote
AQL <- 0.05 # Límite de calidad aceptable (Acceptable Quality Level)
LTPD <- 0.10 # Límite de porcentaje de tolerancia de lote defectuoso (Lot Tolerance Percent Defective)
alpha_des <- 0.05 # Riesgo del productor
beta_des <- 0.10 # Riesgo del consumidor

# Parámetros para la distribución Uniforme de la proporción de defectuosos
a_unif <- 0; b_unif <- 0.25

# Parámetros para la distribución Beta de la proporción de defectuosos
# (Hemos elegido estos para simular una distribución sesgada hacia lotes de buena calidad)
a_beta <- 8; b_beta <- 105

# parametros para la distribucion hypergeometrica

a_hyper <- rbinom(1, 200, AQL); b_hyper <- N - a_hyper


# Cantidad de lotes a inspeccionar en la simulación
k <- 100000

# plan de muestreo clásico con AcceptanceSampling
plan_clasic <- find.plan(PRP = c(AQL, 1 - alpha_des),
                         CRP = c(LTPD, beta_des),
                         N=N, type = "hypergeom")

# n: tamaño de la muestra; c: número de aceptación
n_muestra <- plan_clasic$n
c_aceptacion <- plan_clasic$c
print(paste("Plan de muestreo: n =", n_muestra, ", c =", c_aceptacion))

# Vectores para almacenar los resultados de la simulación con Distribución Uniforme
p_est_unif <- vector("numeric", k)
veredicto_lote_unif <- vector("character", k)
p_lote_unif_real <- vector("numeric", k) # Guardamos la p_lote real para el análisis

# Vectores para almacenar los resultados de la simulación con Distribución Beta
p_est_beta <- vector("numeric", k)
veredicto_lote_beta <- vector("character", k)
p_lote_beta_real <- vector("numeric", k) # Guardamos la p_lote real para el análisis


# Vectores para almacenar los resultados de la simulación con Distribución hipergeomtrica
p_est_beta <- vector("numeric", k)
veredicto_lote_hyper <- vector("character", k)
p_lote_hyper_real <- vector("numeric", k) # Guardamos la p_lote real para el análisis


# Bucle principal de la simulación
for (i in 1:k) {
  
  # --- Simulación con Distribución Uniforme ---
  p_lote_unif_real[i] <- runif(1, a_unif, b_unif)
  sim_lote_binom_unif <- rbinom(N, size = 1, p = p_lote_unif_real[i]) 
  def_muestra_unif <- sim_lote_binom_unif[sample(1:N, n_muestra, replace = FALSE)]
  p_est_unif[i] <- sum(def_muestra_unif) / n_muestra
  
  # --- Simulación con Distribución hipergeometrica ---
  p_lote_hyper_real[i] <- runif(1, a_hyper, b_hyper, plan_clasic$n)
  sim_lote_binom_hyper <- rbinom(N, size = 1, p = p_lote_unif_real[i])  
  def_muestra_unif <- sim_lote_binom_unif[sample(1:N, n_muestra, replace = FALSE)]
  p_est_unif[i] <- sum(def_muestra_unif) / n_muestra
  
  
  if (sum(def_muestra_unif) <= c_aceptacion) {
    veredicto_lote_unif[i] <- "Aceptado"
  } else {
    veredicto_lote_unif[i] <- "Rechazado"
  }
  
  # --- Simulación con Distribución Beta ---
  p_lote_beta_real[i] <- rbeta(1, a_beta, b_beta)
  sim_lote_binom_beta <- rbinom(N, size = 1, p = p_lote_beta_real[i]) 
  def_muestra_beta <- sim_lote_binom_beta[sample(1:N, n_muestra, replace = FALSE)]
  p_est_beta[i] <- sum(def_muestra_beta) / n_muestra
  
  if (sum(def_muestra_beta) <= c_aceptacion) {
    veredicto_lote_beta[i] <- "Aceptado"
  } else {
    veredicto_lote_beta[i] <- "Rechazado"
  }
}

# ---
# Análisis de los resultados
# ---

# Configuración de los gráficos para que se muestren 2x2
par(mfrow = c(2, 2))

# Histograma para la Distribución Uniforme
hist(p_est_unif, main = "Distribución Uniforme",
     xlab = "Proporción estimada", ylab = "Frecuencia", col = "skyblue", border = "black")

# Histograma para la Distribución Beta
hist(p_est_beta, main = "Distribución Beta",
     xlab = "Proporción estimada", ylab = "Frecuencia", col = "salmon", border = "black")

# Gráfico de barras para los resultados de la simulación Uniforme
conteo_unif <- table(veredicto_lote_unif)
barplot(conteo_unif,
        main = "Veredicto - Distribución Uniforme",
        xlab = "Resultado", ylab = "Número de Lotes",
        col = c("tomato", "green"), ylim = c(0, k), border = "black")

# Gráfico de barras para los resultados de la simulación Beta
conteo_beta <- table(veredicto_lote_beta)
barplot(conteo_beta,
        main = "Veredicto - Distribución Beta",
        xlab = "Resultado", ylab = "Número de Lotes",
        col = c("tomato", "green"), ylim = c(0, k), border = "black")

# Resumen de errores para la Distribución Uniforme
falsos_positivos_unif <- sum(p_lote_unif_real <= AQL & veredicto_lote_unif == "Rechazado")
falsos_negativos_unif <- sum(p_lote_unif_real >= LTPD & veredicto_lote_unif == "Aceptado")
verdaderos_positivos_unif <- sum(p_lote_unif_real >= LTPD & veredicto_lote_unif == "Rechazado")
verdaderos_negativos_unif <- sum(p_lote_unif_real <= AQL & veredicto_lote_unif == "Aceptado")
lotes_indefinidos_unif <- sum(p_lote_unif_real > AQL & p_lote_unif_real < LTPD)

print("---")
print("Análisis de errores para la Distribución Uniforme:")
print(paste("Fracción de falsos positivos:", round(falsos_positivos_unif / k, 4)))
print(paste("Fracción de falsos negativos:", round(falsos_negativos_unif / k, 4)))
print("---")

# Resumen de errores para la Distribución Beta
falsos_positivos_beta <- sum(p_lote_beta_real <= AQL & veredicto_lote_beta == "Rechazado")
falsos_negativos_beta <- sum(p_lote_beta_real >= LTPD & veredicto_lote_beta == "Aceptado")
verdaderos_positivos_beta <- sum(p_lote_beta_real >= LTPD & veredicto_lote_beta == "Rechazado")
verdaderos_negativos_beta <- sum(p_lote_beta_real <= AQL & veredicto_lote_beta == "Aceptado")
lotes_indefinidos_beta <- sum(p_lote_beta_real > AQL & p_lote_beta_real < LTPD)

print("---")
print("Análisis de errores para la Distribución Beta:")
print(paste("Fracción de falsos positivos:", round(falsos_positivos_beta / k, 4)))
print(paste("Fracción de falsos negativos:", round(falsos_negativos_beta / k, 4)))
print("---")
