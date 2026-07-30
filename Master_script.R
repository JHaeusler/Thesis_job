# ============================================================================== #
# PANEL DE CONTROL MAESTRO - DISEÑO ADAPTATIVO DE PLANES DE MUESTREO (SOMA)
# ============================================================================== #
rm(list = ls())
setwd("~/Thesis_job") 

# --- 0. CARGA GLOBAL DE PAQUETES ---
paquetes <- c("readxl", "readr", "dplyr", "stats", "boot",
              "coda", "matrixStats", "AcceptanceSampling",
              "sjstats", "tidyr", "pracma")
for (p in paquetes) {
  if(!require(p, character.only = TRUE)){
    install.packages(p); require(p, character.only = TRUE)
  }
}

# --- 1. INTERRUPTORES LÓGICOS ---
EJECUTAR_OBJETIVO_1 <- TRUE   
EJECUTAR_OBJETIVO_2 <- FALSE  
GENERAR_FIGURAS     <- FALSE   

FUENTES_ACTIVAS <- c(1, 2, 3)

# --- 2. CARGA DE DATOS EMPÍRICOS ESTÁTICOS (FUENTES 1 Y 2) ---
data1_excel <- read_excel("Acceptance Sampling MIL-STD 105E for Quality Control.xlsx", sheet = "tab")
vec_data1 <<- data1_excel %>% select(Defect) %>% pull()

data2_excel <- read_excel("Jurnal Rekavasi.xlsx", sheet = "tabla")
vec_data2 <<- data2_excel %>% select(p) %>% pull()

# --- 3. PREÁMBULO: MALLA CONTRACTUAL Y SIMULACIÓN (FUENTE 3) ---
cat("\n>>> Inicializando Malla de Escenarios y Simulando Data 3...\n")

# 3.1 Crear Malla Base
Esce_Global <- expand.grid(N = 1000, alpha = c(0.01, 0.05), beta = c(0.05, 0.10, 0.20), 
                           AQL = c(0.01, 0.02, 0.05), LTPD = c(0.08, 0.10, 0.15, 0.20))
Esce_Global <- Esce_Global[Esce_Global$AQL < Esce_Global$LTPD, ]
rownames(Esce_Global) <- NULL

# 3.2 Adicionar columnas para el plan clásico (Kiermeier)
Esce_Global$n_clasico <- NA
Esce_Global$c_clasico <- NA

# 3.3 Crear Matriz para la simulación (Filas = Escenarios, Columnas = 35 Lotes)
lotes <- 35
Matriz_Data3 <<- matrix(NA, nrow = nrow(Esce_Global), ncol = lotes)
rownames(Matriz_Data3) <- paste0("Esce_", 1:nrow(Esce_Global))
colnames(Matriz_Data3) <- paste0("Lote_", 1:lotes)

# 3.4 Llenar las tablas
set.seed(060722)
for (k in 1:nrow(Esce_Global)) {
  
  # Calcular n y c para el escenario k
  plan_kier <- find.plan(PRP = c(Esce_Global$AQL[k], 1 - Esce_Global$alpha[k]), 
                         CRP = c(Esce_Global$LTPD[k], Esce_Global$beta[k]), 
                         N = Esce_Global$N[k], type = "hypergeom")
  
  Esce_Global$n_clasico[k] <- plan_kier$n
  Esce_Global$c_clasico[k] <- plan_kier$c
  
  # Simular los 35 lotes
  sum_n <- numeric(lotes)
  for (l in 1:lotes) {
    sample_n <- sample(ifelse(runif(Esce_Global$N[k]) <= Esce_Global$AQL[k], 1, 0), plan_kier$n)
    sum_n[l] <- sum(sample_n)
  }
  
  # Guardar proporciones estimadas directamente en la fila k de la matriz
  Matriz_Data3[k, ] <- sum_n / plan_kier$n
}

# 3.5 Identificación de Escenarios Extremos para Simulación (Prov 3A y 3B)
cat(">>> Identificando escenarios de simulación extremos...\n")

Malla_Extremos <- Esce_Global %>% 
  mutate(id_escenario = row_number(),
         dif = LTPD - AQL)

# Prov 3A: Menor esfuerzo muestral
prov3_A <- Malla_Extremos %>% 
  arrange(n_clasico, c_clasico, dif) %>% 
  slice(1) 

# Prov 3B: Mayor esfuerzo muestral
prov3_B <- Malla_Extremos %>% 
  arrange(desc(n_clasico), desc(c_clasico), dif) %>% 
  slice(1) 

# Guardar los IDs en el entorno global para que los otros scripts los vean
idx_min <<- prov3_A$id_escenario
idx_max <<- prov3_B$id_escenario

cat(sprintf("    - Prov 3A (Mínimo): Escenario %d (n = %d, c = %d)\n", idx_min, prov3_A$n_clasico, prov3_A$c_clasico))
cat(sprintf("    - Prov 3B (Máximo): Escenario %d (n = %d, c = %d)\n", idx_max, prov3_B$n_clasico, prov3_B$c_clasico))

# Exportar Esce_Global al entorno global
Esce_Global <<- Esce_Global

cat(sprintf(">>> Preámbulo completado. Esce_Global (%d filas) y Matriz_Data3 listas en memoria.\n", nrow(Esce_Global)))

# ============================================================================== #
# A PARTIR DE AQUÍ VENDRÁN LOS LLAMADOS A LOS DEMÁS SCRIPTS (source...)
# ============================================================================== #

# --- 4. OBJETIVO 1: ELICITACIÓN Y EXTRACCIÓN ---
if (EJECUTAR_OBJETIVO_1) {
  cat("\n=================================================================\n")
  cat(" OBJETIVO 1A: DICCIONARIO COOK (2010)\n")
  cat("=================================================================\n")
  source("elicit_cook_2010.R", encoding = "UTF-8")
  
  cat("\n=================================================================\n")
  cat(" OBJETIVO 1B: EXTRACCIÓN EMPÍRICA MCMC\n")
  cat("=================================================================\n")
  diccionario_empirico <<- data.frame()
  source("MCMC_TTB.R", encoding = "UTF-8")
  
  if (GENERAR_FIGURAS) source("generar_atlas_figuras.R", encoding = "UTF-8")
}

# --- 5. OBJETIVO 2: BÚSQUEDA DEL PLAN ÓPTIMO (RAP) ---
if (EJECUTAR_OBJETIVO_2) {
  cat("\n=================================================================\n")
  cat(" OBJETIVO 2: BÚSQUEDA DETERMINÍSTICA RAP\n")
  cat("=================================================================\n")
  source("RAP_SSP_adaptative.R", encoding = "UTF-8")
}