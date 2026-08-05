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
# EJECUTAR_OBJETIVO_1 <- TRUE   
# EJECUTAR_OBJETIVO_2 <- FALSE  
# GENERAR_FIGURAS     <- FALSE   
# 
# FUENTES_ACTIVAS <- c(1, 2, 3)

# # --- 2. CARGA DE DATOS EMPÍRICOS ESTÁTICOS (FUENTES 1 Y 2) ---
# data1_excel <- read_excel("Acceptance Sampling MIL-STD 105E for Quality Control.xlsx", sheet = "tab")
# vec_data1 <<- data1_excel %>% select(Defect) %>% pull()
# 
# data2_excel <- read_excel("Jurnal Rekavasi.xlsx", sheet = "tabla")
# vec_data2 <<- data2_excel %>% select(p) %>% pull()

N = 1000
alpha_des = c(0.01, 0.05)
beta_des = c(0.05, 0.10, 0.20)
AQL = c(0.01, 0.02, 0.05)
LTPD = c(0.08, 0.10, 0.15, 0.20)

# 3.1 Crear Malla Base
Esce_Global <- expand.grid(N = N, alpha = alpha_des, beta = beta_des, 
                           AQL = AQL, LTPD = LTPD)
Esce_Global <- Esce_Global[Esce_Global$AQL < Esce_Global$LTPD, ]
rownames(Esce_Global) <- NULL

# 3.2 Adicionar columnas para el plan clásico (Kiermeier)
Esce_Global$n_clasico <- NA
Esce_Global$c_clasico <- NA

# 2. Opción A: Base R con mapply (Rápido, limpio y sin librerías extra)
# Pasa las columnas directamente como vectores
res <- mapply(
  FUN = calcular_kiermeier,
  N = Esce_Global$N,
  AQL = Esce_Global$AQL,
  LTPD = Esce_Global$LTPD,
  alpha = Esce_Global$alpha_des,
  beta = Esce_Global$beta_des
)

# 3. Asignar resultados transponiendo la matriz de retorno
Esce_Global$n_clasico <- res["n", ]
Esce_Global$c_clasico <- res["c", ]


# --- 4. OBJETIVO 1: ELICITACIÓN Y EXTRACCIÓN ---
source("elicit_cook_2010.R", encoding = "UTF-8")
source("MCMC_TTB.R", encoding = "UTF-8")
source("RAP_SSP_adaptative.R", encoding = "UTF-8")
