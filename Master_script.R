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



source("Clasic_plans_kiermeier.R", encoding = "UTF-8")



# --- 4. OBJETIVO 1: ELICITACIÓN Y EXTRACCIÓN ---
source("elicit_cook_2010.R", encoding = "UTF-8")
source("MCMC_TTB.R", encoding = "UTF-8")
source("RAP_SSP_adaptative.R", encoding = "UTF-8")
