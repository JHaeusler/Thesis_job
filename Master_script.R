# ============================================================================== #
# PANEL DE CONTROL MAESTRO - DISEÑO ADAPTATIVO DE PLANES DE MUESTREO (SOMA)
# Autor: Juan Sebastián Haeusler
# Tesis de Maestría - Optimización Bayesiana de Calidad
# ============================================================================== #

# Limpiar el entorno de trabajo (Solo se hace aquí en el Master)
rm(list = ls())

# --- 1. CONFIGURACIÓN DEL USUARIO (PANEL DE CONTROL) ---
setwd("~/Thesis_job") 

# A. Selección MÚLTIPLE de Fuentes de Datos Empíricos para el MCMC
# 1 = Defiatri & Damayanti (2023)
# 2 = Isnanto et al. (2019)
# 3 = Simulación Estocástica Propia
# ¡El usuario puede colocar c(1), c(1,2) o c(1,2,3)!
FUENTES_ACTIVAS <- c(1, 2, 3) 

# B. Parámetros del Lote y Malla de Escenarios a evaluar
N_vec     <- 1000                      
alpha_vec <- c(0.01, 0.05)             
beta_vec  <- c(0.05, 0.10, 0.20)       
AQL_vec   <- c(0.01, 0.02, 0.05)       
LTPD_vec  <- c(0.08, 0.10, 0.15, 0.20) 

# Homologación de variables para los scripts esclavos
N <- N_vec; alpha <- alpha_vec; beta <- beta_vec; AQL <- AQL_vec; LTPD <- LTPD_vec

cat("=================================================================\n")
cat(" INICIANDO SISTEMA DE OPTIMIZACIÓN DE MUESTREO ADAPTATIVO (SOMA)\n")
cat("=================================================================\n\n")

# --- 2. ELICITACIÓN DE EXPERTOS (Módulo Independiente de los Datos) ---
cat(">> EJECUTANDO ELICITACIÓN DE EXPERTOS (elicit_cook_2010.R)...\n")
# Corre 1 sola vez porque la grilla AQL/LTPD es la misma para todas las fuentes
source("elicit_cook_2010.R", encoding = "UTF-8")
cat(">> ELICITACIÓN COMPLETADA Y ATLAS PDF GENERADO.\n\n")

# --- 3. BUCLE DE ANÁLISIS POR FUENTE DE DATOS ---
# Lista maestra que guardará los planes óptimos de todos los proveedores
ANALISIS_GLOBAL <- list()

for(fuente in FUENTES_ACTIVAS) {
  
  # La variable FUENTE_ACTIVA informará al MCMC qué Excel debe cargar
  FUENTE_ACTIVA <- fuente 
  
  cat(sprintf("\n#################################################################\n"))
  cat(sprintf("   PROCESANDO FUENTE DE DATOS: %d\n", FUENTE_ACTIVA))
  cat(sprintf("#################################################################\n"))
  
  # MÓDULO I: Estimación Bayesiana (MCMC) para la fuente actual
  cat("\n>> FASE I: ESTIMACIÓN ESTOCÁSTICA (MCMC_TTB.R)...\n")
  source("MCMC_TTB.R", encoding = "UTF-8")
  
  # MÓDULO II: Búsqueda del Plan Óptimo (RAP) para la fuente actual
  cat("\n>> FASE II: BÚSQUEDA DETERMINÍSTICA RAP (RAP_SSP_adaptative.R)...\n")
  source("RAP_SSP_adaptative.R", encoding = "UTF-8")
  
  # Guardar los resultados de los escenarios de esta fuente específica
  nombre_lista <- paste0("Fuente_Datos_", FUENTE_ACTIVA)
  ANALISIS_GLOBAL[[nombre_lista]] <- lista_tablas_resultados
  
  cat(sprintf("\n>> Análisis de la Fuente %d guardado exitosamente.\n", FUENTE_ACTIVA))
}

# ============================================================================== #
# --- 4. CONSOLIDACIÓN FINAL ---
# ============================================================================== #
cat("\n=================================================================\n")
cat(" PROCESO BATCH FINALIZADO CON ÉXITO\n")
cat("=================================================================\n")
cat("Puede explorar el objeto 'ANALISIS_GLOBAL' en su entorno de RStudio\n")
cat("para comparar cómo varían los tamaños de muestra (n) óptimos según\n")
cat("el comportamiento histórico de cada empresa.\n")