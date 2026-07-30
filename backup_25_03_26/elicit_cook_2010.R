# ============================================================================== #
# SCRIPT 2: ELICITACIÓN DE EXPERTOS Y CURVAS CO MULTIESCENARIO
# Metodología: Método de Cook (2010) y Malla de Escenarios
# Dependencia: Variables globales definidas en 'Master_script.R'
# Autor: Juan Sebastián Haeusler
# ============================================================================== #

# --- 1. Carga de Paquetes ---
if(!require(sjstats)) install.packages("sjstats")
if(!require(AcceptanceSampling)) install.packages("AcceptanceSampling")
if(!require(dplyr)) install.packages("dplyr")
if(!require(tidyr)) install.packages("tidyr")

# --- Control de Herencia del Master Script ---

if(!exists("AQL_vec")) {
  cat("\n[!] Ejecución aislada detectada. Cargando malla paramétrica por defecto...\n")
  N_vec     <- 1000                      
  alpha_vec <- c(0.01, 0.05)             
  beta_vec  <- c(0.05, 0.10, 0.20)       
  AQL_vec   <- c(0.01, 0.02, 0.05)       
  LTPD_vec  <- c(0.08, 0.10, 0.15, 0.20)
}

# ============================================================================== #
# --- 2. PARTE A: CREACIÓN DEL DICCIONARIO DE PRIORIS (PARA EL CAPÍTULO 5) ---
# ============================================================================== #

# Malla exclusiva para Cook (Solo requiere AQL y LTPD)
Esce_Cook <- expand.grid(AQL = AQL_vec, LTPD = LTPD_vec)
Esce_Cook <- Esce_Cook[Esce_Cook$AQL < Esce_Cook$LTPD, ]
rownames(Esce_Cook) <- NULL

# Perfiles del Proveedor
p1_val <- c(0.95, 0.75, 0.50, 0.25, 0.10, 0.05, 0.01)
p2_val <- c(1e-4, 0.05, 0.10, 0.25, 0.40, 0.05, 0.01)
nombres_perfiles <- c("Excelente", "Bueno", "Regular", "Malo", "Muy Malo",
                      "ZI_90", "ZI_98")

cat("Generando Diccionario de Prioris para el Capítulo 5...\n")

# Variable temporal para guardar el formato largo
diccionario_temp <- data.frame()

for (k in 1:nrow(Esce_Cook)) {
  aql_c <- Esce_Cook$AQL[k]
  ltpd_c <- Esce_Cook$LTPD[k]
  
  for(i in 1:length(p1_val)) {
    res <- find_beta(x1 = aql_c, p1 = p1_val[i], x2 = ltpd_c, p2 = 1 - p2_val[i])
    fila <- data.frame(
      AQL = aql_c, LTPD = ltpd_c, Perfil = nombres_perfiles[i],
      a = round(res$shape1, 4), b = round(res$shape2, 4)
    )
    diccionario_temp <- rbind(diccionario_temp, fila)
  }
}

# --- 3. CONSOLIDACIÓN DE MEMORIA GLOBAL Y SALIDA LATEX ---

# Sobreescribimos diccionario_cook en formato ANCHO (Requisito indispensable para RAP)
diccionario_cook <<- diccionario_temp %>%
  pivot_wider(names_from = Perfil, values_from = c(a, b), names_glue = "{Perfil}_{.value}") %>%
  select(AQL, LTPD, Excelente_a, Excelente_b, Bueno_a, Bueno_b, 
         Regular_a, Regular_b, Malo_a, Malo_b, `Muy Malo_a`, `Muy Malo_b`,
         ZI_90_a, ZI_90_b, ZI_98_a, ZI_98_b) %>%
  arrange(AQL, LTPD)

cat("\n=== TABLA PARA EL CAPÍTULO 5 (Diccionario Cook) ===\n")
print(as.data.frame(diccionario_cook), row.names = FALSE)