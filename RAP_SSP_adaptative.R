# ========================================================================= #
# SCRIPT 3: MOTOR DE BÚSQUEDA SECUENCIAL RAP (PARALELIZADO)
# ========================================================================= #

cat("\n=================================================================\n")
cat(" INICIANDO BÚSQUEDA SECUENCIAL PARALELIZADA (MODO TURBO)...\n")
cat("=================================================================\n")

# Cargar librerías necesarias
if(!require(doParallel)) install.packages("doParallel")
if(!require(foreach)) install.packages("foreach")

# 1. Preparar el Clúster (Los trabajadores)
nucleos_totales <- parallel::detectCores()
nucleos_a_usar <- nucleos_totales - 1 # Dejamos 1 libre para el sistema operativo
cat(sprintf(">>> Detectados %d núcleos. Usando %d para el proceso...\n", nucleos_totales, nucleos_a_usar))

cl <- parallel::makeCluster(nucleos_a_usar)
doParallel::registerDoParallel(cl)


# Inicializamos el data.frame vacío donde apilaremos todo
Matriz_General_ASSP <- data.frame()

# 1. Extraemos las estimaciones Bivariadas (TTB_a y TTB_b que son columnas 6 y 7)
# de los perfiles: Prov1 (fila 1), Prov2 (fila 2), Prov3_Esce_13 (15) y Prov3_Esce_66 (68)
datos_empiricos <- diccionario_emp_temp[c(1, 2, 15, 68), c(6, 7)]
colnames(datos_empiricos) <- c("a", "b")

# 2. Construimos el diccionario base de forma segura (sin coerción a texto)
diccionario_base <- data.frame(
  Perfil = c("naive", "Prov1", "Prov2", "Prov3_A", "Prov3_B"),
  a      = c(1, datos_empiricos$a),
  b      = c(1, datos_empiricos$b)
)

# Iteramos sobre los 72 escenarios contractuales de Esce_Global
for (esce in 1:nrow(Esce_Global)) {#esce <- 1+esce
  
  # Extraemos la fila del escenario actual
  escenario_actual <- Esce_Global[esce, ]
  
  # Filtramos las 7 prioris de Cook correspondientes al AQL y LTPD actuales
  cook_actual <- subset(diccionario_temp, 
                        AQL == escenario_actual$AQL & 
                          LTPD == escenario_actual$LTPD)
  
  # Nos quedamos solo con las columnas estructurales de Cook
  cook_actual <- cook_actual[, c("Perfil", "a", "b")]
  
  # Ensamblamos los 12 perfiles: 7 dinámicos (Cook) + 5 fijos (Empíricos + Naive)
  perfiles_escenario <- rbind(cook_actual, diccionario_base)
  
  # Cruzamos la información del contrato con los 12 perfiles
  bloque_iteracion <- data.frame(
    Id_Escenario = esce,
    N            = escenario_actual$N,
    alpha_des    = escenario_actual$alpha,
    beta_des     = escenario_actual$beta,
    AQL          = escenario_actual$AQL,
    LTPD         = escenario_actual$LTPD,
    n_clasico    = escenario_actual$n_clasico,
    c_clasico    = escenario_actual$c_clasico,
    Perfil       = perfiles_escenario$Perfil,
    a            = as.numeric(perfiles_escenario$a), # Aseguramos que sea numérico
    b            = as.numeric(perfiles_escenario$b)  # Aseguramos que sea numérico
  )
  
  # Apilamos el bloque en nuestra matriz maestra
  Matriz_General_ASSP <- rbind(Matriz_General_ASSP, bloque_iteracion)
}

# Reiniciamos los nombres de las filas para que vayan del 1 al 864
rownames(Matriz_General_ASSP) <- NULL

cat(">>> Matriz General construida con éxito.\n")
cat(">>> Dimensiones finales:", nrow(Matriz_General_ASSP), "filas y", ncol(Matriz_General_ASSP), "columnas.\n")

# 2. Pre-asignar las nuevas columnas en la matriz (vacías por ahora)
Matriz_General_ASSP$w_0 <- NA
Matriz_General_ASSP$w_1 <- NA
Matriz_General_ASSP$RAP_t_clasico <- NA
Matriz_General_ASSP$RAP_t_objetivo <- NA
Matriz_General_ASSP$n_opt <- NA
Matriz_General_ASSP$c_opt <- NA
Matriz_General_ASSP$RAP_t_opt <- NA

# 3. Ciclo Paralelizado (foreach)
# En lugar de modificar la matriz global directamente, cada núcleo procesa una fila 
# y al final (.combine = rbind) une todas las filas procesadas en una nueva matriz.
Matriz_Resultados <- foreach(i = 1:nrow(Matriz_General_ASSP), 
                             .combine = rbind, 
                             .packages = c("pracma"), # Cada núcleo necesita cargar 'pracma'
                             .export = c("calc_rap", "Pa", "rap_p", "rap_c", "dens_acum") # Funciones que los núcleos deben conocer
) %dopar% {
  
  # Copiamos la fila que este núcleo específico va a procesar
  fila <- Matriz_General_ASSP[i, ]
  
  # --- Extracción de variables ---
  N_lote   <- fila$N
  a_prior  <- fila$a
  b_prior  <- fila$b
  AQL_val  <- fila$AQL
  LTPD_val <- fila$LTPD
  alfa_des <- fila$alpha_des
  beta_des <- fila$beta_des
  
  # --- I. Calcular Masas de Probabilidad ---
  pesos <- dens_acum(a_b = a_prior, b_b = b_prior, AQL = AQL_val, LTPD = LTPD_val)
  w_0 <- pesos["PA_Good"] 
  w_1 <- pesos["PA_Bad"]  
  
  fila$w_0 <- w_0
  fila$w_1 <- w_1
  
  # Riesgo Objetivo
  rap_t_objetivo <- alfa_des * w_0 + beta_des * w_1
  fila$RAP_t_objetivo <- rap_t_objetivo
  
  # --- II. Evaluar el Plan Clásico ---
  riesgos_clasicos <- calc_rap(N_ = N_lote, n = fila$n_clasico, c = fila$c_clasico, 
                               a_b = a_prior, b_b = b_prior, 
                               AQL = AQL_val, LTPD = LTPD_val, w_p = w_0, w_c = w_1)
  
  fila$RAP_t_clasico <- riesgos_clasicos["RAP_t_val"]
  
  # --- III. Búsqueda Exhaustiva del Plan Óptimo ---
  n_opt_found <- NA
  c_opt_found <- NA
  rap_t_found <- NA
  cumple <- FALSE
  
  for (n_ in 1:N_lote) {
    for (c_ in 0:(n_ - 1)) {
      
      riesgos_opt <- calc_rap(N_ = N_lote, n = n_, c = c_, 
                              a_b = a_prior, b_b = b_prior, 
                              AQL = AQL_val, LTPD = LTPD_val, w_p = w_0, w_c = w_1)
      
      RAP_t_current <- riesgos_opt["RAP_t_val"]
      
      if (!is.na(RAP_t_current) && RAP_t_current <= rap_t_objetivo) {
        n_opt_found <- n_
        c_opt_found <- c_
        rap_t_found <- RAP_t_current
        cumple <- TRUE
        break 
      }
    }
    if (cumple) break 
  }
  
  # Guardar resultados en la fila
  fila$n_opt <- n_opt_found
  fila$c_opt <- c_opt_found
  fila$RAP_t_opt <- rap_t_found
  
  # Retornamos la fila completada para que se una al final
  return(fila)
}

# 4. Apagar el Clúster (¡Muy importante para liberar la memoria RAM!)
parallel::stopCluster(cl)

# Sobreescribimos nuestra matriz general con la matriz que ya tiene los resultados
Matriz_General_ASSP <- Matriz_Resultados
write.csv(Matriz_Resultados, "Matriz_Resultados.csv", row.names = FALSE)
cat("\n>>> Búsqueda Secuencial Paralelizada finalizada con éxito.\n")
