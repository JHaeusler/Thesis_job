cases <- c(13:24, 35:46, 59:70, 66:77)
length(cases)
# 1:nrow(Matriz_General_ASSP)

# Se omite el .export para evitar el Warning, doParallel lo detecta automáticamente
Lista_Matrices_RAP <- foreach(i = cases, 
                              .packages = c("pracma")
) %dopar% { # i <- 1 + i
  
  fila <- Matriz_General_ASSP[i, ]
  
  # Forzamos N=200 
  N_lote <- as.numeric(unname(fila$N))
  
  # --- ESCUDO 1: DESINFECCIÓN DIMENSIONAL EXHAUSTIVA ---
  # Garantizamos que todo sea un escalar numérico puro (double) sin atributos
  a_prior  <- as.numeric(unname(fila$a))
  b_prior  <- as.numeric(unname(fila$b))
  AQL_val  <- as.numeric(unname(fila$AQL))
  LTPD_val <- as.numeric(unname(fila$LTPD))
  
  w_0 <- as.numeric(unname(fila$w_0))
  w_1 <- as.numeric(unname(fila$w_1))
  
  n_rap <- as.numeric(unname(fila$n_clasico))
  
  # --- 1. CREACIÓN AUTOMÁTICA DE LA MATRIZ (N x N) ---
  Z_matrix <- matrix(NA, nrow = n_rap, ncol = n_rap)
  rownames(Z_matrix) <- paste0("n=", 1:n_rap)
  colnames(Z_matrix) <- paste0("c=", 0:(n_rap - 1))
  
  
  # --- 2. LLENADO CON AMORTIGUADOR DE FALLOS (tryCatch) ---
  for(n_ in 1:n_rap) {   # n_ <- 1 + n_          
    for(c_ in 0:(n_ - 1)) {         # c_ <- 0 + 1 + c_
      
      # ESCUDO 2: TRATAMIENTO DE ERRORES DE INTEGRACIÓN
      RAP_evaluado <- tryCatch({
        
        # Intenta calcular la integral
        riesgos_opt <- calc_rap(N_ = N_lote, n = n_, c = c_, 
                                a_b = a_prior, b_b = b_prior, 
                                AQL = AQL_val, LTPD = LTPD_val, w_p = w_0, w_c = w_1)
        # Retorna el valor
        as.numeric(riesgos_opt["RAP_t_val"])
        
      }, error = function(e) {
        # Si pracma colapsa en esta celda específica, imprime silencio y retorna NA
        return(NA)
      })
      
      Z_matrix[n_, c_ + 1] <- RAP_evaluado
    }
  }
  
  return(Z_matrix)
}

# Definimos la paleta de colores
paleta_calor <- colorRampPalette(c("forestgreen", "khaki1", "firebrick"))(50)

# Opcional: Ajustar la ventana gráfica (ej. 2x2 para no saturar)
# par(mfrow = c(2, 2))

# Iteramos sobre la LONGITUD de los casos (1 hasta 48)
for(idx in 1:length(cases)) { # idx <- 1 + idex
  
  # Extraemos el ID real de la fila (13, 14, ..., 77)
  fila_real <- Matriz_General_ASSP[cases[idx], ]
  
  # Extraemos la matriz correspondiente a esta iteración (de 1 a 48)
  matriz_raps <- Lista_Matrices_RAP[[idx]]
  
  # Extraemos metadatos para el título dinámico
  nombre_perfil <- Matriz_General_ASSP$Perfil[cases[idx]]
  escenario_id  <- Matriz_General_ASSP$Id_Escenario[cases[idx]]
  
  # --- 1. BÚSQUEDA DEL PLAN ÓPTIMO ---
  # Encontramos el valor RAP mínimo (omitiendo los NA del triángulo superior)
  rap_minimo <- min(matriz_raps, na.rm = TRUE)
  
  # Encontramos las coordenadas en la matriz (arr.ind = TRUE devuelve fila y columna)
  coord_opt <- which(matriz_raps == rap_minimo, arr.ind = TRUE)
  
  # Si hay empates (mismo RAP), tomamos el primero (el de menor n y menor c)
  n_opt <- fila_real$n_opt
  c_opt <- fila_real$c_opt
  
  
  # Definimos las secuencias de los ejes
  n_seq <- 1:nrow(matriz_raps)
  c_seq <- 0:(ncol(matriz_raps) - 1)
  
  # --- 2. CREACIÓN DEL GRÁFICO ---
  # Título compuesto con el plan óptimo detectado
  titulo_grafico <- sprintf("Heatmap RAP - %s (Escenario %s)\nÓptimo: n* = %d, c* = %d (RAP = %.4f)", 
                            nombre_perfil, escenario_id, n_opt, c_opt, rap_minimo)
  
  # Mapa de calor base
  image(x = n_seq, 
        y = c_seq, 
        z = matriz_raps, 
        col = paleta_calor,
        xlab = "Tamaño de Muestra (n)", 
        ylab = "Número de Aceptación (c)",
        main = titulo_grafico,
        cex.main = 0.9) # Reduce un poco el tamaño del título para que quepa bien
  
  # --- 3. CAPAS DE CONTORNO (CONTOUR) ---
  # Contornos de nivel generales (topología)
  contour(x = n_seq, y = c_seq, z = matriz_raps, 
          add = TRUE, col = rgb(1, 1, 1, 0.5), nlevels = 8)
  
  # Frontera óptima gruesa negra (donde z == rap_minimo)
  contour(x = n_seq, y = c_seq, z = matriz_raps, 
          add = TRUE, levels = rap_minimo, 
          col = "black", lwd = 3, lty = 1, drawlabels = FALSE)
  
  # --- 4. PUNTO ÓPTIMO Y LEYENDA ---
  # Dibujamos un punto azul conspicuo en la coordenada exacta del diseño óptimo
  points(x = n_opt, y = c_opt, pch = 16, col = "blue", cex = 1.5)
  
  legend("topleft", 
         legend = c("Frontera Óptima", "Plan Óptimo (n*, c*)"), 
         col = c("black", "blue"), 
         lwd = c(3, NA), pch = c(NA, 19), lty = c(1, NA), 
         bg = "white", cex = 0.8)
}

# Restaurar par() si lo usó
# par(mfrow = c(1, 1))