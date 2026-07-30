orden_perfiles <- c("Excelente", "Malo", "ZI_98", "Prov1", "Prov2", "Prov3_A")

# Escenario 1 - Proveedor Asfixiante AQL = 0.05 | LTPD = 0.08 | alpha = 0.01 | beta = 0.05
esce_asfix <- Matriz_General_ASSP %>% 
  filter(alpha_des == 0.01 & beta_des == 0.05 & AQL == 0.05 & LTPD == 0.08) %>% 
  filter(Perfil %in% orden_perfiles) %>%
  arrange(match(Perfil, orden_perfiles))

# Escenario 2 - Proveedor estandar AQL = 0.05 | LTPD = 0.08 | alpha = 0.01 | beta = 0.05
esce_est <- Matriz_General_ASSP %>% 
  filter(alpha_des == 0.05 & beta_des == 0.10 & AQL == 0.01 & LTPD == 0.10) %>% 
  filter(Perfil %in% orden_perfiles) %>%
  arrange(match(Perfil, orden_perfiles))

# Escenario 3 - Proveedor laxo AQL = 0.05 | LTPD = 0.08 | alpha = 0.01 | beta = 0.05
esce_laxo <- Matriz_General_ASSP %>% 
  filter(alpha_des == 0.05 & beta_des == 0.20 & AQL == 0.01 & LTPD == 0.20) %>% 
  filter(Perfil %in% orden_perfiles) %>%
  arrange(match(Perfil, orden_perfiles))

Esce_tot <- rbind(esce_asfix, esce_est, esce_laxo)

Lista_esce_RAP <- foreach(i = 1:nrow(Esce_tot), 
                          .packages = c("pracma")
) %dopar% { # i <- 1 + i
  
  fila <- Esce_tot[i, ]
  
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

# ============================================================================== #
# FASE DE VISUALIZACIÓN: EXPORTACIÓN AUTOMÁTICA DE PANELES 2x3 (PNG ALTA CALIDAD)
# ============================================================================== #

# Creamos el directorio si no existe
if(!dir.exists("imagenes")) dir.create("imagenes")

# --- 1. FUNCIÓN DE GRAFICACIÓN CON EXPORTACIÓN PNG ---
graficar_panel_escenario <- function(indices_rango, lista_matrices, datos_generales, archivo_salida = NULL) {
  
  # Si el usuario proporciona una ruta, abrimos el dispositivo PNG a 300 dpi (Alta Resolución)
  if(!is.null(archivo_salida)) {
    png(archivo_salida, width = 4400, height = 2700, res = 300)
  }
  
  # Configurar la ventana gráfica (2 filas x 3 columnas) con espacio inferior para la leyenda
  par_original <- par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(4, 0, 3, 0))
  paleta_calor <- colorRampPalette(c("forestgreen", "khaki1", "firebrick"))(50)
  
  # Extraer metadatos globales del escenario
  idx_base <- indices_rango[1]
  alfa_esce <- datos_generales$alpha_des[idx_base]
  beta_esce <- datos_generales$beta_des[idx_base]
  aql_esce  <- datos_generales$AQL[idx_base]
  ltpd_esce <- datos_generales$LTPD[idx_base]
  
  # --- BUCLE DE GRÁFICAS ---
  for(idx in indices_rango) { 
    
    matriz_raps <- lista_matrices[[idx]]
    fila_datos <- datos_generales[idx, ]
    nombre_perfil <- fila_datos$Perfil
    
    # Cálculo del umbral real de decisión (RAP Objetivo)
    alfa_des <- as.numeric(fila_datos$alpha_des)
    beta_des <- as.numeric(fila_datos$beta_des)
    w_0      <- as.numeric(fila_datos$w_0)
    w_1      <- as.numeric(fila_datos$w_1)
    
    rap_t_objetivo <- (alfa_des * w_0) + (beta_des * w_1)
    
    # Búsqueda del plan óptimo (Menor n)
    cumple_condicion <- which(matriz_raps <= rap_t_objetivo, arr.ind = TRUE)
    
    if(nrow(cumple_condicion) > 0) {
      cumple_condicion <- cumple_condicion[order(cumple_condicion[, "row"], cumple_condicion[, "col"]), , drop = FALSE]
      n_optimo <- cumple_condicion[1, "row"]
      c_optimo <- cumple_condicion[1, "col"] - 1 
    } else {
      n_optimo <- NA; c_optimo <- NA
    }
    
    n_seq <- 1:nrow(matriz_raps)
    c_seq <- 0:(ncol(matriz_raps) - 1)
    
    titulo_grafico <- sprintf("%s\nRAP Obj: %.4f | Óptimo (n*=%d, c*=%d)", 
                              nombre_perfil, rap_t_objetivo, n_optimo, c_optimo)
    
    # Renderizado del plot (Mapa, topología y límite estricto)
    image(x = n_seq, y = c_seq, z = matriz_raps, 
          col = paleta_calor, xlab = "Tamaño de Muestra (n)", ylab = "Número de Aceptación (c)",
          main = titulo_grafico,
          cex.main = 1.5,   
          cex.lab = 1.5,    
          cex.axis = 1.5) 
    
    contour(x = n_seq, y = c_seq, z = matriz_raps, 
            add = TRUE, col = rgb(1, 1, 1, 0.5), nlevels = 10)
    
    contour(x = n_seq, y = c_seq, z = matriz_raps, 
            add = TRUE, levels = rap_t_objetivo, 
            col = "black", lwd = 2, lty = 2, drawlabels = FALSE)
    
    if(!is.na(n_optimo)) {
      points(x = n_optimo, y = c_optimo, pch = 18, col = "blue", cex = 1.3)
    }
  }
  
  # --- TÍTULO GLOBAL SUPERIOR ---
  titulo_global <- sprintf("Topología RAP - AQL: %.2f | LTPD: %.2f | \u03b1: %.2f | \u03b2: %.2f", 
                           aql_esce, ltpd_esce, alfa_esce, beta_esce)
  mtext(titulo_global, side = 3, outer = TRUE, cex = 1.3, font = 2)
  
  # --- LEYENDA GLOBAL INFERIOR ---
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", legend = c("Frontera de Cumplimiento (RAP \u2264 Obj)", "Plan Óptimo (Menor n)"), 
         col = c("black", "blue"), lwd = c(4, NA), pch = c(NA, 19), 
         lty = c(2, NA), cex = 1.4, horiz = TRUE, bty = "n", inset = c(0, 0.01))
  
  # Restaurar parámetros gráficos a la normalidad
  par(par_original)
  
  # Si exportamos a PNG, cerramos el lienzo para que el archivo se guarde correctamente
  if(!is.null(archivo_salida)) {
    dev.off()
    cat(sprintf("   -> Exportado exitosamente: %s\n", archivo_salida))
  }
}

# --- 2. GENERACIÓN Y EXPORTACIÓN DE LOS 3 PANELES ---

cat("\n>>> Generando Panel 1: Escenario Asfixiante...\n")
graficar_panel_escenario(indices_rango = 1:6, 
                         lista_matrices = Lista_esce_RAP, 
                         datos_generales = Esce_tot,
                         archivo_salida = "imagenes/Panel_RAP_Esce_Asfixiante.png")

cat(">>> Generando Panel 2: Escenario Estándar...\n")
graficar_panel_escenario(indices_rango = 7:12, 
                         lista_matrices = Lista_esce_RAP, 
                         datos_generales = Esce_tot,
                         archivo_salida = "imagenes/Panel_RAP_Esce_Estandar.png")

cat(">>> Generando Panel 3: Escenario Laxo...\n")
graficar_panel_escenario(indices_rango = 13:18, 
                         lista_matrices = Lista_esce_RAP, 
                         datos_generales = Esce_tot,
                         archivo_salida = "imagenes/Panel_RAP_Esce_Laxo.png")

cat("\n>>> ¡Ciclo de visualización y exportación completado!\n")

