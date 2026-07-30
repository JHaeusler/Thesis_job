# ============================================================================== #
# --- SCRIPT: PANEL MATRICIAL (FILAS=AQL, COLUMNAS=LTPD) EN PNG ---
# ============================================================================== #

# 1. Paletas de colores y nombres
colores <- c(rgb(0.1, 0.6, 0.2, 0.4), rgb(0.2, 0.5, 0.8, 0.4), 
             rgb(0.8, 0.6, 0.1, 0.4), rgb(0.8, 0.4, 0.1, 0.4), 
             rgb(0.8, 0.1, 0.1, 0.4), rgb(238, 18, 137, alpha =102, maxColorValue = 255),
             rgb(155, 205, 155, alpha =102, maxColorValue = 255))

bordes  <- c("forestgreen", "dodgerblue", "darkgoldenrod", "darkorange", 
             "firebrick", "deeppink2", "darkseagreen3")

nombres_perfiles <- c("Excelente", "Bueno", "Regular", "Malo", "Muy Malo", "ZI_90", "ZI_98")

cat("\nGenerando Panel Matricial de Densidades en PNG...\n")

# --- ¡LÍNEAS DESCOMENTADAS! ---
png("imagenes/Panel_dens_cook.png", width = 1400,
    height = 800, res = 100)

# 2. Configuración del panel: 3 filas (AQL) x 4 columnas (LTPD)
par(mfrow = c(3, 4), mar = c(4, 4, 3, 1), oma = c(4, 0, 3, 0))

p_seq <- seq(0, 0.3, length.out = 500)
lim_Y <- 35 

# Ordenamos para asegurar que se llene la matriz por filas correctamente
diccionario_cook <- diccionario_cook[order(diccionario_cook$AQL, diccionario_cook$LTPD), ]

# 3. Iteramos DIRECTAMENTE sobre diccionario_cook
for(k in 1:nrow(diccionario_cook)) {# k <- 1 + k
  
  aql_cur  <- diccionario_cook$AQL[k]
  ltpd_cur <- diccionario_cook$LTPD[k]
  
  # A. Lienzo vacío para las densidades
  plot(NULL, xlim = c(0, 0.3), ylim = c(0, lim_Y), bty = "l",
       cex.main = 1.5,
       cex.lab = 1.5,
       cex.axis = 1.5,
       xlab = "Proporción de Defectuosos (p)", ylab = "Densidad A Priori",
       main = sprintf("AQL = %.2f | LTPD = %.2f", aql_cur, ltpd_cur))
  
  # B. Iteramos sobre los 7 perfiles para dibujar sus polígonos
  for(i in 1:length(nombres_perfiles)) {
    nombre_perf <- nombres_perfiles[i]
    
    # Extraemos 'a' y 'b' asumiendo que el dataframe está en formato "ancho"
    a_val <- diccionario_cook[[paste0(nombre_perf, "_a")]][k]
    b_val <- diccionario_cook[[paste0(nombre_perf, "_b")]][k]
    
    if (length(a_val) > 0 && length(b_val) > 0 && !is.na(a_val) && !is.na(b_val)) {
      
      y_seq <- dbeta(p_seq, shape1 = a_val, shape2 = b_val)
      
      y_seq[!is.finite(y_seq)] <- lim_Y
      y_seq[y_seq > lim_Y] <- lim_Y
      
      polygon(c(p_seq, rev(p_seq)), c(y_seq, rep(0, length(y_seq))), 
              col = colores[i], border = bordes[i], lwd = 2)
    }
  }
  
  # C. Líneas verticales contractuales
  abline(v = aql_cur, col = "forestgreen", lty = 2, lwd = 1.5)
  abline(v = ltpd_cur, col = "firebrick", lty = 2, lwd = 1.5)
}

# 4. Leyenda Global y Título
mtext("Perfiles de Calidad: Densidades A Priori (Filas = AQL, Columnas = LTPD)", 
      outer = TRUE, side = 3, line = 0.5, cex = 1.6, font = 2)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend("bottom", legend = nombres_perfiles, cex = 1.5,
       fill = colores, border = bordes, bty = "n", horiz = TRUE, inset = c(0, 0.02))

# --- ¡LÍNEA DESCOMENTADA! ---
dev.off()

cat("\n¡Panel de Densidades guardado exitosamente como .png!\n")