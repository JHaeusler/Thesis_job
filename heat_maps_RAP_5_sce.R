# ============================================================================== #
# SCRIPT: MAPAS DE CALOR PARA LA SUPERFICIE DE RIESGO (5 PERFILES)
# Objetivo: Evolución topológica del WRT (Excelente -> Muy Malo)
# ============================================================================== #

cat("\nIniciando cálculo de mallas espaciales para los 5 perfiles...\n")

# 1. Selección del Escenario Base a graficar (Ej: Escenario 3: AQL=0.05, LTPD=0.10)
esce_plot <- 3 
N_lote <- Esce[esce_plot, 1]
AQL_val <- Esce[esce_plot, 4]
LTPD_val <- Esce[esce_plot, 5]
alpha_val <- Esce[esce_plot, 2]
beta_val <- Esce[esce_plot, 3]

# 2. Rango de la malla (Ajustado para visualizar bien la transición)
n_seq <- 1:120
c_seq <- 0:8

# 3. Lista para almacenar las matrices de riesgo de los 5 proveedores
Z_list <- list()
des_WRT_list <- numeric(5)
nombres_5 <- c("Excelente", "Bueno", "Regular", "Malo", "Muy Malo")

# 4. Cálculo de las 5 matrices (Esto tomará un par de minutos)
for (k in 1:5) {
  cat(sprintf(">> Procesando Proveedor: %s (%d de 5)...\n", nombres_5[k], k))
  
  a_val <- alpha_beta_params$alpha_b[k]
  b_val <- alpha_beta_params$beta_b[k]
  
  # Límite dinámico (Frontera estocástica)
  mass <- calc_prob_mass(a_val, b_val, AQL_val, LTPD_val)
  des_WRT_list[k] <- alpha_val * mass["P_Good"] + beta_val * mass["P_Bad"]
  
  # Matriz temporal para este proveedor
  Z_temp <- matrix(NA, nrow = length(n_seq), ncol = length(c_seq))
  
  for (i in seq_along(n_seq)) {
    for (j in seq_along(c_seq)) {
      n_ <- n_seq[i]
      c_ <- c_seq[j]
      
      if (c_ < n_) {
        r_val <- calc_wr(N_lote, n_, c_, a_val, b_val, AQL_val, LTPD_val)
        Z_temp[i, j] <- r_val["WRT_val"]
      }
    }
  }
  Z_list[[k]] <- Z_temp
}

# 5. Renderizado del Panel Completo (2 filas x 3 columnas)
cat("\nRenderizando el gráfico consolidado...\n")
file_name <- "Heatmap_5_Perfiles.png"
png(file_name, width = 1500, height = 900, res = 120)

# Configuramos un layout de 2 filas y 3 columnas
par(mfrow = c(2, 3), mar = c(4, 4, 3, 2) + 0.1, oma = c(0, 0, 3, 0))

# Encontrar el valor máximo global para compartir la misma escala de colores
z_max <- max(sapply(Z_list, max, na.rm = TRUE))
paleta_calor <- colorRampPalette(c("forestgreen", "khaki1", "firebrick"))(50)

# Bucle para graficar los 5 paneles
for (k in 1:5) {
  image(x = n_seq, y = c_seq, z = Z_list[[k]], 
        col = paleta_calor, zlim = c(0, z_max),
        xlab = "Tamaño de Muestra (n)", ylab = "Número de Aceptación (c)",
        main = paste("Proveedor:", nombres_5[k]))
  
  # Líneas de contorno generales
  contour(x = n_seq, y = c_seq, z = Z_list[[k]], 
          add = TRUE, col = rgb(1,1,1,0.5), nlevels = 8)
  
  # FRONTERA ÓPTIMA (Línea negra gruesa)
  contour(x = n_seq, y = c_seq, z = Z_list[[k]], 
          add = TRUE, levels = des_WRT_list[k], 
          col = "black", lwd = 4, lty = 1, drawlabels = FALSE)
}

# Sexto panel: Leyenda y anotaciones
plot.new()
legend("center", legend = c("Zona Segura (Riesgo Bajo)", 
                            "Zona Crítica (Riesgo Alto)", 
                            "Frontera de Factibilidad (des_WRT)"),
       fill = c("forestgreen", "firebrick", NA),
       border = c("black", "black", NA),
       lty = c(NA, NA, 1), lwd = c(NA, NA, 4), col = c(NA, NA, "black"),
       cex = 1.2, bty = "n")

mtext(paste("Evolución de la Superficie de Riesgo WRT - Escenario AQL:", AQL_val, "LTPD:", LTPD_val), 
      outer = TRUE, cex = 1.5, font = 2)

dev.off()
cat("¡Gráfico 'Heatmap_5_Perfiles.png' generado exitosamente en su directorio!\n")