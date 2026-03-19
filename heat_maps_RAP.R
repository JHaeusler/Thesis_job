# ============================================================================== #
# SCRIPT: MAPAS DE CALOR PARA LA SUPERFICIE DE RIESGO (WRT)
# Objetivo: Contraste topológico entre Proveedor Excelente vs. Malo
# ============================================================================== #

cat("\nGenerando Mapas de Calor del Riesgo Ponderado...\n")

# 1. Selección del Escenario a graficar (Ej: AQL = 0.05, LTPD = 0.10)
# Asumimos que el Master Script ya dejó 'Esce' y 'alpha_beta_params' en memoria
esce_plot <- 3 # Escoja el escenario que desea analizar en el Capítulo 5
N_lote <- Esce[esce_plot, 1]
AQL_val <- Esce[esce_plot, 4]
LTPD_val <- Esce[esce_plot, 5]
alpha_val <- Esce[esce_plot, 2]
beta_val <- Esce[esce_plot, 3]

# 2. Rango de exploración de la malla (Puede ajustarlo según necesidad)
n_seq <- 1:120
c_seq <- 0:8

# 3. Inicializar matrices de riesgo
Z_excelente <- matrix(NA, nrow = length(n_seq), ncol = length(c_seq))
Z_malo      <- matrix(NA, nrow = length(n_seq), ncol = length(c_seq))

# 4. Extraer parámetros (a, b) de los dos proveedores a contrastar
# Proveedor 1: Excelente | Proveedor 4: Malo (según su lista de Cook)
a_exc <- alpha_beta_params$alpha_b[1]; b_exc <- alpha_beta_params$beta_b[1]
a_mal <- alpha_beta_params$alpha_b[4]; b_mal <- alpha_beta_params$beta_b[4]

# Calcular los límites dinámicos (des_WRT) para trazar la frontera
mass_exc <- calc_prob_mass(a_exc, b_exc, AQL_val, LTPD_val)
des_WRT_exc <- alpha_val * mass_exc["P_Good"] + beta_val * mass_exc["P_Bad"]

mass_mal <- calc_prob_mass(a_mal, b_mal, AQL_val, LTPD_val)
des_WRT_mal <- alpha_val * mass_mal["P_Good"] + beta_val * mass_mal["P_Bad"]

# 5. Llenado de la matriz calculando integrales reales (Esto tomará unos segundos)
cat("Calculando integrales numéricas para la malla espacial...\n")
for (i in seq_along(n_seq)) {
  for (j in seq_along(c_seq)) {
    n_ <- n_seq[i]
    c_ <- c_seq[j]
    
    if (c_ < n_) {
      # Riesgo Proveedor Excelente
      r_exc <- calc_wr(N_lote, n_, c_, a_exc, b_exc, AQL_val, LTPD_val)
      Z_excelente[i, j] <- r_exc["WRT_val"]
      
      # Riesgo Proveedor Malo
      r_mal <- calc_wr(N_lote, n_, c_, a_mal, b_mal, AQL_val, LTPD_val)
      Z_malo[i, j] <- r_mal["WRT_val"]
    }
  }
}

# 6. Renderizado del Gráfico
file_name <- "Heatmap_Contraste_Riesgos.png"
png(file_name, width = 1200, height = 600, res = 120)

par(mfrow = c(1, 2), mar = c(5, 4, 4, 2) + 0.1)
paleta_calor <- colorRampPalette(c("forestgreen", "khaki1", "firebrick"))(50)
z_max <- max(c(Z_excelente, Z_malo), na.rm = TRUE) # Misma escala de color para ambos

# --- PANEL 1: PROVEEDOR EXCELENTE ---
image(x = n_seq, y = c_seq, z = Z_excelente, col = paleta_calor, zlim = c(0, z_max),
      xlab = "Tamaño de Muestra (n)", ylab = "Número de Aceptación (c)",
      main = "Superficie de Riesgo WRT - Proveedor Excelente")
contour(x = n_seq, y = c_seq, z = Z_excelente, add = TRUE, col = rgb(1,1,1,0.5), nlevels = 8)
contour(x = n_seq, y = c_seq, z = Z_excelente, add = TRUE, levels = des_WRT_exc, 
        col = "black", lwd = 4, lty = 1, drawlabels = FALSE)
legend("topleft", legend = "Frontera Óptima", col = "black", lwd = 4, lty = 1, bg="white", cex=0.8)

# --- PANEL 2: PROVEEDOR MALO ---
image(x = n_seq, y = c_seq, z = Z_malo, col = paleta_calor, zlim = c(0, z_max),
      xlab = "Tamaño de Muestra (n)", ylab = "Número de Aceptación (c)",
      main = "Superficie de Riesgo WRT - Proveedor Malo")
contour(x = n_seq, y = c_seq, z = Z_malo, add = TRUE, col = rgb(1,1,1,0.5), nlevels = 8)
contour(x = n_seq, y = c_seq, z = Z_malo, add = TRUE, levels = des_WRT_mal, 
        col = "black", lwd = 4, lty = 1, drawlabels = FALSE)
legend("topleft", legend = "Frontera Óptima", col = "black", lwd = 4, lty = 1, bg="white", cex=0.8)

dev.off()
cat("¡Mapa de calor generado exitosamente!\n")