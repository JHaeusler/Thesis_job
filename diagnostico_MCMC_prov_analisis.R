# ============================================================================== #
# SCRIPT REFINADO: FIGURAS PARA CUERPO DE TESIS (ESTILO ATLAS FINAL)
# Ajuste de Títulos (cex.main) y Nomenclatura Académica (Prov 3A/3B)
# ============================================================================== #

library(coda)

# 1. Mapeo de nombres para títulos académicos
# Relacionamos el nombre técnico del objeto con el nombre que irá en la Tesis
nombres_tesis <- c(
  "Prov1"          = "Proveedor 1",
  "Prov2"          = "Proveedor 2",
  "Prov3_Esce_66"  = "Proveedor 3A",
  "Prov3_Esce_13"  = "Proveedor 3B"
)

escenarios_clave <- names(nombres_tesis)

# 2. Configuración estética
col_priors <- c("Gamma" = rgb(0.5, 0.5, 0.5, 0.4), "Jeffreys" = rgb(0.8, 0.2, 0.2, 0.4), "TTB" = rgb(0.1, 0.5, 0.8, 0.4))
brd_priors <- c("Gamma" = "gray30", "Jeffreys" = "firebrick", "TTB" = "dodgerblue")

for (esce in escenarios_clave) { # esce <- 3 + esce
  if (!esce %in% names(Cadenas_MCMC_List)) next
  
  cadenas <- Cadenas_MCMC_List[[esce]]
  nombre_grafico <- nombres_tesis[esce] # Nombre amigable para el título
  
  # --- FIGURA 1: DIAGNÓSTICO 4 PANELES (Parámetros a y b) ---
  for (par_idx in 1:2) { # par_idx <- 1 + par_idx
    p_name <- ifelse(par_idx == 1, "a", "b")
    chain <- cadenas$TTB[, par_idx]
    met   <- calcular_metricas(chain) 
    
    png(paste0("imagenes/MCMC_Diag_", esce, "_", p_name, ".png"), width = 1000, height = 800, res = 120)
    par(mfrow = c(2, 2), mar = c(4, 5, 4, 1), oma = c(0, 0, 6, 0))
    
    # (a) Histograma y Densidad - AJUSTE cex.main = 1.5
    hist(chain, prob = TRUE, col = "gray90", border = "gray",
         main = "(a) Histograma y Densidad", xlab = p_name, ylab = "Densidad",
         cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
    lines(density(chain), col = brd_priors["TTB"], lwd = 2)
    
    # (b) Traza con límites
    plot(chain, type = "l", col = "gray40", main = "(b) Traza de la Cadena",
         xlab = "Iteración t", ylab = p_name, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
    abline(h = range(chain), col = "firebrick", lty = 2)
    
    # (c) Media Ergódica
    erg_mean <- cumsum(chain) / seq_along(chain)
    N_chain  <- length(chain)
    mcse     <- sd(chain) / sqrt(met$ESS)
    
    plot(erg_mean, type = "l", col = "firebrick", lwd = 2, main = "(c) Media Ergódica",
         xlab = "Iteraciones t", ylab = "Promedio Acumulado", 
         ylim = range(c(erg_mean, erg_mean[N_chain] + 1.96*mcse, erg_mean[N_chain] - 1.96*mcse)),
         cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
    abline(h = erg_mean[N_chain], col = "black", lty = 1, lwd = 1.5)
    abline(h = c(erg_mean[N_chain] + 1.96*mcse, erg_mean[N_chain] - 1.96*mcse), col = "blue", lty = 2)
    
    # (d) ACF
    acf(chain, main = "", ylab = "ACF", cex.lab = 1.5, cex.axis = 1.5)
    title(main = "(d) Autocorrelación (ACF)", cex.main = 1.5)
    
    # Títulos Superiores con Nombre de Tesis
    title_main <- paste("Diagnóstico MCMC:", nombre_grafico, "- Parámetro", p_name)
    diag_stats <- sprintf("ESS: %.1f | Acceptance Rate: %.1f%% | R-hat: %.3f", 
                          met$ESS, met$AR * 100, ifelse(is.na(met$Rhat), 1.000, met$Rhat))
    
    mtext(title_main, outer = TRUE, cex = 1.5, font = 2, line = 2.5)
    mtext(diag_stats, outer = TRUE, cex = 1.2, font = 1, line = 0.5, col = "darkblue")
    
    dev.off()
  }
  
  # --- FIGURA 2: COMPARATIVA DE POSTERIORS (2 PANELES) ---
  png(paste0("imagenes/Posterior_Comparativa_", esce, ".png"), width = 1100, height = 550, res = 120)
  par(mfrow = c(1, 2), oma = c(0, 0, 4, 0), mar = c(4, 5, 4, 2))
  
  for (j in 1:2) {
    p_n <- ifelse(j == 1, "a", "b")
    dg <- density(cadenas$Gamma[,j]); dj <- density(cadenas$Jeffreys[,j]); db <- density(cadenas$TTB[,j])
    
    plot(NULL, xlim = range(c(dg$x, dj$x, db$x)), ylim = c(0, max(c(dg$y, dj$y, db$y)) * 1.2),
         main = paste("Posteriors del Parámetro", p_n), xlab = p_n, ylab = "Densidad", 
         bty = "l", cex.main = 1.5, cex.lab = 1.2)
    
    polygon(dg, col = col_priors["Gamma"], border = brd_priors["Gamma"], lwd = 2)
    polygon(dj, col = col_priors["Jeffreys"], border = brd_priors["Jeffreys"], lwd = 2)
    polygon(db, col = col_priors["TTB"], border = brd_priors["TTB"], lwd = 2)
    
    if(j == 1) legend("topright", legend = names(col_priors), fill = col_priors, border = brd_priors, bty = "n", cex = 1)
  }
  mtext(paste("Análisis de Sensibilidad de Priors -", nombre_grafico), outer = TRUE, cex = 1.6, font = 2)
  dev.off()
}

tabla_latex <- data.frame()

for (esce in escenarios_clave) {
  if (!esce %in% names(Cadenas_MCMC_List)) next
  
  chain_a <- Cadenas_MCMC_List[[esce]]$TTB[, 1]
  chain_b <- Cadenas_MCMC_List[[esce]]$TTB[, 2]
  
  met_a <- calcular_metricas(chain_a)
  met_b <- calcular_metricas(chain_b)
  
  # Extraemos métricas
  a_hat <- mean(chain_a)
  b_hat <- mean(chain_b)
  ess_min <- min(met_a$ESS, met_b$ESS)
  rhat_max <- max(met_a$Rhat, met_b$Rhat, na.rm = TRUE)
  
  # Calculamos el AR mínimo en porcentaje (para ser conservadores)
  ar_min_pct <- min(met_a$AR, met_b$AR) * 100
  
  fila <- data.frame(
    Proveedor = nombres_tesis[esce],
    a_hat     = round(a_hat, 2),
    b_hat     = round(b_hat, 2),
    ESS_min   = round(ess_min, 0),
    AR_pct    = round(ar_min_pct, 1), # Nueva columna
    Rhat_max  = round(rhat_max, 3)
  )
  
  tabla_latex <- rbind(tabla_latex, fila)
}

cat("\n--- DATOS PARA LA TABLA DE LATEX ---\n")
print(tabla_latex, row.names = FALSE)



