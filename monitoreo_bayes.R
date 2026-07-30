# ============================================================================== #
# SCRIPT VISUAL Y TABULAR: DIAGNÓSTICOS MCMC (ESTILO TTB)
# Genera el Atlas en PDF y exporta la tabla resumen sin recalcular
# ============================================================================== #

if(!require(coda)) {
  install.packages("coda", dependencies = TRUE); library(coda)
} else { library(coda) }

if(!exists("Cadenas_MCMC_List")) {
  stop("Error: No se encontró 'Cadenas_MCMC_List' en la memoria.")
}

cat("\n=================================================================\n")
cat(" GENERANDO ATLAS MCMC Y TABLA RESUMEN SIMULTÁNEAMENTE...\n")
cat("=================================================================\n")

# --- 1. FUNCIÓN MATEMÁTICA PURA (No grafica, solo calcula) ---
calcular_metricas <- function(chain) {
  N_chain <- length(chain)
  acc_rate <- sum(diff(chain) != 0) / (N_chain - 1)
  ess_val <- effectiveSize(mcmc(chain))
  
  half <- floor(N_chain / 2)
  c1 <- mcmc(chain[1:half]); c2 <- mcmc(chain[(half + 1):(2 * half)])
  rhat_val <- tryCatch(
    as.numeric(gelman.diag(mcmc.list(c1, c2), autoburnin=FALSE)$psrf[1]), 
    error = function(e) NA
  )
  return(list(ESS = ess_val, AR = acc_rate, Rhat = rhat_val))
}

# --- 2. FUNCIONES GRÁFICAS MODIFICADAS (Reciben las métricas ya calculadas) ---
graficar_4_paneles <- function(chain, param_name, title_text, metricas) {
  N_chain <- length(chain)
  old_par <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 4, 0))
  
  # (a) Histograma
  hist(chain, prob = TRUE, col = rgb(0.8, 0.8, 0.8, 0.5), border = "gray",
       main = "(a) Histograma y Densidad", xlab = param_name, ylab = "Densidad")
  lines(density(chain), col = "dodgerblue", lwd = 2)
  
  # (b) Traza
  min_val <- min(chain); max_val <- max(chain)
  plot(chain, type = "l", col = "gray40", main = paste("(b) Traza de la Cadena (N =", N_chain, ")"),
       xlab = "Iteración t", ylab = param_name)
  abline(h = min_val, col = "firebrick", lty = 2); abline(h = max_val, col = "firebrick", lty = 2) 
  
  # (c) Media Ergódica
  ergodic_mean <- cumsum(chain) / seq_along(chain)
  final_mean <- ergodic_mean[N_chain]
  mcse <- sd(chain) / sqrt(metricas$ESS)
  upper_bound <- final_mean + 1.96 * mcse; lower_bound <- final_mean - 1.96 * mcse
  y_limits <- range(c(ergodic_mean, lower_bound, upper_bound))
  
  plot(ergodic_mean, type = "l", col = "firebrick", lwd = 2, main = "(c) Media Ergódica",
       xlab = "Iteraciones t", ylab = "Promedio Acumulado", ylim = y_limits)
  abline(h = final_mean, col = "black", lty = 1, lwd = 1.5)
  abline(h = upper_bound, col = "blue", lty = 2, lwd = 1.5); abline(h = lower_bound, col = "blue", lty = 2, lwd = 1.5) 
  text(x = N_chain * 0.8, y = final_mean, labels = sprintf("Media: %.4f", final_mean), pos = 3, cex = 0.9, font = 2)
  
  # (d) ACF
  acf(chain, main = "(d) Función de Autocorrelación (ACF)", ylab = "ACF", xlab = "Lag")
  
  # Textos Superiores RECICLANDO los valores
  diag_text <- sprintf("ESS: %.1f | Acceptance Rate: %.1f%% | Split R-hat: %.3f", 
                       metricas$ESS, metricas$AR * 100, ifelse(is.na(metricas$Rhat), 1.000, metricas$Rhat))
  
  mtext(title_text, outer = TRUE, cex = 1.2, font = 2, line = 1.8)
  mtext(diag_text, outer = TRUE, cex = 1, font = 1, line = 0.2, col = "darkblue")
  par(old_par)
}

graficar_posteriores <- function(res_g, res_j, res_b, title_text) {
  old_par <- par(mfrow = c(1, 2), oma = c(0, 0, 3, 0), mar = c(4, 4, 3, 1))
  # Param a
  da_g <- density(res_g[,1]); da_j <- density(res_j[,1]); da_b <- density(res_b[,1])
  plot(NULL, xlim = range(c(da_g$x, da_j$x, da_b$x)), ylim = c(0, max(c(da_g$y, da_j$y, da_b$y)) * 1.1),
       main = "Densidad Posterior Parámetro 'a'", xlab = "a", ylab = "Densidad", bty = "l")
  polygon(da_g, col = rgb(0.5, 0.5, 0.5, 0.4), border = "gray30", lwd = 2)
  polygon(da_j, col = rgb(0.8, 0.2, 0.2, 0.4), border = "firebrick", lwd = 2)
  polygon(da_b, col = rgb(0.1, 0.5, 0.8, 0.4), border = "dodgerblue", lwd = 2)
  legend("topright", legend = c("Gamma", "Jeffreys", "TTB"),
         fill = c(rgb(0.5,0.5,0.5,0.4), rgb(0.8,0.2,0.2,0.4), rgb(0.1,0.5,0.8,0.4)),
         border = c("gray30", "firebrick", "dodgerblue"), bty = "n", cex = 0.8)
  
  # Param b
  db_g <- density(res_g[,2]); db_j <- density(res_j[,2]); db_b <- density(res_b[,2])
  plot(NULL, xlim = range(c(db_g$x, db_j$x, db_b$x)), ylim = c(0, max(c(db_g$y, db_j$y, db_b$y)) * 1.1),
       main = "Densidad Posterior Parámetro 'b'", xlab = "b", ylab = "Densidad", bty = "l")
  polygon(db_g, col = rgb(0.5, 0.5, 0.5, 0.4), border = "gray30", lwd = 2)
  polygon(db_j, col = rgb(0.8, 0.2, 0.2, 0.4), border = "firebrick", lwd = 2)
  polygon(db_b, col = rgb(0.1, 0.5, 0.8, 0.4), border = "dodgerblue", lwd = 2)
  mtext(title_text, outer = TRUE, cex = 1.2, font = 2)
  par(old_par)
}

# --- 3. EJECUCIÓN MAESTRA: DATOS Y PDF ---
tabla_diagnosticos <- data.frame()
nombres_fuentes <- names(Cadenas_MCMC_List)

# timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
pdf("Atlas_MCMC_Diagnosticos.pdf", width = 15, height = 10)

for (fuente in nombres_fuentes) {
  cat(sprintf("\r   Procesando e imprimiendo %s...", fuente))
  cadenas <- Cadenas_MCMC_List[[fuente]]
  
  # Calcular métricas para TODOS los priors y armar la tabla
  for (prior_name in names(cadenas)) {
    cad_a <- cadenas[[prior_name]][, 1]
    cad_b <- cadenas[[prior_name]][, 2]
    
    met_a <- calcular_metricas(cad_a)
    met_b <- calcular_metricas(cad_b)
    
    fila <- data.frame(
      Escenario = fuente, Prior = prior_name,
      ESS_a = round(met_a$ESS, 1), ESS_b = round(met_b$ESS, 1),
      AR_a_pct = round(met_a$AR * 100, 1), AR_b_pct = round(met_b$AR * 100, 1),
      Rhat_a = round(met_a$Rhat, 3), Rhat_b = round(met_b$Rhat, 3)
    )
    tabla_diagnosticos <- rbind(tabla_diagnosticos, fila)
    
    # SI es el modelo TTB, reciclamos las métricas para dibujar los 4 paneles
    if (prior_name == "TTB") {
      graficar_4_paneles(cad_a, "a", paste("Diagnósticos MCMC 'a' (TTB) -", fuente), met_a)
      graficar_4_paneles(cad_b, "b", paste("Diagnósticos MCMC 'b' (TTB) -", fuente), met_b)
    }
  }
  
  # Al final del escenario, comparamos densidades
  graficar_posteriores(cadenas$Gamma, cadenas$Jeffreys, cadenas$TTB, paste("Comparación de Priors -", fuente))
}

dev.off()
write.csv(tabla_diagnosticos, "Resumen_Diagnosticos_MCMC.csv", row.names = FALSE)

cat("\n¡Proceso finalizado con eficiencia!\n")
cat(sprintf("- PDF generado: %s\n", "Atlas_Escenarios_Calidad.pdf"))
cat("- Tabla guardada: Resumen_Diagnosticos_MCMC.csv\n")