# ============================================================================== #
# --- 3. PARTE B: GENERACIÓN DEL ATLAS EN PDF (BASADO EN ESCE_GLOBAL) ---
# ============================================================================== #

colores <- c(rgb(0.1, 0.6, 0.2, 0.4), rgb(0.2, 0.5, 0.8, 0.4), 
             rgb(0.8, 0.6, 0.1, 0.4), rgb(0.8, 0.4, 0.1, 0.4), 
             rgb(0.8, 0.1, 0.1, 0.4), rgb(238, 18, 137, alpha =102, maxColorValue = 255),
             rgb(155, 205, 155, alpha =102, maxColorValue = 255))
bordes  <- c("forestgreen", "dodgerblue", "darkgoldenrod", "darkorange", "firebrick", "deeppink2", "darkseagreen3")

cat("\nGenerando Atlas de Escenarios en PDF...\n")
pdf("Atlas_Escenarios_Calidad.pdf", width = 8, height = 10)

# Iteramos DIRECTAMENTE sobre la tabla maestra que ya tiene los planes calculados
for (k in 1:nrow(Esce_Global)) {
  
  N_cur    <- Esce_Global$N[k]
  AQL_cur  <- Esce_Global$AQL[k]
  LTPD_cur <- Esce_Global$LTPD[k]
  alp_cur  <- Esce_Global$alpha[k]
  bet_cur  <- Esce_Global$beta[k]
  
  # MAGIA AQUÍ: Extraemos el plan exacto que se usó para la simulación
  n_c <- Esce_Global$n_clasico[k]
  c_c <- Esce_Global$c_clasico[k]
  
  # Buscar los parámetros en el diccionario (usando round() para evitar problemas)
  params_cook <- diccionario_cook %>% 
    filter(round(AQL,4) == round(AQL_cur,4), round(LTPD,4) == round(LTPD_cur,4))
  
  # Gráficas
  par(mfrow = c(2, 1), mar = c(4, 4, 3, 2), oma = c(0, 0, 3, 0))
  
  # PANEL SUPERIOR (Curva Clásica real del escenario)
  p_seq <- seq(0, 0.3, length.out = 500)
  Pa_seq <- phyper(c_c, round(N_cur * p_seq), round(N_cur * (1 - p_seq)), n_c)
  
  plot(p_seq, Pa_seq, type = "l", lwd = 2, col = "black", bty = "l",
       xlab = "Proporción de Defectuosos (p)", ylab = "Probabilidad de Aceptación (Pa)",
       main = sprintf("Plan Clásico (n = %d, c = %d)", n_c, c_c))
  polygon(c(0, AQL_cur, AQL_cur, 0), c(0, 0, 1, 1), col = rgb(0.1,0.6,0.2,0.1), border = NA)
  polygon(c(LTPD_cur, 0.3, 0.3, LTPD_cur), c(0, 0, 1, 1), col = rgb(0.8,0.1,0.1,0.1), border = NA)
  abline(v = AQL_cur, col = "forestgreen", lty = 2, lwd = 1.5)
  abline(v = LTPD_cur, col = "firebrick", lty = 2, lwd = 1.5)
  points(AQL_cur, 1 - alp_cur, pch = 19, col = "forestgreen")
  points(LTPD_cur, bet_cur, pch = 19, col = "firebrick")
  
  # PANEL INFERIOR (Densidades Elicitadas de Cook)
  lim_Y <- 35
  plot(NULL, xlim = c(0, 0.3), ylim = c(0, lim_Y), bty = "l",
       xlab = "Proporción de Defectuosos (p)", ylab = "Densidad A Priori",
       main = "Perfiles de Calidad (Método de Cook)")
  
  for(i in 1:length(nombres_perfiles)) {
    nombre_perf <- nombres_perfiles[i]
    
    a_val <- params_cook[[paste0(nombre_perf, "_a")]]
    b_val <- params_cook[[paste0(nombre_perf, "_b")]]
    
    if (length(a_val) > 0 && length(b_val) > 0 && !is.na(a_val) && !is.na(b_val)) {
      y_seq <- dbeta(p_seq, shape1 = a_val, shape2 = b_val)
      y_seq[!is.finite(y_seq)] <- lim_Y
      y_seq[y_seq > lim_Y] <- lim_Y
      
      polygon(c(p_seq, rev(p_seq)), c(y_seq, rep(0, length(y_seq))), 
              col = colores[i], border = bordes[i], lwd = 2)
    }
  }
  
  abline(v = AQL_cur, col = "forestgreen", lty = 2, lwd = 1.5)
  abline(v = LTPD_cur, col = "firebrick", lty = 2, lwd = 1.5)
  legend("topright", legend = nombres_perfiles, fill = colores, border = bordes, bty = "n", cex = 0.8)
  
  mtext(sprintf("Escenario %d: N=%d | AQL=%.2f | LTPD=%.2f | \u03B1=%.2f | \u03B2=%.2f", 
                k, N_cur, AQL_cur, LTPD_cur, alp_cur, bet_cur), outer = TRUE, cex = 1.2, font = 2)
}

dev.off()
par(mfrow = c(1, 1))
cat("\n¡Proceso finalizado! Revise 'Atlas_Escenarios_Calidad.pdf'.\n")