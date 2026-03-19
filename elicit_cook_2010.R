# ============================================================================== #
# SCRIPT 2: ELICITACIÓN DE EXPERTOS Y CURVAS CO MULTIESCENARIO (R BASE)
# Metodología: Método de Cook (2010) y Malla de Escenarios
# Autor: Juan Sebastián Haeusler
# ============================================================================== #

# --- 1. Carga de Paquetes ---
if(!require(stats)) install.packages("stats")
if(!require(AcceptanceSampling)) install.packages("AcceptanceSampling")

# ============================================================================== #
# --- 2. Función del Método de Cook ---
# ============================================================================== #
find_beta <- function(x1, p1, x2, p2) {
  loss_function <- function(params) {
    a <- params[1]; b <- params[2]
    if(a <= 0 || b <= 0) return(Inf)
    error <- (pbeta(x1, a, b) - p1)^2 + (pbeta(x2, a, b) - p2)^2
    return(error)
  }
  opt <- optim(par = c(1, 1), fn = loss_function, method = "Nelder-Mead")
  if(opt$convergence == 0) return(list(alpha = opt$par[1], beta = opt$par[2]))
  else return(list(alpha = NA, beta = NA))
}

# ============================================================================== #
# --- 3. Definición de la Malla de Escenarios ---
# ============================================================================== #
# N_vec     <- 1000
# alpha_vec <- 0.01
# beta_vec  <- 0.05
# AQL_vec   <- c(0.01, 0.02, 0.05)
# LTPD_vec  <- c(0.08, 0.10, 0.15, 0.20)

Esce <- expand.grid(N = N_vec, alpha = alpha_vec, beta = beta_vec, 
                    AQL = AQL_vec, LTPD = LTPD_vec)
# Filtrar combinaciones lógicas (AQL estrictamente menor que LTPD)
Esce <- Esce[Esce$AQL < Esce$LTPD, ]
rownames(Esce) <- NULL

# Perfiles del Proveedor (p1 = PA en AQL, p2 = Riesgo Consumidor en LTPD)
p1_val <- c(0.95, 0.75, 0.50, 0.25, 0.10)
p2_val <- c(1e-4, 0.05, 0.10, 0.25, 0.40)
nombres_perfiles <- c("Excelente", "Bueno", "Regular", "Malo", "Muy Malo")

colores <- c(rgb(0.1, 0.6, 0.2, 0.4), rgb(0.2, 0.5, 0.8, 0.4), 
             rgb(0.8, 0.6, 0.1, 0.4), rgb(0.8, 0.4, 0.1, 0.4), rgb(0.8, 0.1, 0.1, 0.4))
bordes  <- c("forestgreen", "dodgerblue", "darkgoldenrod", "darkorange", "firebrick")

# Lista maestra para guardar los parámetros elicitados de todos los escenarios
lista_parametros_escenarios <- list()

# ============================================================================== #
# --- 4. Ejecución y Generación del Atlas en PDF ---
# ============================================================================== #
cat("Generando Atlas de Escenarios en PDF...\n")

# Se crea un PDF en su directorio de trabajo para almacenar todas las gráficas
pdf("Atlas_Escenarios_Calidad.pdf", width = 8, height = 10)

for (k in 1:nrow(Esce)) {
  
  N_cur    <- Esce$N[k]
  AQL_cur  <- Esce$AQL[k]
  LTPD_cur <- Esce$LTPD[k]
  alp_cur  <- Esce$alpha[k]
  bet_cur  <- Esce$beta[k]
  
  # A. Hallar el Plan Clásico
  plan_clasic <- find.plan(PRP = c(AQL_cur, 1 - alp_cur), 
                           CRP = c(LTPD_cur, bet_cur), 
                           N = N_cur, type = "hypergeom")
  n_c <- plan_clasic$n; c_c <- plan_clasic$c
  
  # B. Elicitación de Parámetros (Cook) para los 5 perfiles
  param_df <- data.frame(Proveedor = nombres_perfiles, alpha_b = NA, beta_b = NA)
  
  for(i in 1:5) {
    res <- find_beta(x1 = AQL_cur, p1 = p1_val[i], x2 = LTPD_cur, p2 = 1 - p2_val[i])
    param_df$alpha_b[i] <- res$alpha
    param_df$beta_b[i]  <- res$beta
  }
  
  lista_parametros_escenarios[[paste0("Escenario_", k)]] <- param_df
  
  # C. Configuración del Panel Dual (2 filas, 1 columna)
  par(mfrow = c(2, 1), mar = c(4, 4, 3, 2), oma = c(0, 0, 3, 0))
  
  # --- PANEL SUPERIOR: CURVA CO CLÁSICA ---
  p_seq <- seq(0, 0.3, length.out = 500)
  Pa_seq <- phyper(c_c, N_cur * p_seq, N_cur * (1 - p_seq), n_c)
  
  plot(p_seq, Pa_seq, type = "l", lwd = 2, col = "black", bty = "l",
       xlab = "Proporción de Defectuosos (p)", ylab = "Probabilidad de Aceptación (Pa)",
       main = sprintf("Plan Clásico (n = %d, c = %d)", n_c, c_c))
  
  # Zonas de riesgo y líneas punteadas
  polygon(c(0, AQL_cur, AQL_cur, 0), c(0, 0, 1, 1), col = rgb(0.1,0.6,0.2,0.1), border = NA)
  polygon(c(LTPD_cur, 0.3, 0.3, LTPD_cur), c(0, 0, 1, 1), col = rgb(0.8,0.1,0.1,0.1), border = NA)
  
  abline(v = AQL_cur, col = "forestgreen", lty = 2, lwd = 1.5)
  abline(v = LTPD_cur, col = "firebrick", lty = 2, lwd = 1.5)
  points(AQL_cur, 1 - alp_cur, pch = 19, col = "forestgreen")
  points(LTPD_cur, bet_cur, pch = 19, col = "firebrick")
  
  # --- PANEL INFERIOR: DENSIDADES ELICITADAS ---
  plot(NULL, xlim = c(0, 0.3), ylim = c(0, 35), bty = "l",
       xlab = "Proporción de Defectuosos (p)", ylab = "Densidad A Priori",
       main = "Perfiles de Calidad del Proveedor (Cook)")
  
  for(i in 1:5) {
    y_seq <- dbeta(p_seq, shape1 = param_df$alpha_b[i], shape2 = param_df$beta_b[i])
    y_seq[!is.finite(y_seq)] <- 0 
    polygon(c(p_seq, rev(p_seq)), c(y_seq, rep(0, length(y_seq))), 
            col = colores[i], border = bordes[i], lwd = 2)
  }
  
  abline(v = AQL_cur, col = "forestgreen", lty = 2, lwd = 1.5)
  abline(v = LTPD_cur, col = "firebrick", lty = 2, lwd = 1.5)
  
  legend("topright", legend = nombres_perfiles, fill = colores, border = bordes, 
         bty = "n", cex = 0.8, title = "Proveedor")
  
  # Título global de la página
  mtext(sprintf("Escenario %d: N=%d | AQL=%.2f | LTPD=%.2f | \u03B1=%.2f | \u03B2=%.2f", 
                k, N_cur, AQL_cur, LTPD_cur, alp_cur, bet_cur), 
        outer = TRUE, cex = 1.2, font = 2)
}

dev.off() # Cierra el PDF
par(mfrow = c(1, 1)) # Restaura la ventana gráfica de RStudio

cat("\n¡Proceso finalizado! Revise el archivo 'Atlas_Escenarios_Calidad.pdf' en su directorio de trabajo.\n")