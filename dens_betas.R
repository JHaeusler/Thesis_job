library(GA)
library(AcceptanceSampling)

# Parámetros AQL y RQL fijos
AQL <- 0.05
RQL <- 0.10

# Vectores de probabilidades para los 5 escenarios
p1 <- c(0.95, 0.75, 0.25, 0.07, 0.02)
p2 <- 1 - c(0.02, 0.07, 0.25, 0.75, 0.95)

# Función Fitness (Ahora acepta p_target1 y p_target2)
fitness_function_iterative <- function(x, p_target1, p_target2) {
  alpha <- x[1]
  beta <- x[2]
  
  if (alpha <= 0 || beta <= 0) {
    return(-Inf)
  }
  
  # La función ahora usa los valores escalares p_target1 y p_target2
  diff1 <- pbeta(AQL, shape1 = alpha, shape2 = beta) - p_target1
  diff2 <- pbeta(RQL, shape1 = alpha, shape2 = beta) - p_target2
  
  SSE <- diff1^2 + diff2^2
  
  return(-SSE) # Maximizamos el NEGATIVO del SSE
}

## Iteración de la Optimización 

# Crear una estructura para almacenar los resultados
results_matrix <- matrix(NA, nrow = 5, ncol = 4, 
                         dimnames = list(
                           paste0("Escenario ", 1:5),
                           c("alpha (a)", "beta (b)", "SSE Mínimo", "Max Fitness")
                         ))

# Bucle para cada uno de los 5 escenarios
for (i in 1:5) {
  cat(sprintf("Procesando Escenario %d: P(X<=0.05)=%g, P(X<=0.10)=%g\n", i, p1[i], p2[i]))
  
  # Ejecución del GA:
  ga_result <- ga(type = "real-valued",
                  # Usamos el argumento '...' para pasar los targets p1[i] y p2[i]
                  fitness = fitness_function_iterative,
                  p_target1 = p1[i],
                  p_target2 = p2[i],
                  
                  lower = c(1e-4, 1e-4),
                  upper = c(100, 10000),
                  popSize = 100,
                  maxiter = 500,
                  run = 100,
                  seed = 42 + i, # Semilla diferente para cada corrida
                  optim = TRUE,
                  monitor = FALSE)
  
  # Almacenar resultados
  results_matrix[i, "alpha (a)"] <- ga_result@solution[1, 1]
  results_matrix[i, "beta (b)"] <- ga_result@solution[1, 2]
  results_matrix[i, "Max Fitness"] <- ga_result@fitnessValue
  results_matrix[i, "SSE Mínimo"] <- -ga_result@fitnessValue
}

## Presentación de Resultados

# Combinar las probabilidades iniciales con los resultados
Final_Results <- data.frame(
  Escenario = 1:5,
  P_Target_0.05 = p1,
  P_Target_0.10 = p2,
  results_matrix
)

# alphas <- Final_Results$alpha..a.
# betas <- Final_Results$beta..b.
# num_scenarios <- length(alphas)
# 
# 
# # 2. Definir los puntos de evaluación y colores
# x_vals <- seq(0, 0.20, length.out = 500) # Enfocamos el gráfico en el rango [0, 0.20]
# colors <- c("darkred", "orange", "darkgreen", "blue", "purple")
# legend_labels <- paste0("Escenario ", 1:num_scenarios, 
#                         " (α=", round(alphas, 2), ", β=", round(betas, 2), ")")
# 
# # 3. Crear el primer gráfico (Escenario 1)
# # Calculamos la altura máxima del gráfico para establecer el límite del eje Y
# y_max <- max(dbeta(x_vals, alphas, betas)) * 1.05 
# 
# plot(x_vals, 
#      dbeta(x_vals, alphas[1], betas[1]), 
#      type = "l", 
#      col = colors[1], 
#      lwd = 2,
#      ylim = c(0, y_max),
#      xlim = c(0, 0.25), # Limitamos el eje X hasta 0.15 para mejor visualización
#      xlab = "Valor (X)",
#      ylab = "Densidad de Probabilidad",
#      main = "Funciones de Densidad Beta para 5 Escenarios"
# )
# 
# # 4. Agregar las densidades restantes al gráfico
# for (i in 2:num_scenarios) {
#   lines(x_vals, dbeta(x_vals, alphas[i], betas[i]), col = colors[i], lwd = 2)
# }
# 
# # 5. Agregar los cuantiles AQL y RQL como líneas verticales
# abline(v = 0.05, lty = 2, col = "grey50") # AQL
# abline(v = 0.10, lty = 2, col = "grey50") # RQL
# 
# # 6. Agregar Leyenda
# legend("topright", 
#        legend = legend_labels, 
#        col = colors, 
#        lwd = 2, 
#        cex = 0.8)
# 
# # Agregar la etiqueta de los cuantiles
# text(0.05, -1, "AQL", pos = 4, cex = 0.7, col = "grey50")
# text(0.10, -1, "RQL", pos = 4, cex = 0.7, col = "grey50")

N <- 200       # Tamaño del lote
alpha <- 0.05 # Riesgo del productor
beta <- 0.10  # Riesgo del consumidor

digs <- 5
delta_p <- 10^-digs # Tolerancia para el paso de p (1e-5)

prop <- seq(0, 1, by = delta_p)

ind_aql <- which(prop <= AQL)          # Zona AQL (Aceptación)
ind_rql <- which(prop >= RQL)          # Zona RQL (Rechazo)
ind_zi <- which(prop > AQL & prop < RQL) # Zona de Indiferencia (ZI)

m_hyp = N * prop
n_hyp = N * (1 - prop)

plan_clasic <- find.plan(PRP = c(AQL, 1 - alpha),
                         CRP = c(RQL, beta),
                         N = N, type = "hypergeom")

n_clasic <- plan_clasic$n 
c_clasic <- plan_clasic$c

esce_beta_ <- Final_Results[, c("alpha..a.", "beta..b.")]
num_escenarios <- nrow(esce_beta_)

PA_clasic <- phyper(c_clasic, m_hyp, n_hyp, n_clasic)

# Inicializar la tabla de resultados para los 5 escenarios Beta
resultados_beta <- data.frame(
  Escenario = 1:num_escenarios,
  alpha_beta = esce_beta_[, 1],
  beta_beta = esce_beta_[, 2],
  n_opt = NA,
  c_opt = NA,
  Risk_AQL_opt = NA,
  Risk_ZI_opt = NA,
  Risk_RQL_opt = NA
)

for (dens_idx in 1:num_escenarios) { # dens_idx <- 1 + dens_idx
  
  shape1_val <- esce_beta_[dens_idx, "alpha..a."]
  shape2_val <- esce_beta_[dens_idx, "beta..b."]
  density_eval <- dbeta(prop, shape1 = shape1_val, shape2 = shape2_val)
  
  Risk_AQL_c <- sum((1 - PA_clasic[ind_aql]) * density_eval[ind_aql] * delta_p)
  Risk_ZI_c <- sum((1 - PA_clasic[ind_zi]) * density_eval[ind_zi] * delta_p)
  Risk_RQL_c <- sum(PA_clasic[ind_rql] * density_eval[ind_rql] * delta_p)
  clasics <- c(Risk_AQL_c, Risk_ZI_c, Risk_RQL_c)
  cumple <- FALSE
  tablaz <- matrix(0, nrow = (n_clasic*(n_clasic - 1))/2, ncol = 5)
  rr <- 1
  for (n_ in 1:n_clasic) { # n_ <- 1 + n_
    for (c_ in 0:(n_ - 1)) { # c_ <- 0 + 1 + c_
      
      PA <- phyper(c_, m_hyp, n_hyp, n_) 
      
      Risk_AQL <- sum((1 - PA[ind_aql]) * density_eval[ind_aql] * delta_p)
      Risk_ZI <- sum((1 - PA[ind_zi]) * density_eval[ind_zi] * delta_p)
      Risk_RQL <- sum(PA[ind_rql] * density_eval[ind_rql] * delta_p)
      
      tablaz[rr, 1] <- n_
      tablaz[rr, 2] <- c_
      tablaz[rr, 3] <- Risk_AQL
      tablaz[rr, 4] <- Risk_ZI
      tablaz[rr, 5] <- Risk_RQL
      
      prueba <- all(Risk_AQL < Risk_AQL_c, Risk_ZI < Risk_ZI_c, Risk_RQL < Risk_RQL_c)
      
      if (prueba) {
        # ¡Plan Óptimo Encontrado!
        n_opt <- n_
        c_opt <- c_
        
        # Almacenar resultados
        resultados_beta[dens_idx, "n_opt"] <- n_opt
        resultados_beta[dens_idx, "c_opt"] <- c_opt
        resultados_beta[dens_idx, "Risk_AQL_opt"] <- Risk_AQL
        resultados_beta[dens_idx, "Risk_ZI_opt"] <- Risk_ZI
        resultados_beta[dens_idx, "Risk_RQL_opt"] <- Risk_RQL
        
        cumple <- TRUE
        break # Sale del bucle 'c_'
      }
      rr <- 1 + rr
    }
   if (cumple) {
     break
   }
    
  }
}
tab_filt <- as.data.frame(tablaz) %>% filter(V5 <= 4e-18)
# Mostrar la tabla de resultados
print(resultados_beta)
