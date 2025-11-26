

lista_de_paquetes <- c("pracma", "GA", "AcceptanceSampling",
                       "scatterplot3d", "tidyverse", "reshape",
                       "sf", "dplyr", "RColorBrewer") # Reemplaza con tus paquetes

for (paquete in lista_de_paquetes) {
  if(!require(paquete, character.only = TRUE)){
    install.packages(paquete)
    require(paquete, character.only = TRUE)
  }
}

# semilla
set.seed(060722)
percen <- 0.05
escen_unif <- data.frame(rbind(sort(runif(2)),
                    c(0.00, 0.03),
                    c(0.03, 0.07),
                    c(0.03, 0.13),
                    c(0.07, 0.13),
                    c(0.13, 0.20)))

colnames(escen_unif) <- c("t_min", "t_max")

esce_beta <- matrix(NA, nrow = nrow(escen_unif), ncol = ncol(escen_unif))
for(i in 1:dim(escen_unif)[1]){
  E_t <- sum(escen_unif[i, ])/2
  sigma_sqrt <- percen*((escen_unif[i, 2] - E_t)^2)
  w <- E_t / (1 - E_t)
  esce_beta[i, 2] <- (w - (((w + 1)^2) * sigma_sqrt))/((w + 1)^3 * sigma_sqrt)
  esce_beta[i, 1] <- esce_beta[i, 2] * w
}

escen_unif_ <- cbind(as.data.frame(escen_unif), density = "uniform")
esce_beta_ <- cbind(as.data.frame(esce_beta), density = "beta")
colnames(esce_beta_)[1:2] <- c("t_min", "t_max")
compl_df <- rbind(escen_unif_, esce_beta_)

N <- 200       # Tamaño del lote
alpha <- 0.05 # Riesgo del productor
beta <- 0.10  # Riesgo del consumidor
AQL <- 0.05   # Nivel de calidad aceptable
RQL <- 0.10   # Nivel de calidad de rechazo

digs <- 5
delta_p <- 10^-digs # Tolerancia para el paso de p (1e-5)

prop <- seq(0, 1, by = delta_p)

ind_aql <- which(prop <= AQL)          # Zona AQL (Aceptación)
ind_rql <- which(prop >= RQL)          # Zona RQL (Rechazo)
ind_zi <- which(prop > AQL & prop < RQL) # Zona de Indiferencia (ZI)

m_hyp = N * prop
n_hyp = N * (1 - prop)


# for (dens_ in 1:nrow(compl_df)) { # dens_ <- 8 + dens_
#   
#   if (compl_df[dens_, 3] == "uniform") {
#     density_eval <- dunif(prop, min = compl_df[dens_, 1], max = compl_df[dens_, 2])
#   } else{
#     density_eval <- dbeta(prop, shape1 = compl_df[dens_, 1], shape2 = compl_df[dens_, 2])
#   }
#   
#   for(j in 1:80){
#     for (k in 0:79) {
#       n_ <- j
#       c_ <- k
#       
#       CO <- phyper(c_, m_hyp, n_hyp, n_)
#     
#       risk[j, 1] <- sum((1 - CO[ind_aql])*density_eval[ind_aql]*delta_p)
#       risk[j, 2] <- sum((1 - CO[ind_zi])*density_eval[ind_zi]*delta_p)
#       risk[j, 3] <- sum((CO[ind_rql])*density_eval[ind_rql]*delta_p)
#       
#    } 
#   }
# 
# }

# Inicializar la lista para almacenar los resultados de los 12 escenarios
resultados_escenarios <- vector("list", nrow(compl_df))
names(resultados_escenarios) <- paste("Escenario", 1:nrow(compl_df), 
                                      compl_df$density, round(compl_df$t_min, 4), round(compl_df$t_max, 4))
n_val <- 30
c_val <- n_val - 1
# Bucle principal: Itera sobre los 12 escenarios (6 Uniformes y 6 Beta)
for (dens_idx in 1:nrow(compl_df)) {
  
  # 1. Calcular la densidad (uniforme o beta) para el vector 'prop'
  if (compl_df[dens_idx, "density"] == "uniform") {
    # Usar nombres de columna para mayor claridad
    min_val <- compl_df[dens_idx, "t_min"]
    max_val <- compl_df[dens_idx, "t_max"]
    density_eval <- dunif(prop, min = min_val, max = max_val)
  } else {
    # Usar nombres de columna para mayor claridad
    shape1_val <- compl_df[dens_idx, "t_min"] # En Beta, t_min/t_max son alfa/beta
    shape2_val <- compl_df[dens_idx, "t_max"]
    # Se añade un manejo de errores por si los parámetros Beta son NA o no válidos
    if (is.na(shape1_val) || is.na(shape2_val) || shape1_val <= 0 || shape2_val <= 0) {
      warning(paste("Parámetros Beta inválidos o NA en Escenario", dens_idx, ". Saltando."))
      resultados_escenarios[[dens_idx]] <- NULL
      next
    }
    density_eval <- dbeta(prop, shape1 = shape1_val, shape2 = shape2_val)
  }
  
  # 2. Inicializar la tabla de riesgos para este escenario (n, c, Riesgo AQL, Riesgo ZI, Riesgo RQL)
  # El número de filas será (longitud de n_values) * (longitud de c_values)
  num_planes <- 30 * 29
  # Creamos un data.frame vacío para pre-asignar y almacenar los resultados
  risk_table <- data.frame(
    n = rep(NA, num_planes),
    c = rep(NA, num_planes),
    Risk_AQL = rep(NA, num_planes),
    Risk_ZI = rep(NA, num_planes),
    Risk_RQL = rep(NA, num_planes)
  )
  
  row_idx <- 1 # Contador de filas para la tabla de riesgos
  
  # 3. Bucle para evaluar los distintos planes (n, c)
  for (n_ in 1:n_val) {
    for (c_ in 0:(n_ - 1)) {
      
      # Calcular la Probabilidad de Aceptación (PA) para cada p en el vector 'prop'
      # PA = P(X <= c | p) = P(X <= c | m, n, n_)
      # phyper(q=c_, m=defectuosos, n=no_defectuosos, k=muestra)
      PA <- phyper(c_, m_hyp, n_hyp, n_) 
      
      # 4. Calcular los riesgos (Riesgo = Sumatoria [P(Rechazo/Acept) * f(p) * dp])
      # NOTA: La probabilidad de rechazo es (1 - PA)
      
      # Riesgo AQL (Riesgo de Rechazo en Zona de Aceptación)
      Risk_AQL <- sum((1 - PA[ind_aql]) * density_eval[ind_aql] * delta_p)
      
      # Riesgo ZI (Riesgo de Rechazo en Zona de Indiferencia)
      Risk_ZI <- sum((1 - PA[ind_zi]) * density_eval[ind_zi] * delta_p)
      
      # Riesgo RQL (Riesgo de Aceptación en Zona de Rechazo)
      Risk_RQL <- sum(PA[ind_rql] * density_eval[ind_rql] * delta_p)
      
      # 5. Almacenar resultados en la tabla
      risk_table[row_idx, ] <- c(n_, c_, Risk_AQL, Risk_ZI, Risk_RQL)
      
      row_idx <- row_idx + 1
    }
  }
  
  # 6. Almacenar la tabla de riesgos completa para el escenario actual en la lista
  resultados_escenarios[[dens_idx]] <- risk_table
}

# Mostrar la estructura de los resultados (los primeros 5 planes del primer escenario)
cat("✅ Proceso de cálculo de riesgos completado.\n\n")
cat("--- Estructura de la Lista de Resultados ---\n")
print(names(resultados_escenarios))
cat("\nPrimeros 5 Planes de Muestreo (n, c) del Escenario 1 (Uniforme):\n")
print(head(resultados_escenarios[[8]], 200))

plan_clasic <- find.plan(PRP = c(AQL, 1 - alpha),
                         CRP = c(RQL, beta),
                         N = N, type = "hypergeom")

n_clasic <- plan_clasic$n 
c_clasic <- plan_clasic$c

PA_clasic <- phyper(c_clasic, m_hyp, n_hyp, n_clasic)

resultados_clasicos <- data.frame(
  Escenario = character(nrow(compl_df)),
  Densidad = character(nrow(compl_df)),
  t_min = numeric(nrow(compl_df)),
  t_max = numeric(nrow(compl_df)),
  Risk_AQL = numeric(nrow(compl_df)),
  Risk_ZI = numeric(nrow(compl_df)),
  Risk_RQL = numeric(nrow(compl_df))
)

# Bucle principal: Itera sobre los 12 escenarios
for (dens_idx in 1:nrow(compl_df)) {
  
  # 1. Calcular la densidad (uniforme o beta) para el vector 'prop'
  if (compl_df[dens_idx, "density"] == "uniform") {
    min_val <- compl_df[dens_idx, "t_min"]
    max_val <- compl_df[dens_idx, "t_max"]
    density_eval <- dunif(prop, min = min_val, max = max_val)
    
  } else { # Densidad Beta
    shape1_val <- compl_df[dens_idx, "t_min"]
    shape2_val <- compl_df[dens_idx, "t_max"]
    
    # Manejo de errores para parámetros Beta inválidos (necesario del paso anterior)
    if (is.na(shape1_val) || is.na(shape2_val) || shape1_val <= 0 || shape2_val <= 0) {
      warning(paste("Parámetros Beta inválidos en Escenario", dens_idx, ". Saltando cálculo."))
      next
    }
    density_eval <- dbeta(prop, shape1 = shape1_val, shape2 = shape2_val)
  }
  
  # 2. Calcular la Probabilidad de Aceptación (PA) para el plan clásico (n_clasic, c_clasic)
  # PA = P(Aceptación | p) usando la Hipergeométrica
  PA_clasic <- phyper(c_clasic, m_hyp, n_hyp, n_clasic)
  
  
  # 3. Calcular los tres riesgos (Integración numérica con delta_p)
  
  # Riesgo AQL (Riesgo de Rechazo en Zona de Aceptación)
  # R_AQL = Integral[ P(Rechazo | p) * f(p) dp ] sobre la Zona AQL
  Risk_AQL_c <- sum((1 - PA_clasic[ind_aql]) * density_eval[ind_aql] * delta_p)
  
  # Riesgo ZI (Riesgo de Rechazo en Zona de Indiferencia)
  # R_ZI = Integral[ P(Rechazo | p) * f(p) dp ] sobre la Zona ZI
  Risk_ZI_c <- sum((1 - PA_clasic[ind_zi]) * density_eval[ind_zi] * delta_p)
  
  # Riesgo RQL (Riesgo de Aceptación en Zona de Rechazo)
  # R_RQL = Integral[ P(Aceptación | p) * f(p) dp ] sobre la Zona RQL
  Risk_RQL_c <- sum(PA_clasic[ind_rql] * density_eval[ind_rql] * delta_p)
  
  # 4. Almacenar los resultados
  resultados_clasicos[dens_idx, "Escenario"] <- paste("Esc", dens_idx)
  resultados_clasicos[dens_idx, "Densidad"] <- compl_df[dens_idx, "density"]
  resultados_clasicos[dens_idx, "t_min"] <- compl_df[dens_idx, "t_min"]
  resultados_clasicos[dens_idx, "t_max"] <- compl_df[dens_idx, "t_max"]
  resultados_clasicos[dens_idx, "Risk_AQL"] <- Risk_AQL_c
  resultados_clasicos[dens_idx, "Risk_ZI"] <- Risk_ZI_c
  resultados_clasicos[dens_idx, "Risk_RQL"] <- Risk_RQL_c
}

# Crear una columna de Riesgo Total para facilitar la comparación
resultados_clasicos$Risk_Total <- (resultados_clasicos$Risk_AQL + 
                                     resultados_clasicos$Risk_ZI + 
                                     resultados_clasicos$Risk_RQL)
