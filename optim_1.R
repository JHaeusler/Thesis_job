

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

plan_clasic <- find.plan(PRP = c(AQL, 1 - alpha),
                         CRP = c(RQL, beta),
                         N = N, type = "hypergeom")

n_clasic <- plan_clasic$n 
c_clasic <- plan_clasic$c

# Inicializar la lista para almacenar los resultados de los 12 escenarios
resultados_escenarios <- vector("list", nrow(compl_df))
names(resultados_escenarios) <- paste("Escenario", 1:nrow(compl_df), 
                                      compl_df$density, round(compl_df$t_min, 4), round(compl_df$t_max, 4))
n_val <- n_clasic
c_val <- n_val - 1
# Bucle principal: Itera sobre los 12 escenarios (6 Uniformes y 6 Beta)
for (dens_idx in 1:nrow(esce_beta_)) {

  shape1_val <- esce_beta_[dens_idx, "t_min"]
  shape2_val <- esce_beta_[dens_idx, "t_max"]
  density_eval <- dbeta(prop, shape1 = shape1_val, shape2 = shape2_val)

  PA_clasic <- phyper(c_clasic, m_hyp, n_hyp, n_clasic)
  
  Risk_AQL_c <- sum((1 - PA_clasic[ind_aql]) * density_eval[ind_aql] * delta_p)
  Risk_ZI_c <- sum((1 - PA_clasic[ind_zi]) * density_eval[ind_zi] * delta_p)
  Risk_RQL_c <- sum(PA_clasic[ind_rql] * density_eval[ind_rql] * delta_p)

  for (n_ in 1:n_val) {
    for (c_ in 0:(n_ - 1)) {
      
      PA <- phyper(c_, m_hyp, n_hyp, n_) 
      
      Risk_AQL <- sum((1 - PA[ind_aql]) * density_eval[ind_aql] * delta_p)

      Risk_ZI <- sum((1 - PA[ind_zi]) * density_eval[ind_zi] * delta_p)

      Risk_RQL <- sum(PA[ind_rql] * density_eval[ind_rql] * delta_p)
      
      prueba <- sum(c(Risk_AQL < Risk_AQL_c, Risk_ZI < Risk_ZI_c, Risk_RQL < Risk_RQL_c))
      
      if (prueba == 3) {
        n_opt <- n_
        c_opt <- c_
        
      }
      
      
    }
  }
  

}
