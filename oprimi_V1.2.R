library(AcceptanceSampling)
library(tidyverse)
library(dplyr)
library(tidyr)
library(gridExtra)

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

SCALING_FACTOR_N <- (AQL / RQL) * (beta / N)

max_p <- 0.13
min_p <- 0.03

dens_eval_p <- dunif(prop, min = min_p, max = max_p)

CO_clasic <- phyper(q = c_clasic, m = m_hyp, n = n_hyp, k = n_clasic, lower.tail = TRUE)

rp_clasic <- sum((1 - CO_clasic[ind_aql]) * dens_eval_p[ind_aql] * delta_p)
rc_clasic <- sum(CO_clasic[ind_rql] * dens_eval_p[ind_rql] * delta_p)

zi_acceptance_risk_clasic <- sum(CO_clasic[ind_zi] * dens_eval_p[ind_zi] * delta_p)

penalizacion_clasic_zi_severidad <- SCALING_FACTOR_N * (n_clasic / (c_clasic + 1)) * zi_acceptance_risk_clasic

penalizacion_clasic_costo <- SCALING_FACTOR_N * ((n_clasic+c_clasic)/N) * (sum(dens_eval_p[ind_aql]*delta_p) + sum(dens_eval_p[ind_rql]*delta_p))

L_clasico <- (rp_clasic + rc_clasic) + penalizacion_clasic_zi_severidad + penalizacion_clasic_costo

RAT_Clasico_Base <- rp_clasic + rc_clasic 

L_actual <- 0

i <- 1 ; n_v <- 0; c_v <- 0;RAT_actual <- 0
penalizacion_actual_costo <- 0

for (n_actual in 1:20) { # n_actual <- 1 + n_actual
  for (c_actual in 0:(n_actual - 1)) { # c_actual <- 0 + 1 + c_actual
   
    CO_actual <- phyper(q = c_actual, m = m_hyp, n = n_hyp, k = n_actual, lower.tail = TRUE)
    
    rp_actual <- sum((1 - CO_actual[ind_aql]) * dens_eval_p[ind_aql] * delta_p)
    rc_actual <- sum(CO_actual[ind_rql] * dens_eval_p[ind_rql] * delta_p)
    
    RAT_actual[i] <- rp_actual + rc_actual

    zi_acceptance_risk_actual <-  sum(CO_actual[ind_zi] * dens_eval_p[ind_zi] * delta_p)

    # Penalización 1: Severidad del Plan Actual * Riesgo ZI Actual * lambda
    
    penalizacion_actual_zi_severidad <- SCALING_FACTOR_N * (n_actual / (c_actual + 1)) * (zi_acceptance_risk_actual)
    
    # Penalización 2: Costo Operacional (Nuevo término)
    
    penalizacion_actual_costo[i] <- SCALING_FACTOR_N * ((n_actual+c_actual)/N) * (sum(dens_eval_p[ind_aql]*delta_p) + sum(dens_eval_p[ind_rql]*delta_p))
 
    L_actual[i] <- RAT_actual[i] + penalizacion_actual_zi_severidad + penalizacion_actual_costo[i]
    n_v[i] <- n_actual 
    c_v[i] <- c_actual
    i <- 1 + i
  }}     

dat <- cbind(n_v, c_v, L_clasico, L_actual, penalizacion_clasic_costo, penalizacion_actual_costo) 
dat
#which(dat[,4] < L_clasico)

