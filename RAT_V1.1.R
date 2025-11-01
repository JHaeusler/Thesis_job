# --- 1. Carga de Librerías ---
  
library(AcceptanceSampling)
library(tidyr) 
library(dplyr) # Necesario para bind_rows al recolectar resultados

# --- 2. Parámetros Fijos del Problema ---

N <- 200     # Tamaño del lote
AQL <- 0.05  # Nivel de Calidad Aceptable
RQL <- 0.10  # Nivel de Calidad Rechazable
alpha_clasico <- 0.05
beta_clasico <- 0.10

digs <- 5
delta_p <- 10^-digs 
prop <- seq(0, 1, delta_p) # Vector de proporciones (p)

# Vectores para la Curva CO (Hipergeométrica)

m_hyp <- N * prop 
n_hyp <- N - m_hyp

# Índices para la Integración (RAT)

ind_aql <- which(prop <= AQL) 
ind_rql <- which(prop >= RQL) 

# --- 4. Plan Clásico de Referencia (Constante) ---

plan_clasico <- find.plan(PRP = c(AQL, 1 - alpha_clasico),
                          CRP = c(RQL, beta_clasico),
                          N = N, type = "hypergeom")

n_clasico <- plan_clasico$n
c_clasico <- plan_clasico$c

# --- 5. Función de Optimización ---

dens_eval_p <- dunif(prop, min = min_p, max = max_p)

# b. Riesgos del Plan Clásico (Referencia)
  
CO_clasic <- phyper(q = c_clasico, m = m_hyp, n = n_hyp, k = n_clasico, lower.tail = TRUE)

rp_eval <- sum((1 - CO_clasic[ind_aql]) * dens_eval_p[ind_aql]) * delta_p
rc_eval <- sum(CO_clasic[ind_rql] * dens_eval_p[ind_rql]) * delta_p

# c. Métrica de Pérdida del Clásico (Distancia al origen en el espacio de riesgos)
  
L_eval_sqrt <- rp_eval^2 + rc_eval^2
  
# d. Inicializar variables de búsqueda
  
n_optimo <- NA
c_optimo <- NA
  
Lsqrt_minimo <- Inf # Usamos Inf para asegurar que el primer plan lo mejore

# e. Búsqueda Secuencial (Mínimo n que cumple la condición)

for (n_indx in 1:n_clasico) {
for(c_indx in 0:(n_indx - 1)){

n_actual <- designs_posib[[i, "n"]]
c_actual <- designs_posib[[i, "c"]]

# 1. Curva CO del plan actual

CO_actual <- phyper(q = c_actual, m = m_hyp, n = n_hyp, k = n_actual, lower.tail = TRUE)

# 2. Riesgos del plan actual

rp_actual <- sum((1 - CO_actual[ind_aql]) * dens_eval_p[ind_aql]) * delta_p
rc_actual <- sum(CO_actual[ind_rql] * dens_eval_p[ind_rql]) * delta_p

# 3. Cálculo de la Desviación Cuadrática (Distancia al plan clásico)

L_cuadratic <- (rp_actual - rp_eval)^2 + (rc_actual - rc_eval)^2

# 4. Condición de Parada: El primer plan con L_cuadratic <= L_eval_sqrt

if (L_cuadratic <= L_eval_sqrt) {

n_optimo <- n_actual
c_optimo <- c_actual

Lsqrt_minimo <- L_cuadratic

# Detenemos la búsqueda para asegurar el n mínimo (objetivo del usuario)

break 

}

}

# f. Retornar resultados en un dataframe para fácil comparación

return(data.frame(

Escenario = escenario,

P_Min = min_p,

P_Max = max_p,

n_Clasico = n_clasico,

c_Clasico = c_clasico,

RP_Clasico = rp_eval,

RC_Clasico = rc_eval,

L_Ref_Cuadratica = L_eval_sqrt,

n_Optimo = n_optimo,

c_Optimo = c_optimo,

L_Opt_Cuadratica = Lsqrt_minimo

))




# --- 6. Definición y Ejecución de los 5 Escenarios ---


escenarios_params <- list(

list(min_p = 0.00, max_p = 0.03, nombre = "1. Proveedor Excelente (p muy bajo)"),

list(min_p = 0.00, max_p = 0.08, nombre = "2. Proveedor Bueno (hasta AQL)"),

list(min_p = 0.03, max_p = 0.13, nombre = "3. Proveedor Regular (cerca de AQL/RQL)"),

list(min_p = 0.08, max_p = 0.13, nombre = "4. Proveedor Malo (cercano a RQL)"),

list(min_p = 0.13, max_p = 0.20, nombre = "5. Proveedor Muy Malo (en RQL)")

)


resultados <- data.frame()


for (esc in escenarios_params) {

res <- optimizar_desviacion_cuadratica(esc$min_p, esc$max_p, esc$nombre)

resultados <- bind_rows(resultados, res)

}

