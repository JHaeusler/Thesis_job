
# excelente proveedor .95 (AQL) - 0.0001 (LTPD)
# buen proveedor      .75 (AQL) - 0.05 (LTPD)     
# regular             .50 (AQL) - .10 (LTPD)
# malo                .25 (AQL) - .25
# muy malo            .10 (AQL) - .40

# install.packages("pracma") # Si no lo tienes instalado
library(pracma)
# install.packages("sjstats")
library(sjstats)
library(AcceptanceSampling)

# --- Variables Globales ---
N <- 200        # Tamaño del lote
alpha <- 0.05   # Riesgo del productor (para plan clasico)
beta <- 0.10    # Riesgo del consumidor (para plan clasico)
AQL=0.05
LTPD=0.10

# Vectores de probabilidades para los 5 escenarios (Para encontrar alpha y beta)
p1 <- c(0.95, 0.75, 0.50, 0.25, 0.10)
p2 <- 1 - c(1e-4, 0.05, 0.10, 0.25, 0.40)

Pa <- function(n, c, p, N=200){

  Pa_val = phyper(c, N*p, N*(1-p), n) 
  return(Pa_val)
}

wr_p <- function(p, n, c, alpha_b, beta_b, N){
  f.p = dbeta(p, alpha_b, beta_b)
  prod_wr <- (1 - Pa(n, c, p, N)) * f.p
  return(prod_wr)
}

wr_c <- function(p, n, c, alpha_b, beta_b, N){
  f.p = dbeta(p, alpha_b, beta_b)
  cons_wr <- Pa(n, c, p, N) * f.p
  return(cons_wr)
}


calc_wr <- function(n, c, alpha_b, beta_b, AQL, LTPD) {
  

  rp_val <- integral(f = function(p) wr_p(p, n_clasic, c_clasic, a, b, N),
                       xmin = 0, xmax = AQL, method = "Kron")
             
  rc_val <- integral(f = function(p) wr_c(p, n_clasic, c_clasic, a, b, N),
                         xmin = LTPD, xmax = 1, method = "Kron")
  
  # La función integral de pracma devuelve solo el valor, no una lista.
  return(c(RP_val = rp_val, RC_val = rc_val))
}


# 1. Determinar el Plan Clásico
plan_clasic <- find.plan(PRP = c(AQL, 1 - alpha),
                         CRP = c(LTPD, beta),
                         N = N, type = "hypergeom")

n_clasic <- plan_clasic$n
c_clasic <- plan_clasic$c

delta_p <- 1e-5
p <- seq(0, 1, delta_p)

Shape <- find_beta(x1=AQL, p1=p1[1], x2=LTPD, p2=p2[1])
a <- Shape$shape1
b <- Shape$shape2

pa <- phyper(c_clasic, N*p, N*(1-p), n_clasic)

denst <- dbeta(p, a, b)


sum((1-pa[p<=AQL])*denst[p<=AQL]*delta_p)

integral(f = function(p) wr_p(p, n_clasic, c_clasic, a, b, N),
         xmin = 0, xmax = AQL, method = "Kron")

sum((pa[p>LTPD])*denst[p>LTPD]*delta_p)

integral(f = function(p) wr_c(p, n_clasic, c_clasic, a, b, N),
         xmin = LTPD, xmax = 1, method = "Kron")

calc_wr(n = n_clasic, c = c_clasic, alpha_b = a,
        beta_b = b)



sum((1-pa[p <= AQL])*dbeta(p, a, b)[p <= AQL])*delta_p



# 2. Obtener los 5 pares de (alpha, beta)
# Almacenaremos los pares en una matriz/data.frame para iterar
alpha_beta_params <- data.frame(alpha_b = rep(NA, 5), beta_b = rep(NA, 5))

for (i in 1:5) {
  Shape <- find_beta(x1=AQL, p1=p1[i], x2=LTPD, p2=p2[i])
  alpha_beta_params[i, "alpha_b"] <- Shape$shape1
  alpha_beta_params[i, "beta_b"] <- Shape$shape2
}

# 3. Inicializar la tabla de resultados para los riesgos
resultados_riesgo <- data.frame(
  Escenario = 1:5,
  alpha_b = alpha_beta_params$alpha_b,
  beta_b = alpha_beta_params$beta_b,
  n_clasic = n_clasic,
  c_clasic = c_clasic,
  RP_clasic = NA,
  RC_clasic = NA
)

# 4. Iterar sobre los 5 escenarios y calcular los riesgos ponderados
for (i in 1:5) {
  alpha_b_val <- resultados_riesgo[i, "alpha_b"]
  beta_b_val <- resultados_riesgo[i, "beta_b"]
  
  # Calcula los riesgos RP y RC para el plan clásico bajo la densidad Beta actual
  risks <- calc_wr(n = n_clasic, c = c_clasic, 
                   alpha_b = alpha_b_val, beta_b = beta_b_val)
  
  resultados_riesgo[i, "RP_clasic"] <- risks["RP_val"]
  resultados_riesgo[i, "RC_clasic"] <- risks["RC_val"]
}


# Asumiendo que el bucle for (i in 1:5) { ... } está activo
# y que las variables resultados_riesgo, alpha_b_val, beta_b_val están definidas.

# 1. Variables para almacenar el plan óptimo encontrado para el escenario actual (i)
n_opt_found <- NA
c_opt_found <- NA
cumple <- FALSE # Bandera para detener la búsqueda

# 2. Búsqueda Secuencial
# El límite superior para n debe ser N (200), no n_clasic, para asegurar que se encuentra la mejor solución.
for (n_ in 1:n_clasic) { 
  
  # Bucle C corregido: c_ debe ir de 0 a n_ - 1
  for (c_ in 0:(n_ - 1)) {
    
    # 3. Calcular los riesgos ponderados para el plan (n_, c_)
    risks_opt <- calc_wr(n = n_, c = c_, 
                         alpha_b = alpha_b_val, beta_b = beta_b_val)
    
    # 4. Condición de Cumplimiento (Mejora en ambos riesgos sobre el plan clásico)
    # Se compara el RP con el RP clásico y el RC con el RC clásico.
    if (risks_opt["RP_val"] < resultados_riesgo[i, "RP_clasic"] && 
        risks_opt["RC_val"] < resultados_riesgo[i, "RC_clasic"]) {
      
      # Plan óptimo (el de n más pequeño) encontrado
      n_opt_found <- n_
      c_opt_found <- c_
      cumple <- TRUE
      
      # El primer plan que cumple el criterio es el de menor 'n' y 'c'
      break # Sale del bucle 'c_'
    }
  }
  
  if (cumple) {
    break # Sale del bucle 'n_'
  }
}

# 5. Almacenar los resultados óptimos en la tabla
resultados_riesgo[i, "n_opt"] <- n_opt_found
resultados_riesgo[i, "c_opt"] <- c_opt_found
# Recalcular y almacenar los riesgos óptimos para el n_opt, c_opt encontrado (si aplica)
if (cumple) {
  risks_opt_final <- calc_wr(n = n_opt_found, c = c_opt_found, 
                             alpha_b = alpha_b_val, beta_b = beta_b_val)
  resultados_riesgo[i, "RP_opt"] <- risks_opt_final["RP_val"]
  resultados_riesgo[i, "RC_opt"] <- risks_opt_final["RC_val"]
}
