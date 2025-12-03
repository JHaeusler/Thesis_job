# Script para comparar planes de muestreo de aceptación clásicos
# con planes óptimos basados en el Riesgo Ponderado (Enfoque Bayesiano).

# Las prioris de calidad de los proveedores se definen por pares (PA en AQL, PA en LTPD)
# Excelente proveedor: PA(AQL)=0.95, PA(LTPD)=1-0.0001
# Buen proveedor:      PA(AQL)=0.75, PA(LTPD)=1-0.05
# Regular:             PA(AQL)=0.50, PA(LTPD)=1-0.10
# Malo:                PA(AQL)=0.25, PA(LTPD)=1-0.25
# Muy malo:            PA(AQL)=0.10, PA(LTPD)=1-0.40

# --- Carga de Librerías ---
# install.packages("pracma") # Si no lo tienes instalado
library(pracma)
# install.packages("sjstats")
library(sjstats)
library(AcceptanceSampling)

# --- Variables Globales y Parámetros del Problema ---
N <- 200        # Tamaño del lote
alpha <- 0.05   # Riesgo del productor (para el plan clásico)
beta <- 0.10    # Riesgo del consumidor (para el plan clásico)
AQL <- 0.05     # Nivel de Calidad Aceptable
LTPD <- 0.10    # Tolerancia de Porcentaje Defectuoso de Lote

# Vectores de probabilidades para los 5 escenarios
# p1 = Probabilidad de Aceptación en el AQL
# p2 = Probabilidad de No Aceptación (o Rechazo) en el LTPD
p1 <- c(0.95, 0.75, 0.50, 0.25, 0.10)
p2 <- c(1e-4, 0.05, 0.10, 0.25, 0.40) # 1 - PA(LTPD) -> Beta (Riesgo del Consumidor)

# --- Funciones de Riesgo Ponderado ---

# Probabilidad de Aceptación (PA) - Curva CO tipo A (Hypergeométrica)
# Usa N (tamaño del lote) definido globalmente
Pa <- function(n, c, p, N){

  # El paquete sjstats redondea X a un entero para la distribución hipergeométrica.
  Pa_val <- phyper(c, N * p, N * (1 - p), n)
  return(Pa_val)
}

# Integrando el Riesgo del Productor (RP): (1 - PA) * f(p)
# RP es la probabilidad de rechazar un lote de buena calidad (p <= AQL).
wr_p <- function(p, n, c, alpha_b, beta_b, N){
  f.p <- dbeta(p, alpha_b, beta_b) # Densidad a priori Beta
  prod_wr <- (1 - Pa(n, c, p, N)) * f.p
  return(prod_wr)
}

# Integrando el Riesgo del Consumidor (RC): PA * f(p)
# RC es la probabilidad de aceptar un lote de mala calidad (p >= LTPD).
wr_c <- function(p, n, c, alpha_b, beta_b, N){
  f.p <- dbeta(p, alpha_b, beta_b) # Densidad a priori Beta
  cons_wr <- Pa(n, c, p, N) * f.p
  return(cons_wr)
}

# Función principal para calcular el Riesgo Ponderado Integrado
calc_wr <- function(n, c, alpha_b, beta_b, AQL, LTPD) {

  # Riesgo del Productor (RP) Integrado [0, AQL]
  rp_val <- integral(f = function(p) wr_p(p, n, c, alpha_b, beta_b, N), 
                     xmin = 0, xmax = AQL, 
                     method = "Kron") # Método de Gauss-Kronrod (preciso)
  
  # Riesgo del Consumidor (RC) Integrado [LTPD, 1]
  rc_val <- integral(f = function(p) wr_c(p, n, c, alpha_b, beta_b, N), 
                     xmin = LTPD, xmax = 1, 
                     method = "Kron")
  
  return(c(RP_val = rp_val, RC_val = rc_val))
}


# 1. Determinar el Plan Clásico (basado en alpha y beta fijos)
plan_clasic <- find.plan(PRP = c(AQL, 1 - alpha),
                         CRP = c(LTPD, beta),
                         N = N, type = "hypergeom")

n_clasic <- plan_clasic$n
c_clasic <- plan_clasic$c

# 2. Obtener los 5 pares de (alpha_b, beta_b) para la distribución Beta
# La distribución Beta debe cumplir: P(p <= AQL) = p1[i] y P(p >= LTPD) = p2[i]
alpha_beta_params <- data.frame(alpha_b = rep(NA, 5), beta_b = rep(NA, 5),
                                acum_aql = rep(NA, 5), acum_ltpd = rep(NA, 5))

for (i in 1:5) {
  # find_beta encuentra los parámetros alpha_b y beta_b que satisfacen la probabilidad acumulada
  Shape <- find_beta(x1 = AQL, p1 = p1[i], x2 = LTPD, p2 = 1 - p2[i]) # p2 es 1-PA(LTPD), find_beta espera P(X<=x2)
  alpha_beta_params[i, "alpha_b"] <- Shape$shape1
  alpha_beta_params[i, "beta_b"] <- Shape$shape2
  alpha_beta_params[i, "acum_aql"] <- pbeta(AQL, Shape$shape1, Shape$shape2) # Debe ser cercano a p1[i]
  # Acumulada del lado derecho (Riesgo del consumidor)
  alpha_beta_params[i, "acum_ltpd"] <- 1 - pbeta(LTPD, Shape$shape1, Shape$shape2) # Debe ser cercano a p2[i]
}

# 3. Inicializar la tabla de resultados
resultados_riesgo <- data.frame(
  Escenario = 1:5,
  alpha_b = alpha_beta_params$alpha_b,
  beta_b = alpha_beta_params$beta_b,
  acum_aql = alpha_beta_params$acum_aql,
  acum_ltpd = alpha_beta_params$acum_ltpd,
  n_clasic = n_clasic,
  c_clasic = c_clasic,
  RP_clasic = NA,
  RC_clasic = NA,
  n_opt = NA,
  c_opt = NA,
  RP_opt = NA,
  RC_opt = NA
)

# 4. Iterar sobre los 5 escenarios y calcular los riesgos ponderados para cada plan
for (i in 1:5) {
  
  alpha_b_val <- resultados_riesgo[i, "alpha_b"]
  beta_b_val <- resultados_riesgo[i, "beta_b"]
  
  # I. Calcular los riesgos RP y RC para el PLAN CLÁSICO bajo la densidad Beta actual
  risks_clasic <- calc_wr(n = n_clasic, c = c_clasic, 
                          alpha_b = alpha_b_val, beta_b = beta_b_val, AQL, LTPD)
  
  resultados_riesgo[i, "RP_clasic"] <- risks_clasic["RP_val"]
  resultados_riesgo[i, "RC_clasic"] <- risks_clasic["RC_val"]
  
  # --- II. Búsqueda del Plan Óptimo para el Escenario i ---
  
  # Reestablecer variables de búsqueda para cada escenario
  n_opt_found <- NA
  c_opt_found <- NA
  cumple <- FALSE 
  
  # Búsqueda Secuencial: Itera n de 1 hasta n_clasic (o N si se quiere ser más exhaustivo)
  for (n_ in 1:n_clasic) { 
    
    # c_ debe ser menor que n_ (el número máximo de defectuosos debe ser menor que la muestra)
    for (c_ in 0:(n_ - 1)) {
      
      # Calcular los riesgos ponderados para el plan (n_, c_)
      risks_opt <- calc_wr(n = n_, c = c_, 
                           alpha_b = alpha_b_val, beta_b = beta_b_val, AQL, LTPD)
      
      # Condición de Cumplimiento (Busca el plan con n más pequeño que mejore a ambos riesgos)
      if (risks_opt["RP_val"] < resultados_riesgo[i, "RP_clasic"] && 
          risks_opt["RC_val"] < resultados_riesgo[i, "RC_clasic"]) {
        
        # Plan óptimo encontrado
        n_opt_found <- n_
        c_opt_found <- c_
        cumple <- TRUE
        
        break # Sale del bucle 'c_' (encontró el c_ más pequeño para ese n_)
      }
    }
    
    if (cumple) {
      break # Sale del bucle 'n_' (encontró el n_ más pequeño que cumple)
    }
  }
  
  # III. Almacenar los resultados óptimos en la tabla
  resultados_riesgo[i, "n_opt"] <- n_opt_found
  resultados_riesgo[i, "c_opt"] <- c_opt_found
  
  # Recalcular y almacenar los riesgos óptimos para el n_opt, c_opt encontrado (si aplica)
  if (cumple) {
    risks_opt_final <- calc_wr(n = n_opt_found, c = c_opt_found, 
                               alpha_b = alpha_b_val, beta_b = beta_b_val, AQL, LTPD)
    resultados_riesgo[i, "RP_opt"] <- risks_opt_final["RP_val"]
    resultados_riesgo[i, "RC_opt"] <- risks_opt_final["RC_val"]
  } else {
    # Si no se encuentra un plan óptimo, se rellenan con NA
    resultados_riesgo[i, "RP_opt"] <- NA
    resultados_riesgo[i, "RC_opt"] <- NA
  }
}

# 5. Mostrar la tabla de resultados final
# Se recomienda usar opciones para mostrar más decimales en los riesgos
options(digits = 6) 
print("Tabla de Comparación de Riesgos Ponderados")
print("------------------------------------------")
print(resultados_riesgo)