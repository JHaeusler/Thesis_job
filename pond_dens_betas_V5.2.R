# --- Carga de Librerías ---
library(pracma)
library(sjstats)
library(AcceptanceSampling)
library(GA)

# --- Configuración de Parámetros ---
N_vec <- c(200, 1000)
alpha_des <- c(0.01, 0.03, 0.05)
beta_des <- c(0.07, 0.10, 0.15, 0.20)
AQL_vec <- c(0.01, 0.03, 0.05)
LTPD_vec <- c(0.05, 0.10, 0.15)

# Generación de la matriz de combinaciones
combinaciones <- expand.grid(N = N_vec, AQL = AQL_vec, LTPD = LTPD_vec, 
                             alpha = alpha_des, beta = beta_des)


# --- Aplicación de la restricción: AQL debe ser estrictamente menor que LTPD ---
combinaciones <- combinaciones[combinaciones$AQL < combinaciones$LTPD, ]

# (Opcional) Reiniciar los índices del dataframe para que sean correlativos
rownames(combinaciones) <- NULL


# Configuración de Perfiles de Proveedor
perfiles <- data.frame(
  Proveedor = c("Excelente", "Bueno", "Regular", "Malo", "Muy Malo"),
  p1 = c(0.95, 0.75, 0.50, 0.25, 0.10),
  p2 = c(1e-4, 0.05, 0.10, 0.25, 0.40)
)

# Lista para almacenar los resultados
lista_resultados <- list()

# --- Funciones Auxiliares ---
calc_prob_mass <- function(alpha_b, beta_b, AQL, LTPD) {
  P_Good <- pbeta(AQL, shape1 = alpha_b, shape2 = beta_b)
  P_Bad <- 1 - pbeta(LTPD, shape1 = alpha_b, shape2 = beta_b)
  return(c(P_Good = P_Good, P_Bad = P_Bad))
}

calc_wr <- function(n, c, N_val, AQL_val, LTPD_val, k_p, k_c) {
  wrp_val <- k_p * (1 - phyper(c, N_val * AQL_val, N_val * (1 - AQL_val), n))
  wrc_val <- k_c * phyper(c, N_val * LTPD_val, N_val * (1 - LTPD_val), n)
  return(c(WRP_val = wrp_val, WRC_val = wrc_val, WR_val = wrp_val + wrc_val))
}

# --- Ciclo Principal ---
for (j in 1:nrow(combinaciones)) {
  
  # Extraer parámetros actuales
  curr <- combinaciones[j, ]
  
  # 1. Determinar Plan Clásico para esta combinación
  plan_clasic <- find.plan(PRP = c(curr$AQL, 1 - curr$alpha),
                           CRP = c(curr$LTPD, curr$beta),
                           N = curr$N, type = "hypergeom")
  
  n_cl <- plan_clasic$n
  c_cl <- plan_clasic$c
  
  # 2. Preparar tabla para los 5 proveedores + Naive
  res_escenario <- data.frame(
    Proveedor = c(perfiles$Proveedor, "Naive"),
    n_cl = n_cl, c_cl = c_cl,
    n_opt = NA, c_opt = NA,
    WR_cl = NA, WR_opt = NA, Ganancia = NA
  )
  
  # 3. Iterar por cada perfil de proveedor dentro de la combinación
  for (i in 1:6) {
    # Determinar parámetros Beta
    if (i <= 5) {
      Shape <- find_beta(x1 = curr$AQL, p1 = perfiles$p1[i], 
                         x2 = curr$LTPD, p2 = 1 - perfiles$p2[i])
      a_b <- Shape$shape1; b_b <- Shape$shape2
    } else {
      a_b <- 1; b_b <- 1 # Caso Naive (Uniforme)
    }
    
    mass <- calc_prob_mass(a_b, b_b, curr$AQL, curr$LTPD)
    kp <- as.numeric(mass["P_Good"]); kc <- as.numeric(mass["P_Bad"])
    
    # Riesgo Plan Clásico
    r_cl <- calc_wr(n_cl, c_cl, curr$N, curr$AQL, curr$LTPD, kp, kc)
    res_escenario$WR_cl[i] <- r_cl["WR_val"]
    
    # Búsqueda Plan Óptimo (Simplificada)
    min_wr <- r_cl["WR_val"]
    n_best <- n_cl; c_best <- c_cl
    
    # Optimización local (puedes reintegrar el GA aquí si prefieres)
    for (n_test in 1:curr$N) {
      for (c_test in 0:min(n_test-1, 20)) { # Límite de búsqueda en c para velocidad
        r_test <- calc_wr(n_test, c_test, curr$N, curr$AQL, curr$LTPD, kp, kc)
        if (r_test["WR_val"] < min_wr) {
          min_wr <- r_test["WR_val"]
          n_best <- n_test; c_best <- c_test
        }
      }
    }
    
    res_escenario$n_opt[i] <- n_best
    res_escenario$c_opt[i] <- c_best
    res_escenario$WR_opt[i] <- min_wr
    res_escenario$Ganancia[i] <- res_escenario$WR_cl[i] - min_wr
  }
  
  # Guardar en la lista con un nombre identificador
  nombre_lista <- paste0("N", curr$N, "_AQL", curr$AQL, "_a", curr$alpha)
  lista_resultados[[nombre_lista]] <- res_escenario
}