# Propuesta para diseño de planes de muestreo simple para atributos
# considerando el histórico de la proporción estimada de unidades
# defectuosas en la inspección de lotes

# Calidad del proveedor según el histórico
# Las prioris de calidad de los proveedores se definen por pares (PA en AQL, PA en LTPD)
# Excelente proveedor: PA(AQL)=0.95, PA(LTPD)=1-0.0001
# Buen proveedor: PA(AQL)=0.75, PA(LTPD)=1-0.05
# Regular: PA(AQL)=0.50, PA(LTPD)=1-0.10
# Malo: PA(AQL)=0.25, PA(LTPD)=1-0.25
# Muy malo: PA(AQL)=0.10, PA(LTPD)=1-0.40

# --- Carga de Librerías ---
# install.packages("pracma") # Si no lo tienes instalado
# install.packages("sjstats")

library(pracma)
library(sjstats)
library(AcceptanceSampling)
library(GA)

# --- Funciones de Riesgo Ponderado y Masa de Probabilidad ---

# Función para calcular la masa de probabilidad (densidad acumulada) en las regiones de riesgo
calc_prob_mass <- function(alpha_b, beta_b, AQL, LTPD) {
  P_Good <- pbeta(AQL, shape1 = alpha_b, shape2 = beta_b)
  P_Bad <- 1 - pbeta(LTPD, shape1 = alpha_b, shape2 = beta_b)
  return(c(P_Good = P_Good, P_Bad = P_Bad))
}

# Función principal para calcular el Riesgo Ponderado
calc_wr <- function(N_, n, c, alpha_b, beta_b, AQL_, LTPD_, k_p, k_c) {

  wrp_val <- k_p*(1-phyper(c, N_ * AQL_, N_ * (1 - AQL_), n))
  wrc_val <- k_c*phyper(c, N_ * LTPD_, N_ * (1 - LTPD_), n)
  wrt_val <- wrp_val + wrc_val

  return(c(WRP_val = wrp_val, WRC_val = wrc_val, WRT_val = wrt_val))
}

decode <- function(string){
  string <- gray2binary(string)
  n <- binary2decimal(string[1:l1])
  c <- min(n, binary2decimal(string[(l1 + 1):(l1 + l2)]))
  
  return(c(n,c))
}

fitness <- function(string){
  par <- decode(string)
  n <- par[1]
  c <- par[2]
  Pa_p <- phyper(c, Esce[esce, 1]*Esce[esce, 4], Esce[esce, 1]*(1 - Esce[esce, 4]), n)
  Pa_c <- phyper(c, Esce[esce, 1]*Esce[esce, 5], Esce[esce, 1]*(1 - Esce[esce, 5]), n)
  Loss <- (Pa_p - (1 - Esce[esce, 2]))^2 + (Pa_c - Esce[esce, 3])^2
  -Loss
}

# --- Variables Globales y Parámetros del Problema ---

N <- 1000 #c(200, 1000)        # Tamaño del lote
alpha <- c(0.01, 0.05)   # Riesgo del productor (para el plan clásico)
beta <- c(0.05, 0.10, 0.20)    # Riesgo del consumidor (para el plan clásico)
AQL <- c(0.01, 0.02, 0.05)     # Nivel de Calidad Aceptable
LTPD <- c(0.08, 0.10, 0.15, 0.20)    # Tolerancia de Porcentaje Defectuoso de Lote

Esce <- expand.grid(N = N, alpha = alpha, beta = beta, AQL = AQL, LTPD = LTPD)
Esce <- Esce[Esce$AQL < Esce$LTPD,]

# Vectores de probabilidades para los 5 escenarios
# p1 = Probabilidad de Aceptación en el AQL
# p2 = Probabilidad de No Aceptación (o Rechazo) en el LTPD
p1 <- c(0.95, 0.75, 0.50, 0.25, 0.10, 0.1, 0.3)
p2 <- c(1e-4, 0.05, 0.10, 0.25, 0.40, 0.89, 0.69) # 1 - PA(LTPD) -> Beta (Riesgo del Consumidor)

# Inicializar las listas maestras
lista_tablas_resultados <- list()
lista_parametros_beta   <- list()

for (esce in 1:dim(Esce)[1]) { # esce <- 1 + esce

  # 1. Determinar el Plan Clásico (basado en alpha y beta fijos)
  plan_clasic <- find.plan(PRP = c(Esce[esce, 4], 1 - Esce[esce, 2]),
               CRP = c(Esce[esce, 5], Esce[esce, 3]),
               N = Esce[esce, 1], type = "hypergeom")
  
  n_clasic <- plan_clasic$n
  c_clasic <- plan_clasic$c
  
  # n_ran <- 2:Esce[esce, 1]
  # c_ran <- 0:(max(n_ran) - 1)
  # 
  # b1 <- decimal2binary(max(n_ran)); l1 <- length(b1)
  # b2 <- decimal2binary(max(c_ran)); l2 <- length(b2)
  # 
  # plan_genetico <- ga(type = "binary", nBits = l1 + l2,
  #                     fitness = fitness, popSize = 200,
  #                     maxiter = 200, run = 100, seed = 060722)
  # 
  # plan_ga <- decode(plan_genetico@solution)
  # 
  # n_ga <- plan_ga[1]
  # c_ga <- plan_ga[2]
  
  # 2. Obtener los pares de (alpha_b, beta_b) para la distribución Beta
  alpha_beta_params <- data.frame(alpha_b = rep(NA, length(p1)), beta_b = rep(NA, length(p1)))
  
  for (i in 1:length(p1)) { # i <- 1 + i
  
  # find_beta encuentra los parámetros de la distribución Beta
  Shape <- find_beta(x1 = Esce[esce, 4], p1 = p1[i], x2 = Esce[esce, 5], p2 = 1 - p2[i]) 
  alpha_beta_params[i, "alpha_b"] <- Shape$shape1
  alpha_beta_params[i, "beta_b"] <- Shape$shape2
  }
  
  # 3. Inicializar la tabla de resultados
  resultados_riesgo <- data.frame(
  Escenario = 1:length(p1),
  Proveedor = c(c("Excelente", "Bueno", "Regular", "Malo", "Muy Malo"), paste("otro", 1:(abs(5-length(p1))))),
  n_clasic = n_clasic,
  c_clasic = c_clasic,
  
  n_ga = NA, #n_ga,
  c_ga = NA, #c_ga,
  
  P_Mass_Good_Beta = NA,                 # P(p < AQL | Beta)
  P_Mass_Bad_Beta = NA,                  # P(p > LTPD | Beta)
  
  WRP_clasic = NA,
  WRC_clasic = NA,
  WRT_clasic = NA,
  # Plan Óptimo (Minimiza TWR bajo Información)
  n_opt = NA,
  c_opt = NA,
  WRP_opt = NA,
  WRC_opt = NA,
  WRT_opt = NA,
  
  # Ganancias (Plan Clásico Naive vs. Plan Óptimo Informado)
  WRP_Ganancia = NA, 
  WRC_Ganancia = NA,
  WRT_Ganancia = NA
  )
  
  # 4. Iterar sobre los escenarios y calcular los riesgos
  
  for (j in 1:length(p1)) { # j <- 1 + j
    
    # Usar los valores calculados de alpha y beta específicos para el proveedor i
    alpha_b_val <- alpha_beta_params[j, "alpha_b"]
    beta_b_val <- alpha_beta_params[j, "beta_b"]
    
    # I. Calcular la masa de probabilidad para el Prior Beta específico (Informado)
    prob_mass_beta <- calc_prob_mass(alpha_b = alpha_b_val, beta_b = beta_b_val, Esce[esce, 4], Esce[esce, 5])
    resultados_riesgo[j, "P_Mass_Good_Beta"] <- prob_mass_beta["P_Good"]
    resultados_riesgo[j, "P_Mass_Bad_Beta"] <- prob_mass_beta["P_Bad"]
    
    k_p_ <- as.numeric(prob_mass_beta["P_Good"])
    k_c_ <- as.numeric((prob_mass_beta["P_Bad"]))
    
    # II. Calcular los riesgos para el PLAN CLÁSICO (Bajo el Prior Beta ESPECÍFICO)
    risks_clasic <- calc_wr(N_ = Esce[esce, 1], n = n_clasic, c = c_clasic, 
                  alpha_b = alpha_b_val, beta_b = beta_b_val, 
                  Esce[esce, 4], Esce[esce, 5], k_p = k_p_, k_c = k_c_)
    
    resultados_riesgo[j, "WRT_clasic"] <- risks_clasic["WRT_val"]
    resultados_riesgo[j, "WRP_clasic"] <- risks_clasic["WRP_val"]
    resultados_riesgo[j, "WRC_clasic"] <- risks_clasic["WRC_val"]
    
    
    # --- III. Búsqueda del Plan Óptimo (Minimizar TWR Puro) para el Escenario i ---
    
    n_opt_found <- NA
    c_opt_found <- NA
    des_WRT <- k_p_ * Esce[esce, 2] + k_c_ * Esce[esce, 3]
    cumple <- FALSE
    # Búsqueda exhaustiva
    for (n_ in 1:Esce[esce, 1]) { # n_ <- 1 + n_
      for (c_ in 0:(n_ - 1)) { # c_ <- 0 + 1 + c_
        
        risks_opt <- calc_wr(N_ = Esce[esce, 1], n = n_, c = c_, 
                       alpha_b = alpha_b_val, beta_b = beta_b_val, 
                       Esce[esce, 4], Esce[esce, 5], k_p = k_p_, k_c = k_c_)
        
        WRP_current <- risks_opt["WRP_val"]
        WRC_current <- risks_opt["WRC_val"]
        WRT_current <- risks_opt["WRT_val"]
        
        # Usando la lógica de tu código: buscar un plan 'mejor' que el plan Naive
        if (WRT_current <= des_WRT) {
          n_opt_found <- n_
          c_opt_found <- c_
          cumple <- TRUE
          break
        }
        }
      if(cumple) break
      } 
      
      
      # IV. Almacenar los resultados del plan óptimo (n_opt, c_opt)
      resultados_riesgo[j, "n_opt"] <- n_opt_found
      resultados_riesgo[j, "c_opt"] <- c_opt_found
      
      # Recalcular los riesgos para el análisis final
      if (!is.na(n_opt_found)) {
      risks_opt_final <- calc_wr(N_ = Esce[esce, 1], n = n_opt_found, c = c_opt_found, 
                       alpha_b = alpha_b_val, beta_b = beta_b_val, 
                       Esce[esce, 4], Esce[esce, 5], k_p = k_p_, k_c = k_c_) 
      
      resultados_riesgo[j, "WRT_opt"] <- risks_opt_final["WRT_val"]
      resultados_riesgo[j, "WRP_opt"] <- risks_opt_final["WRP_val"]
      resultados_riesgo[j, "WRC_opt"] <- risks_opt_final["WRC_val"]
      
      # V. ¡Cálculo de Ganancias Añadido!
      # Ganancia = Riesgo Naive - Riesgo Óptimo
      resultados_riesgo[j, "WRP_Ganancia"] <- resultados_riesgo[j, "WRP_clasic"] - resultados_riesgo[j, "WRP_opt"]
      resultados_riesgo[j, "WRC_Ganancia"] <- resultados_riesgo[j, "WRC_clasic"] - resultados_riesgo[j, "WRC_opt"]
      resultados_riesgo[j, "WRT_Ganancia"] <- resultados_riesgo[j, "WRT_clasic"] - resultados_riesgo[j, "WRT_opt"]
      
    
    } else {
      resultados_riesgo[j, "WRT_opt"] <- NA
      resultados_riesgo[j, "WRP_opt"] <- NA
      resultados_riesgo[j, "WRC_opt"] <- NA
      
      resultados_riesgo[j, "WRP_Ganancia"] <- NA
      resultados_riesgo[j, "WRC_Ganancia"] <- NA
      resultados_riesgo[j, "WRT_Ganancia"] <- NA
    }
  }
  
  
  # Guardar la tabla de resultados del escenario actual en la lista maestra
  lista_tablas_resultados[[esce]] <- resultados_riesgo
  
  # Guardar los parámetros de las densidades Beta calculadas para este escenario
  lista_parametros_beta[[esce]] <- alpha_beta_params

} 



# --- Configuración del Panel ---
# 12 ventanas (3 filas x 4 columnas) + espacio para leyenda
# Definimos los colores solicitados
colores_prov <- c("darkolivegreen4", "dodgerblue1", "gold2", "darkorange1", "tomato2")
nombres_prov <- c("Excelente", "Bueno", "Regular", "Malo", "Muy Malo")

# guardar archivo
file_name <- paste0("densidades_historico.png")
png(file_name, width = 1200, height = 600, res = 100)


# x11()
# Ajustar márgenes: mfrow para la cuadrícula, oma para la leyenda global
par(mfrow = c(4, 3), mar = c(3, 3, 2, 1), oma = c(4, 1, 2, 1))
cases <- seq(1, dim(Esce)[1], 6)


for (k in cases) { # k <- 1 + k
  
  # Extraer parámetros del escenario k
  params <- lista_parametros_beta[[k]][1:5,]
  
  # Configurar el área de dibujo para la ventana k
  # Usamos un xlim hasta 0.25 para cubrir el LTPD máximo común
  plot(NULL, xlim = c(0, 0.25), ylim = c(0, 45), 
       main = paste("AQL:", Esce$AQL[k], "| LTPD:", Esce$LTPD[k]),
       xlab = "", ylab = "", cex.main = 0.9)
  
  # Eje X de densidad
  x_seq <- seq(0, 0.3, length.out = 200)
  
  abline(v = c(Esce$AQL[k], Esce$LTPD[k]),
         col = "bisque3", lty = 2)
  
  # Graficar cada una de las 5 densidades
  for (l in 1:dim(params)[1]) { # l <- 1 + l
    lines(x_seq, dbeta(x_seq, shape1 = params$alpha_b[l],
                       shape2 = params$beta_b[l]), 
          col = colores_prov[l], lwd = 1.5)
  }
  
  # Agregar una rejilla sutil para facilitar la lectura
  #grid(col = "gray90")
}

# --- Leyenda Global ---
# Volvemos al espacio exterior para la leyenda
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend = nombres_prov, col = colores_prov, 
       lwd = 3, horiz = TRUE, bty = "n", inset = c(0, 0.002), cex = 1.2)

# Título general del panel
mtext("Distribución de Probabilidad del Histórico por Tipo de Proveedor", 
      side = 3, line = -2, outer = TRUE, font = 2, cex = 1.2)

dev.off()