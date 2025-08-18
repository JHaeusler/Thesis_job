# Optimización Fernandez - Correa - Pericchi
# Balancing producer and consumer risks in optimal attribute testing: A unified Bayesian/Frequentist design
# 2020

# 1. Parámetros Iniciales y Pesos
# -----------------------------------------------------------------------------
N <- 200 # Tamaño máximo de muestra a considerar
p_0 <- 0.05 # Proporción de defectuosos aceptable (riesgo del productor)
p_1 <- 0.10 # Proporción de defectuosos inaceptable (riesgo del consumidor)
w_0 <- 0.2 # Peso asignado al Riesgo del Productor (Error Tipo I)
w_1 <- 1 - w_0 # Peso asignado al Riesgo del Consumidor (Error Tipo II)

# La variable 'prop' no se utiliza en el código actual, se puede omitir si no es necesaria para futuras adiciones.
# prop <- seq(0, 1, 0.001)

# 2. Generación de Combinaciones de Planes de Muestreo (n, c)
# -----------------------------------------------------------------------------
# Crear todas las combinaciones posibles de n (tamaño de muestra) y c (número de aceptación)
n.opt_all <- 1:N # n puede ir de 1 a N
c.opt_all <- 0:max(n.opt_all) # c puede ir de 0 hasta el n más grande
Esc <- expand.grid(n.opt = n.opt_all, c.opt = c.opt_all)

# Filtrar combinaciones donde n > c. Esto es un requisito lógico: no puedes aceptar más defectos
# de los que hay en la muestra, y c=n implica que siempre aceptas la muestra si todos son no defectuosos.
# Sin embargo, si c puede ser igual a n (por ejemplo, c=0 para n=0), la condición debería ser n >= c.
# Para el muestreo de aceptación típico, n > c es común. Si c=n se considera una aceptación trivial, se mantiene.
Esc <- Esc[Esc$n.opt > Esc$c.opt, ]

# 3. Cálculo de Riesgos y Razón de Verosimilitud para cada Plan
# -----------------------------------------------------------------------------
# Inicializar una matriz para almacenar los resultados para cada plan (n, c).
# Columna 1: n (tamaño de la muestra)
# Columna 2: c (número de aceptación)
# Columna 3: Riesgo total ponderado (w0*RiesgoProductor + w1*RiesgoConsumidor)
# Columna 4: Razón de Verosimilitud (Prob(c|p0) / Prob(c|p1))
wr_t <- matrix(NA,nrow = nrow(Esc), ncol = 4) # 4 columnas inicialmente, la 5ª se añade después

# Asignar n y c a las primeras columnas de wr_t
wr_t[,1] <- Esc$n.opt # Usar nombres de columna es más claro que Esc[,1]
wr_t[,2] <- Esc$c.opt # Usar nombres de columna es más claro que Esc[,2]

# Bucle para calcular los valores para cada combinación de (n, c)
for(j in 1:nrow(Esc)){
  # Calcular el Riesgo del Productor (Alpha): Probabilidad de RECHAZAR un lote BUENO (1 - P(aceptar | p0))
  alpha_riesgo <- 1 - pbinom(Esc$c.opt[j], Esc$n.opt[j], p_0)
  # Calcular el Riesgo del Consumidor (Beta): Probabilidad de ACEPTAR un lote MALO (P(aceptar | p1))
  beta_riesgo <- pbinom(Esc$c.opt[j], Esc$n.opt[j], p_1)
  
  # Columna 3: Suma Ponderada de Riesgos (Función de Costo Objetivo)
  wr_t[j,3] <- w_0 * alpha_riesgo + w_1 * beta_riesgo
  
  # Columna 4: Razón de Verosimilitud (Likelihood Ratio)
  # Probabilidad de observar 'c' bajo p0 / Probabilidad de observar 'c' bajo p1
  # Se usa 'dbinom' porque es la función de masa de probabilidad (PMF)
  # Se añade una pequeña constante al denominador para evitar divisiones por cero
  # si dbinom(Esc$c.opt[j], Esc$n.opt[j], p_1) fuera 0.
  prob_p0 <- dbinom(Esc$c.opt[j], Esc$n.opt[j], p_0)
  prob_p1 <- dbinom(Esc$c.opt[j], Esc$n.opt[j], p_1)
  wr_t[j,4] <- prob_p0 / (prob_p1 + .Machine$double.eps) # Añadir epsilon para estabilidad numérica
}

# Convertir wr_t a un data.frame para facilitar el manejo con nombres de columna
wr_t <- as.data.frame(wr_t)
colnames(wr_t) <- c("n", "c", "RiesgoPonderado", "RazonVerosimilitud")

# 4. Filtrado y Selección del Plan Óptimo
# -----------------------------------------------------------------------------
# Calcular la razón 'r' de los pesos (Umbral para la Razón de Verosimilitud)
r <- w_1 / w_0 # En tu caso, 0.8 / 0.2 = 4

# Filtrar los planes de muestreo donde la Razón de Verosimilitud es menor que 'r'
# Este filtro se basa en un criterio de optimización (ej. Lema de Neyman-Pearson)
wr_t_R_ast <- wr_t[wr_t$RazonVerosimilitud < r, ]

# Encuentra la(s) fila(s) en el subconjunto filtrado donde el Riesgo Ponderado es mínimo.
# Esto identifica el plan(es) (n, c) que minimiza la función objetivo dentro de los que cumplen el criterio.
indices_optimos <- which(wr_t_R_ast$RiesgoPonderado == min(wr_t_R_ast$RiesgoPonderado))

# Obtener el(los) plan(es) óptimo(s)
plan_optimo <- wr_t_R_ast[indices_optimos, ]

# Imprimir el resultado principal
print("---")
print("Plan(es) de muestreo óptimo(s) encontrado(s):")
print(plan_optimo)
print("---")

# 5. Cálculos Adicionales (g_n y c_n) - Contexto y Ajustes
# -----------------------------------------------------------------------------
# La fórmula para 'g_n' se relaciona con la línea de discriminación en gráficos (n, c)
# o en muestreo secuencial. Es una función de los parámetros de probabilidad y los pesos.
g_n <- (log((1-p_0)/(1-p_1)) - log(r))/(log((p_1*(1-p_0))/(p_0*(1-p_1))))

# Añadir una 5ª columna a wr_t_R_ast con un 'c' teórico calculado para cada 'n' filtrado.
# Esta columna te da una idea de cuál sería el 'c' óptimo si siguieras una aproximación continua.
wr_t_R_ast$c_teorico_gn <- floor(wr_t_R_ast$n * g_n)

# Para 'c_n', si tu objetivo es obtener un 'c' general basado en esta fórmula
# y luego usarlo, es importante saber si debe aplicarse a todos los planes
# filtrados o solo al 'plan_optimo' encontrado.
# Asumiendo que quieres el 'c_n' para el *primer* plan óptimo si hay varios, o el único si es uno:
# También aseguramos que c_n no sea menor que 0, ya que c representa un conteo de defectos.
if (nrow(plan_optimo) > 0) {
  n_optimo_single <- plan_optimo$n[1] # Tomar el 'n' del primer plan óptimo
  c_n <- floor(n_optimo_single * g_n) # Calcular c_n basado en el 'n' óptimo
  c_n <- max(0, c_n) # c no puede ser negativo, por lo que lo limitamos a 0 como mínimo
  
  print(paste("Valor de g_n calculado:", round(g_n, 4)))
  print(paste("Valor de c_n (c teórico para el n óptimo) calculado:", c_n))
  
  # Opcional: Comparar el 'c' real del plan óptimo con el 'c' teórico
  print(paste("c real del plan óptimo:", plan_optimo$c[1]))
  if (plan_optimo$c[1] == c_n) {
    print("El c del plan óptimo coincide con el c teórico calculado por g_n.")
  } else {
    print("El c del plan óptimo NO coincide con el c teórico calculado por g_n (debido a redondeo, etc.).")
  }
} else {
  print("No se encontraron planes óptimos que cumplan los criterios.")
}
