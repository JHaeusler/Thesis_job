# ============================================================================== #
# SCRIPT VISUALIZADOR: DENSIDADES A PRIORI ELICITADAS (MÉTODO DE COOK)
# Dependencia: Requiere que 'diccionario_resultados' exista en memoria 
#              (Generado por 'elicit_cook_2010.R' en la Fase 0)
# Autor: Juan Sebastián Haeusler
# ============================================================================== #

# --- 1. Verificación de Dependencias ---
if(!require(dplyr)) install.packages("dplyr")

# Si el diccionario no existe en memoria, llama al script de elicitación para crearlo
if(!exists("diccionario_resultados")) {
  cat("[!] El diccionario de prioris no se encontró. Ejecutando 'elicit_cook_2010.R'...\n")
  source("elicit_cook_2010.R", encoding = "UTF-8")
}

# --- 2. Configuración de Parámetros Gráficos ---
colores_prov <- c("forestgreen", "dodgerblue", "darkgoldenrod", "darkorange", "firebrick")
nombres_prov <- c("Excelente", "Bueno", "Regular", "Malo", "Muy Malo")

# Extraer las 12 combinaciones únicas de AQL y LTPD del diccionario
escenarios_unicos <- diccionario_resultados %>% select(AQL, LTPD) %>% distinct()

# Crear carpeta de imágenes si no existe
if(!dir.exists("imagenes")) dir.create("imagenes")

file_name <- "imagenes/Panel_Densidades_Cook.png"
cat(sprintf("\nGenerando panel gráfico maestro en: %s ...\n", file_name))

# Iniciar el dispositivo gráfico PNG (Alta resolución)
png(file_name, width = 1400, height = 1000, res = 120)

# Configurar el lienzo: 3 filas (AQLs) x 4 columnas (LTPDs)
# mar = márgenes de cada gráfico, oma = margen exterior para título y leyenda
par(mfrow = c(3, 4), mar = c(4, 4, 3, 1), oma = c(4, 1, 4, 1))

# Secuencia del eje X (Proporción de defectuosos p)
x_seq <- seq(0, 0.3, length.out = 500)
limite_Y_max <- 40 # Tope visual para las asíntotas (a < 1)

# --- 3. Ciclo de Graficación por Escenario ---
for (k in 1:nrow(escenarios_unicos)) {
  
  aql_cur <- escenarios_unicos$AQL[k]
  ltpd_cur <- escenarios_unicos$LTPD[k]
  
  # Filtrar los 5 perfiles para este AQL y LTPD específico
  params_esce <- diccionario_resultados %>% filter(AQL == aql_cur, LTPD == ltpd_cur)
  
  # Dibujar lienzo base de la subgráfica
  plot(NULL, xlim = c(0, 0.25), ylim = c(0, limite_Y_max), bty = "l",
       main = sprintf("AQL: %.2f | LTPD: %.2f", aql_cur, ltpd_cur),
       xlab = "Proporción de Defectos (p)", ylab = "Densidad",
       cex.main = 1.2, cex.lab = 0.9, las = 1)
  
  # Líneas de referencia para AQL y LTPD
  abline(v = aql_cur, col = "forestgreen", lty = 2, lwd = 1.5)
  abline(v = ltpd_cur, col = "firebrick", lty = 2, lwd = 1.5)
  
  # Graficar las 5 densidades de los expertos
  for (i in 1:5) {
    a_val <- params_esce$a[i]
    b_val <- params_esce$b[i]
    
    # Manejo de excepciones: Si el parámetro es NA (Ej: Excelente en AQL 0.01 y LTPD 0.10)
    if (is.na(a_val) || is.na(b_val)) next 
    
    # Calcular densidad
    y_seq <- dbeta(x_seq, shape1 = a_val, shape2 = b_val)
    
    # Control de infinitos en p=0 cuando a < 1
    y_seq[!is.finite(y_seq)] <- limite_Y_max
    y_seq[y_seq > limite_Y_max] <- limite_Y_max
    
    # Trazar la curva
    lines(x_seq, y_seq, col = colores_prov[i], lwd = 2.5)
  }
}

# --- 4. Título General y Leyenda Exterior ---
mtext("Diccionario de Prioris Elicitadas (Método de Cook, 2010)", 
      outer = TRUE, side = 3, cex = 1.8, font = 2, line = 1)

# Leyenda global en la parte inferior del panel
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", legend = nombres_prov, col = colores_prov, lwd = 3, 
       horiz = TRUE, bty = "n", cex = 1.2, xpd = TRUE, inset = c(0, 0.02),
       title = "Perfiles de Creencia del Experto")

dev.off() # Cierra y guarda la imagen

cat("¡Panel generado exitosamente! Revise la carpeta 'imagenes'.\n")