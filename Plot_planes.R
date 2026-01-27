

# --- Configuración del Panel ---
# 12 ventanas (3 filas x 4 columnas) + espacio para leyenda
# Definimos los colores solicitados
colores_prov <- c("darkolivegreen4", "dodgerblue1", "gold2", "darkorange1", "tomato2",
                  "chartreuse2", "darkslategray")
nombres_prov <- c("Excelente", "Bueno", "Regular", "Malo", "Muy Malo","Naive", "OTB_est")

# guardar archivo
file_name <- paste0("planes.png")
png(file_name, width = 1200, height = 600, res = 100)


# x11()
# Ajustar márgenes: mfrow para la cuadrícula, oma para la leyenda global
par(mfrow = c(4, 3), mar = c(3, 3, 2, 1), oma = c(4, 1, 2, 1))
cases <- seq(1, dim(Esce)[1], 6)


for (k in cases) { # k <- 1 + k
  
  # Extraer parámetros del escenario k
  planes <- lista_tablas_resultados[[k]]#[1:5,]
  
  n_plans <- planes$n_opt
  c_plans <- planes$c_opt
  
  
  # Configurar el área de dibujo para la ventana k
  # Usamos un xlim hasta 0.25 para cubrir el LTPD máximo común
  plot(NULL, xlim = c(0, 0.25), ylim = c(0, 1), 
       main = paste("AQL:", Esce$AQL[k], "| LTPD:", Esce$LTPD[k]),
       xlab = "", ylab = "", cex.main = 0.9)
  
  # Eje X de densidad
  x_seq <- seq(0, 1, length.out = 200)
  
  abline(v = c(Esce$AQL[k], Esce$LTPD[k]),
         col = "bisque3", lty = 2)
  
  # Graficar cada una de las 5 densidades
  for (l in 1:dim(planes)[1]) { # l <- 1 + l
    Pa_planes <- phyper(c_plans[l], N*x_seq, N*(1 - x_seq), n_plans[l])
    lines(x_seq, Pa_planes, col = colores_prov[l], lwd = 1.5)
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
mtext("OC tipo A  para los planes encontrados según el histórico del Proveedor", 
      side = 3, line = -2, outer = TRUE, font = 2, cex = 1.2)

dev.off()