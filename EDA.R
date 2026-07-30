# ==============================================================================
# SCRIPT CORREGIDO: EDA CON 4 ESCENARIOS (PROV 1, 2, 3A, 3B)
# ==============================================================================

# --- 1. Preparación de vectores ---
# Asegúrate de haber corrido el Panel de Control Maestro antes
idx_A <- which.min(Esce_Global$n_clasico)
idx_B <- which.max(Esce_Global$n_clasico)

lista_datos <- list(
  "Prov. 1" = na.omit(vec_data1),
  "Prov. 2" = na.omit(vec_data2),
  "Prov. 3A (n min)" = na.omit(Matriz_Data3[idx_A, ]),
  "Prov. 3B (n max)" = na.omit(Matriz_Data3[idx_B, ])
)

# --- 2. Tabla de Estadísticos Descriptivos ---
# --- Función Corregida ---
resumen_stats_final <- function(v, nombre, n_plan) {
  # Aseguramos que el paquete esté cargado para esta función
  if(!require(moments)) install.packages("moments")
  
  v_clean <- na.omit(v)
  q <- quantile(v_clean, probs = c(0.25, 0.5, 0.75))
  
  data.frame(
    ID = nombre,
    n_inspeccion = n_plan,
    Media = round(mean(v_clean), 4),
    SD = round(sd(v_clean), 4),
    CV_pct = round((sd(v_clean)/mean(v_clean))*100, 2),
    Q1 = round(q[1], 4),
    Q2 = round(q[2], 4),
    Q3 = round(q[3], 4),
    # Usamos el prefijo del paquete para evitar el error "function not found"
    Asim = round(moments::skewness(v_clean), 4)
  )
}

# n_plan para Prov 1 y 2 son ilustrativos o basados en sus papers, 
# para 3A y 3B vienen de Esce_Global
tabla_resultados <- rbind(
  resumen_stats_final(lista_datos[[1]], "Prov. 1", 27),
  resumen_stats_final(lista_datos[[2]], "Prov. 2", 30),
  resumen_stats_final(lista_datos[[3]], paste0("Prov. 3A (Esc ", idx_A, ")"), Esce_Global$n_clasico[idx_A]),
  resumen_stats_final(lista_datos[[4]], paste0("Prov. 3B (Esc ", idx_B, ")"), Esce_Global$n_clasico[idx_B])
)
print(tabla_resultados)

# --- 3. Gráfica de Densidades (Figura 5.1) ---
png("imagenes/EDA_Densidades_Empiricas.png", width = 900, height = 600, res = 100)
colores <- c("dodgerblue", "firebrick", "forestgreen", "purple")
colores_light <- c(rgb(0.1,0.5,0.8,0.3), rgb(0.8,0.2,0.2,0.3), rgb(0.2,0.8,0.2,0.3), rgb(0.5,0,0.5,0.3))

plot(NULL, xlim=c(0, 0.25), ylim=c(0, 35), bty="l", 
     main="Comparativa de Densidades Empíricas (4 Escenarios)",
     xlab="Proporción de Defectos (p)", ylab="Densidad",
     cex.main = 1.5,   
     cex.lab = 1.5,    
     cex.axis = 1.5)

for(i in 1:4) {
  d <- density(lista_datos[[i]], from=0, to=1)
  polygon(c(min(d$x), d$x, max(d$x)), c(0, d$y, 0), col=colores_light[i], border=colores[i], lwd=2)
}
abline(v=c(0.05, 0.10), col="darkgray", lty=2, lwd=2)
legend("topright", legend=names(lista_datos), fill=colores_light, bty="n", cex = 1.3)
text(x = 0.04, y = 35, labels = "AQL", cex = 1.3, col = "darkgray")
text(x = 0.11, y = 35, labels = "LTPD", cex = 1.3, col = "darkgray")

dev.off()

# --- 4. Boxplots (Figura 5.2) ---
png("imagenes/EDA_Boxplots.png", width = 1100, height = 600, res = 100)
boxplot(lista_datos, col=colores_light, border=colores,
        main="Dispersión de la Calidad por Fuente", 
        ylab="Proporción (p)", las=1,
        cex.main = 1.5,   
        cex.lab = 1.5,    
        cex.axis = 1.5)
abline(h=c(0.05, 0.10), col="darkgray", lty=2)
text(x = 0.5, y = 0.055, labels = "AQL", cex = 1.3, col = "darkgray")
text(x = 0.5, y = 0.105, labels = "LTPD", cex = 1.3, col = "darkgray")

dev.off()