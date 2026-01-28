setwd("D:/Github/Thesis_job/Thesis_job")

N <- 1000
prop <- seq(0, 1, 0.01)

n <- 1:8 
# Definición de colores con transparencia
violeta_trans <- rgb(0.5, 0, 0.5, alpha = 0.3) 
verde_trans   <- rgb(0, 0.5, 0, alpha = 0.3)   

for(n_ in 1:length(n)){
    
    # --- PASO 1: Definir nombre y abrir archivo ---
    # Creamos un nombre descriptivo: Ej: "CCO_n1_c0.png"
    file_name <- paste0("OC_n", n_,".png")
    png(file_name, width = 800, height = 600, res = 100) 
    
  for (c_ in 0:(n_ - 1)) {

    
    # Mantenemos tu lógica de plot
    plot(NULL, xlim = c(0, 1), ylim = c(0, 1),
         main = paste("OC (Tipo A) de un lote de tamaño N=", N, 
                      ", tamaño de muestra \n","n=", n_,
                      "y diferentes valores para c"), type = "l",
         ylab = "Pa", xlab = "p")
    
    for (i in 0:c_) {
      Pa <- phyper(i, N*prop, N*(1 - prop), n_)
      lines(prop, Pa, lty = 1, col = "lightskyblue4")
      
      if (i == 0) {
        # Zona ALPHA
        p_alpha <- seq(0, 0.05, 0.001)
        pa_alpha <- phyper(0, N*p_alpha, N*(1 - p_alpha), n_)
        polygon(c(p_alpha, rev(p_alpha)), 
                c(pa_alpha, rep(1, length(p_alpha))), 
                col = violeta_trans, border = NA)
        
        # Zona BETA
        p_beta <- seq(0.10, 1, 0.001)
        pa_beta <- phyper(0, N*p_beta, N*(1 - p_beta), n_)
        polygon(c(p_beta, rev(p_beta)), 
                c(pa_beta, rep(0, length(p_beta))), 
                col = verde_trans, border = NA)
      }
      
      x_pos <- 0.15 + 0.10*i
      y_pos <- phyper(i, N*x_pos, N*(1 - x_pos), n_)
      
      text(x = x_pos, y = y_pos + 0.02, 
           labels = paste("c =", i), 
           pos = 4, cex = 0.8, col = "lightskyblue4", font = 2)
    }
    
    # Elementos fijos del gráfico
    segments(x0= 0.05, y0= 0, x1= 0.05, y1= 1, lty = 3, col = "burlywood4")
    segments(x0= 0.10, y0= 0, x1= 0.10, y1= 1, lty = 3, col = "burlywood4")
    abline(v= 0, h= c(0, 1), lty = 2, col = "darkgray")
    
    # Etiquetas de riesgo y parámetros
    text(x= 0.05, y= -0.02, label= "AQL", cex = 0.6)
    text(x= 0.10, y= -0.02, label= "LTPD", cex = 0.6)
    text(x= 0.05, y= 1.02, label= "0.05", cex = 0.6)
    text(x= 0.10, y= 1.02, label= "0.10", cex = 0.6)
    text(x= 0.07, y= 0.97 - 0.02*c_, label= expression(alpha), cex = 1, col = rgb(0.5, 0, 0.5))
    text(x= 0.08, y= 0.4 - 0.02*c_, label= expression(beta), cex = 1, col = rgb(0, 0.5, 0))
    
    # --- PASO 2: Cerrar y guardar el archivo ---
    
  }
    dev.off() 
}

# ----------------- grafico convencional de dos puntos --------####
library(pBrackets)
library(AcceptanceSampling)

alpha <- 0.05
beta <- 0.10

AQL <- 0.05
LTPD <- 0.10

plan <- find.plan(PRP = c(AQL, 1 - alpha), CRP = c(LTPD, beta),
                  N = N, type = "hypergeom")
n_kir <- plan$n
c_kir <- plan$c

# guardar archivo
file_name <- paste0("OC_kiermeier.png")
png(file_name, width = 800, height = 600, res = 100)

plot(NULL, xlim = c(0, 0.2), ylim = c(0, 1.01), bty = "l", type = "n",
     xaxs = "i", yaxs = "i", ylab = expression(P[a]), xlab = expression(p),
     main = paste("OC (Tipo A) de un lote de tamaño N=", N, 
                  ", tamaño de muestra \n","n=", n_kir,
                  "y criterio de aceptación c=", c_kir))
lines(c(0,AQL), rep(1 - alpha, 2), lty = 2, col = "gray")
lines(rep(AQL,2), c(1 - alpha, 0), lty = 2, col = "gray")
lines(c(0,LTPD), rep(beta,2), lty = 2, col = "gray")
lines(rep(LTPD, 2), c(beta,0), lty = 2, col = rgb(0, 0.5, 0, alpha = 0.3), lwd = 3)
lines(prop, phyper(c_kir, N*prop, N*(1-prop), n_kir),
      lty = 1, col = "lightskyblue4", lwd = 2)
points(c(AQL, LTPD), c(1 - alpha, beta), pch = 19)

# 1. Llave para ALPHA (Vertical, al lado del eje Y o cerca del AQL)
# brackets(x1, y1, x2, y2, h, ticks, curvature, type, col)
brackets(x1 = AQL, y1 = 1, 
         x2 = AQL, y2 = 1 - alpha, 
         h = 0.005, # Curvatura/ancho de la llave
         ticks = 0.5, curvature = 0.5, type = 1, 
         col = rgb(0.5, 0, 0.5), lwd = 2)

# 2. Llave para BETA (Vertical, sobre la línea del LTPD)
brackets(x1 = LTPD, y1 = beta, 
         x2 = LTPD, y2 = 0, 
         h = - 0.005, 
         ticks = 0.5, curvature = 0.5, type = 1, 
         col = rgb(0, 0.5, 0), lwd = 2)

text(AQL - 0.025, 1 - alpha - 0.04, cex = 0.8,
     labels = expression(paste("(", AQL, ", ", 1 - alpha, ")")), pos = 4)
     
text(LTPD, beta + 0.02, cex = 0.8,
     labels = expression(paste("(", LTPD, ", ", beta, ")")), pos = 4)
text(x= AQL + 0.008, y= 0.97, label= expression(alpha), cex = 1, col = rgb(0.5, 0, 0.5))
text(x= LTPD - 0.008, y= 0.03, label= expression(beta), cex = 1, col = rgb(0, 0.5, 0))
segments(x0= AQL, y0= 0.95, x1= AQL, y1= 1, lty = 3, col = rgb(0.5, 0, 0.5), lwd = 3)
segments(x0= 0, y0= 1, x1= AQL, y1= 1, lty = 2, col = "gray")

dev.off()