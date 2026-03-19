rm(list=ls())
setwd("G:/Mi unidad/Trabajo de Investigación - JSH/TG_Maestría_Juan_Sebastián_Haeusler_Final_V_4_0/")
if(!dir.exists("imagenes")) dir.create("imagenes")
# parámtros dados en # Facchinetti, S. 2022
N <- c(1000)
n <- c(50, 100, 150, 200, 250) 
c <- c(0, 1, 2, 3, 4)
n1 <- c(100) # n fijo
c1 <- c(3)  # c fijo
prop <- seq(0, 1, 0.01)

cols_a <- c("aquamarine4", "aquamarine3", "cadetblue2", "darkorchid", "chartreuse4")
cols_b <- c("aquamarine4", "darkorchid", "cadetblue2", "darkolivegreen3", "chartreuse4")
# darkorchid
lty_a <- c(2, 1, 1, 3, 1)
lty_b <- c(2, 3, 1, 1, 1)

# guardar archivo
file_name <- paste0("imagenes/OC_Facchinetti_n.png")
png(file_name, width = 800, height = 600, res = 100)

plot(prop, phyper(c[1], N[1]*prop, N[1]*(1 - prop), n1), lty = lty_a[1],
     main = paste("CO (Tipo A) de un lote de tamaño N=", N[1], "\n",
     "n= ", n1, " y diferentes valores para c"), type = "l", xlim = c(0, 0.155),
     col = cols_a[1], ylab = "Pa", xlab = "p", lwd = 2,
     cex.main = 1.5,   
     cex.lab = 1.5,    
     cex.axis = 1.5)

c_rest <- c[-1]

for (i in 1:length(c_rest)) {
  lines(prop, phyper(c_rest[i], N[1]*prop,N[1]*(1 - prop), n1), lty = lty_a[i+1],
        col = cols_a[i + 1], lwd = 2)
}

legend(x = "topright", legend = paste("c=",c), col = cols_a, lwd = 2, cex = 1.5, 
       title = "Criterio", lty= c(2, 1, 1, 3, 1), box.lty = 0, bg = NA)
# abline(v=c(0.05, 0.10), lty = 3, col = c("burlywood4"))
segments(x0= 0.05, y0= 0, x1= 0.05, y1= 1.5, lty = 3, col = c("burlywood4"))
segments(x0= 0.10, y0= 0, x1= 0.10, y1= 1.5, lty = 3, col = c("burlywood4"))
abline(v= 0, h= c(0,1), lty = 1, col = c("azure3"))
text(x= 0.05, y= -0.02, label= "AQL", cex = 1.1)
text(x= 0.10, y= -0.02, label= "LTPD", cex = 1.1)
text(x=0.042, y = 0.6, label = "c = 3", col = "darkorchid", cex = 1.5)

dev.off()

# guardar archivo
file_name <- paste0("imagenes/OC_Facchinetti_c.png")
png(file_name, width = 800, height = 600, res = 100)
plot(prop, phyper(c1, N[1]*prop, N[1]*(1 - prop), n[1]), lty = lty_b[1],
     main = paste(paste("CO (Tipo A) de un lote de tamaño N=", N[1], "\n",
                      "c=", c1, " y variando el tamaño de la muestra n")), type = "l", xlim = c(0, 0.155),
     col = cols_b[1], ylab = "Pa", xlab = "p", lwd = 2,
     cex.main = 1.5,   
     cex.lab = 1.5,    
     cex.axis = 1.5)

n_rest <- n[-1]
for (i in 1:length(n_rest)) {
  lines(prop, phyper(c1, N[1]*prop,N[1]*(1 - prop), n_rest[i]), lty= lty_b[i+1],
        col = cols_b[i + 1], lwd = 2)
}

legend(x = "topright", legend = paste("n=",n), col = cols_b, lwd = 2, cex = 1.5, 
       title = "tamaño de muestra", lty= lty_b, box.lty = 0, bg = NA)
# abline(v=c(0.05, 0.10), lty = 3, col = c("burlywood4"))
segments(x0= 0.05, y0= 0, x1= 0.05, y1= 1.5, lty = 3, col = c("burlywood4"))
segments(x0= 0.10, y0= 0, x1= 0.10, y1= 1.5, lty = 3, col = c("burlywood4"))
abline(v= 0, h= c(0,1), lty = 1, col = c("azure3"))
text(x= 0.05, y= -0.02, label= "AQL", cex = 1.1)
text(x= 0.10, y= -0.02, label= "LTPD", cex = 1.1)
text(x=0.042, y = 0.6, label = "n = 100", col = "darkorchid", cex = 1.5)

dev.off()

# Curva CO escalón ideal

# Definimos el punto de corte (ejemplo: AQL = 0.05)

# La probabilidad de aceptación es 1 si p <= AQL, y 0 si p > AQL

AQL <- 0.05; LTPD <- 0.10
pa_ideal <- ifelse(prop < AQL, 1, 0)

file_name <- "imagenes/OC_ideal.png"
png(file_name, width = 800, height = 600, res = 100)

# Usamos type="s" para que R dibuje una escalera (staircase) 
# que representa la discontinuidad del plan ideal
plot(prop, pa_ideal, 
     type = "s", 
     main = "Curva CO Ideal (Discriminación Perfecta)", 
     xlim = c(0, 0.20), 
     ylim = c(0, 1.1), 
     col = "aquamarine4", 
     ylab = "Pa", 
     xlab = "p", 
     lwd = 2,
     bty = "l",
     cex.main = 1.5,   
     cex.lab = 1.5,    
     cex.axis = 1.5)

# 1. Punto de quiebre (opcional, ayuda a la vista)
points(AQL, 1, pch = 21, bg = "red", col = "aquamarine4")

# 2. Etiqueta AQL: 
# x = AQL, y = 1.05 (justo arriba del punto)
# adj = c(0.5, 0) centra horizontalmente y apoya la base del texto en la coordenada
text(x = AQL, y = 1.05, cex = 1.5, 
     labels = bquote(AQL == .(AQL)), 
     col = "aquamarine4", font = 2, adj = c(0.5, 0))

dev.off()
