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
x11()
plot(prop, phyper(c[1], N[1]*prop, N[1]*(1 - prop), n1), lty = lty_a[1],
     main = paste("CCO (Tipo A) de un lote de tamaño N=", N[1], "\n",
     "n= ", n1, " y diferentes valores para c"), type = "l", xlim = c(0, 0.155),
     col = cols_a[1], ylab = "Pa", xlab = "p")

c_rest <- c[-1]

for (i in 1:length(c_rest)) {
  lines(prop, phyper(c_rest[i], N[1]*prop,N[1]*(1 - prop), n1), lty = lty_a[i+1],
        col = cols_a[i + 1])
}

legend(x = "topright", legend = paste("c=",c), col = cols_a, 
       title = "Criterio", lty= c(2, 1, 1, 3, 1), box.lty = 0, bg = NA)
# abline(v=c(0.05, 0.10), lty = 3, col = c("burlywood4"))
segments(x0= 0.05, y0= 0, x1= 0.05, y1= 1.5, lty = 3, col = c("burlywood4"))
segments(x0= 0.10, y0= 0, x1= 0.10, y1= 1.5, lty = 3, col = c("burlywood4"))
abline(v= 0, h= c(0,1), lty = 1, col = c("azure3"))
text(x= 0.05, y= -0.02, label= "AQL", cex = 0.6)
text(x= 0.10, y= -0.02, label= "LTPD", cex = 0.6)
text(x=0.042, y = 0.6, label = "c = 3", col = "darkorchid")


x11()
plot(prop, phyper(c1, N[1]*prop, N[1]*(1 - prop), n[1]), lty = lty_b[1],
     main = paste(paste("CCO (Tipo A) de un lote de tamaño N=", N[1], "\n",
                      "c=", c1, " y variando el tamaño de la muestra n")), type = "l", xlim = c(0, 0.155),
     col = cols_b[1], ylab = "Pa", xlab = "p")

n_rest <- n[-1]
for (i in 1:length(n_rest)) {
  lines(prop, phyper(c1, N[1]*prop,N[1]*(1 - prop), n_rest[i]), lty= lty_b[i+1],
        col = cols_b[i + 1])
}

legend(x = "topright", legend = paste("n=",n), col = cols_b, 
       title = "tamaño de muestra", lty= lty_b, box.lty = 0, bg = NA)
# abline(v=c(0.05, 0.10), lty = 3, col = c("burlywood4"))
segments(x0= 0.05, y0= 0, x1= 0.05, y1= 1.5, lty = 3, col = c("burlywood4"))
segments(x0= 0.10, y0= 0, x1= 0.10, y1= 1.5, lty = 3, col = c("burlywood4"))
abline(v= 0, h= c(0,1), lty = 1, col = c("azure3"))
text(x= 0.05, y= -0.02, label= "AQL", cex = 0.6)
text(x= 0.10, y= -0.02, label= "LTPD", cex = 0.6)
text(x=0.042, y = 0.6, label = "n = 100", col = "darkorchid")
