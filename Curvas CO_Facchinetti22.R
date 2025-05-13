# parámtros dados en # Facchinetti, S. 2022
N <- c(1000)
n <- c(50, 100, 150, 200, 250) 
c <- c(0, 1, 2, 3, 4)
n1 <- c(100) # n fijo
c1 <- c(3)  # c fijo
prop <- seq(0, 1, 0.01)

cols <- c("aquamarine4", "aquamarine3", "aquamarine2", "cadetblue2", "chartreuse4",
          "chartreuse3", "chartreuse2", "darkolivegreen3", "darkorange3", "dodgerblue2")

x11()
plot(prop, phyper(c[1], N[1]*prop, N[1]*(1 - prop), n1),
     main = paste("CCO (Tipo A) de un lote de tamaño N=", N[1], ", tamaño de muestra \n",
                  "n=",n1, "y diferentes valores para c"), type = "l", xlim = c(0, 0.22),
     col = cols[1], ylab = "Pa", xlab = "p")

c_rest <- c[-1]
for (i in 1:length(c_rest)) {
  lines(prop, phyper(c_rest[i], N[1]*prop,N[1]*(1 - prop), n1), lty= 1,
        col = cols[i + 1])
}






for(t_N in 1:length(N)){
  
  for (j in 1:length(n)) {
    x11()
    plot(prop, phyper(c[1], N[t_N ]*prop, N[t_N ]*(1 - prop), n[j]),
         main = paste("CCO (Tipo A) de un lote de tamaño N=", N[t_N ], ", tamaño de muestra \n",
                      "n=",n[j], "y diferentes valores para c"), type = "l",
         col = cols[1], ylab = "Pa", xlab = "p")
    index <- which(c < n[j])
    c_rest <- c[index][-1]
    for (i in 1:length(c_rest)) {
      lines(prop, phyper(c_rest[i], N[t_N ]*prop,N[t_N ]*(1 - prop), n[j]), lty= 1,
            col = cols[i + 1])
    }
    legend(x = "topright", legend = paste("c=",c[index]), col = cols[index], 
           title = "Criterio", lty= 1, box.lty = 0, bg = NA)
    # abline(v=c(0.05, 0.10), lty = 3, col = c("burlywood4"))
    segments(x0= 0.05, y0= 0, x1= 0.05, y1= 1.5, lty = 3, col = c("burlywood4"))
    segments(x0= 0.10, y0= 0, x1= 0.10, y1= 1.5, lty = 3, col = c("burlywood4"))
    abline(v= 0, h= 0, lty = 2, col = c("darkgray"))
    text(x= 0.05, y= -0.02, label= "0.05", cex = 0.6)
    text(x= 0.10, y= -0.02, label= "0.10", cex = 0.6)
  }
  
  for (j in 1:length(c1)) {
    x11()
    index <- which(n1 >  c1[j])
    if(length(index) > 0){ 
      n_rest <- n1[index]
      plot(prop, phyper(c1[j], N[t_N ]*prop, N[t_N ]*(1- prop), n_rest[1]),
           main = paste("CCO (Tipo A) de un lote de tamaño N=", N[t_N ], ", variando el \n tamaño de 
       la muestra para c=", c1[j]), type = "l",
           col = cols[1], ylab = "Pa", xlab = "p")
      n_rest2 <- n_rest[-1] 
      if(length(n_rest2)>0){
        for (i in 1:length(n_rest2)) {
          lines(prop, phyper(c1[j], N[t_N ]*prop,N[t_N ]*(1 - prop), n_rest2[i]), lty= 1,
                col = cols[i + 1])
        }
      }
      legend(x = "topright", legend = paste("n=",n_rest), col = cols[1:length(n_rest)],
             title = "tamaño muestra", lty= 1, box.lty = 0, bg = NA)
      segments(x0= 0.05, y0= 0, x1= 0.05, y1= 1.5, lty = 3, col = c("burlywood4"))
      segments(x0= 0.10, y0= 0, x1= 0.10, y1= 1.5, lty = 3, col = c("burlywood4"))
      abline(v= 0, h= 0, lty = 2, col = c("darkgray"))
      text(x= 0.05, y= -0.02, label= "0.05", cex = 0.6)
      text(x= 0.10, y= -0.02, label= "0.10", cex = 0.6)
    }
  }
}

x11()
plot(prop, phyper(2, 500*prop, 500*(1- prop), 200),
     main = paste("CCO (Tipo A) variando tamaño de lote con tamaño de \n 
       muestra n=", 200,"para c=", 10), type = "l",
     col = cols[1], ylab = "Pa", xlab = "p")
lines(prop, phyper(2, 5000*prop,5000*(1 - prop), 200), lty= 1,
      col = cols[2])
lines(prop, phyper(2, 50000*prop,50000*(1 - prop), 200), lty= 1,
      col = cols[3])
legend(x = "topright", legend = paste("n=",N), col = cols[1:length(N)],
       title = "tamaño muestra", lty= 1, box.lty = 0, bg = NA)
segments(x0= 0.05, y0= 0, x1= 0.05, y1= 1.5, lty = 3, col = c("burlywood4"))
segments(x0= 0.10, y0= 0, x1= 0.10, y1= 1.5, lty = 3, col = c("burlywood4"))
abline(v= 0, h= 0, lty = 2, col = c("darkgray"))
text(x= 0.05, y= -0.02, label= "0.05", cex = 0.6)
text(x= 0.10, y= -0.02, label= "0.10", cex = 0.6)
