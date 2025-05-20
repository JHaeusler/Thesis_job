# Trabajo Investigación ----

lista_de_paquetes <- c("pracma", "GA", "AcceptanceSampling",
                       "scatterplot3d", "tidyverse", "reshape") # Reemplaza con tus paquetes

for (paquete in lista_de_paquetes) {
  if(!require(paquete, character.only = TRUE)){
    install.packages(paquete)
    require(paquete, character.only = TRUE)
  }
}

# semilla
set.seed(123)

#---- Parámetros de entrada ----
N <- 1000; AQL <- 0.05; LTPD <- 0.10; alpha_des <- 0.05; beta_des <- 0.10
par1 <- 0; par2 <- 0.05
dens_n <- parse(text = paste(paste(
  paste("dunif(p",sep=""),
  par1, par2, sep=","),")", sep = "")) 

par11 <- 0.1; par12 <- 1
dens_n1 <- parse(text = paste(paste(
  paste("dunif(p",sep=""),
  par11, par12, sep=","),")", sep = "")) 

par21 <- 0; par22 <- 1
dens_n2 <- parse(text = paste(paste(
  paste("dunif(p",sep=""),
  par1, par2, sep=","),")", sep = "")) 

Pa <- function(n, c, p){
  Pa = phyper(c, N*p, N*(1-p), n) # COO tipo A
  return(Pa)
}

# riesgo ponderado prod
wr_p <- function(n, c, densidad, p){
  f.p = eval(densidad)
  prod_wr <- (1-Pa(n, c, p))*f.p
  return(prod_wr)
}

# riesgo ponderado cons
wr_c <- function(n, c, densidad, p){
  f.p = eval(densidad)
  cons_wr <- Pa(n, c, p)*f.p
  return(cons_wr)
}

calc_wr <- function(n, c, dens_n, N, AQL, LTPD){
  # Riesgo Ponderado para algunas distribuciones dadas
  wr_dist <- (integral(f = function(p) wr_p(n, c, dens_n, p),
                       xmin = 0, xmax = AQL, method = "Kron") +
                integral(f = function(p) wr_c(n, c, dens_n, p),
                         xmin = LTPD, xmax = 1, method = "Kron"))
  return(wr_dist)
} 



n.opt=1:300;c.opt=0:300; Esc=expand.grid(c.opt=c.opt,n.opt=n.opt);rm(n.opt,c.opt);Esc=Esc[Esc$n.opt>Esc$c.opt,]
# View(Esc);dim(Esc)

for(i in 1:nrow(Esc)){

  Esc[i,3]=calc_wr(Esc[i,2],Esc[i,1],dens_n,N,AQL,LTPD)
}

ggplot(Esc, aes(x = n.opt, y = c.opt, fill = V3)) +
  geom_tile() +
  geom_tile(color = "white",
            lwd = 0.0025,
            linetype = 2) +
  coord_fixed()




# Ajuste de la densidad a trabajar

# densidad normal
quant <- 0.9973 # cuantil

# par1 <- AQL/(3.5 + qnorm(quant))  # sigma
# par2<- AQL - qnorm(quant)*par1 # mu



# densidad uniforme
dens_u <- parse(text = paste(paste(
  paste("dunif(p",sep=""),
  0, 1, sep=","),")", sep = ""))

# decode para obtener el n y p propuesto por el GA
decode <- function(string){
  string <- gray2binary(string)
  n <- min(max(Ran.n), binary2decimal(string[1:l1]))
  c <- min(n, binary2decimal(string[(l1+1):(l1+l2)]))
  return(c(n, c))
}

# funciones auxiliares:
# Pa calcula la dist acum de hyper.


# espacios de búsqueda para n y c
Ran.n <- 2:ceiling(N*0.8)
Ran.c <- 0:ceiling(max(Ran.n)*0.8)

b1 <- decimal2binary(max(Ran.n)); l1 <- length(b1)
b2 <- decimal2binary(max(Ran.c)); l2 <- length(b2)

# solución mediante AcceptanceSampling (clásico)
plan <- find.plan(PRP = c(AQL, 1 - alpha_des),
                  CRP = c(LTPD, beta_des),
                  N=N, type = "hypergeom") ; plan$n; plan$c





# fitness function
fit_ <- function(string, dens_n, dens_u, N){
  par <- decode(string)
  n <- par[1]; c <- par[2]
  
  # Riesgo Ponderado para algunas distribuciones dadas
  wr_dist <- (integral(f = function(p) wr_p(n, c, dens_n, p),
                       xmin = 0, xmax = AQL, method = "Kron") +
                integral(f = function(p) wr_c(n, c, dens_n, p),
                         xmin = LTPD, xmax = 1, method = "Kron"))
  
  # Riesgo ponderado, asignando principio de equiporbabilidad
  wr_unif <- (integral(f = function(p) wr_p(n, c, dens_u, p),
                       xmin = 0, xmax = AQL, method = "Kron") +
                integral(f = function(p) wr_c(n, c, dens_u, p),
                         xmin = LTPD, xmax = 1, method = "Kron"))
  
  fit_ <- 
    # fit_ <- 1000 -((c/n)*N + (wr_dist-wr_unif)^2)
    #fit_ <- 00
    #fit_ <- 
}                    

# Ejecución algoritmo genético

Genetic_A <- ga(type = "binary", nBits = l1 + l2,
                fitness = function(string) fit_(string, dens_n, dens_u, N),
                popSize = 300, maxiter = 400, run = 100, pmutation=0.01)

# Solución encontrada
Sol_ga = decode(Genetic_A@solution)#; sols[g, 1] <- Sol_ga[1]; sols[g, 2] <- Sol_ga[2]
Sol_ga[1]; Sol_ga[2]



# fit_ <- (wr_dist - wr_unif)^2
# - fit_
# fit_ <- 10 - ((n/plan$n) + 3*(wr_unif/wr_dist))
#- fit_
#--- por ajustar ----
# plot
x11(); plot(Genetic_A)
prop <- seq(0, LTPD + 0.05, 0.001); prop

lab1 <- paste("clásico (", plan$n, ", ", plan$c,")")
lab2 <- paste("GA (", n_ga, ", ", c_ga,")")
lab3 <- paste("N(", round(a,2), ", ", round(b,2),")")

x11(); layout(matrix(c(rep(1,6), rep(2,3)), nrow=3, byrow=T))

plot(prop,
     phyper(plan$c, N*prop,N - (N*prop), plan$n), type="l",
     ylab= "Pa", xlab="p", lty= 2, col = "green")
lines(prop,phyper(c_ga, N*prop,N - (N*prop), n_ga), lty= 4, col = "red")
abline(v=c(AQL,LTPD), lty=2, col="blue")
abline(v=c(AQL, LTPD), lty=2, col="blue")
legend(x = "topright", legend = c(lab1, lab2), col = c("green", "red"), 
       title = "Método", lty= c(2, 4), box.lty = 0)

plot(prop, dnorm(prop, a, b), type="l", xlab= "p", ylab="fp")
abline(v=c(AQL,LTPD), lty=2, col="blue")
legend(x = "topright", legend = lab3, lty=1,
       title = "densidad", box.lty = 0)

#---- Pruebas ----

N <- 500; AQL <- 0.05; LTPD <- 0.10; alpha_des <- 0.025; beta_des <- 0.075

# Ajuste de la densidad a trabajar
dist_ <- c("normal", "uniforme")
# dist <- "uniforme"
dist <- substr(dist_, start = 1, stop = 4)

quant <- c(0.2, 0.5, 0.8) # cuantil

par1 <- AQL/(4 + qnorm(quant))  # sigma
par2<- AQL - qnorm(quant)*par1 # mu

# par1 <- LTPD/(3.5+qnorm(quant))  # sigma
# par2<- LTPD - qnorm(quant)*par1 # mu

# par2 <- 0.025
# par1 <- sqrt(par2*(1 - par2))

dens_n <- parse(text = paste(paste(
  paste("dnorm(p",sep=""),
  par2, par1, sep=","),")", sep = "")) 

dens_u <- parse(text = paste(paste(
  paste("dunif(p",sep=""),
  0, 1, sep=","),")", sep = ""))

n <- 80; c <- 3

# funciones auxiliares:
# Pa calcula la dist acum de hyper.
Pa <- function(n, c, p){
  Pa = phyper(c, N*p, N*(1-p), n)
  return(Pa)
}

# riesgo ponderado prod
wr_p <- function(n, c, densidad, p){
  f.p = eval(densidad)
  prod_wr <- (1 - Pa(n, c, p))*f.p
  return(prod_wr)
}

# riesgo ponderado cons
wr_c <- function(n, c, densidad, p){
  f.p = eval(densidad)
  cons_wr <- Pa(n, c, p)*f.p
  return(cons_wr)
}

wr_dist <- 0
for (i in 1:length(dens_n)) {
  
  wr_dist[i] <- (integral(f = function(p) wr_p(n, c, dens_n[i], p),
                          xmin = 0, xmax = AQL, method = "Kron") +
                   integral(f = function(p) wr_c(n, c, dens_n[i], p),
                            xmin = LTPD, xmax = 1, method = "Kron"))
  
}

wr_unif <- (integral(f = function(p) wr_p(n, c, dens_u, p),
                     xmin = 0, xmax = AQL, method = "Kron") +
              integral(f = function(p) wr_c(n, c, dens_u, p),
                       xmin = LTPD, xmax = 1, method = "Kron"))
wr_dist; wr_unif

wr_dist/wr_unif
wr_unif/wr_dist

x11()
p <- seq(0,1, by = 0.01)
plot(p, eval(dens_n[1])/n, type = "l", col = "orange",
     ylim = c(0,1), xlim = c(0,0.2))
lines(p, eval(dens_n[2])/n, lty = 2, col = "red")
lines(p, eval(dens_n[3])/n, lty = 2, col = "blue")
lines(p, Pa(n, c, p), lty = 2, col = "coral2")
abline(v = AQL, col = "chocolate4", lty = 3)
abline(v = LTPD, col = "chocolate4", lty = 3)
#----
