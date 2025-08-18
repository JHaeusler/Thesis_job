# Trabajo Investigación ----

lista_de_paquetes <- c("pracma", "GA", "AcceptanceSampling",
                       "scatterplot3d", "tidyverse", "reshape",
                       "sf", "dplyr", "RColorBrewer") # Reemplaza con tus paquetes

for (paquete in lista_de_paquetes) {
  if(!require(paquete, character.only = TRUE)){
    install.packages(paquete)
    require(paquete, character.only = TRUE)
  }
}

# semilla
set.seed(123)

cols_a <- c("chocolate", "chartreuse4", "chartreuse3", "coral3", "coral1",
            "aquamarine4", "antiquewhite4")
#---- Parámetros de entrada ----
N <- 150; AQL <- 0.05; LTPD <- 0.10; alpha_des <- 0.05; beta_des <- 0.10
par11 <- 0; par12 <- 0.05
dens_1 <- parse(text = paste(paste(
  paste("dunif(p",sep=""),
  par11, par12, sep=","),")", sep = "")) 

par21 <- 0.1; par22 <- 1
dens_2 <- parse(text = paste(paste(
  paste("dunif(p",sep=""),
  par21, par22, sep=","),")", sep = "")) 

par31 <- 0; par32 <- 1
dens_3 <- parse(text = paste(paste(
  paste("dunif(p",sep=""),
  par31, par32, sep=","),")", sep = "")) 

### ###
par41 <- 0; par42 <- 0.2
dens_4 <- parse(text = paste(paste(
  paste("dunif(p",sep=""),
  par41, par42, sep=","),")", sep = "")) 

par51 <- 0.2; par52 <- 0.4
dens_5 <- parse(text = paste(paste(
  paste("dunif(p",sep=""),
  par51, par52, sep=","),")", sep = "")) 

par61 <- 0.4; par62 <- 0.6
dens_6 <- parse(text = paste(paste(
  paste("dunif(p",sep=""),
  par21, par22, sep=","),")", sep = "")) 

par71 <- 0.6; par72 <- 0.8
dens_7 <- parse(text = paste(paste(
  paste("dunif(p",sep=""),
  par71, par72, sep=","),")", sep = "")) 

par81 <- 0.8; par82 <- 1
dens_8 <- parse(text = paste(paste(
  paste("dunif(p",sep=""),
  par81, par82, sep=","),")", sep = "")) 

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

# riesgo ponderado comp
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


n.opt=1:N;c.opt=0:N; Esc=expand.grid(c.opt=c.opt,n.opt=n.opt);rm(n.opt,c.opt);Esc=Esc[Esc$n.opt>=Esc$c.opt,]
# View(Esc);dim(Esc)

for(i in 1:nrow(Esc)){
  Esc[i,3] <- calc_wr(Esc[i,2],Esc[i,1],dens_1,N,AQL,LTPD)
  Esc[i,4] <- calc_wr(Esc[i,2],Esc[i,1],dens_2,N,AQL,LTPD)
  Esc[i,5] <- calc_wr(Esc[i,2],Esc[i,1],dens_3,N,AQL,LTPD)
  Esc[i,6] <- calc_wr(Esc[i,2],Esc[i,1],dens_4,N,AQL,LTPD)
  Esc[i,7] <- calc_wr(Esc[i,2],Esc[i,1],dens_5,N,AQL,LTPD)
  Esc[i,8] <- calc_wr(Esc[i,2],Esc[i,1],dens_6,N,AQL,LTPD)
  Esc[i,9] <- calc_wr(Esc[i,2],Esc[i,1],dens_7,N,AQL,LTPD)
  Esc[i,10] <- calc_wr(Esc[i,2],Esc[i,1],dens_8,N,AQL,LTPD)
}
limite_inferior_resaltado <- 0.049
limite_superior_resaltado <- 0.050
Esc_resaltado <- Esc %>%
  mutate(
    # Crea una nueva columna 'es_riesgo_critico'
    es_riesgo_critico6 = (V6 >= limite_inferior_resaltado & V6 <= limite_superior_resaltado)
  )%>%
  mutate(
    # Crea una nueva columna 'es_riesgo_critico'
    es_riesgo_critico7 = (V7 >= limite_inferior_resaltado & V7 <= limite_superior_resaltado)
  )%>%
  mutate(
    # Crea una nueva columna 'es_riesgo_critico'
    es_riesgo_critico8 = (V8 >= limite_inferior_resaltado & V8 <= limite_superior_resaltado)
  )%>%
  mutate(
    # Crea una nueva columna 'es_riesgo_critico'
    es_riesgo_critico9 = (V9 >= limite_inferior_resaltado & V9 <= limite_superior_resaltado)
  )%>%
  mutate(
    # Crea una nueva columna 'es_riesgo_critico'
    es_riesgo_critico10 = (V10 >= limite_inferior_resaltado & V10 <= limite_superior_resaltado)
  )

# print(Esc_resaltado %>% filter(es_riesgo_critico10 == TRUE))
valores_n_opt <- sort(unique(Esc$n.opt))
breaks_n <- c(1, seq(5, max(valores_n_opt, na.rm = TRUE), by = 5))
breaks_n_final <- intersect(breaks_n, valores_n_opt)

col_sel <- Esc$V9

x11()
ggplot(data = Esc_resaltado, aes(x = as.factor(n.opt), y = as.factor(c.opt),
                                 fill = V10)) +
  geom_tile(
    aes(color = es_riesgo_critico10), # Mapea el color del borde a la nueva variable
    linewidth = ifelse(Esc_resaltado$es_riesgo_critico10, 0.5, 0.2) # Grosor del borde: 1.5 si TRUE, 0.5 si FALSE
  ) +
  scale_fill_viridis_c(
    option = "cividis",
    name = "Riesgo ponderado",
    breaks = seq(0, 1, 0.1)
  ) +
  # --- Configura la escala de color para el BORDE ---
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "white"), # Color rojo para el borde de celdas críticas, blanco para el resto
    guide = "none" # Desactiva la leyenda para el color del borde, ya que es auxiliar
  ) +
  scale_x_discrete(breaks = breaks_n_final) + # Muestra solo la primera etiqueta
  # Agrega esto para controlar las etiquetas del eje Y
  scale_y_discrete(breaks = function(x) x[seq(1, length(x), by = 5)]) +
  labs(
    title = "Riesgo ponderado planes de muestreo simples para atributos",
    subtitle = paste0("Planes n, c (Riesgo Crítico menor a ", limite_superior_resaltado, ")"),
    x = "n",
    y = "c"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 1)
  )

# solución mediante AcceptanceSampling (clásico)
plan <- find.plan(PRP = c(AQL, 1 - alpha_des),
                  CRP = c(LTPD, beta_des),
                  N=N[1], type = "hypergeom")#; plan$n; plan$c
# rm(plan)
jujuj <- Esc$V9
cumplen <- which(jujuj >= limite_inferior_resaltado & jujuj <= limite_superior_resaltado)
n_cumplen <- Esc[cumplen, 2]; c_cumplen <- Esc[cumplen, 1] 

prop <- seq(0, 1, 0.01)
# 6 <- 4
x11()
titulo_expresion <- bquote(atop(paste("COO (Tipo A) Lote de Tamaño N=", .(N[1]),
                                      hat(p) %~% (.(par71) ~ "," ~ .(par72)))))
plot(prop, phyper(plan$c, N[1]*prop, N[1]*(1 - prop), plan$n), lty = 2,
     main = titulo_expresion, type = "l", lwd = 1.5,
     col = cols_a[1], ylab = "Pa", xlab = "p")#, xlim = c(0, 0.5))
for (resp in 1:length(cumplen)) {
  lines(prop, phyper(c_cumplen[resp], N[1]*prop,N[1]*(1 - prop), n_cumplen[resp]), lty = 1,
      col = cols_a[resp+1], lwd = 1.5)
}
legend(x = "topright", legend = paste("n=", c(plan$n,n_cumplen), "c=", c(plan$c,c_cumplen)),
       col = cols_a[1:(length(cumplen)+1)], 
       title = "Plan", lty= c(2,rep(1, length(cumplen))), box.lty = 0, bg = NA)
# abline(v=c(0.05, 0.10), lty = 3, col = c("burlywood4"))
segments(x0= 0.05, y0= 0, x1= 0.05, y1= 1.1, lty = 3, col = c("burlywood4"))
segments(x0= 0.10, y0= 0, x1= 0.10, y1= 1.1, lty = 3, col = c("burlywood4"))
segments(x0= 0, y0= 1, x1= .2, y1= 1, lty = 2, col = c("azure3"))
abline(v= 0, h= c(0), lty = 1, col = c("azure3"))
text(x= 0.05, y= -0.02, label= "AQL", cex = 0.6)
text(x= 0.10, y= -0.02, label= "LTPD", cex = 0.6)



# text(x=0.042, y = 0.6, label = "c = 3", col = "darkorchid")



ggplot(data = Esc, aes(x = as.factor(n.opt), y = as.factor(c.opt),
                       fill = col_sel)) +
  geom_tile(color = "white", size = 0.05) + # Dibuja las celdas, con un borde blanco
  #geom_text(aes(label = round(V6, 3)), color = "black", size = 1.7) +
  scale_fill_viridis_c(option = "cividis", name = "Riesgo ponderado",
                       breaks = seq(0,1, 0.1)) +#,
                       #labels = as.character(seq(0,1, 0.01))) + # Escala de color para la medida continua

  labs(
  title = "Riesgo ponderado planes de muestreo simples para atributos",
  subtitle = "planes n, c",
  x = "n",
  y = "c"
  ) +
  theme_minimal() + # Un tema limpio
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 1) # Rota las etiquetas del eje X si son largas
  )
# x11()
# ggplot(data = Esc, aes(x = as.factor(n.opt), y = as.factor(c.opt),
#                        fill = V4)) +
#   geom_tile(color = "white", size = 0.5) + # Dibuja las celdas, con un borde blanco
#   scale_fill_viridis_c(option = "n, c", name = "Riesgo ponderado") + # Escala de color para la medida continua
#   labs(
#     title = "Riesgo ponderado planes de muestreo simples para atributos",
#     subtitle = "plasnes n, c",
#     x = "n",
#     y = "c"
#   ) +
#   theme_minimal() + # Un tema limpio
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold"),
#     plot.subtitle = element_text(hjust = 0.5),
#     axis.text.x = element_text(angle = 45, hjust = 1) # Rota las etiquetas del eje X si son largas
#   )
# x11()
# ggplot(data = Esc, aes(x = as.factor(n.opt), y = as.factor(c.opt),
#                        fill = V5)) +
#   geom_tile(color = "white", size = 0.5) + # Dibuja las celdas, con un borde blanco
#   scale_fill_viridis_c(option = "n, c", name = "Riesgo ponderado") + # Escala de color para la medida continua
#   labs(
#     title = "Riesgo ponderado planes de muestreo simples para atributos",
#     subtitle = "plasnes n, c",
#     x = "n",
#     y = "c"
#   ) +
#   theme_minimal() + # Un tema limpio
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold"),
#     plot.subtitle = element_text(hjust = 0.5),
#     axis.text.x = element_text(angle = 45, hjust = 1) # Rota las etiquetas del eje X si son largas
#   )
#---- mapa de coropletas ----