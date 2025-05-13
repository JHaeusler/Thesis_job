#--- Bayesian Inference ----

setwd("D:/Github/Thesis_job")

lista_de_paquetes <- c("readxl", "ggplot2", "readr", "tidyr") # Reemplaza con tus paquetes

for (paquete in lista_de_paquetes) {
  if(!require(paquete, character.only = TRUE)){
    install.packages(paquete)
    require(paquete, character.only = TRUE)
  }
}

data1 <- read_excel("Acceptance Sampling MIL-STD 105E for Quality Control.xlsx",
                    sheet = "tab")
data2 <- read_excel("Jurnal Rekavasi.xlsx", sheet = "tabla")






