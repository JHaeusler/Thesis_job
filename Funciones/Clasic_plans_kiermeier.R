# función para completar la matriz Esce_Global que será el insumo
# para lograr el objetivo 3; evaluar multiples combinaciones
# generados entre los inputs: N, AQL, LTPD, alpha_des, beta_des
calcular_kiermeier <- function(N, AQL, LTPD, alpha_des, beta_des) {
  tryCatch({
    # PRP = (AQL, 1 - alpha) -> Punto del Productor
    # CRP = (LTPD, beta)     -> Punto del Consumidor
    plan <- find.plan(
      PRP = c(AQL, 1 - alpha_des),
      CRP = c(LTPD, beta_des),
      type = "hypergeom", # O "binomial" según tu marco teórico
      N = N
    )
    return(c(n = plan$n, c = plan$c))
  }, error = function(e) {
    return(c(n = NA, c = NA))
  })
}