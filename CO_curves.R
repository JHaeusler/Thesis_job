set.seed(123)

N <- 300
n <- 1:N
c <- 0:max(n)

alpha <- 0.025
beta <- 0.025

AQL <- 0.05
LTPD <- 0.125

min_unif = 0.025
max_unif = 0.125
delta_prop <- 0.001

prop <- seq(0, 1, by = delta_prop)

combi <- expand.grid(n = n, c = c)

comb_filt <- combi[combi$n > combi$c,]

comb_filt_ordered <- comb_filt[order(comb_filt$n),]

co_curves_bin <- lapply(1:nrow(comb_filt_ordered), function(i) {
  pbinom(comb_filt_ordered[i, 2], comb_filt_ordered[i, 1], prop)
})

df_bin <- as.data.frame(do.call(rbind, co_curves_bin))

dens_unif <- dunif(prop, min = min_unif, max = max_unif)

ind_aql <- which.min(abs(prop - AQL))
ind_ltpd <- which.min(abs(prop - LTPD))

# --- Aplica la función para calcular los riesgos en cada fila ---
# Ahora la función devuelve un vector numérico `c()`
resultados_riesgo <- apply(df_bin, 1, function(pa_curve) {
  pa_curve <- as.numeric(pa_curve)
  pr_curve <- 1 - pa_curve
  
  risk_prod <- sum(pr_curve[1:ind_aql] * dens_unif[1:ind_aql] * delta_prop)
  
  risk_cons <- sum(pa_curve[ind_ltpd:length(prop)] * dens_unif[ind_ltpd:length(prop)] * delta_prop)
  
  weighted_risk <- risk_prod + risk_cons
  
  return(c(risk_prod, risk_cons, weighted_risk))
})

# Transpone la matriz y la convierte en un dataframe
resultados_riesgo <- as.data.frame(t(resultados_riesgo))

# Asigna nombres a las columnas para mayor claridad
colnames(resultados_riesgo) <- c("riesgo_productor", "riesgo_consumidor", "riesgo_total")
options(digits = 5)

# Une los resultados con el dataframe de combinaciones de n y c
comb_filt_ordered <- cbind(comb_filt_ordered, resultados_riesgo)

# Muestra los primeros 10 resultados para verificar
head(comb_filt_ordered, 10)
