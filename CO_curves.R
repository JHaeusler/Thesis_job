set.seed(123)

N <- 100
n <- 1:N
c <- 0:max(n)

alpha <- 0.025
beta <- 0.025

AQL <- 0.05
LTPD <- 0.125

min_unif = 0.025
max_unif = 0.075
delta_prop <- 0.001

prop <- seq(0, 1, by = delta_prop)

combi <- expand.grid(n = n, c = c)

comb_filt <- combi[combi$n > combi$c,]

comb_filt_ordered <- comb_filt[ordered(comb_filt$n),]

co_curves_bin <- lapply(1:nrow(comb_filt_ordered), function(i) {
  pbinom(comb_filt_ordered[i, 2], comb_filt_ordered[i, 1], prop)
})

df_bin <- as.data.frame(do.call(rbind, co_curves_bin))

dens_unif <- dunif(prop, min = min_unif, max = max_unif)

# --- Corrección del cálculo de riesgos ---

# Probabilidad de Aceptación (PA) y de Rechazo (PR) para la primera curva
pa_curve <- as.numeric(df_bin[1, ])
pr_curve <- 1 - pa_curve

# Encuentra los índices de AQL y LTPD en el vector `prop`
ind_aql <- which.min(abs(prop - AQL))
ind_ltpd <- which.min(abs(prop - LTPD))

# 1. Riesgo del Productor (integrado de 0 a AQL)
# Es (probabilidad de rechazo) * densidad
risk_prod <- sum(pr_curve[1:ind_aql] * dens_unif[1:ind_aql] * delta_prop)

# 2. Riesgo del Consumidor (integrado de LTPD a 1)
# Es (probabilidad de aceptación) * densidad
risk_cons <- sum(pa_curve[ind_ltpd:length(prop)] * dens_unif[ind_ltpd:length(prop)] * delta_prop)

# Riesgo Ponderado Total es la suma de ambos
weighted_risk <- risk_prod + risk_cons

print(weighted_risk)