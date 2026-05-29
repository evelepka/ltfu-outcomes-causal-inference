# 21_itt_qba_evalue_analysis.R

# Setup pacotes essenciais de inferencia causal e analise de sensibilidade
if (!require(episensr)) install.packages("episensr", repos="http://cran.rstudio.com/")
if (!require(EValue)) install.packages("EValue", repos="http://cran.rstudio.com/")

library(episensr)
library(EValue)

cat("==================================================================\n")
cat("=== TENTATIVA DE DERROTAR O SICK-TO-STAY BIAS COM QBA ==========\n")
cat("==================================================================\n")

# Pilar II: QBA (Quantitative Bias Analysis) Determinístico 
# Observamos no modelo marginal absoluto G-Formula (Tempo 0): RR de LFTU = 0.94
# Premissas do Viés Oculto U ("Morbidade Terminal Precoce"):
# - RR_UY: Essa morbidade eleva o risco de morte em 4.5 vezes.
# - P(U|A=1): Prevalencia em Abandonadores = 5% (Geralmente estáveis o suficiente pra fugir)
# - P(U|A=0): Prevalencia em Complacentes/Non-LTFU = 30% (Ficam no sistema pq estão graves/hospitalizados/moribundos)

# Fórmula Pura do Risco Relativo Ajustado
RR_obs <- 0.947
RR_UY <- 4.5
Prev_U_Exp <- 0.05  # P(U|A=1) Abandoners
Prev_U_Unexp <- 0.30 # P(U|A=0) Completers

bias_factor <- (RR_UY * Prev_U_Exp + (1 - Prev_U_Exp)) / (RR_UY * Prev_U_Unexp + (1 - Prev_U_Unexp))
RR_adj <- RR_obs / bias_factor

cat("-> RR Observado (Paradoxal):", RR_obs, "\n")
cat("-> Bias Factor (Fator de Correção):", bias_factor, "\n")
cat("-> RR Ajustado (Real):", RR_adj, "\n")


cat("\n==================================================================\n")
cat("=== ROBUSTEZ DO EFEITO REAL PÓS-SOBREVIVÊNCIA (E-VALUE) ========\n")
cat("==================================================================\n")

# Pilar III: E-Value no Landmark Analítico (Aos 180 Dias)
# Nosso modelo Landmark limpo que exigiu sobrevivência de 180 dias encontrou um RR nocivo de 2.15.
# Quão forte precisaria ser um Confounding Oculto para destruir esse achado real?

# Assumindo RR observado = 2.15. 
# Calcularemos o E-value pontual que nulifica (traz o RR real pra 1.0)
e_val <- evalues.RR(est = 2.15, lo = 1.95, hi = 2.37)

print(e_val)

cat("\nDone!\n")
