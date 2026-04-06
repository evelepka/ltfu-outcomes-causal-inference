# Dicionário de Dados Oficial: ITT Cohort (Mortalidade e Abandono)

Este documento foi criado para garantir a reprodutibilidade das análises causais por outros pesquisadores ou membros da equipe. Ele detalha as variáveis da base `itt_cohort.csv`.

## 1. Identificação e Estrutura Principal
* `sinan_clean`: Identificador único do paciente (anonimizado).
* `itt_group`: Grupo de interesse causal ("Loss to follow-up" ou "Non-LTFU"). Definido com base no desfecho do **primeiro episódio Novo** do paciente.
* `tx_month_grp`: Duração categorizada do tratamento índice antes do desfecho (ex: "< 1 month", "≥ 4 months").

## 2. Variáveis de Tempo (Obrigatórias para Análise de Sobrevida)
*Atenção: Para evitar o Viés de Tempo Imortal (Immortal Time Bias), o tempo de início do tratamento é definido por uma lógica Proxy (melhor data disponível).*

* `best_start`: Data "Proxy" de Início. Ordem de precedência: Data de Início do Tratamento (`tx_start`) -> Data do Diagnóstico -> Data da Notificação.
* `end_date`: Data em que o tratamento índice terminou (por cura ou abandono). **Marca o início do Landmark Pós-Tratamento.**
* `death_date`: Data de óbito validada (consolidada entre SIM e SINAN).
* `time_d_tx`: Tempo de sobrevida calculado a partir do `best_start` (em anos). Usada no G-Formula clássico.
* `time_d`: Tempo de sobrevida calculado a partir da `end_date` (em anos). Usada na Análise Landmark (Pós-Tratamento).
* `event_d`: Indicador de evento de morte (1 = Morreu, 0 = Censurado em Dez/2024).

## 3. Variáveis de Confusão Harmonizadas (Tabela 2)
Estas **13 covariáveis** foram rigorosamente selecionadas e recategorizadas para ajuste em todos os modelos causais:

1. `age_group`: Faixa etária (15-24, 25-44, 45-64, ≥65).
2. `sex`: Sexo biológico (Female, Male).
3. `race_clean`: Raça harmonizada (White, Black or Mixed, Other).
4. `edu_clean`: Escolaridade (None, ≤ 7 years, 8 - 11 years, ≥ 12 years).
5. `hiv_aids`: Status de coinfecção HIV (Negative, Positive).
6. `dot_status`: Terapia Diretamente Observada (Yes, No). *Mapeada do campo original tx_administration_type.*
7. `alcohol`: Etilismo (Yes, No).
8. `drug_use`: Uso de drogas ilícitas (Yes, No).
9. `incarcerated`: Situação prisional (Yes, No). *Mapeada do tipo de residência.*
10. `homelessness`: População em situação de rua (Yes, No). *Mapeada do tipo de residência.*
11. `hosp_admission`: Internação hospitalar relacionada à TB (Yes, No).
12. `clinical_clean`: Forma clínica consolidada (Pulmonary, Extrapulmonary, Pulmonary and Extrapulmonary).
13. `diabetes`: Diagnóstico de Diabetes Mellitus (Yes, No).

## 4. Variáveis Secundárias e Competitivas
* `time_rn` / `time_rn_tx`: Tempo até um novo episódio de re-tratamento.
* `event_rn`: Indicador de re-tratamento (Competeting Event).
* `mental_health`, `tobacco_use`, `other_immuno_condition`: Outras morbidades basais (para análises exploratórias em subgrupos).
* `diagnosis_setting`: Onde a doença foi descoberta (Ex: Outpatient, Emergency / Inpatient).

## 5. Como usar na prática (R)
```r
# Exemplo de carregamento e modelo Cox básico considerando Landmark pós-tratamento
df <- read.csv("data/itt_cohort.csv")
cox_model <- coxph(
  Surv(time_d, event_d) ~ itt_group + age_group + sex + race_clean + edu_clean + 
                          hiv_aids + dot_status + alcohol + drug_use + incarcerated + 
                          homelessness + hosp_admission + clinical_clean + diabetes,
  data = df
)
```
---
*Gerado a partir do Protocolo de Harmonização de Cohort (2026).*
