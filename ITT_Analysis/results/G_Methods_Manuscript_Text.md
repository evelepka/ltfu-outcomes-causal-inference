# Suggested Manuscript Text: G-Methods and Evidence-Based Triangulation

*This document contains a structured academic draft for the **Methods** and **Results** sections of your manuscript, focusing on how to report the G-formula analyses, the inherent biases, and the triangulation strategy (QBA and E-Value).*

---

## Methods

### Causal Inference via Parametric G-Formula
To estimate the causal effect of loss to follow-up (LTFU) on long-term mortality, we employed the parametric G-formula. This approach allows for the estimation of standardized absolute risks (cumulative incidence) and risk ratios under hypothetical interventions (e.g., intervening to prevent all LTFU versus allowing observed LTFU rates), standardized to the baseline confounder distribution of the entire cohort. Models were adjusted for baseline covariates including age, sex, race, education, HIV/AIDS status, diabetes, alcohol and drug use, incarceration history, homelessness, hospital admission at diagnosis, clinical presentation, and DOT status.

### Addressing Immortal Time and Selection Biases: A Triangulation Approach
Evaluating the long-term mortality impact of LTFU is methodologically challenging due to two competing dynamic biases: 
1. **The "Sick-to-stay" bias (Selection/Immortal Time Bias):** Patients with the most severe disease presentations are more likely to die early during the intensive phase of treatment before they have the opportunity to abandon. Because LTFU only occurs among those who survive long enough to drop out, simple comparisons from treatment initiation (Time 0) paradoxically bias the LTFU group towards better apparent survival.
2. **The "Healthy-to-quit" bias (Unmeasured Confounding):** Conversely, patients who experience rapid clinical improvement might abandon treatment early because they feel cured, introducing unmeasured confounding related to uncaptured positive clinical severity metrics.

To untangle these biases and establish a robust causal narrative, we employed a multi-pronged **Evidence-Based Triangulation Strategy**, comparing estimates across three distinct analytical timelines:
* **Strategy 1 (Time 0 - Treatment Initiation):** G-formula estimation starting at the date of treatment initiation. We hypothesized this estimate would suffer severely from the sick-to-stay bias, serving as our biased baseline.
* **Strategy 2 (180-Day Landmark):** To mitigate early mortality survival bias and isolate post-treatment effects, we applied a landmark analysis restricted to individuals who survived the intensive phase (at least 180 days). Follow-up time origin was set to day 180.
* **Strategy 3 (Treatment Outcome Landmark):** Follow-up was strictly isolated to the post-treatment period by defining the time origin as the exact date of the treatment outcome (date of cure/completion or date of abandonment).

### Quantitative Bias Analysis (QBA) and E-Value
To formally evaluate the susceptibility of our estimates to these specific biases:
* **QBA for Time 0:** We conducted a multidimensional Quantitative Bias Analysis to determine what magnitude of unmeasured baseline severity (U)—disproportionately affecting treatment completers—would be required to fully explain any paradoxical protective effect observed at Time 0. We modeled $U$ as a binary unmeasured confounder strongly associated with early mortality ($RR_{UY} = 4.5$), with a significantly higher prevalence among completers ($30\%$) compared to abandoners ($5\%$).
* **E-Value for Landmark Analysis:** To assess the robustness of our primary causal effect (Strategy 2) against the "Healthy-to-quit" unmeasured confounding, we calculated the E-value. This quantifies the minimum strength of association that an unmeasured confounder would need to have with both LTFU and mortality to shift the observed risk ratio to the null.

---

## Results

### Evaluating the Methodological Biases (Triangulation)

**Strategy 1 (Time 0 - Treatment Initiation):** 
When starting follow-up at the time of treatment initiation, the crude analysis heavily favored the LTFU group. Parametric G-formula estimation from Time 0 yielded a 12-year cumulative mortality of 8.19% for the non-LTFU group compared to 7.76% for the LTFU group, resulting in an apparent, paradoxical Risk Ratio (RR) of 0.94. 

**Quantitative Bias Analysis (QBA) of the Sick-to-Stay Bias:**
Recognizing this paradoxical finding as an artifact of the sick-to-stay survival bias, we applied a QBA. Assuming the presence of an unmeasured severity factor highly prevalent among completers who die early, mathematical correction of the Time 0 estimate reversed the effect entirely. After adjusting for this hypothesized bias, LTFU was associated with a 65% increased risk of mortality (Adjusted RR: 1.65). 

**Strategy 2: The 180-Day Landmark:**
Restricting the analysis to patients who survived the intensive phase of treatment successfully mitigated the early mortality bias. In this landmark analysis, the true long-term impact of LTFU emerged clearly. The 10-year G-formula estimated cumulative mortality was 3.25% for the non-LTFU group and 7.06% for the LTFU group. At 12 years post-landmark, LTFU was causally associated with a significantly higher risk of death (Risk Difference: 4.02%; RR: 2.15).

**Strategy 3: Outcome-Day Landmark:**
Aligning the follow-up strictly to the exact day of treatment outcome (cure/complete or abandonment) confirmed the Landmark 180-day findings. At 10 years post-outcome, the estimated cumulative mortality was 3.24% for non-LTFU versus 6.67% for LTFU (RR: 2.05). At 12 years, the causal risk ratio remained elevated at 1.94.

**Robustness to Unmeasured Confounding (E-Value):**
To ensure the robust 2.15-fold increased risk observed in our 180-day Landmark model was not driven by the competing "Healthy-to-quit" bias or other unmeasured factors, computational E-values were derived. The E-value for the 180-day estimate was 3.72 (lower confidence bound: 3.31). This indicates that an unmeasured confounder would have to be associated with both a 3.72-fold increase in the likelihood of abandoning treatment and a 3.72-fold increase in the risk of long-term mortality—above and beyond all clinical and demographic factors already included in our model—to explain away the observed effect. As no single measurable clinical factor in our rich TB dataset approached this magnitude, the causal detrimental effect of LTFU is considered highly robust.
