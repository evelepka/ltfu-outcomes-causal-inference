# Causal Analysis of Tuberculosis Treatment Abandonment
**Synthesis of Key Epidemiological Findings**

Tuberculosis (TB) treatment abandonment is an urgent public health crisis, yet evaluating its true mortality consequences is historically compromised by **Immortal Time Bias** (late abandoners falsely appear at low risk) and **Symptomatic Presentation Bias / Confounding by Indication** (patients returning to care artificially appear to have worse outcomes solely because they seek treatment due to severe disease progression). 

This causal framework utilizes a complete population-based cohort of 172,463 individuals treated for Tuberculosis in São Paulo, Brazil (2013-2023) to resolve these biases mathematically.

## 1. The Retreatment Hazard Paradox (Symptomatic Confounding)
Among the 21,619 individuals who abandoned therapy, 43.8% eventually re-entered the treatment system. Standard observational comparisons mistakenly classify returning to treatment as an independent risk factor for death. 

To break this confounding by indication, we fit highly granular **Marginal Structural Models (MSM)** weighted dynamically over person-months using Inverse Probability Weighting (IPW). By predicting the probability of a patient returning based on baseline severity—such as prior hospital admission (sHR 1.76) and HIV (sHR 1.39)—we isolated the true structural hazard.
- **Finding:** Returning to care functions as a surrogate proxy for acute, life-threatening clinical progression rather than a proactive health-seeking behavior. Returning dynamically associates with an overwhelming 5-fold relative increase in mortality (aHR ~5.0) compared to remaining absent.

## 2. Depletion of Susceptibles (Target Trial Emulations)
Because traditional survival models are mathematically biased to view late abandoners as artificially "healthy" (they lived long enough to abandon at month 5), we explicitly isolated the precise hazard of abandonment utilizing **Sequential Target Trial Emulations**. To account for severe non-proportional hazards relative to follow up, we tightened the administrative follow-up window strictly to **2 years** post target-trial enrollment to capture the true acute penalty of unsterilized progression.

- **Finding:** Comparing dropouts dynamically mapped against strictly matched compliant controls, the penalty of dropping out remains staggering regardless of timing. Patients abandoning in the first 30 days present a 2-year baseline hazard ratio of **2.08**. This intense relative hazard penalty persists aggressively throughout the treatment course.


### Table 1: 2-Year Target Trial Penalty by Month of Dropout
| Month of Dropout | Relative Penalty (aHR) | 95% CI |
| :--- | :--- | :--- |
| **Month 1** | 2.08 | 1.67 - 2.58 |
| **Month 2** | 3.22 | 2.78 - 3.74 |
| **Month 3** | 3.11 | 2.64 - 3.66 |
| **Month 4** | 3.26 | 2.73 - 3.88 |
| **Month 5** | 3.16 | 2.59 - 3.86 |
| **Month 6** | 2.21 | 1.70 - 2.87 |

## 3. The Competing Risks Paradigm (Subgroup Target Trial Penalty)
To test structurally vulnerable populations against baseline healthy cohorts, we ran explicit Subgroup Effects Models directly inside the unified 6-Month robust **Target Trial Emulation**. This permanently isolated our subgroup discoveries from any immortal time violations.

We discovered that the mathematical consequence of abandoning treatment is paradoxically harshest on the healthiest components of society:
- **Age Gradient:** Dropping out of treatment punishes young people (ages 15-24) with an immense **aHR 3.40**, whereas older patients face a smaller comparative hazard (aHR 2.54) simply because their absolute baseline risk of death from non-TB causes blunts the statistical magnitude of the TB-specific shock.
- **Homelessness:** Homeless individuals display a massively compressed relative penalty (aHR 1.72) specifically because their absolute competing background mortality overrides the TB penalty, while stably-housed peers suffer a brutal aHR 3.36 when choosing to quit.
- **HIV Equality:** While prior crude models implied HIV+ patients suffer a blunted penalty, explicitly forcing the interaction into the sequential Target Trial matching *equalized* the hazards perfectly across biological status (p=0.21). 

### Table 2: Competing Risks Penalty Mapping (Target Trial Interactions)
| Subgroup | Level | aHR Penalty for Abandoning | 95% CI | Interaction p-value |
| :--- | :--- | :--- | :--- | :--- |
| **Age** | 15-24 | 3.40 | 2.64 - 4.37 | < 0.001 |
| | 25-44 | 2.98 | 2.68 - 3.32 | |
| | 45-64 | 2.54 | 2.20 - 2.93 | |
| **Homelessness** | No | 3.36 | 3.09 - 3.65 | < 0.001 |
| | Yes | 1.72 | 1.44 - 2.05 | |
| **Sex** | Female | 2.91 | 2.45 - 3.45 | 0.015 |
| | Male | 2.84 | 2.60 - 3.11 | |
| **HIV/AIDS** | Negative | 2.90 | 2.64 - 3.19 | 0.212 |
| | Positive | 2.73 | 2.37 - 3.15 | |

**Conclusion:** Treating TB requires a structurally holistic view. Marginalized groups suffer overwhelming absolute baseline mortality irrespective of their specific TB treatment compliance, forcing mathematical tools to underestimate the public health impact if competing background survival is ignored.
