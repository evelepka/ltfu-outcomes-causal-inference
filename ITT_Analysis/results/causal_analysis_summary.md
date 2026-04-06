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

- **Finding:** Comparing dropouts dynamically mapped against strictly matched compliant controls, the penalty of dropping out remains staggering regardless of timing. Patients abandoning in the first 30 days present a 2-year baseline hazard ratio of **2.42**. This intense relative hazard penalty persists aggressively throughout the treatment course, punishing patients abandoning at Month 4 (aHR ~3.12) as extensively as those leaving earlier. 

## 3. The Competing Risks Paradigm (Subgroup Target Trial Penalty)
To test structurally vulnerable populations against baseline healthy cohorts, we ran Explicit Subgroup Effects Models explicitly inside the unified 6-Month robust **Target Trial Emulation**. This permanently isolated our subgroup discoveries from any immortal time violations.

We discovered that the mathematical consequence of abandoning treatment is paradoxically harshest on the healthiest components of society:
- **Age Gradient:** Dropping out of treatment punishes young people (ages 15-24) with an immense **aHR 3.40**. Older patients (45-64) face a smaller comparative hazard (aHR 2.54) simply because their absolute baseline risk of death from non-TB causes blunts the statistical magnitude of the TB-specific shock.
- **Homelessness:** Homeless individuals display a massively compressed relative penalty (aHR 1.72) specifically because their absolute competing background mortality overrides the TB penalty, while stably-housed peers suffer a brutal aHR 3.36 when choosing to quit.
- **HIV Errata:** While crude models implied HIV+ patients suffer a blunted penalty, explicitly forcing the interaction into the Target Trial matching *equalized* the hazards perfectly across biological status (HIV- aHR 2.90 vs HIV+ aHR 2.73, p=0.21). The previously suspected phenomenon was purely an artifact of skewed longitudinal drop-out timing.

**Conclusion:** Treating TB requires a structurally holistic view. Marginalized groups suffer overwhelming absolute baseline mortality irrespective of their specific TB treatment compliance, forcing mathematical tools to underestimate the public health impact if competing background survival is ignored.
