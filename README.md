# COPD Wearables: Sleep and Activity Phenotypes
This repository contains analyses of wrist-worn accelerometer data and clinical characteristics from participants with chronic obstructive pulmonary disease (COPD) enrolled in home-based pulmonary rehabilitation (HBPR).

The project includes two related workflows:
- **Sleep and activity phenotypes:** Identifies participant phenotypes from wearable activity measures and evaluates their associations with program engagement and dropout.
- **Sleep phenotypes and HBPR engagement:** Evaluates whether a wearable-derived Composite Sleep Health Score improves prediction of 12-week HBPR engagement beyond standard clinical variables.

## Repository structure
```text
.
├── sleep_and_activity_phenotypes/
│   ├── phenotypes.html
│   ├── phenotypes_data_dictionary.xlsx
│   └── phenotypes_documentation.docx
└── sleep_phenotypes/
    ├── ml.R
    ├── ml_supplemental.R
    └── ml_documentation.docx
```

## Repository contents

### `sleep_and_activity_phenotypes/`

- **`phenotypes.html`** — HTML export of the Python phenotype analysis. The workflow uses K-means clustering based on peak one-minute cadence and Movement Index, examines cluster stability, summarizes participant characteristics, and fits survival and mixed-effects models.
- **`phenotypes_data_dictionary.xlsx`** — Definitions for variables used in the phenotype analysis.
- **`phenotypes_documentation.docx`** — Methodological documentation covering data processing, feature engineering, phenotype creation, evaluation, and outcomes.

### `sleep_phenotypes/`

- **`ml.R`** — Primary analysis comparing clinical-only models with models that also include a partial least squares (PLS)-derived Composite Sleep Health Score.
- **`ml_supplemental.R`** — Supplemental analysis evaluating the Composite Sleep Health Score without clinical covariates.
- **`ml_documentation.docx`** — Methodological documentation covering data acquisition, feature engineering, outcomes, modeling, evaluation, and sensitivity analyses.

The R analyses evaluate decision tree, logistic regression, linear support vector machine, and Naive Bayes classifiers using cross-validation. Reported metrics include receiver operating characteristic area under the curve (ROC AUC), accuracy, sensitivity, specificity, confidence intervals, DeLong tests, and variable importance.

## Related manuscripts

The data and analyses are described in the following manuscripts:

### Sleep and activity phenotypes

Faust, L., Zawada, S., Grady, M., Winham, S. J., Benzo, R. P., & Fortune, E. (2026, June 1–3). *One week of pre-program actigraphy identifies dropout risk in home-based pulmonary rehabilitation for COPD* [Conference paper]. 14th IEEE International Conference on Healthcare Informatics (ICHI 2026), Minneapolis, Minnesota, United States.

### Sleep phenotypes and HBPR engagement

Zawada, S. J., Faust, L., Enayati, M., Madigan, N. J., Winham, S. J., Benzo, R. P., & Fortune, E. (2026). Wearable sleep measures may improve machine learning prediction of home-based pulmonary rehabilitation engagement among patients with chronic obstructive pulmonary disease: A proof-of-concept study. *Mayo Clinic Proceedings: Digital Health*, 100345.
