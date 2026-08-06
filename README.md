**COPD Wearables: Sleep and Activity Phenotypes**

This repository contains analyses of wrist-worn accelerometer data and clinical characteristics from participants with chronic obstructive pulmonary disease (COPD) enrolled in home-based pulmonary rehabilitation (HBPR).

The project includes two related workflows:

Sleep and activity phenotypes: identifies participant phenotypes from wearable activity measures and evaluates their associations with program engagement and dropout.

Sleep phenotypes and HBPR engagement: evaluates whether a wearable-derived Composite Sleep Health Score improves prediction of 12-week HBPR engagement beyond standard clinical variables.

Repository structure
.
├── sleep_and_activity_phenotypes/
│   ├── phenotypes.html
│   └── phenotypes_data_dictionary.xlsx
└── sleep_phenotypes/
    ├── ml.R
    ├── ml_supplemental.R
    └── ml_documentation.docx

sleep_and_activity_phenotypes/
phenotypes.html — HTML export of the Python phenotype analysis. The workflow uses K-means clustering based on peak one-minute cadence and Movement Index, examines cluster stability, summarizes participant characteristics, and fits survival and mixed-effects models.
phenotypes_data_dictionary.xlsx — definitions for variables used in the phenotype analysis.

sleep_phenotypes/
ml.R — primary analysis comparing clinical-only models with models that also include a PLS-derived Composite Sleep Health Score.
ml_supplemental.R — supplemental analysis evaluating the Composite Sleep Health Score without the clinical covariates.
ml_documentation.docx — methodological documentation covering data acquisition, feature engineering, outcomes, modeling, evaluation, and sensitivity analyses.

The R analyses evaluate decision tree, logistic regression, linear support vector machine, and Naive Bayes classifiers using cross-validation. Reported metrics include ROC AUC, accuracy, sensitivity, specificity, confidence intervals, DeLong tests, and variable importance.

Data are not included in this repository, however, they are detailed in the following manuscripts:

sleep_and_activity_phenotypes: Faust, L., Zawada, S., Grady, M., Winham, S. J., Benzo, R. P., & Fortune, E. (2026, June 1–3). One week of pre-program actigraphy identifies dropout risk in home-based pulmonary rehabilitation for COPD [Conference paper]. 14th IEEE International Conference on Healthcare Informatics (ICHI 2026), Minneapolis, MN, United States.

sleep_phenotypes: Zawada, S. J., Faust, L., Enayati, M., Madigan, N. J., Winham, S. J., Benzo, R. P., & Fortune, E. (2026). Wearable Sleep Measures May Improve Machine Learning Prediction of Home-based Pulmonary Rehabilitation Engagement Among Patients With Chronic Obstructive Pulmonary Disease: A Proof-of-Concept Study. Mayo Clinic Proceedings: Digital Health, 100345.
