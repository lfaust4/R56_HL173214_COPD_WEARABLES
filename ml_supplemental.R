#install.packages("rpart.plot")
library(rpart.plot)
#remove.packages("mice")
#install.packages("mice", dependencies = TRUE)
library(mice)
# Install if not already installed
#install.packages("caret")
#install.packages("pROC")
#install.packages("tidyverse")  # Optional, for dplyr/mutate
#install.packages("rlang")
#install.packages("tidymodels")  # Optional, for dplyr/mutate
library(sandwich)
library(dplyr)
library(lubridate)
library(tidyverse)
library(caret)
library(pROC)
library(tidymodels)
library(caret)
library(yardstick)
library(ggplot2)
library(rpart)
library(workflows)
library(workflowsets)
library(tidymodels)
#install.packages("vip")  # install the package
library(vip)             # load it
library(parsnip)
library(workflows)
library(tune)


df_imputed <- read.csv("imputeddf.csv")
df_interfirst <- read.csv("file_participant_arm_map.csv")
#data_adh <- read.csv("adherence_metrics1.csv")

#pca_data <- df_imputed[, c("mean_waso", "sd_waso", "mean_time", "sd_time", "mean_eff", "mean_noa", "mean_loaim", "mean_ac", "mean_mi", "mean_fi", "mean_sfi")]
#pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
#summary(pca_result)
#pca_result$rotation
#head(pca_result$x)
#pc1_scores <- pca_result$x[, 1]
#df_imputed$PC1 <- pc1_scores

df_filtered_interfirst <- df_interfirst[grepl("^I", df_interfirst$arm), ]
df_f_interfirst <- df_filtered_interfirst[grepl("^b", df_filtered_interfirst$period), ]
df_imputed <- semi_join(df_imputed, df_f_interfirst, by = "Subject.Name")
df_imputed$current_smoker_num_a <- ifelse(df_imputed$current_smoker_a == "Yes", 1, 0)
df_imputed$adh_yes_factor <- as.factor(df_imputed$adh_yes_factor)

write.csv(df_imputed, "df_intervfirst.csv", row.names = FALSE)

df_imputed <- df_imputed[, c("mean_waso", "sd_waso", "mean_time", "sd_time" , "mean_eff", "mean_noa", "mean_loaim", "mean_ac", "mean_mi", "mean_fi", "mean_sfi", "mmrc_a", "age2_a", "sex_numeric_a", "charlson_a", "fev1_updated_a", "current_smoker_num_a", "adh_yes_factor")]
outer_folds <- vfold_cv(df_imputed, v = 10, strata = adh_yes_factor)

#df_imputed <- df_imputed %>%
#  mutate(sex_numeric_a = as.integer(sex_numeric_a))
# Sensor variables to use in PLS
pls_vars <- setdiff(c("mean_waso", "sd_waso", "mean_time", "sd_time" , "mean_eff", "mean_noa", "mean_loaim", "mean_ac", "mean_mi", "mean_fi", "mean_sfi"), "adh_yes_factor")
# Covariates to include without PLS
covariates <- setdiff(c("mmrc_a", "age2_a", "sex_numeric_a", "charlson_a", "fev1_updated_a", "current_smoker_num_a"), "adh_yes_factor")
# Outcome variable
outcome <- "adh_yes_factor"

recipe_no_pls <- recipe(adh_yes_factor ~ ., data = df_imputed) %>%
  update_role(all_of(covariates), new_role = "predictor") %>%
  update_role(all_of(pls_vars), new_role = "ignore") %>%
  step_dummy(all_nominal_predictors()) %>%
  step_normalize(all_numeric_predictors())

recipe_pls <- recipe(adh_yes_factor ~ ., data = df_imputed) %>%
  update_role(all_of(covariates), new_role = "ignore") %>%
  update_role(all_of(pls_vars), new_role = "predictor") %>%
  step_dummy(all_nominal_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_pls(all_of(pls_vars), outcome = "adh_yes_factor", num_comp = 1, prefix = "PLS_") %>%
  step_normalize(starts_with("PLS_"))

#prep(recipe_pls) %>% juice() %>% names()
#library(recipes)
#prepped <- prep(recipe_pls)
#summary(prepped)

tree_spec <- decision_tree() %>%
  set_engine("rpart") %>%
  set_mode("classification")  # or regression

wf_no_pls <- workflow() %>%
  add_model(tree_spec) %>%
  add_recipe(recipe_no_pls)

wf_pls <- workflow() %>%
  add_model(tree_spec) %>%
  add_recipe(recipe_pls)

res_no_pls <- fit_resamples(
  wf_no_pls,
  resamples = outer_folds,
  metrics = metric_set(roc_auc, accuracy),
  control = control_resamples(save_pred = TRUE)
)

res_pls <- fit_resamples(
  wf_pls,
  resamples = outer_folds,
  metrics = metric_set(roc_auc, accuracy),
  control = control_resamples(save_pred = TRUE)
)

results <- bind_rows(
  collect_metrics(res_no_pls) %>% mutate(model = "No PLS"),
  collect_metrics(res_pls) %>% mutate(model = "PLS + Baseline")
)

# Suppose res_pls and res_no_pls were run with save_pred = TRUE
preds_pls <- collect_predictions(res_pls) %>% arrange(.row) 
preds_no_pls <- collect_predictions(res_no_pls) %>% arrange(.row)

roc_pls <- roc(preds_pls$adh_yes_factor, preds_pls$.pred_Yes)
roc_no_pls <- roc(preds_no_pls$adh_yes_factor, preds_no_pls$.pred_Yes)

# Delong's test
roc.test(roc_pls, roc_no_pls, method = "delong")

ci_pls <- ci.auc(roc_pls)
ci_no_pls <- ci.auc(roc_no_pls)
ci_pls
ci_no_pls

# No PLS
conf_mat_no_pls <- conf_mat(preds_no_pls, truth = adh_yes_factor, estimate = .pred_class)
conf_mat_no_pls

# PLS + baseline
conf_mat_pls <- conf_mat(preds_pls, truth = adh_yes_factor, estimate = .pred_class)
conf_mat_pls

# No PLS
accuracy_no_pls <- accuracy(preds_no_pls, truth = adh_yes_factor, estimate = .pred_class)
sens_no_pls     <- sens(preds_no_pls, truth = adh_yes_factor, estimate = .pred_class)
spec_no_pls     <- spec(preds_no_pls, truth = adh_yes_factor, estimate = .pred_class)

# PLS + baseline
accuracy_pls <- accuracy(preds_pls, truth = adh_yes_factor, estimate = .pred_class)
sens_pls     <- sens(preds_pls, truth = adh_yes_factor, estimate = .pred_class)
spec_pls     <- spec(preds_pls, truth = adh_yes_factor, estimate = .pred_class)

tibble(
  model = c("No PLS", "PLS + Baseline"),
  accuracy = c(accuracy_no_pls$.estimate, accuracy_pls$.estimate),
  sensitivity = c(sens_no_pls$.estimate, sens_pls$.estimate),
  specificity = c(spec_no_pls$.estimate, spec_pls$.estimate)
)

##
true <- preds_pls$adh_yes_factor
pred <- preds_pls$.pred_class

accuracy_value <- mean(pred == true)  # accuracy as proportion
n <- length(true)
# 95% CI using exact (Clopper-Pearson) method
accuracy_ci <- binom.confint(sum(pred == true), n, methods = "exact")
accuracy_value
accuracy_ci

# True positives
TP <- sum(pred == "Yes" & true == "Yes")
FN <- sum(pred == "No"  & true == "Yes")
sens_ci <- binom.confint(TP, TP + FN, methods = "exact")
sens_ci

TN <- sum(pred == "No" & true == "No")
FP <- sum(pred == "Yes" & true == "No")
spec_ci <- binom.confint(TN, TN + FP, methods = "exact")
spec_ci

tibble(
  metric = c("Accuracy", "Sensitivity", "Specificity"),
  estimate = c(mean(pred == true), TP / (TP + FN), TN / (TN + FP)),
  lower_95 = c(accuracy_ci$lower, sens_ci$lower, spec_ci$lower),
  upper_95 = c(accuracy_ci$upper, sens_ci$upper, spec_ci$upper)
)

# Example
true <- preds_no_pls$adh_yes_factor
pred <- preds_no_pls$.pred_class

accuracy_value <- mean(pred == true)  # accuracy as proportion
n <- length(true)
# 95% CI using exact (Clopper-Pearson) method
accuracy_ci <- binom.confint(sum(pred == true), n, methods = "exact")
accuracy_value
accuracy_ci

# True positives
TP <- sum(pred == "Yes" & true == "Yes")
FN <- sum(pred == "No"  & true == "Yes")
sens_ci <- binom.confint(TP, TP + FN, methods = "exact")
sens_ci

TN <- sum(pred == "No" & true == "No")
FP <- sum(pred == "Yes" & true == "No")
spec_ci <- binom.confint(TN, TN + FP, methods = "exact")
spec_ci

tibble(
  metric = c("Accuracy", "Sensitivity", "Specificity"),
  estimate = c(mean(pred == true), TP / (TP + FN), TN / (TN + FP)),
  lower_95 = c(accuracy_ci$lower, sens_ci$lower, spec_ci$lower),
  upper_95 = c(accuracy_ci$upper, sens_ci$upper, spec_ci$upper)
)

#######

ggplot(results, aes(x = model, y = mean, fill = model)) +
  geom_col(position = "dodge", alpha = 0.7) +
  geom_errorbar(aes(ymin = mean - std_err, ymax = mean + std_err), 
                width = 0.2, position = position_dodge(0.7)) +
  facet_wrap(~ .metric, scales = "free_y") +
  labs(
    title = "CV Performance Comparison: No PLS vs PLS + Baseline",
    x = "Model",
    y = "Mean CV Metric"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

#####
fit_tree <- fit(wf_no_pls, data = df_imputed)
tree_model <- extract_fit_parsnip(fit_tree)$fit
vip(tree_model)  # plots importance

fit_tree <- fit(wf_pls, data = df_imputed)
tree_model <- extract_fit_parsnip(fit_tree)$fit
vip(tree_model)  # plots importance

###
log_spec <- logistic_reg() %>%
  set_engine("glm") %>%
  set_mode("classification")

# No PLS (baseline predictors only)
wf_no_pls <- workflow() %>%
  add_model(log_spec) %>%
  add_recipe(recipe_no_pls)

# PLS + baseline
wf_pls <- workflow() %>%
  add_model(log_spec) %>%
  add_recipe(recipe_pls)

res_no_pls <- fit_resamples(
  wf_no_pls,
  resamples = outer_folds,
  metrics = metric_set(roc_auc, accuracy),
  control = control_resamples(save_pred = TRUE)
)

res_pls <- fit_resamples(
  wf_pls,
  resamples = outer_folds,
  metrics = metric_set(roc_auc, accuracy),
  control = control_resamples(save_pred = TRUE)
)

results_log <- bind_rows(
  collect_metrics(res_no_pls) %>% mutate(model = "No PLS"),
  collect_metrics(res_pls) %>% mutate(model = "PLS + Baseline")
)

results_log


###
# Suppose res_pls and res_no_pls were run with save_pred = TRUE
preds_pls <- collect_predictions(res_pls) %>% arrange(.row) 
preds_no_pls <- collect_predictions(res_no_pls) %>% arrange(.row)

roc_pls <- roc(preds_pls$adh_yes_factor, preds_pls$.pred_Yes)
roc_no_pls <- roc(preds_no_pls$adh_yes_factor, preds_no_pls$.pred_Yes)
roc.test(roc_pls, roc_no_pls, method = "delong")

# Default is 95% CI
ci_pls <- ci.auc(roc_pls)
ci_no_pls <- ci.auc(roc_no_pls)

ci_pls
ci_no_pls

# No PLS
conf_mat_no_pls <- conf_mat(preds_no_pls, truth = adh_yes_factor, estimate = .pred_class)
conf_mat_no_pls

# PLS + baseline
conf_mat_pls <- conf_mat(preds_pls, truth = adh_yes_factor, estimate = .pred_class)
conf_mat_pls

# No PLS
accuracy_no_pls <- accuracy(preds_no_pls, truth = adh_yes_factor, estimate = .pred_class)
sens_no_pls     <- sens(preds_no_pls, truth = adh_yes_factor, estimate = .pred_class)
spec_no_pls     <- spec(preds_no_pls, truth = adh_yes_factor, estimate = .pred_class)

# PLS + baseline
accuracy_pls <- accuracy(preds_pls, truth = adh_yes_factor, estimate = .pred_class)
sens_pls     <- sens(preds_pls, truth = adh_yes_factor, estimate = .pred_class)
spec_pls     <- spec(preds_pls, truth = adh_yes_factor, estimate = .pred_class)

tibble(
  model = c("No PLS", "PLS + Baseline"),
  accuracy = c(accuracy_no_pls$.estimate, accuracy_pls$.estimate),
  sensitivity = c(sens_no_pls$.estimate, sens_pls$.estimate),
  specificity = c(spec_no_pls$.estimate, spec_pls$.estimate)
)

# Example
true <- preds_pls$adh_yes_factor
pred <- preds_pls$.pred_class

accuracy_value <- mean(pred == true)  # accuracy as proportion
n <- length(true)
# 95% CI using exact (Clopper-Pearson) method
accuracy_ci <- binom.confint(sum(pred == true), n, methods = "exact")
accuracy_value
accuracy_ci

# True positives
TP <- sum(pred == "Yes" & true == "Yes")
FN <- sum(pred == "No"  & true == "Yes")
sens_ci <- binom.confint(TP, TP + FN, methods = "exact")
sens_ci

TN <- sum(pred == "No" & true == "No")
FP <- sum(pred == "Yes" & true == "No")
spec_ci <- binom.confint(TN, TN + FP, methods = "exact")
spec_ci

tibble(
  metric = c("Accuracy", "Sensitivity", "Specificity"),
  estimate = c(mean(pred == true), TP / (TP + FN), TN / (TN + FP)),
  lower_95 = c(accuracy_ci$lower, sens_ci$lower, spec_ci$lower),
  upper_95 = c(accuracy_ci$upper, sens_ci$upper, spec_ci$upper)
)

# Example
true <- preds_no_pls$adh_yes_factor
pred <- preds_no_pls$.pred_class

accuracy_value <- mean(pred == true)  # accuracy as proportion
n <- length(true)
# 95% CI using exact (Clopper-Pearson) method
accuracy_ci <- binom.confint(sum(pred == true), n, methods = "exact")
accuracy_value
accuracy_ci

# True positives
TP <- sum(pred == "Yes" & true == "Yes")
FN <- sum(pred == "No"  & true == "Yes")
sens_ci <- binom.confint(TP, TP + FN, methods = "exact")
sens_ci

TN <- sum(pred == "No" & true == "No")
FP <- sum(pred == "Yes" & true == "No")
spec_ci <- binom.confint(TN, TN + FP, methods = "exact")
spec_ci

tibble(
  metric = c("Accuracy", "Sensitivity", "Specificity"),
  estimate = c(mean(pred == true), TP / (TP + FN), TN / (TN + FP)),
  lower_95 = c(accuracy_ci$lower, sens_ci$lower, spec_ci$lower),
  upper_95 = c(accuracy_ci$upper, sens_ci$upper, spec_ci$upper)
)

fit_tree <- fit(wf_no_pls, data = df_imputed)
tree_model <- extract_fit_parsnip(fit_tree)$fit
vip(tree_model)  # plots importance


fit_tree <- fit(wf_pls, data = df_imputed)
tree_model <- extract_fit_parsnip(fit_tree)$fit
vip(tree_model)  # plots importance


############
svm_spec <- svm_linear() %>%
  set_engine("kernlab") %>%   # linear SVM engine
  set_mode("classification")

wf_no_pls_svm <- workflow() %>%
  add_model(svm_spec) %>%
  add_recipe(recipe_no_pls)

# PLS + baseline
wf_pls_svm <- workflow() %>%
  add_model(svm_spec) %>%
  add_recipe(recipe_pls)

res_no_pls_svm <- fit_resamples(
  wf_no_pls_svm,
  resamples = outer_folds,
  metrics = metric_set(roc_auc, accuracy),
  control = control_resamples(save_pred = TRUE, save_workflow = TRUE)
  # important for later ROC/Delong
)

res_pls_svm <- fit_resamples(
  wf_pls_svm,
  resamples = outer_folds,
  metrics = metric_set(roc_auc, accuracy),
  control = control_resamples(save_pred = TRUE, save_workflow = TRUE)
)

preds_no_pls_svm <- collect_predictions(res_no_pls_svm) %>% arrange(.row)
preds_pls_svm <- collect_predictions(res_pls_svm) %>% arrange(.row)

roc_no_pls_svm <- roc(preds_no_pls_svm$adh_yes_factor, preds_no_pls_svm$.pred_Yes)
roc_pls_svm    <- roc(preds_pls_svm$adh_yes_factor, preds_pls_svm$.pred_Yes)

roc.test(roc_no_pls_svm, roc_pls_svm, method = "delong")

# Default is 95% CI
ci_pls <- ci.auc(roc_pls_svm)
ci_no_pls <- ci.auc(roc_no_pls_svm)

ci_pls
ci_no_pls

# No PLS
conf_mat_no_pls <- conf_mat(preds_no_pls_svm, truth = adh_yes_factor, estimate = .pred_class)
conf_mat_no_pls

# PLS + baseline
conf_mat_pls <- conf_mat(preds_pls_svm, truth = adh_yes_factor, estimate = .pred_class)
conf_mat_pls

# No PLS
accuracy_no_pls <- accuracy(preds_no_pls_svm, truth = adh_yes_factor, estimate = .pred_class)
sens_no_pls     <- sens(preds_no_pls_svm, truth = adh_yes_factor, estimate = .pred_class)
spec_no_pls     <- spec(preds_no_pls_svm, truth = adh_yes_factor, estimate = .pred_class)

# PLS + baseline
accuracy_pls <- accuracy(preds_pls_svm, truth = adh_yes_factor, estimate = .pred_class)
sens_pls     <- sens(preds_pls_svm, truth = adh_yes_factor, estimate = .pred_class)
spec_pls     <- spec(preds_pls_svm, truth = adh_yes_factor, estimate = .pred_class)

tibble(
  model = c("No PLS", "PLS + Baseline"),
  accuracy = c(accuracy_no_pls$.estimate, accuracy_pls$.estimate),
  sensitivity = c(sens_no_pls$.estimate, sens_pls$.estimate),
  specificity = c(spec_no_pls$.estimate, spec_pls$.estimate)
)

# Example
true <- preds_no_pls_svm$adh_yes_factor
pred <- preds_no_pls_svm$.pred_class

accuracy_value <- mean(pred == true)  # accuracy as proportion
n <- length(true)
# 95% CI using exact (Clopper-Pearson) method
accuracy_ci <- binom.confint(sum(pred == true), n, methods = "exact")
accuracy_value
accuracy_ci

# True positives
TP <- sum(pred == "Yes" & true == "Yes")
FN <- sum(pred == "No"  & true == "Yes")
sens_ci <- binom.confint(TP, TP + FN, methods = "exact")
sens_ci

TN <- sum(pred == "No" & true == "No")
FP <- sum(pred == "Yes" & true == "No")
spec_ci <- binom.confint(TN, TN + FP, methods = "exact")
spec_ci

tibble(
  metric = c("Accuracy", "Sensitivity", "Specificity"),
  estimate = c(mean(pred == true), TP / (TP + FN), TN / (TN + FP)),
  lower_95 = c(accuracy_ci$lower, sens_ci$lower, spec_ci$lower),
  upper_95 = c(accuracy_ci$upper, sens_ci$upper, spec_ci$upper)
)


# Example
true <- preds_pls_svm$adh_yes_factor
pred <- preds_pls_svm$.pred_class

accuracy_value <- mean(pred == true)  # accuracy as proportion
n <- length(true)
# 95% CI using exact (Clopper-Pearson) method
accuracy_ci <- binom.confint(sum(pred == true), n, methods = "exact")
accuracy_value
accuracy_ci

# True positives
TP <- sum(pred == "Yes" & true == "Yes")
FN <- sum(pred == "No"  & true == "Yes")
sens_ci <- binom.confint(TP, TP + FN, methods = "exact")
sens_ci

TN <- sum(pred == "No" & true == "No")
FP <- sum(pred == "Yes" & true == "No")
spec_ci <- binom.confint(TN, TN + FP, methods = "exact")
spec_ci

tibble(
  metric = c("Accuracy", "Sensitivity", "Specificity"),
  estimate = c(mean(pred == true), TP / (TP + FN), TN / (TN + FP)),
  lower_95 = c(accuracy_ci$lower, sens_ci$lower, spec_ci$lower),
  upper_95 = c(accuracy_ci$upper, sens_ci$upper, spec_ci$upper)
)

# Fit the workflow on training data (or use a resampled workflow if you want CV-aware importance)
final_fit <- fit(wf_pls_svm, data = df_imputed)
svm_fit <- extract_fit_parsnip(final_fit)
ksvm_model <- svm_fit$fit  # this is an S4 object


library(parsnip)
install.packages("LiblineaR")
library(LiblineaR)
svm_model <- svm_linear(mode = "classification") %>%
  set_engine("LiblineaR")
wf_pls_svm <- wf_pls_svm %>% 
  update_model(svm_model)
final_fit <- fit(wf_pls_svm, data = df_imputed)
svm_fit <- extract_fit_parsnip(final_fit)
importance <- abs(svm_fit$fit$W)   # numeric vector
importance

# Drop the bias column
importance_no_bias <- importance[, colnames(importance) != "Bias", drop = FALSE]

importance_df <- data.frame(
  variable = colnames(importance_no_bias),
  importance = as.numeric(importance_no_bias[1, ])
) %>%
  arrange(desc(importance_no_bias))

ggplot(importance_df, aes(x = reorder(variable, importance), y = importance)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Variable Importance (Linear SVM)",
    x = "Variable",
    y = "Absolute Coefficient"
  ) +
  theme_minimal()


svm_model <- svm_linear(mode = "classification") %>%
  set_engine("LiblineaR")
wf_no_pls_svm <- wf_no_pls_svm %>% 
  update_model(svm_model)
final_fit <- fit(wf_no_pls_svm, data = df_imputed)
svm_fit <- extract_fit_parsnip(final_fit)
importance <- abs(svm_fit$fit$W)   # numeric vector
importance

# Drop the bias column
importance_no_bias <- importance[, colnames(importance) != "Bias", drop = FALSE]

importance_df <- data.frame(
  variable = colnames(importance_no_bias),
  importance = as.numeric(importance_no_bias[1, ])
) %>%
  arrange(desc(importance_no_bias))

ggplot(importance_df, aes(x = reorder(variable, importance), y = importance)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Variable Importance (Linear SVM)",
    x = "Variable",
    y = "Absolute Coefficient"
  ) +
  theme_minimal()

####
library(parsnip)
library(workflows)
library(tidymodels)
install.packages("discrim")  # only if not installed
library(discrim)

nb_spec <- naive_Bayes() %>%
  set_engine("klaR") %>%   # common engine for NB in R
  set_mode("classification")

wf_no_pls_nb <- workflow() %>%
  add_model(nb_spec) %>%
  add_recipe(recipe_no_pls)

# PLS + baseline
wf_pls_nb <- workflow() %>%
  add_model(nb_spec) %>%
  add_recipe(recipe_pls)

res_no_pls_nb <- fit_resamples(
  wf_no_pls_nb,
  resamples = outer_folds,
  metrics = metric_set(roc_auc, accuracy),
  control = control_resamples(save_pred = TRUE)
)

res_pls_nb <- fit_resamples(
  wf_pls_nb,
  resamples = outer_folds,
  metrics = metric_set(roc_auc, accuracy),
  control = control_resamples(save_pred = TRUE)
)

preds_no_pls_nb <- collect_predictions(res_no_pls_nb) %>% arrange(.row)
preds_pls_nb    <- collect_predictions(res_pls_nb) %>% arrange(.row)


roc_no_pls_nb <- roc(preds_no_pls_nb$adh_yes_factor, preds_no_pls_nb$.pred_Yes)
roc_pls_nb    <- roc(preds_pls_nb$adh_yes_factor, preds_pls_nb$.pred_Yes)

# Extract numeric AUCs
auc1 <- auc(roc_pls_nb)
auc2 <- auc(roc_no_pls_nb)

# Delong's test
roc.test(roc_no_pls_nb, roc_pls_nb, method = "delong")

# Default is 95% CI
ci_pls <- ci.auc(roc_pls)
ci_no_pls <- ci.auc(roc_no_pls)

ci_pls
ci_no_pls


# No PLS
conf_mat_no_pls <- conf_mat(preds_no_pls_nb, truth = adh_yes_factor, estimate = .pred_class)
conf_mat_no_pls

# PLS + baseline
conf_mat_pls <- conf_mat(preds_pls_nb, truth = adh_yes_factor, estimate = .pred_class)
conf_mat_pls

# No PLS
accuracy_no_pls <- accuracy(preds_no_pls_nb, truth = adh_yes_factor, estimate = .pred_class)
sens_no_pls     <- sens(preds_no_pls_nb, truth = adh_yes_factor, estimate = .pred_class)
spec_no_pls     <- spec(preds_no_pls_nb, truth = adh_yes_factor, estimate = .pred_class)

# PLS + baseline
accuracy_pls <- accuracy(preds_pls_nb, truth = adh_yes_factor, estimate = .pred_class)
sens_pls     <- sens(preds_pls_nb, truth = adh_yes_factor, estimate = .pred_class)
spec_pls     <- spec(preds_pls_nb, truth = adh_yes_factor, estimate = .pred_class)

tibble(
  model = c("No PLS", "PLS + Baseline"),
  accuracy = c(accuracy_no_pls$.estimate, accuracy_pls$.estimate),
  sensitivity = c(sens_no_pls$.estimate, sens_pls$.estimate),
  specificity = c(spec_no_pls$.estimate, spec_pls$.estimate)
)


#############
install.packages("binom")  # if not installed
library(binom)

# Example
true <- preds_pls_nb$adh_yes_factor
pred <- preds_pls_nb$.pred_class

accuracy_value <- mean(pred == true)  # accuracy as proportion
n <- length(true)
# 95% CI using exact (Clopper-Pearson) method
accuracy_ci <- binom.confint(sum(pred == true), n, methods = "exact")
accuracy_value
accuracy_ci

# True positives
TP <- sum(pred == "Yes" & true == "Yes")
FN <- sum(pred == "No"  & true == "Yes")
sens_ci <- binom.confint(TP, TP + FN, methods = "exact")
sens_ci

TN <- sum(pred == "No" & true == "No")
FP <- sum(pred == "Yes" & true == "No")
spec_ci <- binom.confint(TN, TN + FP, methods = "exact")
spec_ci

tibble(
  metric = c("Accuracy", "Sensitivity", "Specificity"),
  estimate = c(mean(pred == true), TP / (TP + FN), TN / (TN + FP)),
  lower_95 = c(accuracy_ci$lower, sens_ci$lower, spec_ci$lower),
  upper_95 = c(accuracy_ci$upper, sens_ci$upper, spec_ci$upper)
)


# Example
true <- preds_no_pls_nb$adh_yes_factor
pred <- preds_no_pls_nb$.pred_class

accuracy_value <- mean(pred == true)  # accuracy as proportion
n <- length(true)
# 95% CI using exact (Clopper-Pearson) method
accuracy_ci <- binom.confint(sum(pred == true), n, methods = "exact")
accuracy_value
accuracy_ci

# True positives
TP <- sum(pred == "Yes" & true == "Yes")
FN <- sum(pred == "No"  & true == "Yes")
sens_ci <- binom.confint(TP, TP + FN, methods = "exact")
sens_ci

TN <- sum(pred == "No" & true == "No")
FP <- sum(pred == "Yes" & true == "No")
spec_ci <- binom.confint(TN, TN + FP, methods = "exact")
spec_ci

tibble(
  metric = c("Accuracy", "Sensitivity", "Specificity"),
  estimate = c(mean(pred == true), TP / (TP + FN), TN / (TN + FP)),
  lower_95 = c(accuracy_ci$lower, sens_ci$lower, spec_ci$lower),
  upper_95 = c(accuracy_ci$upper, sens_ci$upper, spec_ci$upper)
)

library(vip)

pred_wrapper <- function(object, newdata) {
  predict(object, new_data = newdata, type = "prob")$.pred_Yes
}

final_fit_no_pls_nb <- fit(wf_no_pls_nb, data = df_imputed)
baseline_features <- names(df_imputed)[names(df_imputed) %in% c("mmrc_a", "charlson_a", "current_smoker_num_a", "age2_a", "sex_number_a", "fev1_updated_a")] 
pls_features <- c(baseline_features, "PLS_1")  # only the PLS component, not the original variables

vi_baseline <- vi_permute(
  object = final_fit_no_pls_nb,
  feature_names = baseline_features,
  target = df_imputed$adh_yes_factor,
  metric = "roc_auc",
  pred_wrapper = pred_wrapper,
  nsim = 10,
  train = df_imputed
)

vi_baseline$Importance <- -vi_baseline$Importance  # flip sign so higher = more important
vi_baseline$Importance <- pmax(vi_baseline$Importance, 0)  # optional: remove small negatives

ggplot(vi_baseline, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Permutation Variable Importance (PLS + Baseline)",
       x = "Variable",
       y = "Importance") +
  theme_minimal()

final_fit_pls_nb <- fit(wf_pls_nb, data = df_imputed)

vi_pls <- vi_permute(
  object = final_fit_pls_nb,
  feature_names = pls_features,
  target = df_imputed$adh_yes_factor,
  metric = "roc_auc",
  pred_wrapper = pred_wrapper,
  nsim = 10,
  train = df_imputed
)

vi_pls$Importance <- -vi_pls$Importance  # flip sign so higher = more important
vi_pls$Importance <- pmax(vi_pls$Importance, 0)  # optional: remove small negatives

ggplot(vi_pls, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Permutation Variable Importance (PLS + Baseline)",
       x = "Variable",
       y = "Importance") +
  theme_minimal()
