# ==============================================================================
# 03_evaluate_model.R
# Evaluation tables and figures for the fitted biomod2 models.
#
#   - AUC (ROC) summaries of the single models (cross-validation)
#   - Ensemble model AUC/TSS
#   - Variable importance (single models and ensemble)
#   - Ensemble (EMwmean) response curves
#   - AUC bar chart, ROC curves and combined figure
#
# Inputs : biomod2/<crop>_biomod_models.rds, biomod2/<crop>_biomod_ensemble.rds
#          (written by 02_fit_ensemble_model.R)
# Outputs: tables/ and figures/ under outputs/<crop>/
#
# Author : Rapheal Wangalwa
# ==============================================================================

source(here::here("R", "00_config.R"))
source(here::here("R", "utils.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(pROC)
  library(biomod2)
})

cat("Crop:", CROP_LABEL, "\n")

bm_models   <- readRDS(FILES$model_out)
bm_ensemble <- readRDS(FILES$ensemble_out)

# ==============================================================================
# 1. Single-model AUC (validation folds)
# ==============================================================================
model_eval <- get_evaluations(bm_models)
str(model_eval)

auc_summary <- model_eval %>%
  filter(metric.eval == "ROC") %>%
  group_by(algo) %>%
  summarise(
    Avg_AUC = round(mean(validation, na.rm = TRUE), 3),
    SD_AUC  = round(sd(validation,   na.rm = TRUE), 3),
    Min_AUC = round(min(validation,  na.rm = TRUE), 3),
    Max_AUC = round(max(validation,  na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(desc(Avg_AUC)) %>%
  rename(Model = algo)

print(auc_summary)
write_table(auc_summary, out_table("model_auc_summary.csv"))

# ==============================================================================
# 2. Ensemble model performance
# ==============================================================================
em_eval <- get_evaluations(bm_ensemble)
print(em_eval)

ensemble_scores <- em_eval %>%
  filter(metric.eval %in% c("ROC", "TSS")) %>%
  select(algo, metric.eval, calibration)
print(ensemble_scores)
write_table(ensemble_scores, out_table("ensemble_evaluation.csv"))

# ==============================================================================
# 3. Variable importance
# ==============================================================================
# 3a. Single models
var_imp_summary <- get_variables_importance(bm_models) %>%
  group_by(algo, expl.var) %>%
  summarise(
    Mean_Importance = mean(var.imp, na.rm = TRUE),
    SD_Importance   = sd(var.imp,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean_Importance))
print(var_imp_summary)
write_table(var_imp_summary, out_table("single_model_variable_importance.csv"))

p_var_imp_models <- ggplot(var_imp_summary,
                           aes(x = reorder(expl.var, Mean_Importance),
                               y = Mean_Importance, fill = algo)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(title = "Variable Importance by Algorithm",
       x = "Environmental Variable", y = "Mean Importance") +
  theme_bw()
print(p_var_imp_models)
save_plot(p_var_imp_models, out_figure("single_model_variable_importance"),
          width = 8, height = 6)

# 3b. Ensemble (averaged over EMmean and EMwmean), expressed as % of total
ensemble_var_imp_summary <- get_variables_importance(bm_ensemble) %>%
  group_by(expl.var) %>%
  summarise(
    Mean_Importance = mean(var.imp, na.rm = TRUE),
    SD_Importance   = sd(var.imp,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean_Importance)) %>%
  mutate(Importance_Percent = 100 * Mean_Importance / sum(Mean_Importance, na.rm = TRUE))

print(ensemble_var_imp_summary)
write_table(ensemble_var_imp_summary, out_table("ensemble_variable_importance.csv"))
saveRDS(ensemble_var_imp_summary, out_table("ensemble_variable_importance.rds"))

p_var_imp <- ggplot(ensemble_var_imp_summary,
                    aes(x = reorder(expl.var, Importance_Percent),
                        y = Importance_Percent, fill = Importance_Percent)) +
  geom_bar(stat = "identity", alpha = 1, color = "black", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f", Importance_Percent)),
            hjust = -0.2, color = "black", size = 3.5) +
  scale_fill_gradient(low = "orange", high = "darkblue", name = "Importance (%)") +
  coord_flip() +
  # NOTE: bars above 20 % would be dropped by this axis limit
  scale_y_continuous(limits = c(0, 20), expand = c(0, 0)) +
  labs(x = "Environmental Variables", y = "Importance (%)") +
  theme_bw() +
  theme(
    legend.position    = "none",
    axis.text.x        = element_text(color = "black", face = "bold", size = 9),
    axis.text.y        = element_text(color = "black", face = "bold", size = 9),
    axis.title         = element_text(size = 12, face = "bold"),
    plot.title         = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border       = element_rect(color = "black", linewidth = 1.0),
    axis.ticks         = element_line(color = "black", linewidth = 0.8)
  )
print(p_var_imp)
save_plot(p_var_imp, out_figure("ensemble_variable_importance"),
          width = 8, height = 6, pdf_height = 9)

# ==============================================================================
# 4. Ensemble response curves (EMwmean)
# ==============================================================================
response_data <- bm_PlotResponseCurves(
  bm.out        = bm_ensemble,
  models.chosen = get_built_models(bm_ensemble),
  fixed.var     = "median",
  do.bivariate  = FALSE,
  do.plot       = FALSE
)

emwmean_data <- response_data$tab %>%
  filter(grepl("EMwmean", pred.name))
saveRDS(emwmean_data, out_table("ensemble_response_curves_EMwmean.rds"))

p_response <- ggplot(emwmean_data, aes(x = expl.val, y = pred.val)) +
  geom_line(linewidth = 1, color = "darkred") +
  facet_wrap(~ expl.name, scales = "free") +
  labs(x = "Environmental Variable Value", y = "Predicted Probability",
       title = "Response Curves") +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 10),
    axis.text        = element_text(size = 8),
    axis.title       = element_text(size = 10),
    plot.title       = element_text(face = "bold", hjust = 0.5)
  )
print(p_response)
save_plot(p_response, out_figure("ensemble_response_curves_EMwmean"),
          width = 12, height = 10)

# ==============================================================================
# 5. AUC bar chart (calibration vs validation, mean +/- SD)
# ==============================================================================
auc_long <- model_eval %>%
  filter(metric.eval == "ROC") %>%
  mutate(dataset = ifelse(grepl("allData", full.name), "allData", "subset")) %>%
  pivot_longer(cols = c(calibration, validation),
               names_to = "evaluation_type", values_to = "auc_score") %>%
  filter(!is.na(auc_score))

auc_summary_detailed <- auc_long %>%
  group_by(algo, dataset, evaluation_type) %>%
  summarise(
    mean_auc = mean(auc_score, na.rm = TRUE),
    sd_auc   = sd(auc_score,   na.rm = TRUE),
    n        = n(),
    se_auc   = sd_auc / sqrt(n),
    .groups  = "drop"
  )
print(auc_summary_detailed)
write_table(auc_summary_detailed, out_table("model_auc_summary_detailed.csv"))

auc_bar_plot <- ggplot(auc_summary_detailed,
                       aes(x = algo, y = mean_auc, fill = evaluation_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           alpha = 1, width = 0.7, color = "black", linewidth = 0.5) +
  geom_errorbar(aes(ymin = mean_auc - sd_auc, ymax = mean_auc + sd_auc),
                position = position_dodge(width = 0.8),
                width = 0.25, linewidth = 0.8, color = "black") +
  scale_fill_manual(values = c(calibration = "blue", validation = "orange")) +
  labs(x = "Algorithm", y = "AUCroc Score", fill = "Evaluation Type") +
  theme_bw() +
  theme(
    axis.text.x        = element_text(color = "black", face = "bold"),
    axis.text.y        = element_text(color = "black", face = "bold"),
    legend.position    = "bottom",
    legend.text        = element_text(size = 12),
    legend.title       = element_text(size = 13, face = "bold"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x       = element_text(face = "bold", size = 12),
    axis.title.y       = element_text(face = "bold", size = 12),
    panel.border       = element_rect(color = "black", linewidth = 1.5),
    axis.ticks         = element_line(color = "black", linewidth = 0.8)
  ) +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "red",
             linewidth = 1, alpha = 0.7) +
  annotate("text", x = 0.5, y = 0.72, label = "Good (0.7)",
           color = "red", hjust = 0, size = 4) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "darkblue",
             linewidth = 1, alpha = 0.7) +
  annotate("text", x = 0.5, y = 0.82, label = "Very Good (0.8)",
           color = "darkblue", hjust = 0, size = 4) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "green",
             linewidth = 1, alpha = 0.7) +
  annotate("text", x = 0.5, y = 0.92, label = "Excellent (0.9)",
           color = "green", hjust = 0, size = 4) +
  geom_text(aes(label = round(mean_auc, 3), y = mean_auc + 0.02),
            position = position_dodge(width = 0.8),
            size = 3.5, fontface = "bold")
print(auc_bar_plot)
save_plot(auc_bar_plot, out_figure("auc_performance_bar"),
          width = 10, height = 9, pdf_height = 11, rds = TRUE)

# ==============================================================================
# 6. ROC curves per algorithm (pooled over all runs)
# ==============================================================================
# NOTE: these use the stored predictions on all points (calibration + validation),
# so the pooled AUC is not an independent validation score.
obs   <- get_formal_data(bm_models, "resp.var")  # observed 0/1 response
preds <- get_predictions(bm_models)              # columns: full.name, run, algo, points, pred

overall_aucs <- preds %>%
  group_by(algo) %>%
  summarise(auc_value = round(as.numeric(auc(roc(obs[points], pred, quiet = TRUE))), 3))
print(overall_aucs)
write_table(overall_aucs, out_table("model_overall_aucs.csv"))

# One ROC object per algorithm; names taken from the grouping keys so they
# always match the order used by group_map()
preds_by_algo <- group_by(preds, algo)
roc_list_by_model <- group_map(preds_by_algo,
                               ~ roc(obs[.x$points], .x$pred, quiet = TRUE))
names(roc_list_by_model) <- as.character(group_keys(preds_by_algo)$algo)

roc_to_df <- function(roc_obj, algo, auc_val, delta = 0.02) {
  sens <- rev(roc_obj$sensitivities)
  spec <- rev(roc_obj$specificities)
  data.frame(
    specificity = spec,
    sensitivity = sens,
    ymin  = pmax(sens - delta, 0),
    ymax  = pmin(sens + delta, 1),
    model = paste0(algo, " (AUC=", auc_val, ")")
  )
}

roc_df <- do.call(rbind, lapply(names(roc_list_by_model), function(name) {
  auc_val <- overall_aucs$auc_value[overall_aucs$algo == name]
  roc_to_df(roc_list_by_model[[name]], name, auc_val)
}))

roc_plot <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity,
                               color = model, fill = model)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.3, colour = NA) +
  geom_step(linewidth = 0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("red", "green", "blue")) +
  scale_color_manual(values = c("red", "green", "blue")) +
  labs(x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.6, "cm"),
    axis.title.x    = element_text(face = "bold", size = 13),
    axis.title.y    = element_text(face = "bold", size = 13),
    panel.border    = element_rect(color = "black", linewidth = 1.5),
    axis.ticks      = element_line(color = "black", linewidth = 0.8),
    axis.text.x     = element_text(color = "black", face = "bold"),
    axis.text.y     = element_text(color = "black", face = "bold")
  )
print(roc_plot)
save_plot(roc_plot, out_figure("roc_curves_overall_auc"),
          width = 10, height = 8, rds = TRUE)

# ==============================================================================
# 7. Combined AUC + ROC figure
# ==============================================================================
combined_plot <- auc_bar_plot + roc_plot
print(combined_plot)
save_plot(combined_plot, out_figure("combined_auc_roc"),
          width = 14, height = 8, rds = TRUE)
