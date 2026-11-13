#!/usr/bin/env Rscript
# ==============================================================================
# Unified GAMM Modeling Suite on Final Cohort (N = 4,377 Windows, n = 1,099 Users)
# Ensures 100% consistency across all RDS files, Figures 4-8, and Tables 1-9.
# ==============================================================================

suppressPackageStartupMessages({
  library(mgcv)
  library(dplyr)
})

cat(">>> [Step 1] Loading and verifying final N = 4,377 cohort...\n")
raw_path <- "/Users/speediance/Desktop/paper数据处理/windows_unrestricted_symmetric_y20.csv"
df <- read.csv(raw_path, stringsAsFactors = FALSE) %>%
  filter(!is.na(train_avg_pct_1rm),
         !is.na(y_delta_1rm_pct),
         !is.na(pre_1rm),
         !is.na(weight_kg),
         !is.na(user_age_at_training),
         !is.na(gender_male)) %>%
  mutate(
    train_total_squat_sets = train_num_sessions * train_avg_squat_sets,
    user_id_fac = factor(user_id),
    gender_fac = factor(ifelse(gender_male == 1, "Male", "Female"), levels = c("Male", "Female")),
    age_group = factor(ifelse(user_age_at_training <= 35, "<=35", ">35"), levels = c("<=35", ">35")),
    strength_ratio = pre_1rm / weight_kg,
    strength_level_fac = factor(
      case_when(
        strength_ratio <= 1.0 ~ "Novice",
        strength_ratio <= 1.5 ~ "Intermediate",
        TRUE ~ "Advanced"
      ),
      levels = c("Novice", "Intermediate", "Advanced")
    )
  )

stopifnot(nrow(df) == 4377)
stopifnot(n_distinct(df$user_id) == 1099)
write.csv(df, "/Users/speediance/Desktop/paper数据处理/cleaned_windows_4377.csv", row.names = FALSE)

# Helper to extract clean summary table from a bam model
extract_model_table <- function(m, group_name) {
  s_tab <- summary(m)$s.table
  p_tab <- summary(m)$p.table
  vc <- tryCatch(gam.vcomp(m), error = function(e) NULL)
  re_var <- NA_real_
  if (!is.null(vc) && "s(user_id_fac)" %in% rownames(vc)) {
    re_var <- vc["s(user_id_fac)", "std.dev"]^2
  }
  
  rows <- list()
  for (rn in rownames(s_tab)) {
    rows[[length(rows) + 1]] <- data.frame(
      Group = group_name,
      Type = ifelse(rn == "s(user_id_fac)", "Random Effect", "Smooth Term"),
      Term = rn,
      EDF = s_tab[rn, "edf"],
      Ref_df = s_tab[rn, "Ref.df"],
      Stat_Type = ifelse(rn == "s(user_id_fac)", "sigma2", "F"),
      Stat_Val = ifelse(rn == "s(user_id_fac)", re_var, s_tab[rn, "F"]),
      F_Val = s_tab[rn, "F"],
      p_val = s_tab[rn, "p-value"],
      stringsAsFactors = FALSE
    )
  }
  if ("gender_male" %in% rownames(p_tab)) {
    rows[[length(rows) + 1]] <- data.frame(
      Group = group_name,
      Type = "Parametric Fixed",
      Term = "gender_male",
      EDF = 1.0,
      Ref_df = 1.0,
      Stat_Type = "t",
      Stat_Val = p_tab["gender_male", "t value"],
      F_Val = p_tab["gender_male", "Estimate"],
      p_val = p_tab["gender_male", "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

# ==============================================================================
# 1. Figure 4: Single-Variable Dose-Response GAMMs (Sets & Intensity)
# ==============================================================================
cat(">>> [Step 2] Fitting Figure 4 single-variable GAMM models (N = 4,377)...\n")
m_sets <- bam(
  y_delta_1rm_pct ~ s(train_total_squat_sets, k = 10, bs = "tp") +
    s(pre_1rm, k = 10, bs = "tp") +
    s(weight_kg, k = 8, bs = "tp") +
    s(user_age_at_training, k = 6, bs = "tp") +
    gender_male +
    s(user_id_fac, bs = "re"),
  data = df, method = "fREML", discrete = TRUE
)

m_int <- bam(
  y_delta_1rm_pct ~ s(train_avg_pct_1rm, k = 10, bs = "tp") +
    s(pre_1rm, k = 10, bs = "tp") +
    s(weight_kg, k = 8, bs = "tp") +
    s(user_age_at_training, k = 6, bs = "tp") +
    gender_male +
    s(user_id_fac, bs = "re"),
  data = df, method = "fREML", discrete = TRUE
)

# Predictions for Figure 4 (trimmed to 1st-99th percentile of exposure for visual stability, fitted on full N=4,377)
sets_q <- quantile(df$train_total_squat_sets, c(0.01, 0.99))
sets_seq <- seq(sets_q[1], sets_q[2], length.out = 200)
nd_sets <- data.frame(
  train_total_squat_sets = sets_seq,
  pre_1rm = median(df$pre_1rm),
  weight_kg = median(df$weight_kg),
  user_age_at_training = median(df$user_age_at_training),
  gender_male = 1,
  user_id_fac = df$user_id_fac[1]
)
p_sets <- predict(m_sets, nd_sets, se.fit = TRUE, exclude = "s(user_id_fac)")
nd_sets$pred_fit <- p_sets$fit
nd_sets$pred_se <- p_sets$se.fit
nd_sets$ci_lower <- p_sets$fit - 1.96 * p_sets$se.fit
nd_sets$ci_upper <- p_sets$fit + 1.96 * p_sets$se.fit
nd_sets$dose_value <- sets_seq
nd_sets$dose_var <- "train_total_squat_sets"
nd_sets$dose_label <- "4-Week Total Squat Sets"

int_q <- quantile(df$train_avg_pct_1rm, c(0.01, 0.99))
int_seq_f4 <- seq(int_q[1], int_q[2], length.out = 200)
nd_int <- data.frame(
  train_avg_pct_1rm = int_seq_f4,
  pre_1rm = median(df$pre_1rm),
  weight_kg = median(df$weight_kg),
  user_age_at_training = median(df$user_age_at_training),
  gender_male = 1,
  user_id_fac = df$user_id_fac[1]
)
p_int <- predict(m_int, nd_int, se.fit = TRUE, exclude = "s(user_id_fac)")
nd_int$pred_fit <- p_int$fit
nd_int$pred_se <- p_int$se.fit
nd_int$ci_lower <- p_int$fit - 1.96 * p_int$se.fit
nd_int$ci_upper <- p_int$fit + 1.96 * p_int$se.fit
nd_int$dose_value <- int_seq_f4
nd_int$dose_var <- "train_avg_pct_1rm"
nd_int$dose_label <- "Average Training Intensity (%1RM)"

write.csv(nd_sets, "/Users/speediance/Desktop/paper数据处理/pred_fig4A_sets_4377.csv", row.names = FALSE)
write.csv(nd_int, "/Users/speediance/Desktop/paper数据处理/pred_fig4B_int_4377.csv", row.names = FALSE)

# Update gamm_models_results.rds
models_results <- list(
  total_sets = list(
    model = m_sets,
    predictions = nd_sets,
    raw_data = df[, c("user_id", "train_total_squat_sets", "y_delta_1rm_pct")],
    dose_var = "train_total_squat_sets",
    dose_label = "4-Week Total Squat Sets"
  ),
  intensity = list(
    model = m_int,
    predictions = nd_int,
    raw_data = df[, c("user_id", "train_avg_pct_1rm", "y_delta_1rm_pct")],
    dose_var = "train_avg_pct_1rm",
    dose_label = "Average Training Intensity (%1RM)"
  )
)
saveRDS(models_results, "/Users/speediance/Desktop/one_rm_research/output/gamm_models_results.rds")
saveRDS(models_results, "/Users/speediance/Desktop/paper数据处理/gamm_models_results.rds")

# ==============================================================================
# 2. Table 2 & Figure 5: Main Tensor Product Interaction Model
# ==============================================================================
cat(">>> [Step 3] Fitting Table 2 / Figure 5 Main Interaction GAMM (N = 4,377)...\n")
m_main_te <- bam(
  y_delta_1rm_pct ~ te(train_total_squat_sets, train_avg_pct_1rm, k = c(8, 8)) +
    s(pre_1rm, k = 10, bs = "tp") +
    s(weight_kg, k = 8, bs = "tp") +
    s(user_age_at_training, k = 6, bs = "tp") +
    gender_male +
    s(user_id_fac, bs = "re"),
  data = df, method = "fREML", discrete = TRUE
)

# 2D high-res grid for Figure 5 (v_mesh: 6..32, i_mesh: 52..82, 180x180)
v_mesh <- seq(6, 32, length.out = 180)
i_mesh <- seq(52, 82, length.out = 180)
grid_f5 <- expand.grid(train_total_squat_sets = v_mesh, train_avg_pct_1rm = i_mesh)
grid_f5$pre_1rm <- median(df$pre_1rm)
grid_f5$weight_kg <- median(df$weight_kg)
grid_f5$user_age_at_training <- median(df$user_age_at_training)
grid_f5$gender_male <- 1
grid_f5$user_id_fac <- df$user_id_fac[1]

p_f5 <- predict(m_main_te, grid_f5, se.fit = TRUE, exclude = "s(user_id_fac)")
grid_f5$pred_fit <- p_f5$fit
grid_f5$pred_se <- p_f5$se.fit
write.csv(grid_f5, "/Users/speediance/Desktop/paper数据处理/pred_fig5_surface_4377.csv", row.names = FALSE)

# Exact predictions on the +5.00% iso-gain contour: Anchor A (17, 77.5) and Equivalent B (20, 58.7)
pts_f5 <- data.frame(
  point = c("Anchor_A_17_77.5", "Equiv_B_20_58.7", "Anchor_A_17_72.4", "Equiv_B_21_60.8"),
  train_total_squat_sets = c(17, 20, 17, 21),
  train_avg_pct_1rm = c(77.5, 58.7, 72.4, 60.8),
  pre_1rm = median(df$pre_1rm),
  weight_kg = median(df$weight_kg),
  user_age_at_training = median(df$user_age_at_training),
  gender_male = 1,
  user_id_fac = df$user_id_fac[1]
)
p_pts <- predict(m_main_te, pts_f5, se.fit = TRUE, exclude = "s(user_id_fac)")
pts_f5$pred_fit <- p_pts$fit
pts_f5$pred_se <- p_pts$se.fit
pts_f5$ci_lower <- p_pts$fit - 1.96 * p_pts$se.fit
pts_f5$ci_upper <- p_pts$fit + 1.96 * p_pts$se.fit
write.csv(pts_f5, "/Users/speediance/Desktop/paper数据处理/pred_fig5_points_4377.csv", row.names = FALSE)

saveRDS(m_main_te, "/Users/speediance/Desktop/paper数据处理/model_fig5_te_4377.rds")
interaction_models <- list(vol_int = list(model = m_main_te, grid = grid_f5, raw = df))
saveRDS(interaction_models, "/Users/speediance/Desktop/one_rm_research/output/gamm_interaction_models.rds")
saveRDS(interaction_models, "/Users/speediance/Desktop/paper数据处理/gamm_interaction_models.rds")

# ==============================================================================
# 3. Subgroup Stratified GAMMs (Gender, Age, Training Status)
# ==============================================================================
cat(">>> [Step 4] Fitting Subgroup Stratified GAMM models (Gender, Age, Level)...\n")
int_seq <- seq(50, 85, length.out = 120)
vol_anchors <- c(8, 16, 24)

# --- 3.1 Gender Stratified ---
gender_models <- list()
gender_curves <- list()
gender_tables <- list()

for (g in c("Male", "Female")) {
  sub <- df %>% filter(gender_fac == g) %>% mutate(user_id_fac = factor(user_id))
  m_g <- bam(
    y_delta_1rm_pct ~ te(train_total_squat_sets, train_avg_pct_1rm, k = c(8, 8)) +
      s(pre_1rm, k = 10, bs = "tp") +
      s(weight_kg, k = 8, bs = "tp") +
      s(user_age_at_training, k = 6, bs = "tp") +
      s(user_id_fac, bs = "re"),
    data = sub, method = "fREML", discrete = TRUE
  )
  gender_models[[g]] <- m_g
  gender_tables[[g]] <- extract_model_table(m_g, g)
  
  for (v in vol_anchors) {
    nd <- data.frame(
      train_total_squat_sets = v,
      train_avg_pct_1rm = int_seq,
      pre_1rm = median(sub$pre_1rm),
      weight_kg = median(sub$weight_kg),
      user_age_at_training = median(sub$user_age_at_training),
      user_id_fac = sub$user_id_fac[1]
    )
    pr <- predict(m_g, nd, se.fit = TRUE, exclude = "s(user_id_fac)")
    nd$pred_fit <- pr$fit
    nd$pred_se <- pr$se.fit
    nd$ci_lower <- pr$fit - 1.96 * pr$se.fit
    nd$ci_upper <- pr$fit + 1.96 * pr$se.fit
    nd$total_sets <- v
    nd$gender <- g
    nd$n_windows <- nrow(sub)
    nd$n_users <- n_distinct(sub$user_id)
    gender_curves[[paste(g, v, sep = "_")]] <- nd
  }
}
write.csv(do.call(rbind, gender_curves), "/Users/speediance/Desktop/paper数据处理/curve_gender_4377.csv", row.names = FALSE)

# --- 3.2 Age Stratified ---
age_models <- list()
age_curves <- list()
age_tables <- list()

for (a in c("<=35", ">35")) {
  sub <- df %>% filter(age_group == a) %>% mutate(user_id_fac = factor(user_id))
  m_a <- bam(
    y_delta_1rm_pct ~ te(train_total_squat_sets, train_avg_pct_1rm, k = c(8, 8)) +
      s(pre_1rm, k = 10, bs = "tp") +
      s(weight_kg, k = 8, bs = "tp") +
      gender_male +
      s(user_id_fac, bs = "re"),
    data = sub, method = "fREML", discrete = TRUE
  )
  age_models[[a]] <- m_a
  age_tables[[a]] <- extract_model_table(m_a, a)
  
  for (v in vol_anchors) {
    nd <- data.frame(
      train_total_squat_sets = v,
      train_avg_pct_1rm = int_seq,
      pre_1rm = median(sub$pre_1rm),
      weight_kg = median(sub$weight_kg),
      gender_male = 1,
      user_id_fac = sub$user_id_fac[1]
    )
    pr <- predict(m_a, nd, se.fit = TRUE, exclude = "s(user_id_fac)")
    nd$pred_fit <- pr$fit
    nd$pred_se <- pr$se.fit
    nd$ci_lower <- pr$fit - 1.96 * pr$se.fit
    nd$ci_upper <- pr$fit + 1.96 * pr$se.fit
    nd$total_sets <- v
    nd$age_group <- a
    nd$n_windows <- nrow(sub)
    nd$n_users <- n_distinct(sub$user_id)
    age_curves[[paste(a, v, sep = "_")]] <- nd
  }
}
write.csv(do.call(rbind, age_curves), "/Users/speediance/Desktop/paper数据处理/curve_age_4377.csv", row.names = FALSE)

# --- 3.3 Training Status Stratified ---
level_models <- list()
level_curves <- list()
level_tables <- list()

for (l in c("Novice", "Intermediate", "Advanced")) {
  sub <- df %>% filter(strength_level_fac == l) %>% mutate(user_id_fac = factor(user_id))
  m_l <- bam(
    y_delta_1rm_pct ~ te(train_total_squat_sets, train_avg_pct_1rm, k = c(8, 8)) +
      s(pre_1rm, k = 10, bs = "tp") +
      s(weight_kg, k = 8, bs = "tp") +
      s(user_age_at_training, k = 6, bs = "tp") +
      gender_male +
      s(user_id_fac, bs = "re"),
    data = sub, method = "fREML", discrete = TRUE
  )
  level_models[[l]] <- m_l
  level_tables[[l]] <- extract_model_table(m_l, l)
  
  for (v in vol_anchors) {
    nd <- data.frame(
      train_total_squat_sets = v,
      train_avg_pct_1rm = int_seq,
      pre_1rm = median(sub$pre_1rm),
      weight_kg = median(sub$weight_kg),
      user_age_at_training = median(sub$user_age_at_training),
      gender_male = 1,
      user_id_fac = sub$user_id_fac[1]
    )
    pr <- predict(m_l, nd, se.fit = TRUE, exclude = "s(user_id_fac)")
    nd$pred_fit <- pr$fit
    nd$pred_se <- pr$se.fit
    nd$ci_lower <- pr$fit - 1.96 * pr$se.fit
    nd$ci_upper <- pr$fit + 1.96 * pr$se.fit
    nd$total_sets <- v
    nd$level_fac <- l
    nd$n_windows <- nrow(sub)
    nd$n_users <- n_distinct(sub$user_id)
    level_curves[[paste(l, v, sep = "_")]] <- nd
  }
}
write.csv(do.call(rbind, level_curves), "/Users/speediance/Desktop/paper数据处理/curve_level_4377.csv", row.names = FALSE)

# Combine all model tables
all_tables_df <- rbind(
  extract_model_table(m_sets, "Main_Sets_Fig4A"),
  extract_model_table(m_int, "Main_Intensity_Fig4B"),
  extract_model_table(m_main_te, "Main_Interaction_Table2"),
  do.call(rbind, gender_tables),
  do.call(rbind, age_tables),
  do.call(rbind, level_tables)
)
write.csv(all_tables_df, "/Users/speediance/Desktop/paper数据处理/unified_gamm_summary_tables_4377.csv", row.names = FALSE)

# Save master RDS
unified_rds <- list(
  fig4_sets = m_sets,
  fig4_int = m_int,
  table2_main_te = m_main_te,
  table4_gender = gender_models,
  table6_age = age_models,
  table8_level = level_models
)
saveRDS(unified_rds, "/Users/speediance/Desktop/paper数据处理/unified_gamm_models_4377.rds")
saveRDS(unified_rds, "/Users/speediance/Desktop/paper数据处理/model_gamm_4377.rds")

cat("\n>>> All 9 unified GAMM models fitted and exported successfully!\n")
print(all_tables_df %>% mutate(across(where(is.numeric), ~ round(., 3))))

