{
############################################
# eGFR slope (Slope-adjusted LME) only
# Windows: ≤1 year, ≤3 years, All period
# Figures: (1) 3 windows, (2) ≤1 year only
############################################

# ---- Packages ----
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(purrr)
library(nlme)
library(multcomp)
library(ggplot2)

# =========================================================
# Part 1) Build jin1_inclusion (jin_label, time0, years_from_time0, time0_egfr)
#   - Input: jin1_Eligibile (already in memory) OR read it from CSV
# =========================================================

# If needed:
setwd("X:/R")
jin1_Eligibile <- read_csv("jin1_Eligibile.csv", locale = locale(encoding = "SHIFT-JIS"))

jin1_inclusion <- jin1_Eligibile %>%
  filter(exclude == "include", jin_status %in% c("AKD", "nonAKD")) %>%
  mutate(
    jin_label = case_when(
      jin_status == "nonAKD" ~ "nonAKD",
      jin_status == "AKD" & `150_210recovery` == 1 ~ "Recovery",
      jin_status == "AKD" & `150_210recovery` == 2 ~ "Non-Recovery",
      jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 1 ~ "Recovery",
      jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` %in% c(0, 2) ~ "Non-Recovery",
      TRUE ~ NA_character_
    ),
    jin_label = factor(jin_label, levels = c("nonAKD", "Recovery", "Non-Recovery"))
  ) %>%
  filter(!is.na(jin_label))

# --- (A) max eGFR date within 90-210 days ---
egfr_max_date_210 <- jin1_inclusion %>%
  mutate(days_from_index = as.numeric(date - index_date)) %>%
  filter(days_from_index >= 90, days_from_index <= 210) %>%
  group_by(id) %>%
  filter(egfr == max(egfr, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(id, max_egfr_date_210 = date)

# --- (B) nearest date within 0-90 days (latest) ---
max_date_90 <- jin1_inclusion %>%
  mutate(days_from_index = as.numeric(date - index_date)) %>%
  filter(days_from_index >= 0, days_from_index <= 90) %>%
  group_by(id) %>%
  filter(days_from_index == max(days_from_index, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(id, nearest_date_90 = date)

# --- (C) define time0: max_egfr_date_210 > nearest_date_90 > index_date ---
jin1_inclusion <- jin1_inclusion %>%
  left_join(egfr_max_date_210, by = "id") %>%
  left_join(max_date_90, by = "id") %>%
  mutate(
    time0 = case_when(
      !is.na(max_egfr_date_210) ~ max_egfr_date_210,
      !is.na(nearest_date_90)   ~ nearest_date_90,
      TRUE                      ~ index_date
    ),
    years_from_time0 = as.numeric(date - time0) / 365.25
  )

# --- (D) time0_egfr: if multiple rows on time0, take max ---
time0_egfr_df <- jin1_inclusion %>%
  filter(date == time0) %>%
  group_by(id) %>%
  slice_max(egfr, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(id, time0_egfr = egfr)

jin1_inclusion <- jin1_inclusion %>%
  left_join(time0_egfr_df, by = "id")

# --- Save (optional) ---
# write_csv(jin1_inclusion, "jin1_inclusion.csv")

# =========================================================
# Part 2) Longitudinal dataset for modeling (years_from_time0 >= 0)
# =========================================================

akd_time_m <- jin1_inclusion %>%
  filter(years_from_time0 >= 0)

# CKD_status handling (keep if you need later; not used in slope-adjusted here)
akd_time_m <- akd_time_m %>%
  mutate(CKD_status = na_if(CKD_status, "nd"),
         CKD_status = factor(CKD_status, levels = c("nonCKD", "CKD")))

# Baseline covariates (one row per id; from Eligible include)
baseline_cov <- jin1_Eligibile %>%
  filter(exclude == "include") %>%
  group_by(id) %>%
  arrange(index_date) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(arb_acei_use = if_else(arb == 1 | acei == 1, 1L, 0L)) %>%
  dplyr::select(
    id, age, sex, arb_acei_use,
    dn1, dn3, dn4, dn5, dn6, dn7, dn8, dn9, dn10, dn12, dn13, dn14, dn15
  )

covars <- c("age","sex","arb_acei_use",
            "dn1","dn3","dn4","dn5","dn6","dn7","dn8","dn9","dn10","dn12","dn13","dn14","dn15")

longdat <- akd_time_m %>%
  dplyr::select(-any_of(covars)) %>%
  left_join(baseline_cov, by = "id") %>%
  mutate(
    age_c        = as.numeric(scale(age, center = TRUE, scale = FALSE)),
    time0_egfr_c  = as.numeric(scale(time0_egfr, center = TRUE, scale = FALSE)),
    jin_label    = factor(jin_label, levels = c("nonAKD", "Recovery", "Non-Recovery"))
  )

# =========================================================
# Part 3) Fit Slope-adjusted LME only (≤1y, ≤3y, All period)
# =========================================================

ctrl <- lmeControl(
  maxIter = 1e8, msMaxIter = 1e8,
  opt = "optim", optimMethod = "L-BFGS-B"
)

fit_slope_adj_all <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1 +
    age_c + sex + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15 +
    years_from_time0:(age_c + arb_acei_use +
                        dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15),
  random = list(id = pdSymm(~ 1 + years_from_time0)),
  data = longdat,
  na.action = na.omit,
  method = "REML",
  control = ctrl
)

fit_slope_adj_1y <- update(fit_slope_adj_all,
                           data = dplyr::filter(longdat, years_from_time0 <= 1))

fit_slope_adj_3y <- update(fit_slope_adj_all,
                           data = dplyr::filter(longdat, years_from_time0 <= 3))

fits <- list(
  "≤1 year"    = fit_slope_adj_1y,
  "≤3 years"   = fit_slope_adj_3y,
  "All period" = fit_slope_adj_all
)

# =========================================================
# Part 4) Extract group slopes (Estimate, 95%CI) and contrasts vs nonAKD (Diff, 95%CI, p)
# =========================================================

get_slopes_ci <- function(fit, window_label){
  cf <- names(fixef(fit))
  
  v_slope <- function(g){
    vec <- rep(0, length(cf)); names(vec) <- cf
    if("years_from_time0" %in% cf) vec["years_from_time0"] <- 1
    if(g == "Recovery" && "years_from_time0:jin_labelRecovery" %in% cf)
      vec["years_from_time0:jin_labelRecovery"] <- 1
    if(g == "Non-Recovery" && "years_from_time0:jin_labelNon-Recovery" %in% cf)
      vec["years_from_time0:jin_labelNon-Recovery"] <- 1
    vec
  }
  
  L <- rbind(
    nonAKD         = v_slope("nonAKD"),
    Recovery       = v_slope("Recovery"),
    `Non-Recovery` = v_slope("Non-Recovery")
  )
  
  ci <- suppressMessages(confint(glht(fit, linfct = L)))
  tibble(
    window   = window_label,
    group    = rownames(L),
    estimate = ci$confint[,"Estimate"],
    lower    = ci$confint[,"lwr"],
    upper    = ci$confint[,"upr"]
  )
}

get_contrasts_vs_nonakd <- function(fit, window_label){
  cf <- names(fixef(fit))
  
  v_slope <- function(g){
    vec <- rep(0, length(cf)); names(vec) <- cf
    if("years_from_time0" %in% cf) vec["years_from_time0"] <- 1
    if(g == "Recovery" && "years_from_time0:jin_labelRecovery" %in% cf)
      vec["years_from_time0:jin_labelRecovery"] <- 1
    if(g == "Non-Recovery" && "years_from_time0:jin_labelNon-Recovery" %in% cf)
      vec["years_from_time0:jin_labelNon-Recovery"] <- 1
    vec
  }
  
  b  <- v_slope("nonAKD")
  r  <- v_slope("Recovery")
  nr <- v_slope("Non-Recovery")
  
  K <- rbind(
    `Recovery − nonAKD`     = r  - b,
    `Non-Recovery − nonAKD` = nr - b
  )
  
  gl <- glht(fit, linfct = K)
  ci <- suppressMessages(confint(gl))
  sm <- suppressMessages(summary(gl))
  
  tibble(
    window     = window_label,
    contrast   = rownames(K),
    group      = case_when(
      contrast == "Recovery − nonAKD" ~ "Recovery",
      contrast == "Non-Recovery − nonAKD" ~ "Non-Recovery",
      TRUE ~ contrast
    ),
    diff_value = ci$confint[,"Estimate"],
    diff_lower = ci$confint[,"lwr"],
    diff_upper = ci$confint[,"upr"],
    diff_p     = sm$test$pvalues
  )
}

df_plot <- imap_dfr(fits, ~get_slopes_ci(.x, .y)) %>%
  mutate(
    group  = factor(group, levels = c("nonAKD","Recovery","Non-Recovery")),
    window = factor(window, levels = c("≤1 year","≤3 years","All period"))
  )

df_contrast <- imap_dfr(fits, ~get_contrasts_vs_nonakd(.x, .y))

# sample size: distinct id per window per group
sample_sizes <- bind_rows(
  longdat %>% filter(years_from_time0 <= 1) %>% mutate(window = "≤1 year"),
  longdat %>% filter(years_from_time0 <= 3) %>% mutate(window = "≤3 years"),
  longdat %>% mutate(window = "All period")
) %>%
  distinct(id, jin_label, window) %>%
  count(jin_label, window, name = "n") %>%
  transmute(
    group  = as.character(jin_label),
    window = factor(window, levels = c("≤1 year","≤3 years","All period")),
    n
  )

# add reference rows for nonAKD (Diff, p blank)
diff_ref <- expand.grid(
  window = factor(c("≤1 year","≤3 years","All period"), levels = c("≤1 year","≤3 years","All period")),
  group  = "nonAKD",
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  mutate(
    diff_value = NA_real_, diff_lower = NA_real_, diff_upper = NA_real_, diff_p = NA_real_
  )

differences <- bind_rows(
  diff_ref,
  df_contrast %>% dplyr::select(window, group, diff_value, diff_lower, diff_upper, diff_p) %>%
    mutate(window = factor(window, levels = c("≤1 year","≤3 years","All period")))
)

df_fig_enhanced <- df_plot %>%
  left_join(sample_sizes, by = c("group","window")) %>%
  left_join(differences,  by = c("group","window")) %>%
  mutate(
    group = factor(group, levels = c("nonAKD","Recovery","Non-Recovery"))
  )

# =========================================================
# Part 5) Publication-style plots (Renamed)
#   - Supplement Figure 1: ≤1y, ≤3y, All period (Slope-adjusted)
#   - Figure 4: ≤1y only (Slope-adjusted)
# =========================================================
make_pub_plot <- function(df_in, facet_by_window = TRUE, show_diff = TRUE, show_p = TRUE){
  
  y_min <- min(df_in$lower, na.rm = TRUE)
  y_max <- max(df_in$upper, na.rm = TRUE)
  
  p <- ggplot(df_in, aes(x = group, y = estimate, fill = group)) +
    geom_col(width = 0.7, alpha = 0.85) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, linewidth = 0.7) +
    geom_text(aes(y = pmax(upper, 0) + 0.5, label = paste0("n=", n)),
              size = 3.5, fontface = "bold", color = "grey20") +
    geom_text(aes(y = pmin(lower, 0) - 0.8,
                  label = sprintf("%.2f\n(%.2f, %.2f)", estimate, lower, upper)),
              size = 3, lineheight = 0.95, color = "grey10") +
    labs(
      x = NULL,
      y = expression(paste("Mean change in eGFR (mL/min/1.73 m"^2," per year)")),
      fill = "Group"
    ) +
    scale_fill_manual(
      values = c(nonAKD = "#95A5A6", Recovery = "#2ECC71", `Non-Recovery` = "#E74C3C"),
      labels = c(nonAKD = "No AKD", Recovery = "AKD with Recovery", `Non-Recovery` = "AKD without Recovery")
    ) +
    scale_x_discrete(expand = expansion(add = 0.6)) +   # ←★追加（左右余白）
    coord_cartesian(ylim = c(y_min - 5.8, y_max + 1.5), clip = "off") +
    theme_classic(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.6),
      panel.background = element_rect(fill = "white", color = NA),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 12, margin = margin(r = 10)),
      axis.text.y = element_text(size = 11),
      legend.position = "bottom",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      legend.key.size = unit(1.2, "lines"),
      legend.background = element_blank(),
      plot.margin = margin(t = 5, r = 20, b = 15, l = 25)
    )
  
  # ---- Diff text（棒から追い出す前提なので、必要なときだけ描画）----
  if (show_diff) {
    p <- p +
      geom_text(
        aes(label = ifelse(group == "nonAKD", "",
                           sprintf("Diff: %.2f\n(%.2f, %.2f)", diff_value, diff_lower, diff_upper))),
        y = y_min - 3.6,
        size = 3, lineheight = 0.95,
        fontface = "italic", color = "grey30"
      )
  }
  
  # ---- p-value text ----
  if (show_p) {
    p <- p +
      geom_text(
        aes(label = ifelse(group == "nonAKD", "",
                           dplyr::case_when(
                             is.na(diff_p) ~ "",
                             diff_p < 0.001 ~ "p<0.001",
                             diff_p < 0.01 ~ sprintf("p=%.3f", diff_p),
                             TRUE ~ sprintf("p=%.2f", diff_p)
                           ))),
        y = y_min - 4.7,
        size = 2.9,
        fontface = "bold", color = "grey20"
      )
  }
  
  if (facet_by_window) {
    p <- p + facet_grid(. ~ window)
  }
  
  p
}

# ---- Supplement Figure 1: 3 windows ----
p_supp_fig1 <- make_pub_plot(df_fig_enhanced, facet_by_window = TRUE)

# ---- Figure 4: ≤1 year only ----
p_fig4 <- make_pub_plot(
  df_fig_enhanced %>% filter(window == "≤1 year"),
  facet_by_window = FALSE
) +
  ggtitle("Slope-adjusted eGFR slope (≤1 year)")

# Preview
print(p_supp_fig1)
print(p_fig4)

# =========================================================
# Part 6) Save figures (PDF + TIFF) with new names
# =========================================================

outdir <- "X:/R"   # <- adjust if needed
w_mm <- 200; h_mm <- 140

# ---- Supplement Figure 1 ----
ggsave(file.path(outdir, "SupplementFigure1_eGFR_slope_adj_3windows.pdf"),
       plot = p_supp_fig1,
       width = w_mm, height = h_mm, units = "mm",
       device = cairo_pdf, dpi = 300)

ggsave(file.path(outdir, "SupplementFigure1_eGFR_slope_adj_3windows.tiff"),
       plot = p_supp_fig1,
       width = w_mm, height = h_mm, units = "mm",
       device = "tiff", dpi = 600, compression = "lzw")

# ---- Figure 4 ----
ggsave(file.path(outdir, "Figure4_eGFR_slope_adj_1y_only.pdf"),
       plot = p_fig4,
       width = w_mm, height = h_mm, units = "mm",
       device = cairo_pdf, dpi = 300)

ggsave(file.path(outdir, "Figure4_eGFR_slope_adj_1y_only.tiff"),
       plot = p_fig4,
       width = w_mm, height = h_mm, units = "mm",
       device = "tiff", dpi = 600, compression = "lzw")

} # 本解析
{
#========================
# 0) 前提
#   - fit_slope_adj_1y が存在する
#   - longdat に years_from_time0, jin_label, sex, age_c, time0_egfr_c, arb_acei_use, dn* がある
#========================

#------------------------
# 1) 予測用グリッドを作る（0〜1年）
#------------------------
t_grid <- seq(0, 1, by = 0.05)

# sex の参照レベル（longdat の factor になっている想定）
sex_ref <- if (is.factor(longdat$sex)) levels(longdat$sex)[1] else unique(longdat$sex)[1]

newdat <- expand.grid(
  years_from_time0 = t_grid,
  jin_label = factor(c("nonAKD","Recovery","Non-Recovery"),
                     levels = c("nonAKD","Recovery","Non-Recovery")),
  KEEP.OUT.ATTRS = FALSE
) %>%
  as_tibble() %>%
  mutate(
    age_c = 0,
    time0_egfr_c = 0,
    sex = if (is.factor(longdat$sex))
      factor(sex_ref, levels = levels(longdat$sex))
    else sex_ref,
    arb_acei_use = 0,
    dn1 = 0, dn3 = 0, dn4 = 0, dn5 = 0, dn6 = 0, dn7 = 0,
    dn8 = 0, dn9 = 0, dn10 = 0, dn12 = 0, dn13 = 0, dn14 = 0, dn15 = 0
  )



#------------------------
# 2) 固定効果のみの予測値と 95%CI を計算
#   - model.matrix でデザイン行列 X
#   - Var(beta)=vcov(fit) を使って SE を出す
#------------------------
X <- model.matrix(
  delete.response(terms(fit_slope_adj_1y)),
  newdat
)

beta <- fixef(fit_slope_adj_1y)
V <- vcov(fit_slope_adj_1y)

pred <- as.numeric(X %*% beta)
se   <- sqrt(diag(X %*% V %*% t(X)))

newdat_pred <- newdat %>%
  mutate(
    pred_egfr = pred,
    lwr = pred - 1.96 * se,
    upr = pred + 1.96 * se
  )

#------------------------
# 3-A) 調整済み eGFR 推定曲線（≤1年）
#------------------------
p_adj_egfr_1y <- ggplot(newdat_pred,
                        aes(x = years_from_time0, y = pred_egfr, color = jin_label, fill = jin_label)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, linewidth = 0) +
  geom_line(linewidth = 1.1) +
  labs(
    title = "Adjusted estimated eGFR trajectory.",
    x = "Time from time0 (years)",
    y = expression(paste("Adjusted eGFR (mL/min/1.73 m"^2, ")")),
    color = "Group", fill = "Group"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right")

print(p_adj_egfr_1y)

# 0.25年ごとの点だけ抽出（誤差回避のため丸める）
bar_dat <- newdat_pred %>%
  mutate(t_round = round(years_from_time0, 2)) %>%
  filter(t_round %in% round(seq(0, 1, by = 0.25), 2)) %>%
  dplyr::select(-t_round)

p_adj_egfr_1y_bar <- ggplot(newdat_pred,
                            aes(x = years_from_time0, y = pred_egfr,
                                color = jin_label, group = jin_label)) +
  # 平滑な推定曲線（0.05年刻み）
  geom_line(linewidth = 1.1) +
  
  # 0.25年ごとの点
  geom_point(data = bar_dat, size = 2.4) +
  
  # 0.25年ごとの95%CIバー
  geom_errorbar(
    data = bar_dat,
    aes(ymin = lwr, ymax = upr),
    width = 0.03, linewidth = 0.7
  ) +
  
  labs(
    title = "Adjusted estimated eGFR trajectory",
    x = "Time from time0 (years)",
    y = expression(paste("Adjusted eGFR (mL/min/1.73 m"^2, ")")),
    color = "Group"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right")

print(p_adj_egfr_1y_bar)

}#1年以内の折れ線グラフ
{
# =========================================================
# Figure 3 (final): two-panel (A/B) + two-line x-axis label
#  - A: adjusted eGFR trajectory (≤1 year)  -> p_adj_egfr_1y_bar
#  - B: slope-adjusted bar plot (≤1 year)  -> p_fig4 (from make_pub_plot)
#  - Figure-level title with line break
#  - X-axis label in 2 lines; show only on bottom panel to avoid duplication
#  - A/B tags added
# =========================================================
library(dplyr)
library(ggplot2)
library(patchwork)

# ---- A) top panel (trajectory) ----
p_line <- p_adj_egfr_1y_bar +
  labs(
    tag = "A",
    x = "Time from time0 (year)",
    y = "Adjusted eGFR\n(mL/min/1.73 m²)"
  ) +
  theme(
    plot.margin = margin(b = 5),
    legend.position = "right",
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0, 0.98)
  )

# ---- B) bottom panel (bar plot; relaxed spacing only) ----
df_1y <- df_fig_enhanced %>%
  dplyr::filter(window == "≤1 year") %>%
  mutate(group = factor(group, levels = c("nonAKD","Recovery","Non-Recovery")))

# 下側の基準（Diff / p / 注釈）を“段”として分離
make_pub_plot_F4 <- function(df_in, facet_by_window = TRUE,
                             show_diff = TRUE, show_p = TRUE){
  
  y_min <- min(df_in$lower, na.rm = TRUE)
  y_max <- max(df_in$upper, na.rm = TRUE)
  # ---- positions for bottom annotations (increase spacing) ----
  y_diff_pos <- y_min - 5.2
  y_pval_pos <- y_min - 7.4
  y_note_pos <- y_min - 8.9
  y_n_pos <- y_max + 0.8   # ← バーの一番上より常に上
  
  p <- ggplot(df_in, aes(x = group, y = estimate, fill = group)) +
    geom_col(width = 0.7, alpha = 0.85) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper),
      width = 0.25, linewidth = 0.7
    ) +
    
    # ---- n ----
  geom_text(
    aes(y = y_n_pos,
        label = paste0("n=", n)),
    size = 3.5,
    fontface = "bold",
    color = "grey20"
  ) +
    
    # ---- Estimate + 95%CI（★ここだけ下へ： -0.8 → -1.8）----
  geom_text(
    aes(
      y = pmin(lower, 0) - 1.8,
      label = sprintf("%.2f\n(%.2f, %.2f)",
                      estimate, lower, upper)
    ),
    size = 3,
    lineheight = 0.95,
    color = "grey10"
  ) +
    
    labs(
      tag = "B",
      x = NULL,
      y = expression(
        paste("Mean change in eGFR \n(mL/min/1.73 m² per year)")
      ),
      fill = "Group"
    ) +
    
    scale_fill_manual(
      values = c(
        nonAKD = "#95A5A6",
        Recovery = "#2ECC71",
        `Non-Recovery` = "#E74C3C"
      ),
      labels = c(
        nonAKD = "No AKD",
        Recovery = "AKD with Recovery",
        `Non-Recovery` = "AKD without Recovery"
      )
    ) +
    
    coord_cartesian(
      ylim = c(y_min - 9.8, y_max + 1.5),
      clip = "off"
    ) +
    
    theme_classic(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.6),
      panel.background = element_rect(fill = "white", color = NA),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 12, margin = margin(r = 10)),
      axis.text.y = element_text(size = 11),
      legend.position = "bottom",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      legend.key.size = unit(1.2, "lines"),
      legend.background = element_blank(),
      plot.margin = margin(t = 5, r = 20, b = 15, l = 25)
    )
  
  # ---- Diff ----
  if (show_diff) {
    p <- p +
      geom_text(
        aes(
          label = ifelse(
            group == "nonAKD", "Reference",
            sprintf("Diff: %.2f\n(%.2f, %.2f)",
                    diff_value, diff_lower, diff_upper)
          )
        ),
        y = y_diff_pos,
        size = 3.1,
        lineheight = 0.95,
        fontface = "italic",
        color = "grey30"
      )
  }
  
  # ---- p-value ----
  if (show_p) {
    p <- p +
      geom_text(
        aes(
          label = ifelse(
            group == "nonAKD", "",
            dplyr::case_when(
              is.na(diff_p) ~ "",
              diff_p < 0.001 ~ "p<0.001",
              diff_p < 0.01  ~ sprintf("p=%.3f", diff_p),
              TRUE           ~ sprintf("p=%.2f", diff_p)
            )
          )
        ),
        y = y_pval_pos,
        size = 3.0,
        fontface = "bold",
        color = "grey20"
      )
  }
  
  if (facet_by_window) {
    p <- p + facet_grid(. ~ window)
  }
  
  p
}
# ---- Combine into Figure 3 ----
p_fig4_relaxed <- make_pub_plot_F4(
  df_fig_enhanced %>% dplyr::filter(window == "≤1 year"),
  facet_by_window = FALSE
) +
  ggtitle("Adjusted differences in annual eGFR change within 1 year") +
  scale_x_discrete(expand = expansion(add = 0.8)) +
  theme(
    plot.margin = margin(t = 5, r = 40, b = 20, l = 40)
  )

fig4 <- p_line / p_fig4_relaxed +
  plot_layout(heights = c(2.5, 5.5)) +
  plot_annotation(
    title = "Figure 3. Adjusted eGFR trajectory and slope differences\nwithin 1 year after time0",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0,
                                margin = margin(b = 12))
    )
  )

print(fig4)

# ---- Save ----
ggsave(
  "Figure3_eGFR_trajectory_and_slope_1y.pdf",
  plot  = fig4,
  width = 230, height = 180, units = "mm",
  device = cairo_pdf
)

ggsave(
  "Figure3_eGFR_trajectory_and_slope_1y.tiff",
  plot  = fig4,
  width = 230, height = 180, units = "mm",
  dpi = 600, compression = "lzw"
)

}#1年以内の折れ線グラフと下向き棒グラフの結合
{
  ############################################################
  # Combine SupplementalFigure1A (trajectory) + 1B (bar)
  # into one figure (A/B panels) and save
  ############################################################
  
  library(patchwork)
  library(ggplot2)
  
  # ---- Panel A: trajectory ----
  p_supp1A_tag <- p_supp1A +
    labs(
      tag   = "A",
      title = "Adjusted eGFR trajectory"
    ) +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0),
      plot.tag = element_text(face = "bold", size = 14),
      plot.tag.position = c(0, 0.98)
    )
  
  # ---- Panel B: bar (slope differences) ----
  p_supp1B_tag <- p_supp1B +
    labs(
      tag   = "B",
      title = "Adjusted differences in annual eGFR change"
    ) +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0),
      plot.tag = element_text(face = "bold", size = 14),
      plot.tag.position = c(0, 0.98)
    )
  
  # ---- Combine (vertical) ----
  supp_fig1 <- p_supp1A_tag / p_supp1B_tag +
    plot_layout(heights = c(2.2, 3.2)) +
    plot_annotation(
      title = "Supplemental Figure 1. Adjusted eGFR trajectory and slope differences\nwithin 3 years and all period",
      theme = theme(
        plot.title = element_text(
          size = 14, face = "bold", hjust = 0,
          margin = margin(b = 12)
        )
      )
    )
  
  print(supp_fig1)
  
  # ---- Save ----
  outdir <- "X:/R"
  
  ggsave(
    filename = file.path(outdir, "SupplementalFigure1_main_trajectory_and_bar_3y_all.pdf"),
    plot = supp_fig1,
    width = 230, height = 260, units = "mm",
    device = cairo_pdf
  )
  
  ggsave(
    filename = file.path(outdir, "SupplementalFigure1_main_trajectory_and_bar_3y_all.tiff"),
    plot = supp_fig1,
    width = 230, height = 260, units = "mm",
    dpi = 600, compression = "lzw"
  )
  


}#3年以内と全期間

{

#感度分析(recoveryの定義変更)#####
  ############################################################
  # Sensitivity analysis (recovery definition change)
  # Slope-adjusted models (≤1y / ≤3y / All period)
  ############################################################
  
  library(readr)
  library(dplyr)
  library(tidyr)
  library(nlme)
  library(multcomp)
  library(purrr)
  
  setwd("X:/R")
  
  # -----------------------------
  # 0) Load
  # -----------------------------
  jin1_Eligibile <- read_csv("jin1_Eligibile.csv", locale = locale(encoding = "SHIFT-JIS"))
  jin1_inclusion <- read_csv("jin1_inclusion.csv", locale = locale(encoding = "SHIFT-JIS"))
  
  # -----------------------------
  # 1) Sensitivity AKD label (exclude No-data later)
  # -----------------------------
  jin1_inclusion_sens <- jin1_inclusion %>%
    filter(exclude == "include", jin_status %in% c("AKD","nonAKD")) %>%
    mutate(
      jin_label_sens = case_when(
        jin_status == "nonAKD" ~ "nonAKD",
        jin_status == "AKD" & `150_210recovery` == 1 ~ "Recovery",
        jin_status == "AKD" & `150_210recovery` == 2 ~ "Non-Recovery",
        jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 1 ~ "Recovery",
        jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 2 ~ "Non-Recovery",
        jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 0 ~ "No-data",
        TRUE ~ NA_character_
      ),
      jin_label_sens = factor(
        jin_label_sens,
        levels = c("nonAKD","Recovery","Non-Recovery","No-data")
      )
    )
  
  # -----------------------------
  # 2) Time variable
  # -----------------------------
  akd_time_sens <- jin1_inclusion_sens %>%
    mutate(
      years_from_time0 = as.numeric(date - time0) / 365.25
    ) %>%
    filter(years_from_time0 >= 0)
  
  # -----------------------------
  # 3) Baseline covariates
  # -----------------------------
  baseline_cov <- jin1_Eligibile %>%
    filter(exclude == "include") %>%
    group_by(id) %>%
    arrange(index_date) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(arb_acei_use = if_else(arb == 1 | acei == 1, 1L, 0L)) %>%
    dplyr::select(
      id, age, sex, arb_acei_use,
      dn1, dn3, dn4, dn5, dn6, dn7, dn8,
      dn9, dn10, dn12, dn13, dn14, dn15
    )
  
  # -----------------------------
  # 4) Long data for modeling
  # -----------------------------
  longdat_sens <- akd_time_sens %>%
    filter(jin_label_sens != "No-data") %>%
    left_join(baseline_cov, by = "id") %>%
    mutate(
      age_c        = as.numeric(scale(age.y, center = TRUE, scale = FALSE)),
      time0_egfr_c  = as.numeric(scale(time0_egfr, center = TRUE, scale = FALSE)),
      sex          = sex.y,              # baseline の sex を採用
      dn1  = dn1.y, dn3  = dn3.y, dn4  = dn4.y, dn5  = dn5.y, dn6  = dn6.y,
      dn7  = dn7.y, dn8  = dn8.y, dn9  = dn9.y, dn10 = dn10.y,
      dn12 = dn12.y, dn13 = dn13.y, dn14 = dn14.y, dn15 = dn15.y,
      jin_label_sens = factor(jin_label_sens, levels = c("nonAKD","Recovery","Non-Recovery"))
    )
  longdat_sens <- longdat_sens %>%
    dplyr::select(-dplyr::ends_with(".x"), -dplyr::ends_with(".y"))
  
  # -----------------------------
  # 5) Data by window
  # -----------------------------
  dat_1y  <- longdat_sens %>% filter(years_from_time0 <= 1)
  dat_3y  <- longdat_sens %>% filter(years_from_time0 <= 3)
  dat_all <- longdat_sens

  
  ctrl <- lmeControl(
    maxIter = 1e8,
    msMaxIter = 1e8,
    opt = "optim",
    optimMethod = "L-BFGS-B"
  )
  form_slope_adj <- egfr ~ years_from_time0 * jin_label_sens + time0_egfr_c - 1 +
    age_c + sex + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 +
    dn9 + dn10 + dn12 + dn13 + dn14 + dn15 +
    years_from_time0:(age_c + arb_acei_use +
                        dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 +
                        dn9 + dn10 + dn12 + dn13 + dn14 + dn15)
  # ---- ≤1 year ----
  fit_sens_1y <- lme(
    fixed   = form_slope_adj,
    random  = list(id = pdSymm(~ 1 + years_from_time0)),
    data    = dat_1y,
    na.action = na.omit,
    method  = "REML",
    control = ctrl
  )
  
  # ---- ≤3 years ----
  fit_sens_3y <- lme(
    fixed   = form_slope_adj,
    random  = list(id = pdSymm(~ 1 + years_from_time0)),
    data    = dat_3y,
    na.action = na.omit,
    method  = "REML",
    control = ctrl
  )
  
  # ---- All period ----
  fit_sens_all <- lme(
    fixed   = form_slope_adj,
    random  = list(id = pdSymm(~ 1 + years_from_time0)),
    data    = dat_all,
    na.action = na.omit,
    method  = "REML",
    control = ctrl
  )
  
  
  library(dplyr)
  library(purrr)
  library(multcomp)
  
  # ---------------------------------------------------------
  # 1) 群ごとの slope (Estimate, 95%CI) を返す関数
  # ---------------------------------------------------------
  get_slopes_ci_sens <- function(fit, window_label){
    
    cf <- names(fixef(fit))
    
    v <- function(g){
      vec <- rep(0, length(cf)); names(vec) <- cf
      
      # base slope
      vec["years_from_time0"] <- 1
      
      # group x time interactions
      if(g == "Recovery" && "years_from_time0:jin_label_sensRecovery" %in% cf)
        vec["years_from_time0:jin_label_sensRecovery"] <- 1
      
      if(g == "Non-Recovery" && "years_from_time0:jin_label_sensNon-Recovery" %in% cf)
        vec["years_from_time0:jin_label_sensNon-Recovery"] <- 1
      
      vec
    }
    
    L <- rbind(
      nonAKD         = v("nonAKD"),
      Recovery       = v("Recovery"),
      `Non-Recovery` = v("Non-Recovery")
    )
    
    ci <- suppressMessages(confint(glht(fit, linfct = L)))
    
    tibble(
      window   = window_label,
      group    = rownames(L),
      estimate = ci$confint[,"Estimate"],
      lower    = ci$confint[,"lwr"],
      upper    = ci$confint[,"upr"]
    )
  }
  
  # ---------------------------------------------------------
  # 2) nonAKDとの差 (Diff, 95%CI, p) を返す関数
  # ---------------------------------------------------------
  get_diff_vs_nonakd_sens <- function(fit, window_label){
    
    cf <- names(fixef(fit))
    
    v <- function(g){
      vec <- rep(0, length(cf)); names(vec) <- cf
      vec["years_from_time0"] <- 1
      if(g == "Recovery" && "years_from_time0:jin_label_sensRecovery" %in% cf)
        vec["years_from_time0:jin_label_sensRecovery"] <- 1
      if(g == "Non-Recovery" && "years_from_time0:jin_label_sensNon-Recovery" %in% cf)
        vec["years_from_time0:jin_label_sensNon-Recovery"] <- 1
      vec
    }
    
    b  <- v("nonAKD")
    r  <- v("Recovery")
    nr <- v("Non-Recovery")
    
    K <- rbind(
      `Recovery − nonAKD`     = r  - b,
      `Non-Recovery − nonAKD` = nr - b
    )
    
    g  <- glht(fit, linfct = K)
    ci <- suppressMessages(confint(g))
    sm <- suppressMessages(summary(g))
    
    tibble(
      window     = window_label,
      group      = c("Recovery","Non-Recovery"),
      diff_value = ci$confint[,"Estimate"],
      diff_lower = ci$confint[,"lwr"],
      diff_upper = ci$confint[,"upr"],
      diff_p     = sm$test$pvalues
    )
  }
  
  # ---------------------------------------------------------
  # 3) windowごとの n（distinct id）を作る
  # ---------------------------------------------------------
  sample_sizes_sens <- bind_rows(
    dat_1y  %>% mutate(window = "≤1 year"),
    dat_3y  %>% mutate(window = "≤3 years"),
    dat_all %>% mutate(window = "All period")
  ) %>%
    distinct(id, jin_label_sens, window) %>%
    count(jin_label_sens, window, name = "n") %>%
    transmute(
      group  = as.character(jin_label_sens),
      window = window,
      n
    )
  
  # ---------------------------------------------------------
  # 4) slope と diff を 3 window 分まとめて抽出
  # ---------------------------------------------------------
  fits_sens <- list(
    "≤1 year"    = fit_sens_1y,
    "≤3 years"   = fit_sens_3y,
    "All period" = fit_sens_all
  )
  
  df_plot_sens <- imap_dfr(fits_sens, ~get_slopes_ci_sens(.x, .y))
  
  df_diff_sens <- imap_dfr(fits_sens, ~get_diff_vs_nonakd_sens(.x, .y)) %>%
    # nonAKD 用の reference 行を追加（Diff/pはNA）
    bind_rows(
      tibble(
        window = rep(names(fits_sens), each = 1),
        group  = "nonAKD",
        diff_value = NA_real_,
        diff_lower = NA_real_,
        diff_upper = NA_real_,
        diff_p     = NA_real_
      )
    )
  
  # ---------------------------------------------------------
  # 5) 図用の最終データ df_fig_sens（これが欲しかったもの）
  # ---------------------------------------------------------
  df_fig_sens <- df_plot_sens %>%
    left_join(sample_sizes_sens, by = c("group","window")) %>%
    left_join(df_diff_sens,      by = c("group","window")) %>%
    mutate(
      group  = factor(group, levels = c("nonAKD","Recovery","Non-Recovery")),
      window = factor(window, levels = c("≤1 year","≤3 years","All period"))
    )
  
  # 確認
  df_fig_sens %>% count(window, group) 
  
}#感度分析
{  
  ############################################################
  # Figure 5A / 5B (Sensitivity analysis, ≤1 year)
  #   5A: Adjusted eGFR trajectory (≤1y)
  #   5B: Slope differences (≤1y, downward bars)
  ############################################################
  
  library(dplyr)
  library(ggplot2)
  library(nlme)
  library(patchwork)
  library(tibble)  
  
  # ---- prediction helper (fixed effects only) ----
  predict_fixed_lme_1y <- function(fit, times, sex_ref = 0){
    
    newdat <- expand.grid(
      years_from_time0 = times,
      jin_label_sens = factor(
        c("nonAKD","Recovery","Non-Recovery"),
        levels = c("nonAKD","Recovery","Non-Recovery")
      ),
      KEEP.OUT.ATTRS = FALSE
    ) %>%
      tibble::as_tibble() %>%       # ← tibble を確実に使えるようにした
      dplyr::mutate(
        age_c = 0,
        time0_egfr_c = 0,
        sex = sex_ref,
        arb_acei_use = 0,
        dn1  = 0, dn3  = 0, dn4  = 0, dn5  = 0, dn6  = 0, dn7  = 0,
        dn8  = 0, dn9  = 0, dn10 = 0, dn12 = 0, dn13 = 0, dn14 = 0, dn15 = 0
      )
    
    X <- model.matrix(delete.response(terms(fit)), newdat)
    beta <- fixef(fit)
    V <- vcov(fit)
    
    pred <- as.numeric(X %*% beta)
    se   <- sqrt(diag(X %*% V %*% t(X)))
    
    newdat %>%
      dplyr::mutate(
        pred_egfr = pred,
        lwr = pred - 1.96 * se,
        upr = pred + 1.96 * se
      )
  }
  
  # ---- time grids ----
  t_line <- seq(0, 1, by = 0.05)
  t_bar  <- seq(0, 1, by = 0.25)
  
  # ★ここで pred_1y を作る（エラーがあればここで止まる）
  pred_1y <- predict_fixed_lme_1y(fit_sens_1y, t_line)
  
  pred_1y_bar <- pred_1y %>%
    mutate(t_round = round(years_from_time0, 2)) %>%
    filter(t_round %in% round(t_bar, 2)) %>%
    dplyr::select(-t_round)
  
  # ---- Figure 5A ----
  p_fig6A <- ggplot(pred_1y,
                    aes(x = years_from_time0, y = pred_egfr,
                        color = jin_label_sens, group = jin_label_sens)) +
    geom_line(linewidth = 1.1) +
    geom_point(data = pred_1y_bar, size = 2.3) +
    geom_errorbar(
      data = pred_1y_bar,
      aes(ymin = lwr, ymax = upr),
      width = 0.04, linewidth = 0.7
    ) +
    labs(
      tag = "A",
      x = "Time from time0 (year)",
      y = "Adjusted eGFR\n(mL/min/1.73 m²)",
      color = "Group"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.tag = element_text(face = "bold", size = 14),
      plot.tag.position = c(0, 0.98)
    )
  
  # ---- Figure 5B（df_fig_sens が既に作られている前提） ----
  df_1y <- df_fig_sens %>%
    filter(window == "≤1 year") %>%
    mutate(group = factor(group, levels = c("nonAKD","Recovery","Non-Recovery")))
  
  y_base <- min(df_1y$lower, na.rm = TRUE)
  y_n    <- max(df_1y$upper, na.rm = TRUE) + 0.8
  y_est  <- pmin(df_1y$lower, 0) - 1.8
  y_diff <- y_base - 5.2
  y_pval <- y_base - 7.0
  
  p_fig6B <- ggplot(df_1y, aes(x = group, y = estimate, fill = group)) +
    geom_col(width = 0.7, alpha = 0.85) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  width = 0.25, linewidth = 0.7) +
    geom_text(aes(y = y_n, label = paste0("n=", n)),
              size = 3.5, fontface = "bold", color = "grey20") +
    geom_text(aes(y = y_est,
                  label = sprintf("%.2f\n(%.2f, %.2f)", estimate, lower, upper)),
              size = 3, lineheight = 0.95, color = "grey10") +
    geom_text(aes(y = y_diff,
                  label = ifelse(group == "nonAKD", "Reference",
                                 sprintf("Diff: %.2f\n(%.2f, %.2f)",
                                         diff_value, diff_lower, diff_upper))),
              size = 3.1, lineheight = 0.95,
              fontface = "italic", color = "grey30") +
    geom_text(aes(y = y_pval,
                  label = ifelse(group == "nonAKD", "",
                                 case_when(
                                   diff_p < 0.001 ~ "p<0.001",
                                   diff_p < 0.01  ~ sprintf("p=%.3f", diff_p),
                                   TRUE           ~ sprintf("p=%.2f", diff_p)
                                 ))),
              size = 3.0, fontface = "bold", color = "grey20") +
    labs(
      tag = "B",
      x = NULL,
      y = "Mean change in eGFR\n(mL/min/1.73 m² per year)",
      fill = "Group"
    ) +
    scale_fill_manual(
      values = c(nonAKD = "#95A5A6", Recovery = "#2ECC71", `Non-Recovery` = "#E74C3C"),
      labels = c(nonAKD = "No AKD",
                 Recovery = "AKD with Recovery",
                 `Non-Recovery` = "AKD without Recovery")
    ) +
    scale_x_discrete(expand = expansion(add = 0.8)) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 13) +
    theme(
      panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.6),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "bottom",
      plot.tag = element_text(face = "bold", size = 14),
      plot.tag.position = c(0, 0.98),
      plot.margin = margin(t = 5, r = 25, b = 20, l = 25)
    )
  
  fig6 <- p_fig6A / p_fig6B +
    plot_layout(heights = c(2.5, 5.5)) +
    plot_annotation(
      title = "Figure 5. Sensitivity analysis: adjusted eGFR trajectory and slope differences\nwithin 1 year after time0",
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0,
                                  margin = margin(b = 12))
      )
    )
  
  print(fig6)
  
  # ---- Save ----
  ggsave("Figure5_sensitivity_trajectory_and_slope_1y.pdf",
         plot = fig6, width = 230, height = 180, units = "mm",
         device = cairo_pdf)
  
  ggsave("Figure5_sensitivity_trajectory_and_slope_1y.tiff",
         plot = fig6, width = 230, height = 180, units = "mm",
         dpi = 600, compression = "lzw")
  
}#1年以内
{
  # =========================================================
  # Supplement Figure 2B: Bar plots (≤3 years + All period)
  #   - uses df_fig_sens (derived from fit_sens_3y & fit_sens_all)
  # =========================================================
  
  df_supp2A <- df_fig_sens %>%
    filter(window %in% c("≤3 years","All period")) %>%
    mutate(
      window = factor(window, levels = c("≤3 years","All period")),
      group  = factor(group,  levels = c("nonAKD","Recovery","Non-Recovery"))
    )
  
  # Per-window spacing variables to prevent overlaps
  df_supp2A <- df_supp2A %>%
    group_by(window) %>%
    mutate(
      y_min_w = min(lower, na.rm = TRUE),
      y_max_w = max(upper, na.rm = TRUE),
      y_n_pos    = y_max_w + 0.8,
      y_est_pos  = pmin(lower, 0) - 1.8,
      y_diff_pos = y_min_w - 5.2,
      y_pval_pos = y_min_w - 7.4,
      y_note_pos = y_min_w - 8.9
    ) %>%
    ungroup()
  
  p_supp2A <- ggplot(df_supp2A, aes(x = group, y = estimate, fill = group)) +
    geom_col(width = 0.7, alpha = 0.85) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, linewidth = 0.7) +
    geom_text(aes(y = y_n_pos, label = paste0("n=", n)),
              size = 3.5, fontface = "bold", color = "grey20") +
    geom_text(aes(y = y_est_pos,
                  label = sprintf("%.2f\n(%.2f, %.2f)", estimate, lower, upper)),
              size = 3, lineheight = 0.95, color = "grey10") +
    geom_text(aes(y = y_diff_pos,
                  label = ifelse(group == "nonAKD", "Reference",
                                 sprintf("Diff: %.2f\n(%.2f, %.2f)", diff_value, diff_lower, diff_upper))),
              size = 3.1, lineheight = 0.95, fontface = "italic", color = "grey30") +
    geom_text(aes(y = y_pval_pos,
                  label = ifelse(group == "nonAKD", "",
                                 dplyr::case_when(
                                   is.na(diff_p) ~ "",
                                   diff_p < 0.001 ~ "p<0.001",
                                   diff_p < 0.01  ~ sprintf("p=%.3f", diff_p),
                                   TRUE           ~ sprintf("p=%.2f", diff_p)
                                 ))),
              size = 3.0, fontface = "bold", color = "grey20") +
    facet_grid(. ~ window) +
    labs(
      x = NULL,
      y = "Mean change in eGFR\n(mL/min/1.73 m² per year)",
      fill = "Group",
      title = "Supplemental Figure 2B. Sensitivity analysis: slope differences \n(≤3 years and All period)"
    ) +
    scale_fill_manual(
      values = c(nonAKD = "#95A5A6", Recovery = "#2ECC71", `Non-Recovery` = "#E74C3C"),
      labels = c(nonAKD = "No AKD", Recovery = "AKD with Recovery", `Non-Recovery` = "AKD without Recovery")
    ) +
    scale_x_discrete(expand = expansion(add = 0.8)) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 13) +
    theme(
      panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.6),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "bottom",
      plot.margin = margin(t = 10, r = 20, b = 20, l = 20)
    )
  
  print(p_supp2A)
  
  ggsave("SupplementalFigure2B_sensitivity_bar_3y_all.pdf", plot = p_supp2A,
         width = 230, height = 140, units = "mm", device = cairo_pdf)
  ggsave("SupplementalFigure2B_sensitivity_bar_3y_all.tiff", plot = p_supp2A,
         width = 230, height = 140, units = "mm", dpi = 600, compression = "lzw")
  
  # =========================================================
  # Supplement Figure 2A: Trajectories (≤3 years + All period)
  #   - ≤3 years: uses fit_sens_3y
  #   - All period: uses fit_sens_all
  # =========================================================
  
  # ≤3 years trajectory: 0-3 years
  times_3y_line <- seq(0, 3, by = 0.10)
  times_3y_bar  <- seq(0, 3, by = 0.25)
  
  pred_3y <- predict_fixed_lme_sens(fit_sens_3y, times_3y_line) %>%
    mutate(window = "≤3 years")
  
  pred_3y_bar <- pred_3y %>%
    mutate(t_round = round(years_from_time0, 2)) %>%
    filter(t_round %in% round(times_3y_bar, 2)) %>%
    dplyr::select(-t_round)
  
  # All period trajectory: choose a reasonable x-range for display
  max_all <- max(longdat_sens$years_from_time0, na.rm = TRUE)
  max_all <- min(max_all, 10)  # <- remove this cap if you want the full range
  
  times_all_line <- seq(0, max_all, by = 0.25)
  times_all_bar  <- seq(0, max_all, by = 0.5)
  
  pred_all <- predict_fixed_lme_sens(fit_sens_all, times_all_line) %>%
    mutate(window = "All period")
  
  pred_all_bar <- pred_all %>%
    mutate(t_round = round(years_from_time0, 2)) %>%
    filter(t_round %in% round(times_all_bar, 2)) %>%
    dplyr::select(-t_round)
  
  traj_df <- bind_rows(pred_3y, pred_all) %>%
    mutate(window = factor(window, levels = c("≤3 years","All period")))
  
  traj_bar_df <- bind_rows(pred_3y_bar, pred_all_bar) %>%
    mutate(window = factor(window, levels = c("≤3 years","All period")))
  
  p_supp2B <- ggplot(traj_df,
                     aes(x = years_from_time0, y = pred_egfr,
                         color = jin_label_sens, group = jin_label_sens)) +
    geom_line(linewidth = 1.1) +
    geom_point(data = traj_bar_df, size = 2.2) +
    geom_errorbar(
      data = traj_bar_df,
      aes(ymin = lwr, ymax = upr),
      width = 0.05, linewidth = 0.6
    ) +
    facet_grid(. ~ window, scales = "free_x") +
    labs(
      x = "Time from time0 (year)",
      y = "Adjusted eGFR\n(mL/min/1.73 m²)",
      color = "Group",
      title = "Supplemental Figure 2A. Sensitivity analysis: adjusted eGFR trajectories \n(≤3 years and All period)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "right",
      plot.margin = margin(t = 10, r = 20, b = 10, l = 20)
    )
  
  print(p_supp2B)
  
  ggsave("Supplemental_Figure2A_sensitivity_trajectory_3y_all.pdf", plot = p_supp2B,
         width = 230, height = 120, units = "mm", device = cairo_pdf)
  ggsave("Supplemental_Figure2A_sensitivity_trajectory_3y_all.tiff", plot = p_supp2B,
         width = 230, height = 120, units = "mm", dpi = 600, compression = "lzw")
  
  ############################################################
  # Supplement Figure 2A / 2B (Sensitivity analysis)
  #   2A: Bar plots (≤3 years, All period) using df_fig_sens
  #   2B: Trajectories (≤3 years uses fit_sens_3y,
  #                     All period uses fit_sens_all)
  ############################################################
    
    library(dplyr)
    library(ggplot2)
    library(nlme)
    # patchwork を確実に使える状態にする
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)
  # =========================================================
  # Helper: fixed-effect prediction with Wald 95% CI
  #   (covariates fixed at reference/mean)
  # =========================================================
  predict_fixed_lme_sens <- function(fit, times,
                                     group_levels = c("nonAKD","Recovery","Non-Recovery"),
                                     sex_ref = 0){
    newdat <- expand.grid(
      years_from_time0 = times,
      jin_label_sens   = factor(group_levels, levels = group_levels),
      KEEP.OUT.ATTRS = FALSE
    ) %>%
      tibble::as_tibble() %>%
      mutate(
        age_c = 0,
        time0_egfr_c = 0,
        sex = sex_ref,
        arb_acei_use = 0,
        dn1 = 0, dn3 = 0, dn4 = 0, dn5 = 0, dn6 = 0, dn7 = 0,
        dn8 = 0, dn9 = 0, dn10 = 0, dn12 = 0, dn13 = 0, dn14 = 0, dn15 = 0
      )
    
    X <- model.matrix(delete.response(terms(fit)), newdat)
    beta <- fixef(fit)
    V <- vcov(fit)
    
    pred <- as.numeric(X %*% beta)
    se   <- sqrt(diag(X %*% V %*% t(X)))
    
    newdat %>%
      mutate(
        pred_egfr = pred,
        lwr = pred - 1.96 * se,
        upr = pred + 1.96 * se
      )
  }
  
  # ---- Panel A: trajectory（上） ----
  p_supp2B_tag <- p_supp2B +
    labs(
      tag   = "A",
      title = "Adjusted eGFR trajectory"
    ) +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0),
      plot.tag = element_text(face = "bold", size = 14),
      plot.tag.position = c(0, 0.98)
    )
  
  # ---- Panel B: eGFR change / slope differences（下） ----
  p_supp2A_tag <- p_supp2A +
    labs(
      tag   = "B",
      title = "Adjusted differences in annual eGFR change"
    ) +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0),
      plot.tag = element_text(face = "bold", size = 14),
      plot.tag.position = c(0, 0.98)
    )
  
  library(patchwork)
  
  # ★ 保存先ディレクトリを定義（これが抜けていた）
  outdir <- "X:/R"
  
  # ---- Combine (vertical): A (trajectory) on top, B (bar) bottom ----
  supp_fig2 <- p_supp2B_tag / p_supp2A_tag +
    plot_layout(heights = c(2.2, 3.2)) +
    plot_annotation(
      title = "Supplemental Figure 2. Sensitivity analysis: adjusted eGFR trajectory and slope differences\nwithin 3 years and all period"
    ) &
    theme(
      plot.title = element_text(
        size = 14, face = "bold", hjust = 0,
        margin = margin(b = 12)
      )
    )
  
  # ---- Save ----
  ggsave(
    filename = file.path(outdir, "SupplementalFigure2_sensitivity_bar_and_trajectory_3y_all.pdf"),
    plot = supp_fig2,
    width = 230, height = 260, units = "mm",
    device = cairo_pdf
  )
  
}#3年以内と全期間



