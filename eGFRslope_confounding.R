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
# Part 5) Publication-style plots (Renamed)  ※修正版
#   - "Recovery" 表記 → "recovery"（表示ラベルのみ）
#   - 棒の中（/下）に出していた annual eGFR change（estimate, 95%CI）を削除
#   - Difference vs nonAKD のみ残す
# =========================================================
make_pub_plot <- function(df_in, facet_by_window = TRUE, show_diff = TRUE, show_p = TRUE){
  
  y_min <- min(df_in$lower, na.rm = TRUE)
  y_max <- max(df_in$upper, na.rm = TRUE)
  
  p <- ggplot(df_in, aes(x = group, y = estimate, fill = group)) +
    geom_col(width = 0.7, alpha = 0.85) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, linewidth = 0.7) +
    geom_text(aes(y = pmax(upper, 0) + 0.5, label = paste0("n=", n)),
              size = 3.5, fontface = "bold", color = "grey20") +
    # ★削除：annual eGFR change の Estimate + 95%CI 表示（固定効果の数値）
    # geom_text(aes(y = pmin(lower, 0) - 0.8,
    #               label = sprintf("%.2f\n(%.2f, %.2f)", estimate, lower, upper)),
    #           size = 3, lineheight = 0.95, color = "grey10") +
    labs(
      x = NULL,
      y = expression(paste("Mean change in eGFR (mL/min/1.73 m"^2," per year)")),
      fill = "Group"
    ) +
    scale_fill_manual(
      values = c(nonAKD = "#95A5A6", Recovery = "#2ECC71", `Non-Recovery` = "#E74C3C"),
      # ★ここだけ変更：表示上 “Recovery”→“recovery”
      labels = c(
        nonAKD = "Non AKD",
        Recovery = "AKD with recovery",
        `Non-Recovery` = "AKD without recovery"
      )
    ) +
    scale_x_discrete(expand = expansion(add = 0.6)) +
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
  
  # ---- Difference vs nonAKD（これだけ残す）----
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
  
  # ---- p-value text（必要なら残す）----
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

# ---- Rebuild plots ----
p_supp_fig1 <- make_pub_plot(df_fig_enhanced, facet_by_window = TRUE)

p_fig4 <- make_pub_plot(
  df_fig_enhanced %>% filter(window == "≤1 year"),
  facet_by_window = FALSE
) +
  ggtitle("Slope-adjusted eGFR slope (≤1 year)")

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

} # 本解析

{
library(dplyr)
library(ggplot2)
# years_from_time0 が無い場合は先に作る
akd_time_m <- akd_time_m %>%
  mutate(
    years_from_time0 = as.numeric(date - index_date) / 365.25
  ) %>%
  filter(years_from_time0 >= 0)

# 各ID内の測定間隔（年単位）の中央値を計算
median_interval <- akd_time_m %>%
  arrange(id, years_from_time0) %>%
  group_by(id) %>%
  summarise(
    d = diff(years_from_time0),
    .groups = "drop"
  ) %>%
  pull(d) %>%
  median(na.rm = TRUE)

# window幅を定義
window_width <- median_interval * 0.75

median_interval
window_width

#========================
# 0) 前提（追加）
#   - akd_time_m に id, years_from_time0, jin_label, egfr がある（longdat相当）
#   - target_timepoints がある（例：seq(0,1,by=0.05)など）
#   - window_width がある（観測間隔から決めた幅）
#========================

#------------------------
# 1) 規定時点（0〜1年）
#    ※あなたの t_grid をそのまま流用
#------------------------
t_grid <- seq(0, 1, by = 0.05)
target_timepoints <- t_grid

half_w <- window_width / 2

#------------------------
# 2) window幅ありの「実測」eGFRを規定時点にアライン
#   - 各 id × target_time で窓内の最も近い1点を採用
#------------------------
library(data.table)

dt <- as.data.table(akd_time_m)
dt <- dt[!is.na(egfr)]
dt[, years_from_time0 := as.numeric(date - index_date) / 365.25]
dt <- dt[years_from_time0 >= 0 & years_from_time0 <= 1]

# 規定時点（0.25年刻みならこちら）
target_timepoints <- seq(0, 1, by = 0.25)
# 0.05刻みなら： seq(0, 1, by = 0.05)

half_w <- window_width / 2

# グリッド：id × target_time
ids  <- unique(dt$id)
grid <- CJ(id = ids, target_time = target_timepoints)

setkey(dt, id, years_from_time0)

# ★ joinの「その場」で target_time を列として返す（ここが重要）
# ※ jin_label が無いなら jin_status に置換してください
window_data_obs <- dt[
  grid,
  on = .(id, years_from_time0 = target_time),
  roll = "nearest",
  nomatch = 0L,
  .(id,
    target_time = i.target_time,
    years_from_time0,
    egfr,
    jin_label)   # <- 無ければ jin_status に
]

# 窓内だけ残す
window_data_obs <- window_data_obs[abs(years_from_time0 - target_time) <= half_w]

#------------------------
# 2-b) 集団平均と 95%CI（実測ベース）
#   ※CIは mean ± 1.96*SE（SE=SD/sqrt(n)）
#------------------------
newdat_pred <- window_data_obs %>%
  group_by(jin_label, target_time) %>%
  summarise(
    pred_egfr = mean(egfr, na.rm = TRUE),
    sd       = sd(egfr, na.rm = TRUE),
    n        = sum(!is.na(egfr)),
    se       = sd / sqrt(n),
    lwr      = pred_egfr - 1.96 * se,
    upr      = pred_egfr + 1.96 * se,
    .groups  = "drop"
  )

#------------------------
# 3-A) 実測（window-aligned）eGFR 推定曲線（≤1年）
#   ※あなたの図の体裁を維持
#------------------------
#------------------------
# 3-A) Observed eGFR trajectory (window-aligned) ≤1y
#------------------------

# ---- labels（キーを values_group と完全一致させる）----
lab_group <- c(
  nonAKD         = "non-AKD",
  Recovery       = "AKD with recovery",
  `Non-Recovery` = "AKD without recovery"
)

# ==== Figure2と同じ配色に統一 ====
values_group <- c(
  nonAKD         = "#95A5A6",  # grey
  Recovery       = "#2ECC71",  # green
  `Non-Recovery` = "#E74C3C"   # red
)
# ★重要：jin_label を factor 化（レベル固定）— keys を names(lab_group) に揃える
newdat_pred <- newdat_pred %>%
  mutate(jin_label = factor(as.character(jin_label), levels = names(lab_group)))

# ---- A) ribbon + line ----
p_obs_egfr_1y <- ggplot(newdat_pred,
                        aes(x = target_time, y = pred_egfr,
                            color = jin_label, fill = jin_label)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, linewidth = 0) +
  geom_line(linewidth = 1.1) +
  labs(
    title = "Observed eGFR trajectory (window-aligned)",
    x = "Time from time0 (years)",
    y = expression(paste("Observed eGFR (mL/min/1.73 m"^2, ")")),
    color = "Group", fill = "Group"
  ) +
  scale_color_manual(values = values_group, breaks = names(lab_group), labels = lab_group, drop = FALSE) +
  scale_fill_manual(values  = values_group, breaks = names(lab_group), labels = lab_group, drop = FALSE) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right")

# ---- 0.25年点だけ抽出 ----
bar_dat <- newdat_pred %>%
  mutate(t_round = round(target_time, 2)) %>%
  filter(t_round %in% round(seq(0, 1, by = 0.25), 2)) %>%
  dplyr::select(-t_round)

# ---- B) line + points + errorbar ----
p_obs_egfr_1y_bar <- ggplot(newdat_pred,
                            aes(x = target_time, y = pred_egfr,
                                color = jin_label, group = jin_label)) +
  geom_line(linewidth = 1.1) +
  geom_point(data = bar_dat, size = 2.4) +
  geom_errorbar(data = bar_dat, aes(ymin = lwr, ymax = upr),
                width = 0.03, linewidth = 0.7) +
  labs(
    title = "Observed eGFR trajectory (window-aligned)",
    x = "Time from time0 (years)",
    y = expression(paste("Observed eGFR (mL/min/1.73 m"^2, ")")),
    color = "Group"
  ) +
  scale_color_manual(values = values_group, breaks = names(lab_group), labels = lab_group, drop = FALSE) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right")

# ---- x軸 0.25刻み ----
x_breaks_02 <- seq(0, 1, by = 0.25)

p_obs_egfr_1y <- p_obs_egfr_1y +
  scale_x_continuous(breaks = x_breaks_02, limits = c(0, 1))

p_obs_egfr_1y_bar <- p_obs_egfr_1y_bar +
  scale_x_continuous(breaks = x_breaks_02, limits = c(0, 1))

print(p_obs_egfr_1y)
print(p_obs_egfr_1y_bar)

}#1年以内の折れ線グラフ(実測値)
{
############################################################
# Observed ΔeGFR trajectory (window-aligned) ≤1y
#  - baseline = each id's eGFR at target_time==0 (window-aligned)
#  - delta_egfr = egfr - baseline_egfr
############################################################

library(dplyr)
library(ggplot2)
library(data.table)

# ---------------------------------------------------------
# 0) 前提：あなたのコードで作った window_data_obs がある前提
#    window_data_obs: id, target_time, years_from_time0, egfr, jin_label
# ---------------------------------------------------------
# window_data_obs <- ...（あなたの既存コードで作成済み）

wd <- as.data.table(window_data_obs)

# 念のため型
wd[, target_time := as.numeric(target_time)]
wd[, egfr := as.numeric(egfr)]

# ---------------------------------------------------------
# 1) baseline eGFR（各idの target_time==0 の eGFR）
#    ※ 0が取れないIDは除外（deltaが定義できないため）
# ---------------------------------------------------------
base_dt <- wd[target_time == 0, .(baseline_egfr = egfr[1]), by = id]

# baseline が存在する観測だけ残して結合
wd2 <- merge(wd, base_dt, by = "id", all.x = FALSE, all.y = FALSE)

# ---------------------------------------------------------
# 2) 変化量（ΔeGFR）
# ---------------------------------------------------------
wd2[, delta_egfr := egfr - baseline_egfr]

# ---------------------------------------------------------
# 3) 集団平均と95%CI（ΔeGFRベース）
# ---------------------------------------------------------
newdat_delta <- as.data.frame(wd2) %>%
  group_by(jin_label, target_time) %>%
  summarise(
    mean_delta = mean(delta_egfr, na.rm = TRUE),
    sd         = sd(delta_egfr, na.rm = TRUE),
    n          = sum(!is.na(delta_egfr)),
    se         = sd / sqrt(n),
    lwr        = mean_delta - 1.96 * se,
    upr        = mean_delta + 1.96 * se,
    .groups    = "drop"
  )

# ---------------------------------------------------------
# 4) 図（あなたの体裁を踏襲）
# ---------------------------------------------------------
lab_group <- c(
  nonAKD         = "non-AKD",
  Recovery       = "AKD with recovery",
  `Non-Recovery` = "AKD without recovery"
)

values_group <- c(
  nonAKD         = "#95A5A6",
  Recovery       = "#2ECC71",
  `Non-Recovery` = "#E74C3C"
)

newdat_delta <- newdat_delta %>%
  mutate(jin_label = factor(as.character(jin_label), levels = names(lab_group)))

# ---- ribbon + line（ΔeGFR）----
p_obs_delta_1y <- ggplot(newdat_delta,
                         aes(x = target_time, y = mean_delta,
                             color = jin_label, fill = jin_label)) +
  geom_hline(yintercept = 0, linewidth = 0.5, alpha = 0.5) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, linewidth = 0) +
  geom_line(linewidth = 1.1) +
  labs(
    title = "Observed ΔeGFR trajectory (window-aligned; baseline = time0)",
    x = "Time from time0 (years)",
    y = expression(paste(Delta,"eGFR (mL/min/1.73 m"^2,")")),
    color = "Group", fill = "Group"
  ) +
  scale_color_manual(values = values_group, breaks = names(lab_group), labels = lab_group, drop = FALSE) +
  scale_fill_manual(values  = values_group, breaks = names(lab_group), labels = lab_group, drop = FALSE) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right")

# ---- 0.25年点だけ抽出（point+errorbar版）----
bar_dat_delta <- newdat_delta %>%
  mutate(t_round = round(target_time, 2)) %>%
  filter(t_round %in% round(seq(0, 1, by = 0.25), 2)) %>%
  dplyr::select(-t_round)

p_obs_delta_1y_bar <- ggplot(newdat_delta,
                             aes(x = target_time, y = mean_delta,
                                 color = jin_label, group = jin_label)) +
  geom_hline(yintercept = 0, linewidth = 0.5, alpha = 0.5) +
  geom_line(linewidth = 1.1) +
  geom_point(data = bar_dat_delta, size = 2.4) +
  geom_errorbar(data = bar_dat_delta, aes(ymin = lwr, ymax = upr),
                width = 0.03, linewidth = 0.7) +
  labs(
    title = "Observed ΔeGFR trajectory (window-aligned; baseline = time0)",
    x = "Time from time0 (years)",
    y = expression(paste(Delta,"eGFR (mL/min/1.73 m"^2,")")),
    color = "Group"
  ) +
  scale_color_manual(values = values_group, breaks = names(lab_group), labels = lab_group, drop = FALSE) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right")

# ---- x軸 0.25刻み ----
x_breaks_02 <- seq(0, 1, by = 0.25)
p_obs_delta_1y     <- p_obs_delta_1y     + scale_x_continuous(breaks = x_breaks_02, limits = c(0, 1))
p_obs_delta_1y_bar <- p_obs_delta_1y_bar + scale_x_continuous(breaks = x_breaks_02, limits = c(0, 1))

print(p_obs_delta_1y)
print(p_obs_delta_1y_bar)
} #1年以内変化量
{
############################################################
# Figure 3 (2 versions)
#  (1) Observed eGFR trajectory  + slope differences (bar)
#  (2) Observed ΔeGFR trajectory + slope differences (bar)
# Save to: X:/R/eGFRslope/
############################################################

library(dplyr)
library(ggplot2)
library(patchwork)
library(Cairo)
library(ragg)

# -------------------------
# Paths
# -------------------------
setwd("X:/R")
outdir <- file.path(getwd(), "eGFRslope")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# 0) 前提
#   - p_obs_egfr_1y_bar   : 実測eGFRの折れ線（あなたの既存）
#   - p_obs_delta_1y_bar  : 変化量ΔeGFRの折れ線（先ほど作成したもの）
#   - p_fig4              : 下向き棒グラフ（make_pub_plot()で作成済み）
# =========================================================

# =========================================================
# 1) Bottom panel (共通)：下向き棒グラフ（Figure3用整形）
# =========================================================
p_fig4_for_fig3 <- p_fig4 +
  labs(
    tag   = "B",
    title = "Adjusted differences in annual eGFR change within 1 year",
    y     = expression(paste("Mean change in eGFR \n(mL/min/1.73 m"^2," per year)"))
  ) +
  scale_x_discrete(expand = expansion(add = 0.8)) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0),
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0, 0.98),
    plot.margin = margin(t = 5, r = 40, b = 20, l = 40)
  )

# =========================================================
# 2) Figure 3-1：Observed eGFR（実測） + bar
# =========================================================
p_line_egfr <- p_obs_egfr_1y_bar +
  labs(
    tag = "A",
    x = "Time from time0 (year)",
    y = "Observed eGFR\n(mL/min/1.73 m²)"
  ) +
  theme(
    plot.margin = margin(b = 5),
    legend.position = "right",
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0, 0.98)
  )

fig3_egfr <- p_line_egfr / p_fig4_for_fig3 +
  plot_layout(heights = c(2.5, 5.5)) +
  plot_annotation(
    title = "Figure 3. Observed eGFR trajectory and slope differences\nwithin 1 year after time0",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0,
                                margin = margin(b = 12))
    )
  )

print(fig3_egfr)

ggsave(
  filename = file.path(outdir, "Figure3A_observed_eGFR_trajectory_and_slope_1y.pdf"),
  plot  = fig3_egfr,
  width = 230, height = 180, units = "mm",
  device = cairo_pdf
)

ggsave(
  filename = file.path(outdir, "Figure3A_observed_eGFR_trajectory_and_slope_1y.tiff"),
  plot  = fig3_egfr,
  width = 230, height = 180, units = "mm",
  device = ragg::agg_tiff,
  dpi = 600, compression = "lzw"
)

# =========================================================
# 3) Figure 3-2：Observed ΔeGFR（変化量） + bar
# =========================================================
p_line_delta <- p_obs_delta_1y_bar +
  labs(
    tag = "A",
    x = "Time from time0 (year)",
    y = expression(paste(Delta,"eGFR\n(mL/min/1.73 m"^2,")"))
  ) +
  theme(
    plot.margin = margin(b = 5),
    legend.position = "right",
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0, 0.98)
  )

fig3_delta <- p_line_delta / p_fig4_for_fig3 +
  plot_layout(heights = c(2.5, 5.5)) +
  plot_annotation(
    title = "Figure 3. Observed ΔeGFR trajectory and slope differences\nwithin 1 year after time0",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0,
                                margin = margin(b = 12))
    )
  )

print(fig3_delta)

ggsave(
  filename = file.path(outdir, "Figure3B_observed_delta_eGFR_trajectory_and_slope_1y.pdf"),
  plot  = fig3_delta,
  width = 230, height = 180, units = "mm",
  device = cairo_pdf
)

ggsave(
  filename = file.path(outdir, "Figure3B_observed_delta_eGFR_trajectory_and_slope_1y.tiff"),
  plot  = fig3_delta,
  width = 230, height = 180, units = "mm",
  device = ragg::agg_tiff,
  dpi = 600, compression = "lzw"
)
}#1年以内の折れ線グラフと下向き棒グラフの結合
{
############################################################
# Supplemental Figure 1 (Main analysis) — 2 versions
#  Ver-A) Observed eGFR trajectory (≤3y + All) + Downward bar
#  Ver-B) Observed ΔeGFR trajectory (≤3y + All) + Downward bar
# Save to: X:/R/eGFRslope/
############################################################

library(dplyr)
library(ggplot2)
library(patchwork)
library(data.table)
library(tidyr)
library(Cairo)
library(ragg)

# -----------------------------
# 0) color / legend（指定どおり）
# -----------------------------
col_group <- c(
  nonAKD         = "#95A5A6",
  Recovery       = "#2ECC71",
  `Non-Recovery` = "#E74C3C"
)
lab_group <- c(
  nonAKD         = "Non-AKD",
  Recovery       = "AKD with recovery",
  `Non-Recovery` = "AKD without recovery"
)

# -----------------------------
# Save folder
# -----------------------------
setwd("X:/R")
outdir <- file.path(getwd(), "eGFRslope")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ==========================================================
# 1) Window-aligned observed trajectory function (MAIN)
#    - metric = "egfr"  : 実測eGFR
#    - metric = "delta" : ΔeGFR（baseline = target_time==0 の各id値）
# ==========================================================
make_window_aligned_traj_main <- function(dat_in, end_time,
                                          by_time = 0.25,
                                          width_multiplier = 0.75,
                                          window_label = "≤3 years",
                                          metric = c("egfr","delta")) {
  
  metric <- match.arg(metric)
  
  # --- years_from_time0 を確保（無ければ index_date 起点で作る） ---
  dat_in <- dat_in %>%
    mutate(
      years_from_time0 = if ("years_from_time0" %in% names(dat_in)) years_from_time0
      else as.numeric(date - index_date) / 365.25
    )
  
  # --- range + clean ---
  win_dat <- dat_in %>%
    filter(!is.na(egfr), !is.na(years_from_time0)) %>%
    filter(years_from_time0 >= 0, years_from_time0 <= end_time) %>%
    filter(!is.na(jin_label)) %>%
    mutate(jin_label = as.character(jin_label)) %>%
    filter(jin_label %in% c("nonAKD","Recovery","Non-Recovery"))
  
  # --- window width from within-id median interval ---
  median_interval <- win_dat %>%
    arrange(id, years_from_time0) %>%
    group_by(id) %>%
    summarise(d = diff(years_from_time0), .groups = "drop") %>%
    pull(d) %>%
    median(na.rm = TRUE)
  
  window_width <- median_interval * width_multiplier
  half_w <- window_width / 2
  
  # --- target grid ---
  target_timepoints <- seq(0, end_time, by = by_time)
  
  # --- rolling nearest join ---
  dt <- as.data.table(win_dat)
  ids  <- unique(dt$id)
  grid <- CJ(id = ids, target_time = target_timepoints)
  
  setkey(dt, id, years_from_time0)
  
  window_obs <- dt[
    grid,
    on = .(id, years_from_time0 = target_time),
    roll = "nearest",
    nomatch = 0L,
    .(id,
      target_time = i.target_time,
      years_from_time0,
      egfr,
      jin_label)
  ]
  
  window_obs <- window_obs[abs(years_from_time0 - target_time) <= half_w]
  
  # --- ΔeGFRにする場合：baseline（target_time==0 の各id eGFR）を引く ---
  if (metric == "delta") {
    base_dt <- window_obs[target_time == 0, .(baseline_egfr = egfr[1]), by = id]
    window_obs <- merge(window_obs, base_dt, by = "id", all.x = FALSE, all.y = FALSE)
    window_obs[, value := egfr - baseline_egfr]
  } else {
    window_obs[, value := egfr]
  }
  
  # --- mean ± 95%CI ---
  traj <- as.data.frame(window_obs) %>%
    as_tibble() %>%
    group_by(jin_label, target_time) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      sd         = sd(value, na.rm = TRUE),
      n          = sum(!is.na(value)),
      se         = sd / sqrt(n),
      lwr        = mean_value - 1.96 * se,
      upr        = mean_value + 1.96 * se,
      .groups    = "drop"
    ) %>%
    mutate(
      window = window_label,
      jin_label = factor(jin_label, levels = c("nonAKD","Recovery","Non-Recovery"))
    )
  
  traj
}

# ==========================================================
# 2) trajectory data: ≤3y + All period（egfr / delta の両方）
# ==========================================================
akd_time_m2 <- akd_time_m %>%
  mutate(years_from_time0 = if ("years_from_time0" %in% names(akd_time_m)) years_from_time0
         else as.numeric(date - index_date) / 365.25) %>%
  filter(years_from_time0 >= 0)

end_all_main <- max(akd_time_m2$years_from_time0, na.rm = TRUE)

# ---- Observed eGFR ----
traj_3y_main_egfr <- make_window_aligned_traj_main(
  dat_in = akd_time_m2, end_time = 3, by_time = 0.25,
  window_label = "≤3 years", metric = "egfr"
)
traj_all_main_egfr <- make_window_aligned_traj_main(
  dat_in = akd_time_m2, end_time = end_all_main, by_time = 0.5,
  window_label = "All period", metric = "egfr"
)
traj_supp1A_egfr <- bind_rows(traj_3y_main_egfr, traj_all_main_egfr) %>%
  mutate(window = factor(window, levels = c("≤3 years","All period")))

# ---- Observed ΔeGFR ----
traj_3y_main_delta <- make_window_aligned_traj_main(
  dat_in = akd_time_m2, end_time = 3, by_time = 0.25,
  window_label = "≤3 years", metric = "delta"
)
traj_all_main_delta <- make_window_aligned_traj_main(
  dat_in = akd_time_m2, end_time = end_all_main, by_time = 0.5,
  window_label = "All period", metric = "delta"
)
traj_supp1A_delta <- bind_rows(traj_3y_main_delta, traj_all_main_delta) %>%
  mutate(window = factor(window, levels = c("≤3 years","All period")))

# ==========================================================
# 3) Plot A（Ver-A: eGFR / Ver-B: ΔeGFR）
# ==========================================================
make_Apanel <- function(traj_df, y_lab) {
  ggplot(traj_df,
         aes(x = target_time, y = mean_value,
             color = jin_label, group = jin_label)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.1) +
    geom_errorbar(aes(ymin = lwr, ymax = upr),
                  width = 0.04, linewidth = 0.7) +
    facet_grid(. ~ window, scales = "free_x") +
    labs(
      x = "Time from time0 (years)",
      y = y_lab,
      color = "Group"
    ) +
    scale_color_manual(values = col_group, labels = lab_group, drop = FALSE) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "right")
}

p_supp1A_egfr <- make_Apanel(traj_supp1A_egfr,
                             "Observed eGFR\n(mL/min/1.73 m²)")

p_supp1A_delta <- make_Apanel(traj_supp1A_delta,
                              expression(paste("Observed ", Delta, "eGFR\n(mL/min/1.73 m"^2,")")))

# ==========================================================
# 4) Plot B: downward bars（共通）— Diff only（vs nonAKD）
# ==========================================================
df_bar_supp1B <- df_fig_enhanced %>%
  filter(window %in% c("≤3 years","All period")) %>%
  mutate(
    group  = factor(as.character(group), levels = c("nonAKD","Recovery","Non-Recovery")),
    window = factor(as.character(window), levels = c("≤3 years","All period"))
  )

y_base <- min(df_bar_supp1B$lower, na.rm = TRUE)
y_n    <- max(df_bar_supp1B$upper, na.rm = TRUE) + 0.8
y_diff <- y_base - 5.2
y_pval <- y_base - 7.0

p_supp1B <- ggplot(df_bar_supp1B, aes(x = group, y = estimate, fill = group)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.25, linewidth = 0.7) +
  geom_text(aes(y = y_n, label = paste0("n=", n)),
            size = 3.5, fontface = "bold", color = "grey20") +
  geom_text(aes(y = y_diff,
                label = ifelse(group == "nonAKD", "Reference",
                               sprintf("Diff: %.2f\n(%.2f, %.2f)",
                                       diff_value, diff_lower, diff_upper))),
            size = 3.1, lineheight = 0.95,
            fontface = "italic", color = "grey30") +
  geom_text(aes(y = y_pval,
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
    fill = "Group"
  ) +
  scale_fill_manual(values = col_group, labels = lab_group, drop = FALSE) +
  scale_x_discrete(expand = expansion(add = 0.8)) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.6),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(t = 5, r = 25, b = 20, l = 25)
  )

# ==========================================================
# 5) タグ付け＋結合（2種類）
# ==========================================================
# --- Ver-A: eGFR ---
p_supp1A_egfr_tag <- p_supp1A_egfr +
  labs(tag = "A", title = "Observed eGFR trajectory (window-aligned)") +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0),
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0, 0.98)
  )

p_supp1B_tag <- p_supp1B +
  labs(tag = "B", title = "Adjusted differences in annual eGFR change") +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0),
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0, 0.98)
  )

supp_fig1_main_egfr <- p_supp1A_egfr_tag / p_supp1B_tag +
  plot_layout(heights = c(2.3, 3.3)) +
  plot_annotation(
    title = "Supplemental Figure 1. Observed eGFR trajectory and slope differences\nwithin 3 years and all period after time0",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0,
                                margin = margin(b = 12))
    )
  )

# --- Ver-B: ΔeGFR ---
p_supp1A_delta_tag <- p_supp1A_delta +
  labs(tag = "A", title = "Observed ΔeGFR trajectory (window-aligned; baseline = time0)") +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0),
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0, 0.98)
  )

supp_fig1_main_delta <- p_supp1A_delta_tag / p_supp1B_tag +
  plot_layout(heights = c(2.3, 3.3)) +
  plot_annotation(
    title = "Supplemental Figure 1. Observed ΔeGFR trajectory and slope differences\nwithin 3 years and all period after time0",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0,
                                margin = margin(b = 12))
    )
  )

print(supp_fig1_main_egfr)
print(supp_fig1_main_delta)

# ==========================================================
# 6) Save (PDF + TIFF) — eGFRslope folder
# ==========================================================
ggsave(
  filename = file.path(outdir, "SupplementalFigure1A_main_observed_eGFR_trajectory_and_bar_3y_all.pdf"),
  plot  = supp_fig1_main_egfr,
  width = 230, height = 180, units = "mm",
  device = cairo_pdf
)

ggsave(
  filename = file.path(outdir, "SupplementalFigure1A_main_observed_eGFR_trajectory_and_bar_3y_all.tiff"),
  plot  = supp_fig1_main_egfr,
  width = 230, height = 180, units = "mm",
  device = ragg::agg_tiff,
  dpi = 600, compression = "lzw"
)

ggsave(
  filename = file.path(outdir, "SupplementalFigure1B_main_observed_delta_eGFR_trajectory_and_bar_3y_all.pdf"),
  plot  = supp_fig1_main_delta,
  width = 230, height = 180, units = "mm",
  device = cairo_pdf
)

ggsave(
  filename = file.path(outdir, "SupplementalFigure1B_main_observed_delta_eGFR_trajectory_and_bar_3y_all.tiff"),
  plot  = supp_fig1_main_delta,
  width = 230, height = 180, units = "mm",
  device = ragg::agg_tiff,
  dpi = 600, compression = "lzw"
)

# ---- end ----
}#3年以内と全期間
{
############################################################
# Sensitivity Figure 5 (≤1y): 2 versions
#  Ver-A) 5A: Observed eGFR trajectory (window-aligned) + 5B: slope bar
#  Ver-B) 5A: Observed ΔeGFR trajectory (window-aligned; baseline=time0) + 5B: slope bar
#
# Save to: X:/R/sensitivity_analysis/
#   Figure5A_sens_observed_eGFR_trajectory_and_slope_1y.(pdf/tiff)
#   Figure5B_sens_observed_delta_eGFR_trajectory_and_slope_1y.(pdf/tiff)
############################################################

library(readr)
library(dplyr)
library(tidyr)
library(nlme)
library(multcomp)
library(purrr)
library(ggplot2)
library(patchwork)
library(data.table)
library(Cairo)
library(ragg)

setwd("X:/R")

# -----------------------------
# Save folder
# -----------------------------
outdir <- file.path(getwd(), "sensitivity_analysis")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ==========================================================
# 0) Load
# ==========================================================
jin1_Eligibile <- read_csv("jin1_Eligibile.csv", locale = locale(encoding = "SHIFT-JIS"))
jin1_inclusion <- read_csv("jin1_inclusion.csv", locale = locale(encoding = "SHIFT-JIS"))

# ==========================================================
# 1) Sensitivity label (nonAKD / Recovery / Non-Recovery / No-data)
# ==========================================================
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
    jin_label_sens = factor(jin_label_sens,
                            levels = c("nonAKD","Recovery","Non-Recovery","No-data"))
  )

# ==========================================================
# 2) Time variable
# ==========================================================
akd_time_sens <- jin1_inclusion_sens %>%
  mutate(years_from_time0 = as.numeric(date - time0) / 365.25) %>%
  filter(years_from_time0 >= 0)

# ==========================================================
# 3) Baseline covariates
# ==========================================================
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

# ==========================================================
# 4) Long data for modeling (exclude No-data)
# ==========================================================
longdat_sens <- akd_time_sens %>%
  filter(jin_label_sens != "No-data") %>%
  left_join(baseline_cov, by = "id") %>%
  mutate(
    age_c        = as.numeric(scale(age.y, center = TRUE, scale = FALSE)),
    time0_egfr_c  = as.numeric(scale(time0_egfr, center = TRUE, scale = FALSE)),
    sex          = sex.y,
    dn1  = dn1.y, dn3  = dn3.y, dn4  = dn4.y, dn5  = dn5.y, dn6  = dn6.y,
    dn7  = dn7.y, dn8  = dn8.y, dn9  = dn9.y, dn10 = dn10.y,
    dn12 = dn12.y, dn13 = dn13.y, dn14 = dn14.y, dn15 = dn15.y,
    jin_label_sens = factor(jin_label_sens, levels = c("nonAKD","Recovery","Non-Recovery"))
  ) %>%
  dplyr::select(-dplyr::ends_with(".x"), -dplyr::ends_with(".y"))

dat_1y  <- longdat_sens %>% filter(years_from_time0 <= 1)
dat_3y  <- longdat_sens %>% filter(years_from_time0 <= 3)
dat_all <- longdat_sens

ctrl <- lmeControl(
  maxIter = 1e8, msMaxIter = 1e8,
  opt = "optim", optimMethod = "L-BFGS-B"
)

form_slope_adj <- egfr ~ years_from_time0 * jin_label_sens + time0_egfr_c - 1 +
  age_c + sex + arb_acei_use +
  dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 +
  dn9 + dn10 + dn12 + dn13 + dn14 + dn15 +
  years_from_time0:(age_c + arb_acei_use +
                      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 +
                      dn9 + dn10 + dn12 + dn13 + dn14 + dn15)

fit_sens_1y <- lme(
  fixed   = form_slope_adj,
  random  = list(id = pdSymm(~ 1 + years_from_time0)),
  data    = dat_1y,
  na.action = na.omit,
  method  = "REML",
  control = ctrl
)

fit_sens_3y <- lme(
  fixed   = form_slope_adj,
  random  = list(id = pdSymm(~ 1 + years_from_time0)),
  data    = dat_3y,
  na.action = na.omit,
  method  = "REML",
  control = ctrl
)

fit_sens_all <- lme(
  fixed   = form_slope_adj,
  random  = list(id = pdSymm(~ 1 + years_from_time0)),
  data    = dat_all,
  na.action = na.omit,
  method  = "REML",
  control = ctrl
)

# ==========================================================
# 5) Window-aligned OBSERVED trajectory builder (≤1y)
#    - metric="egfr"  : 実測
#    - metric="delta" : ΔeGFR（baseline=target_time==0の各id eGFR）
# ==========================================================
make_traj_sens_1y <- function(akd_time_sens, metric = c("egfr","delta"),
                              by_time = 0.25, width_multiplier = 0.75) {
  metric <- match.arg(metric)
  
  sens_for_window <- akd_time_sens %>%
    filter(jin_label_sens != "No-data") %>%
    filter(!is.na(egfr), !is.na(years_from_time0)) %>%
    filter(years_from_time0 >= 0, years_from_time0 <= 1) %>%
    mutate(jin_label_sens = as.character(jin_label_sens)) %>%
    filter(jin_label_sens %in% c("nonAKD","Recovery","Non-Recovery"))
  
  median_interval_sens <- sens_for_window %>%
    arrange(id, years_from_time0) %>%
    group_by(id) %>%
    summarise(d = diff(years_from_time0), .groups = "drop") %>%
    pull(d) %>%
    median(na.rm = TRUE)
  
  window_width_sens <- median_interval_sens * width_multiplier
  half_w_sens <- window_width_sens / 2
  target_timepoints <- seq(0, 1, by = by_time)
  
  dt_s <- as.data.table(sens_for_window)
  ids_s  <- unique(dt_s$id)
  grid_s <- CJ(id = ids_s, target_time = target_timepoints)
  
  setkey(dt_s, id, years_from_time0)
  
  window_obs_sens <- dt_s[
    grid_s,
    on = .(id, years_from_time0 = target_time),
    roll = "nearest",
    nomatch = 0L,
    .(id,
      target_time = i.target_time,
      years_from_time0,
      egfr,
      jin_label_sens)
  ]
  window_obs_sens <- window_obs_sens[abs(years_from_time0 - target_time) <= half_w_sens]
  
  # ---- metric 변換 ----
  if (metric == "delta") {
    base_dt <- window_obs_sens[target_time == 0, .(baseline_egfr = egfr[1]), by = id]
    window_obs_sens <- merge(window_obs_sens, base_dt, by = "id", all.x = FALSE, all.y = FALSE)
    window_obs_sens[, value := egfr - baseline_egfr]
  } else {
    window_obs_sens[, value := egfr]
  }
  
  traj <- as.data.frame(window_obs_sens) %>%
    as_tibble() %>%
    group_by(jin_label_sens, target_time) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      sd         = sd(value, na.rm = TRUE),
      n          = sum(!is.na(value)),
      se         = sd / sqrt(n),
      lwr        = mean_value - 1.96 * se,
      upr        = mean_value + 1.96 * se,
      .groups    = "drop"
    ) %>%
    mutate(
      jin_label_sens = factor(as.character(jin_label_sens),
                              levels = c("nonAKD","Recovery","Non-Recovery"))
    )
  
  traj
}

# ---- 色と凡例 ----
col_group <- c(nonAKD="#95A5A6", Recovery="#2ECC71", `Non-Recovery`="#E74C3C")
lab_group <- c(nonAKD="Non-AKD", Recovery="AKD with recovery", `Non-Recovery`="AKD without recovery")

# ---- 5A data（2パターン）----
traj_sens_1y_egfr  <- make_traj_sens_1y(akd_time_sens, metric = "egfr")
traj_sens_1y_delta <- make_traj_sens_1y(akd_time_sens, metric = "delta")

# ---- 5A plot maker ----
make_p_fig5A <- function(traj_df, y_lab) {
  ggplot(traj_df,
         aes(x = target_time, y = mean_value,
             color = jin_label_sens, group = jin_label_sens)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.3) +
    geom_errorbar(aes(ymin = lwr, ymax = upr),
                  width = 0.04, linewidth = 0.7) +
    labs(
      tag = "A",
      x = "Time from time0 (year)",
      y = y_lab,
      color = "Group"
    ) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0, 1)) +
    scale_color_manual(values = col_group, labels = lab_group, drop = FALSE) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.tag = element_text(face = "bold", size = 14),
      plot.tag.position = c(0, 0.98)
    )
}

p_fig5A_obs_egfr <- make_p_fig5A(traj_sens_1y_egfr,
                                 "Observed eGFR\n(mL/min/1.73 m²)")

p_fig5A_obs_delta <- make_p_fig5A(traj_sens_1y_delta,
                                  expression(paste("Observed ", Delta, "eGFR\n(mL/min/1.73 m"^2,")")))

# ==========================================================
# 6) Slopes + Differences (reuse functions)  ※あなたのまま
# ==========================================================
get_slopes_ci_sens <- function(fit, window_label){
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

sample_sizes_sens <- bind_rows(
  dat_1y  %>% mutate(window = "≤1 year"),
  dat_3y  %>% mutate(window = "≤3 years"),
  dat_all %>% mutate(window = "All period")
) %>%
  distinct(id, jin_label_sens, window) %>%
  count(jin_label_sens, window, name = "n") %>%
  transmute(group = as.character(jin_label_sens), window, n)

fits_sens <- list(
  "≤1 year"    = fit_sens_1y,
  "≤3 years"   = fit_sens_3y,
  "All period" = fit_sens_all
)

df_plot_sens <- imap_dfr(fits_sens, ~get_slopes_ci_sens(.x, .y))

df_diff_sens <- imap_dfr(fits_sens, ~get_diff_vs_nonakd_sens(.x, .y)) %>%
  bind_rows(
    tibble(
      window = names(fits_sens),
      group  = "nonAKD",
      diff_value = NA_real_,
      diff_lower = NA_real_,
      diff_upper = NA_real_,
      diff_p     = NA_real_
    )
  )

df_fig_sens <- df_plot_sens %>%
  left_join(sample_sizes_sens, by = c("group","window")) %>%
  left_join(df_diff_sens,      by = c("group","window")) %>%
  mutate(
    group  = factor(group, levels = c("nonAKD","Recovery","Non-Recovery")),
    window = factor(window, levels = c("≤1 year","≤3 years","All period"))
  )

# ==========================================================
# 7) Figure 5B (≤1y): downward bars  ※Differenceのみ表示（共通）
# ==========================================================
df_1y <- df_fig_sens %>%
  filter(window == "≤1 year") %>%
  mutate(group = factor(group, levels = c("nonAKD","Recovery","Non-Recovery")))

y_base <- min(df_1y$lower, na.rm = TRUE)
y_n    <- max(df_1y$upper, na.rm = TRUE) + 0.8
y_diff <- y_base - 5.2
y_pval <- y_base - 7.0

p_fig5B <- ggplot(df_1y, aes(x = group, y = estimate, fill = group)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.25, linewidth = 0.7) +
  geom_text(aes(y = y_n, label = paste0("n=", n)),
            size = 3.5, fontface = "bold", color = "grey20") +
  geom_text(aes(y = y_diff,
                label = ifelse(group == "nonAKD", "Reference",
                               sprintf("Diff: %.2f\n(%.2f, %.2f)",
                                       diff_value, diff_lower, diff_upper))),
            size = 3.1, lineheight = 0.95,
            fontface = "italic", color = "grey30") +
  geom_text(aes(y = y_pval,
                label = ifelse(group == "nonAKD", "",
                               case_when(
                                 is.na(diff_p) ~ "",
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
  scale_fill_manual(values = col_group, labels = lab_group) +
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

# ==========================================================
# 8) Combine (A/B) — 2 versions
# ==========================================================
fig5_egfr <- p_fig5A_obs_egfr / p_fig5B +
  plot_layout(heights = c(2.5, 5.5)) +
  plot_annotation(
    title = "Figure 5. Sensitivity analysis: observed eGFR trajectory and slope differences\nwithin 1 year after time0",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0,
                                margin = margin(b = 12))
    )
  )

fig5_delta <- p_fig5A_obs_delta / p_fig5B +
  plot_layout(heights = c(2.5, 5.5)) +
  plot_annotation(
    title = "Figure 5. Sensitivity analysis: observed ΔeGFR trajectory and slope differences\nwithin 1 year after time0",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0,
                                margin = margin(b = 12))
    )
  )

print(fig5_egfr)
print(fig5_delta)

# ==========================================================
# 9) Save (PDF + TIFF) — sensitivity_analysis folder
# ==========================================================
ggsave(
  filename = file.path(outdir, "Figure5A_sens_observed_eGFR_trajectory_and_slope_1y.pdf"),
  plot  = fig5_egfr,
  width = 230, height = 180, units = "mm",
  device = cairo_pdf
)

ggsave(
  filename = file.path(outdir, "Figure5A_sens_observed_eGFR_trajectory_and_slope_1y.tiff"),
  plot  = fig5_egfr,
  width = 230, height = 180, units = "mm",
  device = ragg::agg_tiff,
  dpi = 600, compression = "lzw"
)

ggsave(
  filename = file.path(outdir, "Figure5B_sens_observed_delta_eGFR_trajectory_and_slope_1y.pdf"),
  plot  = fig5_delta,
  width = 230, height = 180, units = "mm",
  device = cairo_pdf
)

ggsave(
  filename = file.path(outdir, "Figure5B_sens_observed_delta_eGFR_trajectory_and_slope_1y.tiff"),
  plot  = fig5_delta,
  width = 230, height = 180, units = "mm",
  device = ragg::agg_tiff,
  dpi = 600, compression = "lzw"
)

# ---- end ----
}#感度分析　1年以内　折れ線＋下向き棒グラフ
{
############################################################
# Supplemental Figure 2 (Sensitivity) — 2 versions
#  Ver-A) Observed eGFR trajectory (≤3y + All) + downward bar (Diff only)
#  Ver-B) Observed ΔeGFR trajectory (≤3y + All; baseline=time0) + same bar
#
# Save to: X:/R/sensitivity_analysis/
#   SupplementalFigure2A_sens_observed_eGFR_trajectory_and_bar_3y_all.(pdf/tiff)
#   SupplementalFigure2B_sens_observed_delta_eGFR_trajectory_and_bar_3y_all.(pdf/tiff)
############################################################

library(dplyr)
library(ggplot2)
library(patchwork)
library(data.table)
library(tidyr)
library(Cairo)
library(ragg)

# -----------------------------
# 0) color / legend（指定どおり）
# -----------------------------
col_group <- c(
  nonAKD         = "#95A5A6",
  Recovery       = "#2ECC71",
  `Non-Recovery` = "#E74C3C"
)
lab_group <- c(
  nonAKD         = "Non-AKD",
  Recovery       = "AKD with recovery",
  `Non-Recovery` = "AKD without recovery"
)

# -----------------------------
# Save folder（指定）
# -----------------------------
setwd("X:/R")
outdir <- file.path(getwd(), "sensitivity_analysis")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ==========================================================
# 1) A-panel: window-aligned trajectory function（egfr / delta両対応）
#    - dat_in: akd_time_sens（years_from_time0, egfr, id, jin_label_sens を含む）
#    - metric="egfr"  : 実測eGFR
#    - metric="delta" : ΔeGFR（baseline = target_time==0 の各id eGFR）
# ==========================================================
make_window_aligned_traj2 <- function(dat_in, end_time,
                                      by_time = 0.25,
                                      width_multiplier = 0.75,
                                      window_label = "≤3 years",
                                      metric = c("egfr","delta")) {
  metric <- match.arg(metric)
  
  # --- filter range + exclude No-data + keep 3 groups only ---
  sens_for_window <- dat_in %>%
    filter(jin_label_sens != "No-data") %>%
    filter(!is.na(egfr), !is.na(years_from_time0)) %>%
    filter(years_from_time0 >= 0, years_from_time0 <= end_time) %>%
    mutate(jin_label_sens = as.character(jin_label_sens)) %>%
    filter(jin_label_sens %in% c("nonAKD","Recovery","Non-Recovery"))
  
  # --- window width from within-id median interval ---
  median_interval <- sens_for_window %>%
    arrange(id, years_from_time0) %>%
    group_by(id) %>%
    summarise(d = diff(years_from_time0), .groups = "drop") %>%
    pull(d) %>%
    median(na.rm = TRUE)
  
  window_width <- median_interval * width_multiplier
  half_w <- window_width / 2
  
  # --- target time grid ---
  target_timepoints <- seq(0, end_time, by = by_time)
  
  # --- data.table rolling join (nearest) ---
  dt <- as.data.table(sens_for_window)
  ids  <- unique(dt$id)
  grid <- CJ(id = ids, target_time = target_timepoints)
  
  setkey(dt, id, years_from_time0)
  
  window_obs <- dt[
    grid,
    on = .(id, years_from_time0 = target_time),
    roll = "nearest",
    nomatch = 0L,
    .(id,
      target_time = i.target_time,
      years_from_time0,
      egfr,
      jin_label_sens)
  ]
  
  # keep within window band
  window_obs <- window_obs[abs(years_from_time0 - target_time) <= half_w]
  
  # --- metric transform ---
  if (metric == "delta") {
    base_dt <- window_obs[target_time == 0, .(baseline_egfr = egfr[1]), by = id]
    window_obs <- merge(window_obs, base_dt, by = "id", all.x = FALSE, all.y = FALSE)
    window_obs[, value := egfr - baseline_egfr]
  } else {
    window_obs[, value := egfr]
  }
  
  # --- mean ± 95%CI ---
  traj <- as.data.frame(window_obs) %>%
    as_tibble() %>%
    group_by(jin_label_sens, target_time) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      sd         = sd(value, na.rm = TRUE),
      n          = sum(!is.na(value)),
      se         = sd / sqrt(n),
      lwr        = mean_value - 1.96 * se,
      upr        = mean_value + 1.96 * se,
      .groups    = "drop"
    ) %>%
    mutate(
      window = window_label,
      jin_label_sens = factor(jin_label_sens, levels = c("nonAKD","Recovery","Non-Recovery"))
    )
  
  traj
}

# ==========================================================
# 2) Create trajectory data for ≤3y and All period（2パターン）
# ==========================================================
end_all <- max(akd_time_sens$years_from_time0, na.rm = TRUE)

# ---- Ver-A: eGFR ----
traj_3y_egfr <- make_window_aligned_traj2(
  dat_in = akd_time_sens, end_time = 3, by_time = 0.25,
  window_label = "≤3 years", metric = "egfr"
)
traj_all_egfr <- make_window_aligned_traj2(
  dat_in = akd_time_sens, end_time = end_all, by_time = 0.5,
  window_label = "All period", metric = "egfr"
)
traj_supp2_egfr <- bind_rows(traj_3y_egfr, traj_all_egfr) %>%
  mutate(window = factor(window, levels = c("≤3 years","All period")))

# ---- Ver-B: ΔeGFR ----
traj_3y_delta <- make_window_aligned_traj2(
  dat_in = akd_time_sens, end_time = 3, by_time = 0.25,
  window_label = "≤3 years", metric = "delta"
)
traj_all_delta <- make_window_aligned_traj2(
  dat_in = akd_time_sens, end_time = end_all, by_time = 0.5,
  window_label = "All period", metric = "delta"
)
traj_supp2_delta <- bind_rows(traj_3y_delta, traj_all_delta) %>%
  mutate(window = factor(window, levels = c("≤3 years","All period")))

# ==========================================================
# 3) Plot A maker（egfr / delta 共通）
# ==========================================================
make_p_supp2A <- function(traj_df, y_lab) {
  ggplot(traj_df,
         aes(x = target_time, y = mean_value,
             color = jin_label_sens, group = jin_label_sens)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.1) +
    geom_errorbar(aes(ymin = lwr, ymax = upr),
                  width = 0.04, linewidth = 0.7) +
    facet_grid(. ~ window, scales = "free_x") +
    labs(
      tag = "A",
      x = "Time from time0 (years)",
      y = y_lab,
      color = "Group"
    ) +
    scale_color_manual(values = col_group, labels = lab_group, drop = FALSE) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.tag = element_text(face = "bold", size = 14),
      plot.tag.position = c(0, 0.98)
    )
}

p_supp2A_egfr <- make_p_supp2A(traj_supp2_egfr,
                               "Observed eGFR\n(mL/min/1.73 m²)")

p_supp2A_delta <- make_p_supp2A(traj_supp2_delta,
                                expression(paste("Observed ", Delta, "eGFR\n(mL/min/1.73 m"^2,")")))

# ==========================================================
# 4) Plot B: downward bars for ≤3y + All period (Diff only) — 共通
#    - df_fig_sens を使用（あなたの既存成果物）
# ==========================================================
df_bar_supp2 <- df_fig_sens %>%
  filter(window %in% c("≤3 years","All period")) %>%
  mutate(
    group  = factor(as.character(group), levels = c("nonAKD","Recovery","Non-Recovery")),
    window = factor(as.character(window), levels = c("≤3 years","All period"))
  )

y_base <- min(df_bar_supp2$lower, na.rm = TRUE)
y_n    <- max(df_bar_supp2$upper, na.rm = TRUE) + 0.8
y_diff <- y_base - 5.2
y_pval <- y_base - 7.0

p_supp2B <- ggplot(df_bar_supp2, aes(x = group, y = estimate, fill = group)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.25, linewidth = 0.7) +
  geom_text(aes(y = y_n, label = paste0("n=", n)),
            size = 3.5, fontface = "bold", color = "grey20") +
  geom_text(aes(y = y_diff,
                label = ifelse(group == "nonAKD", "Reference",
                               sprintf("Diff: %.2f\n(%.2f, %.2f)",
                                       diff_value, diff_lower, diff_upper))),
            size = 3.1, lineheight = 0.95,
            fontface = "italic", color = "grey30") +
  geom_text(aes(y = y_pval,
                label = ifelse(group == "nonAKD", "",
                               case_when(
                                 is.na(diff_p) ~ "",
                                 diff_p < 0.001 ~ "p<0.001",
                                 diff_p < 0.01  ~ sprintf("p=%.3f", diff_p),
                                 TRUE           ~ sprintf("p=%.2f", diff_p)
                               ))),
            size = 3.0, fontface = "bold", color = "grey20") +
  facet_grid(. ~ window) +
  labs(
    tag = "B",
    x = NULL,
    y = "Mean change in eGFR\n(mL/min/1.73 m² per year)",
    fill = "Group"
  ) +
  scale_fill_manual(values = col_group, labels = lab_group, drop = FALSE) +
  scale_x_discrete(expand = expansion(add = 0.8)) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.6),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0, 0.98),
    plot.margin = margin(t = 5, r = 25, b = 20, l = 25)
  )

# ==========================================================
# 5) Combine + Save (2 versions)
# ==========================================================
supp_fig2_egfr <- p_supp2A_egfr / p_supp2B +
  plot_layout(heights = c(2.6, 5.4)) +
  plot_annotation(
    title = "Supplemental Figure 2. Sensitivity analysis: observed eGFR trajectory and slope differences\nwithin 3 years and all period after time0",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0,
                                margin = margin(b = 12))
    )
  )

supp_fig2_delta <- p_supp2A_delta / p_supp2B +
  plot_layout(heights = c(2.6, 5.4)) +
  plot_annotation(
    title = "Supplemental Figure 2. Sensitivity analysis: observed ΔeGFR trajectory and slope differences\nwithin 3 years and all period after time0",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0,
                                margin = margin(b = 12))
    )
  )

print(supp_fig2_egfr)
print(supp_fig2_delta)

ggsave(
  filename = file.path(outdir, "SupplementalFigure2A_sensitivity_observed_eGFR_trajectory_and_bar_3y_all.pdf"),
  plot  = supp_fig2_egfr,
  width = 230, height = 180, units = "mm",
  device = cairo_pdf
)

ggsave(
  filename = file.path(outdir, "SupplementalFigure2A_sensitivity_observed_eGFR_trajectory_and_bar_3y_all.tiff"),
  plot  = supp_fig2_egfr,
  width = 230, height = 180, units = "mm",
  device = ragg::agg_tiff,
  dpi = 600, compression = "lzw"
)

ggsave(
  filename = file.path(outdir, "SupplementalFigure2B_sensitivity_observed_delta_eGFR_trajectory_and_bar_3y_all.pdf"),
  plot  = supp_fig2_delta,
  width = 230, height = 180, units = "mm",
  device = cairo_pdf
)

ggsave(
  filename = file.path(outdir, "SupplementalFigure2B_sensitivity_observed_delta_eGFR_trajectory_and_bar_3y_all.tiff"),
  plot  = supp_fig2_delta,
  width = 230, height = 180, units = "mm",
  device = ragg::agg_tiff,
  dpi = 600, compression = "lzw"
)

# ---- end ----

}#3年以内と全期間の折れ線グラフ