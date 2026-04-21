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
  
  # ============================
  # ★ Diff と p-value の「縦位置」をここで決める（間隔ノブ）
  # ============================
  diff_y <- y_min - 3.2
  p_y    <- y_min - 5.2   # ← ここを大きくすると、pがさらに下に行ってDiffと離れる
  
  p <- ggplot(df_in, aes(x = group, y = estimate, fill = group)) +
    geom_col(width = 0.7, alpha = 0.85) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, linewidth = 0.7) +
    geom_text(aes(y = pmax(upper, 0) + 0.5, label = paste0("n=", n)),
              size = 3.5, fontface = "bold", color = "grey20") +
    labs(
      x = NULL,
      # ★ expressionをやめて空白問題も回避（m²を使う）
      y = "Mean change in eGFR (mL/min/1.73 m² per year)",
      fill = "Group"
    ) +
    scale_fill_manual(
      values = c(nonAKD = "#95A5A6", Recovery = "#2ECC71", `Non-Recovery` = "#E74C3C"),
      labels = c(
        nonAKD = "Non AKD",
        Recovery = "AKD with recovery",
        `Non-Recovery` = "AKD without recovery"
      )
    ) +
    scale_x_discrete(expand = expansion(add = 0.6)) +
    # ★ ここが「pを下にした分、下端で切れないようにする」修正点
    coord_cartesian(ylim = c(y_min - 6.8, y_max + 1.5), clip = "off") +
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
  
  # ============================
  # ★ Difference vs nonAKD（Referenceもここで入れる）
  # ============================
  if (show_diff) {
    p <- p +
      geom_text(
        aes(label = ifelse(group == "nonAKD",
                           "Reference",
                           sprintf("Diff: %.2f\n(%.2f, %.2f)", diff_value, diff_lower, diff_upper))),
        y = diff_y,                       # ★ここ：y_min - 3.6 → diff_y に置換
        size = 3, lineheight = 0.95,
        fontface = "italic", color = "grey30"
      )
  }
  
  # ============================
  # ★ p-value text（Diffと離すために p_y を使う）
  # ============================
  if (show_p) {
    p <- p +
      geom_text(
        aes(label = ifelse(group == "nonAKD", "",
                           dplyr::case_when(
                             is.na(diff_p) ~ "",
                             diff_p < 0.001 ~ "p<0.001",
                             diff_p < 0.01  ~ sprintf("p=%.3f", diff_p),
                             TRUE           ~ sprintf("p=%.2f", diff_p)
                           ))),
        y = p_y,                          # ★ここ：y_min - 4.7 → p_y に置換
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
ggsave(file.path(outdir, "SupplementFigure2_eGFR_slope_adj_3windows.pdf"),
       plot = p_supp_fig1,
       width = w_mm, height = h_mm, units = "mm",
       device = cairo_pdf, dpi = 300)

ggsave(file.path(outdir, "SupplementFigure2_eGFR_slope_adj_3windows.tiff"),
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
    title = "Observed change in eGFR from time0",
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
    title = "Observed change in eGFR from time0",
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
    title = "Observed change in eGFR from time0",
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
    title = "Observed change in eGFR from time0",
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
    title = "Adjusted differences in annual eGFR change",
    y     = "Mean change in eGFR\n(mL/min/1.73 m² per year)"
  ) +
  scale_x_discrete(expand = expansion(add = 0.8)) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0),
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0, 1.00),
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
    plot.tag.position = c(0, 1.00)
  )

fig3_egfr <- p_line_egfr / p_fig4_for_fig3 +
  plot_layout(heights = c(2.5, 5.5))

print(fig3_egfr)

#ggsave(
#  filename = file.path(outdir, "Figure3A_observed_eGFR_trajectory_and_slope_1y.pdf"),
#  plot  = fig3_egfr,
#  width = 230, height = 180, units = "mm",
#  device = cairo_pdf
#)

#ggsave(
#  filename = file.path(outdir, "Figure3A_observed_eGFR_trajectory_and_slope_1y.tiff"),
#  plot  = fig3_egfr,
#  width = 230, height = 180, units = "mm",
#  device = ragg::agg_tiff,
#  dpi = 600, compression = "lzw"
#)

# =========================================================
# 3) Figure 3-2：Observed ΔeGFR（変化量） + bar
# =========================================================
base_family <- "sans"   # ← Arialにしたいなら "Arial"

common_theme <- theme(
  text = element_text(family = base_family),
  plot.title = element_text(family = base_family, face = "bold"),
  plot.tag   = element_text(family = base_family, face = "bold", size = 14)
)
p_line_delta <- p_obs_delta_1y_bar +
  labs(
    tag = "A",
    x = "Time from time0 (year)",
    y = expression(paste(Delta, "eGFR(mL/min/1.73 m"^2, ")"))
  ) +
  theme(
    plot.margin = margin(b = 5),
    legend.position = "right",
    plot.tag.position = c(0, 1.00)
  ) +
  common_theme

  fig3_delta <- p_line_delta / p_fig4_for_fig3 +
    plot_layout(heights = c(2.5, 5.5))
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

p_supp1A_delta <- make_Apanel(
  traj_supp1A_delta,
  expression(paste(Delta, "eGFR (mL/min/1.73 m"^2, ")"))
)
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
  labs(tag = "A", title = "Observed change in eGFR from time0") +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0),
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0, 1.00)
  )

p_supp1B_tag <- p_supp1B +
  labs(tag = "B", title = "Adjusted differences in annual eGFR change") +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0),
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0, 1.00)
  )

supp_fig1_main_egfr <- p_supp1A_egfr_tag / p_supp1B_tag +
  plot_layout(heights = c(2.3, 3.3)) 

# --- Ver-B: ΔeGFR ---
p_supp1A_delta_tag <- p_supp1A_delta +
  labs(tag = "A", title = "Observed change in eGFR from time0") +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0),
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0, 1.00)
  )

supp_fig1_main_delta <- p_supp1A_delta_tag / p_supp1B_tag +
  plot_layout(heights = c(2.3, 3.3)) 

print(supp_fig1_main_egfr)
print(supp_fig1_main_delta)
# どれでもOK（Windowsなら "sans" が無難。見た目を固定したいなら "Arial" 等）
base_family <- "sans"

title_theme <- theme(
  text = element_text(family = base_family),
  plot.title = element_text(family = base_family, size = 12, face = "bold", hjust = 0),
  plot.tag   = element_text(family = base_family, face = "bold", size = 14)
)

p_supp1A_egfr_tag  <- p_supp1A_egfr_tag  + title_theme
p_supp1A_delta_tag <- p_supp1A_delta_tag + title_theme
p_supp1B_tag       <- p_supp1B_tag       + title_theme
# ==========================================================
# 6) Save (PDF + TIFF) — eGFRslope folder
# ==========================================================
ggsave(
  filename = file.path(outdir, "SupplementalFigure2_main_observed_delta_eGFR_trajectory_and_bar_3y_all.pdf"),
  plot  = supp_fig1_main_delta,
  width = 230, height = 180, units = "mm",
  device = cairo_pdf
)

ggsave(
  filename = file.path(outdir, "SupplementalFigure2_main_observed_delta_eGFR_trajectory_and_bar_3y_all.tiff"),
  plot  = supp_fig1_main_delta,
  width = 230, height = 180, units = "mm",
  device = ragg::agg_tiff,
  dpi = 600, compression = "lzw"
)
# ---- end ----
}#3年以内と全期間
{
  # =========================================================
  # 4) Figure 3B only：poster-friendly version
  # =========================================================
  outdir <- "X:/R/eGFRslope"
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  fig3b_only <- p_fig4 +
    labs(
      tag   = "B",
      title = "Adjusted differences in annual eGFR decline",
      y     = "Mean change in eGFR\n(mL/min/1.73 m² per year)"
    ) +
    scale_x_discrete(expand = expansion(add = 0.8)) +
    theme(
      text = element_text(family = base_family),
      plot.title = element_text(size = 12, face = "bold", hjust = 0, family = base_family),
      plot.tag = element_text(face = "bold", size = 14, family = base_family),
      plot.tag.position = c(0, 1.00),
      legend.position = "bottom",
      plot.margin = margin(t = 8, r = 20, b = 8, l = 20)
    )
  
  print(fig3b_only)
  
  ggsave(
    filename = file.path(outdir, "Figure3B_only_slope_difference_1y_poster.pdf"),
    plot  = fig3b_only,
    width = 180, height = 120, units = "mm",
    device = cairo_pdf
  )
  
  ggsave(
    filename = file.path(outdir, "Figure3B_only_slope_difference_1y_poster.tiff"),
    plot  = fig3b_only,
    width = 180, height = 120, units = "mm",
    device = ragg::agg_tiff,
    dpi = 600,
    compression = "lzw"
  )
}#ポスター用
R.version
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

# -----------------------------
# FONT/Title unify (Figure5)
# -----------------------------
BASE_SIZE_UNIFIED <- 13
common_theme_all <- theme(
  text = element_text(family = "sans"),
  plot.title = element_text(size = 12, face = "bold", hjust = 0),
  
  plot.tag   = element_text(face = "bold", size = 14),
  plot.tag.position = c(0.06, 0.985),   # ★xを0→0.06へ（右へ）＋少し下へ
  
  # 念のため：y軸タイトルも少し左へ
  axis.title.y = element_text(margin = margin(r = 18)),
  
  # 左余白も確保（tagとy軸が詰まるのを防ぐ）
  plot.margin = margin(t = 4, r = 5, b = 5, l = 18)
)

# ---- 5A plot maker ----
make_p_fig5A <- function(traj_df, y_lab, main_title = NULL) {
  ggplot(traj_df,
         aes(x = target_time, y = mean_value,
             color = jin_label_sens, group = jin_label_sens)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.3) +
    geom_errorbar(aes(ymin = lwr, ymax = upr),
                  width = 0.04, linewidth = 0.7) +
    labs(
      tag   = "A",
      title = main_title,
      x     = "Time from time0 (year)",
      y     = y_lab,
      color = "Group"
    ) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0, 1)) +
    scale_color_manual(values = col_group, labels = lab_group, drop = FALSE) +
    theme_bw(base_size = BASE_SIZE_UNIFIED) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "right") +
    common_theme_all
}

p_fig5A_obs_egfr <- make_p_fig5A(
  traj_sens_1y_egfr,
  "Observed eGFR\n(mL/min/1.73 m²)",
  main_title = NULL
)

p_fig5A_obs_delta <- make_p_fig5A(
  traj_sens_1y_delta,
  "\u0394eGFR (mL/min/1.73 m\u00B2)",     # ← 1行のΔeGFR(単位)
  main_title = "Observed change in eGFR from time0"
)
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
    tag   = "B",
    title = "Adjusted differences in annual eGFR change",
    x = NULL,
    y = "Mean change in eGFR\n(mL/min/1.73 m² per year)",
    fill = "Group"
  ) +
  scale_fill_manual(values = col_group, labels = lab_group, drop = FALSE) +
  scale_x_discrete(expand = expansion(add = 0.8)) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = BASE_SIZE_UNIFIED) +
  theme(
    panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.6),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(t = 5, r = 25, b = 20, l = 25)
  ) +
  common_theme_all
fig5_delta <- p_fig5A_obs_delta / p_fig5B +
  plot_layout(heights = c(2.5, 5.5))
# ==========================================================
# 8) Combine (A/B) — 2 versions
# ==========================================================
fig5_egfr <- p_fig5A_obs_egfr / p_fig5B +
  plot_layout(heights = c(2.5, 5.5)) 

fig5_delta <- p_fig5A_obs_delta / p_fig5B +
  plot_layout(heights = c(2.5, 5.5)) 

print(fig5_egfr)
print(fig5_delta)

# ==========================================================
# 9) Save (PDF + TIFF) — sensitivity_analysis folder
# ==========================================================
ggsave(
  filename = file.path(outdir, "Supp_Figure3_sens_observed_delta_eGFR_trajectory_and_slope_1y.pdf"),
  plot  = fig5_delta,
  width = 230, height = 180, units = "mm",
  device = cairo_pdf
)

ggsave(
  filename = file.path(outdir, "Supp_Figure3_sens_observed_delta_eGFR_trajectory_and_slope_1y.tiff"),
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
      
      # --- tag（A）: 左上ギリギリを避ける ---
      plot.tag = element_text(face = "bold", size = 14),
      plot.tag.position = c(0.02, 0.985),   # ★ここがポイント（Figure3寄せ）
      
      # --- y軸タイトル: 左へ逃す（ΔeGFRとAの干渉を防ぐ） ---
      axis.title.y = element_text(margin = margin(r = 18)),
      
      # --- 外側余白: 左と上を少し増やす（見た目をFigure3に寄せる） ---
      plot.margin = margin(t = 4, r = 5, b = 5, l = 18)
    )
}

p_supp2A_egfr <- make_p_supp2A(traj_supp2_egfr,
                               "Observed eGFR\n(mL/min/1.73 m²)")

p_supp2A_delta <- make_p_supp2A(
  traj_supp2_delta,
  "ΔeGFR (mL/min/1.73 m²)"
)

p_supp2A_delta <- p_supp2A_delta +
  labs(title = "Observed change in eGFR from time0")
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
                label = ifelse(group == "nonAKD", "Reference",
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
    title = "Adjusted differences in annual eGFR change",
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
    plot.tag.position = c(0, 1.00),
    plot.margin = margin(t = 5, r = 25, b = 20, l = 25)
  )
common_title_theme <- theme(
  plot.title = element_text(size = 12, face = "bold", hjust = 0)
)

p_supp2A_delta <- p_supp2A_delta + common_title_theme
p_supp2B       <- p_supp2B       + common_title_theme
# ==========================================================
# 5) Combine + Save (2 versions)
# ==========================================================
supp_fig2_egfr <- p_supp2A_egfr / p_supp2B +
  plot_layout(heights = c(2.6, 5.4)) 

supp_fig2_delta <- p_supp2A_delta / p_supp2B +
  plot_layout(heights = c(2.6, 5.4)) 

print(supp_fig2_egfr)
print(supp_fig2_delta)

ggsave(
  filename = file.path(outdir, "SupplementalFigure4_sensitivity_observed_delta_eGFR_trajectory_and_bar_3y_all.pdf"),
  plot  = supp_fig2_delta,
  width = 230, height = 180, units = "mm",
  device = cairo_pdf
)

ggsave(
  filename = file.path(outdir, "SupplementalFigure4_sensitivity_observed_delta_eGFR_trajectory_and_bar_3y_all.tiff"),
  plot  = supp_fig2_delta,
  width = 230, height = 180, units = "mm",
  device = ragg::agg_tiff,
  dpi = 600, compression = "lzw"
)

# ---- end ----

}#3年以内と全期間の折れ線グラフ


#time0をindex_date+210に合わせる
{
  ############################################################
  # Observed delta eGFR trajectory + annual slope
  # Landmark-aligned version
  # time0 = index_date + 210 days
  # Top: A Observed delta eGFR trajectory (≤1 year)
  # Bottom: B Adjusted annual eGFR slope bar plot (≤1 year)
  # Save to: X:/R/eGFRslope
  ############################################################
  
  graphics.off()
  
  # =========================
  # Packages
  # =========================
  pkgs <- c(
    "readr","dplyr","tidyr","stringr","purrr",
    "nlme","multcomp","ggplot2","data.table",
    "grid","gridExtra","gtable","ragg"
  )
  to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
  
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(nlme)
  library(multcomp)
  library(ggplot2)
  library(data.table)
  library(grid)
  library(gridExtra)
  library(gtable)
  library(ragg)
  
  # =========================
  # 0) Read data
  # =========================
  setwd("X:/R")
  jin1_Eligibile <- read_csv("jin1_Eligibile.csv", locale = locale(encoding = "SHIFT-JIS"))
  
  # =========================
  # 1) Build inclusion data
  #    time0 fixed at index_date + 210
  # =========================
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
      jin_label = factor(jin_label, levels = c("nonAKD", "Recovery", "Non-Recovery")),
      index_plus_210 = index_date + 210
    ) %>%
    filter(!is.na(jin_label))
  
  # baseline eGFR nearest to day 210 within 150-210 days
  time0_egfr_df <- jin1_inclusion %>%
    mutate(
      days_from_index = as.numeric(date - index_date),
      dist_to_210 = abs(days_from_index - 210)
    ) %>%
    filter(days_from_index >= 150, days_from_index <= 210) %>%
    group_by(id) %>%
    arrange(dist_to_210, desc(date), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      id,
      time0_egfr = egfr,
      time0_egfr_date = date
    )
  
  jin1_inclusion <- jin1_inclusion %>%
    mutate(
      time0 = index_plus_210,
      years_from_time0 = as.numeric(date - time0) / 365.25
    ) %>%
    left_join(time0_egfr_df, by = "id")
  
  # =========================
  # 2) Longitudinal data
  # =========================
  akd_time_m <- jin1_inclusion %>%
    filter(years_from_time0 >= 0)
  
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
  
  covars <- c(
    "age","sex","arb_acei_use",
    "dn1","dn3","dn4","dn5","dn6","dn7","dn8","dn9","dn10","dn12","dn13","dn14","dn15"
  )
  
  longdat <- akd_time_m %>%
    dplyr::select(-any_of(covars)) %>%
    left_join(baseline_cov, by = "id") %>%
    mutate(
      age_c        = as.numeric(scale(age, center = TRUE, scale = FALSE)),
      time0_egfr_c = as.numeric(scale(time0_egfr, center = TRUE, scale = FALSE)),
      jin_label    = factor(jin_label, levels = c("nonAKD", "Recovery", "Non-Recovery"))
    ) %>%
    filter(!is.na(time0_egfr))
  
  # =========================
  # 3) Slope-adjusted LME (≤1 year only)
  # =========================
  ctrl <- lmeControl(
    maxIter = 1e8, msMaxIter = 1e8,
    opt = "optim", optimMethod = "L-BFGS-B"
  )
  
  fit_slope_adj_1y <- lme(
    egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1 +
      age_c + sex + arb_acei_use +
      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15 +
      years_from_time0:(age_c + arb_acei_use +
                          dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15),
    random = list(id = pdSymm(~ 1 + years_from_time0)),
    data = longdat %>% filter(years_from_time0 <= 1),
    na.action = na.omit,
    method = "REML",
    control = ctrl
  )
  
  # =========================
  # 4) Extract annual slopes + contrast
  # =========================
  cf <- names(fixef(fit_slope_adj_1y))
  
  v_slope <- function(g){
    vec <- rep(0, length(cf)); names(vec) <- cf
    if ("years_from_time0" %in% cf) vec["years_from_time0"] <- 1
    if (g == "Recovery" && "years_from_time0:jin_labelRecovery" %in% cf)
      vec["years_from_time0:jin_labelRecovery"] <- 1
    if (g == "Non-Recovery" && "years_from_time0:jin_labelNon-Recovery" %in% cf)
      vec["years_from_time0:jin_labelNon-Recovery"] <- 1
    vec
  }
  
  L <- rbind(
    nonAKD         = v_slope("nonAKD"),
    Recovery       = v_slope("Recovery"),
    `Non-Recovery` = v_slope("Non-Recovery")
  )
  
  ci_slope <- suppressMessages(confint(glht(fit_slope_adj_1y, linfct = L)))
  
  df_bar <- tibble(
    group    = rownames(L),
    estimate = ci_slope$confint[, "Estimate"],
    lower    = ci_slope$confint[, "lwr"],
    upper    = ci_slope$confint[, "upr"]
  ) %>%
    mutate(group = factor(group, levels = c("nonAKD","Recovery","Non-Recovery")))
  
  # contrasts vs nonAKD
  b  <- v_slope("nonAKD")
  r  <- v_slope("Recovery")
  nr <- v_slope("Non-Recovery")
  
  K <- rbind(
    `Recovery − nonAKD`     = r  - b,
    `Non-Recovery − nonAKD` = nr - b
  )
  
  gl_contrast <- glht(fit_slope_adj_1y, linfct = K)
  ci_contrast <- suppressMessages(confint(gl_contrast))
  sm_contrast <- suppressMessages(summary(gl_contrast))
  
  df_contrast <- tibble(
    group      = c("Recovery","Non-Recovery"),
    diff_value = ci_contrast$confint[, "Estimate"],
    diff_lower = ci_contrast$confint[, "lwr"],
    diff_upper = ci_contrast$confint[, "upr"],
    diff_p     = sm_contrast$test$pvalues
  )
  
  sample_sizes <- longdat %>%
    filter(years_from_time0 <= 1) %>%
    distinct(id, jin_label) %>%
    count(jin_label, name = "n") %>%
    transmute(group = as.character(jin_label), n)
  
  df_bar <- df_bar %>%
    left_join(sample_sizes, by = "group") %>%
    left_join(df_contrast, by = "group")
  
  # =========================
  # 5) Observed delta eGFR data
  # =========================
  dt <- as.data.table(
    akd_time_m %>%
      filter(years_from_time0 >= 0, years_from_time0 <= 1) %>%
      dplyr::select(id, years_from_time0, egfr, jin_label)
  )
  
  dt <- dt[!is.na(egfr)]
  
  median_interval <- dt[order(id, years_from_time0),
                        .(d = diff(years_from_time0)), by = id]$d %>%
    median(na.rm = TRUE)
  
  window_width <- median_interval * 0.75
  half_w <- window_width / 2
  
  target_timepoints <- seq(0, 1, by = 0.25)
  ids <- unique(dt$id)
  grid_dt <- CJ(id = ids, target_time = target_timepoints)
  
  setkey(dt, id, years_from_time0)
  
  window_data_obs <- dt[
    grid_dt,
    on = .(id, years_from_time0 = target_time),
    roll = "nearest",
    nomatch = 0L,
    .(id,
      target_time = i.target_time,
      years_from_time0,
      egfr,
      jin_label)
  ]
  
  window_data_obs <- window_data_obs[abs(years_from_time0 - target_time) <= half_w]
  
  base_dt <- window_data_obs[target_time == 0, .(baseline_egfr = egfr[1]), by = id]
  
  wd2 <- merge(window_data_obs, base_dt, by = "id", all.x = FALSE, all.y = FALSE)
  wd2[, delta_egfr := egfr - baseline_egfr]
  
  plot_dat <- as_tibble(wd2) %>%
    group_by(jin_label, target_time) %>%
    summarise(
      mean_delta = mean(delta_egfr, na.rm = TRUE),
      sd         = sd(delta_egfr, na.rm = TRUE),
      n          = sum(!is.na(delta_egfr)),
      se         = sd / sqrt(n),
      lwr        = mean_delta - 1.96 * se,
      upr        = mean_delta + 1.96 * se,
      .groups    = "drop"
    ) %>%
    mutate(
      jin_label = factor(as.character(jin_label),
                         levels = c("nonAKD","Recovery","Non-Recovery"))
    )
  
  # =========================
  # 6) Colors / labels
  # =========================
  values_group <- c(
    nonAKD         = "#95A5A6",
    Recovery       = "#2ECC71",
    `Non-Recovery` = "#E74C3C"
  )
  
  labels_group <- c(
    nonAKD         = "non-AKD",
    Recovery       = "AKD with recovery",
    `Non-Recovery` = "AKD without recovery"
  )
  
  # =========================
  # 7) Panel A axis range
  #    make sure all CI/lines fit inside panel
  # =========================
  top_ymin <- floor(min(plot_dat$lwr, na.rm = TRUE) - 0.5)
  top_ymax <- ceiling(max(plot_dat$upr, na.rm = TRUE) + 0.5)
  
  # keep a sensible upper bound near reference style
  top_ymax <- max(top_ymax, 0)
  top_ymin <- min(top_ymin, -10)
  
  # =========================
  # 8) Panel A
  # =========================
  p_top <- ggplot(
    plot_dat,
    aes(x = target_time, y = mean_delta, color = jin_label, group = jin_label)
  ) +
    #geom_hline(yintercept = 0, linewidth = 0.5, color = "black") +
    geom_line(linewidth = 1.0) +
    geom_point(size = 2.4) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.03, linewidth = 0.8) +
    scale_color_manual(values = values_group, labels = labels_group, name = "Group") +
    scale_x_continuous(
      breaks = seq(0, 1, by = 0.25),
      limits = c(0, 1),
      labels = sprintf("%.2f", seq(0, 1, by = 0.25)),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      limits = c(top_ymin, top_ymax),
      expand = expansion(mult = c(0.03, 0.05))
    ) +
    labs(
      title = "A  Observed change in eGFR from landmark",
      x = "Time from landmark (year)",
      y = expression(paste(Delta, "eGFR(mL/min/1.73 m"^2, ")"))
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "#D9D9D9", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.title = element_text(size = 12, color = "black"),
      axis.text  = element_text(size = 10, color = "black"),
      plot.title = element_text(size = 13, face = "bold", hjust = 0),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text  = element_text(size = 10),
      legend.key   = element_rect(fill = "white", color = NA),
      plot.margin  = margin(t = 4, r = 4, b = 2, l = 4)
    )
  
  # =========================
  # 9) Panel B positions/range
  #    keep CI and text inside panel
  # =========================
  df_bar <- df_bar %>%
    mutate(
      x = c(1, 2, 3),
      diff_label = case_when(
        group == "nonAKD" ~ "Reference",
        TRUE ~ sprintf("Diff: %.2f\n(%.2f, %.2f)", diff_value, diff_lower, diff_upper)
      ),
      p_label = case_when(
        group == "nonAKD" ~ "",
        diff_p < 0.001 ~ "p<0.001",
        diff_p < 0.01  ~ sprintf("p=%.3f", diff_p),
        TRUE           ~ sprintf("p=%.2f", diff_p)
      )
    )
  
  bar_ymin <- floor(min(df_bar$lower, na.rm = TRUE) - 4.5)
  bar_ymax <- max(1.6, ceiling(max(df_bar$upper, na.rm = TRUE) + 0.8))
  
  n_y    <- 0.55
  diff_y <- bar_ymin + 1.9
  p_y    <- bar_ymin + 0.5
  
  # =========================
  # 10) Panel B
  # =========================
  p_bottom <- ggplot(df_bar, aes(x = x, y = estimate, fill = group)) +
    #geom_hline(yintercept = 0, linewidth = 0.5, color = "black") +
    geom_col(width = 0.62, color = NA) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.16, linewidth = 0.8) +
    geom_text(aes(y = n_y, label = paste0("n=", n)),
              size = 3.8, fontface = "bold", color = "#333333") +
    geom_text(aes(y = diff_y, label = diff_label),
              size = 3.2, fontface = "italic", color = "#4D4D4D", lineheight = 0.95) +
    geom_text(aes(y = p_y, label = p_label),
              size = 3.4, fontface = "bold", color = "#333333") +
    scale_fill_manual(
      values = c(
        nonAKD = "#95A5A6",
        Recovery = "#2ECC71",
        `Non-Recovery` = "#E74C3C"
      ),
      breaks = c("nonAKD", "Recovery", "Non-Recovery"),
      labels = c(
        nonAKD = "Non AKD",
        Recovery = "AKD with recovery",
        `Non-Recovery` = "AKD without recovery"
      ),
      name = "Group"
    ) +
    scale_x_continuous(
      breaks = NULL,
      limits = c(0.4, 3.6),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      limits = c(bar_ymin, bar_ymax),
      expand = expansion(mult = c(0.03, 0.03))
    ) +
    labs(
      title = "B  Adjusted differences in annual eGFR change",
      x = NULL,
      y = "Mean change in eGFR\n(mL/min/1.73 m² per year)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.title = element_text(size = 12, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(size = 13, face = "bold", hjust = 0),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text  = element_text(size = 10),
      legend.key.size = unit(0.9, "lines"),
      legend.key   = element_rect(fill = "white", color = NA),
      plot.margin  = margin(t = 2, r = 4, b = 2, l = 4)
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))
  
  # =========================
  # 11) Combine
  # =========================
  fig_combined <- arrangeGrob(
    p_top,
    p_bottom,
    ncol = 1,
    heights = c(0.9, 1.1)
  )
  
  grid.newpage()
  grid.draw(fig_combined)
  
  # =========================
  # 12) Save
  # =========================
  outdir <- "X:/R/eGFRslope"
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  ragg::agg_tiff(
    filename = file.path(outdir, "observed_delta_and_slope_1y_landmark210_v2.tiff"),
    width = 180, height = 190, units = "mm", res = 600, compression = "lzw"
  )
  grid.newpage()
  grid.draw(fig_combined)
  dev.off()
  
  ggsave(
    filename = file.path(outdir, "observed_delta_and_slope_1y_landmark210_v2.pdf"),
    plot = fig_combined,
    width = 180, height = 190, units = "mm",
    device = cairo_pdf
  )
}


{
  ############################################################
  # Observed delta eGFR trajectory + annual slope
  # Landmark-aligned version
  # time0 = index_date + 210 days
  # Layout: Supplement Figure 4 style
  #   Top    : A Observed change in eGFR from landmark
  #   Bottom : B Adjusted differences in annual eGFR change
  #   Columns: ≤3 years / All period
  # Save to: X:/R/eGFRslope
  ############################################################
  
  graphics.off()
  
  # =========================
  # Packages
  # =========================
  pkgs <- c(
    "readr","dplyr","tidyr","stringr","purrr",
    "nlme","multcomp","ggplot2","data.table",
    "grid","gridExtra","gtable","ragg",
    "patchwork","forcats"
  )
  to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
  
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(nlme)
  library(multcomp)
  library(ggplot2)
  library(data.table)
  library(grid)
  library(gridExtra)
  library(gtable)
  library(ragg)
  library(patchwork)
  library(forcats)
  
  # =========================
  # 0) Read data
  # =========================
  setwd("X:/R")
  jin1_Eligibile <- read_csv(
    "jin1_Eligibile.csv",
    locale = locale(encoding = "SHIFT-JIS")
  )
  
  # =========================
  # 1) Build inclusion data
  #    time0 fixed at index_date + 210
  # =========================
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
      jin_label = factor(jin_label, levels = c("nonAKD", "Recovery", "Non-Recovery")),
      index_plus_210 = index_date + 210
    ) %>%
    filter(!is.na(jin_label))
  
  # baseline eGFR nearest to day 210 within 150-210 days
  time0_egfr_df <- jin1_inclusion %>%
    mutate(
      days_from_index = as.numeric(date - index_date),
      dist_to_210 = abs(days_from_index - 210)
    ) %>%
    filter(days_from_index >= 150, days_from_index <= 210) %>%
    group_by(id) %>%
    arrange(dist_to_210, desc(date), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      id,
      time0_egfr = egfr,
      time0_egfr_date = date
    )
  
  jin1_inclusion <- jin1_inclusion %>%
    mutate(
      time0 = index_plus_210,
      years_from_time0 = as.numeric(date - time0) / 365.25
    ) %>%
    left_join(time0_egfr_df, by = "id")
  
  # =========================
  # 2) Longitudinal data
  # =========================
  akd_time_m <- jin1_inclusion %>%
    filter(years_from_time0 >= 0)
  
  baseline_cov <- jin1_Eligibile %>%
    filter(exclude == "include") %>%
    group_by(id) %>%
    arrange(index_date, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(arb_acei_use = if_else(arb == 1 | acei == 1, 1L, 0L)) %>%
    dplyr::select(
      id, age, sex, arb_acei_use,
      dn1, dn3, dn4, dn5, dn6, dn7, dn8, dn9, dn10, dn12, dn13, dn14, dn15
    )
  
  covars <- c(
    "age","sex","arb_acei_use",
    "dn1","dn3","dn4","dn5","dn6","dn7","dn8","dn9","dn10","dn12","dn13","dn14","dn15"
  )
  
  longdat <- akd_time_m %>%
    dplyr::select(-any_of(covars)) %>%
    left_join(baseline_cov, by = "id") %>%
    mutate(
      age_c        = as.numeric(scale(age, center = TRUE, scale = FALSE)),
      time0_egfr_c = as.numeric(scale(time0_egfr, center = TRUE, scale = FALSE)),
      jin_label    = factor(jin_label, levels = c("nonAKD", "Recovery", "Non-Recovery"))
    ) %>%
    filter(!is.na(time0_egfr))
  
  # =========================
  # 3) Common settings
  # =========================
  ctrl <- lmeControl(
    maxIter   = 1e8,
    msMaxIter = 1e8,
    opt = "optim",
    optimMethod = "L-BFGS-B"
  )
  
  values_group <- c(
    nonAKD         = "#95A5A6",
    Recovery       = "#2ECC71",
    `Non-Recovery` = "#E74C3C"
  )
  
  labels_group <- c(
    nonAKD         = "Non-AKD",
    Recovery       = "AKD with recovery",
    `Non-Recovery` = "AKD without recovery"
  )
  
  # =========================
  # 4) Functions
  # =========================
  
  # ---- fit slope-adjusted model and extract bar data ----
  make_bar_data <- function(dat, years_max = NULL, panel_label = "≤3 years") {
    
    dat_sub <- dat
    if (!is.null(years_max)) {
      dat_sub <- dat_sub %>% filter(years_from_time0 <= years_max)
    }
    
    fit <- lme(
      egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1 +
        age_c + sex + arb_acei_use +
        dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15 +
        years_from_time0:(age_c + arb_acei_use +
                            dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15),
      random = list(id = pdSymm(~ 1 + years_from_time0)),
      data = dat_sub,
      na.action = na.omit,
      method = "REML",
      control = ctrl
    )
    
    cf <- names(fixef(fit))
    
    v_slope <- function(g){
      vec <- rep(0, length(cf)); names(vec) <- cf
      if ("years_from_time0" %in% cf) vec["years_from_time0"] <- 1
      if (g == "Recovery" && "years_from_time0:jin_labelRecovery" %in% cf) {
        vec["years_from_time0:jin_labelRecovery"] <- 1
      }
      if (g == "Non-Recovery" && "years_from_time0:jin_labelNon-Recovery" %in% cf) {
        vec["years_from_time0:jin_labelNon-Recovery"] <- 1
      }
      vec
    }
    
    L <- rbind(
      nonAKD         = v_slope("nonAKD"),
      Recovery       = v_slope("Recovery"),
      `Non-Recovery` = v_slope("Non-Recovery")
    )
    
    ci_slope <- suppressMessages(confint(glht(fit, linfct = L)))
    
    df_bar <- tibble(
      group    = rownames(L),
      estimate = ci_slope$confint[, "Estimate"],
      lower    = ci_slope$confint[, "lwr"],
      upper    = ci_slope$confint[, "upr"]
    ) %>%
      mutate(group = factor(group, levels = c("nonAKD","Recovery","Non-Recovery")))
    
    # contrasts vs nonAKD
    b  <- v_slope("nonAKD")
    r  <- v_slope("Recovery")
    nr <- v_slope("Non-Recovery")
    
    K <- rbind(
      `Recovery − nonAKD`     = r  - b,
      `Non-Recovery − nonAKD` = nr - b
    )
    
    gl_contrast <- glht(fit, linfct = K)
    ci_contrast <- suppressMessages(confint(gl_contrast))
    sm_contrast <- suppressMessages(summary(gl_contrast))
    
    df_contrast <- tibble(
      group      = c("Recovery","Non-Recovery"),
      diff_value = ci_contrast$confint[, "Estimate"],
      diff_lower = ci_contrast$confint[, "lwr"],
      diff_upper = ci_contrast$confint[, "upr"],
      diff_p     = sm_contrast$test$pvalues
    )
    
    sample_sizes <- dat_sub %>%
      distinct(id, jin_label) %>%
      count(jin_label, name = "n") %>%
      transmute(group = as.character(jin_label), n)
    
    df_bar <- df_bar %>%
      left_join(sample_sizes, by = "group") %>%
      left_join(df_contrast, by = "group") %>%
      mutate(
        panel = panel_label,
        group = factor(group, levels = c("nonAKD","Recovery","Non-Recovery")),
        x = c(1, 2, 3),
        diff_label = case_when(
          group == "nonAKD" ~ "Reference",
          TRUE ~ sprintf("Diff: %.2f\n(%.2f, %.2f)", diff_value, diff_lower, diff_upper)
        ),
        p_label = case_when(
          group == "nonAKD" ~ "Reference",
          diff_p < 0.001 ~ "p<0.001",
          diff_p < 0.01  ~ sprintf("p=%.3f", diff_p),
          TRUE           ~ sprintf("p=%.2f", diff_p)
        )
      )
    
    list(fit = fit, df_bar = df_bar)
  }
  
  # ---- observed delta eGFR ----
  make_observed_data <- function(dat, years_max = NULL, panel_label = "≤3 years") {
    
    dt0 <- dat
    if (!is.null(years_max)) {
      dt0 <- dt0 %>% filter(years_from_time0 <= years_max)
    }
    
    dt <- as.data.table(
      dt0 %>%
        dplyr::select(id, years_from_time0, egfr, jin_label)
    )
    
    dt <- dt[!is.na(egfr)]
    
    if (!is.null(years_max)) {
      # 3-year panel: every 0.25 year
      target_timepoints <- seq(0, years_max, by = 0.25)
      window_width <- 0.25 * 0.75
    } else {
      # all-period panel: every 0.5 year
      max_time <- floor(max(dt$years_from_time0, na.rm = TRUE) * 2) / 2
      max_time <- max(max_time, 12)
      target_timepoints <- seq(0, max_time, by = 0.5)
      window_width <- 0.5 * 0.75
    }
    
    half_w <- window_width / 2
    
    ids <- unique(dt$id)
    grid_dt <- CJ(id = ids, target_time = target_timepoints)
    
    setkey(dt, id, years_from_time0)
    
    window_data_obs <- dt[
      grid_dt,
      on = .(id, years_from_time0 = target_time),
      roll = "nearest",
      nomatch = 0L,
      .(id,
        target_time = i.target_time,
        years_from_time0,
        egfr,
        jin_label)
    ]
    
    window_data_obs <- window_data_obs[abs(years_from_time0 - target_time) <= half_w]
    
    base_dt <- window_data_obs[target_time == 0, .(baseline_egfr = egfr[1]), by = id]
    wd2 <- merge(window_data_obs, base_dt, by = "id", all.x = FALSE, all.y = FALSE)
    wd2[, delta_egfr := egfr - baseline_egfr]
    
    as_tibble(wd2) %>%
      group_by(jin_label, target_time) %>%
      summarise(
        mean_delta = mean(delta_egfr, na.rm = TRUE),
        sd         = sd(delta_egfr, na.rm = TRUE),
        n          = sum(!is.na(delta_egfr)),
        se         = sd / sqrt(n),
        lwr        = mean_delta - 1.96 * se,
        upr        = mean_delta + 1.96 * se,
        .groups    = "drop"
      ) %>%
      mutate(
        jin_label = factor(as.character(jin_label),
                           levels = c("nonAKD","Recovery","Non-Recovery")),
        panel = panel_label
      )
  }
  
  # ---- top panel ----
  plot_top_panel <- function(plot_dat, y_limits, show_legend = FALSE) {
    
    max_x <- max(plot_dat$target_time, na.rm = TRUE)
    
    x_breaks <- if (max_x <= 3.01) {
      c(0, 1, 2, 3)
    } else {
      c(0, 3, 6, 9, 12)
    }
    
    ggplot(
      plot_dat,
      aes(x = target_time, y = mean_delta, color = jin_label, group = jin_label)
    ) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 1.8) +
      geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.06, linewidth = 0.6) +
      scale_color_manual(
        values = values_group,
        labels = labels_group,
        name = "Group"
      ) +
      scale_x_continuous(
        breaks = x_breaks,
        limits = c(min(x_breaks), max(x_breaks)),
        expand = expansion(mult = c(0.01, 0.02))
      ) +
      scale_y_continuous(
        limits = y_limits,
        breaks = pretty(y_limits, n = 4),
        expand = expansion(mult = c(0.02, 0.04))
      ) +
      labs(
        x = "Time from time0 (years)",
        y = expression(paste(Delta, "eGFR (mL/min/1.73 m"^2, ")"))
      ) +
      theme_bw(base_size = 12) +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "#D9D9D9", linewidth = 0.35),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.6),
        strip.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 12, color = "black"),
        axis.text  = element_text(size = 10, color = "black"),
        legend.position = if (show_legend) "right" else "none",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text  = element_text(size = 10),
        legend.key = element_rect(fill = "white", color = NA),
        plot.margin = margin(t = 3, r = 4, b = 2, l = 4)
      )
  }
  
  # ---- bottom panel ----
  plot_bottom_panel <- function(df_bar, y_limits, show_legend = FALSE) {
    
    n_y    <- y_limits[2] - 0.10 * diff(y_limits)
    diff_y <- y_limits[1] + 0.18 * diff(y_limits)
    p_y    <- y_limits[1] + 0.04 * diff(y_limits)
    
    ggplot(df_bar, aes(x = x, y = estimate, fill = group)) +
      geom_col(width = 0.62, color = NA) +
      geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.16, linewidth = 0.7) +
      geom_text(aes(y = n_y, label = paste0("n=", n)),
                size = 3.5, fontface = "bold", color = "#333333") +
      geom_text(aes(y = diff_y, label = diff_label),
                size = 3.0, fontface = "italic", color = "#4D4D4D", lineheight = 0.95) +
      geom_text(aes(y = p_y, label = p_label),
                size = 3.1, fontface = "bold", color = "#333333") +
      scale_fill_manual(
        values = values_group,
        breaks = c("nonAKD", "Recovery", "Non-Recovery"),
        labels = labels_group,
        name = "Group"
      ) +
      scale_x_continuous(
        breaks = NULL,
        limits = c(0.4, 3.6),
        expand = expansion(mult = c(0.02, 0.02))
      ) +
      scale_y_continuous(
        limits = y_limits,
        breaks = c(-8, -4, 0, 4),
        expand = expansion(mult = c(0.02, 0.02))
      ) +
      labs(
        x = NULL,
        y = "Mean change in eGFR\n(mL/min/1.73 m² per year)"
      ) +
      theme_bw(base_size = 12) +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.6),
        strip.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = if (show_legend) "bottom" else "none",
        legend.direction = "horizontal",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text  = element_text(size = 10),
        legend.key.size = unit(0.9, "lines"),
        legend.key   = element_rect(fill = "white", color = NA),
        plot.margin  = margin(t = 2, r = 4, b = 2, l = 4)
      ) +
      guides(fill = guide_legend(nrow = 1, byrow = TRUE))
  }
  
  # =========================
  # 5) Create data for ≤3 years / All period
  # =========================
  res_3y  <- make_bar_data(longdat, years_max = 3,    panel_label = "≤3 years")
  res_all <- make_bar_data(longdat, years_max = NULL, panel_label = "All period")
  
  bar_dat <- bind_rows(res_3y$df_bar, res_all$df_bar)
  
  obs_3y  <- make_observed_data(akd_time_m, years_max = 3,    panel_label = "≤3 years")
  obs_all <- make_observed_data(akd_time_m, years_max = NULL, panel_label = "All period")
  
  plot_dat <- bind_rows(obs_3y, obs_all)
  
  # =========================
  # 6) Common axis ranges
  # =========================
  top_ymin <- floor(min(plot_dat$lwr, na.rm = TRUE) - 1)
  top_ymax <- ceiling(max(plot_dat$upr, na.rm = TRUE) + 1)
  
  # Supplement Figure 4 に寄せて少し固定
  top_ymin <- min(top_ymin, -20)
  top_ymax <- max(top_ymax, 0)
  top_ylim <- c(top_ymin, top_ymax)
  
  bar_ymin <- floor(min(bar_dat$lower, na.rm = TRUE) - 3.5)
  bar_ymax <- ceiling(max(bar_dat$upper, na.rm = TRUE) + 2.0)
  
  bar_ymin <- min(bar_ymin, -10)
  bar_ymax <- max(bar_ymax, 4)
  bar_ylim <- c(bar_ymin, bar_ymax)
  
  # =========================
  # 7) Build top row
  # =========================
  p_top_3y <- plot_top_panel(
    plot_dat %>% filter(panel == "≤3 years"),
    y_limits = top_ylim,
    show_legend = FALSE
  ) +
    ggtitle(NULL) +
    facet_wrap(~panel, nrow = 1)
  
  p_top_all <- plot_top_panel(
    plot_dat %>% filter(panel == "All period"),
    y_limits = top_ylim,
    show_legend = TRUE
  ) +
    ggtitle(NULL) +
    facet_wrap(~panel, nrow = 1)
  
  # facet label を残すため 1枚ずつ作る
  g_top_3y  <- p_top_3y
  g_top_all <- p_top_all
  
  top_row <- g_top_3y + g_top_all + plot_layout(widths = c(1, 1.05))
  
  # =========================
  # 8) Build bottom row
  # =========================
  p_bottom_3y <- plot_bottom_panel(
    bar_dat %>% filter(panel == "≤3 years"),
    y_limits = bar_ylim,
    show_legend = FALSE
  ) +
    facet_wrap(~panel, nrow = 1)
  
  p_bottom_all <- plot_bottom_panel(
    bar_dat %>% filter(panel == "All period"),
    y_limits = bar_ylim,
    show_legend = TRUE
  ) +
    facet_wrap(~panel, nrow = 1)
  
  bottom_row <- p_bottom_3y + p_bottom_all + plot_layout(widths = c(1, 1.05))
  
  # =========================
  # 9) Add row titles (A / B)
  # =========================
  title_top <- ggplot() +
    annotate("text", x = 0, y = 1, hjust = 0, vjust = 1,
             label = "A  Observed change in eGFR from time0",
             size = 5, fontface = "bold") +
    theme_void() +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme(plot.margin = margin(t = 0, r = 0, b = -5, l = 2))
  
  title_bottom <- ggplot() +
    annotate("text", x = 0, y = 1, hjust = 0, vjust = 1,
             label = "B  Adjusted differences in annual eGFR change",
             size = 5, fontface = "bold") +
    theme_void() +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme(plot.margin = margin(t = 0, r = 0, b = -5, l = 2))
  
  # =========================
  # 10) Combine final figure
  # =========================
  fig_combined <-
    title_top /
    top_row /
    title_bottom /
    bottom_row +
    plot_layout(heights = c(0.06, 1.0, 0.06, 1.05))
  
  # =========================
  # 11) Draw
  # =========================
  print(fig_combined)
  
  # =========================
  # 12) Save
  # =========================
  outdir <- "X:/R/eGFRslope"
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # TIFF
  ragg::agg_tiff(
    filename = file.path(outdir, "SupplementFigure_landmark210_observed_delta_and_slope_3y_all.tiff"),
    width = 260, height = 180, units = "mm", res = 600, compression = "lzw"
  )
  print(fig_combined)
  dev.off()
  
  # PDF
  ggsave(
    filename = file.path(outdir, "SupplementFigure_landmark210_observed_delta_and_slope_3y_all.pdf"),
    plot = fig_combined,
    width = 260, height = 180, units = "mm",
    device = cairo_pdf
  )
}