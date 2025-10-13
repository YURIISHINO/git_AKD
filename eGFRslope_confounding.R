library(dplyr)
library(readr)
setwd("E:/R")
jin1_Eligibile <- read_csv("jin1_Eligibile.csv", locale = locale(encoding = "SHIFT-JIS"))
jin1_inclusion <- read_csv("jin1_inclusion.csv", locale = locale(encoding = "SHIFT-JIS"))
akd_time_m <- jin1_inclusion %>%
  mutate(
    years_from_time0 = as.numeric(date - time0) / 365.25  # 年単位に変換
  )%>%
  filter(years_from_time0 >= 0)
# CKD_statusを交絡因子に入れるため"nd" を NA に変換して factor 化
akd_time_m <- akd_time_m %>%
  mutate(CKD_status = na_if(CKD_status, "nd"),
         CKD_status = factor(CKD_status, levels = c("nonCKD", "CKD")))


library(dplyr)
library(nlme)

# 1) baseline_cov から index_cre は外す（time0_egfrに統一）
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

# 2) akd_time_m 側に同名列があるなら一旦消してから left_join
covars <- c("age","sex","arb_acei_use",
            "dn1","dn3","dn4","dn5","dn6","dn7","dn8","dn9","dn10","dn12","dn13","dn14","dn15")

longdat <- akd_time_m %>%
#  filter(years_from_time0 <= 1) %>%         #１年以内に言及してしまうと、混合効果モデルでの２年以内や３年以内がワークしないのでは？
  dplyr::select(-any_of(covars)) %>%         # ★ 重複候補を削除
  left_join(baseline_cov, by = "id") %>%
  mutate(
    age_c        = scale(age, center = TRUE, scale = FALSE)[,1],
    time0_egfr_c = scale(time0_egfr, center = TRUE, scale = FALSE)[,1],
    jin_label    = factor(jin_label, levels = c("nonAKD", "Recovery", "Non-Recovery"))
  )

# ①ベースモデル（time0_egfrのみ）
fit_base <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1,
  random = list(id = pdSymm(~ 1 + years_from_time0)),
  na.action = na.omit,
  data = longdat,
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, 
                       opt = "optim", optimMethod = "L-BFGS-B")
)

# ②交絡因子（レベルのみ調整）（Coxの交絡因子の主効果を追加）
fit_level_adj <- update(
  fit_base,
  . ~ . + age_c + sex + arb_acei_use + CKD_status +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + 
    dn10 + dn12 + dn13 + dn14 + dn15
)

# ③交絡因子（スロープまで調整）（years_from_time0との交互作用を追加）
fit_slope_adj <- update(
  fit_level_adj,
  . ~ . + years_from_time0:(age_c + arb_acei_use + CKD_status +
                              dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + 
                              dn10 + dn12 + dn13 + dn14 + dn15)
)

###下向き棒グラフを作成###################################################################
# 必要パッケージ
library(dplyr)
library(tidyr)
library(multcomp)   # glht()
library(ggplot2)
library(stringr)
library(nlme)
#1年以内、2年以内、3年以内のモデル作成
{
# 1年以内
fit_base_1y <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1,
  random = list(id = pdSymm(~ 1 + years_from_time0)),
  na.action = na.omit,
  data = filter(longdat, years_from_time0 <= 1),
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)
# 2年以内
fit_base_2y <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1,
  random = list(id = pdSymm(~ 1 + years_from_time0)),
  na.action = na.omit,
  data = filter(longdat, years_from_time0 <= 2),
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)
# 3年以内
fit_base_3y <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1,
  random = list(id = pdSymm(~ 1 + years_from_time0)),
  na.action = na.omit,
  data = filter(longdat, years_from_time0 <= 3),
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)

# 1年以内
fit_level_adj_1y <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1 +
    age_c + sex + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 +
    dn10 + dn12 + dn13 + dn14 + dn15,
  random = list(id = pdSymm(~ 1 + years_from_time0)),
  na.action = na.omit,
  data = filter(longdat, years_from_time0 <= 1),
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)

# 2年以内
fit_level_adj_2y <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1 +
    age_c + sex + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 +
    dn10 + dn12 + dn13 + dn14 + dn15,
  random = list(id = pdSymm(~ 1 + years_from_time0)),
  na.action = na.omit,
  data = filter(longdat, years_from_time0 <= 2),
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)

# 3年以内
fit_level_adj_3y <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1 +
    age_c + sex + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 +
    dn10 + dn12 + dn13 + dn14 + dn15,
  random = list(id = pdSymm(~ 1 + years_from_time0)),
  na.action = na.omit,
  data = filter(longdat, years_from_time0 <= 3),
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)

# 1年以内
fit_slope_adj_1y <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1 +
    age_c + sex + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 +
    dn10 + dn12 + dn13 + dn14 + dn15 +
    years_from_time0:(age_c + arb_acei_use +
                        dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 +
                        dn10 + dn12 + dn13 + dn14 + dn15),
  random = list(id = pdSymm(~ 1 + years_from_time0)),
  na.action = na.omit,
  data = filter(longdat, years_from_time0 <= 1),
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)

# 2年以内
fit_slope_adj_2y <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1 +
    age_c + sex + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 +
    dn10 + dn12 + dn13 + dn14 + dn15 +
    years_from_time0:(age_c + arb_acei_use +
                        dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 +
                        dn10 + dn12 + dn13 + dn14 + dn15),
  random = list(id = pdSymm(~ 1 + years_from_time0)),
  na.action = na.omit,
  data = filter(longdat, years_from_time0 <= 2),
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)

# 3年以内
fit_slope_adj_3y <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1 +
    age_c + sex + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 +
    dn10 + dn12 + dn13 + dn14 + dn15 +
    years_from_time0:(age_c + arb_acei_use +
                        dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 +
                        dn10 + dn12 + dn13 + dn14 + dn15),
  random = list(id = pdSymm(~ 1 + years_from_time0)),
  na.action = na.omit,
  data = filter(longdat, years_from_time0 <= 3),
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)

}

#---- 1) 各モデルから「群ごとのスロープ推定値＋95%CI」を取る関数 ----#
slope_ci_by_group <- function(fit, model_label, window_label){
  cf <- names(fixef(fit))
  v <- function(g){
    vec <- rep(0, length(cf)); names(vec) <- cf
    # 基本スロープ
    if("years_from_time0" %in% cf) vec["years_from_time0"] <- 1
    # 群×時間の交互作用（あれば加算）
    if(g == "Recovery" && "years_from_time0:jin_labelRecovery" %in% cf)
      vec["years_from_time0:jin_labelRecovery"] <- 1
    if(g == "Non-Recovery" && "years_from_time0:jin_labelNon-Recovery" %in% cf)
      vec["years_from_time0:jin_labelNon-Recovery"] <- 1
    # 交絡因子×時間がモデル③に入っていても、
    # それらの平均（中心化=0）でのスロープを表示する想定
    vec
  }
  L <- rbind(
    nonAKD       = v("nonAKD"),
    Recovery     = v("Recovery"),
    `Non-Recovery` = v("Non-Recovery")
  )
  ci <- suppressMessages(confint(glht(fit, linfct = L)))
  tibble(
    group  = rownames(L),
    estimate = ci$confint[,"Estimate"],
    lower    = ci$confint[,"lwr"],
    upper    = ci$confint[,"upr"],
    model  = model_label,
    window = window_label
  )
}

#---- 2) モデルの入れ物を用意（あなたのオブジェクト名に置き換え） ----#
# 例：各時間窓ごとに ①ベース ②レベル調整 ③スロープ調整 を入れる
fits <- list(
  "≤1 year" = list(
    `① Base`  = fit_base_1y,
    `② Level-adjusted` = fit_level_adj_1y,
    `③ Slope-adjusted` = fit_slope_adj_1y
  ),
  "≤2 years" = list(
    `① Base`  = fit_base_2y,
    `② Level-adjusted` = fit_level_adj_2y,
    `③ Slope-adjusted` = fit_slope_adj_2y
  ),
  "≤3 years" = list(
    `① Base`  = fit_base_3y,
    `② Level-adjusted` = fit_level_adj_3y,
    `③ Slope-adjusted` = fit_slope_adj_3y
  ),
  "All period" = list(
    `① Base`  = fit_base, 
    `② Level-adjusted` = fit_level_adj,
    `③ Slope-adjusted` = fit_slope_adj 
    )
)

#---- 3) ループで推定値を抽出して結合 ----#
df_plot <- purrr::imap_dfr(fits, function(models, win){
  purrr::imap_dfr(models, function(fit, mdl){
    slope_ci_by_group(fit, model_label = mdl, window_label = win)
  })
})

# 因子の順序（色や並びを固定）
df_plot <- df_plot %>%
  mutate(
    group  = factor(group, levels = c("nonAKD","Recovery","Non-Recovery")),
    model  = factor(model, levels = c("① Base","② Level-adjusted","③ Slope-adjusted")),
    window = factor(window, levels = c("≤1 year","≤2 years","≤3 years","All period"))
  )

#---- 4) 作図（モデル×時間窓の2次元ファセット） ----#
ggplot(df_plot, aes(x = group, y = estimate, fill = group)) +
  geom_col(width = 0.65) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 0.5) +
  facet_grid(model ~ window) +
  labs(
    title = "Estimated eGFR Slopes by Group across Models and Time Windows",
    x = NULL,
    y = expression(paste("Slope (mL/min/1.73 m"^2," per year), 95% CI")),
    fill = "Group"
  ) +
  scale_fill_manual(
    values = c(
      nonAKD = "#E41A1C",       # 赤
      Recovery = "#4DAF4A",     # 緑
      `Non-Recovery` = "#377EB8" # 青
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = NA),
    legend.position = "bottom"
  )

#95％信頼区間など######
slope_contrast_by_group <- function(fit, model_label, window_label, level = 0.95){
  cf <- names(fixef(fit))
  v <- function(g){
    vec <- rep(0, length(cf)); names(vec) <- cf
    if("years_from_time0" %in% cf) vec["years_from_time0"] <- 1
    if(g == "Recovery" && "years_from_time0:jin_labelRecovery" %in% cf)
      vec["years_from_time0:jin_labelRecovery"] <- 1
    if(g == "Non-Recovery" && "years_from_time0:jin_labelNon-Recovery" %in% cf)
      vec["years_from_time0:jin_labelNon-Recovery"] <- 1
    vec
  }
  b  <- v("nonAKD")
  r  <- v("Recovery")
  nr <- v("Non-Recovery")
  
  K <- rbind(
    `Recovery − nonAKD`       = r  - b,
    `Non-Recovery − nonAKD`   = nr - b,
    `Non-Recovery − Recovery` = nr - r
  )
  
  gl <- glht(fit, linfct = K)
  ci <- suppressMessages(confint(gl, level = level))
  sm <- suppressMessages(summary(gl))
  
  tibble(
    contrast = rownames(K),
    diff     = ci$confint[, "Estimate"],
    lower    = ci$confint[, "lwr"],
    upper    = ci$confint[, "upr"],
    p_value  = sm$test$pvalues,
    signif   = ifelse(ci$confint[, "lwr"] > 0 | ci$confint[, "upr"] < 0, "Yes", "No"),
    model    = model_label,
    window   = window_label
  )
}
library(purrr)

df_contrast <- purrr::imap_dfr(fits, function(models, win){
  purrr::imap_dfr(models, function(fit, mdl){
    slope_contrast_by_group(fit, model_label = mdl, window_label = win, level = 0.95)
  })
})

# 見やすく整形（丸め、星付けなど）
df_contrast_out <- df_contrast %>%
  mutate(
    diff  = round(diff, 2),
    lower = round(lower, 2),
    upper = round(upper, 2),
    p_value = signif(p_value, 3),
    ci = paste0("[", lower, ", ", upper, "]"),
    star = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ ""
    )
  ) %>%
  dplyr::select(window, model, contrast, diff, ci, p_value, star, signif)

df_contrast_out
print(df_contrast_out, n=Inf)

#抄録用グラフ＋保存####
library(dplyr)
library(ggplot2)

# --- 1) サブセット（≤1 year × 2モデル） ---
df_1y <- df_plot %>%
  filter(window == "≤1 year",
         model %in% c("① Base", "③ Slope-adjusted")) %>%
  mutate(
    group  = factor(group, levels = c("nonAKD", "Recovery", "Non-Recovery")),
    model  = recode(model,
                    "① Base" = "Unadjusted",
                    "③ Slope-adjusted" = "Slope-adjusted")
  )

# --- 2) 図作成 ---
p <- ggplot(df_1y, aes(x = group, y = estimate, fill = group)) +
  geom_col(width = 0.65) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 0.6) +
  facet_wrap(~ model, nrow = 1) +
  labs(
    title = "Comparison of one-year eGFR slopes\nby among recovery groups",
    x = NULL,
    y = expression(paste("Slope (mL/min/1.73 m"^2," per year), 95% CI")),
    fill = "Group"
  ) +
  scale_fill_manual(values = c(
    nonAKD = "#E41A1C",
    Recovery = "#4DAF4A",
    `Non-Recovery` = "#377EB8"
  )) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(size = 10.5, face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_blank(),       # ★ X軸ラベルを非表示
    axis.ticks.x = element_blank(),      # ★ 目盛り線も非表示
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.margin = margin(t = 20, r = 10, b = 20, l = 15),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5, lineheight = 1.1)
  )

# --- 3) 抄録用に保存（JPEG 640×480 px） ---
ggsave(
  filename = "egfr_slope_1y_two_models_no_xlabels.jpeg",
  plot = p,
  width = 640, height = 480, units = "px",
  dpi = 150, device = "jpeg", quality = 95
)
