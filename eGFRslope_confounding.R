{
#eGFR・線形混合効果モデルを走らせる、また解析に必要なcodeのみ####
{
library(nlme)
library(dplyr)
library(ggplot2)
# ① jin_label の定義と必要列の抽出
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
    jin_label = factor(jin_label, levels = c("nonAKD","Recovery","Non-Recovery"))
  ) #%>%
#dplyr::select(id, date, index_date, egfr, jin_status, jin_label)

jin1_inclusion %>%
  group_by(jin_status, jin_label) %>%
  summarise(
    unique_ids = n_distinct(id),
    .groups = "drop"
  )

jin1_inclusion %>%
  filter(!is.na(jin_label)) %>%
  group_by(jin_label) %>%
  summarise(unique_ids = n_distinct(id)) %>%
  arrange(desc(unique_ids))

# ② 90～210日の最大egfr日（max_egfr_date_210）抽出
egfr_max_date_210 <- jin1_inclusion %>%
  mutate(days_from_index = as.numeric(date - index_date)) %>%
  filter(days_from_index >= 90, days_from_index <= 210) %>%
  group_by(id) %>%
  filter(egfr == max(egfr, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(id, max_egfr_date_210 = date)

# ③ 0～90日の最新日（nearest_date_90）抽出
max_date_90 <- jin1_inclusion %>%
  mutate(days_from_index = as.numeric(date - index_date)) %>%
  filter(days_from_index >= 0, days_from_index <= 90) %>%
  group_by(id) %>%
  filter(days_from_index == max(days_from_index, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(id, nearest_date_90 = date)

# ④ 元データに結合
jin1_inclusion <- jin1_inclusion %>%
  left_join(egfr_max_date_210, by = "id") %>%
  left_join(max_date_90, by = "id")

# ⑤ time0 の定義（優先度：max_egfr_date_210 > nearest_date_90 > index_date）
jin1_inclusion <- jin1_inclusion %>%
  mutate(time0 = case_when(
    !is.na(max_egfr_date_210) ~ max_egfr_date_210,
    !is.na(nearest_date_90) ~ nearest_date_90,
    TRUE ~ index_date
  ))

# ⑥ time0 からの年差を計算
jin1_inclusion <- jin1_inclusion %>%
  mutate(
    days_from_time0 = as.numeric(difftime(date, time0, units = "days")),
    years_from_time0 = round(days_from_time0 / 365.25, 3)
  )

# ⑦ time0 に一致する egfr（複数あれば最大値）を抽出
time0_egfr_df <- jin1_inclusion %>%
  filter(date == time0) %>%
  group_by(id, date) %>%
  slice_max(egfr, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(id, time0 = date, time0_egfr = egfr)

# ⑧ id ごとに time0_egfr を付加
jin1_inclusion <- jin1_inclusion %>%
  left_join(time0_egfr_df %>% dplyr::select(id, time0_egfr), by = "id")
colnames(jin1_inclusion)

# CSV保存
library(readr)
write_csv(
  jin1_inclusion,
  "jin1_inclusion.csv"   # 保存ファイル名（作業ディレクトリに保存されます）
) 
}

library(dplyr)
library(readr)
library(nlme)
library(tidyr)
setwd("E:/R")
jin1_Eligibile <- read_csv("jin1_Eligibile.csv", locale = locale(encoding = "SHIFT-JIS"))

# AKD_status × jin_status のクロス集計(確認用)
dat_id <- jin1_Eligibile %>%
  distinct(id, .keep_all = TRUE)
# AKD_status × jin_status のクロス集計
table_id <- dat_id %>%
  count(AKD_status, jin_status) %>%
  pivot_wider(
    names_from = jin_status,
    values_from = n,
    values_fill = 0
  )
table_id

jin1_inclusion <- read_csv("jin1_inclusion.csv", locale = locale(encoding = "SHIFT-JIS"))
colnames(jin1_inclusion)
# AKD_status × jin_status のクロス集計
dat_id_inclusion <- jin1_inclusion %>%
  distinct(id, .keep_all = TRUE)
table_id_inclusion <- dat_id_inclusion %>%
  count(AKD_status, jin_status) %>%
  pivot_wider(
    names_from = jin_status,
    values_from = n,
    values_fill = 0
  )
table_id_inclusion

akd_time_m <- jin1_inclusion %>%
  mutate(
    years_from_time0 = as.numeric(date - time0) / 365.25  # 年単位に変換
  )%>%
  filter(years_from_time0 >= 0)
# CKD_statusを交絡因子に入れるため"nd" を NA に変換して factor 化
akd_time_m <- akd_time_m %>%
  mutate(CKD_status = na_if(CKD_status, "nd"),
         CKD_status = factor(CKD_status, levels = c("nonCKD", "CKD")))

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

ctrl <- lmeControl(
  maxIter  = 50,
  msMaxIter = 50,
  opt = "optim",
  optimMethod = "L-BFGS-B"
)


#下向き棒グラフを作成
# 必要パッケージ
library(dplyr)
library(tidyr)
library(multcomp)   # glht()
library(ggplot2)
library(stringr)
library(nlme)
fit_base <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1,
  random = list(id = pdSymm(~ 1 + years_from_time0)),
  na.action = na.omit,
  data = longdat,
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, 
                       opt = "optim", optimMethod = "L-BFGS-B")
  #control = ctrl
)

fit_level_adj <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1 +
    age_c + sex + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 +
    dn10 + dn12 + dn13 + dn14 + dn15,
  random = list(id = pdSymm(~ 1 + years_from_time0)),
  na.action = na.omit,
  data = longdat,
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, 
                       opt = "optim", optimMethod = "L-BFGS-B")
  #control = ctrl
)
fit_slope_adj <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1 +
    age_c + sex + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 +
    dn10 + dn12 + dn13 + dn14 + dn15 +
    years_from_time0:(age_c + arb_acei_use +
                        dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 +
                        dn10 + dn12 + dn13 + dn14 + dn15),
  random = list(id = pdSymm(~ 1 + years_from_time0)),
  na.action = na.omit,
  data = longdat,
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, 
                       opt = "optim", optimMethod = "L-BFGS-B")
  #control = ctrl
)
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
  #control = ctrl
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
  #control = ctrl
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
  #control = ctrl
)
# 1年以内
fit_base_1y <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1,
  random = list(id = pdSymm(~ 1 + years_from_time0)),
  na.action = na.omit,
  data = filter(longdat, years_from_time0 <= 1),
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, 
                       opt = "optim", optimMethod = "L-BFGS-B")
  #control = ctrl
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
  #control = ctrl
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
  #control = ctrl
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
  #control = ctrl
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
  #control = ctrl
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
  #control = ctrl
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
  #control = ctrl
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
  #control = ctrl
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
  #control = ctrl
)

}

{library(multcomp)
  library(nlme)
  
  
  ## 1) まずは 1年以内のスロープ調整済みモデルを指定
  fit <- fit_slope_adj_1y
  
  ## 2) 係数名を取得
  cf <- names(fixef(fit))
  
  ## 3) ある群 g の「スロープ（years_from_time0 の係数の合成）」を表すベクトルを作る
  v <- function(g){
    vec <- rep(0, length(cf)); names(vec) <- cf
    
    # 基本スロープ
    if("years_from_time0" %in% cf) vec["years_from_time0"] <- 1
    
    # 群×時間の交互作用（あれば加算）
    if(g == "Recovery" && "years_from_time0:jin_labelRecovery" %in% cf)
      vec["years_from_time0:jin_labelRecovery"] <- 1
    
    if(g == "Non-Recovery" && "years_from_time0:jin_labelNon-Recovery" %in% cf)
      vec["years_from_time0:jin_labelNon-Recovery"] <- 1
    
    # 交絡因子×時間は中心化＝0 で評価する前提
    vec
  }
  
  ## 4) 各群のスロープを表す行列 L（あなたの slope_ci_by_group と同じイメージ）
  L <- rbind(
    nonAKD        = v("nonAKD"),
    Recovery      = v("Recovery"),
    `Non-Recovery`= v("Non-Recovery")
  )
  
  ## 5) スロープ差のコントラスト行列 K を作る（3群のペアワイズ比較）
  K <- rbind(
    "Recovery vs nonAKD"        = L["Recovery", ]      - L["nonAKD", ],
    "Non-Recovery vs nonAKD"    = L["Non-Recovery", ]  - L["nonAKD", ],
    "Non-Recovery vs Recovery"  = L["Non-Recovery", ]  - L["Recovery", ]
  )
  
  ## 6) glht で Wald 検定（必要なら多重比較補正も）
  g_slope_diff <- glht(fit, linfct = K)
  
  # 推定値・標準誤差・z値・p値
  summary(g_slope_diff)             
  
  # 95%CI（差分の信頼区間）
  confint(g_slope_diff)
}#1年以内モデルで統計学的評価する

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

library(dplyr)
library(stringr)

# --- ラベル正規化 → 因子化（ここが肝） ---
df_plot <- df_plot %>%
  mutate(
    # 1) model を英語へ統一（丸数字を除去/置換）
    model = case_when(
      model %in% c("① Base", "Base")                         ~ "Base",
      model %in% c("② Level-adjusted", "Level-adjusted")     ~ "Level-adjusted",
      model %in% c("③ Slope-adjusted", "Slope-adjusted")     ~ "Slope-adjusted",
      TRUE ~ NA_character_
    ),
    # 2) window の表記ゆれ吸収（必要なら追加）
    window = case_when(
      window %in% c("≤1 year","<=1 year","≤1year","<=1year")       ~ "≤1 year",
      window %in% c("≤2 years","<=2 years","≤2years","<=2years")   ~ "≤2 years",
      window %in% c("≤3 years","<=3 years","≤3years","<=3years")   ~ "≤3 years",
      window %in% c("All period","All-period","All Period")        ~ "All period",
      TRUE ~ NA_character_
    ),
    # 3) group のダッシュ統一
    group  = trimws(gsub("\u2013|\u2212", "-", as.character(group)))
  ) %>%
  # 4) 欠損や想定外を落とす
  filter(!is.na(model), !is.na(window), group %in% c("nonAKD","Recovery","Non-Recovery")) %>%
  # 5) 最終的に因子化（ここで初めて factor()）
  mutate(
    model  = factor(model,  levels = c("Base","Level-adjusted","Slope-adjusted")),
    window = factor(window, levels = c("≤1 year","≤2 years","≤3 years","All period")),
    group  = factor(group,  levels = c("nonAKD","Recovery","Non-Recovery"))
  ) %>%
  # （もし完全重複があれば念のため除去）
  distinct(window, model, group, .keep_all = TRUE)

#1～3年と全期間
p_all <- 
  ggplot(df_plot, aes(x = group, y = estimate, fill = group)) +
  geom_col(width = 0.65) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 0.5) +
  facet_grid(model ~ window, drop = FALSE) +
  labs(
    title = "Estimated eGFR Slopes \nby Group across Models and Time Windows",
    x = NULL,
    y = expression(paste("Slope (mL/min/1.73 m"^2," per year), 95% CI")),
    fill = "Group"
  ) +
  scale_fill_manual(values = c(nonAKD="#E41A1C", Recovery="#4DAF4A", `Non-Recovery`="#377EB8")) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = NA),
    legend.position = "bottom",
    axis.text.x  = element_blank()  # 軸ラベル非表示
  )

#---- 1年以内と3年以内，Slope-adjustedモデルだけ ----#
p_1y_3y <- 
  ggplot(
    df_plot %>% 
      filter(
        window %in% c("≤1 year", "≤3 years"),
        model  == "Slope-adjusted"          # ★Slope-adjustedだけ抽出
      ),
    aes(x = group, y = estimate, fill = group)
  ) +
  geom_col(width = 0.65) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.2, linewidth = 0.5
  ) +
  facet_grid(
    . ~ window,         # ★model方向のfacetは削除して，windowだけで並べる
    drop = TRUE
  ) +
  labs(
    title = "Estimated eGFR Slopes (≤1 year and ≤3 years)\nSlope-adjusted model",
    x = NULL,
    y = expression(paste("Slope (mL/min/1.73 m"^2," per year), 95% CI")),
    fill = "Group"
  ) +
  scale_fill_manual(values = c(
    nonAKD       = "#E41A1C",
    Recovery     = "#4DAF4A",
    `Non-Recovery` = "#377EB8"
  )) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.placement  = "outside",
    legend.position  = "bottom",
    axis.text.x      = element_blank()
  )


#---- 1年以内と全期間のみ ----#
p_1y_all <- 
  ggplot(df_plot %>% filter(window %in% c("≤1 year", "All period")),
         aes(x = group, y = estimate, fill = group)) +
  geom_col(width = 0.65) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.2, linewidth = 0.5) +
  facet_grid(
    model ~ window,
    drop   = TRUE,
    switch = "y"      # ←★モデル名を左側に表示
  ) +
  labs(
    title = "Estimated eGFR Slopes (≤1 year and All period)",
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
    strip.placement = "outside",   # ★strip を外側に出す
    legend.position = "bottom",
    axis.text.x = element_blank()
  )




#画像保存#####
# 必要パッケージ
install.packages(c("ggplot2","Cairo"))     # 未インストールなら
library(ggplot2)
library(Cairo)
# ── 保存（論文向け推奨：PDF(ベクター) と TIFF(600 dpi)） ─────────
# 仕上がりサイズ：幅 180 mm, 高さ 120 mm（2段組誌に汎用）
w_mm <- 180; h_mm <- 120
outdir <- "E:/R"
## CairoでPDF（フォント埋め込み）; Windowsなら 'device = "cairo_pdf"' が安定
ggsave(file.path(outdir, "Figure4_all_windows.pdf"),
       plot = p_all, width = w_mm, height = h_mm, units = "mm",
       device = "pdf", dpi = 300)

ggsave(file.path(outdir, "Figure4_1y_and_≤3y_period.pdf"),
       plot = p_1y_3y, width = w_mm, height = h_mm, units = "mm",
       device = "pdf", dpi = 300)

## TIFF（600 dpi, 圧縮LZW）
# ragg があれば高品質レンダリング推奨：
use_ragg <- requireNamespace("ragg", quietly = TRUE)

if (use_ragg) {
  ggsave(file.path(outdir, "Figure4_all_windows.tiff"),
         plot = p_all, width = w_mm, height = h_mm, units = "mm",
         dpi = 600, device = ragg::agg_tiff, compression = "lzw")
  ggsave(file.path(outdir, "Figure4_1y_and_≤3y_period.tiff"),
         plot = p_1y_3y, width = w_mm, height = h_mm, units = "mm",
         dpi = 600, device = ragg::agg_tiff, compression = "lzw")
} else {
  # raggが無い場合は標準tiffデバイス
  ggsave(file.path(outdir, "Figure4_all_windows.tiff"),
         plot = p_all, width = w_mm, height = h_mm, units = "mm",
         dpi = 600, device = "tiff", compression = "lzw")
  ggsave(file.path(outdir, "Figure4_1y_and_≤3y_period.tiff"),
         plot = p_1y_3y, width = w_mm, height = h_mm, units = "mm",
         dpi = 600, device = "tiff", compression = "lzw")
}


#95％信頼区間など
{slope_contrast_by_group <- function(fit, model_label, window_label, level = 0.95){
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
print(df_contrast_out, n = Inf, width = Inf)
}
} # 本解析

#感度分析(recoveryの定義変更)#####
{
library(dplyr)
library(nlme)
library(multcomp)  # glhtを使うので明示
# baseline_cov と covars は既存の主解析と同じ定義を使う前提

# --- 1) ラベル定義：No-dataを明示化 ---
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
    )
  )

# ※CSV往復はしない（UTF-8↔SJIS問題回避）

# --- 2) 時間・CKD整備 ---
akd_time_sens <- jin1_inclusion_sens %>%
  mutate(
    years_from_time0 = as.numeric(date - time0) / 365.25,
    CKD_status = factor(na_if(CKD_status, "nd"), levels = c("nonCKD","CKD"))
  ) %>%
  filter(years_from_time0 >= 0)

# --- 3) longデータ作成：No-data除外 → baseline_cov結合 → 中心化 ---
longdat_sens <- akd_time_sens %>%
  dplyr::select(-any_of(covars)) %>%                 # 既存covarsと衝突回避
  filter(jin_label_sens != "No-data") %>%
  left_join(baseline_cov, by = "id") %>%
  mutate(
    age_c        = scale(age, center = TRUE, scale = FALSE)[,1],
    time0_egfr_c = scale(time0_egfr, center = TRUE, scale = FALSE)[,1],
    jin_label_sens = factor(jin_label_sens, levels = c("nonAKD","Recovery","Non-Recovery"))
  )

# --- 4) 収束安定化：各idで時点>=2＆time分散>0、time中心化 ---
longdat_sens <- longdat_sens %>%
  group_by(id) %>%
  mutate(n_time = n_distinct(date),
         var_t  = var(years_from_time0, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(n_time >= 2, !is.na(var_t), var_t > 0) %>%
  mutate(years_c = as.numeric(scale(years_from_time0, center = TRUE, scale = FALSE)))

# --- 5) モデル当て：まずpdDiag（安定）→必要ならpdSymmへ ---
# ① Base（CKDなし）
fit_base_sens <- lme(
  egfr ~ years_c * jin_label_sens + time0_egfr_c - 1,
  random = list(id = pdDiag(~ 1 + years_c)),
  data = longdat_sens, na.action = na.omit, method = "REML",
  control = lmeControl(opt="optim", msMaxIter=1e5, maxIter=1e5, niterEM=30, tolerance=1e-6)
)

# ② Level-adjusted（CKDなし）
fit_level_adj_sens <- update(
  fit_base_sens,
  . ~ . + age_c + sex + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15
)

# ③ Slope-adjusted（CKDなし）
fit_slope_adj_sens <- update(
  fit_level_adj_sens,
  . ~ . + years_c:(age_c + arb_acei_use +
                     dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 +
                     dn10 + dn12 + dn13 + dn14 + dn15)
)

# --- （任意）CKDを含めた感度解析版 ---
fit_level_adj_sens_ckd <- update(
  fit_base_sens,
  . ~ . + age_c + sex + arb_acei_use + CKD_status +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15
)

fit_slope_adj_sens_ckd <- update(
  fit_level_adj_sens_ckd,
  . ~ . + years_c:(age_c + arb_acei_use + CKD_status +
                     dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 +
                     dn10 + dn12 + dn13 + dn14 + dn15)
)

# ==== 1) ≤1年データを作って years_c を再中心化 ====
longdat_sens_1y <- longdat_sens %>%
  filter(years_from_time0 <= 1) %>%                  # 1年以内に限定
  group_by(id) %>%
  mutate(n_time_1y = n_distinct(date),
         var_t_1y  = var(years_from_time0, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(n_time_1y >= 2, !is.na(var_t_1y), var_t_1y > 0) %>%
  mutate(years_c = as.numeric(scale(years_from_time0, center = TRUE, scale = FALSE)))  # 再中心化

# ==== 2) ≤1年の3モデルを推定（式はALLと同じ） ====
fit_base_sens_1y <- lme(
  egfr ~ years_c * jin_label_sens + time0_egfr_c - 1,
  random = list(id = pdDiag(~ 1 + years_c)),
  data = longdat_sens_1y, na.action = na.omit, method = "REML",
  control = lmeControl(opt="optim", msMaxIter=1e5, maxIter=1e5, niterEM=30, tolerance=1e-6)
)

fit_level_adj_sens_1y <- update(
  fit_base_sens_1y,
  . ~ . + age_c + sex + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15
)

fit_slope_adj_sens_1y <- update(
  fit_level_adj_sens_1y,
  . ~ . + years_c:(age_c + arb_acei_use +
                     dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 +
                     dn10 + dn12 + dn13 + dn14 + dn15)
)

# ==== 3) ALL と ≤1年のスロープ+95%CIを抽出して結合 ====
library(dplyr)
library(nlme)

# 1つの lme オブジェクトから
#  nonAKD / Recovery / Non-Recovery のスロープと95%CIを取り出す関数
slope_ci_by_group_sens <- function(fit, model_label) {
  
  cf <- fixef(fit)        # 固定効果係数
  V  <- vcov(fit)         # 固定効果の分散共分散行列
  
  # 線形結合 L'β から推定値と95%CIを出す小さな関数
  extract_slope <- function(group_label) {
    
    # 係数名ベクトルに合わせた L をゼロで初期化
    L <- rep(0, length(cf))
    names(L) <- names(cf)
    
    ## 基本スロープ（nonAKD 群の years_c）
    if ("years_c" %in% names(cf)) {
      L["years_c"] <- 1
    }
    
    ## 交互作用をグループごとに足す
    if (group_label == "Recovery") {
      nm <- "years_c:jin_label_sensRecovery"
      if (nm %in% names(cf)) L[nm] <- 1
    }
    
    if (group_label == "Non-Recovery") {
      nm <- "years_c:jin_label_sensNon-Recovery"
      if (nm %in% names(cf)) L[nm] <- 1
    }
    
    # 推定値とSE, 95%CI
    est <- sum(L * cf)
    se  <- sqrt( as.numeric(t(L) %*% V %*% L) )
    lower <- est - 1.96 * se
    upper <- est + 1.96 * se
    
    tibble(
      model    = model_label,
      group    = group_label,
      estimate = est,
      lower    = lower,
      upper    = upper
    )
  }
  
  # 3群分を bind
  bind_rows(
    extract_slope("nonAKD"),
    extract_slope("Recovery"),
    extract_slope("Non-Recovery")
  )
}

models_all <- list(
  "Base"           = fit_base_sens,
  "Level-adjusted" = fit_level_adj_sens,
  "Slope-adjusted" = fit_slope_adj_sens
)

models_1y <- list(
  "Base"           = fit_base_sens_1y,
  "Level-adjusted" = fit_level_adj_sens_1y,
  "Slope-adjusted" = fit_slope_adj_sens_1y
)

df_all <- purrr::imap_dfr(models_all, ~ slope_ci_by_group_sens(.x, .y)) %>%
  dplyr::mutate(window = "All period")

df_1y  <- purrr::imap_dfr(models_1y,  ~ slope_ci_by_group_sens(.x, .y)) %>%
  dplyr::mutate(window = "\u22641 year")  # "≤1 year"

df_slope_2win <- dplyr::bind_rows(df_1y, df_all) %>%
  dplyr::mutate(
    model  = forcats::fct_relevel(model, c("Base","Level-adjusted","Slope-adjusted")),
    group  = forcats::fct_relevel(group, c("nonAKD","Recovery","Non-Recovery")),
    window = factor(window, levels = c("\u22641 year","All period"))
  )

library(dplyr)
library(multcomp)
library(ggplot2)

get_slopes_3group <- function(fit, window_label,
                              group_levels = c("nonAKD", "Recovery", "Non-Recovery")) {
  
  cf <- names(fixef(fit))
  time_term <- if ("years_c" %in% cf) "years_c" else "years_from_time0"
  if (!time_term %in% cf) stop("時間項（years_c / years_from_time0）が見つかりません。")
  
  int_recovery <- paste0(time_term, ":jin_labelRecovery")
  int_nonrec   <- paste0(time_term, ":jin_labelNon-Recovery")
  
  v <- function(group){
    L <- rep(0, length(cf)); names(L) <- cf
    L[time_term] <- 1
    if (group == "Recovery" && int_recovery %in% cf) L[int_recovery] <- 1
    if (group == "Non-Recovery" && int_nonrec   %in% cf) L[int_nonrec]   <- 1
    L
  }
  
  K <- rbind(
    "nonAKD"       = v("nonAKD"),
    "Recovery"     = v("Recovery"),
    "Non-Recovery" = v("Non-Recovery")
  )
  
  g  <- glht(fit, linfct = K)
  ci <- confint(g)$confint
  
  tibble(
    window   = window_label,
    model    = "Slope-adjusted",
    group    = rownames(ci),
    estimate = ci[, "Estimate"],
    lower    = ci[, "lwr"],
    upper    = ci[, "upr"]
  ) %>%
    mutate(
      window = factor(window, levels = c("≤1 year", "≤3 years")),
      group  = factor(group,  levels = group_levels)
    )
}

# ---- ここで作ったdfを「上書きしない」 ----
df_slope_1y3y_slopeadj <- bind_rows(
  get_slopes_3group(fit_slope_adj_1y, "≤1 year"),
  get_slopes_3group(fit_slope_adj_3y, "≤3 years")
)

# 確認（ここが両方3ずつ出ればOK）
table(df_slope_1y3y_slopeadj$window)

p_1y3y_slopeadj <- ggplot(df_slope_1y3y_slopeadj, aes(x = group, y = estimate, fill = group)) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
  geom_col(width = 0.65) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 0.5) +
  facet_grid(. ~ window)+
  labs(
    title = "Estimated eGFR slopes by group\n(≤1 year vs ≤3 years; Slope-adjusted model; Sensitivity Analysis)",
    x = NULL,
    y = expression(paste("Slope (mL/min/1.73 m"^2," per year), 95% CI")),
    fill = "Group"
  ) +
  scale_fill_manual(values = c(
    nonAKD = "#E41A1C", Recovery = "#4DAF4A", `Non-Recovery` = "#377EB8"
  )) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor  = element_blank(),
    strip.background  = element_rect(fill = "grey95", colour = NA),
    strip.text        = element_text(size = 10, face = "bold", lineheight = 1.05),
    legend.position   = "bottom",
    legend.key.height = unit(0.5, "lines"),
    legend.key.width  = unit(1.2, "lines"),
    axis.title.y      = element_text(margin = margin(r = 6)),
    plot.title        = element_text(hjust = 0, face = "bold", lineheight = 1.1, margin = margin(b = 8)),
    plot.margin       = margin(t = 6, r = 8, b = 6, l = 6),
    panel.spacing.y   = unit(1.2, "lines"),
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank()
  )

# 保存
setwd("E:/R")
w_in <- 180/25.4
h_in <- 120/25.4

ggsave("Figure5b_eGFR_slopes_1y_3y_slopeAdjusted.pdf",
       plot = p_1y3y_slopeadj, device = cairo_pdf, width = w_in, height = h_in, units = "in")

ggsave("Figure5b_eGFR_slopes_1y_3y_slopeAdjusted_600dpi.tiff",
       plot = p_1y3y_slopeadj, device = "tiff", compression = "lzw", dpi = 600,
       width = w_in, height = h_in, units = "in", bg = "white")

}

library(dplyr)
library(purrr)
library(forcats)
library(multcomp)

# --- nonAKD基準の群間差（Recovery−nonAKD, Non-Recovery−nonAKD）を返す関数 ---
contrast_slope_diff <- function(fit, model_label, window_label){
  cf <- names(fixef(fit))
  time_term <- if ("years_c" %in% cf) "years_c" else "years_from_time0"
  term_rec  <- cf[grepl(paste0("^", time_term, ":.*Recovery$"),       cf)]
  term_non  <- cf[grepl(paste0("^", time_term, ":.*Non-Recovery$"),   cf)]
  
  L <- matrix(0, nrow = 2, ncol = length(cf),
              dimnames = list(c("Recovery − nonAKD","Non-Recovery − nonAKD"), cf))
  if (length(term_rec) > 0) L["Recovery − nonAKD", term_rec] <- 1
  if (length(term_non) > 0) L["Non-Recovery − nonAKD", term_non] <- 1
  
  g  <- multcomp::glht(fit, linfct = L)
  s  <- summary(g)
  ci <- confint(g)
  
  tibble(
    contrast = rownames(L),
    diff     = ci$confint[, "Estimate"],
    lower    = ci$confint[, "lwr"],
    upper    = ci$confint[, "upr"],
    p_value  = s$test$pvalues,
    model    = model_label,
    window   = window_label
  ) %>%
    mutate(
      signif = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ "ns"
      )
    )
}

# --- 全期間 / ≤1年 のデータセットを作成（図は作らない） ---
df_contrast <- bind_rows(
  imap_dfr(models_all, ~ contrast_slope_diff(.x, .y, "All period")),
  imap_dfr(models_1y,  ~ contrast_slope_diff(.x, .y, "\u22641 year"))  # ≤1 year
) %>%
  mutate(
    model  = fct_relevel(model, c("Base","Level-adjusted","Slope-adjusted")),
    window = factor(window, levels = c("\u22641 year","All period"))
  ) %>%
  arrange(model, window, match(contrast, c("Recovery − nonAKD","Non-Recovery − nonAKD")))

# --- 表示（丸めたい場合は下の1行を有効化） ---
# df_contrast <- df_contrast %>% mutate(across(c(diff, lower, upper, p_value), ~round(.x, 3)))

df_contrast









#egfr抜き####
{
############################################################
## AKD eGFR slope models (id_under60 を除外した版)
##  - jin1_Eligibile, jin1_inclusion を読み込み
##  - id_under60.csv の id を除外
##  - 混合効果モデル（Base / Level / Slope）
##  - 1年, 2年, 3年, 全期間のスロープ推定 & 棒グラフ
############################################################

## --- パッケージ読み込み -----------------------------------
library(dplyr)
library(readr)
library(nlme)
library(tidyr)
library(multcomp)
library(ggplot2)
library(stringr)
library(purrr)

## --- 作業ディレクトリ -------------------------------------
setwd("E:/R")

## --- 元データ読み込み --------------------------------------
jin1_Eligibile <- read_csv("jin1_Eligibile.csv",
                           locale = locale(encoding = "SHIFT-JIS"))
jin1_inclusion <- read_csv("jin1_inclusion.csv",
                           locale = locale(encoding = "SHIFT-JIS"))

## --- 除外ID（id_under60）読み込み --------------------------
id_under60_tbl <- read_csv("id_under60.csv")  # 少なくとも id 列がある前提

## --- id_under60 に該当する id を除いた版データ -------------
jin1_Eligibile_idrm <- jin1_Eligibile %>%
  filter(!id %in% id_under60_tbl$id)

jin1_inclusion_idrm <- jin1_inclusion %>%
  filter(!id %in% id_under60_tbl$id)

## --- akd_time_m: years_from_time0 >= 0, CKD_status 整形 ----
akd_time_m_idrm <- jin1_inclusion_idrm %>%
  mutate(
    years_from_time0 = as.numeric(date - time0) / 365.25
  ) %>%
  filter(years_from_time0 >= 0) %>%
  mutate(
    CKD_status = na_if(CKD_status, "nd"),
    CKD_status = factor(CKD_status, levels = c("nonCKD", "CKD"))
  )

## --- ベースライン共変量（baseline_cov_idrm）---------------
## index_cre は使わず、time0_egfr と合わせる方針
baseline_cov_idrm <- jin1_Eligibile_idrm %>%
  filter(exclude == "include") %>%
  group_by(id) %>%
  arrange(index_date) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    arb_acei_use = if_else(arb == 1 | acei == 1, 1L, 0L)
  ) %>%
  dplyr::select(
    id, age, sex, arb_acei_use,
    dn1, dn3, dn4, dn5, dn6, dn7, dn8, dn9,
    dn10, dn12, dn13, dn14, dn15
  )

## --- longdat_idrm を作成（解析用ロングデータ）--------------
covars <- c(
  "age","sex","arb_acei_use",
  "dn1","dn3","dn4","dn5","dn6","dn7","dn8","dn9",
  "dn10","dn12","dn13","dn14","dn15"
)

longdat_idrm <- akd_time_m_idrm %>%
  dplyr::select(-any_of(covars)) %>%   # 重複候補を削除
  left_join(baseline_cov_idrm, by = "id") %>%
  mutate(
    age_c        = scale(age, center = TRUE, scale = FALSE)[,1],
    time0_egfr_c = scale(time0_egfr, center = TRUE, scale = FALSE)[,1],
    jin_label    = factor(jin_label,
                          levels = c("nonAKD", "Recovery", "Non-Recovery"))
  )

############################################################
## 混合効果モデルの推定（全期間）
############################################################

## ① ベースモデル（time0_egfr のみ）
fit_base_idrm <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1,
  random = list(id = pdSymm(~ 1 + years_from_time0)),
  na.action = na.omit,
  data = longdat_idrm,
  method = "REML",
  control = lmeControl(
    maxIter = 1e8, msMaxIter = 1e8,
    opt = "optim", optimMethod = "L-BFGS-B"
  )
)

## ② レベル調整モデル（Cox と同じ共変量を主効果で追加）
fit_lvl_idrm <- update(
  fit_base_idrm,
  . ~ . + age_c + sex + arb_acei_use + CKD_status +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 +
    dn10 + dn12 + dn13 + dn14 + dn15
)

## ③ スロープ調整モデル（years_from_time0 との交互作用を追加）
fit_slp_idrm <- update(
  fit_lvl_idrm,
  . ~ . + years_from_time0:(age_c + arb_acei_use + CKD_status +
                              dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 +
                              dn10 + dn12 + dn13 + dn14 + dn15)
)

############################################################
## 時間窓ごとのモデル（≤1年, ≤2年, ≤3年）
############################################################

## ---- ≤1 year ----
fit_base_1y_idrm <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1,
  random = list(id = pdSymm(~ 1 + years_from_time0)),
  na.action = na.omit,
  data = dplyr::filter(longdat_idrm, years_from_time0 <= 1),
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)

fit_lvl_1y_idrm <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1 +
    age_c + sex + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 +
    dn10 + dn12 + dn13 + dn14 + dn15,
  random = list(id = pdSymm(~ 1 + years_from_time0)),
  na.action = na.omit,
  data = dplyr::filter(longdat_idrm, years_from_time0 <= 1),
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)

fit_slp_1y_idrm <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1 +
    age_c + sex + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 +
    dn10 + dn12 + dn13 + dn14 + dn15 +
    years_from_time0:(age_c + arb_acei_use +
                        dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 +
                        dn10 + dn12 + dn13 + dn14 + dn15),
  random = list(id = pdSymm(~ 1 + years_from_time0)),
  na.action = na.omit,
  data = dplyr::filter(longdat_idrm, years_from_time0 <= 1),
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)

## ---- ≤2 years ----
fit_base_2y_idrm <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1,
  random = list(id = pdSymm(~ 1 + years_from_time0)),
  na.action = na.omit,
  data = dplyr::filter(longdat_idrm, years_from_time0 <= 2),
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)

fit_lvl_2y_idrm <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1 +
    age_c + sex + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 +
    dn10 + dn12 + dn13 + dn14 + dn15,
  random = list(id = pdSymm(~ 1 + years_from_time0)),
  na.action = na.omit,
  data = dplyr::filter(longdat_idrm, years_from_time0 <= 2),
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)

fit_slp_2y_idrm <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1 +
    age_c + sex + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 +
    dn10 + dn12 + dn13 + dn14 + dn15 +
    years_from_time0:(age_c + arb_acei_use +
                        dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 +
                        dn10 + dn12 + dn13 + dn14 + dn15),
  random = list(id = pdSymm(~ 1 + years_from_time0)),
  na.action = na.omit,
  data = dplyr::filter(longdat_idrm, years_from_time0 <= 2),
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)

## ---- ≤3 years ----
fit_base_3y_idrm <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1,
  random = list(id = pdSymm(~ 1 + years_from_time0)),
  na.action = na.omit,
  data = dplyr::filter(longdat_idrm, years_from_time0 <= 3),
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)

fit_lvl_3y_idrm <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1 +
    age_c + sex + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 +
    dn10 + dn12 + dn13 + dn14 + dn15,
  random = list(id = pdSymm(~ 1 + years_from_time0)),
  na.action = na.omit,
  data = dplyr::filter(longdat_idrm, years_from_time0 <= 3),
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)

fit_slp_3y_idrm <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1 +
    age_c + sex + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 +
    dn10 + dn12 + dn13 + dn14 + dn15 +
    years_from_time0:(age_c + arb_acei_use +
                        dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 +
                        dn10 + dn12 + dn13 + dn14 + dn15),
  random = list(id = pdSymm(~ 1 + years_from_time0)),
  na.action = na.omit,
  data = dplyr::filter(longdat_idrm, years_from_time0 <= 3),
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)

############################################################
## 群ごとのスロープ推定 + 95%CI を取り出す関数
############################################################

slope_ci_by_group <- function(fit, model_label, window_label){
  cf <- names(fixef(fit))
  v <- function(g){
    vec <- rep(0, length(cf)); names(vec) <- cf
    # 基本スロープ
    if ("years_from_time0" %in% cf) vec["years_from_time0"] <- 1
    # 群×時間の交互作用（あれば加算）
    if (g == "Recovery" &&
        "years_from_time0:jin_labelRecovery" %in% cf) {
      vec["years_from_time0:jin_labelRecovery"] <- 1
    }
    if (g == "Non-Recovery" &&
        "years_from_time0:jin_labelNon-Recovery" %in% cf) {
      vec["years_from_time0:jin_labelNon-Recovery"] <- 1
    }
    vec
  }
  L <- rbind(
    nonAKD        = v("nonAKD"),
    Recovery      = v("Recovery"),
    `Non-Recovery`= v("Non-Recovery")
  )
  ci <- suppressMessages(confint(glht(fit, linfct = L)))
  tibble(
    group    = rownames(L),
    estimate = ci$confint[,"Estimate"],
    lower    = ci$confint[,"lwr"],
    upper    = ci$confint[,"upr"],
    model    = model_label,
    window   = window_label
  )
}

############################################################
## 各時間窓 × 各モデルからスロープを抽出
############################################################

fits_idrm <- list(
  "≤1 year" = list(
    "Base"            = fit_base_1y_idrm,
    "Level-adjusted"  = fit_lvl_1y_idrm,
    "Slope-adjusted"  = fit_slp_1y_idrm
  ),
  "≤2 years" = list(
    "Base"            = fit_base_2y_idrm,
    "Level-adjusted"  = fit_lvl_2y_idrm,
    "Slope-adjusted"  = fit_slp_2y_idrm
  ),
  "≤3 years" = list(
    "Base"            = fit_base_3y_idrm,
    "Level-adjusted"  = fit_lvl_3y_idrm,
    "Slope-adjusted"  = fit_slp_3y_idrm
  ),
  "All period" = list(
    "Base"            = fit_base_idrm,
    "Level-adjusted"  = fit_lvl_idrm,
    "Slope-adjusted"  = fit_slp_idrm
  )
)

df_plot_idrm <- imap_dfr(fits_idrm, function(models, win){
  imap_dfr(models, function(fit, mdl){
    slope_ci_by_group(fit, model_label = mdl, window_label = win)
  })
})

## ラベル正規化と因子化
df_plot_idrm <- df_plot_idrm %>%
  mutate(
    model = case_when(
      model %in% c("Base","① Base")                         ~ "Base",
      model %in% c("Level-adjusted","② Level-adjusted")     ~ "Level-adjusted",
      model %in% c("Slope-adjusted","③ Slope-adjusted")     ~ "Slope-adjusted",
      TRUE ~ NA_character_
    ),
    window = case_when(
      window %in% c("≤1 year","<=1 year","≤1year","<=1year")       ~ "≤1 year",
      window %in% c("≤2 years","<=2 years","≤2years","<=2years")   ~ "≤2 years",
      window %in% c("≤3 years","<=3 years","≤3years","<=3years")   ~ "≤3 years",
      window %in% c("All period","All-period","All Period")        ~ "All period",
      TRUE ~ NA_character_
    ),
    group = trimws(gsub("\u2013|\u2212", "-", as.character(group)))
  ) %>%
  filter(
    !is.na(model),
    !is.na(window),
    group %in% c("nonAKD","Recovery","Non-Recovery")
  ) %>%
  mutate(
    model  = factor(model,
                    levels = c("Base","Level-adjusted","Slope-adjusted")),
    window = factor(window,
                    levels = c("≤1 year","≤2 years","≤3 years","All period")),
    group  = factor(group,
                    levels = c("nonAKD","Recovery","Non-Recovery"))
  ) %>%
  distinct(window, model, group, .keep_all = TRUE)

############################################################
## 作図：全期間＋各時間窓
############################################################

## 全時間窓（1・2・3年 + 全期間）
p_all_idrm <-
  ggplot(df_plot_idrm,
         aes(x = group, y = estimate, fill = group)) +
  geom_col(width = 0.65) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.2, linewidth = 0.5) +
  facet_grid(model ~ window, drop = FALSE) +
  labs(
    title = "Estimated eGFR Slopes\nby Group across Models and Time Windows (id-excluded)",
    x = NULL,
    y = expression(paste("Slope (mL/min/1.73 m"^2," per year), 95% CI")),
    fill = "Group"
  ) +
  scale_fill_manual(values = c(
    nonAKD       = "#E41A1C",
    Recovery     = "#4DAF4A",
    `Non-Recovery` = "#377EB8"
  )) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor  = element_blank(),
    strip.background  = element_rect(fill = "grey95", colour = NA),
    legend.position   = "bottom",
    axis.text.x       = element_blank()
  )

## ≤1年 + 全期間のみ
p_1y_all_idrm <-
  ggplot(df_plot_idrm %>% filter(window %in% c("≤1 year","All period")),
         aes(x = group, y = estimate, fill = group)) +
  geom_col(width = 0.65) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.2, linewidth = 0.5) +
  facet_grid(
    model ~ window,
    drop   = TRUE,
    switch = "y"
  ) +
  labs(
    title = "Estimated eGFR Slopes (≤1 year and All period, id-excluded)",
    x = NULL,
    y = expression(paste("Slope (mL/min/1.73 m"^2," per year), 95% CI")),
    fill = "Group"
  ) +
  scale_fill_manual(values = c(
    nonAKD       = "#E41A1C",
    Recovery     = "#4DAF4A",
    `Non-Recovery` = "#377EB8"
  )) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor   = element_blank(),
    strip.background   = element_rect(fill = "grey95", colour = NA),
    strip.placement    = "outside",
    legend.position    = "bottom",
    axis.text.x        = element_blank()
  )

############################################################
## 図の保存（PDF & TIFF）: *_idrm で出力
############################################################

library(Cairo)

w_mm <- 180
h_mm <- 120
outdir <- "E:/R"

## PDF
ggsave(file.path(outdir, "Figure4_all_windows_idrm.pdf"),
       plot = p_all_idrm, width = w_mm, height = h_mm,
       units = "mm", device = "pdf", dpi = 300)

ggsave(file.path(outdir, "Figure4_1y_and_allperiod_idrm.pdf"),
       plot = p_1y_all_idrm, width = w_mm, height = h_mm,
       units = "mm", device = "pdf", dpi = 300)

## TIFF（600 dpi, ragg があれば優先）
use_ragg <- requireNamespace("ragg", quietly = TRUE)

if (use_ragg) {
  ggsave(file.path(outdir, "Figure4_all_windows_idrm.tiff"),
         plot = p_all_idrm, width = w_mm, height = h_mm,
         units = "mm", dpi = 600,
         device = ragg::agg_tiff, compression = "lzw")
  ggsave(file.path(outdir, "Figure4_1y_and_allperiod_idrm.tiff"),
         plot = p_1y_all_idrm, width = w_mm, height = h_mm,
         units = "mm", dpi = 600,
         device = ragg::agg_tiff, compression = "lzw")
} else {
  ggsave(file.path(outdir, "Figure4_all_windows_idrm.tiff"),
         plot = p_all_idrm, width = w_mm, height = h_mm,
         units = "mm", dpi = 600,
         device = "tiff", compression = "lzw")
  ggsave(file.path(outdir, "Figure4_1y_and_allperiod_idrm.tiff"),
         plot = p_1y_all_idrm, width = w_mm, height = h_mm,
         units = "mm", dpi = 600,
         device = "tiff", compression = "lzw")
}

############################################################
## ここまでで：
## ・id_under60 の id を除外
## ・全期間＋時間窓ごと混合効果モデル
## ・Figure4 相当のスロープ図（id 除外版）を保存
############################################################

}
#感度分析(CKDの要素入れ)#####
{
# 必要パッケージ
library(dplyr)
library(purrr)
library(multcomp)
library(ggplot2)
#感度解析①（CKD の主効果を追加）
fit_level_adj_ckd <- update(fit_base,
                            . ~ . + age_c + sex + arb_acei_use + CKD_status +
                              dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15
)
#感度解析②（CKD の時間交互作用も追加：スロープ差への寄与を評価）
fit_slope_adj_ckd <- update(fit_level_adj_ckd,
                            . ~ . + years_from_time0:CKD_status +
                              years_from_time0:(age_c + arb_acei_use +
                                                  dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15)
)

# ---- 汎用：群別スロープ(固定効果の線形和) + 95%CI を返す ----
slope_ci_by_group_any <- function(fit, model_label) {
  cf <- names(fixef(fit))
  
  # timeの項名を自動検出（years_c が優先、なければ years_from_time0）
  time_term <- if ("years_c" %in% cf) "years_c" else "years_from_time0"
  
  # Recovery / Non-Recovery の相互作用係数名を自動検出（_sens でも無印でもOK）
  term_rec <- cf[grepl(paste0("^", time_term, ":.*Recovery$"), cf)]
  term_non <- cf[grepl(paste0("^", time_term, ":.*Non-Recovery$"), cf)]
  
  mk_vec <- function(g){
    v <- setNames(rep(0, length(cf)), cf)
    if (time_term %in% cf) v[time_term] <- 1
    if (g == "Recovery"     && length(term_rec) == 1) v[term_rec] <- 1
    if (g == "Non-Recovery" && length(term_non) == 1) v[term_non] <- 1
    v
  }
  
  L <- rbind(
    nonAKD        = mk_vec("nonAKD"),
    Recovery      = mk_vec("Recovery"),
    `Non-Recovery`= mk_vec("Non-Recovery")
  )
  
  ci <- suppressMessages(confint(multcomp::glht(fit, linfct = L)))
  tibble(
    group    = rownames(L),
    estimate = ci$confint[, "Estimate"],
    lower    = ci$confint[, "lwr"],
    upper    = ci$confint[, "upr"],
    model    = model_label
  )
}
# ==== 0) データセットの参照（fit_baseに使ったlongデータを想定） ====
dat_all <- if (exists("longdat")) longdat else longdat_sens

# ==== 1) ≤1年データ作成（個体内に時点2つ以上・時間分散>0を担保） ====
dat_1y <- dat_all %>%
  dplyr::filter(years_from_time0 <= 1) %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(n_time_1y = dplyr::n_distinct(date),
                var_t_1y  = stats::var(years_from_time0, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n_time_1y >= 2, !is.na(var_t_1y), var_t_1y > 0)

# ==== 2) ≤1年のCKD感度モデルを再推定（式はALLと同じ、dataだけ差し替え） ====
fit_level_adj_ckd_1y <- update(
  fit_base,
  . ~ . + age_c + sex + arb_acei_use + CKD_status +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
  data = dat_1y
)

fit_slope_adj_ckd_1y <- update(
  fit_level_adj_ckd_1y,
  . ~ . + years_from_time0:CKD_status +
    years_from_time0:(age_c + arb_acei_use +
                        dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15),
  data = dat_1y
)

# ==== 3) スロープ + 95%CI 抽出 → ALL / ≤1年 を結合 ====
models_all_ckd <- list(
  "Level-adjusted + CKD" = fit_level_adj_ckd,
  "Slope-adjusted + CKD" = fit_slope_adj_ckd
)
models_1y_ckd <- list(
  "Level-adjusted + CKD" = fit_level_adj_ckd_1y,
  "Slope-adjusted + CKD" = fit_slope_adj_ckd_1y
)

df_all_ckd <- purrr::imap_dfr(models_all_ckd, ~ slope_ci_by_group_any(.x, .y)) %>%
  dplyr::mutate(window = "All period")
df_1y_ckd  <- purrr::imap_dfr(models_1y_ckd,  ~ slope_ci_by_group_any(.x, .y)) %>%
  dplyr::mutate(window = "\u22641 year")  # ≤1 year

df_ckd_2win <- dplyr::bind_rows(df_1y_ckd, df_all_ckd) %>%
  dplyr::mutate(
    model  = forcats::fct_relevel(model, c("Level-adjusted + CKD","Slope-adjusted + CKD")),
    group  = forcats::fct_relevel(group, c("nonAKD","Recovery","Non-Recovery")),
    window = factor(window, levels = c("\u22641 year","All period"))
  )

# ==== 4) 図：行=window（≤1年/全期間） × 列=model（+CKDの2モデル） ====
{p_ckd_2win <- ggplot(df_ckd_2win, aes(x = group, y = estimate, fill = group)) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
  geom_col(width = 0.65) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 0.5) +
  facet_grid(window ~ model) +
  labs(
    title = "Estimated eGFR slopes by group\n(\u22641 year vs All period; +CKD sensitivity)",
    x = NULL,
    y = expression(paste("Slope (mL/min/1.73 m"^2," per year), 95% CI")),
    fill = "Group"
  ) +
  scale_fill_manual(values = c(
    nonAKD = "#E41A1C", Recovery = "#4DAF4A", `Non-Recovery` = "#377EB8"
  )) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor  = element_blank(),
    strip.background  = element_rect(fill = "grey95", colour = NA),
    strip.text        = element_text(size = 10, face = "bold", lineheight = 1.05),
    legend.position   = "bottom",
    legend.key.height = unit(0.5, "lines"),
    legend.key.width  = unit(1.2, "lines"),
    axis.title.y      = element_text(margin = margin(r = 6)),
    plot.title        = element_text(hjust = 0, face = "bold", lineheight = 1.1, margin = margin(b = 8)),
    plot.margin       = margin(t = 6, r = 8, b = 6, l = 6),
    panel.spacing.x   = unit(1.0, "lines"),
    panel.spacing.y   = unit(1.2, "lines"),
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank()
  )

# ==== 5) 保存（2段×2列に合わせたサイズ） ====
w_in <- 160/25.4  # 6.30 in
h_in <- 150/25.4  # 5.90 in
ggsave("Figure5b_eGFR_slopes_CKD_2windows.pdf",
       plot = p_ckd_2win, device = cairo_pdf, width = w_in, height = h_in, units = "in")
ggsave("Figure5b_eGFR_slopes_CKD_2windows_600dpi.tiff",
       plot = p_ckd_2win, device = "tiff", compression = "lzw", dpi = 600,
       width = w_in, height = h_in, units = "in", bg = "white")


# ==== 本解析とCKDの結果並べる図 ====
library(dplyr)
library(purrr)
library(forcats)
library(ggplot2)
library(multcomp)
library(patchwork)  # install.packages("patchwork")
library(grid)
library(Cairo)

# 1) 傾き+95%CI（years_c / years_from_time0 どちらでもOK）
slope_ci_by_group_any <- function(fit, model_label) {
  cf <- names(fixef(fit))
  time_term <- if ("years_c" %in% cf) "years_c" else "years_from_time0"
  term_rec  <- cf[grepl(paste0("^", time_term, ":.*Recovery$"), cf)]
  term_non  <- cf[grepl(paste0("^", time_term, ":.*Non-Recovery$"), cf)]
  mk_vec <- function(g){
    v <- setNames(rep(0, length(cf)), cf)
    if (time_term %in% cf) v[time_term] <- 1
    if (g == "Recovery"     && length(term_rec) == 1) v[term_rec] <- 1
    if (g == "Non-Recovery" && length(term_non) == 1) v[term_non] <- 1
    v
  }
  L  <- rbind(nonAKD = mk_vec("nonAKD"),
              Recovery = mk_vec("Recovery"),
              `Non-Recovery` = mk_vec("Non-Recovery"))
  ci <- suppressMessages(confint(multcomp::glht(fit, linfct = L)))
  tibble(
    group    = rownames(L),
    estimate = ci$confint[, "Estimate"],
    lower    = ci$confint[, "lwr"],
    upper    = ci$confint[, "upr"],
    model    = model_label
  )
}

# 2) 図作成（行=window, 列=model、X軸ラベル非表示）
make_slope_plot <- function(df, title_text) {
  ggplot(df, aes(x = group, y = estimate, fill = group)) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
    geom_col(width = 0.65) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 0.5) +
    facet_grid(window ~ model) +
    labs(
      title = title_text, x = NULL,
      y = expression(paste("Slope (mL/min/1.73 m"^2," per year), 95% CI")),
      fill = "Group"
    ) +
    scale_fill_manual(values = c(nonAKD="#E41A1C", Recovery="#4DAF4A", `Non-Recovery`="#377EB8")) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text       = element_text(size = 10, face = "bold", lineheight = 1.05),
      legend.position  = "bottom",
      legend.key.height= unit(0.5,"lines"),
      legend.key.width = unit(1.2,"lines"),
      axis.title.y     = element_text(margin = margin(r = 6)),
      plot.title       = element_text(hjust = 0, face = "bold", margin = margin(b = 8)),
      panel.spacing.x  = unit(1.0, "lines"),
      panel.spacing.y  = unit(1.2, "lines"),
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank()
    )
}

# 3) ALL/≤1年の2窓データを一括生成（モデルリスト→df）
build_df_2windows <- function(models_all, models_1y, model_levels, title_models){
  df_all <- imap_dfr(models_all, ~ slope_ci_by_group_any(.x, .y)) %>%
    mutate(window = "All period")
  df_1y  <- imap_dfr(models_1y,  ~ slope_ci_by_group_any(.x, .y)) %>%
    mutate(window = "\u22641 year")  # ≤1 year
  
  bind_rows(df_1y, df_all) %>%
    mutate(
      model  = fct_relevel(model, model_levels),
      group  = fct_relevel(group, c("nonAKD","Recovery","Non-Recovery")),
      window = factor(window, levels = c("\u22641 year","All period")),
      .title = title_models
    )
}

# ===== ① 本解析（CKDを“入れない”想定の3モデル） ============================
# 既存オブジェクト名に合わせる：fit_base / fit_level_adj / fit_slope_adj （全期間）
# と、fit_base_1y / fit_level_adj_1y / fit_slope_adj_1y （≤1年）
models_main_all <- list(
  "Base"           = fit_base,
  "Level-adjusted" = fit_level_adj,
  "Slope-adjusted" = fit_slope_adj
)
models_main_1y <- list(
  "Base"           = fit_base_1y,
  "Level-adjusted" = fit_level_adj_1y,
  "Slope-adjusted" = fit_slope_adj_1y
)
df_main <- build_df_2windows(
  models_all   = models_main_all,
  models_1y    = models_main_1y,
  model_levels = c("Base","Level-adjusted","Slope-adjusted"),
  title_models = "Base / Level-adjusted / Slope-adjusted"
)
p_main <- make_slope_plot(
  df_main,
  "Estimated eGFR slopes by group (≤1 year vs All period)\n— Primary analysis"
)

# ===== ② 感度分析（CKDを入れる2モデル） =====================================
# 既に fit_level_adj_ckd / fit_slope_adj_ckd と _1y があれば使用。
# 無ければ mainモデルから update で作成（安全ネット）。
if (!exists("fit_level_adj_ckd")) {
  fit_level_adj_ckd <- update(
    fit_base, . ~ . + age_c + sex + arb_acei_use + CKD_status +
      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15
  )
}
if (!exists("fit_slope_adj_ckd")) {
  fit_slope_adj_ckd <- update(
    fit_level_adj_ckd, . ~ . + years_from_time0:(age_c + arb_acei_use + CKD_status +
                                                   dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15)
  )
}
if (!exists("fit_level_adj_ckd_1y")) {
  fit_level_adj_ckd_1y <- update(
    fit_base_1y, . ~ . + age_c + sex + arb_acei_use + CKD_status +
      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15
  )
}
if (!exists("fit_slope_adj_ckd_1y")) {
  fit_slope_adj_ckd_1y <- update(
    fit_level_adj_ckd_1y, . ~ . + years_from_time0:(age_c + arb_acei_use + CKD_status +
                                                      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15)
  )
}

models_ckd_all <- list(
  "Level-adjusted + CKD" = fit_level_adj_ckd,
  "Slope-adjusted + CKD" = fit_slope_adj_ckd
)
models_ckd_1y <- list(
  "Level-adjusted + CKD" = fit_level_adj_ckd_1y,
  "Slope-adjusted + CKD" = fit_slope_adj_ckd_1y
)
df_ckd <- build_df_2windows(
  models_all   = models_ckd_all,
  models_1y    = models_ckd_1y,
  model_levels = c("Level-adjusted + CKD","Slope-adjusted + CKD"),
  title_models = "+CKD sensitivity"
)
p_ckd <- make_slope_plot(
  df_ckd,
  "Estimated eGFR slopes by group (≤1 year vs All period)\n— Sensitivity (+CKD)"
)

# ===== ③ 2枚を並べて1枚に（A/Bパネル） =====================================
combined <- p_main / p_ckd + plot_annotation(tag_levels = 'A')
# ===== ④ 保存（2×?列のため少し縦長） ========================================
w_in <- 190/25.4  # 7.48 in
h_in <- 300/25.4  # 11.8 in
ggsave("Figure5_main_vs_CKD_2windows_combined.pdf",
       plot = combined, device = cairo_pdf, width = w_in, height = h_in, units = "in")
ggsave("Figure5_main_vs_CKD_2windows_combined_600dpi.tiff",
       plot = combined, device = "tiff", compression = "lzw", dpi = 600,
       width = w_in, height = h_in, units = "in", bg = "white")







#recovery判定変更とCKDをすべて入れる
{
  # ---- 1) モデルをリスト化（主＋CKD版を並べる）----
  # 既に fit_base_sens, fit_level_adj_sens, fit_slope_adj_sens,
  #       fit_level_adj_ckd, fit_slope_adj_ckd が存在する前提
  models_all <- list(
    "Base"                       = fit_base_sens,
    "Level-adjusted"             = fit_level_adj_sens,
    "Slope-adjusted"             = fit_slope_adj_sens,
    "Level-adjusted + CKD"       = fit_level_adj_ckd,
    "Slope-adjusted + CKD"       = fit_slope_adj_ckd
  )
  
  df_slope <- imap_dfr(models_all, ~ slope_ci_by_group_any(.x, .y)) %>%
    mutate(
      model = factor(model, levels = names(models_all)),               # 並び固定
      group = factor(group, levels = c("nonAKD","Recovery","Non-Recovery"))
    )
  # ---- 2) 図：facet_wrapで2段構成（上3・下2） ----
  p <- ggplot(df_slope, aes(x = group, y = estimate, fill = group)) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
    geom_col(width = 0.65) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 0.5) +
    facet_wrap(~ model, nrow = 2) +                                   # ★ここだけ変更
    labs(
      title = "Estimated eGFR slopes by group \n(main vs +CKD sensitivity)",  # 誤字修正
      x = NULL,
      y = expression(paste("Slope (mL/min/1.73 m"^2," per year), 95% CI")),
      fill = "Group"
    ) +
    scale_fill_manual(values = c(
      nonAKD = "#E41A1C", Recovery = "#4DAF4A", `Non-Recovery` = "#377EB8"
    )) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor  = element_blank(),
      strip.background  = element_rect(fill = "grey95", colour = NA),
      strip.text        = element_text(size = 10, face = "bold"),
      legend.position   = "bottom",
      legend.key.height = unit(0.5, "lines"),
      legend.key.width  = unit(1.2, "lines"),
      axis.title.y      = element_text(margin = margin(r = 6)),
      plot.title        = element_text(hjust = 0, face = "bold", lineheight = 1.1, margin = margin(b = 8)),
      plot.margin       = margin(t = 6, r = 8, b = 6, l = 6),
      panel.spacing.x   = unit(1.0, "lines"),                         # ★横余白を少し広げる
      panel.spacing.y   = unit(1.2, "lines")                          # ★段間の余白
    ) +
    # X軸ラベルを消すならこの2行（不要なら削除）
    scale_x_discrete(labels = NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  # 必要パッケージ
  install.packages(c("ggplot2","Cairo"))  # 未インストールなら
  library(ggplot2)
  install.packages("Cairo")
  library(Cairo)
  
  # ----- 仕上げ：図オブジェクトを作成（あなたのコードを微調整）-----
  # 再現性のための因子順（任意）
  df_slope$group <- factor(df_slope$group, levels = c("nonAKD","Recovery","Non-Recovery"))
  
  p <- ggplot(df_slope, aes(x = group, y = estimate, fill = group)) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
    geom_col(width = 0.65) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 0.5) +
    facet_grid(. ~ model) +
    labs(
      title = "Estimated eGFR slopes by group \n(main vs +CKD sensitivity)",
      x = NULL,
      y = expression(paste("Slope (mL/min/1.73 m"^2," per year), 95% CI")),
      fill = "Group"
    ) +
    scale_fill_manual(values = c(
      nonAKD = "#E41A1C",          # 赤
      Recovery = "#4DAF4A",        # 緑
      `Non-Recovery` = "#377EB8"   # 青
    )) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor  = element_blank(),
      strip.background  = element_rect(fill = "grey95", colour = NA),
      strip.text        = element_text(size = 10, face = "bold"),
      legend.position   = "bottom",
      legend.key.height = unit(0.5, "lines"),
      legend.key.width  = unit(1.2, "lines"),
      axis.title.y      = element_text(margin = margin(r = 6)),
      plot.title        = element_text(hjust = 0, face = "bold", lineheight = 1.1,
                                       margin = margin(b = 8)),
      plot.margin       = margin(t = 6, r = 8, b = 6, l = 6)
    )
  p <- p +
    scale_x_discrete(labels = NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  
  # ----- 保存サイズの基準：主要誌の1/2カラム想定 -----
  # 1カラム ~85 mm、2カラム ~180 mm を inch に変換
  w_1col_in <- 85  / 25.4   # 3.35 in
  w_2col_in <- 180 / 25.4   # 7.09 in
  h_in      <- 120 / 25.4   # 4.72 in（縦横比は誌面に合わせて調整可）
  setwd("E:/R")  
  # ----- 1) ベクター：PDF（フォント埋め込み） -----
  ggsave("Figure5_eGFR_slopes_main_vs_CKD.pdf",
         plot   = p,
         device = cairo_pdf,          # フォントを埋め込み
         width  = w_2col_in, height = h_in, units = "in")
  
  # （必要ならEPS）
  ggsave("Figure5_eGFR_slopes_main_vs_CKD.eps",
         plot   = p,
         device = cairo_ps,           # ベクターEPS
         width  = w_2col_in, height = h_in, units = "in",
         fallback_resolution = 600)
  
  # ----- 2) ラスター：TIFF 600 dpi（投稿規定で頻出） -----
  ggsave("Figure5_eGFR_slopes_main_vs_CKD_600dpi.tiff",
         plot   = p,
         device = "tiff",
         compression = "lzw",         # 可逆圧縮（多くの誌で可）
         dpi    = 600,
         width  = w_2col_in, height = h_in, units = "in",
         bg     = "white")            # 透過を避ける雑誌が多い
  
  # ----- 3) 1カラム幅のバージョンも出力（必要に応じて） -----
  ggsave("Figure5_eGFR_slopes_main_vs_CKD_1col.pdf",
         plot = p, device = cairo_pdf, width = w_1col_in, height = h_in, units = "in")
  
  ggsave("Figure5_eGFR_slopes_main_vs_CKD_1col_600dpi.tiff",
         plot = p, device = "tiff", compression = "lzw", dpi = 600,
         width = w_1col_in, height = h_in, units = "in", bg = "white")
}
{
  # ---- 前提：longdat_sens（全期間データ；years_c済）まで既に作成済み ----
  #   ここから ≤1年データを作って再中心化 → モデル当て → ALL/≤1y を結合して作図
  
  library(dplyr)
  library(nlme)
  library(multcomp)
  library(purrr)
  library(ggplot2)
  
  # ========= 0) 汎用：群別スロープと95%CI =========
  slope_ci_by_group_any <- function(fit, model_label) {
    cf <- names(fixef(fit))
    time_term <- if ("years_c" %in% cf) "years_c" else "years_from_time0"
    term_rec  <- cf[grepl(paste0("^", time_term, ":.*Recovery$"),       cf)]
    term_non  <- cf[grepl(paste0("^", time_term, ":.*Non-Recovery$"),   cf)]
    
    mk_vec <- function(g){
      v <- setNames(rep(0, length(cf)), cf)
      if (time_term %in% cf) v[time_term] <- 1
      if (g == "Recovery"     && length(term_rec) == 1) v[term_rec] <- 1
      if (g == "Non-Recovery" && length(term_non) == 1) v[term_non] <- 1
      v
    }
    
    L <- rbind(nonAKD = mk_vec("nonAKD"),
               Recovery = mk_vec("Recovery"),
               `Non-Recovery` = mk_vec("Non-Recovery"))
    
    ci <- suppressMessages(confint(multcomp::glht(fit, linfct = L)))
    tibble::tibble(
      group    = rownames(L),
      estimate = ci$confint[, "Estimate"],
      lower    = ci$confint[, "lwr"],
      upper    = ci$confint[, "upr"],
      model    = model_label
    )
  }
  
  # ========= 1) ≤1年データ作成（フィルタ → years_c を再中心化 → 安定化） =========
  longdat_sens_1y <- longdat_sens %>%                 # longdat_sens は全期間版（years_c あり）
    filter(years_from_time0 <= 1) %>%                 # ★ 1年以内に制限
    group_by(id) %>%
    mutate(n_time_1y = n_distinct(date),
           var_t_1y  = var(years_from_time0, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(n_time_1y >= 2, !is.na(var_t_1y), var_t_1y > 0) %>%
    mutate(years_c = as.numeric(scale(years_from_time0, center = TRUE, scale = FALSE)))  # ★ 再中心化
  
  # ========= 2) モデル当て（全期間：既存） =========
  # 既に fit_base_sens / fit_level_adj_sens / fit_slope_adj_sens
  #      fit_level_adj_sens_ckd / fit_slope_adj_sens_ckd がある前提
  
  # ========= 3) モデル当て（≤1年：同仕様で再推定） =========
  fit_base_sens_1y <- lme(
    egfr ~ years_c * jin_label_sens + time0_egfr_c - 1,
    random = list(id = pdDiag(~ 1 + years_c)),
    data = longdat_sens_1y, na.action = na.omit, method = "REML",
    control = lmeControl(opt="optim", msMaxIter=1e5, maxIter=1e5, niterEM=30, tolerance=1e-6)
  )
  
  fit_level_adj_sens_1y <- update(
    fit_base_sens_1y,
    . ~ . + age_c + sex + arb_acei_use +
      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15
  )
  
  fit_slope_adj_sens_1y <- update(
    fit_level_adj_sens_1y,
    . ~ . + years_c:(age_c + arb_acei_use +
                       dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 +
                       dn10 + dn12 + dn13 + dn14 + dn15)
  )
  
  # --- CKD含む ≤1年版（主効果 / 交互作用） ---
  fit_level_adj_sens_ckd_1y <- update(
    fit_base_sens_1y,
    . ~ . + age_c + sex + arb_acei_use + CKD_status +
      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15
  )
  
  fit_slope_adj_sens_ckd_1y <- update(
    fit_level_adj_sens_ckd_1y,
    . ~ . + years_c:(age_c + arb_acei_use + CKD_status +
                       dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 +
                       dn10 + dn12 + dn13 + dn14 + dn15)
  )
  
  # ========= 4) スロープ抽出（ALL / ≤1年）を結合 =========
  models_all <- list(
    "Base"                 = fit_base_sens,
    "Level-adjusted"       = fit_level_adj_sens,
    "slope-adjusted"       = fit_slope_adj_sens,
    "Level-adjusted \n+ CKD" = fit_level_adj_sens_ckd,
    "Slope-adjusted \n+ CKD" = fit_slope_adj_sens_ckd
  )
  
  models_1y <- list(
    "Base"                 = fit_base_sens_1y,
    "Level-adjusted"       = fit_level_adj_sens_1y,
    "Slope-adjusted"       = fit_slope_adj_sens_1y,
    "Level-adjusted \n+ CKD" = fit_level_adj_sens_ckd_1y,
    "Slope-adjusted \n+ CKD" = fit_slope_adj_sens_ckd_1y
  )
  
  df_all <- imap_dfr(models_all, ~ slope_ci_by_group_any(.x, .y)) %>%
    mutate(window = "All period")
  df_1y  <- imap_dfr(models_1y,  ~ slope_ci_by_group_any(.x, .y)) %>%
    mutate(window = "\u22641 year")     # "≤1 year"
  
  df_slope <- bind_rows(df_1y, df_all) %>%
    mutate(
      # モデル名の表記ゆれを統一
      model = dplyr::recode(model,
                            `slope-adjusted` = "Slope-adjusted"),
      # レベルも統一（大文字Sで揃える）
      model  = factor(model, levels = c(
        "Base","Level-adjusted","Slope-adjusted","Level-adjusted \n+ CKD","Slope-adjusted \n+ CKD"
      )),
      group  = factor(group,  levels = c("nonAKD","Recovery","Non-Recovery")),
      window = factor(window, levels = c("\u22641 year","All period"))
    )
  df_slope <- df_slope %>%
    mutate(model = forcats::fct_recode(model,
                                       "Slope-adjusted \n+ CKD" = "Slope-adjusted \n+ CKD",
                                       "Level-adjusted \n+ CKD" = "Level-adjusted \n+ CKD"
    ))
  
  p <- p +
    facet_grid(window ~ model) +
    theme(strip.text.x = element_text(size = 9, face = "bold", lineheight = 1.05),
          panel.spacing.x = unit(1.4, "lines"))
  
  # ========= 5) 作図：window（行） × model（列）の2次元ファセット =========
  p <- ggplot(df_slope, aes(x = group, y = estimate, fill = group)) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
    geom_col(width = 0.65) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 0.5) +
    facet_grid(window ~ model) +   # ★ 行=窓, 列=モデル（上段：≤1年／下段：全期間）
    labs(
      title = "Estimated eGFR slopes by group \n(≤1 year vs All period; main and +CKD sensitivity)",
      x = NULL,
      y = expression(paste("Slope (mL/min/1.73 m"^2," per year), 95% CI")),
      fill = "Group"
    ) +
    scale_fill_manual(values = c(
      nonAKD = "#E41A1C", Recovery = "#4DAF4A", `Non-Recovery` = "#377EB8"
    )) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor  = element_blank(),
      strip.background  = element_rect(fill = "grey95", colour = NA),
      strip.text        = element_text(size = 10, face = "bold"),
      legend.position   = "bottom",
      legend.key.height = unit(0.5, "lines"),
      legend.key.width  = unit(1.2, "lines"),
      axis.title.y      = element_text(margin = margin(r = 6)),
      plot.title        = element_text(hjust = 0, face = "bold", lineheight = 1.1,
                                       margin = margin(b = 8)),
      plot.margin       = margin(t = 6, r = 8, b = 6, l = 6),
      panel.spacing.x   = unit(1.0, "lines"),
      panel.spacing.y   = unit(1.2, "lines")
    ) +
    # X軸ラベルを消す場合は以下を有効化
    scale_x_discrete(labels = NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  # ========= 6) 保存（幅は5列ファセットを想定してやや広め） =========
  w_2col_in <- 180 / 25.4  # 7.09 in
  h_in      <- 160 / 25.4  # 行が2段なので少し背を高く（例：160mm）
  setwd("E:/R")  
  ggsave("Figure5_eGFR_slopes_ALL_vs_1year.pdf",
         plot = p, device = cairo_pdf, width = w_2col_in, height = h_in, units = "in")
  
  ggsave("Figure5_eGFR_slopes_ALL_vs_1year_600dpi.tiff",
         plot = p, device = "tiff", compression = "lzw", dpi = 600,
         width = w_2col_in, height = h_in, units = "in", bg = "white")
  
}
}

#95%信頼区間、有意差あるかどうか####
library(multcomp)
library(dplyr)
library(ggplot2)
library(purrr)
library(grid)

# ---- 1) 各モデルに対して contrasts（Recovery−nonAKD, Non-Recovery−nonAKD）を計算 ----
contrast_slope_diff <- function(fit, model_label, window_label) {
  cf <- names(fixef(fit))
  time_term <- if ("years_c" %in% cf) "years_c" else "years_from_time0"
  term_rec <- cf[grepl(paste0("^", time_term, ":.*Recovery$"), cf)]
  term_non <- cf[grepl(paste0("^", time_term, ":.*Non-Recovery$"), cf)]
  
  # 対比行列
  L <- rbind(
    "Recovery − nonAKD"     = c(ifelse(time_term %in% cf, 0, 0), rep(0, length(cf))),
    "Non-Recovery − nonAKD" = c(ifelse(time_term %in% cf, 0, 0), rep(0, length(cf)))
  )
  L <- matrix(0, nrow = 2, ncol = length(cf), dimnames = list(c("Recovery − nonAKD","Non-Recovery − nonAKD"), cf))
  if (length(term_rec) > 0) L[1, term_rec] <- 1
  if (length(term_non) > 0) L[2, term_non] <- 1
  
  # 推定値・CI・p値
  out <- summary(multcomp::glht(fit, linfct = L))
  ci  <- confint(multcomp::glht(fit, linfct = L))
  
  tibble(
    contrast = rownames(L),
    diff     = ci$confint[, "Estimate"],
    lower    = ci$confint[, "lwr"],
    upper    = ci$confint[, "upr"],
    p_value  = out$test$pvalues,
    model    = model_label,
    window   = window_label
  )
}

# ---- 2) すべてのモデル・ウィンドウに対して実行（例：CKD感度解析）----
df_contrast <- purrr::imap_dfr(models_all_ckd, ~ contrast_slope_diff(.x, .y, "All period")) %>%
  bind_rows(purrr::imap_dfr(models_1y_ckd, ~ contrast_slope_diff(.x, .y, "≤1 year"))) %>%
  mutate(
    window = factor(window, levels = c("≤1 year","All period")),
    model  = factor(model,  levels = c("Level-adjusted + CKD","Slope-adjusted + CKD")),
    signif = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )
}
