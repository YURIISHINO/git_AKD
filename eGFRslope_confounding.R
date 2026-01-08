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

# ********参考: 論文スタイルの改良版（サンプルサイズ・差分・P値付き）********

# Step 1: サンプルサイズの取得（longdatから計算）
sample_sizes <- longdat %>%
  mutate(
    window = case_when(
      years_from_time0 <= 1 ~ "≤1 year",
      years_from_time0 <= 3 ~ "≤3 years",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(window)) %>%
  distinct(id, jin_label, window) %>%
  count(jin_label, window, name = "n") %>%
  rename(group = jin_label) %>%
  mutate(group = as.character(group))

# Step 2: 群間差（Difference）の計算
# まず、df_contrastの中身を確認
cat("\n=== df_contrast の window 列のユニーク値 ===\n")
print(unique(df_contrast$window))
cat("\n=== df_contrast の model 列のユニーク値 ===\n")
print(unique(df_contrast$model))
cat("\n=== df_contrast の contrast 列のユニーク値 ===\n")
print(unique(df_contrast$contrast))

# df_contrastから≤1年と≤3年のSlope-adjustedを抽出
# 注意: df_contrastのmodelは"③ Slope-adjusted"（丸数字付き）
# nonAKDとの比較のみを抽出（Non-Recovery − Recoveryは除外）
differences <- df_contrast %>%
  dplyr::filter(
    window %in% c("≤1 year", "≤3 years"),
    model == "③ Slope-adjusted",  # 丸数字付きに修正
    contrast %in% c("Recovery − nonAKD", "Non-Recovery − nonAKD")  # nonAKDとの比較のみ
  ) %>%
  dplyr::mutate(
    group = case_when(
      contrast == "Recovery − nonAKD" ~ "Recovery",
      contrast == "Non-Recovery − nonAKD" ~ "Non-Recovery",
      TRUE ~ contrast
    )
  ) %>%
  dplyr::select(window, group,
                diff_value = diff,
                diff_lower = lower,
                diff_upper = upper,
                diff_p = p_value)

cat("\n=== differences の行数 ===\n")
print(nrow(differences))
cat("\n=== differences の中身 ===\n")
print(differences)

# nonAKD用のReference行を追加
differences <- bind_rows(
  tibble(
    window = rep(c("≤1 year", "≤3 years"), each = 1),
    group = "nonAKD",
    diff_value = NA_real_,
    diff_lower = NA_real_,
    diff_upper = NA_real_,
    diff_p = NA_real_
  ),
  differences
)

# Step 3: データ統合（groupの順序を正しく設定）
# 注意: df_plotのmodelは"Slope-adjusted"（丸数字なし）
df_fig4_enhanced <- df_plot %>%
  filter(
    window %in% c("≤1 year", "≤3 years"),
    model == "Slope-adjusted"  # df_plotでは丸数字なし
  ) %>%
  mutate(group = as.character(group)) %>%
  left_join(sample_sizes, by = c("group", "window")) %>%
  left_join(differences, by = c("group", "window"), suffix = c("", "_diff")) %>%
  # groupの順序を正しく設定（nonAKD, Recovery, Non-Recovery）
  mutate(
    group = factor(group, levels = c("nonAKD", "Recovery", "Non-Recovery"))
  )

cat("\n=== df_fig4_enhanced の結合結果確認 ===\n")
print(df_fig4_enhanced)

# Step 4: プロット作成
# Y軸の最小値と最大値を事前に計算（警告回避のため）
y_min <- min(df_fig4_enhanced$lower, na.rm = TRUE)
y_max <- max(df_fig4_enhanced$upper, na.rm = TRUE)

p_1y_3y_enhanced <- ggplot(
  df_fig4_enhanced,
  aes(x = group, y = estimate, fill = group)
) +
  # 棒グラフ
  geom_col(width = 0.7, alpha = 0.85) +

  # エラーバー
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.25,
    linewidth = 0.7
  ) +

  # サンプルサイズ（棒の上）
  geom_text(
    aes(
      y = pmax(upper, 0) + 0.5,
      label = paste0("n=", n)
    ),
    size = 3.5,
    fontface = "bold",
    color = "grey20"
  ) +

  # Slope値と95%CI（棒の下、Y=0より下）
  geom_text(
    aes(
      y = pmin(lower, 0) - 0.8,
      label = sprintf("%.2f\n(%.2f, %.2f)", estimate, lower, upper)
    ),
    size = 3,
    lineheight = 0.9,
    color = "grey10"
  ) +

  # Difference（グラフ最下部）
  geom_text(
    aes(
      label = ifelse(
        group == "nonAKD",
        "Reference",
        sprintf("Diff: %.2f (%.2f, %.2f)", diff_value, diff_lower, diff_upper)
      )
    ),
    y = y_min - 3.5,
    size = 3,
    fontface = "italic",
    color = "grey30"
  ) +

  # P-value（Differenceの下）
  geom_text(
    aes(
      label = ifelse(
        group == "nonAKD",
        "",
        case_when(
          is.na(diff_p) ~ "",
          diff_p < 0.001 ~ "p<0.001",
          diff_p < 0.01 ~ sprintf("p=%.3f", diff_p),
          TRUE ~ sprintf("p=%.2f", diff_p)
        )
      )
    ),
    y = y_min - 4.5,
    size = 3,
    fontface = "bold",
    color = "grey20"
  ) +

  # ファセット
  facet_grid(. ~ window) +

  # ラベル
  labs(
    x = NULL,
    y = expression(paste("Mean change in eGFR (mL/min/1.73 m"^2," per year)")),
    fill = "Group"
  ) +

  # 色（シンプルで直感的）
  scale_fill_manual(
    values = c(
      nonAKD         = "#95A5A6",
      Recovery       = "#2ECC71",
      `Non-Recovery` = "#E74C3C"
    ),
    labels = c(
      nonAKD = "No AKD",
      Recovery = "AKD with Recovery",
      `Non-Recovery` = "AKD without Recovery"
    )
  ) +

  # Y軸範囲（Difference/P値表示のため下方向に拡張）
  coord_cartesian(
    ylim = c(
      y_min - 5.5,
      max(df_fig4_enhanced$upper, na.rm = TRUE) + 1.5
    ),
    clip = "off"
  ) +

  # テーマ（グリッドなし、シンプル）
  theme_classic(base_size = 13) +
  theme(
    # パネル（上部に境界線を表示）
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.6),
    panel.background = element_rect(fill = "white", color = NA),

    # ファセットラベル（背景色なし）
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold", color = "grey20",
                              margin = margin(t = 2, b = 5)),
    strip.placement = "outside",
    panel.spacing.x = unit(1.5, "lines"),

    # 軸
    axis.line.y = element_blank(),  # panel.borderを使うのでaxis.lineは不要
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    axis.text.y = element_text(size = 11),

    # 凡例
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(1.2, "lines"),
    legend.background = element_blank(),

    # マージン（下部を広げてDifference/P値の表示スペース確保、上部も少し広げる）
    plot.margin = margin(t = 15, r = 15, b = 20, l = 10)
  )

# プレビュー
print(p_1y_3y_enhanced)

# 保存（論文用）
ggsave(
  "Figure4_eGFR_slopes_1y_3y_publication_style.pdf",
  plot = p_1y_3y_enhanced,
  width = 200, height = 140, units = "mm",
  device = cairo_pdf, dpi = 300
)

ggsave(
  "Figure4_eGFR_slopes_1y_3y_publication_style.tiff",
  plot = p_1y_3y_enhanced,
  width = 200, height = 140, units = "mm",
  device = "tiff", dpi = 600, compression = "lzw"
)

# ************************************************************************


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

