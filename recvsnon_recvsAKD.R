#必要な情報をすべて含有するcsvファイル(jin1_Eligible? primay_ESKDやdeath情報や感度分析に必要な情報などすべて含有)を作成し、適当なdfを作成
library(readr)
library(dplyr)
setwd("E:/R")
jin1_Eligibile <- read_csv("jin1_Eligibile.csv", locale = locale(encoding = "SHIFT-JIS"))

#集計#####
jin1_Eligibile %>%
  group_by(jin_status, exclude) %>%
  summarise(n = n(), .groups = "drop")
jin1_Eligibile %>%
  filter(exclude == "include") %>%
  group_by(AKD_status) %>%
  summarise(n_id = n_distinct(id)) %>%
  arrange(desc(n_id))

jin1_Eligibile %>%
  filter(exclude == "include") %>%
  group_by(jin_status) %>%
  summarise(n_unique_ids = n_distinct(id), .groups = "drop")

#観察期間計算####
library(dplyr)
# 観察期間（日数 → 月数へ変換）
follow_summary <- jin1_Eligibile %>%
  mutate(
    follow_days = as.numeric(last_follow_composite - index_date),
    follow_months = follow_days / 30.44  # 日→月換算（平均1か月=30.44日）
  ) %>%
  summarise(
    median_follow = median(follow_months, na.rm = TRUE),
    q1 = quantile(follow_months, 0.25, na.rm = TRUE),
    q3 = quantile(follow_months, 0.75, na.rm = TRUE)
  )

# 出力文の作成
result_text <- with(follow_summary,
                    paste0(
                      "観察期間中央値は", round(median_follow, 1),
                      "ヶ月（四分位範囲：", round(q1, 1), "–", round(q3, 1), "ヶ月）"
                    )
)

# 結果を表示
cat(result_text, "\n")

#eGFR・線形混合効果モデルを走らせる、また解析に必要なcodeのみ####
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

#Time Window Approachの実装
# データ準備
akd_time_m <- jin1_inclusion %>%
  mutate(
    years_from_time0 = as.numeric(date - time0) / 365.25  # 年単位に変換
  )%>%
  filter(years_from_time0 >= 0)

# 測定間隔の中央値を確認
measurement_intervals <- akd_time_m %>%
  arrange(id, years_from_time0) %>%
  group_by(id) %>%
  mutate(interval = years_from_time0 - lag(years_from_time0)) %>%
  filter(!is.na(interval) & interval > 0)
median_interval <- median(measurement_intervals$interval, na.rm = TRUE)
# years_from_time0 の範囲取得
time_range <- range(akd_time_m$years_from_time0, na.rm = TRUE)
# 規定時点の設定
if (time_range[2] <= 2) {
  target_timepoints <- seq(0, 2, by = 0.25)  # 3ヶ月ごと
} else if (time_range[2] <= 5) {
  target_timepoints <- seq(0, ceiling(time_range[2]), by = 0.5)  # 6ヶ月ごと
} else {
  target_timepoints <- seq(0, ceiling(time_range[2]), by = 1)    # 1年ごと
}
print(paste("規定時点:", paste(target_timepoints, collapse = ", ")))
# ウィンドウ幅の設定
window_width <- median_interval * 0.75
# Time Window Approach 関数の定義
create_time_window_data <- function(data, targets, window) {
  result_list <- list()
  for (i in seq_along(targets)) {
    target <- targets[i]
    window_data <- data %>%
      filter(abs(years_from_time0 - target) <= window) %>%
      group_by(id, jin_label) %>%
      slice_min(abs(years_from_time0 - target), n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      mutate(target_time = target)
    result_list[[i]] <- window_data
  }
  bind_rows(result_list)
}

# 実行：Time Window データの生成
window_data <- create_time_window_data(akd_time_m, target_timepoints, window_width)

#slope描出に必要なcodeのみ####
{# 線形混合モデル　1年####
fit_window_within_1year <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr - 1,
  random = list(id = pdSymm(form = ~ 1 + years_from_time0)),
  na.action = na.omit,
  data = filter(akd_time_m, years_from_time0 <= 1),
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)

# 実測値を使ったΔeGFRを計算
window_data_obs <- window_data %>%
  group_by(id) %>%
  mutate(time0_egfr = egfr[target_time == 0][1]) %>%
  ungroup() %>%
  mutate(observed_diff = egfr - time0_egfr)

# 集団平均と95%信頼区間を計算
summary_obs_aligned <- window_data_obs %>%
  group_by(jin_label, target_time) %>%
  summarise(
    mean_diff = mean(observed_diff, na.rm = TRUE),
    sd = sd(observed_diff, na.rm = TRUE),
    n = sum(!is.na(observed_diff)),
    se = sd / sqrt(n),
    ci_lower = mean_diff - 1.96 * se,
    ci_upper = mean_diff + 1.96 * se,
    .groups = "drop"
  )

# ΔeGFRの実測値によるプロット（X軸制限あり）
ggplot(summary_obs_aligned, aes(x = target_time, y = mean_diff, color = jin_label)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.05) +
  scale_x_continuous(limits = c(0, 10)) +  # X軸を0～10年に限定
  coord_cartesian(ylim = c(-30, 0)) +  #y軸を0～-30に限定
  labs(
    title = "Observed ΔeGFR from time0 within 1_year",
    x = "Time from time0 (years)",
    y = expression(Delta~eGFR~"(mL/min/1.73m2)")
  ) +
  theme_minimal()
# 線形混合モデル　2年####
fit_window_within_2year <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr - 1,
  random = list(id = pdSymm(form = ~ 1 + years_from_time0)),
  na.action = na.omit,
  data = filter(akd_time_m, years_from_time0 <= 2),
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)
# 線形混合モデル　3年####
fit_window_within_3year <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr - 1,
  random = list(id = pdSymm(form = ~ 1 + years_from_time0)),
  na.action = na.omit,
  data = filter(akd_time_m, years_from_time0 <= 3),
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)

# 線形混合モデル 全期間#####
fit_window <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr - 1,
  random = list(id = pdSymm(form = ~ 1 + years_from_time0)),
  na.action = na.omit,
  data = akd_time_m,
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B"))



# 実測値を使ったΔeGFRを計算
window_data_obs <- window_data %>%
  group_by(id) %>%
  mutate(time0_egfr = egfr[target_time == 0][1]) %>%
  ungroup() %>%
  mutate(observed_diff = egfr - time0_egfr)

# 集団平均と95%信頼区間を計算
summary_obs_aligned <- window_data_obs %>%
  group_by(jin_label, target_time) %>%
  summarise(
    mean_diff = mean(observed_diff, na.rm = TRUE),
    sd = sd(observed_diff, na.rm = TRUE),
    n = sum(!is.na(observed_diff)),
    se = sd / sqrt(n),
    ci_lower = mean_diff - 1.96 * se,
    ci_upper = mean_diff + 1.96 * se,
    .groups = "drop"
  )

# ΔeGFRの実測値によるプロット（X軸制限あり）
ggplot(summary_obs_aligned, aes(x = target_time, y = mean_diff, color = jin_label)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.05) +
  scale_x_continuous(limits = c(0, 10)) +  # X軸を0～10年に限定
  coord_cartesian(ylim = c(-30, 0)) +  #y軸を0～-30に限定
  labs(
    title = "Observed ΔeGFR from time0 all time",
    x = "Time from time0 (years)",
    y = expression(Delta~eGFR~"(mL/min/1.73m2)")
  ) +
  theme_minimal()

#縮尺をそろえる2####

# 必要パッケージ
library(dplyr)
library(ggplot2)
library(multcomp)
library(tibble)
library(forcats)

# 1) スロープ抽出ヘルパー
extract_slope_df <- function(fit, horizon_label){
  cn <- names(fixef(fit))
  p  <- length(cn)
  
  v_main <- rep(0, p)
  v_main[which(cn == "years_from_time0")] <- 1
  
  # 相互作用（環境で「Non-Recovery」「Non.Recovery」の揺れを許容）
  idx_rec    <- which(cn == "years_from_time0:jin_labelRecovery")
  idx_nonrec <- which(cn %in% c("years_from_time0:jin_labelNon-Recovery",
                                "years_from_time0:jin_labelNon.Recovery"))
  
  v_rec    <- v_main; if (length(idx_rec)==1)    v_rec[idx_rec]       <- 1
  v_nonrec <- v_main; if (length(idx_nonrec)==1) v_nonrec[idx_nonrec] <- 1
  
  ci_nonAKD <- confint(glht(fit, linfct = rbind(nonAKD        = v_main),    vcov = vcov(fit)))$confint
  ci_rec    <- confint(glht(fit, linfct = rbind(Recovery      = v_rec),      vcov = vcov(fit)))$confint
  ci_nonrec <- confint(glht(fit, linfct = rbind(`Non-Recovery`= v_nonrec),   vcov = vcov(fit)))$confint
  
  tibble(
    Group = factor(c("nonAKD","Recovery","Non-Recovery"),
                   levels = c("nonAKD","Recovery","Non-Recovery")),
    Estimate = c(ci_nonAKD[,"Estimate"], ci_rec[,"Estimate"], ci_nonrec[,"Estimate"]),
    CI_Lower = c(ci_nonAKD[,"lwr"],      ci_rec[,"lwr"],      ci_nonrec[,"lwr"]),
    CI_Upper = c(ci_nonAKD[,"upr"],      ci_rec[,"upr"],      ci_nonrec[,"upr"]),
    Horizon  = horizon_label
  )
}

# 2) 4つのモデルからデータ作成（名称はご利用中のオブジェクトに合わせて）
#    - 1年:    fit_window_within_1year
#    - 2年:    fit_window_within_2year
#    - 3年:    fit_window_within_3year
#    - 全期間: fit_window
slope_all <- dplyr::bind_rows(
  extract_slope_df(fit_window_within_1year, "within1 year"),
  extract_slope_df(fit_window_within_2year, "within2 years"),
  extract_slope_df(fit_window_within_3year, "within3 years"),
  extract_slope_df(fit_window,              "All period")
) %>%
  mutate(
    # 表示順（左→右）
    Horizon = factor(Horizon, levels = c("within1 year","within2 years","within3 years","All period"))
  )

# 3) 共通Y軸レンジを決定
y_min <- floor(min(slope_all$CI_Lower, na.rm = TRUE)) - 0.5
y_max <- ceiling(max(slope_all$CI_Upper, na.rm = TRUE)) + 0.5

# 4) 1枚の図にfacetで並べる
fill_cols <- c("nonAKD"="#E41A1C", "Recovery"="#1B9E77", "Non-Recovery"="#377EB8")
p <- ggplot(slope_all, aes(x = Group, y = Estimate, fill = Group)) +
  geom_col(width = 0.7) +  # 棒枠（任意）
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.15, linewidth = 0.6) +
  facet_wrap(~ Horizon, nrow = 1) +
  scale_fill_manual(values = fill_cols) +
  scale_y_continuous(limits = c(y_min, y_max)) +
  labs(
    title = "Estimated eGFR Slopes by Group across Time Horizons",
    x = NULL,
    y = "Slope (mL/min/1.73 m2/year) with 95% CI"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8) # ← 追加
  )
print(p)


# （任意）保存：論文用の高解像度
ggsave("Fig_Slopes_by_Group_TimeHorizons_facet.tiff", plot = p,
       width = 10.5, height = 3.2, units = "in", dpi = 600, compression = "lzw")

}
#アウトカム解析に必要なcodeのみ####
#jin1_Eligibileから一意のidだけ抽出
jin1_Eligibile_unique_id <- jin1_Eligibile %>%
  filter(exclude == "include", jin_status %in% c("AKD", "nonAKD")) %>%
  group_by(id) %>%
  arrange(index_date) %>%
  slice(1) %>%
  ungroup()
#死亡についてのカプランマイヤー####
install.packages("survminer")
library(survival)
library(survminer)
library(dplyr)

#死亡についてのカプランマイヤーデータ整形
jin1_Eligibile_death <- jin1_Eligibile_unique_id %>%
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
    jin_label = factor(jin_label, levels = c("nonAKD", "Recovery", "Non-Recovery")),  # 順序付け
    time_years = as.numeric(last_follow_death - index_plus_210) / 365.25
  ) %>%
  filter(!is.na(jin_label), time_years >= 0)

fit_death <- survfit(Surv(time_years, primary_death) ~ jin_label, data = jin1_Eligibile_death)

ggsurvplot(fit_death, data = jin1_Eligibile_death,
           fun = "event",                     # 1 - survival（累積死亡率）
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE,
           xlim = c(0, 10),                   # 10年まで表示
           ylim = c(0, 0.4),
           title = "death",
           xlab = "year",
           ylab = "event")

colnames(jin1_Eligibile)
#競合エンドポイントについてのカプランマイヤー 2025/8/5 ####
install.packages("survminer")
library(survival)
library(survminer)
library(dplyr)

#競合エンドポイントについてのカプランマイヤーデータ整形
jin1_Eligibile_composite <- jin1_Eligibile_unique_id %>%
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
    time_years = as.numeric(last_follow_composite - index_plus_210) / 365.25　
  ) %>%
  filter(!is.na(jin_label), time_years >= 0)

fit_composite <- survfit(Surv(time_years, primary_composite_event) ~jin_label, data = jin1_Eligibile_composite)

ggsurvplot(fit_composite, data = jin1_Eligibile_composite,
           fun = "event",                     # 1 - survival（累積死亡率）
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE,
           xlim = c(0, 10),                   # 10年まで表示
           ylim = c(0, 0.4),
           title = "composite_event",
           xlab = "year",
           ylab = "event")

#coxphに必要なcodeのみ####
# 必要パッケージ
library(dplyr)
library(survival)
library(broom)
library(ggplot2)
# （Recovery vs Non-Recoveryの直接比較を出すなら）
# library(multcomp)

# -------------------------------
# 1) 解析用データ（3群）の作成
# -------------------------------
# 1人1行（index_date当日レコードのみ・最初のindex_dateを採用）
jin1_Eligibile_unique_id <- jin1_Eligibile %>%
  filter(exclude == "include", jin_status %in% c("AKD", "nonAKD")) %>%
  group_by(id) %>%
  arrange(index_date, date, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

#確認
jin1_Eligibile_unique_id  %>%
  group_by(jin_status) %>%
  summarise(n_unique_ids = n_distinct(id), .groups = "drop")

jin1_Eligibile_cox_3group <- jin1_Eligibile_unique_id %>%
  # nonAKD と AKDのみ（あなたの方針に合わせて必要ならここは保持）
  filter(exclude == "include", jin_status %in% c("AKD", "nonAKD")) %>%
  # index_date当日の行だけ（同日複数行の可能性に備えてdistinct）
  #filter(date == index_date) %>%
  distinct(id, .keep_all = TRUE) %>%
  # 3群ラベル
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
    # 併用フラグ
    arb_acei_use = if_else(coalesce(arb, 0) == 1 | coalesce(acei, 0) == 1, 1L, 0L),
    # 追跡時間（年）：index_date+210日以降に統一
    time_years   = as.numeric(last_follow_death - index_plus_210) / 365.25
  ) %>%
  # 追跡開始以降のみ
  filter(!is.na(jin_label), !is.na(time_years), time_years >= 0)


#確認
jin1_Eligibile_cox_3group  %>%
  group_by(jin_label) %>%
  summarise(n_unique_ids = n_distinct(id), .groups = "drop")

# -------------------------------
# 2) Cox比例ハザード（3群）
# -------------------------------
cox_model_3group <- coxph(
  Surv(time_years, primary_death) ~
    jin_label + age + index_cre + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15 +
    CKD_status,
  data = jin1_Eligibile_cox_3group
)
confint(cox_model_3group)
exp(confint(cox_model_3group))


# RecoveryとNon-Recoveryの差を検定→有意差なし#####
fit <- coxph(Surv(time_years, primary_death) ~ jin_label + age + index_cre + arb_acei_use +
               dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 +
               dn12 + dn13 + dn14 + dn15 + CKD_status,
             data = jin1_Eligibile_cox_3group)
anova(fit)
# または contrast 検定
library(multcomp)
# Recovery － Non-Recovery = 0（＝差がない）を検定
g <- glht(fit, linfct = c("`jin_labelRecovery` - `jin_labelNon-Recovery` = 0"))
summary(g)          # z と p値（線形予測子スケール）
ci <- confint(g)    # 95%CI（線形予測子スケール）

# HRと95%CIに変換して表示
est    <- summary(g)$test$coef[1]
pval   <- summary(g)$test$pvalues[1]
HR     <- exp(est)
HR_CI  <- exp(ci$confint[1, c("lwr","upr")])
c(HR = HR, CI_low = HR_CI[1], CI_high = HR_CI[2], p = pval)



# （任意）PH仮定チェック
# print(cox.zph(cox_model_3group))
# --- 3群Cox : tidy + formula順でのフォレスト -------------------------------
library(broom)
library(dplyr)
library(ggplot2)

tidy3 <- broom::tidy(cox_model_3group, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    term_raw  = term,  # tidy()が吐く元の係数名を保持
    # 読みやすい表示名に置き換え（必要に応じて追加/調整OK）
    term_nice = dplyr::recode(term,
                              "jin_labelRecovery"       = "AKD Recovery vs nonAKD",
                              "jin_labelNon-Recovery"   = "AKD Non-Recovery vs nonAKD",
                              "age"                     = "Age",
                              "index_cre"               = "Index Creatinine",
                              "arb_acei_use"            = "ARB or ACEi Use",
                              "dn1"  = "CHF",           "dn3"  = "Rheumatologic", "dn4"  = "Malignancy",
                              "dn5"  = "Liver Disease", "dn6"  = "Peptic Ulcer",  "dn7"  = "MI",
                              "dn8"  = "Renal Disease", "dn9"  = "Metastatic Cancer",
                              "dn10" = "Diabetes",      "dn12" = "Stroke",        "dn13" = "Hemiplegia",
                              "dn14" = "PVD",           "dn15" = "COPD",
                              .default = term  # 上で指定していないものは元のまま（例: CKD_statusの水準など）
    )
  )

# 1) モデル式の項ラベルを取得（式に書いた順）
form_order <- attr(stats::terms(cox_model_3group), "term.labels")
# 2) その順に一致する係数（ダミーや多水準を含む）を展開
order_raw <- unlist(lapply(form_order, function(x) tidy3$term_raw[startsWith(tidy3$term_raw, x)]))
order_raw <- order_raw[order_raw %in% tidy3$term_raw]  # 念のため整合チェック

# 3) 表示名の順ベクトルを作成（上＝式の先頭項にしたいのでREVする）
display_order <- tidy3 %>%
  filter(term_raw %in% order_raw) %>%
  arrange(match(term_raw, order_raw)) %>%
  pull(term_nice)

tidy3 <- tidy3 %>%
  mutate(term_nice = factor(term_nice, levels = rev(display_order))) 
# ↑ グラフのY軸は下→上にレベル順で描かれるのでrev()で
#    「上から順に：式に書いた順」になります。
#    下から順にしたい場合はrev()を外してください。

# 必要なら見やすい表（この順のまま）
tidy3 %>%
  dplyr::select(term = term_nice, HR = estimate, `CI low` = conf.low, `CI high` = conf.high, p.value) %>%
  arrange(desc(term)) -> tidy3_table  # arrangeは任意。factorレベルが順序を担保

# フォレストプロット（順序固定）
ggplot(tidy3, aes(x = estimate, y = term_nice)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  scale_x_log10() +
  labs(
    x = "Hazard Ratio (log scale)",
    y = NULL,
    title = "Death: Adjusted Hazard Ratios (3-group Cox model)"
  ) +
  theme_minimal(base_size = 14)


#上記解析・描出を感度分析を行うのに必要なcodeのみ####
#①死亡のカプランマイヤー####
#データ整形
jin1_Eligibile_sens <- jin1_Eligibile_unique_id %>%
  filter(exclude == "include", jin_status %in% c("AKD", "nonAKD")) %>%
  mutate(
    group = case_when(
      jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 0 ~ NA_character_,  # 除外対象をNAに
      jin_status == "nonAKD" ~ "nonAKD",
      jin_status == "AKD" & `150_210recovery` == 1 ~ "Recovery",
      jin_status == "AKD" & `150_210recovery` == 2 ~ "Non-Recovery",
      jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 1 ~ "Recovery",
      jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 2 ~ "Non-Recovery",
      TRUE ~ NA_character_
    ),
    group = factor(group, levels = c("nonAKD", "Recovery", "Non-Recovery")), 
    time_years = as.numeric(last_follow_death - index_plus_210) / 365.25
  ) %>%
  filter(!is.na(group), time_years >= 0)

# モデル作成、グラフ（死亡）
fit_death <- survfit(Surv(time_years, primary_death) ~ group, data = jin1_Eligibile_sens)

ggsurvplot(fit_death, data = jin1_Eligibile_sens,
           fun = "event",                     # 1 - survival（累積死亡率）
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE,
           xlim = c(0, 10),                   # 10年まで表示
           ylim = c(0, 0.4),
           title = "Sensitivity Analysis(death)",
           xlab = "year",
           ylab = "event")

#②eGFRslope ####
window_data
#`150_210recovery` == 0,`90_150recovery` == 0を除外する
akd_no_recovery_ids <- jin1_Eligibile %>%
  filter(
    exclude == "include",
    jin_status == "AKD",
    `150_210recovery` == 0,
    `90_150recovery` == 0
  ) %>%
  distinct(id)
akd_time_m_sens <- akd_time_m %>%
  anti_join(akd_no_recovery_ids, by = "id") %>%
  mutate(jin_label = factor(jin_label,
                            levels = c("nonAKD", "Recovery", "Non-Recovery")))

window_data_sens<- window_data %>%
  anti_join(akd_no_recovery_ids, by = "id")

# 線形混合モデル
fit_window_sens <- lme(
  egfr ~ years_from_time0 * jin_label + time0_egfr - 1,
  random = list(id = pdSymm(form = ~ 1 + years_from_time0)),
  na.action = na.omit,
  data = akd_time_m_sens,
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B",
                       msVerbose = TRUE))
# ステップ1：絶対egfrを得る
window_data_sens <- window_data_sens %>%
  mutate(pred_egfr = predict(fit_window_sens, newdata = window_data_sens))


# ステップ2：time0時点の「モデル予測値」を患者ごとに抽出
pred_time0 <- window_data_sens %>%
  filter(target_time == 0) %>%
  dplyr::select(id, pred_egfr_time0 = pred_egfr)

# ステップ3：全体に time0予測値を結合し、ΔeGFRを計算
window_data_sens <- window_data_sens %>%
  group_by(id) %>%
  # idごとに target_time==0 の予測値を拾って列として作り直す
  mutate(pred_egfr_time0 = pred_egfr[target_time == 0][1]) %>%
  ungroup() %>%
  # 列参照を明示（.data）してアライン
  mutate(pred_egfr_aligned = .data$pred_egfr - .data$pred_egfr_time0)

# 集計
summary_aligned <- window_data_sens %>%
  group_by(jin_label, target_time) %>%
  summarise(
    mean_diff = mean(pred_egfr_aligned, na.rm = TRUE),
    sd = sd(pred_egfr_aligned, na.rm = TRUE),
    n = sum(!is.na(pred_egfr_aligned)),
    se = sd / sqrt(n),
    ci_lower = mean_diff - 1.96 * se,
    ci_upper = mean_diff + 1.96 * se,
    .groups = "drop"
  )

# プロット
ggplot(summary_aligned, aes(x = target_time, y = mean_diff, color = jin_label)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.05) +
  labs(
    title = "ΔeGFR from model-based time0 (adjusted by time0_egfr) sensitivity analysis",
    x = "Time from time0 (year
    s)",
    y = expression(Delta~eGFR~"(mL/min/1.73m2)")
  ) +
  scale_x_continuous(limits = c(0, 10)) +  # X軸を0～10年に限定
  theme_minimal()


#③eGFRslope固定効果の棒グラフ####
# 必要パッケージ
library(nlme)
library(dplyr)
library(ggplot2)
library(multcomp)
library(purrr)
library(tibble)

## 0) 前処理：因子水準の固定（再現性のため必須）
#   データは akd_time_m_sens を想定
dat0 <- akd_time_m_sens %>%
  mutate(
    jin_label = factor(jin_label, levels = c("nonAKD", "Recovery", "Non-Recovery"))
  )

## 1) LMEを当てる関数（期間でフィルタ → LME適合）
fit_lme_within <- function(df, horizon_years) {
  df_sub <- if (is.infinite(horizon_years)) df else dplyr::filter(df, years_from_time0 <= horizon_years)
  
  lme(
    egfr ~ years_from_time0 * jin_label + time0_egfr - 1,
    random = list(id = pdSymm(form = ~ 1 + years_from_time0)),
    na.action = na.omit,
    data = df_sub,
    method = "REML",
    control = lmeControl(
      maxIter = 1e8, msMaxIter = 1e8,
      opt = "optim", optimMethod = "L-BFGS-B",
      msVerbose = TRUE
    )
  )
}

## 2) スロープ（年あたり変化量）の推定値と95%CIを glht で取り出す関数
extract_slopes <- function(fit, horizon_label = "") {
  cn <- names(fixef(fit))
  p  <- length(cn)
  
  v_main <- rep(0, p)
  idx_main <- which(cn == "years_from_time0")
  if (length(idx_main) != 1) stop("Could not find 'years_from_time0' in fixed effects.")
  v_main[idx_main] <- 1
  
  # 相互作用係数名（環境差を許容）
  idx_rec    <- which(cn == "years_from_time0:jin_labelRecovery")
  idx_nonrec <- which(cn %in% c("years_from_time0:jin_labelNon-Recovery",
                                "years_from_time0:jin_labelNon.Recovery"))
  
  v_rec <- v_main
  if (length(idx_rec) == 1) v_rec[idx_rec] <- 1 else warning("Recovery interaction term not found.")
  
  v_nonrec <- v_main
  if (length(idx_nonrec) == 1) v_nonrec[idx_nonrec] <- 1 else warning("Non-Recovery interaction term not found.")
  
  ci_nonAKD <- confint(glht(fit, linfct = rbind(nonAKD = v_main),            vcov = vcov(fit)))$confint
  ci_rec    <- confint(glht(fit, linfct = rbind(Recovery = v_rec),            vcov = vcov(fit)))$confint
  ci_nonrec <- confint(glht(fit, linfct = rbind(`Non-Recovery` = v_nonrec),   vcov = vcov(fit)))$confint
  
  bind_rows(
    tibble(Group = factor("nonAKD",       levels = c("nonAKD","Recovery","Non-Recovery")),
           Estimate = ci_nonAKD[,"Estimate"], CI_Lower = ci_nonAKD[,"lwr"], CI_Upper = ci_nonAKD[,"upr"]),
    tibble(Group = factor("Recovery",     levels = c("nonAKD","Recovery","Non-Recovery")),
           Estimate = ci_rec[,"Estimate"],    CI_Lower = ci_rec[,"lwr"],    CI_Upper = ci_rec[,"upr"]),
    tibble(Group = factor("Non-Recovery", levels = c("nonAKD","Recovery","Non-Recovery")),
           Estimate = ci_nonrec[,"Estimate"], CI_Lower = ci_nonrec[,"lwr"], CI_Upper = ci_nonrec[,"upr"])
  ) %>%
    mutate(Horizon = horizon_label)
}

## 3) 全期間＋1/2/3年で一括実行 → 結合
horizons       <- c(1, 2, 3, Inf)
horizon_labels <- c("within1 year", "within2 years", "within3 years","All period")

fits <- map(horizons, ~ fit_lme_within(dat0, .x))
slopes_all <- map2_dfr(fits, horizon_labels, ~ extract_slopes(.x, .y))

order_levels <- c("within1 year", "within2 years", "within3 years", "All period")
slopes_all <- slopes_all %>%
  mutate(Horizon = factor(Horizon, levels = order_levels))


## 4) （任意）各期間の症例数（ユニーク id）をサブタイトルに表示
n_by_horizon <- map_int(horizons, function(h) {
  if (is.infinite(h)) dat0 %>% distinct(id) %>% nrow()
  else dat0 %>% filter(years_from_time0 <= h) %>% distinct(id) %>% nrow()
})
subtitle_text <- paste0(
  paste0(horizon_labels, ": N=", format(n_by_horizon, big.mark=",")),
  collapse = "   "
)

## 5) 図：棒グラフ＋95%CI（facetで4枠）
p_slopes <- ggplot(slopes_all, aes(x = Group, y = Estimate, fill = Group)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.15, linewidth = 0.6) +
  facet_wrap(~ Horizon, nrow = 1) +
  scale_fill_manual(values = c(
    "nonAKD"       = "#E41A1C",  # 赤
    "Recovery"     = "#1B9E77",  # 緑
    "Non-Recovery" = "#377EB8"   # 青
  )) +
  labs(
    title = "Estimated eGFR Slopes by Group across Time Horizons(Sensitivity Analysis)",
    x = NULL,
    y = expression(Slope~"(mL/min/1.73 m"^2*"/year) with 95% CI")
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p_slopes)

## 6) 表の書き出し（補足資料向け）
readr::write_csv(slopes_all, "Table_Slopes_by_Group_TimeHorizons.csv")

## 7) 図の高解像度保存（論文投稿用）
ggsave("Fig_Slopes_by_Group_TimeHorizons.tiff", plot = p_slopes,
       width = 10.5, height = 3.2, units = "in", dpi = 600, compression = "lzw")

#感度分析：CKDあるなし#####
install.packages("survminer")
library(survival)
library(survminer)
library(dplyr)
##①死亡についてのカプランマイヤー
jin1_Eligibile_akd_ckd <- jin1_Eligibile_unique_id %>%
  mutate(
    group = case_when(
      jin_status == "nonAKD" & CKD_status == "CKD"                          ~ "CKD",
      jin_status == "nonAKD" & (CKD_status %in% c("nonCKD", "nd") | is.na(CKD_status))
      ~ "NKD",
      jin_status == "AKD"    & CKD_status == "CKD"                          ~ "AKD with CKD",
      jin_status == "AKD"    & (CKD_status %in% c("nonCKD", "nd") | is.na(CKD_status))
      ~ "AKD without CKD",
      TRUE ~ NA_character_
    ),
    time_years = as.numeric(last_follow_death - index_plus_210) / 365.25
  ) %>%
  filter(!is.na(group), time_years >= 0) %>%
  mutate(
    group = factor(group,
                   levels = c("NKD", "CKD", "AKD without CKD", "AKD with CKD"))
  )


# モデル作成、グラフ（死亡）
fit_death <- survfit(Surv(time_years, primary_death) ~ group, data = jin1_Eligibile_akd_ckd )

ggsurvplot(fit_death, data = jin1_Eligibile_akd_ckd ,
           fun = "event",                     # 1 - survival（累積死亡率）
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE,
           xlim = c(0, 10),                   # 10年まで表示
           ylim = c(0, 0.4),
           title = "Sensitivity Analysis(death)",
           xlab = "year",
           ylab = "event")

#②eGFRslope 
library(dplyr)
library(ggplot2)
library(nlme)

# 0) 4群の定義（順序も固定）
four_levels <- c("NKD", "CKD", "AKD without CKD", "AKD with CKD")

# 1) 対象抽出：include のみ + 4群へ再分類
jin1_inclusion_4g <- jin1_Eligibile_ver2 %>%
  filter(exclude == "include") %>%
  mutate(
    group = case_when(
      jin_status == "nonAKD" & CKD_status == "CKD"                                 ~ "CKD",
      jin_status == "nonAKD" & (CKD_status %in% c("nonCKD", "nd") | is.na(CKD_status))
      ~ "NKD",
      jin_status == "AKD"    & CKD_status == "CKD"                                 ~ "AKD with CKD",
      jin_status == "AKD"    & (CKD_status %in% c("nonCKD", "nd") | is.na(CKD_status))
      ~ "AKD without CKD",
      TRUE ~ NA_character_
    ),
    group = factor(group, levels = four_levels)
  ) %>%
  filter(!is.na(group))

# 2) years_from_time0 を index_plus_90 基準で作成（0年以降のみ）
akd_time_4 <- jin1_inclusion_4g %>%
  mutate(years_from_time0 = as.numeric(date - index_plus_90) / 365.25) %>%
  filter(years_from_time0 >= 0)

# 3) 規定時点（target_timepoints）の自動設定
time_range <- range(akd_time_4$years_from_time0, na.rm = TRUE)
if (time_range[2] <= 2) {
  target_timepoints <- seq(0, 2, by = 0.25)   # 3ヶ月ごと
} else if (time_range[2] <= 5) {
  target_timepoints <- seq(0, ceiling(time_range[2]), by = 0.5) # 6ヶ月ごと
} else {
  target_timepoints <- seq(0, ceiling(time_range[2]), by = 1)   # 1年ごと
}
message("規定時点: ", paste(target_timepoints, collapse = ", "))

# 代表ウインドウ幅（観測間隔の中央値×0.75）を自動計算
if (!exists("median_interval")) {
  median_interval <- akd_time_4 %>%
    arrange(id, years_from_time0) %>%
    group_by(id) %>%
    summarise(d = diff(years_from_time0), .groups = "drop") %>%
    pull(d) %>%
    median(na.rm = TRUE)
}
window_width <- median_interval * 0.75

# 4) Time Window データ生成関数（そのまま利用）
create_time_window_data_2 <- function(data, targets, window) {
  result_list <- vector("list", length(targets))
  for (i in seq_along(targets)) {
    target <- targets[i]
    window_data <- data %>%
      filter(abs(years_from_time0 - target) <= window) %>%
      group_by(id, group) %>%
      slice_min(abs(years_from_time0 - target), n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      mutate(target_time = target)
    result_list[[i]] <- window_data
  }
  bind_rows(result_list)
}
window_data_4 <- create_time_window_data_2(akd_time_4, target_timepoints, window_width)

# 5) 実測ΔeGFR（time0 は各IDで target_time==0 に最も近い測定を基準）
window_data_obs4 <- window_data_4 %>%
  group_by(id) %>%
  mutate(time0_egfr_obs = egfr[which.min(abs(target_time - 0))]) %>%
  ungroup() %>%
  mutate(observed_diff = egfr - time0_egfr_obs)

# 6) 集団平均と95%CI（n=1のときSEはNAに）
summary_obs_4 <- window_data_obs4 %>%
  group_by(group, target_time) %>%
  summarise(
    mean_diff = mean(observed_diff, na.rm = TRUE),
    sd        = sd(observed_diff, na.rm = TRUE),
    n         = sum(!is.na(observed_diff)),
    se        = ifelse(n > 1, sd / sqrt(n), NA_real_),
    ci_lower  = mean_diff - 1.96 * se,
    ci_upper  = mean_diff + 1.96 * se,
    .groups   = "drop"
  )

# 7) プロット（4群）
cols_4g <- c(
  "NKD"             = "#4E79A7",
  "CKD"             = "#F28E2B",
  "AKD without CKD" = "#59A14F",
  "AKD with CKD"    = "#E15759"
)

ggplot(summary_obs_4, aes(x = target_time, y = mean_diff, color = group)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.05, na.rm = TRUE) +
  scale_color_manual(values = cols_4g, name = "Group") +
  scale_x_continuous(limits = c(0, 10)) +
  labs(
    title = "Observed ΔeGFR from time0 (4 groups by CKD×AKD status)",
    x = "Time from time0 (years)",
    y = expression(Delta~eGFR~"(mL/min/1.73m^2)")
  ) +
  theme_minimal()

# ---（任意）LME も4群版に更新する場合 ---
#  time0_egfr を作る：規定時点0年の egfr（ID内で最初のtarget_time==0の値）
window_data_4 <- window_data_4 %>%
  group_by(id) %>%
  mutate(time0_egfr = egfr[target_time == 0][1]) %>%
  ungroup() %>%
  mutate(egfr_aligned = egfr - time0_egfr)

# 1年以内のLME
fit_window_4_within1year <- lme(
  egfr ~ years_from_time0 * group + time0_egfr - 1,
  random = list(id = pdSymm(form = ~ 1 + years_from_time0)),
  na.action = na.omit,
  data = akd_time_4 %>% filter(years_from_time0 <= 1),
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)

# 全期間のLME
fit_window_4 <- lme(
  egfr ~ years_from_time0 * group + time0_egfr - 1,
  random = list(id = pdSymm(form = ~ 1 + years_from_time0)),
  na.action = na.omit,
  data = akd_time_4,
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)

# 実測値ベースの ΔeGFR プロット（4群）
summary_aligned_4 <- window_data_4 %>%
  group_by(group, target_time) %>%
  summarise(
    mean_diff = mean(egfr_aligned, na.rm = TRUE),
    sd        = sd(egfr_aligned, na.rm = TRUE),
    n         = sum(!is.na(egfr_aligned)),
    se        = ifelse(n > 1, sd / sqrt(n), NA_real_),
    ci_lower  = mean_diff - 1.96 * se,
    ci_upper  = mean_diff + 1.96 * se,
    .groups   = "drop"
  )

ggplot(summary_aligned_4, aes(x = target_time, y = mean_diff, color = group)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.05, na.rm = TRUE) +
  scale_color_manual(values = cols_4g, name = "Group") +
  scale_x_continuous(limits = c(0, 13)) +
  labs(
    title = "ΔeGFR from observed time0 (4 groups, raw egfr)",
    x = "Time from time0 (years)",
    y = expression(Delta~eGFR~"(mL/min/1.73m^2)")
  ) +
  theme_minimal()
