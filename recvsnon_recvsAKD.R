#必要な情報をすべて含有するcsvファイル(jin1_Eligible? primay_ESKDやdeath情報や感度分析に必要な情報などすべて含有)を作成し、適当なdfを作成
View(jin1_Eligibile)
library(readr)
setwd("E:/R")
jin1_Eligibile <- read_csv("jin1_Eligibile.csv", locale = locale(encoding = "SHIFT-JIS"))
colnames(jin1_Eligibile)
library(dplyr)
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
{
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

#Time Window Approachの実装{
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
#slope描出に必要なcodeのみ
#slope描出に必要なcodeのみ####
# 線形混合モデル　1年####
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

p_death <- ggsurvplot(
  fit_death,
  data        = jin1_Eligibile_death,
  fun         = "event",
  pval        = TRUE,
  conf.int    = TRUE,
  risk.table  = TRUE,
  xlim        = c(0, 9.5),
  ylim        = c(0, 0.4),
  title       = "Cumulative incidence of all-cause death",
  xlab        = "Years",
  ylab        = "Cumulative incidence",
  break.time.by = 2,
  
  # ★ここを追加（凡例ラベルをシンプルに）
  legend.labs = c("nonAKD", "Recovery", "Non-Recovery"),
  
  # 任意（凡例タイトル変更）
  legend.title = "Group"
)


colnames(jin1_Eligibile)

##論文用の図（デーブル付き）####
# 2カラム用の例：幅 180 mm ≒ 7.1 inch
setwd("E:/R")
tiff(
  filename = "Figure2_death_KM.tiff",
  width    = 7.1,   # inch
  height   = 6.0,   # inch
  units    = "in",
  res      = 600,
  compression = "lzw",
  type     = "cairo"   # ← ここがポイント
)

print(p_death)         # カーブ＋リスクテーブルをまとめて出力

dev.off()

## ---- 図の作成（曲線のみ、抄録仕様） ----
km_plot <- ggsurvplot(
  fit_death,
  data = jin1_Eligibile_death,
  fun = "event",
  pval = TRUE,
  conf.int = TRUE,
  conf.int.alpha = 0.20,          # CIを薄く
  size = 1.1,                     # 線を少し太く
  risk.table = FALSE,             # 抄録では非表示推奨
  xlim = c(0, 9.5),
  ylim = c(0, 0.4),
  title = "Cumulative incidence of death by AKD recovery status",
  xlab = "Years since index date",
  ylab = "Cumulative incidence",
  legend.title = NULL,
  legend.labs = c("non-AKD", "AKD Recovery", "AKD Non-Recovery"),
  palette = c("#E64B35", "#00A087", "#4DBBD5"),
  ggtheme = theme_minimal(base_size = 12)
)

km_plot$plot <- km_plot$plot +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 10)
  )

# （任意）p値の位置を微調整したい場合
# km_plot$plot <- km_plot$plot + annotate("text", x = 0.5, y = 0.36, label = km_plot$pval, hjust = 0)

# ---- JPEGで保存（640×480px, ≤1.5MB）----
jpeg("death_curve_for_abstract.jpeg", width = 640, height = 480, units = "px", quality = 95)
print(km_plot$plot)
dev.off()

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

## ---- 図の作成（曲線のみ、抄録仕様） ----
km_plot <- ggsurvplot(
  fit_composite,
  data = jin1_Eligibile_composite,
  fun = "event",
  pval = TRUE,
  conf.int = TRUE,
  conf.int.alpha = 0.20,          # CIを薄く
  size = 1.1,                     # 線を少し太く
  risk.table = FALSE,             # 抄録では非表示推奨
  xlim = c(0, 9.5),
  ylim = c(0, 0.4),
  title = "Cumulative incidence of composite_event by AKD recovery status",
  xlab = "Years since index date",
  ylab = "Cumulative incidence",
  legend.title = NULL,
  legend.labs = c("non-AKD", "AKD Recovery", "AKD Non-Recovery"),
  palette = c("#E64B35", "#00A087", "#4DBBD5"),
  ggtheme = theme_minimal(base_size = 12)
)

km_plot$plot <- km_plot$plot +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 10)
  )

# （任意）p値の位置を微調整したい場合
# km_plot$plot <- km_plot$plot + annotate("text", x = 0.5, y = 0.36, label = km_plot$pval, hjust = 0)

# ---- JPEGで保存（640×480px, ≤1.5MB）----
jpeg("composite_event_curve_for_abstract.jpeg", width = 640, height = 480, units = "px", quality = 95)
print(km_plot$plot)
dev.off()


#cox比例ハザード#####
{
# ---- 必要パッケージ ----
library(dplyr)
library(survival)
library(broom)
library(ggplot2)

# -------------------------------
# 1) 解析用データ（3群）の作成（CKD_statusの基準を"nonCKD"に固定）
# -------------------------------
jin1_Eligibile_unique_id <- jin1_Eligibile %>%
  filter(exclude == "include", jin_status %in% c("AKD", "nonAKD")) %>%
  group_by(id) %>%
  arrange(index_date, date, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

# 確認：unique id数（任意）
jin1_Eligibile_unique_id %>%
  group_by(jin_status) %>%
  summarise(n_unique_ids = n_distinct(id), .groups = "drop")

# 3群 + 共変量整形（CKD_status を factor にし、"nonCKD" を参照カテゴリへ）
jin1_Eligibile_cox_3group <- jin1_Eligibile_unique_id %>%
  filter(exclude == "include", jin_status %in% c("AKD", "nonAKD")) %>%
  distinct(id, .keep_all = TRUE) %>%
  mutate(
    jin_label = case_when(
      jin_status == "nonAKD" ~ "nonAKD",
      jin_status == "AKD" & `150_210recovery` == 1 ~ "Recovery",
      jin_status == "AKD" & `150_210recovery` == 2 ~ "Non-Recovery",
      jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 1 ~ "Recovery",
      jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` %in% c(0, 2) ~ "Non-Recovery",
      TRUE ~ NA_character_
    ),
    jin_label    = factor(jin_label, levels = c("nonAKD", "Recovery", "Non-Recovery")),
    arb_acei_use = if_else(coalesce(arb, 0) == 1 | coalesce(acei, 0) == 1, 1L, 0L),
    time_years   = as.numeric(last_follow_death - index_plus_210) / 365.25,
    
    # ★ ここがポイント：CKD_status を文字列→factor化し参照を "nonCKD" に固定
    CKD_status = case_when(
      as.character(CKD_status) %in% c("CKD","1")     ~ "CKD",
      as.character(CKD_status) %in% c("nonCKD","0") ~ "nonCKD",
      TRUE ~ as.character(CKD_status)
    ),
    CKD_status = factor(CKD_status, levels = c("nonCKD","CKD"))
  ) %>%
  filter(!is.na(jin_label), !is.na(time_years), time_years >= 0, !is.na(CKD_status))

# 確認：群別症例数（任意）
jin1_Eligibile_cox_3group %>%
  group_by(jin_label) %>%
  summarise(n_unique_ids = n_distinct(id), .groups = "drop")

# -------------------------------
# 2) Cox比例ハザード（3群）— CKD_statusは「CKD vs nonCKD」で推定
# -------------------------------
cox_model_3group <- coxph(
  Surv(time_years, primary_death) ~
    jin_label + age + index_cre + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15 +
    CKD_status,   # ← 参照=nonCKD なので、出力係数は CKD vs nonCKD
  data = jin1_Eligibile_cox_3group
)

# 係数の95%CI（HRスケール）
exp(confint(cox_model_3group))

coef(cox_model_3group)

# AKDrecovery vs non-recoveryを比較
library(multcomp)

summary(
  glht(
    cox_model_3group,
    linfct = c("`jin_labelNon-Recovery` - jin_labelRecovery = 0")
  )
)



#AKDの比較とテーブル作成
{
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(Cairo)
  library(broom)
  library(dplyr)
  
  # ---- 1) HRテーブル（上で作ったものと同じ）----
  hr_table <- tidy(
    cox_model_3group,
    exponentiate = TRUE,
    conf.int    = TRUE
  ) %>%
    filter(term %in% c("jin_labelRecovery", "jin_labelNon-Recovery")) %>%
    mutate(
      Contrast = dplyr::recode(
        term,
        "jin_labelRecovery"     = "AKD Recovery vs nonAKD",
        "jin_labelNon-Recovery" = "AKD Non-Recovery vs nonAKD"
      ),
      HR        = round(estimate, 2),
      CI_lower  = round(conf.low, 2),
      CI_upper  = round(conf.high, 2),
      p_raw     = p.value,
      `p-value` = ifelse(p_raw < 0.01,
                         "<0.01",
                         sprintf("%.2f", p_raw))
    ) %>%
    dplyr::select(
      Contrast,
      HR,
      `Lower 95% CI` = CI_lower,
      `Upper 95% CI` = CI_upper,
      `p-value`
    )
  
  # ---- 2) フォレストプロット用データ・図 ----
  hr_plot_df <- tidy(
    cox_model_3group,
    exponentiate = TRUE,
    conf.int    = TRUE
  ) %>%
    filter(term %in% c("jin_labelRecovery", "jin_labelNon-Recovery")) %>%
    mutate(
      Contrast = dplyr::recode(
        term,
        "jin_labelRecovery"     = "AKD Recovery vs nonAKD",
        "jin_labelNon-Recovery" = "AKD Non-Recovery vs nonAKD"
      )
    )
  
  hr_plot_df$Contrast <- factor(
    hr_plot_df$Contrast,
    levels = c("AKD Non-Recovery vs nonAKD", "AKD Recovery vs nonAKD")
  )
  
  p_forest <- ggplot(hr_plot_df, aes(x = Contrast, y = estimate)) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                  width = 0.1, linewidth = 0.7) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    coord_flip() +
    scale_y_log10(
      name   = "Hazard ratio for all-cause mortality (log scale)",
      breaks = c(0.5, 1, 2, 3, 4, 7),
      limits = c(0.5, 7)
    ) +
    labs(
      x = NULL,
      title = "Adjusted hazard ratios for all-cause mortality\n(AKD groups vs nonAKD)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.y      = element_text(size = 10),
      axis.text.x      = element_text(size = 9),
      plot.title       = element_text(hjust = 0.5, face = "bold"),
      plot.margin      = margin(t = 5.5, r = 10, b = 5.5, l = 7, unit = "pt")
    )
  
  # ---- 3) テーブルを tableGrob に変換 ----
  table_theme <- gridExtra::ttheme_minimal(
    core = list(
      fontsize = 9,
      padding  = unit(c(2, 2), "mm")
    ),
    colhead = list(
      fontsize = 9,
      fontface = "bold",
      padding  = unit(c(2, 2), "mm")
    )
  )
  
  table_grob <- gridExtra::tableGrob(
    hr_table,
    rows  = NULL,
    theme = table_theme
  )
  
  # ---- 4) テーブル＋フォレストを縦に並べる ----
  combined <- gridExtra::arrangeGrob(
    table_grob,
    p_forest,
    ncol    = 1,
    heights = c(1, 2)   # テーブル：図 = 1:2 の高さ
  )
  setwd("E:/R")  
  # ---- 5) 一枚の図としてTIFF保存（横切れ対策で width 少し広め）----
  CairoTIFF(
    filename = "Figure3_forest_with_table_death_akd_3groups.tif",
    width  = 160,   # mm（必要なら 170, 180 などに調整）
    height = 140,   # mm
    units  = "mm",
    res    = 600
  )
  grid::grid.newpage()
  grid::grid.draw(combined)
  dev.off()
  
}

# （任意）PH仮定チェック
# print(cox.zph(cox_model_3group))

# -------------------------------
# 3) Recovery vs Non-Recovery の直接比較（任意：変更なし）
# -------------------------------
 library(multcomp)
 fit <- cox_model_3group
 g <- multcomp::glht(fit, linfct = c("`jin_labelRecovery` - `jin_labelNon-Recovery` = 0"))
 summary(g); ci <- confint(g)
 est  <- summary(g)$test$coef[1]; pval <- summary(g)$test$pvalues[1]
 HR   <- exp(est); HR_CI <- exp(ci$confint[1, c("lwr","upr")])
 c(HR = HR, CI_low = HR_CI[1], CI_high = HR_CI[2], p = pval)

# -------------------------------
# 4) tidy & フォレスト（CKDを "CKD vs nonCKD" と明示）
# -------------------------------
tidy3 <- broom::tidy(cox_model_3group, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    term_raw  = term,
    term_nice = dplyr::recode(
      term,
      "jin_labelRecovery"     = "AKD Recovery vs nonAKD",
      "jin_labelNon-Recovery" = "AKD Non-Recovery vs nonAKD",
      "age"                   = "Age (per 1y)",
      "index_cre"             = "Index Creatinine (per 1 mg/dL)",
      "arb_acei_use"          = "ARB or ACEi Use",
      # Charlson等のラベル
      "dn1"  = "CHF",           "dn3"  = "Rheumatologic", "dn4"  = "Malignancy",
      "dn5"  = "Liver Disease", "dn6"  = "Peptic Ulcer",  "dn7"  = "MI",
      "dn8"  = "Renal Disease", "dn9"  = "Metastatic Cancer",
      "dn10" = "Diabetes",      "dn12" = "Stroke",        "dn13" = "Hemiplegia",
      "dn14" = "PVD",           "dn15" = "COPD",
      
      # ★ CKDの表示名を明示：係数名は "CKD_statusCKD" になる想定
      "CKD_statusCKD"         = "CKD vs nonCKD",
      
      .default = term
    )
  )

# モデル式の順に並べる（上から：式先頭→末尾）
form_order <- attr(stats::terms(cox_model_3group), "term.labels")
order_raw <- unlist(lapply(form_order, function(x) tidy3$term_raw[startsWith(tidy3$term_raw, x)]))
order_raw <- order_raw[order_raw %in% tidy3$term_raw]
display_order <- tidy3 %>%
  filter(term_raw %in% order_raw) %>%
  arrange(match(term_raw, order_raw)) %>%
  pull(term_nice)

tidy3 <- tidy3 %>%
  mutate(term_nice = factor(term_nice, levels = rev(display_order)))

# 表（任意で出力）
tidy3_table <- tidy3 %>%
  dplyr::select(term = term_nice, HR = estimate, `CI low` = conf.low, `CI high` = conf.high, p.value) %>%
  arrange(desc(term))

# -------------------------------
# 5) フォレストプロット作成＆論文用に保存（TIFF/PDFの2形式）
# -------------------------------
p_forest <- ggplot(tidy3, aes(x = estimate, y = term_nice)) +
  geom_point(size = 2.8) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.18, linewidth = 0.6) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  labs(
    x = "Hazard Ratio (log scale)",
    y = NULL,
    title = "All-cause death: Adjusted hazard ratios (3-group model)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0, face = "bold"),
    axis.text.y = element_text(size = 9)
  )

# 画面確認
print(p_forest)

# 保存：雑誌で求められやすいTIFF(600 dpi)とPDF（埋め込みフォント）を作成
# ※ Windows などでアンチエイリアス/フォント埋め込みを安定させるために Cairo を使用
#   （必要なら install.packages("Cairo")）
#   出力サイズは投稿規定に合わせ調整（例：幅 7 inch × 高さ 6 inch）

# TIFF
ggsave(
  filename = "forest_cox3_death.tiff",
  plot = p_forest,
  device = "tiff",
  dpi = 600,
  width = 7, height = 6, units = "in",
  compression = "lzw"
)

# PDF（フォント埋め込み用に cairo_pdf を明示）
ggsave(
  filename = "forest_cox3_death.pdf",
  plot = p_forest,
  device = cairo_pdf,
  width = 7, height = 6, units = "in"
)
}

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
library(survival)
library(survminer)
library(ggplot2)

fit_death_sens <- survfit(Surv(time_years, primary_death) ~ group,
                          data = jin1_Eligibile_sens)

p_death_sens <- ggsurvplot(
  fit_death_sens,
  data        = jin1_Eligibile_sens,
  fun         = "event",                     # 1 - survival（累積死亡率）
  pval        = TRUE,
  conf.int    = TRUE,
  risk.table  = TRUE,
  xlim        = c(0, 9.5),                    # 10年まで表示
  ylim        = c(0, 0.4),
  title       = "Cumulative incidence of all-cause death\n(Sensitivity analysis)",
  xlab        = "Years",
  ylab        = "Cumulative incidence",
  break.time.by = 2,
  legend.title  = "Group",
  legend.labs   = c("nonAKD", "Recovery", "Non-Recovery")  # 必要に応じて
)
# 必要なら作業ディレクトリを指定
setwd("E:/R")
tiff(
  filename   = "Figure5a_death_KM_sensitivity.tiff",
  width      = 7.1,   # inch（約180 mm）
  height     = 6.0,   # inch（適宜調整）
  units      = "in",
  res        = 600,
  compression = "lzw",
  type       = "cairo"   # ← これが CairoTIFF の代わり
)

print(p_death_sens)      # カーブ + リスクテーブルをまとめて出力

dev.off()
