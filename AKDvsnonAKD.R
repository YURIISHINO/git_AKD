#non_AKDvsAKDにする8/8 ####
library(dplyr)
library(lubridate)
setwd("E:/R")
#データ整形####
####除外要件を0/1でラベルし除外
# Eligible patients
# 浜松医科大学にて18歳以上かつ2013年に一回以上血清クレアチニン測定をしている患者
jin1_index_cre_egfr
#### follow
jin1_index_follow_ver2 <- jin1_index_cre_egfr %>%
  group_by(id) %>%
  mutate(no_follow = if_else(any(date > index_date & date <= index_date + years(1)), 0, 1)) %>%
  ungroup()

### cre < 0.1 または egfr > 300 を除外 ###
#jin1_index_follow_ver2$cre_egfr_over <- ifelse(jin1_index_follow_ver2$cre < 0.1 | jin1_index_cre_egfr$egfr > 300, 1, 0)
#回復判定までの範囲（index_date+210日までに制限）
jin1_index_follow_ver2 <- jin1_index_follow_ver2 %>%
  mutate(
    cre_egfr_over = if_else(
      (cre < 0.1 | egfr > 300) & date <= index_date + 90,
      1,
      0
    )
  )

### 維持透析 ###
hd_summary_index <- hd_summary %>%
  left_join(jin1_index_cre_egfr %>% dplyr::select(id, index_date) %>% distinct(), by = "id")

hd_summary_delete <- hd_summary_index %>%
  filter(
    !is.na(index_date) & (
      start < index_date |
        (start >= index_date & start <= index_date + 90)
    )
  )
jin1_index_follow_ver2$hd <- ifelse(jin1_index_follow_ver2$id %in% hd_summary_delete$id, 1, 0)

### 腎移植 ###
rr_index <- rr %>%
  left_join(jin1_index_cre_egfr %>% dplyr::select(id, index_date) %>% distinct(), by = "id")

rr_delete <- rr_index %>%
  filter(
    !is.na(index_date) & (
      rr_date < index_date | 
        (rr_date >= index_date & rr_date <= index_date + 90)
    )
  )

jin1_index_follow_ver2$rr <- ifelse(jin1_index_follow_ver2$id %in% rr_delete$id, 1, 0)

### 死亡除外 ###
death_index <- death %>%
  left_join(jin1_index_cre_egfr %>% dplyr::select(id, index_date) %>% distinct(), by = "id")

death_delete <- death_index %>%
  filter(
    !is.na(index_date) & (
      death_date < index_date | 
        (death_date >= index_date & death_date <= index_date + 90)
    )
  )

jin1_index_follow_ver2$death <- ifelse(jin1_index_follow_ver2$id %in% death_delete$id, 1, 0)

### 除外ラベルの付与（優先順位に基づく） ###
exclude_priority <- c("hd", "rr", "death", "no_follow", "cre_egfr_over", "include")

jin1_exclude_by_id <- jin1_index_follow_ver2 %>%
  mutate(
    exclude = case_when(
      hd == 1 ~ "hd",
      rr == 1 ~ "rr",
      death == 1 ~ "death",
      no_follow == 1 ~ "no_follow",
      cre_egfr_over == 1 ~ "cre_egfr_over",
      TRUE ~ "include"
    )
  ) %>%
  group_by(id) %>%
  summarise(
    exclude = factor(
      exclude_priority[
        min(match(exclude, exclude_priority), na.rm = TRUE)
      ],
      levels = exclude_priority
    )
  ) %>%
  ungroup()

# 除外ステータスごとの一意のID数をカウント
exclude_summary <- jin1_exclude_by_id %>%
  group_by(exclude) %>%
  summarise(unique_ids = n_distinct(id)) %>%
  arrange(desc(unique_ids))
print(exclude_summary)

jin1_index_follow_ver2 <- jin1_index_follow_ver2 %>%
  left_join(jin1_exclude_by_id, by = "id")
View(jin1_index_follow_ver2)

#一意のidとAKDstatus/AKIstatus/CKDstatusの結果を抜き出す
id_status_unique <- jin1_status %>%
  dplyr::select(id, AKD_status,AKI_status,CKD_status) %>%
  distinct()
print(id_status_unique)

# すべて結合
jin1_Eligibile_ver2 <- jin1_index_follow_ver2 %>%
  left_join(jin1_med, by = "id") %>%
  left_join(dn_complete, by = "id") %>%
  left_join(id_status_unique, by = "id")

#　index_date時の年齢を追加
jin1_Eligibile_ver2 <- jin1_Eligibile_ver2 %>%
  mutate(age = floor(interval(birth, index_date) / years(1)))

library(dplyr)
library(lubridate)
setwd("E:/R")
#AKD/nonAKD/non_AKDを定義(変更点)####
#AKD:AKDstatusが"AKD_cre""AKD_egfr"もしくはAKI_statusが"AKI"
#non_AKD:AKD_statusが"nd"もしくは"nonAKD"　かつAKI_statusが "nd"もしくは"nonAKI" かつCKD_statusが"nd"かつ"nonCKD"
#上記以外がother
# ステップ1：jin_status を仮で付与（各行に付ける）
jin1_Eligibile_ver2 <- jin1_Eligibile_ver2 %>%
  mutate(
    jin_status_tmp = case_when(
      AKD_status %in% c("AKD_cre", "AKD_egfr") | AKI_status == "AKI" ~ "AKD",
      AKD_status %in% c("nd", "nonAKD") &
        AKI_status %in% c("nd", "nonAKI") ~ "nonAKD",
      TRUE ~ "other"
    ),
    jin_status_tmp = factor(jin_status_tmp, levels = c("AKD", "nonAKD", "other"))
  )


# ステップ2：idごとの優先順位つき代表ステータス取得
id_status_priority <- jin1_Eligibile_ver2 %>%
  group_by(id) %>%
  summarise(jin_status_rep = first(jin_status_tmp[order(jin_status_tmp)]), .groups = "drop")

# ステップ3：元データに結合して因子化
jin1_Eligibile_ver2 <- jin1_Eligibile_ver2 %>%
  left_join(id_status_priority, by = "id") %>%
  dplyr::select(-jin_status_tmp) %>%
  mutate(jin_status = factor(jin_status_rep, levels = c("AKD", "nonAKD", "other"))) %>%
  dplyr::select(-jin_status_rep)

jin1_Eligibile_ver2 %>%
  group_by(jin_status) %>%
  summarise(n_unique_ids = n_distinct(id), .groups = "drop")


#アウトカム関連まとめ
#death,ESKD(hd_new,rr)日にちをつける
death
hd_new
rr
jin1_Eligibile_ver2 <- jin1_Eligibile_ver2 %>%
  left_join(death, by = "id") %>%
  left_join(hd_new, by = "id") %>%
  left_join(rr, by = "id")
#hd_new,rrはESKD_dateにまとめ
jin1_Eligibile_ver2 <- jin1_Eligibile_ver2 %>%
  mutate(
    ESKD_date = pmin(hd_new_start, rr_date, na.rm = TRUE)
  )
#最後のデータが存在する日をlast_dataとする
jin1_Eligibile_ver2 <- jin1_Eligibile_ver2 %>%
  group_by(id) %>%
  mutate(
    last_data = max(date, na.rm = TRUE)
  ) %>%
  ungroup()

#フォロー開始日、アウトカムの0.1データ
library(dplyr)
library(lubridate)
jin1_Eligibile_ver2 <- jin1_Eligibile_ver2 %>%
  mutate(
    # 90日後の基準日を作成
    index_plus_90 = index_date + 90,
    
    # 死亡とESKDの発症が90日以降かどうかを判定
    primary_death = if_else(!is.na(death_date) & death_date > index_plus_90, 1, 0),
    primary_ESKD = if_else(!is.na(ESKD_date) & ESKD_date > index_plus_90, 1, 0)
  )

jin1_Eligibile_ver2 <- jin1_Eligibile_ver2 %>%
  mutate(
    last_follow_death = if_else(primary_death == 1, as.Date(death_date), as.Date(last_data))
  )

jin1_Eligibile_ver2 <- jin1_Eligibile_ver2 %>%
  mutate(
    # composite endpoint が起きたか判定（primary_death または primary_ESKD のどちらかが1）
    primary_composite_event = (primary_death == 1 | primary_ESKD == 1),
    # event_date は death_date と ESKD_date の早い方（Date型で）
    event_date = pmin(as.Date(death_date), as.Date(ESKD_date), na.rm = TRUE),
    # composite_event が TRUE（=1）なら event_date、それ以外は last_data
    last_follow_composite = if_else(primary_composite_event == 1, event_date, as.Date(last_data))
  )

#time0=index_plus_90 からの年差を計算
jin1_Eligibile_ver2 <- jin1_Eligibile_ver2 %>%
  mutate(
    days_from_time0 = as.numeric(difftime(date, index_plus_90, units = "days")),
    years_from_time0 = round(days_from_time0 / 365.25, 3)
  )

#AKD,AKI発生日を入れる
# CSVファイルの読み込み
library(readr)
library(lubridate)
jin1_AKI_date_nonNA_unique <- read_csv(
  "jin1_AKI_date_nonNA_unique.csv",
  col_types = cols(
    id = col_double(),
    AKI_date = col_date(format = "%Y-%m-%d")  # Date型として読み込み
  )
)
jin1_AKD_date_nonNA_unique <- read_csv(
  "jin1_AKD_date_nonNA_unique.csv",
  col_types = cols(
    id = col_double(),
    AKD_date = col_date(format = "%Y-%m-%d")  # Date型として読み込み
  )
)


jin1_Eligibile_ver2 <- jin1_Eligibile_ver2 %>%
  left_join(jin1_AKI_date_nonNA_unique, by = "id") %>%
  left_join(jin1_AKD_date_nonNA_unique, by = "id")
colnames(jin1_Eligibile_ver2 )

library(dplyr)
library(lubridate)

#time0=index_plus_90
egfr_time0 <- jin1_Eligibile_ver2 %>%
  # 期間境界の準備
  mutate(
    upper_90_180 = index_plus_90 + days(180),
    # AKI_date / AKD_date のうち早い方。両方NAなら index_date を使う
    lower_pre90 = case_when(
      !is.na(AKI_date) & !is.na(AKD_date) ~ pmin(AKI_date, AKD_date),
      !is.na(AKI_date)                    ~ AKI_date,
      !is.na(AKD_date)                    ~ AKD_date,
      TRUE                                ~ index_date
    )
  ) %>%
  group_by(id, index_date, index_plus_90, upper_90_180, lower_pre90) %>%
  summarise(
    # ① index_plus_90以降180日以内の最大eGFR
    egfr_90_180 = max(egfr[date >= index_plus_90 & date <= upper_90_180], na.rm = TRUE),
    # ② ①が無い場合の代替：AKI/AKDの早い方～index_plus_90
    egfr_pre90  = max(egfr[date >= lower_pre90 & date <= index_plus_90], na.rm = TRUE),
    # ③（保険）index_date～index_plus_90 も計算しておく（両日不明のとき用）
    egfr_0_90   = max(egfr[date >= index_date & date <= index_plus_90], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # max(..., na.rm=TRUE) で該当が無いと -Inf になるので処理
  mutate(
    egfr_90_180 = ifelse(is.finite(egfr_90_180), egfr_90_180, NA_real_),
    egfr_pre90  = ifelse(is.finite(egfr_pre90),  egfr_pre90,  NA_real_),
    egfr_0_90   = ifelse(is.finite(egfr_0_90),   egfr_0_90,   NA_real_),
    time0_egfr  = coalesce(egfr_90_180, egfr_pre90, egfr_0_90)
  ) %>%
  dplyr::select(id, index_date, time0_egfr)

# ステップ2: 元のデータに統合
jin1_Eligibile_ver2 <- jin1_Eligibile_ver2 %>%
  left_join(egfr_time0, by = c("id", "index_date"))
colnames(jin1_Eligibile_ver2)
View(jin1_Eligibile_ver2)

#jin1_EligibileをCSVファイルに書き出し
getwd()  # 現在の作業ディレクトリを確認
write.csv(jin1_Eligibile_ver2, file = "E:/R/jin1_Eligibile_ver2.csv", row.names = FALSE)

#eGFR・線形混合効果モデルを走らせる、また解析に必要なcodeのみ####
library(dplyr)
library(ggplot2)
library(nlme)

# 1) 対象抽出：AKD と nonAKD のみ
jin1_inclusion_ver2 <- jin1_Eligibile_ver2 %>%
  filter(exclude == "include", jin_status %in% c("AKD", "nonAKD")) %>%
  mutate(
    # 因子順を固定（凡例順・モデル基準）
    jin_status = factor(jin_status, levels = c("AKD", "nonAKD"))
  )

library(readr)
write_excel_csv(jin1_inclusion_ver2, "jin1_inclusion_ver2.csv", na = "")
#CSVファイルを読み込む
jin1_inclusion_ver2 <- read_csv("jin1_inclusion_ver2.csv")

# 2) years_from_time0 を index_plus_90 基準で作成（0年以降のみ）
akd_time_2 <- jin1_inclusion_ver2 %>%
  mutate(
    years_from_time0 = as.numeric(date - index_plus_90) / 365.25
  ) %>%
  filter(years_from_time0 >= 0)

# 3) 規定時点の設定
time_range <- range(akd_time_2$years_from_time0, na.rm = TRUE)
if (time_range[2] <= 2) {
  target_timepoints <- seq(0, 2, by = 0.25)   # 3ヶ月ごと
} else if (time_range[2] <= 5) {
  target_timepoints <- seq(0, ceiling(time_range[2]), by = 0.5) # 6ヶ月ごと
} else {
  target_timepoints <- seq(0, ceiling(time_range[2]), by = 1)   # 1年ごと
}
message("規定時点: ", paste(target_timepoints, collapse = ", "))

# ※ median_interval が未定義なら、代表値を自動計算（観測間隔の年単位の中央値）
if (!exists("median_interval")) {
  median_interval <- akd_time_2 %>%
    arrange(id, years_from_time0) %>%
    group_by(id) %>%
    summarise(d = diff(years_from_time0), .groups = "drop") %>%
    pull(d) %>%
    median(na.rm = TRUE)
}
window_width <- median_interval * 0.75

# 4) Time Window データ生成関数
create_time_window_data_2 <- function(data, targets, window) {
  result_list <- vector("list", length(targets))
  for (i in seq_along(targets)) {
    target <- targets[i]
    window_data <- data %>%
      filter(abs(years_from_time0 - target) <= window) %>%
      group_by(id, jin_status) %>%
      slice_min(abs(years_from_time0 - target), n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      mutate(target_time = target)
    result_list[[i]] <- window_data
  }
  bind_rows(result_list)
}
window_data_2 <- create_time_window_data_2(akd_time_2, target_timepoints, window_width)

# 5) 実測ΔeGFR（time0 が欠けるIDに備え、最も0に近い測定を基準にする）
window_data_obs2 <- window_data_2 %>%
  group_by(id) %>%
  mutate(time0_egfr_obs = egfr[which.min(abs(target_time - 0))]) %>%
  ungroup() %>%
  mutate(observed_diff = egfr - time0_egfr_obs)

# 6) 集団平均と95%CI
summary_obs_aligned <- window_data_obs2 %>%
  group_by(jin_status, target_time) %>%
  summarise(
    mean_diff = mean(observed_diff, na.rm = TRUE),
    sd       = sd(observed_diff, na.rm = TRUE),
    n        = sum(!is.na(observed_diff)),
    se       = sd / sqrt(n),
    ci_lower = mean_diff - 1.96 * se,
    ci_upper = mean_diff + 1.96 * se,
    .groups  = "drop"
  )

# 7) プロット（AKD と nonAKD の2本が出る）
ggplot(summary_obs_aligned, aes(x = target_time, y = mean_diff, color = jin_status)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.05) +
  scale_x_continuous(limits = c(0, 10)) +
  labs(
    title = "Observed ΔeGFR from time0 (AKD vs nonAKD)",
    x = "Time from time0 (years)",
    y = expression(Delta~eGFR~"(mL/min/1.73m^2)"),
    color = "Group"
  ) +
  theme_minimal()

# --- 参考：1年以内のLME（データ名と因子を修正） ---
fit_window_2_within1year <- lme(
  egfr ~ years_from_time0 * jin_status + time0_egfr - 1,
  random = list(id = pdSymm(form = ~ 1 + years_from_time0)),
  na.action = na.omit,
  data = akd_time_2 %>% filter(years_from_time0 <= 1),
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)

# 線形混合モデル　1年####
fit_window_2_within_1year <- lme(
  egfr ~ years_from_time0 * jin_status + time0_egfr - 1,
  random = list(id = pdSymm(form = ~ 1 + years_from_time0)),
  na.action = na.omit,
  data = filter(akd_time_2, years_from_time0 <= 1),
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)

# 線形混合モデル 全期間####
fit_window_2 <- lme(
  egfr ~ years_from_time0 * jin_status + time0_egfr - 1,
  random = list(id = pdSymm(form = ~ 1 + years_from_time0)),
  na.action = na.omit,
  data = akd_time_2,
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)

# 実測値でΔeGFR（time0との差）を計算
window_data_2 <- window_data_2 %>%
  group_by(id) %>%
  mutate(egfr_time0 = egfr[target_time == 0][1]) %>%
  ungroup() %>%
  mutate(egfr_aligned = egfr - egfr_time0)

# 集計：各グループ・時間ごとに平均・CIを計算
summary_aligned_2 <- window_data_2 %>%
  group_by(jin_status, target_time) %>%
  summarise(
    mean_diff = mean(egfr_aligned, na.rm = TRUE),
    sd = sd(egfr_aligned, na.rm = TRUE),
    n = sum(!is.na(egfr_aligned)),
    se = sd / sqrt(n),
    ci_lower = mean_diff - 1.96 * se,
    ci_upper = mean_diff + 1.96 * se,
    .groups = "drop"
  )

# 実測値ベースの ΔeGFR プロット
ggplot(summary_aligned_2, aes(x = target_time, y = mean_diff, color = jin_status)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.05) +
  scale_x_continuous(limits = c(0, 13)) +
  labs(
    title = "ΔeGFR from observed time0 (based on raw egfr)",
    x = "Time from time0 (years)",
    y = expression(Delta~eGFR~"(mL/min/1.73m2)")
  ) +
  theme_minimal()

#4つのモデルの固定効果について下向き棒グラフ#####
# 必要パッケージ
library(dplyr)
library(ggplot2)
library(multcomp)
library(tibble)
library(forcats)
library(nlme)

# 0) 前処理：jin_statusを明示的に2水準に固定（順序も指定）
akd_time_2 <- akd_time_2 %>%
  mutate(
    jin_status = factor(jin_status, levels = c("nonAKD", "AKD"))
  )

# 1) 4つの期間で同一仕様のLMEをフィット
#    egfr ~ years_from_time0 * jin_status + time0_egfr - 1
#    ランダム効果: ~ 1 + years_from_time0 | id（対称分散共分散）
fit_window_2_within_1year <- lme(
  egfr ~ years_from_time0 * jin_status + time0_egfr - 1,
  random = list(id = pdSymm(form = ~ 1 + years_from_time0)),
  na.action = na.omit,
  data = dplyr::filter(akd_time_2, years_from_time0 <= 1),
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)
summary(fit_window_2_within_1year)#p=0.0179で有意差はある


fit_window_2_within_2year <- lme(
  egfr ~ years_from_time0 * jin_status + time0_egfr - 1,
  random = list(id = pdSymm(form = ~ 1 + years_from_time0)),
  na.action = na.omit,
  data = dplyr::filter(akd_time_2, years_from_time0 <= 2),
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)
summary(fit_window_2_within_2year)#p=0.9618で有意差なし

fit_window_2_within_3year <- lme(
  egfr ~ years_from_time0 * jin_status + time0_egfr - 1,
  random = list(id = pdSymm(form = ~ 1 + years_from_time0)),
  na.action = na.omit,
  data = dplyr::filter(akd_time_2, years_from_time0 <= 3),
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)

fit_window_2 <- lme(  # 全期間
  egfr ~ years_from_time0 * jin_status + time0_egfr - 1,
  random = list(id = pdSymm(form = ~ 1 + years_from_time0)),
  na.action = na.omit,
  data = akd_time_2,
  method = "REML",
  control = lmeControl(maxIter = 1e8, msMaxIter = 1e8,
                       opt = "optim", optimMethod = "L-BFGS-B")
)

# 2) スロープ抽出ヘルパー（AKD / nonAKD）
#    傾き = years_from_time0 + years_from_time0:jin_status<level>
extract_slope_df_2groups <- function(fit, horizon_label){
  cn <- names(fixef(fit))
  p  <- length(cn)
  
  # years_from_time0 の係数位置
  idx_main <- which(cn == "years_from_time0")
  
  # 相互作用の係数位置（環境差異による命名の揺れに少し耐性を持たせる）
  idx_int_nonAKD <- which(cn %in% c("years_from_time0:jin_statusnonAKD",
                                    "jin_statusnonAKD:years_from_time0"))
  idx_int_AKD    <- which(cn %in% c("years_from_time0:jin_statusAKD",
                                    "jin_statusAKD:years_from_time0"))
  
  # ベクトル定義
  v_nonAKD <- rep(0, p)
  v_AKD    <- rep(0, p)
  
  if (length(idx_main) == 1) {
    v_nonAKD[idx_main] <- 1
    v_AKD[idx_main]    <- 1
  } else {
    stop("Could not find 'years_from_time0' in fixed effects.")
  }
  if (length(idx_int_nonAKD) == 1) v_nonAKD[idx_int_nonAKD] <- 1
  if (length(idx_int_AKD)    == 1) v_AKD[idx_int_AKD]       <- 1
  
  # glhtで推定値と95%CI
  ci_nonAKD <- confint(glht(fit, linfct = rbind(nonAKD = v_nonAKD), vcov = vcov(fit)))$confint
  ci_AKD    <- confint(glht(fit, linfct = rbind(AKD    = v_AKD),    vcov = vcov(fit)))$confint
  
  tibble(
    Group    = factor(c("nonAKD","AKD"), levels = c("nonAKD","AKD")),
    Estimate = c(ci_nonAKD[,"Estimate"], ci_AKD[,"Estimate"]),
    CI_Lower = c(ci_nonAKD[,"lwr"],      ci_AKD[,"lwr"]),
    CI_Upper = c(ci_nonAKD[,"upr"],      ci_AKD[,"upr"]),
    Horizon  = horizon_label
  )
}

slope_all_2groups <- dplyr::bind_rows(
  extract_slope_df_2groups(fit_window_2_within_1year, "within1 year"),
  extract_slope_df_2groups(fit_window_2_within_2year, "within2 years"),
  extract_slope_df_2groups(fit_window_2_within_3year, "within3 years"),
  extract_slope_df_2groups(fit_window_2,              "All period")
) %>%
  mutate(
    Horizon = factor(Horizon, levels = c("within1 year","within2 years","within3 years","All period"))
    # ← PlotValue/PlotLow/PlotHigh は作らない
  )

y_min <- floor(min(slope_all_2groups$CI_Lower, na.rm = TRUE)) - 0.5
y_max <- ceiling(max(slope_all_2groups$CI_Upper, na.rm = TRUE)) + 0.5

p2 <- ggplot(slope_all_2groups, aes(x = Group, y = Estimate, fill = Group)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.15, linewidth = 0.6) +
  facet_wrap(~ Horizon, nrow = 1) +
  scale_y_continuous(limits = c(y_min, y_max)) +
  labs(
    title = "Estimated eGFR slope by group across time horizons",
    x = NULL,
    y = "Slope (mL/min/1.73 m2 per year, 95% CI)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

print(p2)



# 論文用に高解像度保存（下向き棒グラフ）
ggsave("Fig_Slope_AKD_vs_nonAKD_downward.tiff", plot = p2,
       width = 10.5, height = 3.2, units = "in", dpi = 600, compression = "lzw")



#アウトカム解析に必要なcodeのみ####
#jin1_Eligibileから一意のidだけ抽出
jin1_Eligibile_unique_id_ver2 <- jin1_Eligibile_ver2 %>%
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
jin1_Eligibile_unique_death_ver2 <- jin1_Eligibile_unique_id_ver2 %>%
  filter(jin_status %in% c("nonAKD", "AKD")) %>%  # nonAKDとAKDだけ抽出
  mutate(
    group = jin_status,  # "nonAKD" または "AKD"をそのままgroupに
    time_years = as.numeric(last_follow_death - index_plus_90) / 365.25
  ) %>%
  filter(time_years >= 0)
#モデル作成、グラフ（死亡）
fit_death_2 <- survfit(Surv(time_years, primary_death) ~ group, data = jin1_Eligibile_unique_death_ver2)

ggsurvplot(fit_death_2, data = jin1_Eligibile_unique_death_ver2,
           fun = "event",                     # 1 - survival（累積死亡率）
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE,
           xlim = c(0, 10),                   # 10年まで表示
           ylim = c(0, 0.4),
           title = "death",
           xlab = "year",
           ylab = "event")
#競合エンドポイントについてのカプランマイヤー2025/8/5####
install.packages("survminer")
library(survival)
library(survminer)
library(dplyr)

#競合エンドポイントについてのカプランマイヤーデータ整形
jin1_Eligibile_unique_composite_ver2 <- jin1_Eligibile_unique_id_ver2  %>%
  mutate(
    group = case_when(
      jin_status == "nonAKD" ~ "nonAKD",
      jin_status == "AKD" ~  "AKD",
      TRUE ~ NA_character_
    ),
    time_years = as.numeric(last_follow_composite - index_plus_90) / 365.25　
  ) %>%
  filter(!is.na(group), time_years >= 0)

fit_composite <- survfit(Surv(time_years, primary_composite_event) ~ group, data = jin1_Eligibile_unique_composite_ver2)

ggsurvplot(fit_composite, data = jin1_Eligibile_unique_composite_ver2,
           fun = "event",                     # 1 - survival（累積死亡率）
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE,
           xlim = c(0, 10),                   # 10年まで表示
           ylim = c(0, 0.4),
           title = "composite_event",
           xlab = "year",
           ylab = "event")

#coxphに必要なcodeのみ####
jin1_Eligibile_cox <- jin1_Eligibile_unique_death_ver2 %>% #arb/aceiを変数にする
  filter(jin_status %in% c("nonAKD", "AKD")) %>% 
  mutate(arb_acei_use = if_else(arb == 1 | acei == 1, 1, 0))
cox_model_death_v3 <- coxph(
  Surv(time_years, primary_death) ~ 
    group + age + index_cre + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + 
    dn12 + dn13 + dn14 + dn15,
  data = jin1_Eligibile_cox)



#感度分析：CKDあるなし#####
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
