setwd("X:/R")
library(readr)
library(dplyr)
jin1_index_cre_egfr_3 <- read_csv("jin1_index_cre_egfr_3.csv", locale = locale(encoding = "SHIFT-JIS"))
jin1_index_cre_egfr_3
# index_dateから90日以内のdateの行を抽出
jin1_index_90 <- jin1_index_cre_egfr_3 %>%
  filter(date >= index_date - 90 & date <= index_date + 90)
# AKD判定の途中まで################################
jin1_AKD_intermediate <- jin1_index_90 %>%
  group_by(id) %>%
  mutate(
    AKD_status = case_when(
      # index_date以外のdateが存在しない場合
      n() == 1 ~ "nd",
      
      # index_dateよりも前のdateのcreがindex_creの2/3以下の場合
      date < index_date & cre <= index_cre * (2/3) ~ "AKD_cre1",
      
      # index_dateよりも前のdateのegfrがindex_egfrの0.74倍以下の場合
      #date < index_date & egfr <= index_egfr * 0.74 ~ "AKD_egfr1",
      #修正
      #index_dateよりも前のdateのegfrがindex_egfrの1/0.65≒1.54倍以上の場合
      date < index_date & egfr >= index_egfr /0.65 ~ "AKD_egfr1",
      
      # index_dateよりも前のdateのegfrが60以上かつindex_egfrが60未満の場合
      date < index_date & egfr >= 60 & index_egfr < 60 ~ "AKD_egfr2",
      
      # 一旦nonAKDと判定
      TRUE ~ "nonAKD"
    )
  ) %>%
  # 後続の条件を再判定
  mutate(
    AKD_status = case_when(
      # index_dateよりも後のdateのcreがindex_creの1.5倍以上の場合
      date > index_date & cre >= index_cre * 1.5 ~ "AKD_cre2",
      
      # index_dateよりも後のdateのegfrがindex_egfrの1.35倍以上の場合
      # date > index_date & egfr >= index_egfr * 1.35 ~ "AKD_egfr3",
      #修正
      # index_dateよりも後のdateのegfrがindex_egfrの0.65倍以下の場合
      date > index_date & egfr <= index_egfr * 0.65 ~ "AKD_egfr3",
      
      # それ以外はnonAKD
      TRUE ~ AKD_status
    )
  )
View(jin1_AKD_intermediate)
library(dplyr)

jin1_AKD_intermediate %>%
  group_by(AKD_status) %>%
  summarise(n_id = n_distinct(id)) %>%
  arrange(desc(n_id))

library(dplyr)

akd_combo <- jin1_AKD_intermediate %>%
  group_by(id) %>%
  summarise(
    combo = paste(sort(unique(AKD_status)), collapse = " & "),
    .groups = "drop"
  ) %>%
  count(combo, name = "n_id") %>%
  arrange(desc(n_id))

akd_combo

#最終判定###################################################
#AKD_cre1 or AKD_cre2 → "AKD_cre"
#AKD_egfr1 or AKD_egfr3 → "AKD_egfr"
#AKD_egfr2 → "AKD_egfr_under60"

library(dplyr)

jin1_AKD_with_details <- jin1_AKD_intermediate %>%
  group_by(id) %>%
  mutate(
    final_AKD_status = case_when(
      any(AKD_status %in% c("AKD_cre1", "AKD_cre2"), na.rm = TRUE) ~ "AKD_cre",             # 最優先
      any(AKD_status %in% c("AKD_egfr1", "AKD_egfr3"), na.rm = TRUE) ~ "AKD_egfr",           # 次点
      any(AKD_status == "AKD_egfr2", na.rm = TRUE)                      ~ "AKD_egfr_under60", # 最後
      all(is.na(AKD_status) | AKD_status == "nd")                       ~ "nd",
      TRUE                                                              ~ "nonAKD"
    )
  ) %>%
  ungroup()

View(jin1_AKD_with_details)

jin1_AKD_with_details %>%
  group_by(final_AKD_status) %>%
  summarise(n_id = n_distinct(id)) %>%
  arrange(desc(n_id))

jin1_AKD_mod <- jin1_AKD_with_details %>%
  dplyr::select(id, AKD_status = final_AKD_status)

View(jin1_AKD_with_details)
View(jin1_AKD_mod)

setwd("X:/R")
write.csv(jin1_AKD_mod, file = "jin1_AKD_mod.csv", row.names = FALSE)
library(readr)
setwd("X:/R")
jin1_AKD_mod <- read_csv("jin1_AKD_mod.csv", locale = locale(encoding = "SHIFT-JIS"))

#集計
library(dplyr)
AKD_status_counts <- jin1_AKD_mod %>%
  group_by(AKD_status) %>%
  summarise(n_unique_id = n_distinct(id)) %>%
  arrange(AKD_status)
AKD_status_counts
# "AKD_egfr_under60"だけ抜き出し
library(dplyr)
id_under60 <- jin1_AKD_mod %>%
  filter(AKD_status == "AKD_egfr_under60") %>%
  distinct(id)
id_under60
write.csv(id_under60, "E:/R/id_under60.csv", row.names = FALSE)

#statusまとめ####
library(readr)
setwd("X:/R")
jin1_AKI <- read_csv("jin1_AKI.csv", locale = locale(encoding = "SHIFT-JIS"))
jin1_AKD_mod <- read_csv("jin1_AKD_mod.csv", locale = locale(encoding = "SHIFT-JIS"))
jin1_CKD_all <- read_csv("jin1_CKD_all.csv", locale = locale(encoding = "SHIFT-JIS"))
# AKD_status, AKI_status, CKD_statusを連結するコード
jin1_status <- jin1_index_cre_egfr_3 %>%
  # AKD_statusの連結
  left_join(jin1_AKD_mod %>%
              group_by(id) %>%
              summarise(AKD_status = paste(unique(AKD_status), collapse = ", ")), 
            by = "id") %>%
  mutate(AKD_status = if_else(is.na(AKD_status), "nd", AKD_status)) %>%
  # AKI_statusの連結
  left_join(jin1_AKI %>%
              group_by(id) %>%
              summarise(AKI_status = paste(unique(AKI_status), collapse = ", ")), 
            by = "id") %>%
  mutate(AKI_status = if_else(is.na(AKI_status), "nd", AKI_status)) %>%
  # CKD_statusの連結
  left_join(jin1_CKD_all %>%
              group_by(id) %>%
              summarise(CKD_status = paste(unique(CKD_status), collapse = ", ")), 
            by = "id") %>%
  mutate(CKD_status = if_else(is.na(CKD_status), "nd", CKD_status))

# 結果を表示
print(jin1_status)
View(jin1_status)
colnames(jin1_status)

jin1_status %>%
  group_by(AKD_status) %>%
  summarise(n_id = n_distinct(id)) %>%
  arrange(desc(n_id))

setwd("X:/R")
write.csv(jin1_status, file = "jin1_status.csv", row.names = FALSE)

library(readr)
setwd("X:/R")
jin1_status <- read_csv("jin1_status.csv", locale = locale(encoding = "SHIFT-JIS"))
jin1_status %>%
  group_by(AKD_status) %>%
  summarise(n_id = n_distinct(id)) %>%
  arrange(desc(n_id))

#recovery#####
library(dplyr)
setwd("X:/R")
jin1_status <- read_csv("jin1_status.csv", locale = locale(encoding = "SHIFT-JIS"))
#AKDの行だけ抜き出し
jin1_status_AKD <- jin1_status %>%
  filter(AKD_status %in% c("AKD_cre", "AKD_egfr") | AKI_status == "AKI") %>%
 # filter(AKD_status %in% c("AKD_cre", "AKD_egfr", "AKD_egfr_under60") | AKI_status == "AKI") %>%
  dplyr::select(-AKI_status, -AKD_status, -CKD_status)

# baseline_cre と baseline_egfr を計算して追加
jin1_status_AKD_baseline <- jin1_status_AKD %>%
  group_by(id) %>%
  mutate(
    baseline_cre = min(cre[date >= (index_date - 365) & date <= index_date], na.rm = TRUE),
    baseline_egfr = max(egfr[date >= (index_date - 365) & date <= index_date], na.rm = TRUE)
  ) %>%
  ungroup()

#index_dateからの日数を計算し列に追加
jin1_status_AKD_baseline <- jin1_status_AKD_baseline %>%
  mutate(days_from_index = as.numeric(difftime(date, index_date, units = "days")))

# 結果を確認
View(jin1_status_AKD_baseline)

# baseline_cre/egfrにNAがないか確認→なし
jin1_missing_baseline <- jin1_status_AKD_baseline %>%
  filter(is.na(baseline_cre) | is.na(baseline_egfr))
jin1_missing_baseline 

{#150~210日データで
library(dplyr)

jin1_status_AKD_recovery <- jin1_status_AKD_baseline %>%
  group_by(id) %>%
  mutate(
    has_150_210 = any(between(days_from_index, 150, 210)), #150~210にデータがあるかどうか
    recovery_flag = case_when(
      days_from_index >= 150 & days_from_index <= 210 &
        (cre <= baseline_cre * 1.25 | egfr >= baseline_egfr * 0.8) ~ 1,　#150～210日以内」かつ「回復基準を満たす（creが基準の1.25倍以下 or egfrが0.8倍以上）」なら 1: recovery
      days_from_index >= 150 & days_from_index <= 210 ~ 2, #150～210日以内だが、上記基準を満たしていない」なら 2:non-recovery
      TRUE ~ NA_real_ #それ以外の行（期間外の行）は NA にしておく
    ),
    `150_210recovery` = case_when(
      has_150_210 == FALSE ~ 0,                          # 範囲内データなし → 0:no data
      any(recovery_flag == 1, na.rm = TRUE) ~ 1,         # 条件を満たす → 1: recovery
      has_150_210 == TRUE ~ 2,                           # 範囲内あり、条件満たさず → 2:non-recovery
      TRUE ~ NA_real_
    )
  ) %>%
  ungroup() %>%
  dplyr::select(-has_150_210, -recovery_flag)
View(jin1_status_AKD_recovery)

# 要素ごとの一意のid数を数える
jin1_status_AKD_recovery %>%
  group_by(`150_210recovery`) %>%
  summarise(unique_id_count = n_distinct(id))}#150~210日データで
{#90~150日データで

jin1_status_AKD_recovery <- jin1_status_AKD_recovery %>%
  group_by(id) %>%
  mutate(
    has_90_150 = any(between(days_from_index, 90, 150)),
    recovery_flag_90 = case_when(
      days_from_index >= 90 & days_from_index <= 150 &
        (cre <= baseline_cre * 1.25 | egfr >= baseline_egfr * 0.8) ~ 1,
      days_from_index >= 90 & days_from_index <= 150 ~ 2,
      TRUE ~ NA_real_
    ),
    `90_150recovery` = case_when(
      has_90_150 == FALSE ~ 0,                          # 範囲内データなし → 0
      any(recovery_flag_90 == 1, na.rm = TRUE) ~ 1,     # 条件を満たす → 1
      has_90_150 == TRUE ~ 2,                           # 範囲内あり、条件満たさず → 2
      TRUE ~ NA_real_
    )
  ) %>%
  ungroup() %>%
  dplyr::select(-has_90_150, -recovery_flag_90)
View(jin1_status_AKD_recovery)}##90~150日データで
setwd("X:/R")
write.csv(jin1_status_AKD_recovery, file = "jin1_status_AKD_recovery.csv", row.names = FALSE)
jin1_status_AKD_recovery %>%
  group_by(`150_210recovery`, `90_150recovery`) %>%
  summarise(unique_ids = n_distinct(id), .groups = "drop")

#jin1_Eligibileにまとめ####
library(dplyr)
library(lubridate)
library(readr)
setwd("X:/R")
jin1_status_AKD_recovery <- read_csv("jin1_status_AKD_recovery.csv", locale = locale(encoding = "SHIFT-JIS"))
jin1_med <- read_csv("jin1_med.csv", locale = locale(encoding = "SHIFT-JIS"))
dn_complete　<- read_csv("dn_complete.csv", locale = locale(encoding = "SHIFT-JIS"))
jin1_status <- read_csv("jin1_status.csv", locale = locale(encoding = "SHIFT-JIS"))
jin1_index_cre_egfr_3 <- read_csv("jin1_index_cre_egfr_3.csv", locale = locale(encoding = "SHIFT-JIS"))
#一意のidとrecoveryの結果を抜き出す
id_recovery_unique <- jin1_status_AKD_recovery %>%
  dplyr::select(id, `150_210recovery`, `90_150recovery`) %>%
  distinct()
#一意のidとAKDstatus/AKIstatus/CKDstatusの結果を抜き出す
id_status_unique <- jin1_status %>%
  dplyr::select(id, AKD_status,AKI_status,CKD_status) %>%
  distinct()
print(id_status_unique)

# 基準: jin1_index_follow の index_date を採用
jin1_Eligibile <- jin1_index_cre_egfr_3 %>%
  left_join(id_recovery_unique %>% dplyr::select(-any_of("index_date")), by = "id") %>%
  left_join(jin1_med           %>% dplyr::select(-any_of("index_date")), by = "id") %>%
  left_join(dn_complete        %>% dplyr::select(-any_of("index_date")), by = "id") %>%
  left_join(id_status_unique   %>% dplyr::select(-any_of("index_date")), by = "id") %>%
  mutate(age = floor(lubridate::interval(birth, index_date) / lubridate::years(1)))
#　index_date時の年齢を追加
jin1_Eligibile <- jin1_Eligibile %>%
  mutate(age = floor(interval(birth, index_date) / years(1)))
colnames(jin1_Eligibile)
#AKD/nonAKD/NKDを定義####
#AKD:AKDstatusが"AKD_cre""AKD_egfr","AKD_egfr_under60"もしくはAKI_statusが"AKI"
#上記と反対だがinclusion
#上記以外がother（すべてのステータスがNA）
# ステップ1：jin_status_tmp 作成
jin1_Eligibile <- jin1_Eligibile %>%
  mutate(
    jin_status_tmp = case_when(
      AKD_status %in% c("AKD_cre", "AKD_egfr") | AKI_status == "AKI" ~ "AKD",
      #AKD_status %in% c("AKD_cre", "AKD_egfr","AKD_egfr_under60") | AKI_status == "AKI" ~ "AKD",
      AKD_status %in% c("nd", "nonAKD","AKD_egfr_under60") &
        AKI_status %in% c("nd", "nonAKI") ~ "nonAKD",
      TRUE ~ "other"
    ),
    jin_status_tmp = factor(jin_status_tmp, levels = c("AKD", "nonAKD", "other"))
  )
# ステップ2：idごとの優先順位つき代表ステータス取得
id_status_priority <- jin1_Eligibile %>%
  group_by(id) %>%
  summarise(jin_status_rep = first(jin_status_tmp[order(jin_status_tmp)]), .groups = "drop")

# ステップ3：元データに結合して因子化
jin1_Eligibile <- jin1_Eligibile %>%
  left_join(id_status_priority, by = "id") %>%
  dplyr::select(-jin_status_tmp) %>%
  mutate(jin_status = factor(jin_status_rep, levels = c("AKD", "nonAKD", "other"))) %>%
  dplyr::select(-jin_status_rep)
#集計
jin1_Eligibile %>%
  group_by(jin_status) %>%
  summarise(n_unique_ids = n_distinct(id), .groups = "drop")
jin1_Eligibile %>%
  group_by(exclude,jin_status) %>%
  summarise(n_id = n_distinct(id)) %>%
  arrange(desc(n_id)) # otherは存在しない



#アウトカム関連まとめ####
#death,ESKD(hd_new,rr)日にちをつける
death
hd_new
rr
jin1_Eligibile <- jin1_Eligibile %>%
  left_join(death, by = "id") %>%
  left_join(hd_new, by = "id") %>%
  left_join(rr, by = "id")
#hd_new,rrはESKD_dateにまとめ
jin1_Eligibile <- jin1_Eligibile %>%
  mutate(
    ESKD_date = pmin(hd_new_start, rr_date, na.rm = TRUE)
  )
#最後のデータが存在する日をlast_dataとする
jin1_Eligibile <- jin1_Eligibile %>%
  group_by(id) %>%
  mutate(
    last_data = max(date, na.rm = TRUE)
  ) %>%
  ungroup()

#フォロー開始日、アウトカムの0.1データ
library(dplyr)
library(lubridate)
jin1_Eligibile <- jin1_Eligibile %>%
  mutate(
    # 210日後の基準日を作成
    index_plus_210 = index_date + 210,
    
    # 死亡とESKDの発症が210日以降かどうかを判定
    primary_death = if_else(!is.na(death_date) & death_date > index_plus_210, 1, 0),
    primary_ESKD = if_else(!is.na(ESKD_date) & ESKD_date > index_plus_210, 1, 0)
  )

jin1_Eligibile <- jin1_Eligibile %>%
  mutate(
    last_follow_death = if_else(primary_death == 1, as.Date(death_date), as.Date(last_data))
  )

jin1_Eligibile <- jin1_Eligibile %>%
  mutate(
    # composite endpoint が起きたか判定（primary_death または primary_ESKD のどちらかが1）
    primary_composite_event = (primary_death == 1 | primary_ESKD == 1),
    # event_date は death_date と ESKD_date の早い方（Date型で）
    event_date = pmin(as.Date(death_date), as.Date(ESKD_date), na.rm = TRUE),
    # composite_event が TRUE（=1）なら event_date、それ以外は last_data
    last_follow_composite = if_else(primary_composite_event == 1, event_date, as.Date(last_data))
  )


View(jin1_Eligibile)
jin1_Eligibile %>%
  group_by(jin_status) %>%
  summarise(unique_ids = n_distinct(id)) %>%
  arrange(desc(unique_ids))  # オプション：多い順に並べる

colnames(jin1_Eligibile)

#AKD,AKI発生日を入れる####
jin1_Eligibile <- jin1_Eligibile %>%
  left_join(jin1_AKI_date_nonNA_unique, by = "id") %>%
  left_join(jin1_AKD_date_nonNA_unique, by = "id")

#腎摘出######
library(readr)
nephrectomy <- read_csv("E:/nephrectomy.csv", locale = locale(encoding = "SHIFT-JIS"), skip = 1)
print(nephrectomy)
View(nephrectomy)
colnames(nephrectomy)

library(dplyr)
library(stringr)
library(lubridate)
library(stringi)  # 全角→半角

nephrectomy_sub <- nephrectomy %>%
  dplyr::select(ID, 手術日, 術式...4) %>%
  filter(術式...4 %in% c("腎摘", "腎摘＋部切", "腎摘出術", "腎摘出術 K772-00")) %>%
  mutate(
    # 日付文字列を正規化（全角→半角、余分な空白除去、非数字→"-"）
    手術日_norm = stri_trans_general(手術日, "Fullwidth-Halfwidth"),
    手術日_norm = str_squish(手術日_norm),
    手術日_norm = str_replace_all(手術日_norm, "[^0-9]", "-"),
    手術日_norm = str_replace_all(手術日_norm, "-+", "-"),
    手術日_norm = str_replace_all(手術日_norm, "^-|-$", ""),
    
    # 日付に変換（"YYYY-M-D" でもOK）
    ope_date = ymd(手術日_norm, quiet = TRUE),
    
    # ID→数値id
    id = suppressWarnings(as.numeric(ID))
  ) %>%
 dplyr::select(id, ope_date, 術式...4)

# 変換できなかった行があるか軽くチェック（必要なら）
sum(is.na(nephrectomy_sub$ope_date))

nephrectomy_sub %>%
  group_by(術式...4) %>%
  summarise(unique_ids = n_distinct(id)) %>%
  arrange(desc(unique_ids))  # オプション：多い順に並べる

library(dplyr)
library(lubridate)

# ① idごとの代表 index_date（最も早い日）を用意
id_index <- jin1_Eligibile %>%
  filter(!is.na(index_date)) %>%
  arrange(id, index_date) %>%
  group_by(id) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(id, index_date)

# ② nephrectomy_sub に index_date を突合し、-90～210日以内の手術か判定
nephrectomy_flag <- nephrectomy_sub %>%
  left_join(id_index, by = "id") %>%
  mutate(
    nephrectomy = if_else(
      !is.na(index_date) & !is.na(ope_date) &
        #ope_date >= index_date & ope_date <= index_date + days(210), #14
        ope_date >= index_date - days(90) & ope_date <= index_date + days(210), #19
      1L, 0L
    )
  )

# 確認
print(nephrectomy_flag, n = 10)

nephrectomy_1 <- nephrectomy_flag %>%
  filter(nephrectomy == 1)
print(nephrectomy_1, n = Inf)  # すべて表示したい場合

#出産#####
library(readr)
pregnancy <- read_csv("E:/pregnancy.csv", locale = locale(encoding = "SHIFT-JIS"))
print(pregnancy)
library(dplyr)
library(stringr)
library(stringi)  # 全角→半角

pregnancy_sub <- pregnancy %>%
  transmute(
    # ID列：全角→半角、前後の引用符/空白/改行を除去、数字以外を除去（先頭ゼロは保持）
    id = `あ母入院番号１` |>
      stri_trans_general("Fullwidth-Halfwidth") |>
      str_squish() |>
      str_remove_all('^"|"$') |>
      str_replace_all("[^0-9]", "") |>
      as.numeric(),
    # 分娩日：日付型へ（もともと<date>でも明示的に統一）
    date_of_birth = as.Date(`い分娩日`)
  )
sum(is.na(pregnancy_sub$date_of_birth))
print(pregnancy_sub, n = 10)

# ② pregnancy_subに index_date を突合し、-90～210日以内の出産か判定
pregnancy_flag <- pregnancy_sub %>%
  left_join(id_index, by = "id") %>%
  mutate(
    pregnancy = if_else(
      !is.na(index_date) & !is.na(date_of_birth) &
        #date_of_birth >= index_date & date_of_birth <= index_date + days(210),           #380 
        date_of_birth >= index_date - days(90) & date_of_birth <= index_date + days(210), #404
      1L, 0L
    )
  )
pregnancy_1 <- pregnancy_flag %>%
  filter(pregnancy == 1)
print(pregnancy_1, n = Inf)
# 重複しているidの行すべて表示→すべて同じidと出産日のため重複は削除
pregnancy_dup <- pregnancy_1 %>%
  filter(duplicated(id) | duplicated(id, fromLast = TRUE))
print(pregnancy_dup, n = Inf)  
pregnancy_unique <- pregnancy_1 %>%
  filter(!(duplicated(id) | duplicated(id, fromLast = TRUE)))
print(pregnancy_unique, n = 10)

library(dplyr)

# 1) idごとのフラグを用意（重複idがあっても安全）
neph_flag <- nephrectomy_1 %>%
  distinct(id) %>%
  mutate(nephrectomy_flag = 1L)

preg_flag <- pregnancy_unique %>%
  distinct(id) %>%
  mutate(pregnancy_flag = 1L)

# 2) 突合してフラグ付け → excludeを書き換え（includeのみ）
jin1_Eligibile_removed <- jin1_Eligibile %>%
  left_join(neph_flag, by = "id") %>%
  left_join(preg_flag, by = "id") %>%
  mutate(
    nephrectomy_flag = coalesce(nephrectomy_flag, 0L),
    pregnancy_flag   = coalesce(pregnancy_flag,   0L),
    exclude_before   = exclude,  # 監査用に保存
    exclude = case_when(
      # 両方に該当 → nephrectomyを優先
      exclude_before == "include" & nephrectomy_flag == 1L              ~ "nephrectomy",
      exclude_before == "include" & nephrectomy_flag == 0L & pregnancy_flag == 1L ~ "pregnancy",
      TRUE ~ exclude_before
    )
  )

# 変更件数のサマリ
library(tidyr)

summary_changes <- jin1_Eligibile_removed %>%
  transmute(
    changed = exclude_before == "include" & exclude %in% c("nephrectomy","pregnancy"),
    reason  = exclude
  ) %>%
  count(changed, reason, name = "n")

print(summary_changes)

# 重複（両フラグ1）の件数確認
both_flag_n <- jin1_Eligibile_removed %>%
  filter(nephrectomy_flag == 1L, pregnancy_flag == 1L) %>%
  nrow()
both_flag_n
library(dplyr)
jin1_Eligibile_removed %>%
  group_by(jin_status, exclude) %>%
  summarise(n = n(), .groups = "drop")
jin1_Eligibile_removed %>%
  filter(exclude == "include") %>%
  group_by(jin_status) %>%
  summarise(n_unique_ids = n_distinct(id), .groups = "drop")


####jin1_EligibileをCSVファイルに書き出し####
setwd("X:/R")
write.csv(jin1_Eligibile_removed, file = "jin1_Eligibile.csv", row.names = FALSE)
library(readr)
setwd("X:/R")
jin1_Eligibile <- read_csv("jin1_Eligibile.csv", locale = locale(encoding = "SHIFT-JIS"))

jin1_Eligibile %>%
  filter(exclude == "include") %>%
  group_by(AKD_status) %>%
  summarise(n_id = n_distinct(id)) %>%
  arrange(desc(n_id))

jin1_Eligibile %>%
  filter(exclude == "include") %>%
  group_by(jin_status) %>%
  summarise(n_unique_ids = n_distinct(id), .groups = "drop")
colnames(jin1_Eligibile)
id_under60 <- read.csv("E:/R/id_under60.csv", fileEncoding = "CP932")
library(dplyr)
jin1_Eligibile_under60 <- jin1_Eligibile %>%
  inner_join(id_under60, by = "id")
library(dplyr)
AKD_status_count_id <- jin1_Eligibile_under60 %>%
  group_by(jin_status) %>%
  summarise(n_unique_id = n_distinct(id)) %>%
  arrange(jin_status)
AKD_status_count_id

{#eGFR・線形混合効果モチE��を走らせる、また解析に忁E��なcodeのみ####
  library(nlme)
  library(dplyr)
  library(ggplot2)
  # ① jin_label の定義と忁E���Eの抽出
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
  
  # ② 90�E�E10日の最大egfr日�E�Eax_egfr_date_210�E�抽出
  egfr_max_date_210 <- jin1_inclusion %>%
    mutate(days_from_index = as.numeric(date - index_date)) %>%
    filter(days_from_index >= 90, days_from_index <= 210) %>%
    group_by(id) %>%
    filter(egfr == max(egfr, na.rm = TRUE)) %>%
    slice(1) %>%
    ungroup() %>%
    dplyr::select(id, max_egfr_date_210 = date)
  
  # ③ 0�E�E0日の最新日�E�Eearest_date_90�E�抽出
  max_date_90 <- jin1_inclusion %>%
    mutate(days_from_index = as.numeric(date - index_date)) %>%
    filter(days_from_index >= 0, days_from_index <= 90) %>%
    group_by(id) %>%
    filter(days_from_index == max(days_from_index, na.rm = TRUE)) %>%
    slice(1) %>%
    ungroup() %>%
    dplyr::select(id, nearest_date_90 = date)
  
  # ④ 允E��ータに結合
  jin1_inclusion <- jin1_inclusion %>%
    left_join(egfr_max_date_210, by = "id") %>%
    left_join(max_date_90, by = "id")
  
  # ⑤ time0 の定義�E�優先度�E�max_egfr_date_210 > nearest_date_90 > index_date�E�E
  jin1_inclusion <- jin1_inclusion %>%
    mutate(time0 = case_when(
      !is.na(max_egfr_date_210) ~ max_egfr_date_210,
      !is.na(nearest_date_90) ~ nearest_date_90,
      TRUE ~ index_date
    ))
  
  # ⑥ time0 からの年差を計箁E
  jin1_inclusion <- jin1_inclusion %>%
    mutate(
      days_from_time0 = as.numeric(difftime(date, time0, units = "days")),
      years_from_time0 = round(days_from_time0 / 365.25, 3)
    )
  
  # ⑦ time0 に一致する egfr�E�褁E��あれば最大値�E�を抽出
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
  
  # CSV保孁E
  library(readr)
  write_csv(
    jin1_inclusion,
    "jin1_inclusion.csv"   # 保存ファイル名（作業チE��レクトリに保存されます！E
  ) 
}
