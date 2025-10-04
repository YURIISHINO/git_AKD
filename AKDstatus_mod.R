setwd("E:/R")
library(readr)
jin1_index_cre_egfr_3 <- read_csv("jin1_index_cre_egfr_3.csv", locale = locale(encoding = "SHIFT-JIS"))
jin1_index_cre_egfr_3

#AKD判定_90日以内データを抽出

# index_dateから90日以内のdateの行を抽出
jin1_index_90 <- jin1_index_cre_egfr_3 %>%
  filter(date >= index_date - 90 & date <= index_date + 90)

# 結果を表示
print(jin1_index_90)
View(jin1_index_90)

###################################################
# AKD判定の途中まで
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

#####################################################
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

jin1_AKD <- jin1_AKD_with_details %>%
  dplyr::select(id, AKD_status = final_AKD_status)

View(jin1_AKD_with_details)
View(jin1_AKD)