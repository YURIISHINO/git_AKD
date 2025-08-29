library(dplyr)
library(tibble)
library(readr)
library(lubridate)
#jin_label作成####
library(readr)
library(dplyr)

# 読み込み設定
loc <- locale(encoding = "SHIFT-JIS")

# 必要な列だけ読み込み（.default = "c" は他列を文字列扱いにして後で選択）
cre_2012 <- read_csv("jin/cre_over18/cre_2012_over18.csv",
                     locale = loc, skip = 3,
                     col_types = cols(.default = "c")) %>%
  select(患者ID, 科ｺｰﾄﾞ, 検査日)
cre_2013 <- read_csv("jin/cre_over18/cre_2013_over18.csv",
                     locale = loc, skip = 3,
                     col_types = cols(.default = "c")) %>%
  select(患者ID, 科ｺｰﾄﾞ, 検査日)
# 行結合
cre_2012_2013 <- bind_rows(cre_2012, cre_2013)
# 検査日を日付型に変換（必要に応じて）
cre_2012_2013 <- cre_2012_2013 %>%
  mutate(検査日 = as.Date(as.character(検査日), format = "%Y%m%d"))

cre_2012_2013_sub <- cre_2012_2013 %>%
  transmute(
    id   = as.numeric(患者ID),
    code = 科ｺｰﾄﾞ,
    date = 検査日
  ) %>%
  group_by(id, date) %>%
  summarise(
    code = paste(unique(code), collapse = "&"),  # A&B形式にまとめる
    .groups = "drop"
  )
cre_2012_2013_sub
#判定期間にデータがないものはnon-recoveryに#####
jin1_Eligibile_include_code <- jin1_Eligibile %>%
  left_join(cre_2012_2013_sub, by = c("id", "date")) %>%
  filter(exclude == "include") %>%
  mutate(
    jin_label = case_when(
      jin_status == "nonAKD" ~ "nonAKD",
      jin_status == "AKD" & `150_210recovery` == 1 ~ "Recovery",
      jin_status == "AKD" & `150_210recovery` == 2 ~ "Non-Recovery",
      jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 1 ~ "Recovery",
      jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` %in% c(0, 2) ~ "Non-Recovery",
      TRUE ~ NA_character_
    )
  )
jin1_Eligibile_include_code <- jin1_Eligibile_include_code %>%
  mutate(jin_label = factor(jin_label, levels = c("nonAKD", "Non-Recovery", "Recovery")))

#判定期間にデータがないものを除外#####
jin1_Eligibile_include_code <- jin1_Eligibile %>%
  left_join(cre_2012_2013_sub, by = c("id", "date")) %>%
  filter(exclude == "include") %>%
  mutate(
    jin_label = case_when(
      jin_status == "nonAKD" ~ "nonAKD",
      jin_status == "AKD" & `150_210recovery` == 1 ~ "Recovery",
      jin_status == "AKD" & `150_210recovery` == 2 ~ "Non-Recovery",
      jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 1 ~ "Recovery",
      jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 0 ~ "No-data",
      jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 2 ~ "Non-Recovery",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(jin_label = factor(jin_label, levels = c("nonAKD", "No-data", "Non-Recovery", "Recovery")))


#####出現回数をみる####
main_code_df_n <- jin1_Eligibile_include_code %>%
  filter(!is.na(code)) %>%
  filter(date <= index_date) %>%
  group_by(id, code) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(id) %>%
  filter(n == max(n)) %>%
  summarise(
    main_code = paste(sort(unique(code)), collapse = " & "),
    count_main_code = max(n),  # 再頻出codeの出現回数を保持
    .groups = "drop"
  )
jin1_Eligibile_include_code_n <- jin1_Eligibile_include_code %>%
  left_join(main_code_df_n, by = "id")
colnames(jin1_Eligibile_include_code_n)

main_code_stats <- jin1_Eligibile_include_code_n %>%
  distinct(id, jin_label, main_code, count_main_code) %>%
  group_by(jin_label, main_code) %>%
  summarise(
    n_patients = n(),  # そのmain_codeを持つ患者数
    mean_n = mean(count_main_code),
    median_n = median(count_main_code),
    max_n = max(count_main_code),
    .groups = "drop"
  ) %>%
  arrange(desc(n_patients))
print(main_code_stats, n = Inf)

#消化器外科まとめ####
library(dplyr)
library(stringr)
library(tidyr)

# Digestive Surgery に該当するコード
digestive_codes <- c("DK", "DM", "DL", "DR")

# main_code を再構成する関数
normalize_main_code <- function(code_string) {
  codes <- str_split(code_string, " & ", simplify = TRUE)
  codes <- ifelse(codes %in% digestive_codes, "DigestiveSurgery", codes)
  codes <- sort(unique(codes))
  paste(codes, collapse = " & ")
}

# 再頻出科のまとめ（index_date以前に限定）
main_code_df <- jin1_Eligibile_include_code %>%
  filter(!is.na(code)) %>%
  filter(date <= index_date) %>%  # index_date以前のみに限定
  group_by(id, code) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(id) %>%
  filter(n == max(n)) %>%
  summarise(
    main_code = paste(sort(unique(code)), collapse = " & "),
    .groups = "drop"
  ) %>%
  # 消化器外科コードをまとめる
  mutate(main_code = sapply(main_code, normalize_main_code))

# main_code を元データに結合
jin1_Eligibile_include_code <- jin1_Eligibile_include_code %>%
  left_join(main_code_df, by = "id")

# jin_label 毎に再頻出科を集計
main_code_summary <- jin1_Eligibile_include_code %>%
  distinct(id, jin_label, main_code) %>%
  count(jin_label, main_code, name = "n") %>%
  arrange(desc(n))

# 上位10位を選択
main_code_by_label_top10 <- main_code_summary %>%
  group_by(jin_label) %>%
  slice_max(order_by = n, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  rename(n_patients = n)
print(main_code_by_label_top10, n = Inf)

#科コードと科名の対応表
department_lookup_jp <- tibble(
  main_code = c(
    "AL", "AQ", "HH", "AP", "AM", "AG", "AN", "AH", "FF",
    "DE", "LL", "AR", "MM", "DQ", "GG", "EE", "DN",
    "DH", "DF", "BB", "PS", "RR", "AK", "CC", "KK", "YY",
    "PP", "NN", "AS", "FS", "DG", "DigestiveSurgery"
  ),
  clinical_department = c(
    "内分泌", "血液内", "泌尿", "循環内", "呼吸内", "消化内", "肝臓内", "腎臓内", "整形",
    "心臓外", "耳鼻", "免疫内", "産婦", "救急部", "皮膚", "脳外", "血管外",
    "乳腺外", "呼吸外", "精神", "形成", "歯口", "脳神内", "小児", "眼科", "担当部",
    "麻酔", "放射治", "臨薬", "リハ科", "小児外", "消化器外科"
  )
)
#英語####
department_lookup_en <- tibble(
  main_code = c(
    "AL", "AQ", "HH", "AP", "AM", "AG", "AN", "AH", "FF",
    "DE", "LL", "AR", "MM", "DQ", "GG", "EE", "DN",
    "DH", "DF", "BB", "PS", "RR", "AK", "CC", "KK", "YY",
    "PP", "NN", "AS", "FS", "DG", "DigestiveSurgery"
  ),
  clinical_department_en = c(
    "Endocrinology",           # 内分泌
    "Hematology",              # 血液内
    "Urology",                 # 泌尿
    "Cardiology",              # 循環内
    "Respiratory Medicine",    # 呼吸内
    "Gastroenterology",        # 消化内
    "Hepatology",              # 肝臓内
    "Nephrology",              # 腎臓内
    "Orthopedics",             # 整形
    "Cardiovascular Surgery",  # 心臓外
    "Otorhinolaryngology",     # 耳鼻
    "Immunology",              # 免疫内
    "Obstetrics and Gynecology", # 産婦
    "Emergency Medicine",      # 救急部
    "Dermatology",             # 皮膚
    "Neurosurgery",            # 脳外
    "Vascular Surgery",        # 血管外
    "Breast Surgery",          # 乳腺外
    "Thoracic Surgery",        # 呼吸外
    "Psychiatry",              # 精神
    "Plastic Surgery",         # 形成
    "Dentistry and Oral Surgery", # 歯口
    "Neurology",               # 脳神内
    "Pediatrics",              # 小児
    "Ophthalmology",           # 眼科
    "Administrative Dept.",    # 担当部
    "Anesthesiology",          # 麻酔
    "Radiation Therapy",       # 放射治
    "Clinical Pharmacy",       # 臨薬
    "Rehabilitation",          # リハ科
    "Pediatric Surgery",       # 小児外
    "Digestive Surgery"        # 消化器外科
  )
)


#　上位10位表に対応表をつける
main_code_by_label_top10_named_full <- main_code_by_label_top10 %>%
  rowwise() %>%
  mutate(
    clinical_department = {
      code_parts <- strsplit(main_code, " & ")[[1]]
      dept_names <- department_lookup_jp %>%
        filter(main_code %in% code_parts) %>%
        arrange(match(main_code, code_parts)) %>%
        pull(clinical_department)
      paste(dept_names, collapse = " & ")
    }
  ) %>%
  ungroup()
print(main_code_by_label_top10_named_full, n = Inf)

# jin_labelごとの総患者数（指定値）
label_totals <- tibble(
  jin_label = c("nonAKD", "Non-Recovery", "Recovery"),
  total_patients = c(14406, 227, 114)
)

# 構成比を計算
main_code_ratio <- main_code_by_label_top10_named_full %>%
  left_join(label_totals, by = "jin_label") %>%
  mutate(
    pct = round(100 * n_patients / total_patients, 1)
  )
print(main_code_ratio, n = Inf)

library(dplyr)

# 英語名を main_code_ratio に追加
main_code_ratio_en <- main_code_ratio %>%
  left_join(department_lookup_en, by = "main_code")
View(main_code_ratio_en)

# Excelファイルに書き出す
library(writexl)
write_xlsx(main_code_ratio_en, path = "main_code_ratio.xlsx")
#積み上げ縦棒####
library(dplyr)
library(ggplot2)

# 1) 上位5位（同率はすべて採用）＋その他
main_code_ratio_top5 <- main_code_ratio_en %>%
  filter(!is.na(jin_label), !is.na(clinical_department_en), !is.na(pct)) %>%
  group_by(jin_label) %>%
  mutate(rank = dense_rank(desc(pct))) %>%
  mutate(clinical_department_en = if_else(rank <= 5, clinical_department_en, "Other")) %>%
  ungroup() %>%  # ★ いったん解除してから .by を使う
  summarise(pct = sum(pct), .by = c(jin_label, clinical_department_en)) %>%
  group_by(jin_label) %>%
  mutate(pct = 100 * pct / sum(pct)) %>%
  ungroup()


# 2) 色パレット（指定どおり）
color_palette_filtered <- c(
  "Digestive Surgery"         = "#1B9E77",  # 緑系
  "Endocrinology"             = "#E41A1C",  # 赤系
  "Hematology"                = "#377EB8",  # 青系
  "Obstetrics and Gynecology" = "#E78AC3",  # ピンク系
  "Other"                     = "#BEBEBE",  # グレー（その他）
  "Urology"                   = "#FFA07A",  # 薄い赤系（ライトサーモン）
  "Cardiology"                = "#984EA3",  # 紫
  "Cardiovascular Surgery"    = "#66C2A5",  # ティール
  "Emergency Medicine"        = "#FF7F00",  # オレンジ
  "Otorhinolaryngology"       = "#FFD700"   # 黄色
)
library(dplyr)
library(ggplot2)
library(forcats)

# main_code_ratio_top5 : すでに作成済み（jin_label, clinical_department_en, pct がある想定）

# --- 1) まずはプロット用の基礎データ（df_plot2）を作る ---
df_plot2 <- main_code_ratio_top5 %>%
  transmute(
    jin_label,
    dept = clinical_department_en,
    pct  = pct
  ) %>%
  # dept をファクター化（この後でレベル順を上書き）
  mutate(dept_fac = factor(dept))

# --- 2) 凡例レベルの順序（Other を最後に） ---
legend_lvls <- levels(fct_unique(df_plot2$dept_fac))
legend_lvls <- c(setdiff(legend_lvls, "Other"), "Other")

# パレットを凡例レベルに合わせて並べ替え（存在しない科は落とす）
local_palette <- color_palette_filtered[legend_lvls]
# 万一パレットにないキーがあればグレーで埋める
need_fill <- is.na(local_palette)
local_palette[need_fill] <- "#BEBEBE"  # フォールバック
df_rect_top <- df_plot2 %>%
  group_by(jin_label) %>%
  # Other を必ず最初に並べて → 上から積み上げの一番上になる
  arrange(dept == "Other", desc(pct), .by_group = TRUE) %>%
  mutate(
    cum_above = cumsum(pct),
    ymax = 100 - lag(cum_above, default = 0),  # 上から積み上げ開始位置
    ymin = ymax - pct
  ) %>%
  ungroup() %>%
  mutate(
    grp_id = as.integer(forcats::fct_inorder(jin_label)),
    width  = 0.45,
    xmin   = grp_id - width,
    xmax   = grp_id + width,
    x      = grp_id,
    ymid   = (ymin + ymax) / 2,
    dept_fac = factor(dept, levels = legend_lvls)
  )

ggplot(df_rect_top) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = dept_fac)) +
  geom_text(aes(x = x, y = ymid, label = ifelse(pct >= 5, sprintf("%.1f", pct), "")),
            size = 3, color = "black") +
  scale_x_continuous(
    breaks = sort(unique(df_rect_top$x)),
    labels = levels(forcats::fct_inorder(df_rect_top$jin_label))
  ) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(
    values = local_palette,
    breaks = legend_lvls,
    name   = "Clinical Department"
  ) +
  labs(x = NULL, y = "Ratio (%)",
       title = "Department composition by group (top stacking, Other always on top)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

