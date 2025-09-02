library(dplyr)
library(tibble)
library(readr)
library(lubridate)
#jin_label作成####
library(readr)
library(dplyr)

# 読み込み設定
loc <- locale(encoding = "SHIFT-JIS")
setwd("E:/R")
# 必要な列だけ読み込み（.default = "c" は他列を文字列扱いにして後で選択）
cre_2012 <- read_csv("jin/cre_over18/cre_2012_over18.csv",
                     locale = loc, skip = 3,
                     col_types = cols(.default = "c")) %>%
  select(患者ID, 科ｺｰﾄﾞ, 検査日)
cre_2013 <- read_csv("jin/cre_over18/cre_2013_over18.csv",
                     locale = loc, skip = 3,
                     col_types = cols(.default = "c")) %>%
  select(患者ID, 科ｺｰﾄﾞ, 検査日)
cre_2014 <- read_csv("jin/cre_over18/cre_2014_over18.csv",
                     locale = loc, skip = 3,
                     col_types = cols(.default = "c")) %>%
  select(患者ID, 科ｺｰﾄﾞ, 検査日)

# 行結合
cre_2012_2014 <- bind_rows(cre_2012, cre_2013, cre_2014)
# 検査日を日付型に変換（必要に応じて）
cre_2012_2014 <- cre_2012_2014 %>%
  mutate(検査日 = as.Date(as.character(検査日), format = "%Y%m%d"))

cre_2012_2014_sub <- cre_2012_2014 %>%
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
cre_2012_2014_sub
library(readr)
setwd("E:/R")
jin1_Eligibile <- read_csv("jin1_Eligibile.csv", locale = locale(encoding = "SHIFT-JIS"))
#判定期間にデータがないものを除外#####
jin1_Eligibile_include_code <- jin1_Eligibile %>%
  left_join(cre_2012_2014_sub, by = c("id", "date")) %>%
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
  mutate(jin_label = factor(jin_label, levels = c("nonAKD", "Non-Recovery", "Recovery", "No-data")))
library(dplyr)
jin1_Eligibile_include_code %>%
  group_by(jin_label) %>%
  summarise(unique_ids = n_distinct(id), .groups = "drop")
colnames(jin1_Eligibile_include_code)

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
  jin_label = c("nonAKD", "Non-Recovery", "Recovery", "No-data"),
  total_patients = c(14406, 100, 114, 127)
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
library(forcats)

# --- 1) 上位5位（同率はすべて採用）＋Other ---
main_code_ratio_top5 <- main_code_ratio_en %>%
  filter(!is.na(jin_label), !is.na(clinical_department_en), !is.na(pct)) %>%
  group_by(jin_label) %>%
  mutate(rank = dense_rank(desc(pct))) %>%
  mutate(clinical_department_en = if_else(rank <= 5, clinical_department_en, "Other")) %>%
  ungroup() %>%
  summarise(pct = sum(pct), .by = c(jin_label, clinical_department_en)) %>%
  group_by(jin_label) %>%
  mutate(pct = 100 * pct / sum(pct)) %>%
  ungroup()

# --- 2) プロット用データと凡例レベル（Otherを最後） ---
df_plot2 <- main_code_ratio_top5 %>%
  transmute(jin_label, dept = clinical_department_en, pct)

dept_levels_all <- df_plot2 %>% distinct(dept) %>% pull(dept)
dept_levels_all <- c(setdiff(dept_levels_all, "Other"), "Other")

# --- 3) 既定パレット + 未定義科を自動着色（Otherは常にグレー） ---
base_palette <- c(
  "Digestive Surgery"         = "#1B9E77",
  "Endocrinology"             = "#E41A1C",
  "Hematology"                = "#377EB8",
  "Obstetrics and Gynecology" = "#E78AC3",
  "Urology"                   = "#FFA07A",
  "Cardiology"                = "#984EA3",
  "Cardiovascular Surgery"    = "#66C2A5",
  "Emergency Medicine"        = "#FF7F00",
  "Otorhinolaryngology"       = "#FFD700",
  "Other"                     = "#BEBEBE"  # ←固定
)

# 未定義キー（Other以外）を抽出
missing_keys <- setdiff(dept_levels_all, names(base_palette))
missing_keys <- setdiff(missing_keys, "Other")
n_missing <- length(missing_keys)

# 均等色相の自動配色（再現性あり・見分けやすいトーン）
auto_cols <- if (n_missing > 0) {
  grDevices::hcl(
    h = seq(15, 375, length.out = n_missing + 1)[1:n_missing],
    c = 60, l = 65
  )
} else character(0)
names(auto_cols) <- missing_keys

# マージして凡例順に並べ替え（Otherは必ずグレーで上書き）
local_palette <- c(base_palette, auto_cols)
local_palette["Other"] <- "#BEBEBE"
local_palette <- local_palette[dept_levels_all]

# --- 4) 上から積み上げ座標計算（Otherを最上段） ---
df_rect_top <- df_plot2 %>%
  mutate(dept_fac = factor(dept, levels = dept_levels_all)) %>%
  group_by(jin_label) %>%
  arrange(dept == "Other", desc(pct), .by_group = TRUE) %>%
  mutate(
    cum_above = cumsum(pct),
    ymax = 100 - dplyr::lag(cum_above, default = 0),
    ymin = ymax - pct
  ) %>%
  ungroup() %>%
  mutate(
    grp_id = as.integer(forcats::fct_inorder(jin_label)),
    width  = 0.45,
    xmin   = grp_id - width,
    xmax   = grp_id + width,
    x      = grp_id,
    ymid   = (ymin + ymax) / 2
  )

# --- 5) 作図（凡例は出現科のみ・drop=FALSEで固定表示） ---
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
    limits = dept_levels_all,
    drop   = FALSE,
    name   = "Clinical Department"
  ) +
  labs(x = NULL, y = "Ratio (%)",
       title = "Department composition by group (Top 5 + Other, stacked from top)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


#産婦人科、泌尿器科にあたるidを抜き出し#####
library(dplyr)
library(stringr)
library(tidyr)
library(readr)

# --- id×科（英語名）へ展開 ---
id_dept_en <- jin1_Eligibile_include_code %>%
  distinct(id, jin_label, main_code) %>%
  filter(!is.na(main_code)) %>%
  mutate(code_part = str_split(main_code, " & ")) %>%  # "A & B" → c("A","B")
  unnest(code_part) %>%
  left_join(department_lookup_en, by = c("code_part" = "main_code")) %>%
  rename(dept_en = clinical_department_en)
View(jin1_Eligibile_include_code)

# --- ① No-data かつ Obstetrics and Gynecology のID一覧 ---
no_data_obgyn_ids <- id_dept_en %>%
  filter(jin_label == "No-data", dept_en == "Obstetrics and Gynecology") %>%
  distinct(id) %>%
  arrange(id)

# --- ② Non-Recovery かつ Urology のID一覧 ---
nonrec_urology_ids <- id_dept_en %>%
  filter(jin_label == "Non-Recovery", dept_en == "Urology") %>%
  distinct(id) %>%
  arrange(id)

# 確認出力
cat("# No-data × Obstetrics and Gynecology: ", nrow(no_data_obgyn_ids), "人\n", sep = "")
print(no_data_obgyn_ids, n = Inf)

cat("\n# Non-Recovery × Urology: ", nrow(nonrec_urology_ids), "人\n", sep = "")
print(nonrec_urology_ids, n = Inf)

# 必要ならCSVに保存
write_csv(no_data_obgyn_ids, "idlist_NoData_ObstetricsGynecology.csv")
write_csv(nonrec_urology_ids, "idlist_NonRecovery_Urology.csv")


library(dplyr)

ids_no_data_MM <- jin1_Eligibile_include_code %>%
  filter(jin_label == "No-data", main_code == "MM") %>%
  distinct(id) %>%
  arrange(id)

# 結果を確認
print(ids_no_data_MM, n = Inf)

# 必要ならCSVに保存
# write.csv(ids_no_data_MM, "ids_NoData_MM.csv", row.names = FALSE)

#AKI/AKD日の日にち#####
library(readr)
library(dplyr)
library(tibble)
library(lubridate)
jin1_AKI_date_nonNA_unique <- read_csv("E:/R/jin1_AKI_date_nonNA_unique.csv")
jin1_AKD_date_nonNA_unique <- read_csv("E:/R/jin1_AKD_date_nonNA_unique.csv")
jin1_Eligibile_AKD_code <- jin1_Eligibile %>%
  # 既存のコード（cre_2012_2013_sub を id, date で突合）
  left_join(cre_2012_2014_sub, by = c("id", "date")) %>%
  # ★ AKI / AKD の発生日を id で突合（必要列だけに絞るのが安全）
  #left_join(select(jin1_AKI_date_nonNA_unique, id, AKI_date), by = "id") %>%
  #left_join(select(jin1_AKD_date_nonNA_unique, id, AKD_date), by = "id") %>%
  # ★ 型合わせ（文字なら Date に変換）
  mutate(
    AKI_date = as.Date(AKI_date),
    AKD_date = as.Date(AKD_date)
    # もし 20130715 のような数値/文字なら lubridate::ymd(AKI_date) 等に変更
  ) %>%
  # ★ event_date を最小（より早い日付）で作成（片方NAも考慮）
  mutate(
    event_date = if_else(
      !is.na(AKI_date) & !is.na(AKD_date),
      pmin(AKI_date, AKD_date),
      coalesce(AKI_date, AKD_date)
    )
  ) %>%
  # 以降は元のラベル付け
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
  mutate(jin_label = factor(jin_label, levels = c("nonAKD", "Non-Recovery", "Recovery", "No-data")))

jin1_Eligibile_AKD_code_unique <- jin1_Eligibile_AKD_code %>%
  filter(jin_label != "nonAKD") %>%
  distinct(id, date, jin_label, event_date, code)
jin1_Eligibile_AKD_code_event_only <- jin1_Eligibile_AKD_code_unique %>%
  filter(date == event_date)
# 1) 突合して科名を付与
jin1_Eligibile_AKD_code_event_only <- jin1_Eligibile_AKD_code_event_only %>%
  left_join(department_lookup_jp,  by = c("code" = "main_code")) %>%
  left_join(department_lookup_en,  by = c("code" = "main_code"))
# 2) jin_statusごとにcodeを集計
code_summary <- jin1_Eligibile_AKD_code_event_only %>%
  count(jin_label, code, clinical_department, clinical_department_en, name = "n") %>%
  arrange(jin_label, desc(n))
code_summary_top10 <- code_summary %>%
  group_by(jin_label) %>%
  slice_max(order_by = n, n = 10) %>%
  ungroup()
print(code_summary_top10, n=Inf)

library(dplyr)
library(ggplot2)
library(forcats)

# --- 1) Top 5 + Other に分類 ---
df_plot_top5_other <- code_summary %>%
  group_by(jin_label) %>%
  mutate(rank = dense_rank(desc(n))) %>%
  mutate(
    dept = if_else(
      rank <= 5,
      dplyr::coalesce(clinical_department_en, "Unknown"),
      "Other"
    )
  ) %>%
  ungroup() %>%                                  # ★ ここで解除
  summarise(n = sum(n), .by = c(jin_label, dept)) %>%  # ★ .by は非グループ化で
  group_by(jin_label) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

# --- 2) 凡例レベル（Otherは最後） ---
dept_levels_all <- df_plot_top5_other %>% distinct(dept) %>% pull(dept)
dept_levels_all <- c(setdiff(dept_levels_all, "Other"), "Other")

# --- 3) パレット（Otherはグレー固定） ---
base_palette <- c(
  "Digestive Surgery"         = "#1B9E77",
  "Endocrinology"             = "#E41A1C",
  "Hematology"                = "#377EB8",
  "Obstetrics and Gynecology" = "#E78AC3",
  "Urology"                   = "#FFA07A",
  "Cardiology"                = "#984EA3",
  "Cardiovascular Surgery"    = "#66C2A5",
  "Emergency Medicine"        = "#FF7F00",
  "Otorhinolaryngology"       = "#FFD700",
  "Neurosurgery"              = "#8DA0CB",
  "Gastroenterology"          = "#A6D854",
  "Respiratory Medicine"      = "#BC80BD",
  "Nephrology"                = "#80B1D3",
  "Orthopedics"               = "#FDB462",
  "Vascular Surgery"          = "#FB9A99",
  "Immunology"                = "#B3DE69",
  "Other"                     = "#BEBEBE"
)

missing_keys <- setdiff(dept_levels_all, names(base_palette))
missing_keys <- setdiff(missing_keys, "Other")
n_missing <- length(missing_keys)
auto_cols <- if (n_missing > 0) {
  grDevices::hcl(
    h = seq(15, 375, length.out = n_missing + 1)[1:n_missing],
    c = 60, l = 65
  )
} else character(0)
names(auto_cols) <- missing_keys

local_palette <- c(base_palette, auto_cols)
local_palette["Other"] <- "#BEBEBE"
local_palette <- local_palette[dept_levels_all]

# --- 4) 上から積み上げ（Otherを最上段） ---
df_rect_top <- df_plot_top5_other %>%
  mutate(dept_fac = factor(dept, levels = dept_levels_all)) %>%
  group_by(jin_label) %>%
  arrange(dept == "Other", desc(pct), .by_group = TRUE) %>%
  mutate(
    cum_above = cumsum(pct),
    ymax = 100 - dplyr::lag(cum_above, default = 0),
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
    label_pct = ifelse(pct >= 5, sprintf("%.1f%%", pct), "")
  )

# --- 5) 作図 ---
x_lut <- df_rect_top %>% distinct(grp_id = x, jin_label) %>% arrange(grp_id)

ggplot(df_rect_top) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = dept_fac)) +
  geom_text(aes(x = x, y = ymid, label = label_pct), size = 3, color = "black") +
  scale_x_continuous(
    breaks = x_lut$grp_id,
    labels = as.character(x_lut$jin_label)
  ) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(
    values = local_palette,
    limits = dept_levels_all,
    drop   = FALSE,
    name   = "Clinical Department"
  ) +
  labs(
    x = NULL, y = "Ratio within group (%)",
    title = "Top-5 departments by group (AKI/AKD event_date)"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
