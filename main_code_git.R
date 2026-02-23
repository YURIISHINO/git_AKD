library(dplyr)
library(tibble)
library(readr)
library(lubridate)
#jin_label作成####
library(readr)
library(dplyr)



# 読み込み設定
loc <- locale(encoding = "SHIFT-JIS")
setwd("X:/R")
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
library(readr)
setwd("X:/R")
jin1_Eligibile <- read_csv("jin1_Eligibile.csv", locale = locale(encoding = "SHIFT-JIS"))
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

# ---- 推奨: 患者数の自動計算版（再現性・透明性向上のため） ----
# 注意: 上記のハードコーディングされた値は、データ更新時に齟齬が生じるリスクがあります。
#       以下のコードで実データから自動計算することを推奨します。
#       理由: (1) 再現性: データ更新時も自動的に正しい数値が計算される
#             (2) 透明性: 患者数の根拠が明確(実データから計算)
#             (3) エラー防止: 手入力による数値ミスを防ぐ
#             (4) 監査対応: 論文査読時に数値の根拠を明示できる
#
# label_totals <- jin1_Eligibile_include_code %>%
#   distinct(id, jin_label) %>%          # 1人1行に集約
#   count(jin_label, name = "total_patients") %>%
#   arrange(desc(total_patients))
# print("Calculated patient counts by jin_label:")
# print(label_totals)
# ---- 推奨コード終わり ----

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


{
  # ============================================================
  # Standardized residual heatmap (AKD only) among top 10 departments
  # One-shot paste-and-run
  # Requires objects in memory:
  #   - jin1_Eligibile_include_code
  #   - department_lookup_en (main_code <-> clinical_department_en)
  # ============================================================
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(scales)
  })
  
  # --------------------------
  # 0) Safety checks
  # --------------------------
  if (!exists("jin1_Eligibile_include_code")) {
    stop("Object 'jin1_Eligibile_include_code' not found in environment.")
  }
  if (!exists("department_lookup_en")) {
    stop("Object 'department_lookup_en' not found in environment.")
  }
  
  required_cols <- c("id", "exclude", "jin_status")
  missing_cols <- setdiff(required_cols, colnames(jin1_Eligibile_include_code))
  if (length(missing_cols) > 0) {
    stop(paste0("Missing required columns in jin1_Eligibile_include_code: ",
                paste(missing_cols, collapse = ", ")))
  }
  
  # --------------------------
  # 1) Build 1-row-per-id dataset + AKD flag
  # --------------------------
  dat_id <- jin1_Eligibile_include_code %>%
    filter(exclude == "include") %>%
    distinct(id, .keep_all = TRUE) %>%
    mutate(akd_flag = if_else(jin_status == "AKD", 1L, 0L))
  
  # --------------------------
  # 2) Create 'dept' safely
  #   Priority:
  #     (a) if dept already exists -> use it
  #     (b) else if clinical_department_en exists -> use it
  #     (c) else if main_code exists -> join lookup table to create dept
  # --------------------------
  if ("dept" %in% colnames(dat_id)) {
    dat_id <- dat_id %>% mutate(dept = as.character(dept))
  } else if ("clinical_department_en" %in% colnames(dat_id)) {
    dat_id <- dat_id %>% mutate(dept = as.character(clinical_department_en))
  } else {
    # Need main_code to join lookup
    if (!"main_code" %in% colnames(dat_id)) {
      stop("No 'dept' or 'clinical_department_en' column found, and 'main_code' is also missing. Cannot define department.")
    }
    # Join lookup
    if (!all(c("main_code", "clinical_department_en") %in% colnames(department_lookup_en))) {
      stop("department_lookup_en must have columns: main_code, clinical_department_en")
    }
    dat_id <- dat_id %>%
      left_join(department_lookup_en, by = "main_code") %>%
      mutate(dept = as.character(clinical_department_en))
  }
  
  # Remove NA/blank departments
  dat_id <- dat_id %>%
    filter(!is.na(dept), dept != "")
  
  # --------------------------
  # 3) Select top 10 departments by AKD count
  # --------------------------
  top10_dept <- dat_id %>%
    group_by(dept) %>%
    summarise(akd_n = sum(akd_flag == 1), .groups = "drop") %>%
    arrange(desc(akd_n)) %>%
    slice_head(n = 10) %>%
    pull(dept)
  
  if (length(top10_dept) < 2) {
    stop("Too few departments after filtering. Check dept definition / missing values.")
  }
  
  dat_top10 <- dat_id %>%
    filter(dept %in% top10_dept) %>%
    mutate(dept = factor(dept, levels = top10_dept))  # top1 at top
  
  # --------------------------
  # 4) Contingency table -> chi-square -> standardized residuals
  # --------------------------
  tab <- with(dat_top10, table(dept, akd_flag))
  # Force column names to nonAKD/AKD
  if (ncol(tab) != 2) {
    stop("Contingency table does not have exactly 2 columns. Check akd_flag coding.")
  }
  colnames(tab) <- c("nonAKD", "AKD")
  
  chi <- suppressWarnings(chisq.test(tab))
  stdres <- as.matrix(chi$stdres)
  
  # Ensure AKD column exists
  if (!"AKD" %in% colnames(stdres)) {
    stop("Standardized residual matrix has no 'AKD' column. Something is wrong with tab column naming.")
  }
  
  # --------------------------
  # 5) Plot data: AKD column only
  # --------------------------
  df_plot <- tibble(
    dept = rownames(stdres),
    std_resid = as.numeric(stdres[, "AKD"])
  ) %>%
    mutate(dept = factor(dept, levels = rev(top10_dept)))
  df_plot <- df_plot %>%
    mutate(dept = forcats::fct_reorder(dept, std_resid))
  
  # --------------------------
  # 6) Plot (match your target style)
  # --------------------------
  p <- ggplot(df_plot, aes(x = "AKD", y = dept, fill = std_resid)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0,
      limits = c(-10, 10),
      oob = scales::squish,
      breaks = c(-10, -5, 0, 5, 10),
      name = "Std. residual\n(AKD)"
    ) +
    labs(
      title = "Standardized residuals\nfor AKD among top 10 departments",
      x = "AKD",
      y = "Clinical department (AKD top 10)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.text.y  = element_text(size = 13),
      axis.text.x  = element_text(size = 13),
      plot.title   = element_text(size = 18, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 13),
      legend.text  = element_text(size = 12)
    )
  
  print(p)
  
  # --------------------------
  # 7) (Optional) show chi-square warning context
  # --------------------------
  cat("\n[Note] chisq.test warning about approximation may occur when expected counts are small.\n")
  cat("       You can inspect expected counts via: chi$expected\n")
  
}




#AKI/AKD日の日にち#####

{
  ############################################################
  # AKD 3群（Recovery / Non-Recovery / No-data）別
  # AKD_date当日の採血オーダー科：縦の積み上げ棒（Top5 + Other）
  #  - 色/凡例は参照コード（base_palette + auto_cols、Other固定）を踏襲
  ############################################################
  
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(purrr)
  
  # -----------------------------
  # 0) Load
  # -----------------------------
  loc <- locale(encoding = "SHIFT-JIS")
  setwd("X:/R")
  
  # ---- cre 2012/2013 -> id×date で code を "A&B" にまとめる ----
  cre_2012 <- read_csv("jin/cre_over18/cre_2012_over18.csv",
                       locale = loc, skip = 3,
                       col_types = cols(.default = "c")) %>%
    select(患者ID, 科ｺｰﾄﾞ, 検査日)
  
  cre_2013 <- read_csv("jin/cre_over18/cre_2013_over18.csv",
                       locale = loc, skip = 3,
                       col_types = cols(.default = "c")) %>%
    select(患者ID, 科ｺｰﾄﾞ, 検査日)
  
  cre_2012_2013_sub <- bind_rows(cre_2012, cre_2013) %>%
    mutate(検査日 = as.Date(as.character(検査日), format = "%Y%m%d")) %>%
    transmute(
      id   = as.numeric(患者ID),
      date = 検査日,
      code = 科ｺｰﾄﾞ
    ) %>%
    group_by(id, date) %>%
    summarise(code = paste(unique(code), collapse = "&"), .groups = "drop")
  
  # ---- jin1_Eligibile 読み込み（ラベル作成用）----
  jin1_Eligibile <- read_csv("jin1_Eligibile.csv", locale = loc)
  
  # -----------------------------
  # 1) AKD 3群ラベル（Recovery / Non-Recovery / No-data）を作る
  #    ※あなたの定義を踏襲（nonAKDはここでは不要なので落とす）
  # -----------------------------
  akd_label_3 <- jin1_Eligibile %>%
    filter(exclude == "include") %>%
    distinct(id, .keep_all = TRUE) %>%
    mutate(
      jin_label = case_when(
        jin_status == "AKD" & `150_210recovery` == 1 ~ "Recovery",
        jin_status == "AKD" & `150_210recovery` == 2 ~ "Non-Recovery",
        jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 1 ~ "Recovery",
        jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 0 ~ "No-data",
        jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 2 ~ "Non-Recovery",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(jin_label)) %>%
    mutate(jin_label = factor(jin_label, levels = c("Non-Recovery", "Recovery", "No-data"))) %>%
    select(id, jin_label)
  
  # -----------------------------
  # 2) AKD_date当日の採血オーダー科を突合し、3群ラベルを付与
  # -----------------------------
  # 消化器外科コードをまとめる
  digestive_codes <- c("DK","DM","DL","DR")
  
  normalize_codes_vec <- function(code_string) {
    if (is.na(code_string) || code_string == "") return(character(0))
    codes <- unlist(str_split(code_string, "&"))
    codes <- str_trim(codes)
    codes[codes %in% digestive_codes] <- "DigestiveSurgery"
    unique(codes)
  }
  
  # 科コード→科名（英語）lookup（あなたの表を踏襲）
  department_lookup_en <- tibble(
    main_code = c(
      "AL", "AQ", "HH", "AP", "AM", "AG", "AN", "AH", "FF",
      "DE", "LL", "AR", "MM", "DQ", "GG", "EE", "DN",
      "DH", "DF", "BB", "PS", "RR", "AK", "CC", "KK", "YY",
      "PP", "NN", "AS", "FS", "DG", "DigestiveSurgery"
    ),
    clinical_department_en = c(
      "Endocrinology","Hematology","Urology","Cardiology","Respiratory Medicine",
      "Gastroenterology","Hepatology","Nephrology","Orthopedics",
      "Cardiovascular Surgery","Otorhinolaryngology","Immunology",
      "Obstetrics and Gynecology","Emergency Medicine","Dermatology","Neurosurgery",
      "Vascular Surgery","Breast Surgery","Thoracic Surgery","Psychiatry",
      "Plastic Surgery","Dentistry and Oral Surgery","Neurology","Pediatrics",
      "Ophthalmology","Administrative Dept.","Anesthesiology","Radiation Therapy",
      "Clinical Pharmacy","Rehabilitation","Pediatric Surgery","Digestive Surgery"
    )
  )
  
  akd_day3 <- jin1_AKD_date_nonNA_unique %>%
    mutate(AKD_date = as.Date(AKD_date)) %>%
    left_join(akd_label_3, by = "id") %>%                       # 3群ラベル付与
    filter(!is.na(jin_label)) %>%
    left_join(cre_2012_2013_sub, by = c("id" = "id", "AKD_date" = "date"))
  
  # 突合状況チェック
  akd_day3 %>%
    group_by(jin_label) %>%
    summarise(
      n_patients = n_distinct(id),
      n_with_code = sum(!is.na(code)),
      n_missing_code = sum(is.na(code)),
      .groups = "drop"
    ) %>% print(n = Inf)
  
  # -----------------------------
  # 3) code分解→科名→患者ベース構成（群別）
  # -----------------------------
  akd_dept_patient3 <- akd_day3 %>%
    mutate(code_vec = purrr::map(code, normalize_codes_vec)) %>%
    unnest(code_vec, keep_empty = TRUE) %>%
    rename(main_code = code_vec) %>%
    left_join(department_lookup_en, by = "main_code") %>%
    mutate(
      clinical_department_en = if_else(is.na(clinical_department_en), "Other", clinical_department_en)
    ) %>%
    distinct(id, jin_label, AKD_date, clinical_department_en)  # 患者×群×当日×科
  
  # 群別に全体比率（%）
  overall3 <- akd_dept_patient3 %>%
    count(jin_label, clinical_department_en, name = "n") %>%
    group_by(jin_label) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup()
  
  # -----------------------------
  # 4) Top5 + Other（群ごとにTop5）
  # -----------------------------
  plot_df <- overall3 %>%
    group_by(jin_label) %>%
    mutate(rank = dense_rank(desc(pct))) %>%
    mutate(dept = if_else(rank <= 5, clinical_department_en, "Other")) %>%
    ungroup() %>%
    group_by(jin_label, dept) %>%
    summarise(pct = sum(pct), .groups = "drop") %>%
    group_by(jin_label) %>%
    mutate(pct = 100 * pct / sum(pct)) %>%
    ungroup()
  
  # ---- 凡例レベル（出現科のみ。Otherを最後）----
  dept_levels_all <- plot_df %>% distinct(dept) %>% pull(dept)
  dept_levels_all <- c(setdiff(dept_levels_all, "Other"), "Other")
  
  plot_df <- plot_df %>%
    mutate(
      jin_label = factor(jin_label, levels = c("Non-Recovery", "Recovery", "No-data")),
      dept_fac  = factor(dept, levels = dept_levels_all)
    )
  
  # -----------------------------
  # 5) 参照コード準拠の配色（base_palette + auto_cols、Other固定グレー）
  # -----------------------------
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
  
  # -----------------------------
  # 6) Plot：縦の積み上げ棒（3本：3群）
  # -----------------------------
  p <- ggplot(plot_df, aes(x = jin_label, y = pct, fill = dept_fac)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = ifelse(pct >= 5, sprintf("%.1f", pct), "")),
              position = position_stack(vjust = 0.5),
              size = 3, color = "black") +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(
      values = local_palette,
      limits = dept_levels_all,
      drop   = FALSE,
      name   = "Clinical Department"
    ) +
    labs(
      x = NULL,
      y = "Ratio (%)",
      title = "Department composition on the AKD date by recovery status (Top 5 + Other)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 15, hjust = 1),
      panel.grid.minor = element_blank()
    )
  
  print(p)
  
  # 保存（論文用）
  ggsave("Figure_AKDdate_department_stacked_by_recovery.png", p, width = 7.0, height = 4.8, dpi = 600)
  ggsave("Figure_AKDdate_department_stacked_by_recovery.tiff", p, width = 7.0, height = 4.8, dpi = 600, compression = "lzw")
  
  
}  

############################################################
# 7) All（3群をまとめて1本）：縦の積み上げ棒（Top5 + Other）
#    ※母集団は「3群ラベルが付与できた患者」に限定（= akd_dept_patient3）
############################################################

# -----------------------------
# A) 全体（All）で科別割合（%）を作る
#    ※akd_dept_patient3 は 3群コードで既に作成済み
#      （id, jin_label, AKD_date, clinical_department_en）
# -----------------------------
overall_all <- akd_dept_patient3 %>%
  count(clinical_department_en, name = "n") %>%
  mutate(pct = 100 * n / sum(n)) %>%
  arrange(desc(pct)) %>%
  mutate(rank = dense_rank(desc(pct))) %>%
  mutate(dept = if_else(rank <= 5, clinical_department_en, "Other")) %>%
  group_by(dept) %>%
  summarise(pct = sum(pct), .groups = "drop") %>%
  mutate(pct = 100 * pct / sum(pct)) %>%   # 念のため100%に再正規化
  arrange(desc(pct))

# ---- 凡例順（Otherを最後）----
dept_levels_all2 <- overall_all %>% pull(dept) %>% unique()
dept_levels_all2 <- c(setdiff(dept_levels_all2, "Other"), "Other")

overall_all <- overall_all %>%
  mutate(dept_fac = factor(dept, levels = dept_levels_all2))

# -----------------------------
# B) 配色（参照コード踏襲：base_palette + auto_cols、Other固定）
# -----------------------------
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
  "Other"                     = "#BEBEBE"
)

missing_keys2 <- setdiff(dept_levels_all2, names(base_palette))
missing_keys2 <- setdiff(missing_keys2, "Other")
n_missing2 <- length(missing_keys2)

auto_cols2 <- if (n_missing2 > 0) {
  grDevices::hcl(
    h = seq(15, 375, length.out = n_missing2 + 1)[1:n_missing2],
    c = 60, l = 65
  )
} else character(0)
names(auto_cols2) <- missing_keys2

local_palette2 <- c(base_palette, auto_cols2)
local_palette2["Other"] <- "#BEBEBE"
local_palette2 <- local_palette2[dept_levels_all2]

# -----------------------------
# C) Plot：縦の積み上げ棒（1本）
# -----------------------------
p_all <- ggplot(overall_all, aes(x = "All", y = pct, fill = dept_fac)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = ifelse(pct >= 5, sprintf("%.1f", pct), "")),
            position = position_stack(vjust = 0.5),
            size = 3, color = "black") +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(
    values = local_palette2,
    limits = dept_levels_all2,
    drop   = FALSE,
    name   = "Clinical Department"
  ) +
  labs(
    x = NULL,
    y = "Ratio (%)",
    title = "Department composition on the AKD date (All AKD; Top 5 + Other)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.minor = element_blank()
  )

print(p_all)

ggsave("Figure_AKDdate_department_stacked_allAKD.png", p_all,
       width = 4.0, height = 5.0, dpi = 600)
ggsave("Figure_AKDdate_department_stacked_allAKD.tiff", p_all,
       width = 4.0, height = 5.0, dpi = 600, compression = "lzw")



{
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  top5_by_group <- akd_dept_patient3 %>%
    count(jin_label, clinical_department_en, name = "n_patients") %>%
    group_by(jin_label) %>%
    arrange(desc(n_patients)) %>%
    slice_head(n = 5) %>%
    ungroup()
  
  # ---- AKD 3群ラベル（あなたの定義を踏襲）----
  akd_label_3 <- jin1_Eligibile %>%
    filter(exclude == "include") %>%
    distinct(id, .keep_all = TRUE) %>%
    mutate(
      jin_label = case_when(
        jin_status == "AKD" & `150_210recovery` == 1 ~ "Recovery",
        jin_status == "AKD" & `150_210recovery` == 2 ~ "Non-Recovery",
        jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 1 ~ "Recovery",
        jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 0 ~ "No-data",
        jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 2 ~ "Non-Recovery",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(jin_label)) %>%
    mutate(jin_label = factor(jin_label,
                              levels = c("Non-Recovery", "Recovery", "No-data"))) %>%
    select(id, jin_label)
  
  # ---- 消化器外科コードまとめ ----
  digestive_codes <- c("DK","DM","DL","DR")
  
  normalize_codes_vec <- function(code_string) {
    if (is.na(code_string) || code_string == "") return(character(0))
    codes <- unlist(str_split(code_string, "&"))
    codes <- str_trim(codes)
    codes[codes %in% digestive_codes] <- "DigestiveSurgery"
    unique(codes)
  }
  
  # ---- AKD_date当日の科（患者×群×科）----
  akd_dept_patient3 <- jin1_AKD_date_nonNA_unique %>%
    mutate(AKD_date = as.Date(AKD_date)) %>%
    left_join(akd_label_3, by = "id") %>%
    filter(!is.na(jin_label)) %>%
    left_join(cre_2012_2013_sub, by = c("id" = "id", "AKD_date" = "date")) %>%
    mutate(code_vec = map(code, normalize_codes_vec)) %>%
    unnest(code_vec, keep_empty = TRUE) %>%
    rename(main_code = code_vec) %>%
    left_join(department_lookup_en, by = "main_code") %>%
    mutate(
      clinical_department_en =
        if_else(is.na(clinical_department_en), "Other", clinical_department_en)
    ) %>%
    distinct(id, jin_label, clinical_department_en)
  group_totals <- akd_dept_patient3 %>%
    distinct(id, jin_label) %>%
    count(jin_label, name = "total_patients")
  
  print(group_totals)
  final_table <- top5_by_group %>%
    left_join(group_totals, by = "jin_label") %>%
    arrange(jin_label, desc(n_patients))
  
  print(final_table, n = Inf)
  
  
  top5_by_group <- akd_dept_patient3 %>%
    count(jin_label, clinical_department_en, name = "n_patients") %>%
    group_by(jin_label) %>%
    arrange(desc(n_patients)) %>%
    slice_head(n = 5) %>%
    ungroup()
  
  
}