library(readr)
library(dplyr)
library(survival)
library(broom)
library(ggplot2)
library(forestploter)
library(grid)

# CSVファイルをtibbleとして読み込む_藤倉用
setwd("X:/R")
jin1_Eligibile <- read_csv("jin1_Eligibile.csv", locale = locale(encoding = "SHIFT-JIS"))
jin1_Eligibile_cox_3group <- read_csv("jin1_Eligibile_cox_3group.csv", locale = locale(encoding = "SHIFT-JIS"))
colnames(jin1_Eligibile)
# 1人1行（index_date当日レコードのみ・最初のindex_dateを採用）
jin1_Eligibile_unique_sub <- jin1_Eligibile %>%
  group_by(id) %>%
  arrange(index_date, date, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  # nonAKD と AKDのみ（あなたの方針に合わせて必要ならここは保持）
  filter(jin_status %in% c("nonAKD", "AKD")) %>%
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

# --- CKDの2値化（不明"nd"は除外。含めたい場合は0/1へ再符号化を検討） ---
dat <- jin1_Eligibile_unique_sub %>%
  mutate(
    ckd_bin = case_when(
      CKD_status == "CKD"    ~ 1L,
      CKD_status == "nonCKD" ~ 0L,
      TRUE ~ NA_integer_
    )
  ) %>%
  filter(!is.na(ckd_bin), !is.na(jin_label), !is.na(time_years), time_years >= 0)

# サブグループのn・イベント数を把握
dat %>% count(ckd_bin, jin_label, name = "n") %>% print()
dat %>% group_by(ckd_bin, jin_label) %>% summarise(events = sum(primary_death), n = n(), .groups="drop") %>% print()

# --- サブグループ別のCox（CKDはモデルから外す：サブセット内で不変だから） ---
cox_formula_sub <- as.formula(
  Surv(time_years, primary_death) ~
    jin_label + age + index_cre + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15
)

fit_ckd1 <- coxph(cox_formula_sub, data = filter(dat, ckd_bin == 1))  # CKDあり
fit_ckd0 <- coxph(cox_formula_sub, data = filter(dat, ckd_bin == 0))  # CKDなし

# PH仮定チェック（必要に応じて）
# cox.zph(fit_ckd1); cox.zph(fit_ckd0)

# --- 群間HR（Recovery/Non-Recovery vs nonAKD）だけを取り出して整形 ---
tidy_sub <- function(fit, label){
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% c("jin_labelRecovery", "jin_labelNon-Recovery")) %>%
    transmute(
      Subgroup = label,
      Comparison = recode(term,
                          "jin_labelRecovery" = "Recovery vs nonAKD",
                          "jin_labelNon-Recovery" = "Non-Recovery vs nonAKD"),
      HR = estimate,
      CI_low = conf.low,
      CI_high = conf.high,
      p = p.value
    )
}

tbl_ckd <- bind_rows(
  tidy_sub(fit_ckd1, "CKD: yes"),
  tidy_sub(fit_ckd0, "CKD: no")
)

print(tbl_ckd, n = Inf)


cox_interaction <- coxph(
  Surv(time_years, primary_death) ~ 
    jin_label * CKD_status +
    age + index_cre + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
  data = jin1_Eligibile_cox_3group
)

summary(cox_interaction)

# --- サブグループ・フォレスト ---
{
# 並び順を明示（上：CKDあり → 下：CKDなし）
plot_df <- tbl_ckd %>%
  mutate(
    Comparison = factor(Comparison,
                        levels = c("Recovery vs nonAKD","Non-Recovery vs nonAKD")),
    Subgroup   = factor(Subgroup, levels = c("CKD: yes","CKD: no"))
  )

# CKDありの2本（Recovery/Non-Recovery）を上段、
# CKDなしの2本を下段に配置
ggplot(plot_df,
       aes(x = HR, y = Comparison, xmin = CI_low, xmax = CI_high)) +
  geom_point(size = 3) +
  geom_errorbarh(height = 0.18) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  facet_grid(rows = vars(Subgroup), switch = "y") +  # ← ここで上下に分ける
  labs(
    title = "Death (Subgroup Cox): AKD groups within CKD strata",
    x = "Hazard Ratio (log scale)", y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold") # 見出しを左に
  )

fit_total <- coxph(cox_formula_sub, data = dat)  
# 基準群（nonAKD）行を追加してデータを再構築
total_header <- tibble(
  Subgroup = "Total", 
  Comparison = "", 
  HR = NA_real_, 
  CI_low = NA_real_, 
  CI_high = NA_real_, 
  p = NA_real_
)

total_data <- bind_rows(
  tibble(Subgroup = " ", Comparison = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, p = NA_real_),
  tidy_sub(fit_total, "Total") %>%
    mutate(Subgroup = " ")
)

ckd_header <- tibble(
  Subgroup = "CKD", 
  Comparison = "", 
  HR = NA_real_, 
  CI_low = NA_real_, 
  CI_high = NA_real_, 
  p = NA_real_
)

# CKD yesグループ
ckd1_data <- bind_rows(
  tibble(Subgroup = "  yes", Comparison = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, p = NA_real_),
  tidy_sub(fit_ckd1, "CKD: yes") %>%
    mutate(Subgroup = case_when(
      row_number() == 0 ~ "  yes",  # 最初の行のみ"yes"
      TRUE ~ "  "                   # 2行目以降は空白
    ))
)

# CKD noグループ
ckd0_data <- bind_rows(
  tibble(Subgroup = "  no", Comparison = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, p = NA_real_),
  tidy_sub(fit_ckd0, "CKD: no") %>%
    mutate(Subgroup = case_when(
      row_number() == 0 ~ "  no",   # 最初の行のみ"no"
      TRUE ~ "  "                   # 2行目以降は空白
    ))
)

# データを結合
tbl_total <- bind_rows(
  total_header,
  total_data,
  ckd_header, 
  ckd1_data,
  ckd0_data
) %>%
  mutate(
    Comparison = case_when(
      Comparison == "" ~ "",
      Comparison == "Recovery vs nonAKD" ~ "Recovery",
      Comparison == "Non-Recovery vs nonAKD" ~ "Non-Recovery",
      Comparison == "nonAKD" ~ "nonAKD",
      TRUE ~ Comparison
    ),
    is_header = Subgroup %in% c("Total", "CKD"),
    is_reference = Comparison == "nonAKD"
  )
# プロット用データフレームの作成
plot_df_forest <- tbl_total %>%
  mutate(
    HR_with_CI = case_when(
      is_header ~ "",  # ヘッダー行は空白
      is_reference ~ "Reference",
      is.na(HR) ~ "",
      TRUE ~ paste0(
        sprintf("%.2f", HR), " (",
        sprintf("%.2f", CI_low), "–", sprintf("%.2f", CI_high), ")"
      )
    ),
    P_value = case_when(
      is_header ~ "",  # ヘッダー行は空白
      is_reference ~ "-",
      is.na(p) ~ "",
      p < 0.001 ~ "<0.001",
      p < 0.01  ~ formatC(p, format = "fg", digits = 1),
      TRUE      ~ formatC(p, format = "fg", digits = 2)
    )
  ) %>%
  dplyr::select(Subgroup, Comparison, HR_with_CI, P_value)

plot_df_forest$hazard <- paste(rep(" ", 20), collapse = " ")
plot_df_forest <- plot_df_forest %>%
  dplyr::select('Subgroup', 'Comparison', 'hazard', 'HR_with_CI', 'P_value')
colnames(plot_df_forest) <- c("Subgroup", "Comparison", " ", "Hazard Ratio (95% CI)", "P-value")


# フォレストプロット作成
forest_plot <- forestploter::forest(
  data = plot_df_forest,
  est = ifelse(tbl_total$is_header, NA, tbl_total$HR),  # ヘッダー行はNA
  lower = ifelse(tbl_total$is_reference | tbl_total$is_header, 
                 ifelse( tbl_total$is_header, NA, tbl_total$HR), 
                 tbl_total$CI_low),
  upper = ifelse(tbl_total$is_reference | tbl_total$is_header, 
                 ifelse(tbl_total$is_header, NA, tbl_total$HR), 
                 tbl_total$CI_high),
  sizes = ifelse(tbl_total$is_header, 0.1, 0.6),
  ci_column = 3,
  is_summary = tbl_total$is_header,
  ref_line = 1,
  x_trans = "log",
  xlim = c(0.5, 10),
  ticks_at = c(0.5, 1, 2, 4, 8),
  arrow_lab = c("Lower", "Higher")
)
# 適宜行間を調整する
convertHeight(forest_plot$heights, "mm", valueOnly = TRUE) 
forest_plot$heights <- rep(unit(8, "mm"), nrow(forest_plot))
}

## =========================================================
## Supplemental Figure 3: Interaction model forest plot
##  - Title + footnotes
##  - Rename group labels: Recovery/Non-Recovery -> AKD with/without Recovery
##  - Save TIFF (600 dpi) + PDF
##  前提: forest_plot がすでに作成済み
## =========================================================

library(dplyr)
library(grid)
library(gridExtra)

# ---- 1) forest_plot 内で表示している表データ（plot_df_forest）も名称統一 ----
# ※ plot_df_forest を作った後に forest_plot を作っているなら、先にこの置換をして forest_plot を作り直すのが理想。
# すでに forest_plot を作った後なら、この置換→forest_plot再作成 を推奨します。
plot_df_forest <- plot_df_forest %>%
  mutate(
    Comparison = case_when(
      Comparison == "Recovery"      ~ "AKD with Recovery",
      Comparison == "Non-Recovery"  ~ "AKD without Recovery",
      TRUE ~ Comparison
    )
  )

# 列名はそのままでOK（必要ならここで変更可能）
# colnames(plot_df_forest) <- c("Subgroup", "Group", " ", "Hazard Ratio (95% CI)", "P-value")

# ---- 2) HRの数値側（tbl_total）も名称統一（行対応をズラさないため）----
tbl_total <- tbl_total %>%
  mutate(
    Comparison = case_when(
      Comparison == "Recovery"      ~ "AKD with Recovery",
      Comparison == "Non-Recovery"  ~ "AKD without Recovery",
      TRUE ~ Comparison
    )
  )

# ★ forest_plot を「置換後データ」で作り直す（これが一番確実）
forest_plot <- forestploter::forest(
  data = plot_df_forest,
  est = ifelse(tbl_total$is_header, NA, tbl_total$HR),
  lower = ifelse(tbl_total$is_reference | tbl_total$is_header,
                 ifelse(tbl_total$is_header, NA, tbl_total$HR),
                 tbl_total$CI_low),
  upper = ifelse(tbl_total$is_reference | tbl_total$is_header,
                 ifelse(tbl_total$is_header, NA, tbl_total$HR),
                 tbl_total$CI_high),
  sizes = ifelse(tbl_total$is_header, 0.1, 0.6),
  ci_column = 3,
  is_summary = tbl_total$is_header,
  ref_line = 1,
  x_trans = "log",
  xlim = c(0.5, 10),
  ticks_at = c(0.5, 1, 2, 4, 8),
  arrow_lab = c("Lower", "Higher")
)

# ---- 3) タイトル + 注釈（footnote）を追加して1枚にまとめる ----
title_grob <- grid::textGrob(
  "Supplemental Figure 3. Adjusted Hazard Ratios by CKD Status (Interaction Model)",
  x = unit(0, "npc"), y = unit(1, "npc"),
  just = c("left", "top"),
  gp = grid::gpar(fontsize = 14, fontface = "bold")
)

footnote_grob <- grid::textGrob(
  "AKD, acute kidney disease. HR, hazard ratio; CI, confidence interval. Two-sided p value <0.05 was considered statistically significant.",
  x = unit(0, "npc"), y = unit(1, "npc"),
  just = c("left", "top"),
  gp = grid::gpar(fontsize = 9)
)

supp_fig3 <- gridExtra::arrangeGrob(
  title_grob,
  forest_plot,
  footnote_grob,
  ncol = 1,
  heights = c(0.08, 1, 0.12)
)

grid::grid.newpage()
grid::grid.draw(supp_fig3)

# ---- helper: grobを少し縮小して描画（横切れ対策）----
draw_scaled <- function(grob, scale = 0.92){
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(width = grid::unit(scale, "npc"),
                                    height = grid::unit(1, "npc"),
                                    just = c("center", "center")))
  grid::grid.draw(grob)
  grid::popViewport()
}


# ---- 4) 保存（TIFF 600dpi + PDF）----
out_tif <- "X:/R/SupplementalFigure3_interaction_forest_CKD.tiff"
out_pdf <- "X:/R/SupplementalFigure3_interaction_forest_CKD_check.pdf"

if (!requireNamespace("ragg", quietly = TRUE)) install.packages("ragg")

# ★ TIFF：横を少し広く＋縮小描画
ragg::agg_tiff(
  filename = out_tif,
  width = 200, height = 160, units = "mm",   # ← 180→200mmに拡張
  res = 600, compression = "lzw"
)
draw_scaled(supp_fig3, scale = 0.92)         # ← 0.90～0.95で調整可
dev.off()

# ★ PDF：横を少し広く＋縮小描画
pdf(out_pdf, width = 8.3, height = 6.3)      # ← 7.1→8.3 inch に拡張
draw_scaled(supp_fig3, scale = 0.92)
dev.off()

# ---- 5) 交互作用モデルの表（hr_int）も名称統一して保存する場合 ----
# ※ hr_int がある場合のみ
if (exists("hr_int")) {
  hr_int2 <- hr_int %>%
    mutate(
      Comparison = case_when(
        Comparison == "Recovery vs nonAKD"     ~ "AKD with Recovery vs nonAKD",
        Comparison == "Non-Recovery vs nonAKD" ~ "AKD without Recovery vs nonAKD",
        Comparison == "Recovery"              ~ "AKD with Recovery",
        Comparison == "Non-Recovery"          ~ "AKD without Recovery",
        TRUE ~ Comparison
      )
    )
  write.csv(hr_int2, file = "X:/R/SupplementalFigure3_interaction_forest_CKD_hr_by_CKD.xlsx_ready.csv", row.names = FALSE)
}









##fit_ckd1のcoxを走らせるとモデル不安定のアラートが出現するので下記対応がお勧めと、geminiの回答
#報告すべき主要な結果: 交互作用が有意でないため、cox_base（交互作用なしのモデル）の結果を主たる結果として報告します。
#これは「AKDの3群分類が死亡に与える影響は、CKDの有無によって統計的に有意に異なるとは言えない」という検定結果を反映した、最も頑健な結論です。
#hr_int_table の扱い: この結果は「記述的な推定値」または「参考値」として扱います。
#交互作用が有意でない以上、CKDあり/なしでハザード比が異なるという積極的な主張はできません。
#論文の補足資料（Supplement）に載せるか、本文中で「交互作用は有意ではなかったが、参考までに各層でのハザード比を示すと...」と記述する程度に留めるのが一般的です。
#サブグループ別解析（fit_ckd1）の扱い: fit_ckd1で出た警告は、データを分割したことでモデルが不安定になった明確な証拠です。
#したがって、このサブグループ別解析の結果（特にfit_ckd1）は信頼性が低く、主要な結果として採用すべきではありません。
#交互作用モデルのアプローチが、この不安定性を回避するためのより優れた方法であったことが示された、と解釈できます。

# ========================================
# Overall + CKD Subgroup Forest Plot
# ========================================1

# Step 1: Overall解析用のCoxモデル（recvsnon_recvsAKD.Rのモデルを再現）
library(readr)
jin1_Eligibile <- read_csv("/Users/tfuji/Dropbox/臨床研究/石野先生/石野先生_practice/rstudio-export_25.12.18/jin1_Eligibile.csv", locale = locale(encoding = "SHIFT-JIS"))

# Overall用データ作成（recvsnon_recvsAKD.Rと同じロジック）
jin1_Eligibile_cox_overall <- jin1_Eligibile %>%
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
    jin_label = factor(jin_label, levels = c("nonAKD", "Recovery", "Non-Recovery")),
    arb_acei_use = if_else(coalesce(arb, 0) == 1 | coalesce(acei, 0) == 1, 1L, 0L),
    time_years = as.numeric(last_follow_death - index_plus_210) / 365.25,
    CKD_status = case_when(
      as.character(CKD_status) %in% c("CKD","1") ~ "CKD",
      as.character(CKD_status) %in% c("nonCKD","0") ~ "nonCKD",
      TRUE ~ as.character(CKD_status)
    ),
    CKD_status = factor(CKD_status, levels = c("nonCKD","CKD"))
  ) %>%
  filter(!is.na(jin_label), !is.na(time_years), time_years >= 0, !is.na(CKD_status))

# Overallモデル（recvsnon_recvsAKD.Rのcox_model_3groupと同じ仕様）
cox_overall <- coxph(
  Surv(time_years, primary_death) ~
    jin_label + age + index_cre + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15 +
    CKD_status,
  data = jin1_Eligibile_cox_overall
)

# Step 2: 患者数・死亡数集計関数
calculate_summary_stats <- function(data, group_var = NULL) {
  if (is.null(group_var)) {
    # Overall
    data %>%
      group_by(jin_label) %>%
      summarise(
        Number = n(),
        Deaths = sum(primary_death),
        Death_pct = sprintf("%.1f", 100 * Deaths / Number),
        .groups = "drop"
      ) %>%
      mutate(
        Subgroup = "Overall",
        Group = as.character(jin_label)
      ) %>%
      select(Subgroup, Group, Number, Deaths, Death_pct)
  } else {
    # Subgroup
    data %>%
      group_by(!!sym(group_var), jin_label) %>%
      summarise(
        Number = n(),
        Deaths = sum(primary_death),
        Death_pct = sprintf("%.1f", 100 * Deaths / Number),
        .groups = "drop"
      ) %>%
      mutate(
        Subgroup = paste("CKD", if_else(!!sym(group_var) == 0, "no", "yes")),
        Group = as.character(jin_label)
      ) %>%
      select(Subgroup, Group, Number, Deaths, Death_pct)
  }
}

# 集計実行
stats_overall <- calculate_summary_stats(jin1_Eligibile_cox_overall)
stats_ckd <- calculate_summary_stats(dat, "ckd_bin")

# Step 3: HR・CI・P値抽出関数
extract_hr_ci_p <- function(fit, ref_level = "nonAKD") {
  tidy_res <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)

  # jin_label関連の係数のみ抽出
  hr_data <- tidy_res %>%
    filter(grepl("^jin_label", term)) %>%
    mutate(
      Group = case_when(
        term == "jin_labelRecovery" ~ "Recovery",
        term == "jin_labelNon-Recovery" ~ "Non-Recovery",
        TRUE ~ NA_character_
      ),
      HR = estimate,
      CI_low = conf.low,
      CI_high = conf.high,
      P_value = p.value
    ) %>%
    select(Group, HR, CI_low, CI_high, P_value)

  # nonAKDのReference行を追加
  bind_rows(
    tibble(Group = ref_level, HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
    hr_data
  )
}

# Overall HR抽出
hr_overall <- extract_hr_ci_p(cox_overall) %>%
  mutate(Subgroup = "Overall")

# CKD Subgroup HR抽出（既存のhr_int_tableの結果を使用、P値も含む）
hr_ckd <- hr_int %>%
  mutate(
    Subgroup = case_when(
      Subgroup == "CKD: no" ~ "CKD no",
      Subgroup == "CKD: yes" ~ "CKD yes",
      TRUE ~ Subgroup
    ),
    Group = case_when(
      Comparison == "Recovery vs nonAKD" ~ "Recovery",
      Comparison == "Non-Recovery vs nonAKD" ~ "Non-Recovery",
      TRUE ~ NA_character_
    )
  ) %>%
  select(Subgroup, Group, HR, CI_low, CI_high, P_value) %>%
  # nonAKD Reference行を各サブグループに追加
  bind_rows(
    tibble(Subgroup = "CKD no", Group = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
    filter(., Subgroup == "CKD no"),
    tibble(Subgroup = "CKD yes", Group = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
    filter(., Subgroup == "CKD yes")
  ) %>%
  distinct()

# Step 4: 統合データフレーム作成
# 統計量とHRを結合
forest_data_full <- bind_rows(
  stats_overall %>% left_join(hr_overall, by = c("Subgroup", "Group")),
  stats_ckd %>% left_join(hr_ckd, by = c("Subgroup", "Group"))
) %>%
  mutate(
    # Subgroupの順序を固定
    Subgroup = factor(Subgroup, levels = c("Overall", "CKD no", "CKD yes")),
    Group = factor(Group, levels = c("nonAKD", "Recovery", "Non-Recovery"))
  ) %>%
  arrange(Subgroup, Group) %>%
  mutate(
    # 第1列目の表示形式変更（Overall/CKD/no/yes のインデント付き）
    `Subgroup` = case_when(
      Subgroup == "Overall" & Group == "nonAKD" ~ "Overall",
      Subgroup == "CKD no" & Group == "nonAKD" ~ "CKD",
      Subgroup == "CKD no" & Group == "Recovery" ~ "  no",
      Subgroup == "CKD yes" & Group == "nonAKD" ~ "CKD",
      Subgroup == "CKD yes" & Group == "Recovery" ~ "  yes",
      TRUE ~ ""
    ),

    # 患者数と死亡数の列（右揃え用にスペース追加は後で）
    `Number` = as.character(Number),
    `No. of death (%)` = paste0(Deaths, " (", Death_pct, ")"),

    # HR (95%CI)の列
    `HR (95%CI)` = if_else(
      Group == "nonAKD",
      "Reference",
      sprintf("%.2f (%.2f–%.2f)", HR, CI_low, CI_high)
    ),

    # P値の列（有意水準で表記）
    `P-value` = case_when(
      Group == "nonAKD" ~ "",
      is.na(P_value) ~ "",
      P_value < 0.001 ~ "<0.001",
      P_value < 0.01 ~ sprintf("%.3f", P_value),
      TRUE ~ sprintf("%.2f", P_value)
    )
  ) %>%
  select(Subgroup, Group, Number, `No. of death (%)`, `HR (95%CI)`, `P-value`, HR, CI_low, CI_high)

# フォレストプロット用の空白列追加（Forest Plotの位置）
forest_data_full$` ` <- paste(rep(" ", 20), collapse = " ")
forest_data_plot <- forest_data_full %>%
  select(Subgroup, Group, Number, `No. of death (%)`, ` `, `HR (95%CI)`, `P-value`)  # 列順: Forest Plotが5列目

# Step 5: forestploterでプロット作成
library(forestploter)
library(grid)

# フォレストプロット作成
forest_plot_combined <- forestploter::forest(
  data = forest_data_plot,
  est = forest_data_full$HR,
  lower = forest_data_full$CI_low,
  upper = forest_data_full$CI_high,
  sizes = 0.6,
  ci_column = 5,  # 空白列の位置
  ref_line = 1,
  x_trans = "log",
  xlim = c(0.5, 10),
  ticks_at = c(0.5, 1, 2, 4, 8),
  arrow_lab = c("Favors AKD", "Favors nonAKD")
)

# テキストフォーマット設定
# Number列とNo. of death列を右揃えに
forest_plot_combined <- edit_plot(
  forest_plot_combined,
  col = c(3, 4),  # Number列とNo. of death列
  which = "text",
  hjust = unit(1, "npc"),  # 右揃え
  x = unit(1, "npc")
)

# ヘッダー下に横線を追加
forest_plot_combined <- add_border(
  forest_plot_combined,
  row = 0,  # ヘッダー行
  where = "bottom",
  gp = gpar(lwd = 1)
)

# 行間調整
forest_plot_combined$heights <- rep(unit(8, "mm"), nrow(forest_plot_combined))

# 表示
plot(forest_plot_combined)

# Step 6: 論文用保存（TIFF 600dpi）
setwd("E:/R")

# TIFFで保存（余白を最小化）
tiff(
  filename = "Supplementary_Figure_ForestPlot_Overall_CKD_Subgroup.tiff",
  width = 240,   # mm
  height = 110,  # mm
  units = "mm",
  res = 600,
  compression = "lzw"
)
# 余白を最小化（左右0.5%、上1%、下0.5%）
grid.newpage()
pushViewport(viewport(x = unit(0.005, "npc"), y = unit(0.005, "npc"),
                      width = unit(0.99, "npc"), height = unit(0.985, "npc"),
                      just = c("left", "bottom")))
grid.draw(forest_plot_combined)
popViewport()
dev.off()

# PDFでも保存（同じ設定で）
cairo_pdf(
  filename = "Supplementary_Figure_ForestPlot_Overall_CKD_Subgroup.pdf",
  width = 240 / 25.4,  # インチ変換（240mm）
  height = 110 / 25.4  # インチ変換（110mm）
)
# 余白を最小化（左右0.5%、上1%、下0.5%）
grid.newpage()
pushViewport(viewport(x = unit(0.005, "npc"), y = unit(0.005, "npc"),
                      width = unit(0.99, "npc"), height = unit(0.985, "npc"),
                      just = c("left", "bottom")))
grid.draw(forest_plot_combined)
popViewport()
dev.off()