#必要な情報をすべて含有するcsvファイル(jin1_Eligible? primay_ESKDやdeath情報や感度分析に必要な情報などすべて含有)を作成し、適当なdfを作成
View(jin1_Eligibile)
library(readr)
setwd("X:/R")
jin1_Eligibile <- read_csv("jin1_Eligibile.csv", locale = locale(encoding = "SHIFT-JIS"))
colnames(jin1_Eligibile)
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

#死亡者人数を集計
library(dplyr)

jin1_Eligibile_death %>%
  group_by(jin_label) %>%
  summarise(
    n_primary_death = sum(primary_death == 1, na.rm = TRUE)
  )


##論文用の図（デーブル付き）####
# 2カラム用の例：幅 180 mm ≒ 7.1 inch
setwd("X:/R")
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


# CSV保存
library(readr)
write_csv(
  jin1_Eligibile_cox_3group,
  "jin1_Eligibile_cox_3group.csv"   # 保存ファイル名（作業ディレクトリに保存されます）
) 
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
  setwd("X:/R")  
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
    arb_acei_use = if_else(coalesce(arb, 0) == 1 | coalesce(acei, 0) == 1, 1L, 0L),
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
setwd("X:/R")
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

#cox比例ハザード###
library(survival)
library(dplyr)
library(broom)

# 参照群を nonAKD に設定
jin1_Eligibile_sens <- jin1_Eligibile_sens %>%
  mutate(group = relevel(group, ref = "nonAKD"))

# 本解析と同一の共変量セット
cox_model_3group_sens <- coxph(
  Surv(time_years, primary_death) ~
    group + age + index_cre + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15 +
    CKD_status,   # ref = nonCKD
  data = jin1_Eligibile_sens
)

summary(cox_model_3group_sens)

anova(cox_model_3group_sens, test = "LRT")

ph_test <- cox.zph(cox_model_3group_sens)
print(ph_test)

#AKD群を比較する
jin1_Eligibile_sens2 <- jin1_Eligibile_sens %>%
  mutate(group = relevel(group, ref = "Recovery"))

cox_model_recovery_ref <- coxph(
  Surv(time_years, primary_death) ~
    group + age + index_cre + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 +
    dn10 + dn12 + dn12 + dn13 + dn14 + dn15 +
    CKD_status,
  data = jin1_Eligibile_sens2
)

summary(cox_model_recovery_ref)


##併存疾患表の作成#####
# ============================================
# Supplementary Table 1
# 元Excelの順序を保持したまま term を collapse
# → Word（docx）に論文掲載体裁で出力
# ============================================

library(readxl)
library(dplyr)
library(stringr)
library(flextable)
library(officer)
library(openxlsx)

# ---- 1) read ----
file <- "Supplementary_Table_1_combined.xlsx"

dat <- readxl::read_xlsx(file, sheet = "Supplementary_Table_1") %>%
  rename(
    term = `Drug class or comorbidity`,
    definition = `Definitions (generic drug names or ICD-10 codes)`
  ) %>%
  mutate(
    term = str_squish(term),
    definition = str_squish(definition),
    row_id = row_number()   # ★ 元Excelでの出現順を保存
  ) %>%
  filter(!is.na(term), !is.na(definition), term != "", definition != "") %>%
  filter(!(term %in% c("—","-","--") & definition %in% c("—","-","--")))

# ---- 2) termごとにまとめる（最初に出た順序を保持） ----
dat_collapsed <- dat %>%
  group_by(term) %>%
  summarise(
    definition = paste(unique(definition), collapse = "; "),
    first_row  = min(row_id),     # ★ 最初に出た位置
    .groups = "drop"
  ) %>%
  arrange(first_row) %>%          # ★ 元Excel順に並び替え
  mutate(
    definition = stringr::str_replace_all(definition, ";\\s*", ";\n"),
    definition = stringr::str_replace_all(definition, ",\\s*", ", ")
  ) %>%   
  dplyr::select(term, definition)
 


# ---- 3) flextable（論文掲載体裁） ----
caption_txt <- "Supplementary Table 1. Definitions of medications and comorbidities"

ft <- flextable(dat_collapsed) %>%
  set_caption(caption_txt) %>%
  set_header_labels(
    term = "Drug class or comorbidity",
    definition = "Definitions (generic drug names or ICD-10 codes)"
  ) %>%
  bold(part = "header") %>%
  align(align = "left", part = "all") %>%
  valign(valign = "top", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  width(j = "term", width = 2.3) %>%
  width(j = "definition", width = 5.0) %>%
  border_remove() %>%
  hline_top(border = fp_border(width = 1)) %>%
  hline(border = fp_border(width = 0.6), part = "header") %>%
  hline_bottom(border = fp_border(width = 1)) %>%
  autofit() %>%
# ---- 罫線設定 ----
border_remove() %>%
  
  # 表の一番上（太線）
  hline_top(border = fp_border(width = 1)) %>%
  
  # ヘッダ下（中太線）
  hline(border = fp_border(width = 0.8), part = "header") %>%
  
  # ★ term（各行）の下に細い罫線を引く
  hline(
    i = seq_len(nrow(dat_collapsed)),
    border = fp_border(width = 0.4),
    part = "body"
  ) %>%
  
  # 表の一番下（太線）
  hline_bottom(border = fp_border(width = 1))
# ft を作ったあとに追加する（重要）
ft <- ft %>%
  autofit() %>%
  flextable::set_table_properties(
    layout = "autofit",
    width  = 1        # ← Wordページ幅に強制フィット
  ) %>%
  flextable::valign(valign = "top", part = "all")


# ---- 4) Wordに出力 ----
doc <- read_docx()
doc <- body_add_flextable(doc, value = ft)
print(doc, target = "Supplementary_Table_1_collapsed.docx")

# ---- 5) 掲載用Excelも保存 ----
openxlsx::write.xlsx(
  dat_collapsed,
  "Supplementary_Table_1_collapsed.xlsx"
)

# ---- 6) Viewer確認 ----
ft



# 本解析Recovery ID
id_main_rec <- jin1_Eligibile_death %>%
  filter(jin_label == "nonAKD") %>%
  pull(id) %>% unique()

# 感度分析Recovery ID
id_sens_rec <- jin1_Eligibile_sens %>%
  filter(group == "nonAKD") %>%
  pull(id) %>% unique()

# 差分確認
setdiff(id_main_rec, id_sens_rec)   # 本解析にいて感度にいない
setdiff(id_sens_rec, id_main_rec)   # 感度にいて本解析にいない
length(id_main_rec); length(id_sens_rec)
