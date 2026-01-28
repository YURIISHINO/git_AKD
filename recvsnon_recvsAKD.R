{
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


##論文用の図（デーブル付き）
# 2カラム用の例：幅 180 mm ≒ 7.1 inch
setwd("X:/R")
tiff(
  filename = "Figure2A_death_KM.tiff",
  width    = 7.1,   # inch
  height   = 6.0,   # inch
  units    = "in",
  res      = 600,
  compression = "lzw",
  type     = "cairo"   # ← ここがポイント
)

print(p_death)         # カーブ＋リスクテーブルをまとめて出力

dev.off()

## ---- 図の作成（曲線のみ、抄録仕様）
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

# ---- JPEGで保存（640×480px, ≤1.5MB)
jpeg("death_curve_for_abstract.jpeg", width = 640, height = 480, units = "px", quality = 95)
print(km_plot$plot)
dev.off()

#cox比例ハザード
{
# =========================
# Figure 3: HR table + Forest (AKD groups vs nonAKD)
# Output: Figure3_forest_with_table_death_akd_3groups.tif (+ check PDF)
# =========================

library(dplyr)
library(survival)
library(broom)
library(ggplot2)
library(gridExtra)
library(grid)

# ---- 1) 1人1行データ作成（CKD_status参照=nonCKD）----
dat_cox <- jin1_Eligibile %>%
  filter(exclude == "include", jin_status %in% c("AKD", "nonAKD")) %>%
  group_by(id) %>%
  arrange(index_date, date, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
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
    CKD_status = case_when(
      as.character(CKD_status) %in% c("CKD", "1")     ~ "CKD",
      as.character(CKD_status) %in% c("nonCKD", "0") ~ "nonCKD",
      TRUE ~ as.character(CKD_status)
    ),
    CKD_status = factor(CKD_status, levels = c("nonCKD", "CKD"))
  ) %>%
  filter(!is.na(jin_label), !is.na(time_years), time_years >= 0, !is.na(CKD_status))

# ---- 2) Cox（3群）----
fit_main <- coxph(
  Surv(time_years, primary_death) ~
    jin_label + age + index_cre + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15 +
    CKD_status,
  data = dat_cox
)

# ---- 3) 2比較（vs nonAKD）の tidy（1回だけ）----
td <- broom::tidy(fit_main, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term %in% c("jin_labelRecovery", "jin_labelNon-Recovery")) %>%
  mutate(
    Contrast = dplyr::recode(
      term,
      "jin_labelRecovery"     = "AKD Recovery vs nonAKD",
      "jin_labelNon-Recovery" = "AKD Non-Recovery vs nonAKD"
    ),
    Contrast = factor(Contrast, levels = c("AKD Non-Recovery vs nonAKD", "AKD Recovery vs nonAKD"))
  )

# ---- 4) HRテーブル（上段）----
hr_table <- td %>%
  mutate(
    HR        = round(estimate, 2),
    CI_lower  = round(conf.low, 2),
    CI_upper  = round(conf.high, 2),
    p_raw     = p.value,
    `p-value` = ifelse(p_raw < 0.01, "<0.01", sprintf("%.2f", p_raw))
  ) %>%
  dplyr::select(
    Contrast,
    HR,
    `Lower 95% CI` = CI_lower,
    `Upper 95% CI` = CI_upper,
    `p-value`
  )

# ---- 5) フォレスト（下段）----
p_forest <- ggplot(td, aes(x = Contrast, y = estimate)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, linewidth = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  scale_y_log10(
    name   = "Hazard ratio for all-cause mortality (log scale)",
    breaks = c(0.5, 1, 2, 3, 4, 7),
    limits = c(0.5, 7)
  ) +
  labs(
    x = NULL,
    title = NULL 
  )+
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y      = element_text(size = 10),
    axis.text.x      = element_text(size = 9),
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    plot.margin      = margin(t = 5.5, r = 10, b = 5.5, l = 7, unit = "pt")
  )

# ---- 6) 表を grob 化して「表タイトル + テーブル」にする ----
table_theme <- gridExtra::ttheme_minimal(
  core = list(fontsize = 12, padding = unit(c(2, 2), "mm")),
  colhead = list(fontsize = 12, fontface = "bold", padding = unit(c(2, 2), "mm"))
)

table_grob <- gridExtra::tableGrob(
  hr_table,
  rows  = NULL,
  theme = table_theme
)

# ★ 追加：テーブルタイトル（表の上）
table_title <- grid::textGrob(
  "Adjusted hazard ratios for all-cause mortality (Cox proportional hazards model)",
  x = unit(0, "npc"),  # 左寄せ
  just = c("left", "top"),
  gp = grid::gpar(fontsize = 12, fontface = "bold")
)

# タイトル＋テーブルを縦に結合（タイトルは小さめの高さ）
table_block <- gridExtra::arrangeGrob(
  table_title,
  table_grob,
  ncol = 1,
  heights = c(0.20, 1)
)

# ---- 7) 「表ブロック + フォレスト」を縦に並べる ----
combined <- gridExtra::arrangeGrob(
  table_block,   # ← ここが table_grob から置き換わる
  p_forest,
  ncol    = 1,
  heights = c(1, 2)
)

# ---- 8) 保存（TIFF + check PDF）----
out_tif <- "X:/R/Figure2B_forest_with_table_death_akd_3groups.tif"
out_pdf <- "X:/R/Figure2B_forest_with_table_death_akd_3groups_check.pdf"

# TIFF（Cairoの不安定さ回避のため ragg 推奨）
if (!requireNamespace("ragg", quietly = TRUE)) install.packages("ragg")
ragg::agg_tiff(filename = out_tif, width = 160, height = 140, units = "mm", res = 600, compression = "lzw")
grid::grid.newpage()
grid::grid.draw(combined)
dev.off()

# PDF（確認用）
pdf(out_pdf, width = 7, height = 6)
grid::grid.draw(combined)
dev.off()

}
#本解析primaryまとめ
{
# =========================
# Figure 2: remove original plot titles, add (A)/(B) labels + subtitles
# p_death (ggsurvplot) and combined (Cox table+forest grob) are assumed to exist
# =========================

library(grid)
library(gridExtra)
library(ggplot2)

# ---- 0) Figure2A：KM曲線側の「元タイトル」を削除 ----
# ggsurvplotは list なので $plot に対して theme() を当てる
p_death2 <- p_death
p_death2$plot <- p_death2$plot + theme(plot.title = element_blank())

# （任意）risk table 側も余計なタイトル等がある場合に備えて
if (!is.null(p_death2$table)) {
  p_death2$table <- p_death2$table + theme(plot.title = element_blank())
}

# ---- 1) Figure 2A（KM + risk table）に (A) ラベル＋サブタイトル ----
km_with_label <- gridExtra::arrangeGrob(
  grid::textGrob(
    "(A) Cumulative incidence of all-cause death",
    x = unit(0, "npc"), y = unit(1, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 12, fontface = "bold")
  ),
  gridExtra::arrangeGrob(
    p_death2$plot,
    p_death2$table,
    ncol = 1,
    heights = c(3, 1)
  ),
  ncol = 1,
  heights = c(0.08, 1)
)

# ---- 2) Figure 2B（Cox：table + forest）に (B) ラベル＋サブタイトル ----
# ※ Coxの下段フォレストの元タイトルを消した combined を使うのが前提です
#    もし combined の中身が古くてタイトルが残る場合は、p_forest を title=NULL で作り直して combined も作り直してください
# ---- 7) 「表ブロック + フォレスト」を縦に並べる ----
combined_without_title <- gridExtra::arrangeGrob(
  table_grob,   # ← ここが table_grob から置き換わる
  p_forest,
  ncol    = 1,
  heights = c(1, 2)
)

cox_with_label <- gridExtra::arrangeGrob(
  grid::textGrob(
    "(B) Adjusted hazard ratios for all-cause mortality",
    x = unit(0, "npc"), y = unit(1, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 12, fontface = "bold")
  ),
  combined_without_title,
  ncol = 1,
  heights = c(0.08, 1)
)
# ---- 3') Figure全体タイトルを追加して上下に結合 ----
fig2_AB_labeled <- gridExtra::arrangeGrob(
  # ★ Figure全体タイトル
  grid::textGrob(
    "Figure 2. All-cause mortality",
    x = unit(0, "npc"), y = unit(1, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 14, fontface = "bold")
  ),
  
  # (A) + (B) 本体
  gridExtra::arrangeGrob(
    km_with_label,
    cox_with_label,
    ncol = 1,
    heights = c(1.2, 1.0)
  ),
  
  ncol = 1,
  heights = c(0.06, 1)   # ← タイトル分の高さ
)

# ---- 4) 保存 ----
out_tif <- "X:/R/Figure2_combined_vertical_labeled.tif"
out_pdf <- "X:/R/Figure2_combined_vertical_labeled_check.pdf"

if (!requireNamespace("ragg", quietly = TRUE)) install.packages("ragg")
ragg::agg_tiff(
  filename = out_tif,
  width = 180,
  height = 300,
  units = "mm",
  res = 600,
  compression = "lzw"
)
grid::grid.newpage()
grid::grid.draw(fig2_AB_labeled)
dev.off()

pdf(out_pdf, width = 7.1, height = 11.5)
grid::grid.draw(fig2_AB_labeled)
dev.off()
}

#感度分析を行うのに必要なcodeのみ####
#感度分析死亡のカプランマイヤー####
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
  filename   = "Figure4A_death_KM_sensitivity.tiff",
  width      = 7.1,   # inch（約180 mm）
  height     = 6.0,   # inch（適宜調整）
  units      = "in",
  res        = 600,
  compression = "lzw",
  type       = "cairo"   # ← これが CairoTIFF の代わり
)

print(p_death_sens)      # カーブ + リスクテーブルをまとめて出力

dev.off()

}
# ============================================================
# Figure 2 (ONE-PASTE): Primary Outcome: All-Cause Mortality
#  - (A) KM (event) + risk table
#  - (B) Cox HR table + forest
#  - Labels: nonAKD / AKD with Recovery / AKD without Recovery
#  - Figure title: "Figure 2. Primary Outcome: All-Cause Mortality"
#  - Assumes: jin1_Eligibile.csv exists in X:/R
# ============================================================

library(readr)
library(dplyr)
library(survival)
library(survminer)
library(broom)
library(ggplot2)
library(grid)
library(gridExtra)

setwd("X:/R")

# ----------------------------
# 0) Load data
# ----------------------------
jin1_Eligibile <- read_csv("jin1_Eligibile.csv", locale = locale(encoding = "SHIFT-JIS"))

# ----------------------------
# 1) Make 1-row-per-id dataset (analysis cohort)
# ----------------------------
dat1 <- jin1_Eligibile %>%
  filter(exclude == "include", jin_status %in% c("AKD", "nonAKD")) %>%
  group_by(id) %>%
  arrange(index_date, date, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    # 3-group label for plotting / models
    jin_label = case_when(
      jin_status == "nonAKD" ~ "nonAKD",
      jin_status == "AKD" & `150_210recovery` == 1 ~ "Recovery",
      jin_status == "AKD" & `150_210recovery` == 2 ~ "Non-Recovery",
      jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 1 ~ "Recovery",
      jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` %in% c(0, 2) ~ "Non-Recovery",
      TRUE ~ NA_character_
    ),
    jin_label = factor(jin_label, levels = c("nonAKD", "Recovery", "Non-Recovery")),
    
    # follow-up time (years)
    time_years = as.numeric(last_follow_death - index_plus_210) / 365.25,
    
    # CKD_status harmonization (if needed)
    CKD_status = case_when(
      as.character(CKD_status) %in% c("CKD", "1")     ~ "CKD",
      as.character(CKD_status) %in% c("nonCKD", "0") ~ "nonCKD",
      TRUE ~ as.character(CKD_status)
    ),
    CKD_status = factor(CKD_status, levels = c("nonCKD", "CKD")),
    
    # RASi (ARB/ACEi) combined
    arb_acei_use = if_else(coalesce(arb, 0) == 1 | coalesce(acei, 0) == 1, 1L, 0L)
  ) %>%
  filter(!is.na(jin_label), !is.na(time_years), time_years >= 0)

# ----------------------------
# 2) (A) Kaplan–Meier (event = cumulative incidence)
# ----------------------------
fit_death <- survfit(Surv(time_years, primary_death) ~ jin_label, data = dat1)

p_death <- ggsurvplot(
  fit_death,
  data         = dat1,
  fun          = "event",
  pval         = TRUE,
  conf.int     = TRUE,
  risk.table   = TRUE,
  xlim         = c(0, 9.5),
  ylim         = c(0, 0.4),
  title        = NULL,
  xlab         = "Years",
  ylab         = "Cumulative incidence",
  break.time.by = 2,
  legend.title = "Group",
  legend.labs  = c("nonAKD", "AKD with Recovery", "AKD without Recovery")
)

# remove any internal plot titles (safety)
p_death2 <- p_death
p_death2$plot  <- p_death2$plot  + theme(plot.title = element_blank())
if (!is.null(p_death2$table)) p_death2$table <- p_death2$table + theme(plot.title = element_blank())

km_with_label <- gridExtra::arrangeGrob(
  grid::textGrob(
    "(A) Cumulative incidence of all-cause death",
    x = unit(0, "npc"), y = unit(1, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 12, fontface = "bold")
  ),
  gridExtra::arrangeGrob(
    p_death2$plot,
    p_death2$table,
    ncol = 1,
    heights = c(3, 1)
  ),
  ncol = 1,
  heights = c(0.08, 1)
)

# ----------------------------
# 3) (B) Cox model + HR table + forest
# ----------------------------
fit_main <- coxph(
  Surv(time_years, primary_death) ~
    jin_label + age + index_cre + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15 +
    CKD_status,
  data = dat1 %>% filter(!is.na(CKD_status))
)

td <- broom::tidy(fit_main, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term %in% c("jin_labelRecovery", "jin_labelNon-Recovery")) %>%
  mutate(
    Contrast = dplyr::recode(
      term,
      "jin_labelRecovery"     = "AKD with Recovery vs nonAKD",
      "jin_labelNon-Recovery" = "AKD without Recovery vs nonAKD"
    ),
    Contrast = factor(
      Contrast,
      levels = c("AKD without Recovery vs nonAKD", "AKD with Recovery vs nonAKD")
    )
  )

hr_table <- td %>%
  mutate(
    HR        = round(estimate, 2),
    CI_lower  = round(conf.low, 2),
    CI_upper  = round(conf.high, 2),
    p_raw     = p.value,
    `p-value` = ifelse(p_raw < 0.01, "<0.01", sprintf("%.2f", p_raw))
  ) %>%
  select(
    Contrast,
    HR,
    `Lower 95% CI` = CI_lower,
    `Upper 95% CI` = CI_upper,
    `p-value`
  )

p_forest <- ggplot(td, aes(x = Contrast, y = estimate)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, linewidth = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  scale_y_log10(
    name   = "Hazard ratio for all-cause mortality (log scale)",
    breaks = c(0.5, 1, 2, 3, 4, 7),
    limits = c(0.5, 7)
  ) +
  labs(x = NULL, title = NULL) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y      = element_text(size = 10),
    axis.text.x      = element_text(size = 9),
    plot.margin      = margin(t = 5.5, r = 10, b = 5.5, l = 7, unit = "pt")
  )

table_theme <- gridExtra::ttheme_minimal(
  core = list(fontsize = 12, padding = unit(c(2, 2), "mm")),
  colhead = list(fontsize = 12, fontface = "bold", padding = unit(c(2, 2), "mm"))
)

table_grob <- gridExtra::tableGrob(
  hr_table,
  rows  = NULL,
  theme = table_theme
)

# (B) block (table + forest)
combined_B <- gridExtra::arrangeGrob(
  table_grob,
  p_forest,
  ncol    = 1,
  heights = c(1, 2)
)

cox_with_label <- gridExtra::arrangeGrob(
  grid::textGrob(
    "(B) Adjusted hazard ratios for all-cause mortality",
    x = unit(0, "npc"), y = unit(1, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 12, fontface = "bold")
  ),
  combined_B,
  ncol = 1,
  heights = c(0.08, 1)
)

# ----------------------------
# 4) Combine A + B with FIGURE title
# ----------------------------
fig2_AB_labeled <- gridExtra::arrangeGrob(
  grid::textGrob(
    "Figure 2. Primary Outcome: All-Cause Mortality",
    x = unit(0, "npc"), y = unit(1, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 14, fontface = "bold")
  ),
  gridExtra::arrangeGrob(
    km_with_label,
    cox_with_label,
    ncol = 1,
    heights = c(1.2, 1.0)
  ),
  ncol = 1,
  heights = c(0.06, 1)
)

# ----------------------------
# 5) Save (TIFF + PDF check)
# ----------------------------
out_tif <- "X:/R/Figure2_combined_vertical_labeled.tif"
out_pdf <- "X:/R/Figure2_combined_vertical_labeled_check.pdf"

if (!requireNamespace("ragg", quietly = TRUE)) install.packages("ragg")

ragg::agg_tiff(
  filename = out_tif,
  width = 180,
  height = 300,
  units = "mm",
  res = 600,
  compression = "lzw"
)
grid::grid.newpage()
grid::grid.draw(fig2_AB_labeled)
dev.off()

pdf(out_pdf, width = 7.1, height = 11.5)
grid::grid.draw(fig2_AB_labeled)
dev.off()

# optional preview
grid::grid.newpage()
grid::grid.draw(fig2_AB_labeled)



#感度分析cox比例ハザード#####
{
# =========================
# Sensitivity analysis (Death): Cox PH + Forest with HR table (group 3-level)
# Output: Figure4B_forest_with_table_death_akd_3groups_sensitivity.tif (+ check PDF)
# =========================

library(dplyr)
library(survival)
library(broom)
library(ggplot2)
library(gridExtra)
library(grid)

# ---- 0) Safety: CKD_status ref を nonCKD に固定（本解析と揃える）----
# ※ すでに整形済みならこのmutateは同じ結果になります
jin1_Eligibile_sens <- jin1_Eligibile_sens %>%
  mutate(
    CKD_status = case_when(
      as.character(CKD_status) %in% c("CKD", "1")     ~ "CKD",
      as.character(CKD_status) %in% c("nonCKD", "0") ~ "nonCKD",
      TRUE ~ as.character(CKD_status)
    ),
    CKD_status = factor(CKD_status, levels = c("nonCKD", "CKD"))
  ) %>%
  filter(!is.na(CKD_status))

# ---- 1) Cox model（感度分析：group, ref=nonAKD）----
dat_sens <- jin1_Eligibile_sens %>%
  mutate(group = relevel(group, ref = "nonAKD")) %>%
  droplevels()

fit_sens <- coxph(
  Surv(time_years, primary_death) ~
    group + age + index_cre + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15 +
    CKD_status,
  data = dat_sens
)

# ---- 2) Figure作成関数（引数名は fit で統一）----
make_forest_with_table_2contrast_group <- function(fit,
                                                   out_tif,
                                                   out_pdf,
                                                   width_mm = 160,
                                                   height_mm = 140,
                                                   res_dpi = 600,
                                                   hr_limits = c(0.5, 7),
                                                   hr_breaks = c(0.5, 1, 2, 3, 4, 7),
                                                   title_txt = "Sensitivity analysis: Adjusted hazard ratios for all-cause mortality\n(AKD groups vs nonAKD)") {
  
  # group係数が推定されているかチェック
  cf <- names(coef(fit))
  if (!any(grepl("^group", cf))) {
    stop("This Cox model does not contain estimated terms starting with 'group'.\nCheck that the model formula includes 'group' and that group has >=2 levels.")
  }
  
  # tidy（1回だけ）
  td <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  
  # nonAKD基準の2比較
  target_terms <- c("groupRecovery", "groupNon-Recovery")
  target_terms <- target_terms[target_terms %in% td$term]
  if (length(target_terms) == 0) {
    stop("No target terms (groupRecovery / groupNon-Recovery) found in broom::tidy(fit).\nPossible causes: only one group level remains or no events in a group.")
  }
  
  td2 <- td %>%
    filter(term %in% target_terms) %>%
    mutate(
      Contrast = dplyr::recode(
        term,
        "groupRecovery"     = "AKD Recovery vs nonAKD",
        "groupNon-Recovery" = "AKD Non-Recovery vs nonAKD",
        .default = term
      ),
      Contrast = factor(Contrast, levels = c("AKD Non-Recovery vs nonAKD", "AKD Recovery vs nonAKD"))
    )
  
  # HRテーブル
  hr_table <- td2 %>%
    mutate(
      HR        = round(estimate, 2),
      CI_lower  = round(conf.low, 2),
      CI_upper  = round(conf.high, 2),
      p_raw     = p.value,
      `p-value` = ifelse(p_raw < 0.01, "<0.01", sprintf("%.2f", p_raw))
    ) %>%
    dplyr::select(
      Contrast,
      HR,
      `Lower 95% CI` = CI_lower,
      `Upper 95% CI` = CI_upper,
      `p-value`
    )
  
  # フォレスト
  p_forest <- ggplot(td2, aes(x = Contrast, y = estimate)) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, linewidth = 0.7) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    coord_flip() +
    scale_y_log10(
      name   = "Hazard ratio for all-cause mortality (log scale)",
      breaks = hr_breaks,
      limits = hr_limits
    ) +
    labs(x = NULL, title = title_txt) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.y      = element_text(size = 10),
      axis.text.x      = element_text(size = 9),
      plot.title       = element_text(hjust = 0.5, face = "bold"),
      plot.margin      = margin(t = 5.5, r = 10, b = 5.5, l = 7, unit = "pt")
    )
  
  # tableGrob + 結合
  table_theme <- gridExtra::ttheme_minimal(
    core = list(fontsize = 9, padding = unit(c(2, 2), "mm")),
    colhead = list(fontsize = 9, fontface = "bold", padding = unit(c(2, 2), "mm"))
  )
  table_grob <- gridExtra::tableGrob(hr_table, rows = NULL, theme = table_theme)
  
  combined <- gridExtra::arrangeGrob(
    table_grob,
    p_forest,
    ncol = 1,
    heights = c(1, 2)
  )
  
  # 保存：TIFFは ragg（安定）
  if (!requireNamespace("ragg", quietly = TRUE)) install.packages("ragg")
  ragg::agg_tiff(
    filename = out_tif,
    width = width_mm, height = height_mm, units = "mm",
    res = res_dpi,
    compression = "lzw"
  )
  grid::grid.newpage()
  grid::grid.draw(combined)
  dev.off()
  
  # PDF（確認用）
  pdf(out_pdf, width = 7, height = 6)
  grid::grid.draw(combined)
  dev.off()
  
  invisible(list(hr_table = hr_table, forest = p_forest, combined = combined))
}

# ---- 3) 実行：Figure4B を出力（fitの衝突なし）----
make_forest_with_table_2contrast_group(
  fit = fit_sens,
  out_tif = "X:/R/Figure4B_forest_with_table_death_akd_3groups_sensitivity.tif",
  out_pdf = "X:/R/Figure4B_forest_with_table_death_akd_3groups_sensitivity_check.pdf",
  title_txt = "Sensitivity analysis: Adjusted hazard ratios for all-cause mortality\n(AKD groups vs nonAKD)"
)
}
#感度分析まとめ
  # ============================================================
  # Figure 4 (Sensitivity: All-Cause Mortality) - ONE PASTE
  #  (A) KM + risk table
  #  (B) Cox HR table + forest
  #  - Labels: nonAKD / AKD with Recovery / AKD without Recovery
  #  - Figure title: "Figure 4. Sensitivity Analyses: All-Cause Mortality"
  #  - Keep (A)/(B); remove "(sensitivity analysis)" from subtitles
  #
  # Assumes these objects already exist in your session:
  #   - jin1_Eligibile_sens   (sensitivity dataset; must include: group, time_years, primary_death, age, index_cre, arb_acei_use, dn*, CKD_status)
  #   - p_death_sens          (ggsurvplot object for KM in sensitivity analysis)
  # ============================================================
  
  library(dplyr)
  library(survival)
  library(broom)
  library(ggplot2)
  library(gridExtra)
  library(grid)

library(grid)
library(gridExtra)
library(ggplot2)
  
#感度分析を行うのに必要なcodeのみ####
 #感度分析死亡のカプランマイヤー####
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
    filename   = "Figure4A_death_KM_sensitivity.tiff",
    width      = 7.1,   # inch（約180 mm）
    height     = 6.0,   # inch（適宜調整）
    units      = "in",
    res        = 600,
    compression = "lzw",
    type       = "cairo"   # ← これが CairoTIFF の代わり
  )
  
  print(p_death_sens)      # カーブ + リスクテーブルをまとめて出力
  
  dev.off()

  
  # ---- 0) Safety: CKD_status ref fixed to nonCKD ----
  jin1_Eligibile_sens <- jin1_Eligibile_sens %>%
    mutate(
      CKD_status = case_when(
        as.character(CKD_status) %in% c("CKD", "1")     ~ "CKD",
        as.character(CKD_status) %in% c("nonCKD", "0") ~ "nonCKD",
        TRUE ~ as.character(CKD_status)
      ),
      CKD_status = factor(CKD_status, levels = c("nonCKD", "CKD"))
    ) %>%
    filter(!is.na(CKD_status))
  
  # ---- 1) Cox model（Sensitivity: ref=nonAKD）----
  dat_sens <- jin1_Eligibile_sens %>%
    mutate(group = relevel(as.factor(group), ref = "nonAKD")) %>%
    droplevels()
  
  fit_sens <- coxph(
    Surv(time_years, primary_death) ~
      group + age + index_cre + arb_acei_use +
      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15 +
      CKD_status,
    data = dat_sens
  )
  
  # ---- 2) Define function (IMPORTANT: define BEFORE calling) ----
  make_forest_with_table_2contrast_group <- function(fit,
                                                     out_tif,
                                                     out_pdf,
                                                     width_mm = 160,
                                                     height_mm = 140,
                                                     res_dpi = 600,
                                                     hr_limits = c(0.5, 7),
                                                     hr_breaks = c(0.5, 1, 2, 3, 4, 7),
                                                     title_txt = NULL) {
    
    td <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
    
    # expected terms from group factor (default: groupRecovery / groupNon-Recovery)
    target_terms <- c("groupRecovery", "groupNon-Recovery")
    target_terms <- target_terms[target_terms %in% td$term]
    if (length(target_terms) == 0) {
      stop("No target terms found. Available terms include:\n",
           paste(td$term, collapse = ", "),
           "\nCheck group levels and hyphen/spaces in level names.")
    }
    
    td2 <- td %>%
      filter(term %in% target_terms) %>%
      mutate(
        Contrast = dplyr::recode(
          term,
          "groupRecovery"     = "AKD with Recovery vs nonAKD",
          "groupNon-Recovery" = "AKD without Recovery vs nonAKD",
          .default = term
        ),
        Contrast = factor(
          Contrast,
          levels = c("AKD without Recovery vs nonAKD", "AKD with Recovery vs nonAKD")
        )
      )
    
    hr_table <- td2 %>%
      mutate(
        HR        = round(estimate, 2),
        CI_lower  = round(conf.low, 2),
        CI_upper  = round(conf.high, 2),
        p_raw     = p.value,
        `p-value` = ifelse(p_raw < 0.01, "<0.01", sprintf("%.2f", p_raw))
      ) %>%
      dplyr::select(
        Contrast,
        HR,
        `Lower 95% CI` = CI_lower,
        `Upper 95% CI` = CI_upper,
        `p-value`
      )
    
    p_forest <- ggplot(td2, aes(x = Contrast, y = estimate)) +
      geom_point(size = 2.5) +
      geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, linewidth = 0.7) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      coord_flip() +
      scale_y_log10(
        name   = "Hazard ratio for all-cause mortality (log scale)",
        breaks = hr_breaks,
        limits = hr_limits
      ) +
      labs(x = NULL, title = title_txt) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        axis.text.y      = element_text(size = 10),
        axis.text.x      = element_text(size = 9),
        plot.title       = element_text(hjust = 0.5, face = "bold"),
        plot.margin      = margin(t = 5.5, r = 10, b = 5.5, l = 7, unit = "pt")
      )
    
    table_theme <- gridExtra::ttheme_minimal(
      core = list(fontsize = 9, padding = unit(c(2, 2), "mm")),
      colhead = list(fontsize = 9, fontface = "bold", padding = unit(c(2, 2), "mm"))
    )
    table_grob <- gridExtra::tableGrob(hr_table, rows = NULL, theme = table_theme)
    
    combined <- gridExtra::arrangeGrob(
      table_grob,
      p_forest,
      ncol = 1,
      heights = c(1, 2)
    )
    
    # Save standalone panel (optional but keeps your original behavior)
    if (!requireNamespace("ragg", quietly = TRUE)) install.packages("ragg")
    ragg::agg_tiff(filename = out_tif, width = width_mm, height = height_mm, units = "mm",
                   res = res_dpi, compression = "lzw")
    grid::grid.newpage(); grid::grid.draw(combined); dev.off()
    
    pdf(out_pdf, width = 7, height = 6)
    grid::grid.draw(combined)
    dev.off()
    
    invisible(list(hr_table = hr_table, forest = p_forest, combined = combined))
  }
  
  # ---- 3) Create Figure4B grob (no forest title) ----
  res4B <- make_forest_with_table_2contrast_group(
    fit = fit_sens,
    out_tif = "X:/R/Figure4B_forest_with_table_death_akd_3groups_sensitivity.tif",
    out_pdf = "X:/R/Figure4B_forest_with_table_death_akd_3groups_sensitivity_check.pdf",
    title_txt = NULL
  )
  cox4B_grob <- res4B$combined
  
  # ---- 4) Figure4A: remove internal titles + fix legend labels ----
  p_death_sens2 <- p_death_sens
  p_death_sens2$plot  <- p_death_sens2$plot  + theme(plot.title = element_blank())
  if (!is.null(p_death_sens2$table)) p_death_sens2$table <- p_death_sens2$table + theme(plot.title = element_blank())
  
  # Update legend labels (works for survminer ggplot objects)
  p_death_sens2$plot <- p_death_sens2$plot +
    scale_color_discrete(labels = c("nonAKD", "AKD with Recovery", "AKD without Recovery")) +
    scale_fill_discrete(labels  = c("nonAKD", "AKD with Recovery", "AKD without Recovery"))
  
  # ---- 5) Panel labels ----
  km4A_with_label <- gridExtra::arrangeGrob(
    grid::textGrob(
      "(A) Cumulative incidence of all-cause death",
      x = unit(0, "npc"), y = unit(1, "npc"),
      just = c("left", "top"),
      gp = grid::gpar(fontsize = 12, fontface = "bold")
    ),
    gridExtra::arrangeGrob(
      p_death_sens2$plot,
      p_death_sens2$table,
      ncol = 1,
      heights = c(3, 1)
    ),
    ncol = 1,
    heights = c(0.08, 1)
  )
  
  cox4B_with_label <- gridExtra::arrangeGrob(
    grid::textGrob(
      "(B) Adjusted hazard ratios for all-cause mortality",
      x = unit(0, "npc"), y = unit(1, "npc"),
      just = c("left", "top"),
      gp = grid::gpar(fontsize = 12, fontface = "bold")
    ),
    cox4B_grob,
    ncol = 1,
    heights = c(0.08, 1)
  )
  
  # ---- 6) Figure title + combine A/B ----
  fig4_AB_labeled <- gridExtra::arrangeGrob(
    grid::textGrob(
      "Figure 4. Sensitivity Analyses: All-Cause Mortality",
      x = unit(0, "npc"), y = unit(1, "npc"),
      just = c("left", "top"),
      gp = grid::gpar(fontsize = 14, fontface = "bold")
    ),
    gridExtra::arrangeGrob(
      km4A_with_label,
      cox4B_with_label,
      ncol = 1,
      heights = c(1.2, 1.0)
    ),
    ncol = 1,
    heights = c(0.06, 1)
  )
  
  # ---- 7) Save combined figure ----
  out_tif <- "X:/R/Figure4_combined_sensitivity_A_KM_B_Cox_labeled.tif"
  out_pdf <- "X:/R/Figure4_combined_sensitivity_A_KM_B_Cox_labeled_check.pdf"
  
  if (!requireNamespace("ragg", quietly = TRUE)) install.packages("ragg")
  ragg::agg_tiff(filename = out_tif, width = 180, height = 300, units = "mm", res = 600, compression = "lzw")
  grid::grid.newpage(); grid::grid.draw(fig4_AB_labeled); dev.off()
  
  pdf(out_pdf, width = 7.1, height = 11.5)
  grid::grid.draw(fig4_AB_labeled)
  dev.off()
  
  # optional preview
  grid::grid.newpage()
  grid::grid.draw(fig4_AB_labeled)
  
  
  # ============================================================
  # Figure 4 (Sensitivity: All-Cause Mortality) - ONE PASTE
  #  (A) KM + risk table
  #  (B) Cox HR table + forest
  #  - Labels: nonAKD / AKD with Recovery / AKD without Recovery
  #  - Figure title: "Figure 4. Sensitivity Analyses: All-Cause Mortality"
  #  - Keep (A)/(B); remove "(sensitivity analysis)" from subtitles
  #
  # Assumes these objects already exist in your session:
  #   - jin1_Eligibile_sens   (sensitivity dataset; must include: group, time_years, primary_death, age, index_cre, arb_acei_use, dn*, CKD_status)
  #   - p_death_sens          (ggsurvplot object for KM in sensitivity analysis)
  # ============================================================
  
  library(dplyr)
  library(survival)
  library(broom)
  library(ggplot2)
  library(gridExtra)
  library(grid)
#感度分析死亡のカプランマイヤー####
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
# ---- 0) Safety: CKD_status ref を nonCKD に固定（本解析と揃える）----
jin1_Eligibile_sens <- jin1_Eligibile_sens %>%
  mutate(
    CKD_status = case_when(
      as.character(CKD_status) %in% c("CKD", "1")     ~ "CKD",
      as.character(CKD_status) %in% c("nonCKD", "0") ~ "nonCKD",
      TRUE ~ as.character(CKD_status)
    ),
    CKD_status = factor(CKD_status, levels = c("nonCKD", "CKD"))
  ) %>%
  filter(!is.na(CKD_status))

# ---- 1) Cox model（感度分析：group, ref=nonAKD）----
dat_sens <- jin1_Eligibile_sens %>%
  mutate(group = relevel(group, ref = "nonAKD")) %>%
  droplevels()

fit_sens <- coxph(
  Surv(time_years, primary_death) ~
    group + age + index_cre + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15 +
    CKD_status,
  data = dat_sens
)

# ---- 2) Forest + HR table function (2 contrasts vs nonAKD) ----
make_forest_with_table_2contrast_group <- function(fit,
                                                   out_tif,
                                                   out_pdf,
                                                   width_mm = 160,
                                                   height_mm = 140,
                                                   res_dpi = 600,
                                                   hr_limits = c(0.5, 7),
                                                   hr_breaks = c(0.5, 1, 2, 3, 4, 7),
                                                   title_txt = NULL) {
  
  cf <- names(coef(fit))
  if (!any(grepl("^group", cf))) {
    stop("This Cox model does not contain estimated terms starting with 'group'.")
  }
  
  td <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  
  target_terms <- c("groupRecovery", "groupNon-Recovery")
  target_terms <- target_terms[target_terms %in% td$term]
  if (length(target_terms) == 0) {
    stop("No target terms (groupRecovery / groupNon-Recovery) found in broom::tidy(fit).")
  }
  
  td2 <- td %>%
    filter(term %in% target_terms) %>%
    mutate(
      Contrast = dplyr::recode(
        term,
        "groupRecovery"     = "AKD with Recovery vs nonAKD",
        "groupNon-Recovery" = "AKD without Recovery vs nonAKD",
        .default = term
      ),
      Contrast = factor(
        Contrast,
        levels = c("AKD without Recovery vs nonAKD", "AKD with Recovery vs nonAKD")
      )
    )
  
  hr_table <- td2 %>%
    mutate(
      HR        = round(estimate, 2),
      CI_lower  = round(conf.low, 2),
      CI_upper  = round(conf.high, 2),
      p_raw     = p.value,
      `p-value` = ifelse(p_raw < 0.01, "<0.01", sprintf("%.2f", p_raw))
    ) %>%
    dplyr::select(
      Contrast,
      HR,
      `Lower 95% CI` = CI_lower,
      `Upper 95% CI` = CI_upper,
      `p-value`
    )
  
  p_forest <- ggplot(td2, aes(x = Contrast, y = estimate)) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, linewidth = 0.7) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    coord_flip() +
    scale_y_log10(
      name   = "Hazard ratio for all-cause mortality (log scale)",
      breaks = hr_breaks,
      limits = hr_limits
    ) +
    labs(x = NULL, title = title_txt) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.y      = element_text(size = 10),
      axis.text.x      = element_text(size = 9),
      plot.title       = element_text(hjust = 0.5, face = "bold"),
      plot.margin      = margin(t = 5.5, r = 10, b = 5.5, l = 7, unit = "pt")
    )
  
  table_theme <- gridExtra::ttheme_minimal(
    core = list(fontsize = 9, padding = unit(c(2, 2), "mm")),
    colhead = list(fontsize = 9, fontface = "bold", padding = unit(c(2, 2), "mm"))
  )
  table_grob <- gridExtra::tableGrob(hr_table, rows = NULL, theme = table_theme)
  
  combined <- gridExtra::arrangeGrob(
    table_grob,
    p_forest,
    ncol = 1,
    heights = c(1, 2)
  )
  
  if (!requireNamespace("ragg", quietly = TRUE)) install.packages("ragg")
  ragg::agg_tiff(
    filename = out_tif,
    width = width_mm, height = height_mm, units = "mm",
    res = res_dpi,
    compression = "lzw"
  )
  grid::grid.newpage()
  grid::grid.draw(combined)
  dev.off()
  
  pdf(out_pdf, width = 7, height = 6)
  grid::grid.draw(combined)
  dev.off()
  
  invisible(list(hr_table = hr_table, forest = p_forest, combined = combined))
}

# ---- 3) Build Figure4B grob (no forest title) ----
res4B <- make_forest_with_table_2contrast_group(
  fit = fit_sens,
  out_tif = "X:/R/Figure4B_forest_with_table_death_akd_3groups_sensitivity.tif",
  out_pdf = "X:/R/Figure4B_forest_with_table_death_akd_3groups_sensitivity_check.pdf",
  title_txt = NULL
)
cox4B_grob <- res4B$combined

# ---- 4) Figure4A: remove internal titles + fix legend labels ----
p_death_sens2 <- p_death_sens
p_death_sens2$plot  <- p_death_sens2$plot  + theme(plot.title = element_blank())
if (!is.null(p_death_sens2$table)) p_death_sens2$table <- p_death_sens2$table + theme(plot.title = element_blank())

# ★ Legend labels (display only)
p_death_sens2$plot <- p_death_sens2$plot +
  scale_color_discrete(labels = c("nonAKD", "AKD with Recovery", "AKD without Recovery")) +
  scale_fill_discrete(labels  = c("nonAKD", "AKD with Recovery", "AKD without Recovery"))

# ---- 5) Panel labels (keep A/B; remove '(sensitivity analysis)') ----
km4A_with_label <- gridExtra::arrangeGrob(
  grid::textGrob(
    "(A) Cumulative incidence of all-cause death",
    x = unit(0, "npc"), y = unit(1, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 12, fontface = "bold")
  ),
  gridExtra::arrangeGrob(
    p_death_sens2$plot,
    p_death_sens2$table,
    ncol = 1,
    heights = c(3, 1)
  ),
  ncol = 1,
  heights = c(0.08, 1)
)

cox4B_with_label <- gridExtra::arrangeGrob(
  grid::textGrob(
    "(B) Adjusted hazard ratios for all-cause mortality",
    x = unit(0, "npc"), y = unit(1, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 12, fontface = "bold")
  ),
  cox4B_grob,
  ncol = 1,
  heights = c(0.08, 1)
)

# ---- 6) Figure title + combine A/B ----
fig4_AB_labeled <- gridExtra::arrangeGrob(
  grid::textGrob(
    "Figure 4. Sensitivity Analyses: All-Cause Mortality",
    x = unit(0, "npc"), y = unit(1, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 14, fontface = "bold")
  ),
  gridExtra::arrangeGrob(
    km4A_with_label,
    cox4B_with_label,
    ncol = 1,
    heights = c(1.2, 1.0)
  ),
  ncol = 1,
  heights = c(0.06, 1)
)

# ---- 7) Save combined figure ----
out_tif <- "X:/R/Figure4_combined_sensitivity_A_KM_B_Cox_labeled.tif"
out_pdf <- "X:/R/Figure4_combined_sensitivity_A_KM_B_Cox_labeled_check.pdf"

if (!requireNamespace("ragg", quietly = TRUE)) install.packages("ragg")
ragg::agg_tiff(
  filename = out_tif,
  width = 180,
  height = 300,
  units = "mm",
  res = 600,
  compression = "lzw"
)
grid::grid.newpage()
grid::grid.draw(fig4_AB_labeled)
dev.off()

pdf(out_pdf, width = 7.1, height = 11.5)
grid::grid.draw(fig4_AB_labeled)
dev.off()

# optional preview
grid::grid.newpage()
grid::grid.draw(fig4_AB_labeled)





##併存疾患表の作成
{
#####
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

}
#競合エンドポイントについてのカプランマイヤー 2025/8/5
{
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
}

setwd("X:/R")
install.packages("qpdf")
library(qpdf)
pdfs <- c(
  "Figure1_flowchart_AKD.pdf",
  "Figure2_combined_vertical_labeled_check.pdf",
  "Figure3_eGFR_trajectory_and_slope_1y.pdf",
  "Figure4_combined_sensitivity_A_KM_B_Cox_labeled_check.pdf",
  "Figure5_sensitivity_trajectory_and_slope_1y.pdf",
  "SupplementalFigure1_main_trajectory_and_bar_3y_all.pdf",
  "SupplementalFigure2_sensitivity_bar_and_trajectory_3y_all.pdf",
  "SupplementalFigure3_interaction_forest_CKD_check.pdf",
  "Supplementary_Table_1.pdf",
  "Table1_Baseline_with_CKD.pdf"
)

pdf_combine(pdfs, output = "All_Figures_and_Tables.pdf")
