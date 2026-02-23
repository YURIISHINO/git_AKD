{

############################################################
# Figure 2 (KM only: Legend + KM + Number at risk)
#  + Table 2 (HR only)  --- ONE-PASS COPY-PASTE ---
#  - Figure2 is packed into the upper half; large blank space at bottom
#  - Risk panel height is compressed (reduces "spaced-out" impression)
#  - Risk row spacing controlled by fixed y_map (stable)
#  - Table2 saved separately (bordered) + CSV
############################################################

graphics.off()

# --------------------------
# Packages
# --------------------------
pkgs <- c("readr","dplyr","survival","broom","ggplot2",
          "grid","gridExtra","gtable","tibble","ragg")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

library(readr); library(dplyr); library(survival); library(broom); library(ggplot2)
library(grid); library(gridExtra); library(gtable); library(tibble)

# ==========================================================
# TUNING (★ここだけ調整すればOK)
# ==========================================================
base_fs <- 10
risk_fs <- 11
hr_fs   <- 11

km_line_lwd <- 1.05
km_ci_alpha <- 0.18

num_size <- 3.6            # risk numbers
leg_text_size <- 3.2       # legend text size

# ---- Risk row spacing (fixed y) ----
# closer values => tighter rows
y_map <- c(
  "non-AKD"              = 0.50,
  "AKD with recovery"    = 0.40,
  "AKD without recovery" = 0.30
)

# ---- Risk title spacing ----
risk_title_mb <- 2         # title -> table (downward)
risk_margin_t <- 10        # top margin of risk panel (more => more space above title)
risk_margin_b <- 2

# ---- Risk panel view window (controls title-to-table & bottom blank inside risk) ----
# smaller upper => title closer to 1st row; larger lower => less bottom blank
risk_ylim <- c(0.28, 0.56)  # ★あなたの現状に合わせた推奨

# ---- Make the whole Figure2 packed to upper half ----
# 1) KM vs Risk height ratio (smaller risk => tighter look)
km_vs_risk_heights <- c(2.00, 1.30)   # (KM, risk) 例: (3.8,0.75)でさらにrisk圧縮
# 2) Add large blank at bottom of final Figure2
fig2_heights <- c(0.32, 4.40, 1.20)  # (legend, KM+risk, bottom blank)

# ==========================================================
# Paths
# ==========================================================
setwd("X:/R")
in_csv <- "jin1_Eligibile.csv"

out_dir <- file.path("X:/R","primary")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

out_fig2_tif <- file.path(out_dir,"Figure2_primary_KM.tif")
out_fig2_pdf <- file.path(out_dir,"Figure2_primary_KM.pdf")
out_tab2_tif <- file.path(out_dir,"Table2_primary_HR.tif")
out_tab2_pdf <- file.path(out_dir,"Table2_primary_HR.pdf")
out_tab2_csv <- file.path(out_dir,"Table2_primary_HR.csv")

# ==========================================================
# Load & build 1 row per patient
# ==========================================================
jin1_Eligibile <- read_csv(in_csv, locale = locale(encoding = "SHIFT-JIS"))

dat1 <- jin1_Eligibile %>%
  filter(exclude == "include", jin_status %in% c("AKD","nonAKD")) %>%
  group_by(id) %>%
  arrange(index_date, date, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

levels_full <- c("non-AKD","AKD with recovery","AKD without recovery")

dat_km <- dat1 %>%
  mutate(
    group = case_when(
      jin_status == "nonAKD" ~ "non-AKD",
      jin_status == "AKD" & `150_210recovery` == 1 ~ "AKD with recovery",
      jin_status == "AKD" & `150_210recovery` == 2 ~ "AKD without recovery",
      jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 1 ~ "AKD with recovery",
      jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` %in% c(0,2) ~ "AKD without recovery",
      TRUE ~ NA_character_
    ),
    group = factor(group, levels = levels_full),
    time_years = as.numeric(last_follow_death - index_plus_210)/365.25
  ) %>%
  filter(!is.na(group), !is.na(time_years), time_years >= 0)

# ==========================================================
# Common axis / colors
# ==========================================================
ticks_show <- seq(0, 8, by = 2)
x_right <- 9.5
x_right_risk <- 8.2

pal <- c(
  "non-AKD"              = "#95A5A6",
  "AKD with recovery"    = "#2ECC71",
  "AKD without recovery" = "#E74C3C"
)

common_left_margin_pt  <- 80
common_right_margin_pt <- 80

# ==========================================================
# KM fit
# ==========================================================
fit <- survfit(Surv(time_years, primary_death) ~ group, data = dat_km)

# ==========================================================
# KM panel (manual CIF)
# ==========================================================
s <- summary(fit)

km_df <- data.frame(
  time   = s$time,
  surv   = s$surv,
  lower  = s$lower,
  upper  = s$upper,
  strata = s$strata
) %>%
  mutate(
    group = factor(sub("^group=","", strata), levels = levels_full),
    cif   = 1 - surv,
    cif_l = 1 - upper,
    cif_u = 1 - lower
  ) %>%
  arrange(group, time)

p_km <- ggplot(km_df, aes(x = time, y = cif, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = cif_l, ymax = cif_u),
              alpha = km_ci_alpha, linewidth = 0, show.legend = FALSE) +
  geom_step(linewidth = km_line_lwd, show.legend = FALSE) +
  scale_x_continuous(breaks = ticks_show, limits = c(0, x_right), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.40), breaks = seq(0, 0.4, 0.1), expand = c(0, 0)) +
  scale_colour_manual(values = pal, breaks = levels_full) +
  scale_fill_manual(values = pal, breaks = levels_full) +
  labs(x = NULL, y = "Cumulative incidence") +
  theme_classic(base_size = base_fs) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin  = margin(t = 4, r = common_right_margin_pt, b = 10,   # ★ここ(b)を増やす
                          l = common_left_margin_pt, unit = "pt")
  )

# ==========================================================
# Legend panel (full width above KM)
# ==========================================================
leg_df <- tibble(
  group = factor(levels_full, levels = levels_full),
  x0 = c(1.2, 10.0, 20.0),
  y  = 1.10
) %>%
  mutate(x1 = x0 + 0.85)

p_leg <- ggplot(leg_df) +
  geom_rect(aes(xmin = x0, xmax = x1, ymin = y - 0.18, ymax = y + 0.18, fill = group),
            alpha = 0.25, colour = NA) +
  geom_segment(aes(x = x0, xend = x1, y = y, yend = y, colour = group),
               linewidth = 1.1) +
  geom_text(aes(x = x1 + 0.55, y = y, label = group),
            hjust = 0, size = leg_text_size) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_colour_manual(values = pal, guide = "none") +
  coord_cartesian(xlim = c(0.5, 30.0), ylim = c(0.6, 1.7), clip = "off") +
  theme_void(base_size = base_fs) +
  theme(
    plot.margin = margin(t = 0, r = common_right_margin_pt, b = 0,
                         l = common_left_margin_pt, unit = "pt")
  )

# ==========================================================
# Number at risk (0-align version)
#  - x-axis starts at 0 (align with KM)
#  - group labels drawn OUTSIDE panel (annotation_custom)
#  - x-axis line moved downward by lowering risk_ylim[1]
# ==========================================================
sfit <- summary(fit, times = ticks_show, extend = TRUE)

# ★ 0年の数字は現状維持（動かしたくないならこのまま）
time0_shift <- 0.22

# ★ riskパネルの下を深くするほど、Yearsの横棒（x軸線）が下に降り、ラベルと重ならない
risk_ylim <- c(0.18, 0.56)   # ← 下限を 0.28→0.18 に（まずこれ）

# ★ Group名を外に出す量（mm）
label_pad_mm <- 10           # 8〜14で調整（大きいほど左へ）

risk_df <- data.frame(
  time   = sfit$time,
  strata = sfit$strata,
  n_risk = sfit$n.risk
) %>%
  mutate(
    group = factor(sub("^group=","", strata), levels = levels_full),
    y = unname(y_map[as.character(group)]),
    time_plot = ifelse(time == 0, time0_shift, time),
    group_disp = case_when(
      as.character(group) == "AKD with recovery"    ~ "AKD with\nrecovery",
      as.character(group) == "AKD without recovery" ~ "AKD without\nrecovery",
      TRUE ~ as.character(group)
    )
  )

p_risk <- ggplot() +
  # ---- numbers ----
geom_text(
  data = risk_df,
  aes(x = time_plot, y = y, label = n_risk),
  size = num_size
) +
  scale_x_continuous(
    breaks = ticks_show,
    limits = c(0, x_right_risk),     # ★0開始（KMと揃える本体）
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = sort(unique(risk_df$y)),
    labels = rep("", length(unique(risk_df$y))),
    expand = c(0, 0)
  ) +
  labs(x = "Years", title = "Number at risk") +
  theme_classic(base_size = risk_fs) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    plot.margin  = margin(t = risk_margin_t, r = common_right_margin_pt,
                          b = risk_margin_b, l = common_left_margin_pt, unit = "pt"),
    plot.title.position = "plot",
    plot.title = element_text(margin = margin(b = risk_title_mb), vjust = 0)
  ) +
  coord_cartesian(ylim = risk_ylim, clip = "off")

# ---- add GROUP LABELS outside panel (x=0より左に固定表示) ----
for (grp in levels_full) {
  yv  <- unique(risk_df$y[risk_df$group == grp])[1]
  lab <- unique(risk_df$group_disp[risk_df$group == grp])[1]
  col <- pal[grp]
  
  p_risk <- p_risk +
    annotation_custom(
      grob = textGrob(
        lab,
        x = unit(0, "npc") - unit(label_pad_mm, "mm"),
        just = "right",
        gp = gpar(col = col, fontsize = risk_fs)
      ),
      xmin = -Inf, xmax = -Inf, ymin = yv, ymax = yv
    )
}

# ==========================================================
# Bind KM + risk with controlled vertical ratio (compress risk)
# ==========================================================
g_km   <- ggplotGrob(p_km)
g_risk <- ggplotGrob(p_risk)

# hard align widths
g_risk$widths <- g_km$widths

km_risk_block <- arrangeGrob(
  g_km, g_risk,
  ncol = 1,
  heights = km_vs_risk_heights
)

# ==========================================================
# Figure 2 object (KM only) + bottom blank (pack to upper half)
# ==========================================================
fig2_onlyKM <- arrangeGrob(
  textGrob(
    "Figure 2. Primary Outcome: All-Cause Mortality",
    x = unit(0.02, "npc"),   # ← 左余白（0.02〜0.05で調整）
    y = unit(0.95, "npc"),    # ← 上余白（0.9〜0.95で調整）
    just = c("left", "top"),
    gp = gpar(fontsize = 14, fontface = "bold")
  ),
  p_leg,
  km_risk_block,
  nullGrob(),
  ncol = 1,
  heights = c(0.10, fig2_heights)  # ← タイトル分を先頭に追加
)

# ==========================================================
# Cox → HR table (Table 2 only)
# ==========================================================
dat_cox <- dat_km %>%
  mutate(arb_acei_use = if_else(coalesce(arb,0)==1 | coalesce(acei,0)==1, 1L, 0L))

fit_main <- coxph(
  Surv(time_years, primary_death) ~
    group + age + index_cre + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
  data = dat_cox
)

hr_table <- broom::tidy(fit_main, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term %in% c("groupAKD with recovery", "groupAKD without recovery")) %>%
  mutate(
    Contrast = c("AKD with recovery vs non-AKD",
                 "AKD without recovery vs non-AKD"),
    HR = sprintf("%.2f", estimate),
    `95% CI` = sprintf("%.2f–%.2f", conf.low, conf.high),
    `p-value` = ifelse(p.value < 0.01, "<0.01", sprintf("%.2f", p.value))
  ) %>%
  select(Contrast, HR, `95% CI`, `p-value`)

write_csv(hr_table, out_tab2_csv)

hr_grob <- tableGrob(
  hr_table, rows = NULL,
  theme = ttheme_default(
    base_size = hr_fs,
    core = list(bg_params = list(col = "black", lwd = 0.4)),
    colhead = list(bg_params = list(col = "black", lwd = 0.6),
                   fg_params = list(fontface = "bold"))
  )
)

# outer border
hr_grob <- gtable_add_grob(
  hr_grob,
  rectGrob(gp = gpar(fill = NA, col = "black", lwd = 1)),
  t = 1, l = 1, b = nrow(hr_grob), r = ncol(hr_grob)
)

tab2_onlyHR <- arrangeGrob(
  textGrob("Table 2. Adjusted hazard ratios for all-cause mortality",
           x = unit(0, "npc"), just = "left",
           gp = gpar(fontsize = 12, fontface = "bold")),
  hr_grob, ncol = 1, heights = c(0.18, 1)
)

# ==========================================================
# Draw (optional)
# ==========================================================
grid.newpage(); grid.draw(fig2_onlyKM)
grid.newpage(); grid.draw(tab2_onlyHR)

# ==========================================================
# SAVE
# ==========================================================
ragg::agg_tiff(out_fig2_tif, width = 180, height = 220, units = "mm", res = 600, compression = "lzw")
grid.newpage(); grid.draw(fig2_onlyKM); dev.off()

pdf(out_fig2_pdf, width = 7.8, height = 8.4)
grid.newpage(); grid.draw(fig2_onlyKM); dev.off()

ragg::agg_tiff(out_tab2_tif, width = 180, height = 120, units = "mm", res = 600, compression = "lzw")
grid.newpage(); grid.draw(tab2_onlyHR); dev.off()

pdf(out_tab2_pdf, width = 7.8, height = 4.6)
grid.newpage(); grid.draw(tab2_onlyHR); dev.off()

############################################################
# Quick knobs:
#  - pack figure up: fig2_heights = c(legend, km+risk, blank) -> increase blank
#  - shrink risk panel: km_vs_risk_heights -> make 2nd smaller
#  - risk row spacing: y_map values closer together
#  - title-to-table: risk_ylim upper (smaller -> closer)
############################################################

# ==========================================================
# Table 2 を Word（docx）でも保存（officer + flextable）
#  - タイトル行 + 罫線付きテーブル
#  - 1ページに収まりやすいように余白/フォント/幅を調整
# ==========================================================
pkgs2 <- c("officer","flextable")
to_install2 <- pkgs2[!vapply(pkgs2, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install2) > 0) install.packages(to_install2, dependencies = TRUE)

library(officer)
library(flextable)

out_tab2_docx <- file.path(out_dir, "Table2_primary_HR.docx")

# ---- flextable化（列幅/罫線/文字などを整える）----
ft <- flextable(hr_table)

ft <- ft %>%
  theme_booktabs() %>%               # すっきりした罫線（好みで theme_vanilla() でもOK）
  bold(part = "header") %>%
  align(align = "left", part = "all") %>%
  align(j = c("HR","95% CI","p-value"), align = "center", part = "all") %>%
  autofit()

# 列幅を固定（A4縦で読みやすい目安。必要なら微調整）
ft <- width(ft, j = "Contrast", width = 3.6)
ft <- width(ft, j = "HR",       width = 1.0)
ft <- width(ft, j = "95% CI",   width = 1.6)
ft <- width(ft, j = "p-value",  width = 1.0)

# フォントサイズ
ft <- fontsize(ft, size = 11, part = "all")

# 外枠を強めに（TableGrobの“外枠”相当）
outer <- fp_border(color = "black", width = 1)
inner <- fp_border(color = "black", width = 0.5)
ft <- border_remove(ft)
ft <- border_outer(ft, border = outer)
ft <- border_inner_h(ft, border = inner)
ft <- border_inner_v(ft, border = inner)

# ---- Word作成：ページ余白を少し詰めて縦1枚に収めやすく ----
doc <- read_docx()

# セクション（A4縦・余白調整）
sec <- prop_section(
  page_size = page_size(width = 8.27, height = 11.69),     # A4 (inch)
  page_margins = page_mar(
    top = 0.6, bottom = 0.6, left = 0.7, right = 0.7       # inch
  )
)

doc <- doc %>%
  body_add_par("Table 2. Adjusted hazard ratios for all-cause mortality",
               style = "heading 2") %>%
  body_add_flextable(ft) %>%
  body_add_par("", style = "Normal") %>%
  body_end_section_continuous() %>%
  body_set_default_section(sec)

print(doc, target = out_tab2_docx)

} #本解析
{
############################################################
# Figure 4 (Sensitivity)  --- ONE-PASS COPY-PASTE ---
# Layout is IDENTICAL to Primary Figure 2 code:
#  - Figure4: Legend + KM(CIF+CI) + Number at risk (packed to upper half)
#  - HR table saved separately (bordered) + CSV  (like Table 2 in primary)
# Output folder: X:/R/sensitivity_analysis/
############################################################

graphics.off()

# --------------------------
# Packages
# --------------------------
pkgs <- c("readr","dplyr","survival","broom","ggplot2",
          "grid","gridExtra","gtable","tibble","ragg")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

library(readr); library(dplyr); library(survival); library(broom); library(ggplot2)
library(grid); library(gridExtra); library(gtable); library(tibble)

# ==========================================================
# TUNING (★primaryと同じ思想：ここだけ調整すればOK)
# ==========================================================
base_fs <- 10
risk_fs <- 11
hr_fs   <- 11

km_line_lwd <- 1.05
km_ci_alpha <- 0.18

num_size <- 3.6
leg_text_size <- 3.2

# ---- Risk row spacing (fixed y: primaryと同じ) ----
y_map <- c(
  "nonAKD"        = 0.50,
  "Recovery"      = 0.40,
  "Non-Recovery"  = 0.30
)

# ---- Risk title spacing ----
risk_title_mb <- 2
risk_margin_t <- 10
risk_margin_b <- 2

# ---- Risk panel view window ----
risk_ylim <- c(0.18, 0.56)

# ---- KM vs Risk height ratio ----
km_vs_risk_heights <- c(2.00, 1.30)

# ---- Pack whole Figure4 to upper half + big blank bottom ----
fig4_heights <- c(0.32, 4.40, 1.20)  # (legend, KM+risk, bottom blank)

# ==========================================================
# Paths (★sensitivity_analysis folder)
# ==========================================================
setwd("X:/R")
in_csv <- "jin1_Eligibile.csv"

out_dir <- file.path("X:/R","sensitivity_analysis")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

out_fig4_tif <- file.path(out_dir,"Figure4_sensitivity_KM.tif")
out_fig4_pdf <- file.path(out_dir,"Figure4_sensitivity_KM.pdf")

out_tab4_tif <- file.path(out_dir,"Table_3_sensitivity_HR.tif")
out_tab4_pdf <- file.path(out_dir,"Table_3_sensitivity_HR.pdf")
out_tab4_csv <- file.path(out_dir,"Table_3_sensitivity_HR.csv")

# ==========================================================
# Load & build 1 row per patient (primaryと同じ)
# ==========================================================
jin1_Eligibile <- read_csv(in_csv, locale = locale(encoding = "SHIFT-JIS"))

dat1 <- jin1_Eligibile %>%
  filter(exclude == "include", jin_status %in% c("AKD","nonAKD")) %>%
  group_by(id) %>%
  arrange(index_date, date, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

levels_full <- c("nonAKD", "Recovery", "Non-Recovery")

dat_sens <- dat1 %>%
  mutate(
    group = case_when(
      jin_status == "nonAKD" ~ "nonAKD",
      jin_status == "AKD" & `150_210recovery` == 1 ~ "Recovery",
      jin_status == "AKD" & `150_210recovery` == 2 ~ "Non-Recovery",
      jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 1 ~ "Recovery",
      jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` %in% c(0,2) ~ "Non-Recovery",
      TRUE ~ NA_character_
    ),
    group = factor(group, levels = levels_full),
    arb_acei_use = if_else(coalesce(arb,0)==1 | coalesce(acei,0)==1, 1L, 0L),
    time_years = as.numeric(last_follow_death - index_plus_210)/365.25
  ) %>%
  filter(!is.na(group), !is.na(time_years), time_years >= 0)

# ==========================================================
# Common axis / colors (primaryと同じ)
# ==========================================================
ticks_show <- seq(0, 8, by = 2)
x_right <- 9.5
x_right_risk <- 8.2

pal <- c(
  "nonAKD"       = "#95A5A6",
  "Recovery"     = "#2ECC71",
  "Non-Recovery" = "#E74C3C"
)

common_left_margin_pt  <- 80
common_right_margin_pt <- 80

# ==========================================================
# KM fit
# ==========================================================
fit <- survfit(Surv(time_years, primary_death) ~ group, data = dat_sens)

# ==========================================================
# KM panel (manual CIF)  ※primaryと同じ
# ==========================================================
s <- summary(fit)

km_df <- data.frame(
  time   = s$time,
  surv   = s$surv,
  lower  = s$lower,
  upper  = s$upper,
  strata = s$strata
) %>%
  mutate(
    group = factor(sub("^group=","", strata), levels = levels_full),
    cif   = 1 - surv,
    cif_l = 1 - upper,
    cif_u = 1 - lower
  ) %>%
  arrange(group, time)

p_km <- ggplot(km_df, aes(x = time, y = cif, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = cif_l, ymax = cif_u),
              alpha = km_ci_alpha, linewidth = 0, show.legend = FALSE) +
  geom_step(linewidth = km_line_lwd, show.legend = FALSE) +
  scale_x_continuous(breaks = ticks_show, limits = c(0, x_right), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.40), breaks = seq(0, 0.4, 0.1), expand = c(0, 0)) +
  scale_colour_manual(values = pal, breaks = levels_full) +
  scale_fill_manual(values = pal, breaks = levels_full) +
  labs(x = NULL, y = "Cumulative incidence") +
  theme_classic(base_size = base_fs) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin  = margin(t = 4, r = common_right_margin_pt, b = 10,
                          l = common_left_margin_pt, unit = "pt")
  )

# ==========================================================
# Legend panel (full width above KM) ※primaryと同じ構造
#  表示名（改行なし）をここで作る
# ==========================================================
leg_levels_disp <- c(
  "nonAKD"       = "nonAKD",
  "Recovery"     = "AKD with recovery",
  "Non-Recovery" = "AKD without recovery"
)

leg_df <- tibble(
  group = factor(levels_full, levels = levels_full),
  x0 = c(1.2, 10.0, 20.0),
  y  = 1.10,
  label = unname(leg_levels_disp[levels_full])
) %>%
  mutate(x1 = x0 + 0.85)

p_leg <- ggplot(leg_df) +
  geom_rect(aes(xmin = x0, xmax = x1, ymin = y - 0.18, ymax = y + 0.18, fill = group),
            alpha = 0.25, colour = NA) +
  geom_segment(aes(x = x0, xend = x1, y = y, yend = y, colour = group),
               linewidth = 1.1) +
  geom_text(aes(x = x1 + 0.55, y = y, label = label),
            hjust = 0, size = leg_text_size) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_colour_manual(values = pal, guide = "none") +
  coord_cartesian(xlim = c(0.5, 30.0), ylim = c(0.6, 1.7), clip = "off") +
  theme_void(base_size = base_fs) +
  theme(
    plot.margin = margin(t = 0, r = common_right_margin_pt, b = 0,
                         l = common_left_margin_pt, unit = "pt")
  )

# ==========================================================
# Number at risk (0-align + fixed y_map) ※primaryと同じ
# ==========================================================
sfit <- summary(fit, times = ticks_show, extend = TRUE)

time0_shift <- 0.22
label_pad_mm <- 10

risk_df <- data.frame(
  time   = sfit$time,
  strata = sfit$strata,
  n_risk = sfit$n.risk
) %>%
  mutate(
    group = factor(sub("^group=","", strata), levels = levels_full),
    y = unname(y_map[as.character(group)]),
    time_plot = ifelse(time == 0, time0_shift, time),
    group_disp = case_when(
      as.character(group) == "Recovery"     ~ "AKD with\nrecovery",
      as.character(group) == "Non-Recovery" ~ "AKD without\nrecovery",
      TRUE ~ as.character(group)
    )
  )

p_risk <- ggplot() +
  geom_text(
    data = risk_df,
    aes(x = time_plot, y = y, label = n_risk),
    size = num_size
  ) +
  scale_x_continuous(
    breaks = ticks_show,
    limits = c(0, x_right_risk),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = sort(unique(risk_df$y)),
    labels = rep("", length(unique(risk_df$y))),
    expand = c(0, 0)
  ) +
  labs(x = "Years", title = "Number at risk") +
  theme_classic(base_size = risk_fs) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    plot.margin  = margin(t = risk_margin_t, r = common_right_margin_pt,
                          b = risk_margin_b, l = common_left_margin_pt, unit = "pt"),
    plot.title.position = "plot",
    plot.title = element_text(margin = margin(b = risk_title_mb), vjust = 0)
  ) +
  coord_cartesian(ylim = risk_ylim, clip = "off")

# ---- add GROUP LABELS outside panel ----
for (grp in levels_full) {
  yv  <- unique(risk_df$y[risk_df$group == grp])[1]
  lab <- unique(risk_df$group_disp[risk_df$group == grp])[1]
  col <- pal[grp]
  
  p_risk <- p_risk +
    annotation_custom(
      grob = textGrob(
        lab,
        x = unit(0, "npc") - unit(label_pad_mm, "mm"),
        just = "right",
        gp = gpar(col = col, fontsize = risk_fs)
      ),
      xmin = -Inf, xmax = -Inf, ymin = yv, ymax = yv
    )
}

# ==========================================================
# Bind KM + risk with controlled vertical ratio (compress risk)
# ==========================================================
g_km   <- ggplotGrob(p_km)
g_risk <- ggplotGrob(p_risk)

# hard align widths
g_risk$widths <- g_km$widths

km_risk_block <- arrangeGrob(
  g_km, g_risk,
  ncol = 1,
  heights = km_vs_risk_heights
)

# ==========================================================
# Figure 4 object (KM only) + bottom blank (pack to upper half)
# ==========================================================
fig4_onlyKM <- arrangeGrob(
  textGrob(
    "Figure 4. Sensitivity Analyses: All-Cause Mortality",
    x = unit(0.02, "npc"),
    y = unit(0.95, "npc"),
    just = c("left", "top"),
    gp = gpar(fontsize = 14, fontface = "bold")
  ),
  p_leg,
  km_risk_block,
  nullGrob(),
  ncol = 1,
  heights = c(0.10, fig4_heights)
)

# ==========================================================
# Cox → HR table (separate file, like primary Table2)
# ==========================================================
fit_cox <- coxph(
  Surv(time_years, primary_death) ~
    group + age + index_cre + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
  data = dat_sens
)

hr_table <- broom::tidy(fit_cox, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term %in% c("groupRecovery", "groupNon-Recovery")) %>%
  mutate(
    Contrast = c("AKD with recovery vs nonAKD",
                 "AKD without recovery vs nonAKD"),
    HR = sprintf("%.2f", estimate),
    `95% CI` = sprintf("%.2f–%.2f", conf.low, conf.high),
    `p-value` = ifelse(p.value < 0.01, "<0.01", sprintf("%.2f", p.value))
  ) %>%
  select(Contrast, HR, `95% CI`, `p-value`)

write_csv(hr_table, out_tab4_csv)

hr_grob <- tableGrob(
  hr_table, rows = NULL,
  theme = ttheme_default(
    base_size = hr_fs,
    core = list(bg_params = list(col = "black", lwd = 0.4)),
    colhead = list(bg_params = list(col = "black", lwd = 0.6),
                   fg_params = list(fontface = "bold"))
  )
)

# outer border
hr_grob <- gtable_add_grob(
  hr_grob,
  rectGrob(gp = gpar(fill = NA, col = "black", lwd = 1)),
  t = 1, l = 1, b = nrow(hr_grob), r = ncol(hr_grob)
)

tab4_onlyHR <- arrangeGrob(
  textGrob("Table 3. Adjusted hazard ratios for all-cause mortality",
           x = unit(0, "npc"), just = "left",
           gp = gpar(fontsize = 12, fontface = "bold")),
  hr_grob, ncol = 1, heights = c(0.18, 1)
)

# ==========================================================
# Draw (optional check)
# ==========================================================
grid.newpage(); grid.draw(fig4_onlyKM)
grid.newpage(); grid.draw(tab4_onlyHR)

# ==========================================================
# SAVE (★sensitivity_analysisへ)
# ==========================================================
ragg::agg_tiff(out_fig4_tif, width = 180, height = 220, units = "mm", res = 600, compression = "lzw")
grid.newpage(); grid.draw(fig4_onlyKM); dev.off()

pdf(out_fig4_pdf, width = 7.8, height = 8.4)
grid.newpage(); grid.draw(fig4_onlyKM); dev.off()

ragg::agg_tiff(out_tab4_tif, width = 180, height = 120, units = "mm", res = 600, compression = "lzw")
grid.newpage(); grid.draw(tab4_onlyHR); dev.off()

pdf(out_tab4_pdf, width = 7.8, height = 4.6)
grid.newpage(); grid.draw(tab4_onlyHR); dev.off()

############################################################
# Quick knobs (same as primary):
#  - pack figure up: fig4_heights = c(legend, km+risk, blank) -> increase blank
#  - shrink risk panel: km_vs_risk_heights -> make 2nd smaller
#  - risk row spacing: y_map values closer together
#  - title-to-table: risk_ylim upper (smaller -> closer)
############################################################
# ==========================================================
# Table (Figure 4) を Word（docx）でも保存（officer + flextable）
#  - タイトル + 罫線付きテーブル
#  - 1ページに収まりやすい余白/幅/フォント
# ==========================================================
pkgs_w <- c("officer","flextable")
to_install_w <- pkgs_w[!vapply(pkgs_w, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install_w) > 0) install.packages(to_install_w, dependencies = TRUE)

library(officer)
library(flextable)

out_tab4_docx <- file.path(out_dir, "Table_3_sensitivity_HR.docx")

# ---- flextable化（列幅/罫線/文字などを整える）----
ft4 <- flextable(hr_table)

ft4 <- ft4 %>%
  theme_booktabs() %>%                 # すっきり（罫線強めが良ければ theme_vanilla() に変更）
  bold(part = "header") %>%
  align(align = "left", part = "all") %>%
  align(j = c("HR","95% CI","p-value"), align = "center", part = "all") %>%
  autofit()

# 列幅（A4縦で読みやすい目安。必要なら微調整）
ft4 <- width(ft4, j = "Contrast", width = 3.8)
ft4 <- width(ft4, j = "HR",       width = 1.0)
ft4 <- width(ft4, j = "95% CI",   width = 1.6)
ft4 <- width(ft4, j = "p-value",  width = 1.0)

# フォントサイズ
ft4 <- fontsize(ft4, size = 11, part = "all")

# 外枠＋内枠（TableGrobの“外枠”相当）
outer <- fp_border(color = "black", width = 1)
inner <- fp_border(color = "black", width = 0.5)
ft4 <- border_remove(ft4)
ft4 <- border_outer(ft4, border = outer)
ft4 <- border_inner_h(ft4, border = inner)
ft4 <- border_inner_v(ft4, border = inner)

# ---- Word作成（A4縦・余白を詰めて1枚に収めやすく）----
doc4 <- read_docx()

sec4 <- prop_section(
  page_size = page_size(width = 8.27, height = 11.69),   # A4 (inch)
  page_margins = page_mar(top = 0.6, bottom = 0.6, left = 0.7, right = 0.7)
)

doc4 <- doc4 %>%
  body_add_par("Table 3. Adjusted hazard ratios for all-cause mortality",
               style = "heading 2") %>%
  body_add_flextable(ft4) %>%
  body_add_par("", style = "Normal") %>%
  body_end_section_continuous() %>%
  body_set_default_section(sec4)

print(doc4, target = out_tab4_docx)
} #感度分析

##併存疾患表の作成(SuppleT1,2)
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
setwd("X:/R")
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
} #データセットをエクセルにまとめ
{
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(officer)
library(flextable)

setwd("X:/R")
file_path <- "Supplementary_Table_1_collapsed.xlsx"
out_dir   <- "X:/R/word_supp_tables"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_doc1 <- file.path(out_dir, "Supplementary_Table_1_Drugs_portrait_readable.docx")
out_doc2 <- file.path(out_dir, "Supplementary_Table_2_ICD10_portrait_readable.docx")

# ---- Load ----
raw <- read_excel(file_path)
if (!all(c("term","definition") %in% names(raw))) raw <- raw %>% rename(term = 1, definition = 2)

dat <- raw %>%
  mutate(
    term = str_squish(as.character(term)),
    definition = str_squish(as.character(definition))
  ) %>%
  filter(!is.na(term), !is.na(definition), term != "", definition != "")

# ---- Drug terms ----
drug_terms <- c(
  "SGLT2 inhibitor",
  "Angiotensin II receptor blocker (ARB)",
  "Angiotensin receptor–neprilysin inhibitor (ARNI)",
  "Angiotensin-converting enzyme inhibitor (ACE inhibitor)",
  "SGLT2 inhibitor combination",
  "ARB + calcium channel blocker",
  "ARB + diuretic"
)

dat2 <- dat %>% mutate(is_drug = term %in% drug_terms)

drug_dat <- dat2 %>% filter(is_drug) %>% select(term, definition)
icd_dat  <- dat2 %>% filter(!is_drug) %>% select(term, definition)

# ---- Table 1: collapse by Drug class ----
tbl1 <- drug_dat %>%
  mutate(definition = str_replace_all(definition, ";", ",")) %>%
  separate_rows(definition, sep = "[,/]+") %>%
  mutate(Generic_name = str_squish(definition)) %>%
  filter(Generic_name != "") %>%
  transmute(`Drug class` = term, `Generic name` = Generic_name) %>%
  distinct() %>%
  group_by(`Drug class`) %>%
  summarise(`Generic name` = paste(sort(unique(`Generic name`)), collapse = ", "), .groups = "drop") %>%
  arrange(`Drug class`)

# ---- Table 2: collapse by Disease name (NO line breaks; keep original text) ----
tbl2 <- icd_dat %>%
  group_by(`Disease name` = term) %>%
  summarise(`ICD-10 code` = paste(unique(definition), collapse = "; "), .groups = "drop") %>%
  mutate(`ICD-10 code` = str_replace_all(`ICD-10 code`, "\\s*;\\s*", "; ")) %>%
  arrange(`Disease name`)

# ---- A4 portrait, margins slightly tight but not extreme ----
ps <- prop_section(
  page_size = page_size(width = 21.0, height = 29.7),
  page_margins = page_mar(top = 0.6, bottom = 0.6, left = 0.6, right = 0.6,
                          header = 0.2, footer = 0.2)
)

make_ft_readable <- function(df, caption_txt, col_widths, font_size = 10) {
  ft <- flextable(df) %>%
    set_caption(caption_txt) %>%
    bold(part = "header") %>%
    align(align = "left", part = "all") %>%
    valign(valign = "top", part = "all") %>%
    fontsize(size = font_size, part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    border_remove() %>%
    hline_top(border = fp_border(width = 1)) %>%
    hline(border = fp_border(width = 0.8), part = "header") %>%
    hline(i = seq_len(nrow(df)), border = fp_border(width = 0.35), part = "body") %>%
    hline_bottom(border = fp_border(width = 1)) %>%
    set_table_properties(layout = "fixed", width = 1)
  
  # 列幅固定（縦1枚で重要）
  for (nm in names(col_widths)) {
    if (nm %in% colnames(df)) ft <- width(ft, j = nm, width = col_widths[[nm]])
  }
  
  ft
}

# Table 1: 10ptで十分
ft1 <- make_ft_readable(
  tbl1,
  "Supplementary Table 1. Definitions of medications (drug class and generic names)",
  col_widths = c("Drug class" = 7.0, "Generic name" = 9.0),
  font_size = 10
)

doc1 <- read_docx() %>%
  body_set_default_section(ps) %>%
  body_add_flextable(ft1)
print(doc1, target = out_doc1)

# Table 2: Disease名は短め、ICD列を広く
ft2 <- make_ft_readable(
  tbl2,
  "Supplementary Table 2. Disease definitions using ICD-10 codes",
  col_widths = c("Disease name" = 6.0, "ICD-10 code" = 10.0),
  font_size = 10
)

doc2 <- read_docx() %>%
  body_set_default_section(ps) %>%
  body_add_flextable(ft2)
print(doc2, target = out_doc2)

message("Saved:")
message(out_doc1)
message(out_doc2)

library(RDCOMClient)

convert_docx_to_pdf <- function(docx_path, pdf_path) {
  docx_path <- normalizePath(docx_path, winslash = "\\", mustWork = TRUE)
  pdf_path  <- normalizePath(pdf_path,  winslash = "\\", mustWork = FALSE)
  
  # --- Word起動（既存があれば掴む、なければ作る） ---
  word <- NULL
  word <- tryCatch(COMGetActiveObject("Word.Application"), error = function(e) NULL)
  if (is.null(word)) {
    word <- COMCreate("Word.Application")
  }
  
  # Visible は環境によって無いことがあるので触らない（←今回の回避点）
  # word[["Visible"]] <- FALSE
  
  # --- docxを開く（ReadOnly, AddToRecentFiles=FALSE） ---
  docs <- word$Documents()
  doc  <- docs$Open(docx_path, ReadOnly = TRUE, AddToRecentFiles = FALSE)
  
  # --- PDF保存：ExportAsFixedFormat が最も安定 ---
  ok <- FALSE
  try({
    # 17 = wdExportFormatPDF
    # 0 = wdExportOptimizeForPrint
    doc$ExportAsFixedFormat(
      OutputFileName = pdf_path,
      ExportFormat   = 17,
      OpenAfterExport = FALSE,
      OptimizeFor     = 0
    )
    ok <- TRUE
  }, silent = TRUE)
  
  # --- だめなら SaveAs2 にフォールバック ---
  if (!ok) {
    try({
      # 17 = wdFormatPDF
      doc$SaveAs2(pdf_path, FileFormat = 17)
      ok <- TRUE
    }, silent = TRUE)
  }
  
  # --- 後片付け ---
  doc$Close(FALSE)
  
  # ここは「自分で起動したWordだけ閉じたい」けど判定が難しいので、
  # いったん Quit しない（Wordが勝手に閉じるのが嫌な場合）
  # 必要なら次行を有効化
  # word$Quit()
  
  if (!ok) stop("PDF conversion failed. (ExportAsFixedFormat / SaveAs2 both failed)")
  invisible(TRUE)
}

# ---- ファイル指定 ----
docx1 <- "X:/R/word_supp_tables/Supplementary_Table_1_Drugs_portrait_readable.docx"
docx2 <- "X:/R/word_supp_tables/Supplementary_Table_2_ICD10_portrait_readable.docx"

pdf1  <- "X:/R/word_supp_tables/Supplementary_Table_1_Drugs_portrait_readable.pdf"
pdf2  <- "X:/R/word_supp_tables/Supplementary_Table_2_ICD10_portrait_readable.pdf"

# ---- 実行 ----
convert_docx_to_pdf(docx1, pdf1)
convert_docx_to_pdf(docx2, pdf2)

message("PDF conversion completed.")

} #word,PDF変換


#すべての図表をまとめる
{
############################################################
# Merge selected PDFs (exact filenames from screenshot)
############################################################

# ---- 必要パッケージ ----
if (!requireNamespace("pdftools", quietly = TRUE)) {
  install.packages("pdftools")
}
library(pdftools)

# ---- フォルダ ----
dir_in  <- "X:/R/figure_table"
out_pdf <- file.path(dir_in, "Merged_Figures_All.pdf")

# ---- スクショ通りの正確なPDF名 ----
pdf_files <- c(
  "Figure1_flowchart_AKD.pdf",
  "Figure2_primary_KM.pdf",
  "Figure3A_observed_eGFR_trajectory_and_slope_1y.pdf",
  "Figure4_sensitivity_KM.pdf",
  "Figure5A_sens_observed_eGFR_trajectory_and_slope_1y.pdf",
  "SupplementalFigure1A_main_observed_eGFR_trajectory_and_bar_3y_all.pdf",
  "SupplementalFigure2A_sensitivity_observed_eGFR_trajectory_and_bar_3y_all.pdf",
  "SupplementalFigure3_interaction_forest_CKD_check.pdf",
  "Supplementary_Figure4_ForestPlot_Overall_AgeSubgroup.pdf",
  "Supplementary_Figure5_AgeContinuousInteraction.pdf"
)

# ---- フルパス化 ----
pdf_paths <- file.path(dir_in, pdf_files)

# ---- 存在チェック ----
missing <- pdf_paths[!file.exists(pdf_paths)]
if (length(missing)) {
  stop("These PDF files were not found:\n", paste(missing, collapse = "\n"))
}

# ---- 出力ファイルが開いていれば削除 ----
if (file.exists(out_pdf)) {
  file.remove(out_pdf)
}

# ---- 結合 ----
pdf_combine(pdf_paths, output = out_pdf)

cat("Merged successfully:\n", out_pdf)
} #FigureをPDFまとめ
{
############################################################
# Merge selected Word files into ONE document
#  - Switch section orientation (portrait/landscape) per file
#  - Add page break between documents
############################################################

pkgs <- c("officer","stringr","fs")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)

library(officer)
library(stringr)
library(fs)

# ---- folder ----
dir_in <- "X:/R/figure_table"   # ←あなたのフォルダ
stopifnot(dir_exists(dir_in))

out_docx <- file.path(dir_in, "Merged_Tables_Selected_fit.docx")

files <- c(
  "Supplementary_Table_1_Drugs_portrait_readable.docx",
  "Supplementary_Table_2_ICD10_portrait_readable.docx",
  "Table_3_sensitivity_HR.docx",
  "Table1_Baseline_with_CKD_1page.docx",
  "Table2_primary_HR.docx"
)

paths <- file.path(dir_in, files)
missing <- paths[!file.exists(paths)]
if (length(missing)) stop("Not found:\n", paste(missing, collapse = "\n"))

if (file.exists(out_docx)) {
  ok <- tryCatch(file.remove(out_docx), error = function(e) FALSE)
  if (!isTRUE(ok)) stop("Close the output docx first: ", out_docx)
}

# ==========================================================
# Section settings
#  - A4 portrait/landscape, slightly tight margins
# ==========================================================
sec_portrait <- prop_section(
  page_size = page_size(width = 8.27, height = 11.69),     # A4 portrait (inch)
  page_margins = page_mar(top = 0.55, bottom = 0.55, left = 0.55, right = 0.55,
                          header = 0.25, footer = 0.25)
)

sec_landscape <- prop_section(
  page_size = page_size(width = 11.69, height = 8.27),     # A4 landscape (inch)
  page_margins = page_mar(top = 0.45, bottom = 0.45, left = 0.45, right = 0.45,
                          header = 0.25, footer = 0.25)
)

# ---- “横向きにしたいファイル”を指定（広い表があるものだけ）----
# まずは Baseline table がはみ出すことが多いので landscape 推奨
landscape_files <- c(
  "Table1_Baseline_with_CKD_1page.docx"
)
# もし他もはみ出すならここに追加してください
# landscape_files <- c("Table1_Baseline_with_CKD_1page.docx","...")

# ==========================================================
# Merge
# ==========================================================
doc <- read_docx() |> body_set_default_section(sec_portrait)

for (i in seq_along(paths)) {
  
  f <- basename(paths[i])
  is_land <- f %in% landscape_files
  
  # ---- セクション切替（次に入れる文書に合わせる）----
  # ※「continuous section」で切替（ページは続くが向きだけ変わる）
  if (is_land) {
    doc <- doc |> body_end_section_continuous() |> body_set_default_section(sec_landscape)
  } else {
    doc <- doc |> body_end_section_continuous() |> body_set_default_section(sec_portrait)
  }
  
  # (optional) 見出し（不要ならコメントアウトOK）
  heading_txt <- str_replace(f, "\\.docx$", "")
  doc <- doc |>
    body_add_par(heading_txt, style = "heading 1") |>
    body_add_par("", style = "Normal")
  
  # ---- docx を挿入 ----
  doc <- body_add_docx(doc, src = paths[i])
  
  # ---- 次の文書との区切り：改ページ（最後以外）----
  if (i < length(paths)) {
    doc <- body_add_break(doc, pos = "after")
  }
}

print(doc, target = out_docx)
message("Saved: ", out_docx)
} #Tableをwordまとめ