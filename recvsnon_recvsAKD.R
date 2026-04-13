{
############################################################
# Figure 2 (KM only: Legend + KM + Number at risk)
#  + Table 2 (HR only, NO title)
#  --- ONE-PASS COPY-PASTE COMPLETE VERSION ---
############################################################

graphics.off()

# --------------------------
# Packages
# --------------------------
pkgs <- c("readr","dplyr","survival","broom","ggplot2",
          "grid","gridExtra","gtable","tibble","ragg",
          "officer","flextable")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

library(readr)
library(dplyr)
library(survival)
library(broom)
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
library(tibble)
library(ragg)
library(officer)
library(flextable)

# ==========================================================
# TUNING
# ==========================================================
base_fs <- 10
risk_fs <- 11
hr_fs   <- 11

km_line_lwd <- 1.05
km_ci_alpha <- 0.18

num_size <- 3.6
leg_text_size <- 3.2

# ---- Risk row spacing ----
y_map <- c(
  "non-AKD"              = 0.50,
  "AKD with recovery"    = 0.40,
  "AKD without recovery" = 0.30
)

# ---- Risk title spacing ----
risk_title_mb <- 2
risk_margin_t <- 10
risk_margin_b <- 2

# ---- Risk panel view window ----
risk_ylim <- c(0.18, 0.56)

# ---- Figure2 layout ----
km_vs_risk_heights <- c(2.00, 1.30)   # KM, risk
fig2_heights       <- c(0.32, 4.40, 1.20)  # legend, KM+risk, bottom blank

# ==========================================================
# Paths
# ==========================================================
setwd("X:/R")
in_csv <- "jin1_Eligibile.csv"

out_dir <- file.path("X:/R", "word_supp_tables")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

out_fig2_tif  <- file.path(out_dir, "Figure2_primary_KM.tif")
out_fig2_pdf  <- file.path(out_dir, "Figure2_primary_KM.pdf")
out_tab2_tif  <- file.path(out_dir, "Table2_primary_HR.tif")
out_tab2_pdf  <- file.path(out_dir, "Table2_primary_HR.pdf")
out_tab2_csv  <- file.path(out_dir, "Table2_primary_HR.csv")
out_tab2_docx <- file.path(out_dir, "Table2_primary_HR.docx")

# ==========================================================
# Load data
# ==========================================================
jin1_Eligibile <- read_csv(
  in_csv,
  locale = locale(encoding = "SHIFT-JIS"),
  show_col_types = FALSE
)
problems(jin1_Eligibile)

# 必要列チェック
req_cols <- c(
  "id","exclude","jin_status","index_date","date",
  "150_210recovery","90_150recovery",
  "last_follow_death","index_plus_210","primary_death",
  "age","index_cre","arb","acei",
  "dn1","dn3","dn4","dn5","dn6","dn7","dn8","dn9","dn10","dn12","dn13","dn14","dn15"
)
missing_cols <- setdiff(req_cols, names(jin1_Eligibile))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# ==========================================================
# Build analysis dataset (1 row per patient)
# ==========================================================
dat1 <- jin1_Eligibile %>%
  filter(exclude == "include", jin_status %in% c("AKD", "nonAKD")) %>%
  group_by(id) %>%
  arrange(index_date, date, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

levels_full <- c(
  "non-AKD",
  "AKD with recovery",
  "AKD without recovery"
)

dat_km <- dat1 %>%
  mutate(
    group = case_when(
      jin_status == "nonAKD" ~ "non-AKD",
      jin_status == "AKD" & `150_210recovery` == 1 ~ "AKD with recovery",
      jin_status == "AKD" & `150_210recovery` == 2 ~ "AKD without recovery",
      jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 1 ~ "AKD with recovery",
      jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` %in% c(0, 2) ~ "AKD without recovery",
      TRUE ~ NA_character_
    ),
    group = factor(group, levels = levels_full),
    time_years = as.numeric(last_follow_death - index_plus_210) / 365.25,
    primary_death = as.numeric(primary_death)
  ) %>%
  filter(
    !is.na(group),
    !is.na(time_years),
    time_years >= 0,
    !is.na(primary_death)
  )

# 念のため確認
if (nrow(dat_km) == 0) stop("dat_km has 0 rows after filtering.")
if (all(is.na(dat_km$group))) stop("All group values are NA.")
if (all(is.na(dat_km$time_years))) stop("All time_years are NA.")
if (all(is.na(dat_km$primary_death))) stop("All primary_death values are NA.")

print(table(dat_km$group, useNA = "ifany"))
print(summary(dat_km$time_years))
print(table(dat_km$primary_death, useNA = "ifany"))

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
print(fit)

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
    group = factor(sub("^group=", "", strata), levels = levels_full),
    cif   = 1 - surv,
    cif_l = 1 - upper,
    cif_u = 1 - lower
  ) %>%
  arrange(group, time)

p_km <- ggplot(km_df, aes(x = time, y = cif, colour = group, fill = group)) +
  geom_ribbon(
    aes(ymin = cif_l, ymax = cif_u),
    alpha = km_ci_alpha, linewidth = 0, show.legend = FALSE
  ) +
  geom_step(linewidth = km_line_lwd, show.legend = FALSE) +
  scale_x_continuous(
    breaks = ticks_show,
    limits = c(0, x_right),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 0.40),
    breaks = seq(0, 0.4, 0.1),
    expand = c(0, 0)
  ) +
  scale_colour_manual(values = pal, breaks = levels_full) +
  scale_fill_manual(values = pal, breaks = levels_full) +
  labs(x = NULL, y = "Cumulative incidence") +
  theme_classic(base_size = base_fs) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin  = margin(
      t = 4, r = common_right_margin_pt, b = 10,
      l = common_left_margin_pt, unit = "pt"
    )
  )

# ==========================================================
# Legend panel
# ==========================================================
leg_df <- tibble(
  group = factor(levels_full, levels = levels_full),
  x0 = c(1.2, 10.0, 20.0),
  y  = 1.10
) %>%
  mutate(x1 = x0 + 0.85)

p_leg <- ggplot(leg_df) +
  geom_rect(
    aes(xmin = x0, xmax = x1, ymin = y - 0.18, ymax = y + 0.18, fill = group),
    alpha = 0.25, colour = NA
  ) +
  geom_segment(
    aes(x = x0, xend = x1, y = y, yend = y, colour = group),
    linewidth = 1.1
  ) +
  geom_text(
    aes(x = x1 + 0.55, y = y, label = group),
    hjust = 0, size = leg_text_size
  ) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_colour_manual(values = pal, guide = "none") +
  coord_cartesian(xlim = c(0.5, 30.0), ylim = c(0.6, 1.7), clip = "off") +
  theme_void(base_size = base_fs) +
  theme(
    plot.margin = margin(
      t = 0, r = common_right_margin_pt, b = 0,
      l = common_left_margin_pt, unit = "pt"
    )
  )

# ==========================================================
# Number at risk
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
    group = factor(sub("^group=", "", strata), levels = levels_full),
    y = unname(y_map[as.character(group)]),
    time_plot = ifelse(time == 0, time0_shift, time),
    group_disp = case_when(
      as.character(group) == "AKD with recovery"    ~ "AKD with\nrecovery",
      as.character(group) == "AKD without recovery" ~ "AKD without\nrecovery",
      TRUE ~ as.character(group)
    )
  ) %>%
  filter(!is.na(group), !is.na(y))

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
    breaks = unname(y_map[levels_full]),
    labels = rep("", length(levels_full)),
    expand = c(0, 0)
  ) +
  labs(x = "Years", title = "Number at risk") +
  theme_classic(base_size = risk_fs) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    plot.margin  = margin(
      t = risk_margin_t, r = common_right_margin_pt,
      b = risk_margin_b, l = common_left_margin_pt, unit = "pt"
    ),
    plot.title.position = "plot",
    plot.title = element_text(margin = margin(b = risk_title_mb), vjust = 0)
  ) +
  coord_cartesian(ylim = risk_ylim, clip = "off")

for (grp in levels_full) {
  yv  <- unique(risk_df$y[risk_df$group == grp])[1]
  lab <- unique(risk_df$group_disp[risk_df$group == grp])[1]
  col <- if (grp == "non-AKD") "black" else pal[grp]
  
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
# Bind KM + risk
# ==========================================================
g_km   <- ggplotGrob(p_km)
g_risk <- ggplotGrob(p_risk)

g_risk$widths <- g_km$widths

km_risk_block <- arrangeGrob(
  g_km, g_risk,
  ncol = 1,
  heights = km_vs_risk_heights
)

# ==========================================================
# Figure 2 object
# ==========================================================
fig2_onlyKM <- arrangeGrob(
  p_leg,
  km_risk_block,
  nullGrob(),
  ncol = 1,
  heights = fig2_heights
)

{
  library(dplyr)
  library(readr)
  
  jin1_Eligibile <- read_csv(
    "X:/R/jin1_Eligibile.csv",
    locale = locale(encoding = "SHIFT-JIS"),
    show_col_types = FALSE
  )
  
  # -----------------------------
  # Figure 1 相当
  # -----------------------------
  flow_df <- jin1_Eligibile %>%
    distinct(id, .keep_all = TRUE) %>%
    filter(exclude == "include") %>%
    mutate(
      flow_group = case_when(
        jin_status == "nonAKD" ~ "non-AKD",
        jin_status == "AKD" & `150_210recovery` == 1 ~ "AKD with recovery",
        jin_status == "AKD" & `150_210recovery` == 2 ~ "AKD without recovery",
        jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 1 ~ "AKD with recovery",
        jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` %in% c(0, 2) ~ "AKD without recovery",
        TRUE ~ NA_character_
      )
    )
  
  table(flow_df$flow_group, useNA = "ifany")
  
  # -----------------------------
  # Figure 2 相当
  # -----------------------------
  levels_full <- c("non-AKD", "AKD with recovery", "AKD without recovery")
  
  dat1 <- jin1_Eligibile %>%
    filter(exclude == "include", jin_status %in% c("AKD", "nonAKD")) %>%
    group_by(id) %>%
    arrange(index_date, date, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()
  
  dat_km_pre <- dat1 %>%
    mutate(
      group = case_when(
        jin_status == "nonAKD" ~ "non-AKD",
        jin_status == "AKD" & `150_210recovery` == 1 ~ "AKD with recovery",
        jin_status == "AKD" & `150_210recovery` == 2 ~ "AKD without recovery",
        jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 1 ~ "AKD with recovery",
        jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` %in% c(0, 2) ~ "AKD without recovery",
        TRUE ~ NA_character_
      ),
      time_years = as.numeric(last_follow_death - index_plus_210) / 365.25,
      primary_death = as.numeric(primary_death)
    )
  
  # どこで落ちているか確認
  dat_km_pre %>%
    mutate(
      reason = case_when(
        is.na(group) ~ "group missing",
        is.na(time_years) ~ "time_years missing",
        time_years < 0 ~ "time_years < 0",
        is.na(primary_death) ~ "primary_death missing",
        TRUE ~ "kept"
      )
    ) %>%
    count(group, reason)
  
  # 除外されたID一覧
  excluded_from_km <- dat_km_pre %>%
    mutate(
      reason = case_when(
        is.na(group) ~ "group missing",
        is.na(time_years) ~ "time_years missing",
        time_years < 0 ~ "time_years < 0",
        is.na(primary_death) ~ "primary_death missing",
        TRUE ~ "kept"
      )
    ) %>%
    filter(reason != "kept") %>%
    select(id, group, last_follow_death, index_plus_210, time_years, primary_death, reason)
  
  excluded_from_km
  excluded_from_km %>%
    group_by(group) %>%
    summarise(
      n_excluded = n_distinct(id),
      .groups = "drop"
    )
}#Number at riskで人が減ることの確認


# ==========================================================
# Cox model -> HR table
# ==========================================================
dat_cox <- dat_km %>%
  mutate(
    arb_acei_use = if_else(coalesce(arb, 0) == 1 | coalesce(acei, 0) == 1, 1L, 0L)
  )

fit_main <- coxph(
  Surv(time_years, primary_death) ~
    group + age + index_cre + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
  data = dat_cox
)

hr_table <- broom::tidy(fit_main, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term %in% c("groupAKD with recovery", "groupAKD without recovery")) %>%
  mutate(
    Contrast = c(
      "AKD with recovery vs non-AKD",
      "AKD without recovery vs non-AKD"
    ),
    HR = sprintf("%.2f", estimate),
    `95% CI` = sprintf("%.2f–%.2f", conf.low, conf.high),
    `p-value` = ifelse(p.value < 0.01, "<0.01", sprintf("%.2f", p.value))
  ) %>%
  select(Contrast, HR, `95% CI`, `p-value`)

write_csv(hr_table, out_tab2_csv)

# ==========================================================
# Table 2 grob (NO TITLE)
# ==========================================================
hr_grob <- tableGrob(
  hr_table,
  rows = NULL,
  theme = ttheme_default(
    base_size = hr_fs,
    core = list(bg_params = list(col = "black", lwd = 0.4)),
    colhead = list(
      bg_params = list(col = "black", lwd = 0.6),
      fg_params = list(fontface = "bold")
    )
  )
)

hr_grob <- gtable_add_grob(
  hr_grob,
  rectGrob(gp = gpar(fill = NA, col = "black", lwd = 1)),
  t = 1, l = 1, b = nrow(hr_grob), r = ncol(hr_grob)
)

# ★ タイトルなし
tab2_onlyHR <- arrangeGrob(
  hr_grob,
  ncol = 1
)

# ==========================================================
# Draw (optional)
# ==========================================================
grid.newpage()
grid.draw(fig2_onlyKM)

grid.newpage()
grid.draw(tab2_onlyHR)

# ==========================================================
# Save Figure 2 / Table 2
# ==========================================================
ragg::agg_tiff(
  out_fig2_tif,
  width = 180, height = 220, units = "mm",
  res = 600, compression = "lzw"
)
grid.newpage()
grid.draw(fig2_onlyKM)
dev.off()

pdf(out_fig2_pdf, width = 7.8, height = 8.4)
grid.newpage()
grid.draw(fig2_onlyKM)
dev.off()

ragg::agg_tiff(
  out_tab2_tif,
  width = 180, height = 120, units = "mm",
  res = 600, compression = "lzw"
)
grid.newpage()
grid.draw(tab2_onlyHR)
dev.off()

pdf(out_tab2_pdf, width = 7.8, height = 4.6)
grid.newpage()
grid.draw(tab2_onlyHR)
dev.off()

# ==========================================================
# Table 2 -> Word (NO TITLE)
# ==========================================================
ft <- flextable(hr_table)

ft <- ft %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  align(align = "left", part = "all") %>%
  align(j = c("HR", "95% CI", "p-value"), align = "center", part = "all") %>%
  autofit()

ft <- width(ft, j = "Contrast", width = 3.6)
ft <- width(ft, j = "HR",       width = 1.0)
ft <- width(ft, j = "95% CI",   width = 1.6)
ft <- width(ft, j = "p-value",  width = 1.0)

ft <- fontsize(ft, size = 11, part = "all")

outer <- fp_border(color = "black", width = 1)
inner <- fp_border(color = "black", width = 0.5)

ft <- border_remove(ft)
ft <- border_outer(ft, border = outer)
ft <- border_inner_h(ft, border = inner)
ft <- border_inner_v(ft, border = inner)

doc <- read_docx()

sec <- prop_section(
  page_size = page_size(width = 8.27, height = 11.69),
  page_margins = page_mar(
    top = 0.6, bottom = 0.6, left = 0.7, right = 0.7
  )
)

doc <- doc %>%
  body_add_flextable(ft) %>%
  body_add_par("", style = "Normal") %>%
  body_end_section_continuous() %>%
  body_set_default_section(sec)

print(doc, target = out_tab2_docx)

############################################################
# Quick knobs
#  - pack figure up: fig2_heights = c(legend, km+risk, blank)
#  - shrink risk panel: km_vs_risk_heights
#  - risk row spacing: y_map values
#  - title-to-table: risk_ylim upper/lower
############################################################
# ==========================================================
# AKD群内：recovery vs non-recovery の直接比較
# ==========================================================
dat_akd_cox <- dat_cox %>%
  filter(group %in% c("AKD with recovery", "AKD without recovery")) %>%
  mutate(
    group_akd = factor(
      group,
      levels = c("AKD without recovery", "AKD with recovery")
    )
  )

fit_akd_internal <- coxph(
  Surv(time_years, primary_death) ~
    group_akd + age + index_cre + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
  data = dat_akd_cox
)

summary(fit_akd_internal)

# HR table
hr_table_akd_internal <- broom::tidy(
  fit_akd_internal,
  exponentiate = TRUE,
  conf.int = TRUE
) %>%
  filter(term == "group_akdAKD with recovery") %>%
  mutate(
    Contrast = "AKD with recovery vs AKD without recovery",
    HR = sprintf("%.2f", estimate),
    `95% CI` = sprintf("%.2f–%.2f", conf.low, conf.high),
    `p-value` = ifelse(p.value < 0.01, "<0.01", sprintf("%.2f", p.value))
  ) %>%
  select(Contrast, HR, `95% CI`, `p-value`)

print(hr_table_akd_internal)  


} #本解析


{
  ############################################################
  # Figure 4 (Sensitivity)  --- ONE-PASS COPY-PASTE COMPLETE VERSION ---
  # Revised:
  #  - Sensitivity grouping includes "No-data"
  #  - BUT "No-data" is excluded from KM figure and Cox model
  #  - Figure layout is identical to Primary Figure 2 code
  #  - Figure4: Legend + KM(CIF+CI) + Number at risk
  #  - HR table saved separately + CSV + Word
  #  - Table 3 title REMOVED from TIFF/PDF/Word
  # Output folder: X:/R/sensitivity_analysis/
  ############################################################
  
  graphics.off()
  
  # --------------------------
  # Packages
  # --------------------------
  pkgs <- c("readr","dplyr","survival","broom","ggplot2",
            "grid","gridExtra","gtable","tibble","ragg",
            "officer","flextable")
  to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
  
  library(readr)
  library(dplyr)
  library(survival)
  library(broom)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(gtable)
  library(tibble)
  library(ragg)
  library(officer)
  library(flextable)
  
  # ==========================================================
  # TUNING
  # ==========================================================
  base_fs <- 10
  risk_fs <- 11
  hr_fs   <- 11
  
  km_line_lwd <- 1.05
  km_ci_alpha <- 0.18
  
  num_size <- 3.6
  leg_text_size <- 3.2
  
  # ---- Risk row spacing ----
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
  # Paths
  # ==========================================================
  setwd("X:/R")
  in_csv <- "jin1_Eligibile.csv"

  out_dir <- file.path("X:/R", "word_supp_tables")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  out_fig4_tif  <- file.path(out_dir, "Supplementary Figure1 Sensitivity Analysis KM.tif")
  out_fig4_pdf  <- file.path(out_dir, "Supplementary Figure1 Sensitivity Analysis KM.pdf")
  
  out_tab4_tif  <- file.path(out_dir, "Supplementary Table_3_sensitivity_HR.tif")
  out_tab4_pdf  <- file.path(out_dir, "Supplementary Table_3_sensitivity_HR.pdf")
  out_tab4_csv  <- file.path(out_dir, "Supplementary Table_3_sensitivity_HR.csv")
  out_tab4_docx <- file.path(out_dir, "Table_3_sensitivity_HR.docx")
  
  # ==========================================================
  # Load data
  # ==========================================================
  jin1_Eligibile <- read_csv(
    in_csv,
    locale = locale(encoding = "SHIFT-JIS"),
    show_col_types = FALSE
  )
  
  # 必要列チェック
  req_cols <- c(
    "id","exclude","jin_status","index_date","date",
    "150_210recovery","90_150recovery",
    "last_follow_death","index_plus_210","primary_death",
    "age","index_cre","arb","acei",
    "dn1","dn3","dn4","dn5","dn6","dn7","dn8","dn9","dn10","dn12","dn13","dn14","dn15"
  )
  missing_cols <- setdiff(req_cols, names(jin1_Eligibile))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # ==========================================================
  # Build 1 row per patient
  # ==========================================================
  dat1 <- jin1_Eligibile %>%
    filter(exclude == "include", jin_status %in% c("AKD", "nonAKD")) %>%
    group_by(id) %>%
    arrange(index_date, date, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()
  
  # ==========================================================
  # Sensitivity grouping
  #   - Define No-data as a separate category
  #   - But exclude No-data from figure / Cox model
  # ==========================================================
  dat_sens_all <- dat1 %>%
    mutate(
      jin_label_sens = case_when(
        jin_status == "nonAKD" ~ "nonAKD",
        jin_status == "AKD" & `150_210recovery` == 1 ~ "Recovery",
        jin_status == "AKD" & `150_210recovery` == 2 ~ "Non-Recovery",
        jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 1 ~ "Recovery",
        jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 2 ~ "Non-Recovery",
        jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 0 ~ "No-data",
        TRUE ~ NA_character_
      ),
      jin_label_sens = factor(
        jin_label_sens,
        levels = c("nonAKD", "Recovery", "Non-Recovery", "No-data")
      ),
      arb_acei_use = if_else(coalesce(arb, 0) == 1 | coalesce(acei, 0) == 1, 1L, 0L),
      time_years = as.numeric(last_follow_death - index_plus_210) / 365.25,
      primary_death = as.numeric(primary_death)
    ) %>%
    filter(
      !is.na(jin_label_sens),
      !is.na(time_years),
      time_years >= 0,
      !is.na(primary_death)
    )
  
  # ---- Figure / Cox use only 3 groups ----
  levels_full <- c("nonAKD", "Recovery", "Non-Recovery")
  
  dat_sens <- dat_sens_all %>%
    filter(jin_label_sens != "No-data") %>%
    mutate(
      group = factor(as.character(jin_label_sens), levels = levels_full)
    ) %>%
    filter(!is.na(group))
  
  # 確認
  if (nrow(dat_sens) == 0) stop("dat_sens has 0 rows after filtering.")
  if (all(is.na(dat_sens$group))) stop("All group values are NA.")
  if (all(is.na(dat_sens$time_years))) stop("All time_years are NA.")
  if (all(is.na(dat_sens$primary_death))) stop("All primary_death values are NA.")
  
  print(table(dat_sens$group, useNA = "ifany"))
  print(summary(dat_sens$time_years))
  print(table(dat_sens$primary_death, useNA = "ifany"))
  
  # ==========================================================
  # Common axis / colors
  # ==========================================================
  ticks_show   <- seq(0, 8, by = 2)
  x_right      <- 9.5
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
  print(fit)
  
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
      group = factor(sub("^group=", "", strata), levels = levels_full),
      cif   = 1 - surv,
      cif_l = 1 - upper,
      cif_u = 1 - lower
    ) %>%
    arrange(group, time)
  
  p_km <- ggplot(km_df, aes(x = time, y = cif, colour = group, fill = group)) +
    geom_ribbon(
      aes(ymin = cif_l, ymax = cif_u),
      alpha = km_ci_alpha, linewidth = 0, show.legend = FALSE
    ) +
    geom_step(linewidth = km_line_lwd, show.legend = FALSE) +
    scale_x_continuous(
      breaks = ticks_show,
      limits = c(0, x_right),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 0.40),
      breaks = seq(0, 0.4, 0.1),
      expand = c(0, 0)
    ) +
    scale_colour_manual(values = pal, breaks = levels_full) +
    scale_fill_manual(values = pal, breaks = levels_full) +
    labs(x = NULL, y = "Cumulative incidence") +
    theme_classic(base_size = base_fs) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin  = margin(
        t = 4, r = common_right_margin_pt, b = 10,
        l = common_left_margin_pt, unit = "pt"
      )
    )
  
  # ==========================================================
  # Legend panel
  # ==========================================================
  leg_levels_disp <- c(
    "nonAKD"       = "nonAKD",
    "Recovery"     = "AKD with recovery",
    "Non-Recovery" = "AKD without recovery"
  )
  
  leg_df <- tibble(
    group = factor(levels_full, levels = levels_full),
    x0    = c(1.2, 10.0, 20.0),
    y     = 1.10,
    label = unname(leg_levels_disp[levels_full])
  ) %>%
    mutate(x1 = x0 + 0.85)
  
  p_leg <- ggplot(leg_df) +
    geom_rect(
      aes(xmin = x0, xmax = x1, ymin = y - 0.18, ymax = y + 0.18, fill = group),
      alpha = 0.25, colour = NA
    ) +
    geom_segment(
      aes(x = x0, xend = x1, y = y, yend = y, colour = group),
      linewidth = 1.1
    ) +
    geom_text(
      aes(x = x1 + 0.55, y = y, label = label),
      hjust = 0, size = leg_text_size
    ) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_colour_manual(values = pal, guide = "none") +
    coord_cartesian(xlim = c(0.5, 30.0), ylim = c(0.6, 1.7), clip = "off") +
    theme_void(base_size = base_fs) +
    theme(
      plot.margin = margin(
        t = 0, r = common_right_margin_pt, b = 0,
        l = common_left_margin_pt, unit = "pt"
      )
    )
  
  # ==========================================================
  # Number at risk
  # ==========================================================
  sfit <- summary(fit, times = ticks_show, extend = TRUE)
  
  time0_shift  <- 0.22
  label_pad_mm <- 10
  
  risk_df <- data.frame(
    time   = sfit$time,
    strata = sfit$strata,
    n_risk = sfit$n.risk
  ) %>%
    mutate(
      group = factor(sub("^group=", "", strata), levels = levels_full),
      y = unname(y_map[as.character(group)]),
      time_plot = ifelse(time == 0, time0_shift, time),
      group_disp = case_when(
        as.character(group) == "Recovery"     ~ "AKD with\nrecovery",
        as.character(group) == "Non-Recovery" ~ "AKD without\nrecovery",
        TRUE ~ as.character(group)
      )
    ) %>%
    filter(!is.na(group), !is.na(y))
  
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
      breaks = unname(y_map[levels_full]),
      labels = rep("", length(levels_full)),
      expand = c(0, 0)
    ) +
    labs(x = "Years", title = "Number at risk") +
    theme_classic(base_size = risk_fs) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y  = element_blank(),
      plot.margin  = margin(
        t = risk_margin_t, r = common_right_margin_pt,
        b = risk_margin_b, l = common_left_margin_pt, unit = "pt"
      ),
      plot.title.position = "plot",
      plot.title = element_text(margin = margin(b = risk_title_mb), vjust = 0)
    ) +
    coord_cartesian(ylim = risk_ylim, clip = "off")
  
  for (grp in levels_full) {
    yv  <- unique(risk_df$y[risk_df$group == grp])[1]
    lab <- unique(risk_df$group_disp[risk_df$group == grp])[1]
    col <- if (grp == "nonAKD") "black" else pal[grp]
    
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
  # Bind KM + risk with controlled vertical ratio
  # ==========================================================
  g_km   <- ggplotGrob(p_km)
  g_risk <- ggplotGrob(p_risk)
  
  g_risk$widths <- g_km$widths
  
  km_risk_block <- arrangeGrob(
    g_km, g_risk,
    ncol = 1,
    heights = km_vs_risk_heights
  )
  
  # ==========================================================
  # Figure 4 object
  # ==========================================================
  fig4_onlyKM <- arrangeGrob(
    p_leg,
    km_risk_block,
    nullGrob(),
    ncol = 1,
    heights = fig4_heights
  )
  
  # ==========================================================
  # Cox -> HR table
  #   - No-data excluded
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
      Contrast = c(
        "AKD with recovery vs nonAKD",
        "AKD without recovery vs nonAKD"
      ),
      HR = sprintf("%.2f", estimate),
      `95% CI` = sprintf("%.2f–%.2f", conf.low, conf.high),
      `p-value` = ifelse(p.value < 0.01, "<0.01", sprintf("%.2f", p.value))
    ) %>%
    select(Contrast, HR, `95% CI`, `p-value`)
  
  write_csv(hr_table, out_tab4_csv)
  
  # ==========================================================
  # Table 3 grob (NO TITLE)
  # ==========================================================
  hr_grob <- tableGrob(
    hr_table,
    rows = NULL,
    theme = ttheme_default(
      base_size = hr_fs,
      core = list(bg_params = list(col = "black", lwd = 0.4)),
      colhead = list(
        bg_params = list(col = "black", lwd = 0.6),
        fg_params = list(fontface = "bold")
      )
    )
  )
  
  hr_grob <- gtable_add_grob(
    hr_grob,
    rectGrob(gp = gpar(fill = NA, col = "black", lwd = 1)),
    t = 1, l = 1, b = nrow(hr_grob), r = ncol(hr_grob)
  )
  
  # ★ タイトルなし
  tab4_onlyHR <- arrangeGrob(
    hr_grob,
    ncol = 1
  )
  
  # ==========================================================
  # Draw (optional check)
  # ==========================================================
  grid.newpage()
  grid.draw(fig4_onlyKM)
  
  grid.newpage()
  grid.draw(tab4_onlyHR)
  
  # ==========================================================
  # SAVE
  # ==========================================================
  ragg::agg_tiff(
    out_fig4_tif,
    width = 180, height = 220, units = "mm",
    res = 600, compression = "lzw"
  )
  grid.newpage()
  grid.draw(fig4_onlyKM)
  dev.off()
  
  pdf(out_fig4_pdf, width = 7.8, height = 8.4)
  grid.newpage()
  grid.draw(fig4_onlyKM)
  dev.off()
  
  ragg::agg_tiff(
    out_tab4_tif,
    width = 180, height = 120, units = "mm",
    res = 600, compression = "lzw"
  )
  grid.newpage()
  grid.draw(tab4_onlyHR)
  dev.off()
  
  pdf(out_tab4_pdf, width = 7.8, height = 4.6)
  grid.newpage()
  grid.draw(tab4_onlyHR)
  dev.off()
  
  # ==========================================================
  # Table 3 also save as Word (docx) --- NO TITLE
  # ==========================================================
  ft4 <- flextable(hr_table)
  
  ft4 <- ft4 %>%
    theme_booktabs() %>%
    bold(part = "header") %>%
    align(align = "left", part = "all") %>%
    align(j = c("HR", "95% CI", "p-value"), align = "center", part = "all") %>%
    autofit()
  
  ft4 <- width(ft4, j = "Contrast", width = 3.8)
  ft4 <- width(ft4, j = "HR",       width = 1.0)
  ft4 <- width(ft4, j = "95% CI",   width = 1.6)
  ft4 <- width(ft4, j = "p-value",  width = 1.0)
  
  ft4 <- fontsize(ft4, size = 11, part = "all")
  
  outer <- fp_border(color = "black", width = 1)
  inner <- fp_border(color = "black", width = 0.5)
  
  ft4 <- border_remove(ft4)
  ft4 <- border_outer(ft4, border = outer)
  ft4 <- border_inner_h(ft4, border = inner)
  ft4 <- border_inner_v(ft4, border = inner)
  
  doc4 <- read_docx()
  
  sec4 <- prop_section(
    page_size = page_size(width = 8.27, height = 11.69),
    page_margins = page_mar(top = 0.6, bottom = 0.6, left = 0.7, right = 0.7)
  )
  
  # ★ タイトル行なし
  doc4 <- doc4 %>%
    body_add_flextable(ft4) %>%
    body_add_par("", style = "Normal") %>%
    body_end_section_continuous() %>%
    body_set_default_section(sec4)
  
  print(doc4, target = out_tab4_docx)
  
  ############################################################
  # Quick knobs
  #  - pack figure up: fig4_heights = c(legend, km+risk, blank)
  #  - shrink risk panel: km_vs_risk_heights
  #  - risk row spacing: y_map values
  #  - title-to-table: risk_ylim upper/lower
  ############################################################
  # ==========================================================
  # Sensitivity analysis:
  # Direct comparison within AKD
  # Recovery vs Non-Recovery
  # ==========================================================
  dat_sens_akd <- dat_sens %>%
    filter(group %in% c("Recovery", "Non-Recovery")) %>%
    mutate(
      group_akd = factor(group, levels = c("Non-Recovery", "Recovery"))
    )
  
  fit_cox_akd_internal_sens <- coxph(
    Surv(time_years, primary_death) ~
      group_akd + age + index_cre + arb_acei_use +
      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
    data = dat_sens_akd
  )
  
  summary(fit_cox_akd_internal_sens)
  
  hr_table_akd_internal_sens <- broom::tidy(
    fit_cox_akd_internal_sens,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    filter(term == "group_akdRecovery") %>%
    mutate(
      Contrast = "AKD with recovery vs AKD without recovery",
      HR = sprintf("%.2f", estimate),
      `95% CI` = sprintf("%.2f–%.2f", conf.low, conf.high),
      `p-value` = ifelse(p.value < 0.01, "<0.01", sprintf("%.2f", p.value))
    ) %>%
    select(Contrast, HR, `95% CI`, `p-value`)
  
  print(hr_table_akd_internal_sens)
  
  
  
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
    #set_caption(caption_txt) %>%
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
  # =========================================================
  # Supplementary Table 2 (Drugs) - final submission style
  # - Wider Drug class column
  # - Reduced left padding in Generic name column
  # - Minimized wrapping
  # - Save as DOCX and PDF
  # =========================================================
  
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(officer)
  library(flextable)
  library(RDCOMClient)
  
  setwd("X:/R")
  
  # =========================================================
  # Paths
  # =========================================================
  file_path <- "Supplementary_Table_1_collapsed.xlsx"
  out_dir   <- file.path("X:/R", "word_supp_tables")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  out_docx <- file.path(out_dir, "Supplementary_Table_2_Drugs_final.docx")
  out_pdf  <- file.path(out_dir, "Supplementary_Table_2_Drugs_final.pdf")
  
  # =========================================================
  # Load
  # =========================================================
  raw <- read_excel(file_path)
  
  if (!all(c("term","definition") %in% names(raw))) {
    raw <- raw %>% rename(term = 1, definition = 2)
  }
  
  dat <- raw %>%
    mutate(
      term       = str_squish(as.character(term)),
      definition = str_squish(as.character(definition))
    ) %>%
    filter(!is.na(term), !is.na(definition), term != "", definition != "")
  
  # =========================================================
  # Drug terms
  # =========================================================
  drug_terms <- c(
    "SGLT2 inhibitor",
    "Angiotensin II receptor blocker (ARB)",
    "Angiotensin receptor–neprilysin inhibitor (ARNI)",
    "Angiotensin-converting enzyme inhibitor (ACE inhibitor)",
    "SGLT2 inhibitor combination",
    "ARB + calcium channel blocker",
    "ARB + diuretic"
  )
  
  drug_dat <- dat %>%
    filter(term %in% drug_terms) %>%
    select(term, definition)
  
  # =========================================================
  # Build Supplementary Table 2
  # =========================================================
  tbl2 <- drug_dat %>%
    mutate(definition = str_replace_all(definition, ";", ",")) %>%
    separate_rows(definition, sep = "[,/]+") %>%
    mutate(Generic_name = str_squish(definition)) %>%
    filter(Generic_name != "") %>%
    transmute(
      `Drug class`   = term,
      `Generic name` = Generic_name
    ) %>%
    distinct() %>%
    group_by(`Drug class`) %>%
    summarise(
      `Generic name` = paste(sort(unique(`Generic name`)), collapse = ", "),
      .groups = "drop"
    ) %>%
    arrange(`Drug class`)
  
  # =========================================================
  # Optional: prevent awkward line breaks around hyphens
  # =========================================================
  # =========================================================
  # Build Supplementary Table 2
  # 配合剤は「1つの合剤単位」で残す
  # =========================================================
  combo_terms <- c(
    "ARB + calcium channel blocker",
    "ARB + diuretic",
    "Angiotensin receptor–neprilysin inhibitor (ARNI)",
    "SGLT2 inhibitor combination"
  )
  
  tbl2 <- drug_dat %>%
    mutate(
      definition = str_squish(definition),
      definition = str_replace_all(definition, "；", ";"),
      definition = str_replace_all(definition, "\\s*;\\s*", "; ")
    ) %>%
    rowwise() %>%
    mutate(
      generic_clean = if (`term` %in% combo_terms) {
        # 合剤は / や , で壊さない
        # まず ; 単位で候補を分ける
        parts <- unlist(str_split(definition, "\\s*;\\s*"))
        
        # もし ; がなく、, で複数製品が並んでいそうならそのまま1塊として残す
        parts <- parts[parts != ""]
        
        # 表示を整える
        parts <- str_replace_all(parts, "\\s*/\\s*", " + ")
        parts <- str_replace_all(parts, "\\s*,\\s*", " + ")
        
        # 重複除去して並べる
        paste(unique(parts), collapse = "; ")
      } else {
        # 単剤は従来通り成分を展開
        parts <- unlist(str_split(definition, "\\s*[,/]\\s*"))
        parts <- str_squish(parts)
        parts <- parts[parts != ""]
        paste(sort(unique(parts)), collapse = ", ")
      }
    ) %>%
    ungroup() %>%
    transmute(
      `Drug class`   = term,
      `Generic name` = generic_clean
    ) %>%
    distinct() %>%
    group_by(`Drug class`) %>%
    summarise(
      `Generic name` = paste(unique(`Generic name`), collapse = "; "),
      .groups = "drop"
    ) %>%
    arrange(`Drug class`) %>%
    mutate(
      `Drug class` = str_replace_all(`Drug class`, "-", "\u2011")
    )
  # =========================================================
  # Flextable helper
  # =========================================================
  make_ft_supp_drug_final <- function(df) {
    ft <- flextable(df)
    
    ft <- ft %>%
      bold(part = "header") %>%
      align(align = "left", part = "all") %>%
      valign(valign = "top", part = "all") %>%
      font(fontname = "Times New Roman", part = "all") %>%
      fontsize(size = 10, part = "body") %>%
      fontsize(size = 10.5, part = "header") %>%
      line_spacing(space = 1.0, part = "all") %>%
      border_remove() %>%
      hline_top(border = fp_border(width = 1.0)) %>%
      hline(part = "header", border = fp_border(width = 0.8)) %>%
      hline(i = seq_len(nrow(df)), part = "body", border = fp_border(width = 0.35)) %>%
      hline_bottom(border = fp_border(width = 1.0)) %>%
      set_table_properties(layout = "fixed", width = 1)
    
    # ---- Column widths ----
    # Drug class を広めに確保
    # Generic name は最大限活かす
    ft <- ft %>%
      width(j = "Drug class",   width = 3.3) %>%
      width(j = "Generic name", width = 7.7)
    
    # ---- Padding ----
    # Generic name の左余白を小さく
    ft <- ft %>%
      padding(part = "all", padding.top = 1.5, padding.bottom = 1.5,
              padding.left = 2.5, padding.right = 2.5) %>%
      padding(j = "Drug class",   part = "all", padding.left = 2.2, padding.right = 2.0) %>%
      padding(j = "Generic name", part = "all", padding.left = 0.8, padding.right = 1.5)
    
    # ---- Header style ----
    ft <- ft %>%
      bg(part = "header", bg = "white")
    
    # ---- Row height ----
    ft <- ft %>%
      height_all(height = 0.24)
    
    ft
  }
  # =========================================================
  # Build Supplementary Table 1
  # =========================================================
  icd_dat <- dat %>%
    filter(!term %in% drug_terms) %>%
    select(term, definition)
  
  tbl1 <- icd_dat %>%
    group_by(`Disease name` = term) %>%
    summarise(
      `ICD-10 code` = paste(unique(definition), collapse = "; "),
      .groups = "drop"
    ) %>%
    mutate(
      `ICD-10 code` = str_replace_all(`ICD-10 code`, "\\s*;\\s*", "; ")
    ) %>%
    arrange(`Disease name`)
  
  # =========================================================
  # Flextable helper for Supplementary Table 1
  # =========================================================
  make_ft_supp_icd_final <- function(df) {
    ft <- flextable(df)
    
    ft <- ft %>%
      bold(part = "header") %>%
      align(align = "left", part = "all") %>%
      valign(valign = "top", part = "all") %>%
      font(fontname = "Times New Roman", part = "all") %>%
      fontsize(size = 10, part = "body") %>%
      fontsize(size = 10.5, part = "header") %>%
      line_spacing(space = 1.0, part = "all") %>%
      border_remove() %>%
      hline_top(border = fp_border(width = 1.0)) %>%
      hline(part = "header", border = fp_border(width = 0.8)) %>%
      hline(i = seq_len(nrow(df)), part = "body", border = fp_border(width = 0.35)) %>%
      hline_bottom(border = fp_border(width = 1.0)) %>%
      set_table_properties(layout = "fixed", width = 1)
    
    # ---- Column widths ----
    ft <- ft %>%
      width(j = "Disease name", width = 3.8) %>%
      width(j = "ICD-10 code",  width = 7.2)
    
    # ---- Padding ----
    ft <- ft %>%
      padding(part = "all", padding.top = 1.5, padding.bottom = 1.5,
              padding.left = 1.5, padding.right = 1.5) %>%
      padding(j = "Disease name", part = "all", padding.left = 1.2, padding.right = 1.2) %>%
      padding(j = "ICD-10 code",  part = "all", padding.left = 0.6, padding.right = 1.0)
    
    ft <- ft %>%
      bg(part = "header", bg = "white") %>%
      height_all(height = 0.24)
    
    ft
  }
  
  # =========================================================
  # Create flextables
  # =========================================================
  ft1 <- make_ft_supp_icd_final(tbl1)
  ft2 <- make_ft_supp_drug_final(tbl2)
  
  # =========================================================
  # Output paths
  # =========================================================
  out_docx1 <- file.path(out_dir, "Supplementary_Table_1_ICD10_final.docx")
  out_pdf1  <- file.path(out_dir, "Supplementary_Table_1_ICD10_final.pdf")
  
  out_docx2 <- file.path(out_dir, "Supplementary_Table_2_Drugs_final.docx")
  out_pdf2  <- file.path(out_dir, "Supplementary_Table_2_Drugs_final.pdf")
  
  # =========================================================
  # A4 portrait section
  # 余白はかなり狭めて、横切れを防ぐ
  # =========================================================
  sec <- prop_section(
    page_size = page_size(orient = "portrait", width = 8.27, height = 11.69),
    page_margins = page_mar(
      top = 0.35, bottom = 0.35,
      left = 0.35, right = 0.35,
      header = 0.2, footer = 0.2
    )
  )
  
  # =========================================================
  # Save DOCX (titleなし)
  # =========================================================
  save_ft_docx <- function(ft, path, section_def) {
    doc <- read_docx() %>%
      body_set_default_section(section_def) %>%
      body_add_flextable(ft)
    
    print(doc, target = path)
    invisible(path)
  }
  
  save_ft_docx(ft1, out_docx1, sec)
  save_ft_docx(ft2, out_docx2, sec)
 
  } #word変換
# =========================================================
# Supplementary Tables 1 and 2 - final submission style
# - Table 2: combination drugs kept as combinations
# - Table 2: custom Drug class order
# - Save as DOCX
# =========================================================

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(officer)
library(flextable)

setwd("X:/R")

# =========================================================
# Paths
# =========================================================
file_path <- "Supplementary_Table_1_collapsed.xlsx"
out_dir   <- file.path("X:/R", "word_supp_tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_docx1 <- file.path(out_dir, "Supplementary_Table_1_ICD10_final.docx")
out_docx2 <- file.path(out_dir, "Supplementary_Table_2_Drugs_final.docx")

# =========================================================
# Load
# =========================================================
raw <- read_excel(file_path)

if (!all(c("term", "definition") %in% names(raw))) {
  raw <- raw %>% rename(term = 1, definition = 2)
}

dat <- raw %>%
  mutate(
    term       = str_squish(as.character(term)),
    definition = str_squish(as.character(definition))
  ) %>%
  filter(!is.na(term), !is.na(definition), term != "", definition != "")

# =========================================================
# Drug terms
# =========================================================
drug_terms <- c(
  "SGLT2 inhibitor",
  "Angiotensin II receptor blocker (ARB)",
  "Angiotensin receptor–neprilysin inhibitor (ARNI)",
  "Angiotensin-converting enzyme inhibitor (ACE inhibitor)",
  "SGLT2 inhibitor combination",
  "ARB + calcium channel blocker",
  "ARB + diuretic"
)

combo_terms <- c(
  "ARB + calcium channel blocker",
  "ARB + diuretic",
  "Angiotensin receptor–neprilysin inhibitor (ARNI)",
  "SGLT2 inhibitor combination"
)

drug_class_order <- c(
  "Angiotensin II receptor blocker (ARB)",
  "Angiotensin-converting enzyme inhibitor (ACE inhibitor)",
  "Angiotensin receptor–neprilysin inhibitor (ARNI)",
  "ARB + calcium channel blocker",
  "ARB + diuretic",
  "SGLT2 inhibitor",
  "SGLT2 inhibitor combination"
)

# =========================================================
# Build Supplementary Table 2
# 配合剤は「1つの合剤単位」で残す
# =========================================================
drug_dat <- dat %>%
  filter(term %in% drug_terms) %>%
  select(term, definition)

tbl2 <- drug_dat %>%
  mutate(
    definition = str_squish(definition),
    definition = str_replace_all(definition, "；", ";"),
    definition = str_replace_all(definition, "\\s*;\\s*", "; ")
  ) %>%
  rowwise() %>%
  mutate(
    generic_clean = if (term %in% combo_terms) {
      parts <- unlist(str_split(definition, "\\s*;\\s*"))
      parts <- str_squish(parts)
      parts <- parts[parts != ""]
      parts <- str_replace_all(parts, "\\s*/\\s*", " + ")
      parts <- str_replace_all(parts, "\\s*,\\s*", " + ")
      paste(unique(parts), collapse = "; ")
    } else {
      parts <- unlist(str_split(definition, "\\s*[,/]\\s*"))
      parts <- str_squish(parts)
      parts <- parts[parts != ""]
      paste(sort(unique(parts)), collapse = "; ")
    }
  ) %>%
  ungroup() %>%
  transmute(
    `Drug class`   = term,
    `Generic name` = generic_clean
  ) %>%
  distinct() %>%
  group_by(`Drug class`) %>%
  summarise(
    `Generic name` = paste(unique(`Generic name`), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    `Drug class` = str_replace_all(`Drug class`, "\u2011", "-")
  ) %>%
  mutate(
    `Drug class` = factor(`Drug class`, levels = drug_class_order)
  ) %>%
  arrange(`Drug class`) %>%
  mutate(
    `Drug class` = as.character(`Drug class`),
    `Drug class` = str_replace_all(`Drug class`, "-", "\u2011")
  )

# =========================================================
# Build Supplementary Table 1
# =========================================================
icd_dat <- dat %>%
  filter(!term %in% drug_terms) %>%
  select(term, definition)

tbl1 <- icd_dat %>%
  group_by(`Disease name` = term) %>%
  summarise(
    `ICD-10 code` = paste(unique(definition), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    `ICD-10 code` = str_replace_all(`ICD-10 code`, "\\s*;\\s*", "; ")
  ) %>%
  arrange(`Disease name`)

# =========================================================
# Flextable helper for Supplementary Table 2
# =========================================================
make_ft_supp_drug_final <- function(df) {
  ft <- flextable(df)
  
  ft <- ft %>%
    bold(part = "header") %>%
    align(align = "left", part = "all") %>%
    valign(valign = "top", part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 10, part = "body") %>%
    fontsize(size = 10.5, part = "header") %>%
    line_spacing(space = 1.0, part = "all") %>%
    border_remove() %>%
    hline_top(border = fp_border(width = 1.0)) %>%
    hline(part = "header", border = fp_border(width = 0.8)) %>%
    hline(i = seq_len(nrow(df)), part = "body", border = fp_border(width = 0.35)) %>%
    hline_bottom(border = fp_border(width = 1.0)) %>%
    set_table_properties(layout = "fixed", width = 1)
  
  ft <- ft %>%
    width(j = "Drug class",   width = 3.3) %>%
    width(j = "Generic name", width = 7.7)
  
  ft <- ft %>%
    padding(
      part = "all",
      padding.top = 1.5, padding.bottom = 1.5,
      padding.left = 2.5, padding.right = 2.5
    ) %>%
    padding(j = "Drug class", part = "all", padding.left = 2.2, padding.right = 2.0) %>%
    padding(j = "Generic name", part = "all", padding.left = 0.8, padding.right = 1.5)
  
  ft <- ft %>%
    bg(part = "header", bg = "white") %>%
    height_all(height = 0.24)
  
  ft
}

# =========================================================
# Flextable helper for Supplementary Table 1
# =========================================================
make_ft_supp_icd_final <- function(df) {
  ft <- flextable(df)
  
  ft <- ft %>%
    bold(part = "header") %>%
    align(align = "left", part = "all") %>%
    valign(valign = "top", part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 10, part = "body") %>%
    fontsize(size = 10.5, part = "header") %>%
    line_spacing(space = 1.0, part = "all") %>%
    border_remove() %>%
    hline_top(border = fp_border(width = 1.0)) %>%
    hline(part = "header", border = fp_border(width = 0.8)) %>%
    hline(i = seq_len(nrow(df)), part = "body", border = fp_border(width = 0.35)) %>%
    hline_bottom(border = fp_border(width = 1.0)) %>%
    set_table_properties(layout = "fixed", width = 1)
  
  ft <- ft %>%
    width(j = "Disease name", width = 3.8) %>%
    width(j = "ICD-10 code",  width = 7.2)
  
  ft <- ft %>%
    padding(
      part = "all",
      padding.top = 1.5, padding.bottom = 1.5,
      padding.left = 1.5, padding.right = 1.5
    ) %>%
    padding(j = "Disease name", part = "all", padding.left = 1.2, padding.right = 1.2) %>%
    padding(j = "ICD-10 code",  part = "all", padding.left = 0.6, padding.right = 1.0)
  
  ft <- ft %>%
    bg(part = "header", bg = "white") %>%
    height_all(height = 0.24)
  
  ft
}

# =========================================================
# Create flextables
# =========================================================
ft1 <- make_ft_supp_icd_final(tbl1)
ft2 <- make_ft_supp_drug_final(tbl2)

# =========================================================
# A4 portrait section
# =========================================================
sec <- prop_section(
  page_size = page_size(orient = "portrait", width = 8.27, height = 11.69),
  page_margins = page_mar(
    top = 0.35, bottom = 0.35,
    left = 0.35, right = 0.35,
    header = 0.2, footer = 0.2
  )
)

# =========================================================
# Save DOCX (titleなし)
# =========================================================
save_ft_docx <- function(ft, path, section_def) {
  doc <- read_docx() %>%
    body_set_default_section(section_def) %>%
    body_add_flextable(ft)
  
  print(doc, target = path)
  invisible(path)
}

save_ft_docx(ft1, out_docx1, sec)
save_ft_docx(ft2, out_docx2, sec)

# =========================================================
# Message
# =========================================================
message("Saved DOCX:")
message(out_docx1)
message(out_docx2)
