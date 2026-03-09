############################################################
# Figure 4. Subgroup Forest Plot (Sex / Age / CKD)
# - Overall is recalculated using the SAME population/model as Table 2
# - Figure4-3 is subgroup-only version (Overall removed)
# - RR / RD included
############################################################

graphics.off()

# ==========================
# Packages
# ==========================
pkgs <- c("readr","dplyr","survival","broom","tibble",
          "forestploter","grid","gridExtra","ragg")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)

library(readr)
library(dplyr)
library(survival)
library(broom)
library(tibble)
library(forestploter)
library(grid)
library(gridExtra)
library(ragg)

# ==========================
# Paths
# ==========================
setwd("X:/R")
in_csv <- "jin1_Eligibile.csv"
outdir <- file.path(getwd(), "figure_table")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

out_pdf_4   <- file.path(outdir, "Figure4_ForestPlot_Subgroups_Sex_Age_CKD.pdf")
out_tiff_4  <- file.path(outdir, "Figure4_ForestPlot_Subgroups_Sex_Age_CKD.tiff")
out_pdf_43  <- file.path(outdir, "Figure4-3_ForestPlot_Subgroups_Sex_Age_CKD.pdf")
out_tiff_43 <- file.path(outdir, "Figure4-3_ForestPlot_Subgroups_Sex_Age_CKD.tiff")

# ==========================
# Font
# ==========================
base_family <- {
  f <- c("Yu Gothic", "MS Gothic", "Meiryo", "Arial Unicode MS", "Arial")
  ok <- f[f %in% names(grDevices::windowsFonts())]
  if (length(ok) == 0) "sans" else ok[1]
}

# ==========================
# Helpers
# ==========================
disp_group <- function(x){
  dplyr::recode(
    as.character(x),
    "nonAKD"       = "non-AKD",
    "Recovery"     = "AKD with recovery",
    "Non-Recovery" = "AKD without recovery"
  )
}

disp_sex <- function(x){
  xx <- as.character(x)
  dplyr::case_when(
    xx %in% c("F","Female","0","2") ~ "Female",
    xx %in% c("M","Male","1")       ~ "Male",
    TRUE ~ xx
  )
}

disp_agegrp <- function(x){
  dplyr::recode(as.character(x),
                "lt75" = "<75",
                "ge75" = "≥75")
}

disp_ckdgrp <- function(x){
  dplyr::recode(as.character(x),
                "0" = "CKD no",
                "1" = "CKD yes")
}

save_pdf_onepage_safe <- function(grob, filename,
                                  page_w_mm = 360,
                                  page_h_mm = 190,
                                  margin_left_mm = 18,
                                  margin_right_mm = 22,
                                  margin_top_mm = 12,
                                  margin_bottom_mm = 40,
                                  scale = 0.68,
                                  family = "sans"){
  pdf(filename,
      width  = page_w_mm/25.4,
      height = page_h_mm/25.4,
      family = family)
  
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    x = grid::unit(margin_left_mm, "mm"),
    y = grid::unit(margin_bottom_mm, "mm"),
    just = c("left","bottom"),
    width  = grid::unit(1, "npc") - grid::unit(margin_left_mm + margin_right_mm, "mm"),
    height = grid::unit(1, "npc") - grid::unit(margin_top_mm + margin_bottom_mm, "mm"),
    clip = "on"
  ))
  grid::pushViewport(grid::viewport(
    x = 0.5, y = 0.5,
    width  = grid::unit(scale, "npc"),
    height = grid::unit(scale, "npc"),
    just = c("center","center"),
    clip = "on"
  ))
  grid::grid.draw(grob)
  grid::popViewport(2)
  dev.off()
  cat("Saved PDF:", normalizePath(filename), "\n")
}

save_tiff_onepage_safe <- function(grob, filename,
                                   page_w_mm = 360,
                                   page_h_mm = 190,
                                   margin_left_mm = 18,
                                   margin_right_mm = 22,
                                   margin_top_mm = 12,
                                   margin_bottom_mm = 40,
                                   scale = 0.68,
                                   dpi = 600,
                                   family = "sans"){
  ragg::agg_tiff(filename,
                 width  = page_w_mm/25.4,
                 height = page_h_mm/25.4,
                 units = "in",
                 res = dpi,
                 compression = "lzw")
  
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    x = grid::unit(margin_left_mm, "mm"),
    y = grid::unit(margin_bottom_mm, "mm"),
    just = c("left","bottom"),
    width  = grid::unit(1, "npc") - grid::unit(margin_left_mm + margin_right_mm, "mm"),
    height = grid::unit(1, "npc") - grid::unit(margin_top_mm + margin_bottom_mm, "mm"),
    clip = "on"
  ))
  grid::pushViewport(grid::viewport(
    x = 0.5, y = 0.5,
    width  = grid::unit(scale, "npc"),
    height = grid::unit(scale, "npc"),
    just = c("center","center"),
    clip = "on"
  ))
  grid::grid.draw(grob)
  grid::popViewport(2)
  dev.off()
  cat("Saved TIFF:", normalizePath(filename), "\n")
}

# ==========================
# 0) Load
# ==========================
loc <- locale(encoding = "SHIFT-JIS")
jin1_Eligibile <- read_csv(in_csv, locale = loc)

# ==========================
# 1) Build base dataset
# ==========================
dat_base <- jin1_Eligibile %>%
  group_by(id) %>%
  arrange(index_date, date, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  filter(exclude == "include", jin_status %in% c("nonAKD","AKD")) %>%
  mutate(
    group = case_when(
      jin_status == "nonAKD" ~ "nonAKD",
      jin_status == "AKD" & `150_210recovery` == 1 ~ "Recovery",
      jin_status == "AKD" & `150_210recovery` == 2 ~ "Non-Recovery",
      jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 1 ~ "Recovery",
      jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` %in% c(0,2) ~ "Non-Recovery",
      TRUE ~ NA_character_
    ),
    group = factor(group, levels = c("nonAKD","Recovery","Non-Recovery")),
    arb_acei_use = if_else(coalesce(arb,0)==1 | coalesce(acei,0)==1, 1L, 0L),
    time_years   = as.numeric(last_follow_death - index_plus_210) / 365.25
  ) %>%
  filter(!is.na(group),
         !is.na(time_years), time_years >= 0,
         !is.na(primary_death),
         !is.na(age),
         !is.na(index_cre))

# ==========================
# 2) Overall = same as Table 2
# ==========================
dat_overall <- dat_base %>%
  mutate(group = factor(group, levels = c("nonAKD","Recovery","Non-Recovery")))

fit_overall <- coxph(
  Surv(time_years, primary_death) ~
    group + age + index_cre + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
  data = dat_overall
)

# ==========================
# 3) Subgroup datasets
# ==========================
dat_sub <- dat_base %>%
  mutate(
    sex_chr = disp_sex(sex),
    sex_f = case_when(
      sex_chr %in% c("Female","F") ~ "F",
      sex_chr %in% c("Male","M")   ~ "M",
      TRUE ~ NA_character_
    ),
    sex_f = factor(sex_f, levels = c("F","M")),
    age75 = if_else(age >= 75, 1L, 0L),
    age75_f = factor(age75, levels = c(0,1), labels = c("lt75","ge75")),
    CKD_status2 = case_when(
      as.character(CKD_status) %in% c("CKD","1")    ~ "CKD",
      as.character(CKD_status) %in% c("nonCKD","0") ~ "nonCKD",
      TRUE ~ NA_character_
    ),
    CKD_status2 = factor(CKD_status2, levels = c("nonCKD","CKD")),
    ckd_bin = case_when(
      CKD_status2 == "CKD"    ~ 1L,
      CKD_status2 == "nonCKD" ~ 0L,
      TRUE ~ NA_integer_
    )
  )

dat_sex <- dat_sub %>% filter(!is.na(sex_f))
dat_age <- dat_sub %>% filter(!is.na(age75_f))
dat_ckd <- dat_sub %>% filter(!is.na(ckd_bin))

# ==========================
# 4) Interaction models
# ==========================
fit_sex <- coxph(
  Surv(time_years, primary_death) ~
    group * sex_f + age + index_cre + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
  data = dat_sex
)

fit_age <- coxph(
  Surv(time_years, primary_death) ~
    group * age75_f + age + index_cre + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
  data = dat_age
)

fit_ckd <- coxph(
  Surv(time_years, primary_death) ~
    group * ckd_bin + age + index_cre + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
  data = dat_ckd
)

# ==========================
# 5) Summary/stat helpers
# ==========================
calc_stats_overall <- function(data){
  data %>%
    group_by(group) %>%
    summarise(
      Number = n(),
      Deaths = sum(primary_death),
      Death_pct = sprintf("%.1f", 100 * Deaths / Number),
      risk = Deaths / Number,
      .groups = "drop"
    ) %>%
    mutate(Subgroup = "Overall",
           Group = as.character(group)) %>%
    select(Subgroup, Group, Number, Deaths, Death_pct, risk)
}

calc_stats_sub <- function(data, subgroup_var, subgroup_lab_fun = NULL){
  out <- data %>%
    group_by(!!sym(subgroup_var), group) %>%
    summarise(
      Number = n(),
      Deaths = sum(primary_death),
      Death_pct = sprintf("%.1f", 100 * Deaths / Number),
      risk = Deaths / Number,
      .groups = "drop"
    ) %>%
    mutate(
      Subgroup = as.character(!!sym(subgroup_var)),
      Group = as.character(group)
    )
  if (!is.null(subgroup_lab_fun)) out$Subgroup <- subgroup_lab_fun(out$Subgroup)
  out %>% select(Subgroup, Group, Number, Deaths, Death_pct, risk)
}

extract_hr_overall <- function(fit){
  td <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  hr <- td %>%
    filter(term %in% c("groupRecovery", "groupNon-Recovery")) %>%
    mutate(
      Group = case_when(
        term == "groupRecovery"     ~ "Recovery",
        term == "groupNon-Recovery" ~ "Non-Recovery",
        TRUE ~ NA_character_
      ),
      HR = estimate,
      CI_low = conf.low,
      CI_high = conf.high,
      P_value = p.value
    ) %>%
    select(Group, HR, CI_low, CI_high, P_value)
  
  bind_rows(
    tibble(Group = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
    hr
  ) %>%
    mutate(Subgroup = "Overall") %>%
    select(Subgroup, everything())
}

hr_from_interaction <- function(fit, subgroup_levels,
                                b_rec = "groupRecovery",
                                b_nrec = "groupNon-Recovery",
                                int_suffix_other){
  
  cf <- coef(fit)
  vc <- vcov(fit)
  nm <- names(cf)
  
  find_name <- function(patterns){
    hit <- nm[Reduce(`|`, lapply(patterns, function(p) grepl(p, nm)))]
    if (length(hit) == 0) return(NA_character_)
    hit[1]
  }
  
  int_rec <- find_name(c(
    paste0("^", b_rec, ":", int_suffix_other, "$"),
    paste0("^", int_suffix_other, ":", b_rec, "$")
  ))
  int_nrec <- find_name(c(
    paste0("^", b_nrec, ":", int_suffix_other, "$"),
    paste0("^", int_suffix_other, ":", b_nrec, "$")
  ))
  
  comp_one <- function(b, int, is_other){
    if (!is_other){
      est <- cf[b]
      v   <- vc[b, b]
    } else {
      est <- cf[b] + cf[int]
      v   <- vc[b,b] + vc[int,int] + 2 * vc[b,int]
    }
    se <- sqrt(v)
    HR <- exp(est)
    lo <- exp(est - 1.96 * se)
    hi <- exp(est + 1.96 * se)
    z  <- est / se
    p  <- 2 * pnorm(-abs(z))
    c(HR = HR, CI_low = lo, CI_high = hi, P_value = p)
  }
  
  out <- rbind(
    c(subgroup_levels[1], "Recovery",     comp_one(b_rec,  int_rec,  FALSE)),
    c(subgroup_levels[1], "Non-Recovery", comp_one(b_nrec, int_nrec, FALSE)),
    c(subgroup_levels[2], "Recovery",     comp_one(b_rec,  int_rec,  TRUE)),
    c(subgroup_levels[2], "Non-Recovery", comp_one(b_nrec, int_nrec, TRUE))
  )
  
  tibble(
    Subgroup = out[,1],
    Group    = out[,2],
    HR       = as.numeric(out[,3]),
    CI_low   = as.numeric(out[,4]),
    CI_high  = as.numeric(out[,5]),
    P_value  = as.numeric(out[,6])
  )
}

# ==========================
# 6) Stats + HR
# ==========================
stats_overall <- calc_stats_overall(dat_overall)
stats_sex <- calc_stats_sub(dat_sex, "sex_f", function(x) dplyr::recode(x, "F" = "Female", "M" = "Male"))
stats_age <- calc_stats_sub(dat_age, "age75_f", disp_agegrp)
stats_ckd <- calc_stats_sub(dat_ckd, "ckd_bin", disp_ckdgrp)

hr_overall <- extract_hr_overall(fit_overall)

hr_sex0 <- hr_from_interaction(fit_sex, subgroup_levels = c("Female","Male"), int_suffix_other = "sex_fM")
hr_sex <- bind_rows(
  tibble(Subgroup = "Female", Group = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
  filter(hr_sex0, Subgroup == "Female"),
  tibble(Subgroup = "Male", Group = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
  filter(hr_sex0, Subgroup == "Male")
)

hr_age0 <- hr_from_interaction(fit_age, subgroup_levels = c("<75","≥75"), int_suffix_other = "age75_fge75")
hr_age <- bind_rows(
  tibble(Subgroup = "<75", Group = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
  filter(hr_age0, Subgroup == "<75"),
  tibble(Subgroup = "≥75", Group = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
  filter(hr_age0, Subgroup == "≥75")
)

hr_ckd0 <- hr_from_interaction(fit_ckd, subgroup_levels = c("CKD no","CKD yes"), int_suffix_other = "ckd_bin")
hr_ckd <- bind_rows(
  tibble(Subgroup = "CKD no", Group = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
  filter(hr_ckd0, Subgroup == "CKD no"),
  tibble(Subgroup = "CKD yes", Group = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
  filter(hr_ckd0, Subgroup == "CKD yes")
)

# ==========================
# 7) Merge
# ==========================
build_block <- function(stats_df, hr_df){
  stats_df %>% left_join(hr_df, by = c("Subgroup","Group"))
}

forest_full <- bind_rows(
  build_block(stats_overall, hr_overall),
  build_block(stats_sex,     hr_sex),
  build_block(stats_age,     hr_age),
  build_block(stats_ckd,     hr_ckd)
) %>%
  mutate(
    Subgroup = factor(Subgroup,
                      levels = c("Overall","Female","Male","<75","≥75","CKD no","CKD yes")),
    Group = factor(Group, levels = c("nonAKD","Recovery","Non-Recovery"))
  ) %>%
  arrange(Subgroup, Group) %>%
  group_by(Subgroup) %>%
  mutate(
    risk_ref = risk[Group == "nonAKD"][1],
    RR = if_else(Group == "nonAKD", NA_real_, risk / risk_ref),
    RD = if_else(Group == "nonAKD", NA_real_, risk - risk_ref),
    row_in_block = row_number()
  ) %>%
  ungroup()

# ==========================
# 8) Display builders
# ==========================
make_display_df <- function(df, include_overall = TRUE){
  
  if (!include_overall) {
    df <- df %>% filter(Subgroup != "Overall")
  }
  
  df <- df %>%
    mutate(
      block_name = case_when(
        Subgroup == "Overall" ~ "Overall",
        Subgroup %in% c("Female","Male") ~ "Sex subgroup",
        Subgroup %in% c("<75","≥75") ~ "Age subgroup",
        Subgroup %in% c("CKD no","CKD yes") ~ "CKD subgroup",
        TRUE ~ ""
      ),
      subgroup_value = case_when(
        Subgroup == "Overall" ~ "",
        TRUE ~ as.character(Subgroup)
      ),
      Subgroup_disp = if_else(row_in_block == 1, block_name, ""),
      Group_disp = case_when(
        Subgroup == "Overall" & Group == "nonAKD"       ~ "  non-AKD",
        Subgroup == "Overall" & Group == "Recovery"     ~ "  AKD with recovery",
        Subgroup == "Overall" & Group == "Non-Recovery" ~ "  AKD without recovery",
        
        Group == "nonAKD"       ~ paste0("  ", subgroup_value, ": non-AKD"),
        Group == "Recovery"     ~ paste0("  ", subgroup_value, ": AKD with recovery"),
        Group == "Non-Recovery" ~ paste0("  ", subgroup_value, ": AKD without recovery"),
        TRUE ~ ""
      ),
      Number_chr = as.character(Number),
      Deaths_chr = paste0(Deaths, " (", Death_pct, ")"),
      Risk_chr   = sprintf("%.1f", 100 * risk),
      RR_chr = case_when(
        Group == "nonAKD" ~ "Reference",
        is.na(RR) ~ "",
        TRUE ~ sprintf("%.2f", RR)
      ),
      RD_chr = case_when(
        Group == "nonAKD" ~ "",
        is.na(RD) ~ "",
        TRUE ~ sprintf("%+.1f", 100 * RD)
      ),
      HR_chr = case_when(
        Group == "nonAKD" ~ "Reference",
        TRUE ~ sprintf("%.2f (%.2f–%.2f)", HR, CI_low, CI_high)
      ),
      P_chr = case_when(
        Group == "nonAKD" ~ "",
        is.na(P_value) ~ "",
        P_value < 0.001 ~ "<0.001",
        P_value < 0.01  ~ sprintf("%.3f", P_value),
        TRUE            ~ sprintf("%.2f", P_value)
      ),
      blank_ci = paste(rep(" ", 20), collapse = " ")
    )
  
  plot_df <- df %>%
    select(
      Subgroup = Subgroup_disp,
      Group = Group_disp,
      Number = Number_chr,
      `Deaths (%)` = Deaths_chr,
      `Risk (%)` = Risk_chr,
      `RR vs non-AKD` = RR_chr,
      `RD vs non-AKD (pp)` = RD_chr,
      ` ` = blank_ci, 
      `HR (95%CI)` = HR_chr,
      `P-value` = P_chr
    )
  
  list(full = df, plot = plot_df)
}

# ==========================
# 9) Figure builders
# ==========================
build_forest_plot <- function(df_full, df_plot){
  fp <- forestploter::forest(
    data  = df_plot,
    est   = df_full$HR,
    lower = df_full$CI_low,
    upper = df_full$CI_high,
    sizes = 0.6,
    ci_column = 8,
    ref_line  = 1,
    x_trans   = "log",
    xlim      = c(0.5, 12),
    ticks_at  = c(0.5, 1, 2, 4, 8),
    arrow_lab = c("", ""),
    xlab = "Higher risk                                 Lower risk"
  )
  
  fp <- edit_plot(
    fp,
    col = c(3,4,5,6,7,9,10),
    which = "text",
    hjust = unit(1, "npc"),
    x = unit(1, "npc")
  )
  
  fp <- add_border(fp, row = 0, where = "bottom", gp = gpar(lwd = 1))
  fp$heights <- rep(unit(6.4, "mm"), nrow(fp))
  fp
}

# ==========================
# 10) Figure 4
# ==========================
disp4 <- make_display_df(forest_full, include_overall = TRUE)
fp4 <- build_forest_plot(disp4$full, disp4$plot)

grid.newpage()
grid.draw(fp4)

save_pdf_onepage_safe(fp4, out_pdf_4, family = base_family)
save_tiff_onepage_safe(fp4, out_tiff_4, family = base_family)

# ==========================
# 11) Figure 4-3
# ==========================
disp43 <- make_display_df(forest_full, include_overall = FALSE)
fp43 <- build_forest_plot(disp43$full, disp43$plot)

grid.newpage()
grid.draw(fp43)

save_pdf_onepage_safe(fp43, out_pdf_43, family = base_family)
save_tiff_onepage_safe(fp43, out_tiff_43, family = base_family)


############################################################
# Figure 4-1 / Figure 4-2
# Based on original Figure 4 (Sex / Age / CKD)
#  - Figure4-1: keep Risk / RR / RD + add P for interaction
#  - Figure4-2: HR-focused version
############################################################

graphics.off()

# ==========================
# Packages
# ==========================
pkgs <- c("readr","dplyr","survival","broom","tibble",
          "forestploter","grid","gridExtra","ragg")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)

library(readr)
library(dplyr)
library(survival)
library(broom)
library(tibble)
library(forestploter)
library(grid)
library(gridExtra)
library(ragg)

# ==========================
# Paths
# ==========================
setwd("X:/R")
in_csv  <- "jin1_Eligibile.csv"
outdir  <- file.path(getwd(), "figure_table")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

out_pdf_41  <- file.path(outdir, "Figure4-1_ForestPlot_Subgroups_Sex_Age_CKD_with_Pinteraction.pdf")
out_tiff_41 <- file.path(outdir, "Figure4-1_ForestPlot_Subgroups_Sex_Age_CKD_with_Pinteraction.tiff")

out_pdf_42  <- file.path(outdir, "Figure4-2_ForestPlot_Subgroups_Sex_Age_CKD_HRonly_with_Pinteraction.pdf")
out_tiff_42 <- file.path(outdir, "Figure4-2_ForestPlot_Subgroups_Sex_Age_CKD_HRonly_with_Pinteraction.tiff")

# ==========================
# Font
# ==========================
base_family <- {
  f <- c("Yu Gothic", "MS Gothic", "Meiryo", "Arial Unicode MS", "Arial")
  ok <- f[f %in% names(grDevices::windowsFonts())]
  if (length(ok) == 0) "sans" else ok[1]
}

# ==========================
# Helpers
# ==========================
disp_group <- function(x){
  dplyr::recode(
    as.character(x),
    "nonAKD"       = "non-AKD",
    "Recovery"     = "AKD with recovery",
    "Non-Recovery" = "AKD without recovery"
  )
}

disp_sex <- function(x){
  xx <- as.character(x)
  dplyr::case_when(
    xx %in% c("F","Female","0","2") ~ "Female",
    xx %in% c("M","Male","1")       ~ "Male",
    TRUE ~ xx
  )
}

disp_ckd <- function(x){
  dplyr::case_when(
    as.character(x) %in% c("nonCKD","0") ~ "CKD no",
    as.character(x) %in% c("CKD","1")    ~ "CKD yes",
    TRUE ~ as.character(x)
  )
}

disp_agegrp <- function(x){
  dplyr::recode(as.character(x),
                "lt75" = "<75",
                "ge75" = "≥75")
}

fmt_p <- function(p){
  dplyr::case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "<0.001",
    p < 0.01  ~ sprintf("%.3f", p),
    TRUE      ~ sprintf("%.2f", p)
  )
}

get_lrt_p <- function(anova_obj){
  cn <- colnames(anova_obj)
  cand <- c("P(>|Chi|)", "Pr(>|Chi|)", "P(>|Chisq|)", "Pr(>|Chisq|)")
  hit <- cand[cand %in% cn]
  
  if (length(hit) > 0) {
    return(as.numeric(anova_obj[2, hit[1], drop = TRUE]))
  }
  
  suppressWarnings({
    p_last <- as.numeric(anova_obj[2, ncol(anova_obj), drop = TRUE])
  })
  
  if (length(p_last) == 1 && is.finite(p_last)) return(p_last)
  return(NA_real_)
}

# ==========================
# Safe save helpers
# ==========================
save_pdf_onepage_safe <- function(grob, filename,
                                  page_w_mm = 360,
                                  page_h_mm = 190,
                                  margin_left_mm = 18,
                                  margin_right_mm = 22,
                                  margin_top_mm = 12,
                                  margin_bottom_mm = 34,
                                  scale = 0.70,
                                  family = "sans"){
  
  pdf(filename,
      width = page_w_mm/25.4,
      height = page_h_mm/25.4,
      family = family,
      useDingbats = FALSE)
  
  grid::grid.newpage()
  
  grid::pushViewport(grid::viewport(
    x = unit(margin_left_mm, "mm"),
    y = unit(margin_bottom_mm, "mm"),
    just = c("left","bottom"),
    width  = unit(1, "npc") - unit(margin_left_mm + margin_right_mm, "mm"),
    height = unit(1, "npc") - unit(margin_top_mm + margin_bottom_mm, "mm"),
    clip = "on"
  ))
  
  grid::pushViewport(grid::viewport(
    x = 0.5, y = 0.5,
    width  = unit(scale, "npc"),
    height = unit(scale, "npc"),
    just = c("center","center"),
    clip = "on"
  ))
  
  grid::grid.draw(grob)
  grid::popViewport(2)
  dev.off()
  
  cat("Saved PDF :", normalizePath(filename), "\n")
}

save_tiff_onepage_safe <- function(grob, filename,
                                   page_w_mm = 360,
                                   page_h_mm = 190,
                                   margin_left_mm = 18,
                                   margin_right_mm = 22,
                                   margin_top_mm = 12,
                                   margin_bottom_mm = 34,
                                   scale = 0.70,
                                   dpi = 600,
                                   family = "sans"){
  
  ragg::agg_tiff(filename,
                 width  = page_w_mm/25.4,
                 height = page_h_mm/25.4,
                 units = "in",
                 res = dpi,
                 compression = "lzw")
  
  grid::grid.newpage()
  
  grid::pushViewport(grid::viewport(
    x = unit(margin_left_mm, "mm"),
    y = unit(margin_bottom_mm, "mm"),
    just = c("left","bottom"),
    width  = unit(1, "npc") - unit(margin_left_mm + margin_right_mm, "mm"),
    height = unit(1, "npc") - unit(margin_top_mm + margin_bottom_mm, "mm"),
    clip = "on"
  ))
  
  grid::pushViewport(grid::viewport(
    x = 0.5, y = 0.5,
    width  = unit(scale, "npc"),
    height = unit(scale, "npc"),
    just = c("center","center"),
    clip = "on"
  ))
  
  grid::grid.draw(grob)
  grid::popViewport(2)
  dev.off()
  
  cat("Saved TIFF:", normalizePath(filename), "\n")
}

# ==========================
# Load data
# ==========================
loc <- locale(encoding = "SHIFT-JIS")
jin1_Eligibile <- read_csv(in_csv, locale = loc)

# ==========================
# Build analysis dataset
# ==========================
dat0 <- jin1_Eligibile %>%
  group_by(id) %>%
  arrange(index_date, date, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  filter(exclude == "include", jin_status %in% c("nonAKD","AKD")) %>%
  mutate(
    jin_label = case_when(
      jin_status == "nonAKD" ~ "nonAKD",
      jin_status == "AKD" & `150_210recovery` == 1 ~ "Recovery",
      jin_status == "AKD" & `150_210recovery` == 2 ~ "Non-Recovery",
      jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 1 ~ "Recovery",
      jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` %in% c(0,2) ~ "Non-Recovery",
      TRUE ~ NA_character_
    ),
    jin_label = factor(jin_label, levels = c("nonAKD","Recovery","Non-Recovery")),
    arb_acei_use = if_else(coalesce(arb, 0) == 1 | coalesce(acei, 0) == 1, 1L, 0L),
    time_years   = as.numeric(last_follow_death - index_plus_210) / 365.25
  ) %>%
  filter(!is.na(jin_label),
         !is.na(time_years), time_years >= 0,
         !is.na(primary_death),
         !is.na(age),
         !is.na(index_cre))

dat0 <- dat0 %>%
  mutate(
    sex_chr = disp_sex(sex),
    sex_f = case_when(
      sex_chr %in% c("Female","F") ~ "F",
      sex_chr %in% c("Male","M")   ~ "M",
      TRUE ~ NA_character_
    ),
    sex_f = factor(sex_f, levels = c("F","M"))
  ) %>%
  filter(!is.na(sex_f))

dat0 <- dat0 %>%
  mutate(
    age75 = if_else(age >= 75, 1L, 0L),
    age75_f = factor(age75, levels = c(0,1), labels = c("lt75","ge75"))
  ) %>%
  filter(!is.na(age75_f))

dat0 <- dat0 %>%
  mutate(
    CKD_status = case_when(
      as.character(CKD_status) %in% c("CKD","1")    ~ "CKD",
      as.character(CKD_status) %in% c("nonCKD","0") ~ "nonCKD",
      TRUE ~ as.character(CKD_status)
    ),
    CKD_status = factor(CKD_status, levels = c("nonCKD","CKD")),
    ckd_bin = case_when(
      CKD_status == "CKD"    ~ 1L,
      CKD_status == "nonCKD" ~ 0L,
      TRUE ~ NA_integer_
    )
  ) %>%
  filter(!is.na(CKD_status), !is.na(ckd_bin))

dat <- dat0

# ==========================
# Models
# ==========================
cox_core_cov <- as.formula(
  Surv(time_years, primary_death) ~
    jin_label + age + index_cre + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15
)

cox_overall <- coxph(
  update(cox_core_cov, . ~ . + CKD_status + sex_f),
  data = dat
)

cox_base_sex <- coxph(update(cox_core_cov, . ~ . + sex_f), data = dat)
cox_base_age <- coxph(update(cox_core_cov, . ~ . + age75_f), data = dat)
cox_base_ckd <- coxph(update(cox_core_cov, . ~ . + ckd_bin), data = dat)

cox_int_sex <- coxph(update(cox_core_cov, . ~ jin_label * sex_f + age + index_cre + arb_acei_use +
                              dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15),
                     data = dat)

cox_int_age <- coxph(update(cox_core_cov, . ~ jin_label * age75_f + age + index_cre + arb_acei_use +
                              dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15),
                     data = dat)

cox_int_ckd <- coxph(update(cox_core_cov, . ~ jin_label * ckd_bin + age + index_cre + arb_acei_use +
                              dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15),
                     data = dat)

# P for interaction (LRT)
p_int_sex <- get_lrt_p(anova(cox_base_sex, cox_int_sex, test = "LRT"))
p_int_age <- get_lrt_p(anova(cox_base_age, cox_int_age, test = "LRT"))
p_int_ckd <- get_lrt_p(anova(cox_base_ckd, cox_int_ckd, test = "LRT"))

p_int_sex_txt <- fmt_p(p_int_sex)
p_int_age_txt <- fmt_p(p_int_age)
p_int_ckd_txt <- fmt_p(p_int_ckd)

# ==========================
# Summary stats
# ==========================
calc_stats <- function(data, subgroup_var = NULL, subgroup_lab_fun = NULL){
  if (is.null(subgroup_var)) {
    out <- data %>%
      group_by(jin_label) %>%
      summarise(
        Number = n(),
        Deaths = sum(primary_death),
        Death_pct = sprintf("%.1f", 100 * Deaths / Number),
        risk = Deaths / Number,
        .groups = "drop"
      ) %>%
      mutate(Subgroup = "Overall", Group = as.character(jin_label))
  } else {
    out <- data %>%
      group_by(!!sym(subgroup_var), jin_label) %>%
      summarise(
        Number = n(),
        Deaths = sum(primary_death),
        Death_pct = sprintf("%.1f", 100 * Deaths / Number),
        risk = Deaths / Number,
        .groups = "drop"
      ) %>%
      mutate(
        Subgroup = as.character(!!sym(subgroup_var)),
        Group = as.character(jin_label)
      )
    if (!is.null(subgroup_lab_fun)) out$Subgroup <- subgroup_lab_fun(out$Subgroup)
  }
  out %>% dplyr::select(Subgroup, Group, Number, Deaths, Death_pct, risk)
}

stats_overall <- calc_stats(dat)
stats_sex <- calc_stats(dat, "sex_f", function(x) dplyr::recode(x, "F"="Female", "M"="Male"))
stats_age <- calc_stats(dat, "age75_f", disp_agegrp)
stats_ckd <- calc_stats(dat, "ckd_bin", function(x) dplyr::recode(x, "0"="CKD no", "1"="CKD yes"))

# ==========================
# HR extraction
# ==========================
extract_hr_overall <- function(fit, ref_level = "nonAKD"){
  td <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  hr <- td %>%
    filter(grepl("^jin_label", term)) %>%
    mutate(
      Group = case_when(
        term == "jin_labelRecovery" ~ "Recovery",
        term == "jin_labelNon-Recovery" ~ "Non-Recovery",
        TRUE ~ NA_character_
      ),
      HR = estimate, CI_low = conf.low, CI_high = conf.high, P_value = p.value
    ) %>%
    dplyr::select(Group, HR, CI_low, CI_high, P_value)
  
  bind_rows(
    tibble(Group = ref_level, HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
    hr
  )
}

hr_overall <- extract_hr_overall(cox_overall) %>% mutate(Subgroup = "Overall")

hr_from_interaction <- function(fit,
                                subgroup_levels = c("ref","other"),
                                b_rec = "jin_labelRecovery",
                                b_nrec = "jin_labelNon-Recovery",
                                int_suffix_other){
  cf <- coef(fit); vc <- vcov(fit); nm <- names(cf)
  
  find_name <- function(patterns){
    hit <- nm[Reduce(`|`, lapply(patterns, function(p) grepl(p, nm)))]
    if (length(hit) == 0) return(NA_character_)
    hit[1]
  }
  
  int_rec  <- find_name(c(paste0("^", b_rec, ":", int_suffix_other, "$"),
                          paste0("^", int_suffix_other, ":", b_rec, "$")))
  int_nrec <- find_name(c(paste0("^", b_nrec, ":", int_suffix_other, "$"),
                          paste0("^", int_suffix_other, ":", b_nrec, "$")))
  
  comp_one <- function(b, int, is_other){
    if (!is_other){
      est <- cf[b]; v <- vc[b,b]
    } else {
      est <- cf[b] + cf[int]
      v   <- vc[b,b] + vc[int,int] + 2*vc[b,int]
    }
    se <- sqrt(v)
    HR <- exp(est)
    lo <- exp(est - 1.96*se)
    hi <- exp(est + 1.96*se)
    z  <- est / se
    p_value <- 2 * pnorm(-abs(z))
    c(HR = HR, CI_low = lo, CI_high = hi, P_value = p_value)
  }
  
  out <- rbind(
    c(subgroup_levels[1], "Recovery",     comp_one(b_rec,  int_rec,  FALSE)),
    c(subgroup_levels[1], "Non-Recovery", comp_one(b_nrec, int_nrec, FALSE)),
    c(subgroup_levels[2], "Recovery",     comp_one(b_rec,  int_rec,  TRUE)),
    c(subgroup_levels[2], "Non-Recovery", comp_one(b_nrec, int_nrec, TRUE))
  )
  
  tibble(
    Subgroup = out[,1],
    Group    = out[,2],
    HR       = as.numeric(out[,3]),
    CI_low   = as.numeric(out[,4]),
    CI_high  = as.numeric(out[,5]),
    P_value  = as.numeric(out[,6])
  )
}

hr_sex0 <- hr_from_interaction(
  cox_int_sex,
  subgroup_levels = c("Female","Male"),
  int_suffix_other = "sex_fM"
)
hr_sex <- bind_rows(
  tibble(Subgroup="Female", Group="nonAKD", HR=1, CI_low=1, CI_high=1, P_value=NA_real_),
  filter(hr_sex0, Subgroup=="Female"),
  tibble(Subgroup="Male",   Group="nonAKD", HR=1, CI_low=1, CI_high=1, P_value=NA_real_),
  filter(hr_sex0, Subgroup=="Male")
) %>% distinct()

hr_age0 <- hr_from_interaction(
  cox_int_age,
  subgroup_levels = c("<75","≥75"),
  int_suffix_other = "age75_fge75"
)
hr_age <- bind_rows(
  tibble(Subgroup="<75", Group="nonAKD", HR=1, CI_low=1, CI_high=1, P_value=NA_real_),
  filter(hr_age0, Subgroup=="<75"),
  tibble(Subgroup="≥75", Group="nonAKD", HR=1, CI_low=1, CI_high=1, P_value=NA_real_),
  filter(hr_age0, Subgroup=="≥75")
) %>% distinct()

hr_ckd0 <- hr_from_interaction(
  cox_int_ckd,
  subgroup_levels = c("CKD no","CKD yes"),
  int_suffix_other = "ckd_bin"
)
hr_ckd <- bind_rows(
  tibble(Subgroup="CKD no", Group="nonAKD", HR=1, CI_low=1, CI_high=1, P_value=NA_real_),
  filter(hr_ckd0, Subgroup=="CKD no"),
  tibble(Subgroup="CKD yes", Group="nonAKD", HR=1, CI_low=1, CI_high=1, P_value=NA_real_),
  filter(hr_ckd0, Subgroup=="CKD yes")
) %>% distinct()

# ==========================
# Merge
# ==========================
build_block <- function(stats_df, hr_df){
  stats_df %>% left_join(hr_df, by = c("Subgroup","Group"))
}

forest_full <- bind_rows(
  build_block(stats_overall, hr_overall),
  build_block(stats_sex,     hr_sex),
  build_block(stats_age,     hr_age),
  build_block(stats_ckd,     hr_ckd)
) %>%
  mutate(
    Subgroup = factor(Subgroup,
                      levels = c("Overall",
                                 "Female","Male",
                                 "<75","≥75",
                                 "CKD no","CKD yes")),
    Group = factor(Group, levels = c("nonAKD","Recovery","Non-Recovery"))
  ) %>%
  arrange(Subgroup, Group) %>%
  group_by(Subgroup) %>%
  mutate(
    risk_ref = risk[Group=="nonAKD"][1],
    RR = if_else(Group=="nonAKD", NA_real_, risk / risk_ref),
    RD = if_else(Group=="nonAKD", NA_real_, risk - risk_ref)
  ) %>%
  ungroup()

# ==========================
# Build display table
# ==========================
forest_full <- forest_full %>%
  mutate(
    Subgroup_disp = case_when(
      Subgroup == "Overall" & Group=="nonAKD" ~ "Overall",
      
      Subgroup == "Female" & Group=="nonAKD"  ~ "Sex subgroup",
      Subgroup == "Female" & Group=="Recovery" ~ "  Female",
      Subgroup == "Male"   & Group=="nonAKD"  ~ "",
      Subgroup == "Male"   & Group=="Recovery" ~ "  Male",
      
      Subgroup == "<75" & Group=="nonAKD"     ~ "Age subgroup",
      Subgroup == "<75" & Group=="Recovery"   ~ "  <75",
      Subgroup == "≥75" & Group=="nonAKD"     ~ "",
      Subgroup == "≥75" & Group=="Recovery"   ~ "  ≥75",
      
      Subgroup == "CKD no" & Group=="nonAKD"  ~ "CKD subgroup",
      Subgroup == "CKD no" & Group=="Recovery"~ "  CKD no",
      Subgroup == "CKD yes"& Group=="nonAKD"  ~ "",
      Subgroup == "CKD yes"& Group=="Recovery"~ "  CKD yes",
      
      TRUE ~ ""
    ),
    Group_disp = factor(
      disp_group(Group),
      levels = c("non-AKD","AKD with recovery","AKD without recovery")
    ),
    `Number`     = as.character(Number),
    `Deaths (%)` = paste0(Deaths, " (", Death_pct, ")"),
    `Risk (%)`   = sprintf("%.1f", 100 * risk),
    `RR vs non-AKD` = case_when(
      Group=="nonAKD" ~ "Reference",
      is.na(RR) ~ "",
      TRUE ~ sprintf("%.2f", RR)
    ),
    `RD vs non-AKD (pp)` = case_when(
      Group=="nonAKD" ~ "",
      is.na(RD) ~ "",
      TRUE ~ sprintf("%+.1f", 100 * RD)
    ),
    `HR (95% CI)` = if_else(
      Group=="nonAKD",
      "Reference",
      sprintf("%.2f (%.2f–%.2f)", HR, CI_low, CI_high)
    ),
    `P-value` = case_when(
      Group=="nonAKD" ~ "",
      is.na(P_value)  ~ "",
      P_value < 0.001 ~ "<0.001",
      P_value < 0.01  ~ sprintf("%.3f", P_value),
      TRUE            ~ sprintf("%.2f", P_value)
    ),
    `P for interaction` = case_when(
      Subgroup == "Female" & Group == "Recovery" ~ p_int_sex_txt,
      Subgroup == "<75"    & Group == "Recovery" ~ p_int_age_txt,
      Subgroup == "CKD no" & Group == "Recovery" ~ p_int_ckd_txt,
      TRUE ~ ""
    )
  )

# ==========================
# Figure4-1 (keep RR/RD)
# ==========================
forest_full_41 <- forest_full
forest_full_41$` ` <- paste(rep(" ", 18), collapse = " ")

forest_plot_df_41 <- forest_full_41 %>%
  dplyr::select(
    Subgroup = Subgroup_disp,
    Group = Group_disp,
    Number,
    `Deaths (%)`,
    `Risk (%)`,
    `RR vs non-AKD`,
    `RD vs non-AKD (pp)`,
    ` `,
    `HR (95% CI)`,
    `P-value`,
    `P for interaction`
  )

fp41 <- forestploter::forest(
  data  = forest_plot_df_41,
  est   = forest_full_41$HR,
  lower = forest_full_41$CI_low,
  upper = forest_full_41$CI_high,
  sizes = 0.6,
  ci_column = 8,
  ref_line  = 1,
  x_trans   = "log",
  xlim      = c(0.5, 12),
  ticks_at  = c(0.5, 1, 2, 4, 8),
  arrow_lab = c("", "")
)

fp41 <- edit_plot(
  fp41,
  col = c(3,4,5,6,7,9,10,11),
  which = "text",
  hjust = unit(1, "npc"),
  x = unit(1, "npc")
)

fp41 <- add_border(fp41, row = 0, where = "bottom", gp = gpar(lwd = 1))
fp41$heights <- rep(unit(6.2, "mm"), nrow(fp41))

fp41_annot <- gridExtra::arrangeGrob(
  fp41,
  bottom = grid::textGrob("Higher risk (left)    Lower risk (right)",
                          x = 0.72, hjust = 0.5,
                          gp = grid::gpar(fontsize = 9, fontfamily = base_family))
)

save_pdf_onepage_safe(fp41_annot, out_pdf_41,
                      page_w_mm = 380, page_h_mm = 195,
                      margin_left_mm = 16, margin_right_mm = 22,
                      margin_top_mm = 12, margin_bottom_mm = 34,
                      scale = 0.72,
                      family = base_family)

save_tiff_onepage_safe(fp41_annot, out_tiff_41,
                       page_w_mm = 380, page_h_mm = 195,
                       margin_left_mm = 16, margin_right_mm = 22,
                       margin_top_mm = 12, margin_bottom_mm = 34,
                       scale = 0.72,
                       dpi = 600,
                       family = base_family)

# ==========================
# Figure4-2 (HR-focused)
# ==========================
forest_full_42 <- forest_full
forest_full_42$` ` <- paste(rep(" ", 18), collapse = " ")

forest_plot_df_42 <- forest_full_42 %>%
  dplyr::select(
    Subgroup = Subgroup_disp,
    Group = Group_disp,
    Number,
    `Deaths (%)`,
    ` `,
    `HR (95% CI)`,
    `P-value`,
    `P for interaction`
  )

fp42 <- forestploter::forest(
  data  = forest_plot_df_42,
  est   = forest_full_42$HR,
  lower = forest_full_42$CI_low,
  upper = forest_full_42$CI_high,
  sizes = 0.6,
  ci_column = 5,
  ref_line  = 1,
  x_trans   = "log",
  xlim      = c(0.5, 12),
  ticks_at  = c(0.5, 1, 2, 4, 8),
  arrow_lab = c("", "")
)

fp42 <- edit_plot(
  fp42,
  col = c(3,4,6,7,8),
  which = "text",
  hjust = unit(1, "npc"),
  x = unit(1, "npc")
)

fp42 <- add_border(fp42, row = 0, where = "bottom", gp = gpar(lwd = 1))
fp42$heights <- rep(unit(6.2, "mm"), nrow(fp42))

fp42_annot <- gridExtra::arrangeGrob(
  fp42,
  bottom = grid::textGrob("Higher risk (left)    Lower risk (right)",
                          x = 0.66, hjust = 0.5,
                          gp = grid::gpar(fontsize = 9, fontfamily = base_family))
)

save_pdf_onepage_safe(fp42_annot, out_pdf_42,
                      page_w_mm = 330, page_h_mm = 185,
                      margin_left_mm = 16, margin_right_mm = 20,
                      margin_top_mm = 12, margin_bottom_mm = 34,
                      scale = 0.76,
                      family = base_family)

save_tiff_onepage_safe(fp42_annot, out_tiff_42,
                       page_w_mm = 330, page_h_mm = 185,
                       margin_left_mm = 16, margin_right_mm = 20,
                       margin_top_mm = 12, margin_bottom_mm = 34,
                       scale = 0.76,
                       dpi = 600,
                       family = base_family)

cat("\nDone.\n")
cat("Figure4-1:\n", out_pdf_41, "\n", out_tiff_41, "\n")
cat("Figure4-2:\n", out_pdf_42, "\n", out_tiff_42, "\n")

# ==========================
# Remove Overall (Figure4-3)
# ==========================
forest_full <- forest_full %>%
  filter(Subgroup != "Overall")
out_tiff <- file.path(outdir, "Figure4-3_ForestPlot_Subgroups_Sex_Age_CKD.tiff")
out_pdf  <- file.path(outdir, "Figure4-3_ForestPlot_Subgroups_Sex_Age_CKD.pdf")
forest_full <- forest_full %>%
  mutate(
    Subgroup = factor(Subgroup,
                      levels = c(
                        "Female","Male",
                        "<75","≧75",
                        "CKD no","CKD yes"
                      ))
  ) %>%
  arrange(Subgroup, Group)

