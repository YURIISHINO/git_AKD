#年齢曲線#####
{
  ############################################################
  # Supplement Figure 5. Continuous age × AKD interaction (Cox)
  # - HR(age) curves for AKD with recovery / AKD without recovery
  #   versus non-AKD
  # - Analysis population aligned with Figure 4
  # - Save as PDF + TIFF in figure_table
  ############################################################
  
  graphics.off()
  
  # ==========================
  # Packages
  # ==========================
  pkgs <- c("readr", "dplyr", "survival", "ggplot2", "tidyr", "ragg")
  to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(to_install)) install.packages(to_install, dependencies = TRUE)
  
  library(readr)
  library(dplyr)
  library(survival)
  library(ggplot2)
  library(tidyr)
  library(ragg)
  
  # ==========================
  # Paths
  # ==========================
  setwd("X:/R")
  in_csv <- "jin1_Eligibile.csv"
  outdir <- file.path(getwd(), "word_supp_tables")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  out_pdf  <- file.path(outdir, "Supplement Figure5_AgeContinuousInteraction.pdf")
  out_tiff <- file.path(outdir, "Supplement Figure5_AgeContinuousInteraction.tiff")
  
  # ==========================
  # Font
  # ==========================
  base_family <- {
    f <- c("Yu Gothic", "MS Gothic", "Meiryo", "Arial Unicode MS", "Arial")
    ok <- f[f %in% names(grDevices::windowsFonts())]
    if (length(ok) == 0) "sans" else ok[1]
  }
  
  # ==========================
  # 0) Load
  # ==========================
  loc <- locale(encoding = "SHIFT-JIS")
  jin1_Eligibile <- read_csv(in_csv, locale = loc)
  
  # ==========================
  # 1) Build analysis dataset
  #    (aligned with Figure 4)
  # ==========================
  dat <- jin1_Eligibile %>%
    group_by(id) %>%
    arrange(index_date, date, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    filter(exclude == "include",
           jin_status %in% c("nonAKD", "AKD")) %>%
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
      time_years = as.numeric(last_follow_death - index_plus_210) / 365.25
    ) %>%
    filter(!is.na(jin_label),
           !is.na(time_years), time_years >= 0,
           !is.na(primary_death),
           !is.na(age),
           !is.na(index_cre))
  
  # ==========================
  # 2) Center age
  # ==========================
  age_mean <- mean(dat$age, na.rm = TRUE)
  dat <- dat %>%
    mutate(age_c = age - age_mean)
  
  # ==========================
  # 3) Cox models
  #    - base model
  #    - interaction model
  # ==========================
  cox_base_age_cont <- coxph(
    Surv(time_years, primary_death) ~
      jin_label + age_c + index_cre + arb_acei_use +
      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 +
      dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
    data = dat
  )
  
  cox_int_age_cont <- coxph(
    Surv(time_years, primary_death) ~
      jin_label * age_c + index_cre + arb_acei_use +
      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 +
      dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
    data = dat
  )
  
  cat("\nLikelihood ratio test for age interaction:\n")
  print(anova(cox_base_age_cont, cox_int_age_cont, test = "LRT"))
  
  # ==========================
  # 4) Extract coef / vcov
  # ==========================
  cf <- coef(cox_int_age_cont)
  vc <- vcov(cox_int_age_cont)
  
  cat("\nRelevant coefficient names:\n")
  print(names(cf)[grepl("^jin_label|age_c", names(cf))])
  
  # ==========================
  # 5) Function to calculate HR(age)
  # ==========================
  calc_hr_age <- function(age_val, b_main, b_int, age_mean, cf, vc) {
    age_c_val <- as.numeric(age_val - age_mean)
    
    est <- as.numeric(cf[b_main] + cf[b_int] * age_c_val)
    
    v1 <- as.numeric(vc[b_main, b_main])
    v2 <- as.numeric(vc[b_int,  b_int])
    v3 <- as.numeric(vc[b_main, b_int])
    
    var <- as.numeric(v1 + (age_c_val^2) * v2 + 2 * age_c_val * v3)
    var <- pmax(var, 0)
    se  <- sqrt(var)
    
    HR <- exp(est)
    lo <- exp(est - 1.96 * se)
    hi <- exp(est + 1.96 * se)
    
    c(HR = HR, lo = lo, hi = hi)
  }
  
  # ==========================
  # 6) Build plot dataset
  # ==========================
  age_min <- floor(quantile(dat$age, 0.01, na.rm = TRUE))
  age_max <- ceiling(quantile(dat$age, 0.99, na.rm = TRUE))
  age_seq <- seq(age_min, age_max, by = 1)
  
  plot_df <- do.call(rbind, lapply(age_seq, function(a) {
    rec  <- calc_hr_age(a, "jin_labelRecovery", "jin_labelRecovery:age_c",
                        age_mean = age_mean, cf = cf, vc = vc)
    nonr <- calc_hr_age(a, "jin_labelNon-Recovery", "jin_labelNon-Recovery:age_c",
                        age_mean = age_mean, cf = cf, vc = vc)
    
    data.frame(
      age = a,
      Recovery_HR = rec["HR"],
      Recovery_lo = rec["lo"],
      Recovery_hi = rec["hi"],
      NonRecovery_HR = nonr["HR"],
      NonRecovery_lo = nonr["lo"],
      NonRecovery_hi = nonr["hi"]
    )
  })) %>%
    as.data.frame()
  
  plot_long <- plot_df %>%
    pivot_longer(
      cols = -age,
      names_to = c("Group", ".value"),
      names_pattern = "(Recovery|NonRecovery)_(HR|lo|hi)"
    ) %>%
    mutate(
      Group = factor(
        Group,
        levels = c("Recovery", "NonRecovery"),
        labels = c("AKD with recovery", "AKD without recovery")
      ),
      line_type = if_else(Group == "AKD with recovery", "solid", "dashed")
    ) %>%
    filter(is.finite(HR), is.finite(lo), is.finite(hi))
  
  stopifnot(nrow(plot_long) > 0)
  
  # ==========================
  # 7) Axis range
  # ==========================
  ymin <- min(plot_long$lo, 1, na.rm = TRUE)
  ymax <- max(plot_long$hi, 1, na.rm = TRUE)
  
  # 見やすさのため少し余裕を持たせる
  ymin <- max(0, ymin * 0.95)
  ymax <- ymax * 1.05
  
  # ==========================
  # 9) Plot
  # ==========================
  col_rec  <- "#1B9E77"
  col_nonr <- "#D95F02"
  
  p_age_int <- ggplot(plot_long, aes(x = age, y = HR, color = Group, linetype = line_type)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = Group),
                alpha = 0.15, color = NA, show.legend = FALSE) +
    geom_line(linewidth = 1.15) +
    geom_hline(yintercept = 1, linetype = "dotted", linewidth = 0.7) +
    coord_cartesian(ylim = c(ymin, ymax)) +
    scale_color_manual(values = c(col_rec, col_nonr)) +
    scale_fill_manual(values = c(col_rec, col_nonr)) +
    scale_linetype_identity() +
    labs(
      x = "Age (years)",
      y = "Hazard ratio vs non-AKD",
      color = "AKD group"
    ) +
    guides(
      color = guide_legend(override.aes = list(linewidth = 1.8))
    ) +
    theme_minimal(base_family = base_family, base_size = 13) +
    theme(
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.caption = element_text(
        size = 10.5,
        hjust = 0,
        lineheight = 1.12,
        margin = margin(t = 8)
      ),
      plot.margin = margin(t = 10, r = 10, b = 18, l = 10)
    )
  
  print(p_age_int)
  
  # ==========================
  # 10) Per +1 year multiplicative change in HR
  # ==========================
  mult_rec  <- exp(cf["jin_labelRecovery:age_c"])
  mult_nonr <- exp(cf["jin_labelNon-Recovery:age_c"])
  
  cat("\nPer +1 year multiplicative change in HR (vs non-AKD):\n")
  cat("  AKD with recovery    :", mult_rec,  "\n")
  cat("  AKD without recovery :", mult_nonr, "\n")
  
  # ==========================
  # 11) Save
  # ==========================
  pdf(out_pdf, width = 180/25.4, height = 120/25.4, family = base_family)
  print(p_age_int)
  dev.off()
  
  ragg::agg_tiff(
    out_tiff,
    width = 180/25.4,
    height = 120/25.4,
    units = "in",
    res = 600,
    compression = "lzw"
  )
  print(p_age_int)
  dev.off()
  
  cat("\nSaved files:\n", out_pdf, "\n", out_tiff, "\n")
  
  ###proportional hazards（PH仮定）
  cox.zph(cox_int_age_cont)
  plot(cox.zph(cox_int_age_cont))
  
}
#年齢ごとのフォレストプロットのみ
{
  ############################################################
  # Supplementary Figure 5
  # Forest-style plot of age-specific HRs from continuous
  # age × AKD-group interaction Cox model
  # Displayed as:
  #   <40 / 40-49 / 50-59 / 60-69 / 70-79 / ≥80
  ############################################################
  
  graphics.off()
  
  # ==========================
  # Packages
  # ==========================
  pkgs <- c("readr", "dplyr", "survival", "ggplot2", "tibble", "ragg")
  to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
  
  library(readr)
  library(dplyr)
  library(survival)
  library(ggplot2)
  library(tibble)
  library(ragg)
  
  # ==========================
  # Paths
  # ==========================
  setwd("X:/R")
  in_csv <- "jin1_Eligibile.csv"
  outdir <- file.path(getwd(), "word_supp_tables")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  out_pdf  <- file.path(outdir, "Supplementary_Figure5_Forest_AgeBands.pdf")
  out_tiff <- file.path(outdir, "Supplementary_Figure5_Forest_AgeBands.tif")
  out_csv  <- file.path(outdir, "Supplementary_Figure5_Forest_AgeBands_values.csv")
  
  # ==========================
  # Font
  # ==========================
  base_family <- {
    f <- c("Yu Gothic", "MS Gothic", "Meiryo", "Arial Unicode MS", "Arial")
    ok <- f[f %in% names(grDevices::windowsFonts())]
    if (length(ok) == 0) "sans" else ok[1]
  }
  
  # ==========================
  # 1) Load
  # ==========================
  loc <- locale(encoding = "SHIFT-JIS")
  jin1_Eligibile <- read_csv(in_csv, locale = loc)
  
  # ==========================
  # 2) Build analysis dataset
  # ==========================
  dat <- jin1_Eligibile %>%
    group_by(id) %>%
    arrange(index_date, date, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    filter(
      exclude == "include",
      jin_status %in% c("AKD", "nonAKD")
    ) %>%
    mutate(
      group = case_when(
        jin_status == "nonAKD" ~ "non-AKD",
        jin_status == "AKD" & `150_210recovery` == 1 ~ "AKD with recovery",
        jin_status == "AKD" & `150_210recovery` == 2 ~ "AKD without recovery",
        jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 1 ~ "AKD with recovery",
        jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` %in% c(0, 2) ~ "AKD without recovery",
        TRUE ~ NA_character_
      ),
      group = factor(
        group,
        levels = c("non-AKD", "AKD with recovery", "AKD without recovery")
      ),
      arb_acei_use = if_else(coalesce(arb, 0) == 1 | coalesce(acei, 0) == 1, 1L, 0L),
      time_years = as.numeric(last_follow_death - index_plus_210) / 365.25
    ) %>%
    filter(
      !is.na(group),
      !is.na(age),
      !is.na(index_cre),
      !is.na(time_years), time_years >= 0,
      !is.na(primary_death)
    )
  
  cat("\nGroup counts:\n")
  print(table(dat$group, useNA = "ifany"))
  
  # ==========================
  # 3) Center age
  # ==========================
  age_mean <- mean(dat$age, na.rm = TRUE)
  dat <- dat %>%
    mutate(age_c = age - age_mean)
  
  # ==========================
  # 4) Cox model with age interaction
  # ==========================
  cox_model <- coxph(
    Surv(time_years, primary_death) ~
      group * age_c +
      index_cre + arb_acei_use +
      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 +
      dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
    data = dat
  )
  
  cf <- coef(cox_model)
  vc <- vcov(cox_model)
  
  cat("\nCoefficient names:\n")
  print(names(cf))
  
  # ==========================
  # 5) Age-band display labels and representative ages
  # ==========================
  rep_age_lt40 <- dat %>%
    filter(age < 40) %>%
    summarise(rep_age = median(age, na.rm = TRUE)) %>%
    pull(rep_age)
  
  rep_age_ge80 <- dat %>%
    filter(age >= 80) %>%
    summarise(rep_age = median(age, na.rm = TRUE)) %>%
    pull(rep_age)
  
  if (length(rep_age_lt40) == 0 || is.na(rep_age_lt40)) rep_age_lt40 <- 35
  if (length(rep_age_ge80) == 0 || is.na(rep_age_ge80)) rep_age_ge80 <- 82
  
  age_band_df <- tibble(
    band_id    = c("lt40", "40_49", "50_59", "60_69", "70_79", "ge80"),
    band_label = c("<40", "40-49", "50-59", "60-69", "70-79", ">=80"),
    rep_age    = c(rep_age_lt40, 45, 55, 65, 75, rep_age_ge80)
  )
  
  cat("\nRepresentative ages used:\n")
  print(age_band_df)
  
  # ==========================
  # 6) Function to calculate age-specific HR
  # ==========================
  calc_hr_at_age <- function(group_label, age_value, age_mean, cf, vc) {
    coef_main <- paste0("group", group_label)
    coef_int  <- paste0("group", group_label, ":age_c")
    
    beta_main <- unname(cf[coef_main])
    beta_int  <- unname(cf[coef_int])
    
    var_main <- vc[coef_main, coef_main]
    var_int  <- vc[coef_int,  coef_int]
    cov_mi   <- vc[coef_main, coef_int]
    
    age_c_val <- age_value - age_mean
    
    log_hr <- beta_main + beta_int * age_c_val
    se_log_hr <- sqrt(pmax(0, var_main + (age_c_val^2) * var_int + 2 * age_c_val * cov_mi))
    
    tibble(
      rep_age = age_value,
      Group   = group_label,
      HR      = exp(log_hr),
      lower   = exp(log_hr - 1.96 * se_log_hr),
      upper   = exp(log_hr + 1.96 * se_log_hr)
    )
  }
  
  # ==========================
  # 7) Build plotting table
  # ==========================
  plot_df <- bind_rows(
    lapply(seq_len(nrow(age_band_df)), function(i) {
      bind_rows(
        calc_hr_at_age("AKD with recovery",    age_band_df$rep_age[i], age_mean, cf, vc),
        calc_hr_at_age("AKD without recovery", age_band_df$rep_age[i], age_mean, cf, vc)
      ) %>%
        mutate(
          age_band_id    = age_band_df$band_id[i],
          age_band_label = age_band_df$band_label[i]
        )
    })
  ) %>%
    mutate(
      age_band_id = factor(age_band_id, levels = rev(age_band_df$band_id)),
      Group = factor(
        Group,
        levels = c("AKD with recovery", "AKD without recovery")
      ),
      label_ci = sprintf("%.2f (%.2f-%.2f)", HR, lower, upper)
    )
  print(plot_df)
  
  # ==========================
  # 8) Save values
  # ==========================
  write_csv(plot_df, out_csv)
  
  # ==========================
  # 9) Forest-style plot
  # ==========================
  plot_df <- plot_df %>%
    mutate(
      y = case_when(
        Group == "AKD with recovery"    ~ as.numeric(age_band_id) + 0.16,
        Group == "AKD without recovery" ~ as.numeric(age_band_id) - 0.16
      )
    )
  
  x_min <- min(plot_df$lower, na.rm = TRUE)
  x_max <- max(plot_df$upper, na.rm = TRUE)
  
  # 対数軸
  x_lower <- min(0.5, floor(x_min * 10) / 10)
  x_upper <- max(4.0, ceiling(x_max * 10) / 10)
  label_x <- x_upper * 1.08
  
  p <- ggplot(plot_df, aes(x = HR, y = y, color = Group)) +
    geom_vline(xintercept = 1, linetype = "dotted", linewidth = 0.7, color = "black") +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, linewidth = 0.8) +
    geom_point(aes(shape = Group), size = 2.8, stroke = 0.9, fill = "white") +
    geom_text(
      aes(x = label_x, label = label_ci),
      hjust = 0,
      size = 3.3,
      family = base_family,
      color = "black"
    ) +
    scale_color_manual(
      values = c(
        "AKD with recovery"    = "#2ECC71",
        "AKD without recovery" = "#E74C3C"
      )
    ) +
    scale_shape_manual(
      values = c(
        "AKD with recovery"    = 21,
        "AKD without recovery" = 24
      )
    ) +
    scale_x_log10(
      limits = c(x_lower, label_x * 1.08),
      breaks = c(0.5, 1, 2, 4),
      labels = c("0.5", "1", "2", "4")
    ) +
    scale_y_continuous(
      breaks = seq_along(levels(plot_df$age_band_id)),
      labels = rev(age_band_df$band_label),
      expand = expansion(mult = c(0.08, 0.08))
    )+
    labs(
      x = "Hazard ratio vs non-AKD",
      y = "Age group (years)",
      color = "AKD group",
      shape = "AKD group"
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_family = base_family, base_size = 12) +
    theme(
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(t = 10, r = 90, b = 10, l = 10)
    )
  
  print(p)
  
  # ==========================
  # 10) Save
  # ==========================
  pdf(out_pdf, width = 190/25.4, height = 140/25.4, family = base_family)
  print(p)
  dev.off()
  
  ragg::agg_tiff(
    out_tiff,
    width = 190/25.4,
    height = 140/25.4,
    units = "in",
    res = 600,
    compression = "lzw"
  )
  print(p)
  dev.off()
  
  cat("\nSaved files:\n")
  cat(" PDF : ", out_pdf, "\n")
  cat(" TIFF: ", out_tiff, "\n")
  cat(" CSV : ", out_csv, "\n")} 


