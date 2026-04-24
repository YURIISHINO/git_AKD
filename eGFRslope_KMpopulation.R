#KMはいった患者だけでeGFrslope
{
  ############################################################
  # eGFR slope among patients included in KM analysis
  # time0 = index_date + 210 days
  ############################################################
  
  graphics.off()
  
  # ==========================
  # Packages
  # ==========================
  pkgs <- c(
    "readr", "dplyr", "tidyr", "stringr", "purrr",
    "nlme", "multcomp", "ggplot2", "data.table", "ragg"
  )
  to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(to_install)) install.packages(to_install, dependencies = TRUE)
  
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(nlme)
  library(multcomp)
  library(ggplot2)
  library(data.table)
  library(ragg)
  
  # ==========================
  # Paths
  # ==========================
  setwd("X:/R")
  in_csv <- "jin1_Eligibile.csv"
  
  outdir <- "X:/R/eGFRslope"
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # ==========================
  # Load data
  # ==========================
  jin1_Eligibile <- read_csv(
    in_csv,
    locale = locale(encoding = "SHIFT-JIS"),
    show_col_types = FALSE
  )
  
  # ==========================
  # 1) Recreate KM analysis population
  # ==========================
  levels_full <- c(
    "non-AKD",
    "AKD with recovery",
    "AKD without recovery"
  )
  
  dat1 <- jin1_Eligibile %>%
    filter(exclude == "include", jin_status %in% c("AKD", "nonAKD")) %>%
    group_by(id) %>%
    arrange(index_date, date, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()
  
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
  
  km_ids <- dat_km %>%
    distinct(id, group, index_plus_210)
  
  print(table(dat_km$group, useNA = "ifany"))
  
  # ==========================
  # 2) Build eGFR longitudinal dataset restricted to KM patients
  # ==========================
  jin1_slope_km <- jin1_Eligibile %>%
    semi_join(km_ids, by = "id") %>%
    left_join(
      km_ids %>%
        transmute(id, group_km = group, time0 = index_plus_210),
      by = "id"
    ) %>%
    mutate(
      jin_label = case_when(
        group_km == "non-AKD" ~ "nonAKD",
        group_km == "AKD with recovery" ~ "Recovery",
        group_km == "AKD without recovery" ~ "Non-Recovery",
        TRUE ~ NA_character_
      ),
      jin_label = factor(jin_label, levels = c("nonAKD", "Recovery", "Non-Recovery")),
      years_from_time0 = as.numeric(date - time0) / 365.25
    ) %>%
    filter(
      !is.na(jin_label),
      !is.na(egfr),
      years_from_time0 >= 0
    )
  
  # time0 eGFR: index_plus_210当日の値があればそれを使う
  # なければ time0以降で最も近いeGFRをbaselineとして採用
  time0_egfr_df <- jin1_slope_km %>%
    mutate(abs_days_from_time0 = abs(as.numeric(date - time0))) %>%
    group_by(id) %>%
    arrange(abs_days_from_time0, date, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(id, time0_egfr = egfr)
  
  jin1_slope_km <- jin1_slope_km %>%
    left_join(time0_egfr_df, by = "id")
  
  # ==========================
  # 3) Baseline covariates
  # ==========================
  baseline_cov <- dat_km %>%
    mutate(
      arb_acei_use = if_else(coalesce(arb, 0) == 1 | coalesce(acei, 0) == 1, 1L, 0L)
    ) %>%
    dplyr::select(
      id, age, sex, arb_acei_use,
      dn1, dn3, dn4, dn5, dn6, dn7, dn8, dn9, dn10, dn12, dn13, dn14, dn15
    )
  
  covars <- c(
    "age", "sex", "arb_acei_use",
    "dn1", "dn3", "dn4", "dn5", "dn6", "dn7", "dn8", "dn9",
    "dn10", "dn12", "dn13", "dn14", "dn15"
  )
  
  longdat <- jin1_slope_km %>%
    dplyr::select(-any_of(covars)) %>%
    left_join(baseline_cov, by = "id") %>%
    mutate(
      age_c = as.numeric(scale(age, center = TRUE, scale = FALSE)),
      time0_egfr_c = as.numeric(scale(time0_egfr, center = TRUE, scale = FALSE)),
      jin_label = factor(jin_label, levels = c("nonAKD", "Recovery", "Non-Recovery"))
    ) %>%
    filter(!is.na(time0_egfr))
  
  # eGFR slope解析に入る人数確認
  n_by_group <- longdat %>%
    distinct(id, jin_label) %>%
    count(jin_label, name = "n")
  
  print(n_by_group)
  
  # ==========================
  # 4) LME model: slope-adjusted
  # ==========================
  ctrl <- lmeControl(
    maxIter = 1e8,
    msMaxIter = 1e8,
    opt = "optim",
    optimMethod = "L-BFGS-B"
  )
  
  fit_slope_adj_all <- lme(
    egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1 +
      age_c + sex + arb_acei_use +
      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15 +
      years_from_time0:(
        age_c + arb_acei_use +
          dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15
      ),
    random = list(id = pdSymm(~ 1 + years_from_time0)),
    data = longdat,
    na.action = na.omit,
    method = "REML",
    control = ctrl
  )
  
  fit_slope_adj_1y <- update(
    fit_slope_adj_all,
    data = filter(longdat, years_from_time0 <= 1)
  )
  
  fit_slope_adj_3y <- update(
    fit_slope_adj_all,
    data = filter(longdat, years_from_time0 <= 3)
  )
  
  fits <- list(
    "≤1 year" = fit_slope_adj_1y,
    "≤3 years" = fit_slope_adj_3y,
    "All period" = fit_slope_adj_all
  )
  
  # ==========================
  # 5) Extract slopes and CIs
  # ==========================
  get_slopes_ci <- function(fit, window_label) {
    cf <- names(fixef(fit))
    
    v_slope <- function(g) {
      vec <- rep(0, length(cf))
      names(vec) <- cf
      
      if ("years_from_time0" %in% cf) {
        vec["years_from_time0"] <- 1
      }
      if (g == "Recovery" && "years_from_time0:jin_labelRecovery" %in% cf) {
        vec["years_from_time0:jin_labelRecovery"] <- 1
      }
      if (g == "Non-Recovery" && "years_from_time0:jin_labelNon-Recovery" %in% cf) {
        vec["years_from_time0:jin_labelNon-Recovery"] <- 1
      }
      vec
    }
    
    L <- rbind(
      nonAKD = v_slope("nonAKD"),
      Recovery = v_slope("Recovery"),
      `Non-Recovery` = v_slope("Non-Recovery")
    )
    
    ci <- suppressMessages(confint(glht(fit, linfct = L)))
    
    tibble(
      window = window_label,
      group = rownames(L),
      estimate = ci$confint[, "Estimate"],
      lower = ci$confint[, "lwr"],
      upper = ci$confint[, "upr"]
    )
  }
  
  get_contrasts_vs_nonakd <- function(fit, window_label) {
    cf <- names(fixef(fit))
    
    v_slope <- function(g) {
      vec <- rep(0, length(cf))
      names(vec) <- cf
      
      if ("years_from_time0" %in% cf) {
        vec["years_from_time0"] <- 1
      }
      if (g == "Recovery" && "years_from_time0:jin_labelRecovery" %in% cf) {
        vec["years_from_time0:jin_labelRecovery"] <- 1
      }
      if (g == "Non-Recovery" && "years_from_time0:jin_labelNon-Recovery" %in% cf) {
        vec["years_from_time0:jin_labelNon-Recovery"] <- 1
      }
      vec
    }
    
    b <- v_slope("nonAKD")
    r <- v_slope("Recovery")
    nr <- v_slope("Non-Recovery")
    
    K <- rbind(
      `Recovery − nonAKD` = r - b,
      `Non-Recovery − nonAKD` = nr - b
    )
    
    gl <- glht(fit, linfct = K)
    ci <- suppressMessages(confint(gl))
    sm <- suppressMessages(summary(gl))
    
    tibble(
      window = window_label,
      contrast = rownames(K),
      group = case_when(
        contrast == "Recovery − nonAKD" ~ "Recovery",
        contrast == "Non-Recovery − nonAKD" ~ "Non-Recovery",
        TRUE ~ contrast
      ),
      diff_value = ci$confint[, "Estimate"],
      diff_lower = ci$confint[, "lwr"],
      diff_upper = ci$confint[, "upr"],
      diff_p = sm$test$pvalues
    )
  }
  
  df_plot <- imap_dfr(fits, ~get_slopes_ci(.x, .y)) %>%
    mutate(
      group = factor(group, levels = c("nonAKD", "Recovery", "Non-Recovery")),
      window = factor(window, levels = c("≤1 year", "≤3 years", "All period"))
    )
  
  df_contrast <- imap_dfr(fits, ~get_contrasts_vs_nonakd(.x, .y))
  
  sample_sizes <- bind_rows(
    longdat %>% filter(years_from_time0 <= 1) %>% mutate(window = "≤1 year"),
    longdat %>% filter(years_from_time0 <= 3) %>% mutate(window = "≤3 years"),
    longdat %>% mutate(window = "All period")
  ) %>%
    distinct(id, jin_label, window) %>%
    count(jin_label, window, name = "n") %>%
    transmute(
      group = as.character(jin_label),
      window = factor(window, levels = c("≤1 year", "≤3 years", "All period")),
      n
    )
  
  diff_ref <- expand.grid(
    window = factor(c("≤1 year", "≤3 years", "All period"),
                    levels = c("≤1 year", "≤3 years", "All period")),
    group = "nonAKD",
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    mutate(
      diff_value = NA_real_,
      diff_lower = NA_real_,
      diff_upper = NA_real_,
      diff_p = NA_real_
    )
  
  differences <- bind_rows(
    diff_ref,
    df_contrast %>%
      dplyr::select(window, group, diff_value, diff_lower, diff_upper, diff_p) %>%
      mutate(window = factor(window, levels = c("≤1 year", "≤3 years", "All period")))
  )
  
  df_fig <- df_plot %>%
    left_join(sample_sizes, by = c("group", "window")) %>%
    left_join(differences, by = c("group", "window")) %>%
    mutate(
      group = factor(group, levels = c("nonAKD", "Recovery", "Non-Recovery"))
    )
  
  write_csv(df_fig, file.path(outdir, "eGFR_slope_KM_population_results.csv"))
  write_csv(df_contrast, file.path(outdir, "eGFR_slope_KM_population_contrasts.csv"))
  
  # ==========================
  # 6) Plot: slope-adjusted eGFR slope
  # ==========================
  make_pub_plot <- function(df_in, facet_by_window = TRUE) {
    
    y_min <- min(df_in$lower, na.rm = TRUE)
    y_max <- max(df_in$upper, na.rm = TRUE)
    
    diff_y <- y_min - 3.2
    p_y <- y_min - 5.2
    
    p <- ggplot(df_in, aes(x = group, y = estimate, fill = group)) +
      geom_col(width = 0.7, alpha = 0.85) +
      geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, linewidth = 0.7) +
      geom_text(
        aes(y = pmax(upper, 0) + 0.5, label = paste0("n=", n)),
        size = 3.5,
        fontface = "bold",
        color = "grey20"
      ) +
      geom_text(
        aes(
          y = diff_y,
          label = ifelse(
            group == "nonAKD",
            "Reference",
            sprintf("Diff: %.2f\n(%.2f, %.2f)", diff_value, diff_lower, diff_upper)
          )
        ),
        size = 3,
        lineheight = 0.95,
        fontface = "italic",
        color = "grey30"
      ) +
      geom_text(
        aes(
          y = p_y,
          label = ifelse(
            group == "nonAKD",
            "",
            case_when(
              is.na(diff_p) ~ "",
              diff_p < 0.001 ~ "p<0.001",
              diff_p < 0.01 ~ sprintf("p=%.3f", diff_p),
              TRUE ~ sprintf("p=%.2f", diff_p)
            )
          )
        ),
        size = 2.9,
        fontface = "bold",
        color = "grey20"
      ) +
      scale_fill_manual(
        values = c(
          nonAKD = "#95A5A6",
          Recovery = "#2ECC71",
          `Non-Recovery` = "#E74C3C"
        ),
        labels = c(
          nonAKD = "non-AKD",
          Recovery = "AKD with recovery",
          `Non-Recovery` = "AKD without recovery"
        )
      ) +
      scale_x_discrete(expand = expansion(add = 0.6)) +
      coord_cartesian(ylim = c(y_min - 6.8, y_max + 1.5), clip = "off") +
      labs(
        x = NULL,
        y = "Mean change in eGFR (mL/min/1.73 m² per year)",
        fill = "Group"
      ) +
      theme_classic(base_size = 13) +
      theme(
        panel.grid = element_blank(),
        panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.6),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 12, margin = margin(r = 10)),
        axis.text.y = element_text(size = 11),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1.2, "lines"),
        legend.background = element_blank(),
        plot.margin = margin(t = 5, r = 20, b = 18, l = 25)
      )
    
    if (facet_by_window) {
      p <- p + facet_grid(. ~ window)
    }
    
    p
  }
  
  p_slope_3windows <- make_pub_plot(df_fig, facet_by_window = TRUE)
  
  p_slope_1y <- make_pub_plot(
    df_fig %>% filter(window == "≤1 year"),
    facet_by_window = FALSE
  ) +
    ggtitle("Slope-adjusted eGFR slope among KM population (≤1 year)")
  
  print(p_slope_3windows)
  print(p_slope_1y)
  
  ggsave(
    file.path(outdir, "eGFR_slope_KM_population_3windows.pdf"),
    plot = p_slope_3windows,
    width = 200,
    height = 140,
    units = "mm",
    device = cairo_pdf,
    dpi = 300
  )
  
  ggsave(
    file.path(outdir, "eGFR_slope_KM_population_3windows.tiff"),
    plot = p_slope_3windows,
    width = 200,
    height = 140,
    units = "mm",
    device = "tiff",
    dpi = 600,
    compression = "lzw"
  )
  
  ggsave(
    file.path(outdir, "eGFR_slope_KM_population_1year.pdf"),
    plot = p_slope_1y,
    width = 120,
    height = 120,
    units = "mm",
    device = cairo_pdf,
    dpi = 300
  )
  
  ggsave(
    file.path(outdir, "eGFR_slope_KM_population_1year.tiff"),
    plot = p_slope_1y,
    width = 120,
    height = 120,
    units = "mm",
    device = "tiff",
    dpi = 600,
    compression = "lzw"
  )
}

############################################################
# KMに入ったが slopeで除外された患者の抽出
############################################################

# ==========================
# 1) KMに入ったID
# ==========================
km_id_df <- dat_km %>%
  distinct(id)

# ==========================
# 2) slopeに入ったID（LMEに実際に入ったもの）
# ==========================
slope_id_df <- longdat %>%
  distinct(id)

# ==========================
# 3) KMにはいるが slopeにはいないID
# ==========================
excluded_ids <- km_id_df %>%
  anti_join(slope_id_df, by = "id")

# 確認
n_excluded <- nrow(excluded_ids)
print(paste0("Excluded from slope: ", n_excluded))

# ==========================
# 4) death_date / last_data を取得
# ==========================
excluded_detail <- jin1_Eligibile %>%
  semi_join(excluded_ids, by = "id") %>%
  group_by(id) %>%
  summarise(
    index_date = first(index_date),
    death_date = first(death_date),
    last_data  = first(last_data),
    .groups = "drop"
  )


# ==========================
# 5) 出力
# ==========================
print(head(excluded_detail, 20))