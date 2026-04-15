#adjusted cox:original版→おもい####
{
  ############################################################
  # Figure 4. Subgroup Forest Plot (Sex / Age / CKD)
  # Final revised version:
  # - RR column removed
  # - Adjusted RD at 5 years added by G-computation
  # - Bootstrap CI uses bootstrapped refitted Cox model correctly
  # - n_boot_rd = 1000
  # - HR: original full follow-up Cox model
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
  
  out_pdf_4  <- file.path(outdir, "Figure4_ForestPlot_Subgroups_Sex_Age_CKD_adjRD.pdf")
  out_tiff_4 <- file.path(outdir, "Figure4_ForestPlot_Subgroups_Sex_Age_CKD_adjRD.tiff")
  
  # ==========================
  # Font
  # ==========================
  base_family <- {
    f <- c("Yu Gothic", "MS Gothic", "Meiryo", "Arial Unicode MS", "Arial")
    ok <- f[f %in% names(grDevices::windowsFonts())]
    if (length(ok) == 0) "sans" else ok[1]
  }
  
  # ==========================
  # Parameters
  # ==========================
  rd_horizon <- 5       # 最終版では 5-year RD として固定
  n_boot_rd  <- 5    # 最終版
  
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
                  "ge75" = ">=75")
  }
  
  disp_ckdgrp <- function(x){
    dplyr::recode(as.character(x),
                  "0" = "CKD no",
                  "1" = "CKD yes")
  }
  
  save_pdf_onepage_safe <- function(grob, filename,
                                    page_w_mm = 380,
                                    page_h_mm = 235,
                                    margin_left_mm = 24,
                                    margin_right_mm = 24,
                                    margin_top_mm = 24,
                                    margin_bottom_mm = 42,
                                    scale = 0.58,
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
                                     page_w_mm = 380,
                                     page_h_mm = 235,
                                     margin_left_mm = 24,
                                     margin_right_mm = 24,
                                     margin_top_mm = 24,
                                     margin_bottom_mm = 42,
                                     scale = 0.58,
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
  # G-computation for adjusted RD
  # ==========================
  gcomp_rd <- function(fit, data, horizon, group_var,
                       ref_group = "nonAKD",
                       target_groups = c("Recovery","Non-Recovery"),
                       n_boot = 1000){
    
    marginal_risk <- function(fit_obj, dat, grp){
      dat_cf <- dat
      dat_cf[[group_var]] <- factor(grp, levels = levels(dat[[group_var]]))
      
      sf <- survfit(fit_obj, newdata = dat_cf)
      s  <- summary(sf, times = horizon, extend = TRUE)
      
      mean(1 - s$surv, na.rm = TRUE)
    }
    
    # point estimate
    risk_ref <- marginal_risk(fit, data, ref_group)
    
    rd_point <- sapply(target_groups, function(g){
      marginal_risk(fit, data, g) - risk_ref
    })
    
    # bootstrap CI
    rd_boot <- replicate(n_boot, {
      idx <- sample(seq_len(nrow(data)), replace = TRUE)
      boot_dat <- data[idx, , drop = FALSE]
      
      # factor levels を元データに合わせる
      boot_dat[[group_var]] <- factor(boot_dat[[group_var]], levels = levels(data[[group_var]]))
      
      boot_fit <- update(fit, data = boot_dat)
      
      risk_ref_b <- marginal_risk(boot_fit, boot_dat, ref_group)
      
      sapply(target_groups, function(g){
        marginal_risk(boot_fit, boot_dat, g) - risk_ref_b
      })
    })
    
    if (is.null(dim(rd_boot))) {
      rd_boot <- matrix(rd_boot, nrow = length(target_groups))
    }
    
    tibble(
      Group   = target_groups,
      RD      = as.numeric(rd_point),
      RD_low  = apply(rd_boot, 1, quantile, probs = 0.025, na.rm = TRUE),
      RD_high = apply(rd_boot, 1, quantile, probs = 0.975, na.rm = TRUE)
    )
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
  # 2) Overall = same as original Figure 4
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
  # 4) Interaction models for HR
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
        .groups = "drop"
      ) %>%
      mutate(Subgroup = "Overall",
             Group = as.character(group)) %>%
      select(Subgroup, Group, Number, Deaths, Death_pct)
  }
  
  calc_stats_sub <- function(data, subgroup_var, subgroup_lab_fun = NULL){
    out <- data %>%
      group_by(!!sym(subgroup_var), group) %>%
      summarise(
        Number = n(),
        Deaths = sum(primary_death),
        Death_pct = sprintf("%.1f", 100 * Deaths / Number),
        .groups = "drop"
      ) %>%
      mutate(
        Subgroup = as.character(!!sym(subgroup_var)),
        Group = as.character(group)
      )
    if (!is.null(subgroup_lab_fun)) out$Subgroup <- subgroup_lab_fun(out$Subgroup)
    out %>% select(Subgroup, Group, Number, Deaths, Death_pct)
  }
  
  km_risk5_overall <- function(data, t0 = 5){
    sf <- survfit(Surv(time_years, primary_death) ~ group, data = data)
    s  <- summary(sf, times = t0, extend = TRUE)
    
    tibble(
      strata = s$strata,
      surv   = s$surv
    ) %>%
      mutate(
        Group = sub("^group=", "", strata),
        Subgroup = "Overall",
        risk5 = 1 - surv
      ) %>%
      select(Subgroup, Group, risk5)
  }
  
  km_risk5_sub <- function(data, subgroup_var, subgroup_lab_fun = NULL, t0 = 5){
    subs <- sort(unique(as.character(data[[subgroup_var]])))
    
    bind_rows(lapply(subs, function(sv){
      df_sub <- data %>% filter(as.character(.data[[subgroup_var]]) == sv)
      sf <- survfit(Surv(time_years, primary_death) ~ group, data = df_sub)
      s  <- summary(sf, times = t0, extend = TRUE)
      
      sublab <- if (is.null(subgroup_lab_fun)) sv else subgroup_lab_fun(sv)
      
      tibble(
        strata = s$strata,
        surv   = s$surv
      ) %>%
        mutate(
          Group = sub("^group=", "", strata),
          Subgroup = sublab,
          risk5 = 1 - surv
        ) %>%
        select(Subgroup, Group, risk5)
    }))
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
  # 6) Stats + crude 5y risk + HR
  # ==========================
  stats_overall <- calc_stats_overall(dat_overall)
  stats_sex <- calc_stats_sub(dat_sex, "sex_f", function(x) dplyr::recode(x, "F" = "Female", "M" = "Male"))
  stats_age <- calc_stats_sub(dat_age, "age75_f", disp_agegrp)
  stats_ckd <- calc_stats_sub(dat_ckd, "ckd_bin", disp_ckdgrp)
  
  risk5_overall <- km_risk5_overall(dat_overall, t0 = rd_horizon)
  risk5_sex <- km_risk5_sub(dat_sex, "sex_f", function(x) dplyr::recode(x, "F" = "Female", "M" = "Male"), t0 = rd_horizon)
  risk5_age <- km_risk5_sub(dat_age, "age75_f", disp_agegrp, t0 = rd_horizon)
  risk5_ckd <- km_risk5_sub(dat_ckd, "ckd_bin", disp_ckdgrp, t0 = rd_horizon)
  
  hr_overall <- extract_hr_overall(fit_overall)
  
  hr_sex0 <- hr_from_interaction(fit_sex, subgroup_levels = c("Female","Male"), int_suffix_other = "sex_fM")
  hr_sex <- bind_rows(
    tibble(Subgroup = "Female", Group = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
    filter(hr_sex0, Subgroup == "Female"),
    tibble(Subgroup = "Male", Group = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
    filter(hr_sex0, Subgroup == "Male")
  )
  
  hr_age0 <- hr_from_interaction(fit_age, subgroup_levels = c("<75",">=75"), int_suffix_other = "age75_fge75")
  hr_age <- bind_rows(
    tibble(Subgroup = "<75", Group = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
    filter(hr_age0, Subgroup == "<75"),
    tibble(Subgroup = ">=75", Group = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
    filter(hr_age0, Subgroup == ">=75")
  )
  
  hr_ckd0 <- hr_from_interaction(fit_ckd, subgroup_levels = c("CKD no","CKD yes"), int_suffix_other = "ckd_bin")
  hr_ckd <- bind_rows(
    tibble(Subgroup = "CKD no", Group = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
    filter(hr_ckd0, Subgroup == "CKD no"),
    tibble(Subgroup = "CKD yes", Group = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
    filter(hr_ckd0, Subgroup == "CKD yes")
  )
  
  # ==========================
  # 7) Adjusted RD by subgroup
  # ==========================
  rd_overall <- gcomp_rd(
    fit_overall,
    dat_overall,
    horizon = rd_horizon,
    group_var = "group",
    n_boot = n_boot_rd
  ) %>% mutate(Subgroup = "Overall")
  
  rd_female <- gcomp_rd(
    update(fit_overall, data = dat_sex %>% filter(sex_f == "F")),
    dat_sex %>% filter(sex_f == "F"),
    horizon = rd_horizon,
    group_var = "group",
    n_boot = n_boot_rd
  ) %>% mutate(Subgroup = "Female")
  
  rd_male <- gcomp_rd(
    update(fit_overall, data = dat_sex %>% filter(sex_f == "M")),
    dat_sex %>% filter(sex_f == "M"),
    horizon = rd_horizon,
    group_var = "group",
    n_boot = n_boot_rd
  ) %>% mutate(Subgroup = "Male")
  
  rd_lt75 <- gcomp_rd(
    update(fit_overall, data = dat_age %>% filter(age75_f == "lt75")),
    dat_age %>% filter(age75_f == "lt75"),
    horizon = rd_horizon,
    group_var = "group",
    n_boot = n_boot_rd
  ) %>% mutate(Subgroup = "<75")
  
  rd_ge75 <- gcomp_rd(
    update(fit_overall, data = dat_age %>% filter(age75_f == "ge75")),
    dat_age %>% filter(age75_f == "ge75"),
    horizon = rd_horizon,
    group_var = "group",
    n_boot = n_boot_rd
  ) %>% mutate(Subgroup = ">=75")
  
  rd_ckd_no <- gcomp_rd(
    update(fit_overall, data = dat_ckd %>% filter(ckd_bin == 0)),
    dat_ckd %>% filter(ckd_bin == 0),
    horizon = rd_horizon,
    group_var = "group",
    n_boot = n_boot_rd
  ) %>% mutate(Subgroup = "CKD no")
  
  rd_ckd_yes <- gcomp_rd(
    update(fit_overall, data = dat_ckd %>% filter(ckd_bin == 1)),
    dat_ckd %>% filter(ckd_bin == 1),
    horizon = rd_horizon,
    group_var = "group",
    n_boot = n_boot_rd
  ) %>% mutate(Subgroup = "CKD yes")
  
  rd_all <- bind_rows(
    rd_overall, rd_female, rd_male, rd_lt75, rd_ge75, rd_ckd_no, rd_ckd_yes
  )
  
  # nonAKD row を追加
  rd_ref <- tibble(
    Subgroup = c("Overall","Female","Male","<75",">=75","CKD no","CKD yes"),
    Group = "nonAKD",
    RD = NA_real_,
    RD_low = NA_real_,
    RD_high = NA_real_
  )
  
  rd_all <- bind_rows(rd_ref, rd_all)
  
  # ==========================
  # 8) Merge
  # ==========================
  build_block <- function(stats_df, risk5_df, hr_df, rd_df){
    stats_df %>%
      left_join(risk5_df, by = c("Subgroup","Group")) %>%
      left_join(hr_df,    by = c("Subgroup","Group")) %>%
      left_join(rd_df,    by = c("Subgroup","Group"))
  }
  
  forest_full <- bind_rows(
    build_block(stats_overall, risk5_overall, hr_overall, filter(rd_all, Subgroup == "Overall")),
    build_block(stats_sex,     risk5_sex,     hr_sex,     filter(rd_all, Subgroup %in% c("Female","Male"))),
    build_block(stats_age,     risk5_age,     hr_age,     filter(rd_all, Subgroup %in% c("<75",">=75"))),
    build_block(stats_ckd,     risk5_ckd,     hr_ckd,     filter(rd_all, Subgroup %in% c("CKD no","CKD yes")))
  ) %>%
    mutate(
      Subgroup = factor(Subgroup,
                        levels = c("Overall","Female","Male","<75",">=75","CKD no","CKD yes")),
      Group = factor(Group, levels = c("nonAKD","Recovery","Non-Recovery"))
    ) %>%
    arrange(Subgroup, Group) %>%
    mutate(
      section = case_when(
        Subgroup == "Overall" ~ "Overall",
        Subgroup %in% c("Female","Male") ~ "Sex subgroup",
        Subgroup %in% c("<75",">=75") ~ "Age subgroup",
        Subgroup %in% c("CKD no","CKD yes") ~ "CKD subgroup",
        TRUE ~ ""
      )
    ) %>%
    group_by(Subgroup) %>%
    mutate(row_in_subgroup = row_number()) %>%
    ungroup() %>%
    group_by(section) %>%
    mutate(row_in_section = row_number()) %>%
    ungroup()
  
  # ==========================
  # 9) Display builders
  # ==========================
  make_display_df <- function(df){
    
    df <- df %>%
      mutate(
        Section_disp = if_else(row_in_section == 1, section, ""),
        Subgroup_disp = case_when(
          section == "Overall" ~ if_else(row_in_subgroup == 1, " ", ""),
          TRUE ~ if_else(row_in_subgroup == 1, as.character(Subgroup), "")
        ),
        Group_disp = case_when(
          Group == "nonAKD"       ~ "non-AKD",
          Group == "Recovery"     ~ "AKD with recovery",
          Group == "Non-Recovery" ~ "AKD without recovery",
          TRUE ~ ""
        ),
        Number_chr = as.character(Number),
        Deaths_chr = paste0(Deaths, " (", Death_pct, ")"),
        Risk_chr   = sprintf("%.1f", 100 * risk5),
        RD_chr = case_when(
          Group == "nonAKD" ~ "",
          is.na(RD) ~ "",
          TRUE ~ sprintf("%+.1f (%+.1f, %+.1f)",
                         RD * 100,
                         RD_low * 100,
                         RD_high * 100)
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
        blank_ci = paste(rep(" ", 18), collapse = " ")
      )
    
    plot_df <- df %>%
      select(
        Section = Section_disp,
        Subgroup = Subgroup_disp,
        Group = Group_disp,
        Number = Number_chr,
        `Deaths (%)` = Deaths_chr,
        `Risk (5y, %)` = Risk_chr,
        `Adjusted RD (5y, pp)` = RD_chr,
        ` ` = blank_ci,
        `HR (95%CI)` = HR_chr,
        `P-value` = P_chr
      )
    
    list(full = df, plot = plot_df)
  }
  
  # ==========================
  # 10) Figure builder
  # ==========================
  build_forest_plot <- function(df_full, df_plot){
    fp <- forestploter::forest(
      data  = df_plot,
      est   = df_full$HR,
      lower = df_full$CI_low,
      upper = df_full$CI_high,
      sizes = 0.55,
      ci_column = 8,
      ref_line  = 1,
      x_trans   = "log",
      xlim      = c(0.5, 12),
      ticks_at  = c(0.5, 1, 2, 4, 8),
      arrow_lab = c("", ""),
      xlab = "Lower risk                                 Higher risk"
    )
    
    fp <- edit_plot(
      fp,
      col = c(4,5,6,7,9,10),
      which = "text",
      hjust = unit(1, "npc"),
      x = unit(1, "npc")
    )
    
    fp <- add_border(fp, row = 0, where = "bottom", gp = gpar(lwd = 1))
    fp$heights <- rep(unit(6.5, "mm"), nrow(fp))
    fp
  }
  
  # ==========================
  # 11) Draw and save
  # ==========================
  disp4 <- make_display_df(forest_full)
  fp4 <- build_forest_plot(disp4$full, disp4$plot)
  
  grid.newpage()
  grid.draw(fp4)
  
  save_pdf_onepage_safe(fp4, out_pdf_4, family = base_family)
  save_tiff_onepage_safe(fp4, out_tiff_4, family = base_family)
  
  cat("\nSaved:\n", out_pdf_4, "\n", out_tiff_4, "\n")
}
#adjusted cox:軽量版(bootstrap部分削除　点推定値のみ)####
{
  ############################################################
  # Figure 4. Subgroup Forest Plot (Sex / Age / CKD)
  # LIGHT VERSION:
  # - RR column removed
  # - Adjusted 5-year RD added by LIGHT G-computation
  # - NO bootstrap CI for RD  -> much faster
  # - HR: adjusted Cox model
  ############################################################
  
  graphics.off()
  
  # ==========================
  # Packages
  # ==========================
  pkgs <- c("readr","dplyr","survival","broom","tibble",
            "forestploter","grid","gridExtra","ragg","stringr")
  to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
  
  library(readr)
  library(dplyr)
  library(survival)
  library(broom)
  library(tibble)
  library(forestploter)
  library(grid)
  library(gridExtra)
  library(ragg)
  library(stringr)
  
  # ==========================
  # Paths
  # ==========================
  setwd("X:/R")
  in_csv <- "jin1_Eligibile.csv"
  outdir <- file.path(getwd(), "figure_table")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  out_pdf_4  <- file.path(outdir, "Figure4_ForestPlot_Subgroups_Sex_Age_CKD_adjRD_LIGHT.pdf")
  out_tiff_4 <- file.path(outdir, "Figure4_ForestPlot_Subgroups_Sex_Age_CKD_adjRD_LIGHT.tiff")
  out_csv_4  <- file.path(outdir, "Figure4_ForestPlot_Subgroups_Sex_Age_CKD_adjRD_LIGHT.csv")
  
  # ==========================
  # Font
  # ==========================
  base_family <- {
    f <- c("Yu Gothic", "MS Gothic", "Meiryo", "Arial Unicode MS", "Arial")
    ok <- f[f %in% names(grDevices::windowsFonts())]
    if (length(ok) == 0) "sans" else ok[1]
  }
  
  # ==========================
  # Parameters
  # ==========================
  rd_horizon <- 5
  
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
                  "ge75" = ">=75")
  }
  
  disp_ckdgrp <- function(x){
    dplyr::recode(as.character(x),
                  "0" = "CKD no",
                  "1" = "CKD yes")
  }
  
  fmt_p <- function(p){
    case_when(
      is.na(p)  ~ "",
      p < 0.001 ~ "<0.001",
      TRUE      ~ sprintf("%.3f", p)
    )
  }
  
  fmt_hr <- function(hr, lo, hi){
    ifelse(
      is.na(hr), "",
      sprintf("%.2f (%.2f to %.2f)", hr, lo, hi)
    )
  }
  
  fmt_risk <- function(x){
    ifelse(is.na(x), "", sprintf("%.1f", 100 * x))
  }
  
  fmt_rd <- function(x){
    ifelse(is.na(x), "", sprintf("%+.1f", 100 * x))
  }
  
  save_pdf_onepage_safe <- function(grob, filename,
                                    page_w_mm = 380,
                                    page_h_mm = 235,
                                    margin_left_mm = 24,
                                    margin_right_mm = 24,
                                    margin_top_mm = 24,
                                    margin_bottom_mm = 42,
                                    scale = 0.58,
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
                                     page_w_mm = 380,
                                     page_h_mm = 235,
                                     margin_left_mm = 24,
                                     margin_right_mm = 24,
                                     margin_top_mm = 24,
                                     margin_bottom_mm = 42,
                                     scale = 0.58,
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
  # LIGHT G-computation for adjusted risk / RD
  # ==========================
  gcomp_rd_light <- function(fit, data, horizon, group_var,
                             ref_group = "nonAKD",
                             target_groups = c("Recovery","Non-Recovery")) {
    
    marginal_risk <- function(fit_obj, dat, grp){
      dat_cf <- dat
      dat_cf[[group_var]] <- factor(grp, levels = levels(dat[[group_var]]))
      
      sf <- survfit(fit_obj, newdata = dat_cf)
      s  <- summary(sf, times = horizon, extend = TRUE)
      
      mean(1 - s$surv, na.rm = TRUE)
    }
    
    risk_ref <- marginal_risk(fit, data, ref_group)
    
    bind_rows(lapply(target_groups, function(g){
      risk_g <- marginal_risk(fit, data, g)
      tibble(
        Group = g,
        risk5_adj = risk_g,
        RD = risk_g - risk_ref
      )
    }))
  }
  
  # ==========================
  # Summary/stat helpers
  # ==========================
  calc_stats_overall <- function(data){
    data %>%
      group_by(group) %>%
      summarise(
        Number = n(),
        Deaths = sum(primary_death),
        Death_pct = sprintf("%.1f", 100 * Deaths / Number),
        .groups = "drop"
      ) %>%
      mutate(Subgroup = "Overall",
             Group = as.character(group)) %>%
      select(Subgroup, Group, Number, Deaths, Death_pct)
  }
  
  calc_stats_sub <- function(data, subgroup_var, subgroup_lab_fun = NULL){
    out <- data %>%
      group_by(!!sym(subgroup_var), group) %>%
      summarise(
        Number = n(),
        Deaths = sum(primary_death),
        Death_pct = sprintf("%.1f", 100 * Deaths / Number),
        .groups = "drop"
      ) %>%
      mutate(
        Subgroup = as.character(!!sym(subgroup_var)),
        Group = as.character(group)
      )
    if (!is.null(subgroup_lab_fun)) out$Subgroup <- subgroup_lab_fun(out$Subgroup)
    out %>% select(Subgroup, Group, Number, Deaths, Death_pct)
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
  # 0) Load
  # ==========================
  loc <- locale(encoding = "SHIFT-JIS")
  jin1_Eligibile <- read_csv(in_csv, locale = loc, show_col_types = FALSE)
  
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
  # 2) Overall model
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
  # 4) Interaction models for HR
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
  # 5) Stats
  # ==========================
  stats_overall <- calc_stats_overall(dat_overall)
  stats_sex <- calc_stats_sub(dat_sex, "sex_f", function(x) dplyr::recode(x, "F" = "Female", "M" = "Male"))
  stats_age <- calc_stats_sub(dat_age, "age75_f", disp_agegrp)
  stats_ckd <- calc_stats_sub(dat_ckd, "ckd_bin", disp_ckdgrp)
  
  # ==========================
  # 6) Adjusted HR
  # ==========================
  hr_overall <- extract_hr_overall(fit_overall)
  
  hr_sex0 <- hr_from_interaction(fit_sex, subgroup_levels = c("Female","Male"), int_suffix_other = "sex_fM")
  hr_sex <- bind_rows(
    tibble(Subgroup = "Female", Group = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
    filter(hr_sex0, Subgroup == "Female"),
    tibble(Subgroup = "Male", Group = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
    filter(hr_sex0, Subgroup == "Male")
  )
  
  hr_age0 <- hr_from_interaction(fit_age, subgroup_levels = c("<75",">=75"), int_suffix_other = "age75_fge75")
  hr_age <- bind_rows(
    tibble(Subgroup = "<75", Group = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
    filter(hr_age0, Subgroup == "<75"),
    tibble(Subgroup = ">=75", Group = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
    filter(hr_age0, Subgroup == ">=75")
  )
  
  hr_ckd0 <- hr_from_interaction(fit_ckd, subgroup_levels = c("CKD no","CKD yes"), int_suffix_other = "ckd_bin")
  hr_ckd <- bind_rows(
    tibble(Subgroup = "CKD no", Group = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
    filter(hr_ckd0, Subgroup == "CKD no"),
    tibble(Subgroup = "CKD yes", Group = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
    filter(hr_ckd0, Subgroup == "CKD yes")
  )
  
  # ==========================
  # 7) Adjusted RD (LIGHT, point estimate only)
  # ==========================
  # Overall
  rd_overall <- gcomp_rd_light(
    fit = fit_overall,
    data = dat_overall,
    horizon = rd_horizon,
    group_var = "group"
  ) %>%
    mutate(Subgroup = "Overall")
  
  # Sex-specific
  fit_female <- coxph(
    Surv(time_years, primary_death) ~
      group + age + index_cre + arb_acei_use +
      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
    data = dat_sex %>% filter(sex_f == "F")
  )
  
  fit_male <- coxph(
    Surv(time_years, primary_death) ~
      group + age + index_cre + arb_acei_use +
      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
    data = dat_sex %>% filter(sex_f == "M")
  )
  
  rd_female <- gcomp_rd_light(
    fit = fit_female,
    data = dat_sex %>% filter(sex_f == "F"),
    horizon = rd_horizon,
    group_var = "group"
  ) %>%
    mutate(Subgroup = "Female")
  
  rd_male <- gcomp_rd_light(
    fit = fit_male,
    data = dat_sex %>% filter(sex_f == "M"),
    horizon = rd_horizon,
    group_var = "group"
  ) %>%
    mutate(Subgroup = "Male")
  
  # Age-specific
  fit_lt75 <- coxph(
    Surv(time_years, primary_death) ~
      group + age + index_cre + arb_acei_use +
      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
    data = dat_age %>% filter(age75_f == "lt75")
  )
  
  fit_ge75 <- coxph(
    Surv(time_years, primary_death) ~
      group + age + index_cre + arb_acei_use +
      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
    data = dat_age %>% filter(age75_f == "ge75")
  )
  
  rd_lt75 <- gcomp_rd_light(
    fit = fit_lt75,
    data = dat_age %>% filter(age75_f == "lt75"),
    horizon = rd_horizon,
    group_var = "group"
  ) %>%
    mutate(Subgroup = "<75")
  
  rd_ge75 <- gcomp_rd_light(
    fit = fit_ge75,
    data = dat_age %>% filter(age75_f == "ge75"),
    horizon = rd_horizon,
    group_var = "group"
  ) %>%
    mutate(Subgroup = ">=75")
  
  # CKD-specific
  fit_ckd_no <- coxph(
    Surv(time_years, primary_death) ~
      group + age + index_cre + arb_acei_use +
      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
    data = dat_ckd %>% filter(ckd_bin == 0)
  )
  
  fit_ckd_yes <- coxph(
    Surv(time_years, primary_death) ~
      group + age + index_cre + arb_acei_use +
      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
    data = dat_ckd %>% filter(ckd_bin == 1)
  )
  
  rd_ckd_no <- gcomp_rd_light(
    fit = fit_ckd_no,
    data = dat_ckd %>% filter(ckd_bin == 0),
    horizon = rd_horizon,
    group_var = "group"
  ) %>%
    mutate(Subgroup = "CKD no")
  
  rd_ckd_yes <- gcomp_rd_light(
    fit = fit_ckd_yes,
    data = dat_ckd %>% filter(ckd_bin == 1),
    horizon = rd_horizon,
    group_var = "group"
  ) %>%
    mutate(Subgroup = "CKD yes")
  
  rd_all <- bind_rows(
    rd_overall,
    rd_female,
    rd_male,
    rd_lt75,
    rd_ge75,
    rd_ckd_no,
    rd_ckd_yes
  ) %>%
    select(Subgroup, Group, risk5_adj, RD)
  
  # nonAKD row を追加
  rd_ref_rows <- tibble(
    Subgroup = c("Overall","Female","Male","<75",">=75","CKD no","CKD yes"),
    Group = "nonAKD",
    risk5_adj = NA_real_,
    RD = 0
  )
  
  rd_all2 <- bind_rows(rd_ref_rows, rd_all)
  # ==========================
  # 8) Merge all blocks
  # ==========================
  plotdat <- bind_rows(
    stats_overall,
    stats_sex,
    stats_age,
    stats_ckd
  ) %>%
    left_join(
      bind_rows(hr_overall, hr_sex, hr_age, hr_ckd),
      by = c("Subgroup", "Group")
    ) %>%
    left_join(
      rd_all2,
      by = c("Subgroup", "Group")
    ) %>%
    mutate(
      GroupLabel = dplyr::recode(
        as.character(Group),
        "nonAKD"       = "non-AKD",
        "Recovery"     = "AKD with recovery",
        "Non-Recovery" = "AKD without recovery"
      ),
      Risk5y_adj = fmt_risk(risk5_adj),
      RD5y_adj   = fmt_rd(RD),
      HR_CI      = fmt_hr(HR, CI_low, CI_high),
      P_txt      = fmt_p(P_value)
    )
  
  subgroup_order <- c("Overall", "Female", "Male", "<75", ">=75", "CKD no", "CKD yes")
  group_order    <- c("nonAKD", "Recovery", "Non-Recovery")
  
  plotdat <- plotdat %>%
    mutate(
      Subgroup = factor(Subgroup, levels = subgroup_order),
      Group    = factor(Group, levels = group_order)
    ) %>%
    arrange(Subgroup, Group)
  
  # ==========================
  # 9) Section / Subgroup / Group を分けて表示
  # ==========================
  plotdat2 <- plotdat %>%
    mutate(
      Section = case_when(
        Subgroup == "Overall" ~ "Overall",
        Subgroup %in% c("Female", "Male") ~ "Sex subgroup",
        Subgroup %in% c("<75", ">=75") ~ "Age subgroup",
        Subgroup %in% c("CKD no", "CKD yes") ~ "CKD subgroup",
        TRUE ~ ""
      )
    ) %>%
    group_by(Section) %>%
    mutate(
      Section_disp = if_else(row_number() == 1, Section, ""),
      Subgroup_disp = case_when(
        Section == "Overall" ~ if_else(row_number() == 1, "Overall", ""),
        TRUE ~ if_else(Group == "nonAKD", as.character(Subgroup), "")
      ),
      Group_disp  = GroupLabel,
      Number_disp = as.character(Number),
      Deaths_disp = ifelse(is.na(Deaths), "", sprintf("%d (%.1f)", Deaths, as.numeric(Death_pct))),
      HR_CI_disp  = if_else(Group == "nonAKD", "Reference", HR_CI),
      P_txt_disp  = if_else(Group == "nonAKD", "", P_txt)
    ) %>%
    ungroup()
  
  # forest positions
  plotdat2$est   <- plotdat2$HR
  plotdat2$lower <- plotdat2$CI_low
  plotdat2$upper <- plotdat2$CI_high
  
  # ==========================
  # 10) Forest plot table text
  # 1 Section
  # 2 Subgroup
  # 3 Group
  # 4 Number
  # 5 Deaths
  # 6 Risk
  # 7 RD
  # 8 spacer
  # 9 Forest
  # 10 HR
  # 11 P
  # ==========================
  tabletext <- plotdat2 %>%
    transmute(
      Section = Section_disp,
      Subgroup = Subgroup_disp,
      Group = Group_disp,
      Number = Number_disp,
      `Deaths (%)` = Deaths_disp,
      `Risk (5y, %)` = Risk5y_adj,
      `RD (5y, pp)` = RD5y_adj,
      ` ` = "",                  # spacer between RD and forest
      `  ` = "",                 # forest column
      `HR (95% CI)` = HR_CI_disp,
      `P-value` = P_txt_disp
    )
  
  # ==========================
  # 11) Draw forest plot
  # ==========================
  theme1 <- forest_theme(
    base_size = 10.0,
    base_family = base_family,
    refline_col = "grey40",
    footnote_col = "grey30",
    ci_pch = 15,
    ci_lwd = 1.8,
    ci_Theight = 0.30,
    summary_col = "black",
    core = list(
      padding = unit(c(0.8, 0.8), "mm"),
      fg_params = list(
        hjust = c(0, 0, 0, 1, 1, 1, 1, 0.5, 0.5, 1, 1),
        x     = c(0, 0, 0, 1, 1, 1, 1, 0.5, 0.5, 1, 1),
        fontsize = 10.0
      ),
      bg_params = list(fill = c("#F7F7F7", "white"))
    ),
    colhead = list(
      fg_params = list(
        fontface = "bold",
        fontsize = 10.2,
        hjust = c(0, 0, 0, 1, 1, 1, 1, 0.5, 0.5, 1, 1),
        x     = c(0, 0, 0, 1, 1, 1, 1, 0.5, 0.5, 1, 1)
      )
    )
  )
  
  p <- forest(
    tabletext,
    est = plotdat2$est,
    lower = plotdat2$lower,
    upper = plotdat2$upper,
    ci_column = 9,                 # forest column
    ref_line = 1,
    xlim = c(0.5, 8.0),
    ticks_at = c(0.5, 1, 2, 4, 8),
    arrow_lab = c("Lower risk", "Higher risk"),
    xlab = "Adjusted Hazard Ratio",
    theme = theme1
  )
  
  # --------------------------
  # 本文の列位置
  # --------------------------
  p <- edit_plot(p, col = 1,  which = "text", x = unit(0, "npc"),    hjust = 0) # Section
  p <- edit_plot(p, col = 2,  which = "text", x = unit(0.06, "npc"), hjust = 0) # Subgroup ← 少し右
  p <- edit_plot(p, col = 3,  which = "text", x = unit(0, "npc"),    hjust = 0) # Group
  p <- edit_plot(p, col = 4,  which = "text", x = unit(1, "npc"),    hjust = 1) # Number
  p <- edit_plot(p, col = 5,  which = "text", x = unit(1, "npc"),    hjust = 1) # Deaths
  p <- edit_plot(p, col = 6,  which = "text", x = unit(1, "npc"),    hjust = 1) # Risk
  p <- edit_plot(p, col = 7,  which = "text", x = unit(0.72, "npc"), hjust = 1) # RD ← 左へ
  p <- edit_plot(p, col = 10, which = "text", x = unit(1, "npc"),    hjust = 1) # HR
  p <- edit_plot(p, col = 11, which = "text", x = unit(1, "npc"),    hjust = 1) # P
  
  # --------------------------
  # 見出し位置
  # --------------------------
  p <- edit_plot(p, row = 0, col = 1,  which = "text", x = unit(0, "npc"),    hjust = 0)
  p <- edit_plot(p, row = 0, col = 2,  which = "text", x = unit(0.06, "npc"), hjust = 0)
  p <- edit_plot(p, row = 0, col = 3,  which = "text", x = unit(0, "npc"),    hjust = 0)
  p <- edit_plot(p, row = 0, col = 4,  which = "text", x = unit(1, "npc"),    hjust = 1)
  p <- edit_plot(p, row = 0, col = 5,  which = "text", x = unit(1, "npc"),    hjust = 1)
  p <- edit_plot(p, row = 0, col = 6,  which = "text", x = unit(1, "npc"),    hjust = 1)
  p <- edit_plot(p, row = 0, col = 7,  which = "text", x = unit(0.72, "npc"), hjust = 1)
  p <- edit_plot(p, row = 0, col = 10, which = "text", x = unit(1, "npc"),    hjust = 1)
  p <- edit_plot(p, row = 0, col = 11, which = "text", x = unit(1, "npc"),    hjust = 1)
  
  # ==========================
  # 12) 列幅
  # ==========================
  p$widths[1]  <- unit(16, "mm")   # Section
  p$widths[2]  <- unit(24, "mm")   # Subgroup
  p$widths[3]  <- unit(34, "mm")   # Group
  p$widths[4]  <- unit(18, "mm")   # Number
  p$widths[5]  <- unit(28, "mm")   # Deaths
  p$widths[6]  <- unit(30, "mm")   # Risk ← 広げる
  p$widths[7]  <- unit(30, "mm")   # RD   ← 広げる
  p$widths[8]  <- unit(10, "mm")   # spacer between RD and forest
  p$widths[9]  <- unit(110, "mm")  # Forest ← かなり広げる
  p$widths[10] <- unit(36, "mm")   # HR
  p$widths[11] <- unit(14, "mm")   # P
  
  # --------------------------
  # 下余白
  # --------------------------
  if ("xlab-b" %in% p$layout$name) {
    idx <- which(p$layout$name == "xlab-b")
    p$heights[idx] <- unit(3.0, "mm")
  }
  if ("axis-b" %in% p$layout$name) {
    idx <- which(p$layout$name == "axis-b")
    p$heights[idx] <- unit(4.0, "mm")
  }
  
  # ==========================
  # 13) Save
  # scale を無理に上げず、ページ自体を広げる
  # ==========================
  save_pdf_onepage_safe(
    grob = p,
    filename = out_pdf_4,
    page_w_mm = 540,
    page_h_mm = 235,
    margin_left_mm = 0.2,
    margin_right_mm = 0.2,
    margin_top_mm = 0.2,
    margin_bottom_mm = 0.8,
    scale = 1.00,
    family = base_family
  )
  
  save_tiff_onepage_safe(
    grob = p,
    filename = out_tiff_4,
    page_w_mm = 540,
    page_h_mm = 235,
    margin_left_mm = 0.2,
    margin_right_mm = 0.2,
    margin_top_mm = 0.2,
    margin_bottom_mm = 0.8,
    scale = 1.00,
    family = base_family
  )
}
#adjusted cox:ハイブリット版→rd-stage1で2時間たっても実行できず####
{
  ############################################################
  # Figure 4. Subgroup Forest Plot (Sex / Age / CKD)
  # Adjusted RD at 5 years by G-computation with staged bootstrap
  #
  # rd_stage = 1 :
  #   - Overall only bootstrap 50
  #   - subgroup = point estimate only
  #
  # rd_stage = 2 :
  #   - Overall bootstrap 50
  #   - subgroup bootstrap 30
  #
  # rd_stage = 3 :
  #   - Overall/subgroup bootstrap 100
  ############################################################
  
  graphics.off()
  
  # ==========================
  # Packages
  # ==========================
  pkgs <- c(
    "readr", "dplyr", "survival", "broom", "tibble",
    "forestploter", "grid", "gridExtra", "ragg"
  )
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
  
  out_pdf_4  <- file.path(outdir, "Figure4_ForestPlot_Subgroups_Sex_Age_CKD_adjRD_stage.pdf")
  out_tiff_4 <- file.path(outdir, "Figure4_ForestPlot_Subgroups_Sex_Age_CKD_adjRD_stage.tiff")
  out_csv_4  <- file.path(outdir, "Figure4_ForestPlot_Subgroups_Sex_Age_CKD_adjRD_stage.csv")
  
  # ==========================
  # Font
  # ==========================
  base_family <- {
    f <- c("Yu Gothic", "MS Gothic", "Meiryo", "Arial Unicode MS", "Arial")
    ok <- f[f %in% names(grDevices::windowsFonts())]
    if (length(ok) == 0) "sans" else ok[1]
  }
  
  # ==========================
  # Parameters
  # ==========================
  rd_horizon <- 5
  
  # 1 = Overall only bootstrap 50; subgroup point estimate only
  # 2 = Overall 50 + subgroup 30
  # 3 = Final version: overall/subgroup 100
  rd_stage <- 1
  
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
                  "ge75" = ">=75")
  }
  
  disp_ckdgrp <- function(x){
    dplyr::recode(as.character(x),
                  "0" = "CKD no",
                  "1" = "CKD yes")
  }
  
  fmt_p <- function(p){
    dplyr::case_when(
      is.na(p)  ~ "",
      p < 0.001 ~ "<0.001",
      p < 0.01  ~ sprintf("%.3f", p),
      TRUE      ~ sprintf("%.2f", p)
    )
  }
  
  fmt_hr <- function(hr, lo, hi){
    ifelse(
      is.na(hr), "",
      sprintf("%.2f (%.2f–%.2f)", hr, lo, hi)
    )
  }
  
  fmt_risk <- function(x){
    ifelse(is.na(x), "", sprintf("%.1f", 100 * x))
  }
  
  fmt_rd <- function(x){
    ifelse(is.na(x), "", sprintf("%+.1f", 100 * x))
  }
  
  save_pdf_onepage_safe <- function(grob, filename,
                                    page_w_mm = 500,
                                    page_h_mm = 246,
                                    margin_left_mm = 2,
                                    margin_right_mm = 2,
                                    margin_top_mm = 2,
                                    margin_bottom_mm = 4,
                                    scale = 1.00,
                                    family = "sans"){
    pdf(
      filename,
      width = page_w_mm / 25.4,
      height = page_h_mm / 25.4,
      family = family
    )
    
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(
      x = grid::unit(margin_left_mm, "mm"),
      y = grid::unit(margin_bottom_mm, "mm"),
      just = c("left", "bottom"),
      width = grid::unit(1, "npc") - grid::unit(margin_left_mm + margin_right_mm, "mm"),
      height = grid::unit(1, "npc") - grid::unit(margin_top_mm + margin_bottom_mm, "mm"),
      clip = "on"
    ))
    grid::pushViewport(grid::viewport(
      x = 0.5, y = 0.5,
      width = grid::unit(scale, "npc"),
      height = grid::unit(scale, "npc"),
      just = c("center", "center"),
      clip = "on"
    ))
    grid::grid.draw(grob)
    grid::popViewport(2)
    dev.off()
    cat("Saved PDF:", normalizePath(filename), "\n")
  }
  
  save_tiff_onepage_safe <- function(grob, filename,
                                     page_w_mm = 500,
                                     page_h_mm = 246,
                                     margin_left_mm = 2,
                                     margin_right_mm = 2,
                                     margin_top_mm = 2,
                                     margin_bottom_mm = 4,
                                     scale = 1.00,
                                     dpi = 600,
                                     family = "sans"){
    ragg::agg_tiff(
      filename,
      width = page_w_mm / 25.4,
      height = page_h_mm / 25.4,
      units = "in",
      res = dpi,
      compression = "lzw"
    )
    
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(
      x = grid::unit(margin_left_mm, "mm"),
      y = grid::unit(margin_bottom_mm, "mm"),
      just = c("left", "bottom"),
      width = grid::unit(1, "npc") - grid::unit(margin_left_mm + margin_right_mm, "mm"),
      height = grid::unit(1, "npc") - grid::unit(margin_top_mm + margin_bottom_mm, "mm"),
      clip = "on"
    ))
    grid::pushViewport(grid::viewport(
      x = 0.5, y = 0.5,
      width = grid::unit(scale, "npc"),
      height = grid::unit(scale, "npc"),
      just = c("center", "center"),
      clip = "on"
    ))
    grid::grid.draw(grob)
    grid::popViewport(2)
    dev.off()
    cat("Saved TIFF:", normalizePath(filename), "\n")
  }
  
  # ==========================
  # RD bootstrap plan by stage
  # ==========================
  get_rd_boot_n <- function(stage, subgroup_name){
    if (stage == 1) {
      if (subgroup_name == "Overall") return(50)
      return(0)
    }
    if (stage == 2) {
      if (subgroup_name == "Overall") return(50)
      return(30)
    }
    if (stage == 3) {
      return(100)
    }
    stop("rd_stage must be 1, 2, or 3.")
  }
  
  # ==========================
  # G-computation for adjusted RD
  # - bootstrap optional
  # - robust to failed refits
  # ==========================
  gcomp_rd <- function(fit, data, horizon, group_var,
                       ref_group = "nonAKD",
                       target_groups = c("Recovery","Non-Recovery"),
                       n_boot = 0,
                       seed = 1234) {
    
    marginal_risk <- function(fit_obj, dat, grp){
      dat_cf <- dat
      dat_cf[[group_var]] <- factor(grp, levels = levels(dat[[group_var]]))
      
      sf <- survfit(fit_obj, newdata = dat_cf)
      s  <- summary(sf, times = horizon, extend = TRUE)
      
      mean(1 - s$surv, na.rm = TRUE)
    }
    
    # point estimate
    risk_ref <- marginal_risk(fit, data, ref_group)
    
    rd_point <- sapply(target_groups, function(g){
      marginal_risk(fit, data, g) - risk_ref
    })
    
    out <- tibble(
      Group = target_groups,
      RD = as.numeric(rd_point),
      RD_low = NA_real_,
      RD_high = NA_real_,
      boot_n_success = NA_integer_
    )
    
    if (is.null(n_boot) || n_boot <= 0) {
      return(out)
    }
    
    set.seed(seed)
    boot_list <- vector("list", n_boot)
    
    for (b in seq_len(n_boot)) {
      idx <- sample.int(nrow(data), size = nrow(data), replace = TRUE)
      boot_dat <- data[idx, , drop = FALSE]
      
      boot_dat[[group_var]] <- factor(boot_dat[[group_var]], levels = levels(data[[group_var]]))
      
      boot_fit <- try(update(fit, data = boot_dat), silent = TRUE)
      if (inherits(boot_fit, "try-error")) {
        boot_list[[b]] <- rep(NA_real_, length(target_groups))
        next
      }
      
      risk_ref_b <- try(marginal_risk(boot_fit, boot_dat, ref_group), silent = TRUE)
      if (inherits(risk_ref_b, "try-error") || is.na(risk_ref_b)) {
        boot_list[[b]] <- rep(NA_real_, length(target_groups))
        next
      }
      
      vals <- sapply(target_groups, function(g){
        tmp <- try(marginal_risk(boot_fit, boot_dat, g), silent = TRUE)
        if (inherits(tmp, "try-error") || is.na(tmp)) return(NA_real_)
        tmp - risk_ref_b
      })
      
      boot_list[[b]] <- vals
    }
    
    rd_boot <- do.call(cbind, boot_list)
    
    if (is.null(rd_boot) || ncol(rd_boot) == 0) {
      return(out)
    }
    
    success_vec <- apply(rd_boot, 2, function(x) all(is.finite(x)))
    n_success <- sum(success_vec)
    
    if (n_success < 5) {
      warning("Too few successful bootstrap replicates: ", n_success)
      out$boot_n_success <- n_success
      return(out)
    }
    
    rd_boot_ok <- rd_boot[, success_vec, drop = FALSE]
    
    out$RD_low <- apply(rd_boot_ok, 1, quantile, probs = 0.025, na.rm = TRUE)
    out$RD_high <- apply(rd_boot_ok, 1, quantile, probs = 0.975, na.rm = TRUE)
    out$boot_n_success <- n_success
    
    out
  }
  
  # ==========================
  # 0) Load
  # ==========================
  loc <- locale(encoding = "SHIFT-JIS")
  jin1_Eligibile <- read_csv(in_csv, locale = loc, show_col_types = FALSE)
  
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
        jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` %in% c(0, 2) ~ "Non-Recovery",
        TRUE ~ NA_character_
      ),
      group = factor(group, levels = c("nonAKD", "Recovery", "Non-Recovery")),
      arb_acei_use = if_else(coalesce(arb, 0) == 1 | coalesce(acei, 0) == 1, 1L, 0L),
      time_years   = as.numeric(last_follow_death - index_plus_210) / 365.25
    ) %>%
    filter(
      !is.na(group),
      !is.na(time_years), time_years >= 0,
      !is.na(primary_death),
      !is.na(age),
      !is.na(index_cre)
    )
  
  # ==========================
  # 2) Overall model
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
  # 4) Interaction models for HR
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
        .groups = "drop"
      ) %>%
      mutate(
        Subgroup = "Overall",
        Group = as.character(group)
      ) %>%
      select(Subgroup, Group, Number, Deaths, Death_pct)
  }
  
  calc_stats_sub <- function(data, subgroup_var, subgroup_lab_fun = NULL){
    out <- data %>%
      group_by(!!sym(subgroup_var), group) %>%
      summarise(
        Number = n(),
        Deaths = sum(primary_death),
        Death_pct = sprintf("%.1f", 100 * Deaths / Number),
        .groups = "drop"
      ) %>%
      mutate(
        Subgroup = as.character(!!sym(subgroup_var)),
        Group = as.character(group)
      )
    
    if (!is.null(subgroup_lab_fun)) out$Subgroup <- subgroup_lab_fun(out$Subgroup)
    out %>% select(Subgroup, Group, Number, Deaths, Death_pct)
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
  # 6) Summary + HR
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
  
  hr_age0 <- hr_from_interaction(fit_age, subgroup_levels = c("<75",">=75"), int_suffix_other = "age75_fge75")
  hr_age <- bind_rows(
    tibble(Subgroup = "<75", Group = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
    filter(hr_age0, Subgroup == "<75"),
    tibble(Subgroup = ">=75", Group = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
    filter(hr_age0, Subgroup == ">=75")
  )
  
  hr_ckd0 <- hr_from_interaction(fit_ckd, subgroup_levels = c("CKD no","CKD yes"), int_suffix_other = "ckd_bin")
  hr_ckd <- bind_rows(
    tibble(Subgroup = "CKD no", Group = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
    filter(hr_ckd0, Subgroup == "CKD no"),
    tibble(Subgroup = "CKD yes", Group = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
    filter(hr_ckd0, Subgroup == "CKD yes")
  )
  
  # ==========================
  # 7) Adjusted RD by subgroup (stage-based)
  # ==========================
  rd_overall <- gcomp_rd(
    fit = fit_overall,
    data = dat_overall,
    horizon = rd_horizon,
    group_var = "group",
    n_boot = get_rd_boot_n(rd_stage, "Overall"),
    seed = 1001
  ) %>%
    mutate(Subgroup = "Overall")
  
  fit_female <- coxph(
    Surv(time_years, primary_death) ~
      group + age + index_cre + arb_acei_use +
      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
    data = dat_sex %>% filter(sex_f == "F")
  )
  
  rd_female <- gcomp_rd(
    fit = fit_female,
    data = dat_sex %>% filter(sex_f == "F"),
    horizon = rd_horizon,
    group_var = "group",
    n_boot = get_rd_boot_n(rd_stage, "Female"),
    seed = 1002
  ) %>%
    mutate(Subgroup = "Female")
  
  fit_male <- coxph(
    Surv(time_years, primary_death) ~
      group + age + index_cre + arb_acei_use +
      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
    data = dat_sex %>% filter(sex_f == "M")
  )
  
  rd_male <- gcomp_rd(
    fit = fit_male,
    data = dat_sex %>% filter(sex_f == "M"),
    horizon = rd_horizon,
    group_var = "group",
    n_boot = get_rd_boot_n(rd_stage, "Male"),
    seed = 1003
  ) %>%
    mutate(Subgroup = "Male")
  
  fit_lt75 <- coxph(
    Surv(time_years, primary_death) ~
      group + age + index_cre + arb_acei_use +
      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
    data = dat_age %>% filter(age75_f == "lt75")
  )
  
  rd_lt75 <- gcomp_rd(
    fit = fit_lt75,
    data = dat_age %>% filter(age75_f == "lt75"),
    horizon = rd_horizon,
    group_var = "group",
    n_boot = get_rd_boot_n(rd_stage, "<75"),
    seed = 1004
  ) %>%
    mutate(Subgroup = "<75")
  
  fit_ge75 <- coxph(
    Surv(time_years, primary_death) ~
      group + age + index_cre + arb_acei_use +
      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
    data = dat_age %>% filter(age75_f == "ge75")
  )
  
  rd_ge75 <- gcomp_rd(
    fit = fit_ge75,
    data = dat_age %>% filter(age75_f == "ge75"),
    horizon = rd_horizon,
    group_var = "group",
    n_boot = get_rd_boot_n(rd_stage, ">=75"),
    seed = 1005
  ) %>%
    mutate(Subgroup = ">=75")
  
  fit_ckd_no <- coxph(
    Surv(time_years, primary_death) ~
      group + age + index_cre + arb_acei_use +
      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
    data = dat_ckd %>% filter(ckd_bin == 0)
  )
  
  rd_ckd_no <- gcomp_rd(
    fit = fit_ckd_no,
    data = dat_ckd %>% filter(ckd_bin == 0),
    horizon = rd_horizon,
    group_var = "group",
    n_boot = get_rd_boot_n(rd_stage, "CKD no"),
    seed = 1006
  ) %>%
    mutate(Subgroup = "CKD no")
  
  fit_ckd_yes <- coxph(
    Surv(time_years, primary_death) ~
      group + age + index_cre + arb_acei_use +
      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
    data = dat_ckd %>% filter(ckd_bin == 1)
  )
  
  rd_ckd_yes <- gcomp_rd(
    fit = fit_ckd_yes,
    data = dat_ckd %>% filter(ckd_bin == 1),
    horizon = rd_horizon,
    group_var = "group",
    n_boot = get_rd_boot_n(rd_stage, "CKD yes"),
    seed = 1007
  ) %>%
    mutate(Subgroup = "CKD yes")
  
  rd_all <- bind_rows(
    rd_overall,
    rd_female,
    rd_male,
    rd_lt75,
    rd_ge75,
    rd_ckd_no,
    rd_ckd_yes
  )
  
  rd_ref <- tibble(
    Subgroup = c("Overall","Female","Male","<75",">=75","CKD no","CKD yes"),
    Group = "nonAKD",
    RD = NA_real_,
    RD_low = NA_real_,
    RD_high = NA_real_,
    boot_n_success = NA_integer_
  )
  
  rd_all <- bind_rows(rd_ref, rd_all)
  
  # ==========================
  # 8) Merge
  # ==========================
  build_block <- function(stats_df, hr_df, rd_df){
    stats_df %>%
      left_join(hr_df, by = c("Subgroup", "Group")) %>%
      left_join(rd_df, by = c("Subgroup", "Group"))
  }
  
  forest_full <- bind_rows(
    build_block(stats_overall, hr_overall, filter(rd_all, Subgroup == "Overall")),
    build_block(stats_sex,     hr_sex,     filter(rd_all, Subgroup %in% c("Female","Male"))),
    build_block(stats_age,     hr_age,     filter(rd_all, Subgroup %in% c("<75",">=75"))),
    build_block(stats_ckd,     hr_ckd,     filter(rd_all, Subgroup %in% c("CKD no","CKD yes")))
  ) %>%
    mutate(
      Subgroup = factor(Subgroup, levels = c("Overall","Female","Male","<75",">=75","CKD no","CKD yes")),
      Group    = factor(Group, levels = c("nonAKD","Recovery","Non-Recovery"))
    ) %>%
    arrange(Subgroup, Group) %>%
    mutate(
      section = case_when(
        Subgroup == "Overall" ~ "Overall",
        Subgroup %in% c("Female","Male") ~ "Sex subgroup",
        Subgroup %in% c("<75",">=75") ~ "Age subgroup",
        Subgroup %in% c("CKD no","CKD yes") ~ "CKD subgroup",
        TRUE ~ ""
      )
    ) %>%
    group_by(Subgroup) %>%
    mutate(row_in_subgroup = row_number()) %>%
    ungroup() %>%
    group_by(section) %>%
    mutate(row_in_section = row_number()) %>%
    ungroup()
  
  # ==========================
  # 9) Display builders
  # ==========================
  make_display_df <- function(df){
    
    df <- df %>%
      mutate(
        Section_disp = if_else(row_in_section == 1, section, ""),
        Subgroup_disp = case_when(
          section == "Overall" ~ if_else(row_in_subgroup == 1, "Overall", ""),
          TRUE ~ if_else(row_in_subgroup == 1, as.character(Subgroup), "")
        ),
        Group_disp = case_when(
          Group == "nonAKD"       ~ "non-AKD",
          Group == "Recovery"     ~ "AKD with recovery",
          Group == "Non-Recovery" ~ "AKD without recovery",
          TRUE ~ ""
        ),
        Number_chr = as.character(Number),
        Deaths_chr = paste0(Deaths, " (", Death_pct, ")"),
        Risk_chr   = "",
        RD_chr = case_when(
          Group == "nonAKD" ~ "",
          is.na(RD) ~ "",
          is.na(RD_low) | is.na(RD_high) ~ sprintf("%+.1f", RD * 100),
          TRUE ~ sprintf("%+.1f (%+.1f, %+.1f)",
                         RD * 100,
                         RD_low * 100,
                         RD_high * 100)
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
        blank_ci = paste(rep(" ", 18), collapse = " ")
      )
    
    plot_df <- df %>%
      select(
        Section = Section_disp,
        Subgroup = Subgroup_disp,
        Group = Group_disp,
        Number = Number_chr,
        `Deaths (%)` = Deaths_chr,
        `Adjusted RD (5y, pp)` = RD_chr,
        ` ` = blank_ci,
        `HR (95%CI)` = HR_chr,
        `P-value` = P_chr
      )
    
    list(full = df, plot = plot_df)
  }
  
  # ==========================
  # 10) Figure builder
  # ==========================
  build_forest_plot <- function(df_full, df_plot){
    fp <- forestploter::forest(
      data  = df_plot,
      est   = df_full$HR,
      lower = df_full$CI_low,
      upper = df_full$CI_high,
      sizes = 0.65,
      ci_column = 7,
      ref_line  = 1,
      x_trans   = "log",
      xlim      = c(0.5, 8),
      ticks_at  = c(0.5, 1, 2, 4, 8),
      arrow_lab = c("Lower risk", "Higher risk"),
      xlab = "Adjusted Hazard Ratio"
    )
    
    fp <- edit_plot(
      fp,
      col = c(4,5,6,8,9),
      which = "text",
      hjust = unit(1, "npc"),
      x = unit(1, "npc")
    )
    
    fp <- add_border(fp, row = 0, where = "bottom", gp = gpar(lwd = 1))
    fp$heights <- rep(unit(6.8, "mm"), nrow(fp))
    fp
  }
  
  # ==========================
  # 11) Draw and save
  # ==========================
  disp4 <- make_display_df(forest_full)
  fp4 <- build_forest_plot(disp4$full, disp4$plot)
  
  # CSV export
  write_csv(
    forest_full %>%
      mutate(
        Group = as.character(Group),
        Subgroup = as.character(Subgroup)
      ),
    out_csv_4
  )
  
  grid.newpage()
  grid.draw(fp4)
  
  save_pdf_onepage_safe(
    fp4,
    out_pdf_4,
    page_w_mm = 430,
    page_h_mm = 240,
    margin_left_mm = 1,
    margin_right_mm = 1,
    margin_top_mm = 1,
    margin_bottom_mm = 2,
    scale = 1.00,
    family = base_family
  )
  
  save_tiff_onepage_safe(
    fp4,
    out_tiff_4,
    page_w_mm = 430,
    page_h_mm = 240,
    margin_left_mm = 1,
    margin_right_mm = 1,
    margin_top_mm = 1,
    margin_bottom_mm = 2,
    scale = 1.00,
    family = base_family
  )
  
  cat("\nSaved:\n", out_pdf_4, "\n", out_tiff_4, "\n", out_csv_4, "\n")
  cat("\nBootstrap stage =", rd_stage, "\n")
  print(rd_all %>% select(Subgroup, Group, boot_n_success))
}


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
