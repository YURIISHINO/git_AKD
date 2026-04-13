#adjusted cox:おもい####
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
  n_boot_rd  <- 10    # 最終版
  
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
#adjusted cox:軽量版####
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
#若年と高齢で二つのグラフ####
{
  ############################################################
  # Figure 5: Continuous age × AKD interaction
  #  - Main figure: age 40-85 years
  #  - Supplementary figure: full age range
  #  - Grouping: 150_210recovery priority, then 90_150recovery
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
  library(grid)
  
  # ==========================
  # Paths
  # ==========================
  setwd("X:/R")
  in_csv <- "jin1_Eligibile.csv"
  outdir <- file.path(getwd(), "figure_table")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  out_pdf_main  <- file.path(outdir, "Figure5_main_40to85.pdf")
  out_tiff_main <- file.path(outdir, "Figure5_main_40to85.tif")
  out_pdf_full  <- file.path(outdir, "Supplementary_Figure5_fullrange.pdf")
  out_tiff_full <- file.path(outdir, "Supplementary_Figure5_fullrange.tif")
  
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
  
  # 確認
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
  # Reference = non-AKD
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
  # 5) Function to compute HR(age) vs non-AKD
  # ==========================
  get_hr_curve <- function(group_label, ages, age_mean, cf, vc) {
    
    coef_main <- paste0("group", group_label)
    coef_int  <- paste0("group", group_label, ":age_c")
    
    beta_main <- unname(cf[coef_main])
    beta_int  <- unname(cf[coef_int])
    
    var_main <- vc[coef_main, coef_main]
    var_int  <- vc[coef_int,  coef_int]
    cov_mi   <- vc[coef_main, coef_int]
    
    age_c_vals <- ages - age_mean
    log_hr <- beta_main + beta_int * age_c_vals
    se_log_hr <- sqrt(pmax(0, var_main + (age_c_vals^2) * var_int + 2 * age_c_vals * cov_mi))
    
    tibble(
      age   = ages,
      HR    = exp(log_hr),
      lower = exp(log_hr - 1.96 * se_log_hr),
      upper = exp(log_hr + 1.96 * se_log_hr),
      Group = dplyr::case_when(
        group_label == "AKD with recovery"    ~ "AKD with recovery",
        group_label == "AKD without recovery" ~ "AKD without recovery"
      )
    )
  }
  
  # ==========================
  # 6) Create curves
  # full range: use observed 1st-99th percentile
  # main figure: 40-85 years
  # ==========================
  age_min_full <- floor(quantile(dat$age, 0.01, na.rm = TRUE))
  age_max_full <- ceiling(quantile(dat$age, 0.99, na.rm = TRUE))
  
  age_seq_full <- seq(age_min_full, age_max_full, by = 1)
  age_seq_main <- seq(40, 85, by = 1)
  
  curve_full <- bind_rows(
    get_hr_curve("AKD with recovery", age_seq_full, age_mean, cf, vc),
    get_hr_curve("AKD without recovery", age_seq_full, age_mean, cf, vc)
  )
  
  curve_main <- bind_rows(
    get_hr_curve("AKD with recovery", age_seq_main, age_mean, cf, vc),
    get_hr_curve("AKD without recovery", age_seq_main, age_mean, cf, vc)
  )
  
  # ==========================
  # 7) Truncation settings
  # ==========================
  # Main figure: emphasize convergence near HR = 1
  ymin_main <- 0.5
  ymax_main <- 3.0
  
  # Full supplementary figure: allow wider range
  ymin_full <- max(0.3, min(curve_full$lower, 1, na.rm = TRUE) * 0.95)
  ymax_full <- max(curve_full$upper, 1, na.rm = TRUE) * 1.05
  
  prepare_plot <- function(dat, ymin, ymax) {
    dat %>%
      mutate(
        lower_plot = pmax(lower, ymin),
        upper_plot = pmin(upper, ymax),
        lo_trunc   = lower < ymin,
        hi_trunc   = upper > ymax
      )
  }
  
  plot_main <- prepare_plot(curve_main, ymin_main, ymax_main)
  plot_full <- prepare_plot(curve_full, ymin_full, ymax_full)
  
  # ==========================
  # 8) Arrow data (safe version)
  # ==========================
  arrow_main_upper <- plot_main %>%
    filter(hi_trunc) %>%
    group_by(Group) %>%
    filter((row_number() - 1) %% 4 == 0) %>%
    ungroup()
  
  arrow_main_lower <- plot_main %>%
    filter(lo_trunc) %>%
    group_by(Group) %>%
    filter((row_number() - 1) %% 4 == 0) %>%
    ungroup()
  
  arrow_full_upper <- plot_full %>%
    filter(hi_trunc) %>%
    group_by(Group) %>%
    filter((row_number() - 1) %% 4 == 0) %>%
    ungroup()
  
  arrow_full_lower <- plot_full %>%
    filter(lo_trunc) %>%
    group_by(Group) %>%
    filter((row_number() - 1) %% 4 == 0) %>%
    ungroup()
  
  # ==========================
  # 9) Plot function
  # ==========================
  make_plot <- function(dat, arrow_up, arrow_low, ymin, ymax, xlim = NULL, main_plot = FALSE) {
    
    p <- ggplot(dat, aes(x = age, y = HR, color = Group, fill = Group)) +
      geom_ribbon(
        aes(ymin = lower_plot, ymax = upper_plot),
        alpha = if (main_plot) 0.10 else 0.15,
        color = NA,
        show.legend = FALSE
      ) +
      geom_line(aes(linetype = Group), linewidth = 1.15) +
      geom_hline(yintercept = 1, linetype = "dotted", linewidth = 0.7) +
      scale_color_manual(
        values = c(
          "AKD with recovery"    = "#2ECC71",
          "AKD without recovery" = "#E74C3C"
        )
      ) +
      scale_fill_manual(
        values = c(
          "AKD with recovery"    = "#2ECC71",
          "AKD without recovery" = "#E74C3C"
        )
      ) +
      scale_linetype_manual(
        values = c(
          "AKD with recovery"    = "solid",
          "AKD without recovery" = "dashed"
        )
      ) +
      labs(
        x = "Age (years)",
        y = "Hazard ratio vs non-AKD",
        color = "AKD group",
        linetype = "AKD group"
      ) +
      theme_classic(base_family = base_family, base_size = 12) +
      theme(
        legend.position = "top",
        legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold")
      )
    
    if (main_plot) {
      p <- p +
        geom_segment(
          data = arrow_up,
          aes(x = age, xend = age, y = ymax - 0.10, yend = ymax),
          inherit.aes = FALSE,
          color = "black",
          linewidth = 0.35,
          arrow = arrow(length = unit(0.10, "inches"), type = "closed")
        ) +
        geom_segment(
          data = arrow_low,
          aes(x = age, xend = age, y = ymin + 0.10, yend = ymin),
          inherit.aes = FALSE,
          color = "black",
          linewidth = 0.35,
          arrow = arrow(length = unit(0.10, "inches"), type = "closed")
        ) +
        coord_cartesian(
          xlim = xlim,
          ylim = c(ymin, ymax),
          expand = FALSE
        ) +
        scale_x_continuous(breaks = seq(xlim[1], xlim[2], by = 5)) +
        scale_y_continuous(breaks = c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0))
    } else {
      p <- p +
        geom_segment(
          data = arrow_up,
          aes(x = age, xend = age, y = ymax * 0.96, yend = ymax),
          inherit.aes = FALSE,
          color = "black",
          linewidth = 0.35,
          arrow = arrow(length = unit(0.10, "inches"), type = "closed")
        ) +
        geom_segment(
          data = arrow_low,
          aes(x = age, xend = age, y = ymin * 1.04, yend = ymin),
          inherit.aes = FALSE,
          color = "black",
          linewidth = 0.35,
          arrow = arrow(length = unit(0.10, "inches"), type = "closed")
        ) +
        coord_cartesian(
          ylim = c(ymin, ymax),
          expand = TRUE
        )
    }
    
    p
  }
  
  # ==========================
  # 10) Create plots
  # ==========================
  p_main <- make_plot(
    dat = plot_main,
    arrow_up = arrow_main_upper,
    arrow_low = arrow_main_lower,
    ymin = ymin_main,
    ymax = ymax_main,
    xlim = c(40, 85),
    main_plot = TRUE
  )
  
  p_full <- make_plot(
    dat = plot_full,
    arrow_up = arrow_full_upper,
    arrow_low = arrow_full_lower,
    ymin = ymin_full,
    ymax = ymax_full,
    main_plot = FALSE
  )
  
  # ==========================
  # 11) Print check
  # ==========================
  print(p_main)
  print(p_full)
  
  # ==========================
  # 12) Save files
  # ==========================
  pdf(out_pdf_main, width = 180/25.4, height = 120/25.4, family = base_family)
  print(p_main)
  dev.off()
  
  ragg::agg_tiff(
    out_tiff_main,
    width = 180/25.4,
    height = 120/25.4,
    units = "in",
    res = 600,
    compression = "lzw"
  )
  print(p_main)
  dev.off()
  
  pdf(out_pdf_full, width = 180/25.4, height = 120/25.4, family = base_family)
  print(p_full)
  dev.off()
  
  ragg::agg_tiff(
    out_tiff_full,
    width = 180/25.4,
    height = 120/25.4,
    units = "in",
    res = 600,
    compression = "lzw"
  )
  print(p_full)
  dev.off()
  
  cat("\nSaved files:\n")
  cat(" Main PDF : ", out_pdf_main, "\n")
  cat(" Main TIFF: ", out_tiff_main, "\n")
  cat(" Full PDF : ", out_pdf_full, "\n")
  cat(" Full TIFF: ", out_tiff_full, "\n")
}
#年齢ごとのフォレストプロットのみ####
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
  outdir <- file.path(getwd(), "figure_table")
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
  # <40 は実データ内の <40 歳の中央値を使う
  # ≥80 は実データ内の >=80 歳の中央値を使う
  rep_age_lt40 <- dat %>%
    filter(age < 40) %>%
    summarise(rep_age = median(age, na.rm = TRUE)) %>%
    pull(rep_age)
  
  rep_age_ge80 <- dat %>%
    filter(age >= 80) %>%
    summarise(rep_age = median(age, na.rm = TRUE)) %>%
    pull(rep_age)
  
  # データがない場合の保険
  if (length(rep_age_lt40) == 0 || is.na(rep_age_lt40)) rep_age_lt40 <- 35
  if (length(rep_age_ge80) == 0 || is.na(rep_age_ge80)) rep_age_ge80 <- 82
  
  age_band_df <- tibble(
    band_label = c("<40", "40-49", "50-59", "60-69", "70-79", "≥80"),
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
          age_band = age_band_df$band_label[i]
        )
    })
  ) %>%
    mutate(
      age_band = factor(age_band, levels = rev(age_band_df$band_label)),
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
        Group == "AKD with recovery"    ~ as.numeric(age_band) + 0.16,
        Group == "AKD without recovery" ~ as.numeric(age_band) - 0.16
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
      breaks = seq_along(levels(plot_df$age_band)),
      labels = levels(plot_df$age_band),
      expand = expansion(mult = c(0.08, 0.08))
    ) +
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
