{
  ############################################################
  # Supplement Figure 4
  # Overall + Age subgroup (<75 / ≧75)
  # Table columns: N, Deaths(%), Risk(%), RR, RD(pp), HR(95%CI), P
  # Group labels: non-AKD / AKD with recovery / AKD without recovery
  # + Footnote/legend explaining RR/RD/HR
  ############################################################
  
  graphics.off()
  
  # ==========================
  # Packages
  # ==========================
  pkgs <- c("readr","dplyr","survival","broom","tibble",
            "forestploter","grid","gridExtra")
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
  
  # ==========================
  # Paths
  # ==========================
  setwd("X:/R")
  in_csv <- "jin1_Eligibile.csv"
  out_tiff <- "Supplementary_Figure4_ForestPlot_Overall_AgeSubgroup.tiff"
  out_pdf  <- "Supplementary_Figure4_ForestPlot_Overall_AgeSubgroup.pdf"
  
  # ==========================
  # Font (≧ を確実に出す)
  #  - Windowsなら "Yu Gothic" or "MS Gothic" がだいたい使えます
  #  - 使える方を自動で選びます
  # ==========================
  base_family <- {
    f <- c("Yu Gothic", "MS Gothic", "Meiryo", "Arial Unicode MS", "Arial")
    ok <- f[f %in% names(grDevices::windowsFonts())]
    if (length(ok) == 0) "sans" else ok[1]
  }
  
  # ==========================
  # Helper: display labels
  # ==========================
  disp_group <- function(x){
    dplyr::recode(
      as.character(x),
      "nonAKD"       = "non-AKD",
      "Recovery"     = "AKD with recovery",
      "Non-Recovery" = "AKD without recovery"
    )
  }
  
  disp_agegrp <- function(x){
    dplyr::recode(
      as.character(x),
      "lt75" = "<75",
      "ge75" = "≧75"
    )
  }
  
  # ==========================
  # 0) Load
  # ==========================
  jin1_Eligibile <- read_csv(in_csv, locale = locale(encoding = "SHIFT-JIS"))
  
  # ==========================
  # 1) 1人1行 + 3群 + follow-up + age75
  #   ※ age75_f はモデル内部では ASCII (lt75/ge75) にしておく（係数名が安定）
  # ==========================
  dat_age <- jin1_Eligibile %>%
    group_by(id) %>%
    arrange(index_date, date, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    filter(jin_status %in% c("nonAKD","AKD")) %>%
    mutate(
      jin_label = case_when(
        jin_status == "nonAKD" ~ "nonAKD",
        jin_status == "AKD" & `150_210recovery` == 1 ~ "Recovery",
        jin_status == "AKD" & `150_210recovery` == 2 ~ "Non-Recovery",
        jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` == 1 ~ "Recovery",
        jin_status == "AKD" & `150_210recovery` == 0 & `90_150recovery` %in% c(0, 2) ~ "Non-Recovery",
        TRUE ~ NA_character_
      ),
      jin_label = factor(jin_label, levels = c("nonAKD","Recovery","Non-Recovery")),
      arb_acei_use = if_else(coalesce(arb, 0) == 1 | coalesce(acei, 0) == 1, 1L, 0L),
      time_years   = as.numeric(last_follow_death - index_plus_210) / 365.25,
      age75        = case_when(is.na(age) ~ NA_integer_, age >= 75 ~ 1L, TRUE ~ 0L),
      age75_f      = factor(age75, levels = c(0,1), labels = c("lt75","ge75"))  # ←ASCII
    ) %>%
    filter(
      !is.na(jin_label),
      !is.na(time_years), time_years >= 0,
      !is.na(primary_death),
      !is.na(age), !is.na(index_cre),
      !is.na(age75_f)
    )
  
  # （ログ用）nとイベント数
  dat_age %>% count(age75_f, jin_label, name = "n") %>% print(n = Inf)
  dat_age %>% group_by(age75_f, jin_label) %>%
    summarise(deaths = sum(primary_death), n = n(), .groups = "drop") %>% print(n = Inf)
  
  # ==========================
  # 2) Cox models (overall + interaction)
  # ==========================
  cox_overall <- coxph(
    Surv(time_years, primary_death) ~
      jin_label + age + index_cre + arb_acei_use +
      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
    data = dat_age
  )
  
  cox_base_age <- coxph(
    Surv(time_years, primary_death) ~
      jin_label + age75_f + age + index_cre + arb_acei_use +
      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
    data = dat_age
  )
  
  cox_int_age <- coxph(
    Surv(time_years, primary_death) ~
      jin_label * age75_f + age + index_cre + arb_acei_use +
      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
    data = dat_age
  )
  
  cat("\n--- LRT for interaction (age75_f) ---\n")
  print(anova(cox_base_age, cox_int_age, test = "LRT"))
  
  # ==========================
  # 3) Summary stats: N / Deaths / Risk
  # ==========================
  calculate_summary_stats <- function(data, subgroup_var = NULL) {
    if (is.null(subgroup_var)) {
      data %>%
        group_by(jin_label) %>%
        summarise(
          Number = n(),
          Deaths = sum(primary_death),
          Death_pct = sprintf("%.1f", 100 * Deaths / Number),
          risk = Deaths / Number,
          .groups = "drop"
        ) %>%
        mutate(Subgroup = "Overall", Group = as.character(jin_label)) %>%
        select(Subgroup, Group, Number, Deaths, Death_pct, risk)
    } else {
      data %>%
        group_by(!!sym(subgroup_var), jin_label) %>%
        summarise(
          Number = n(),
          Deaths = sum(primary_death),
          Death_pct = sprintf("%.1f", 100 * Deaths / Number),
          risk = Deaths / Number,
          .groups = "drop"
        ) %>%
        mutate(Subgroup = as.character(!!sym(subgroup_var)),
               Group = as.character(jin_label)) %>%
        select(Subgroup, Group, Number, Deaths, Death_pct, risk)
    }
  }
  
  stats_overall <- calculate_summary_stats(dat_age)
  stats_age <- dat_age %>% calculate_summary_stats("age75_f")  # lt75/ge75
  
  # ==========================
  # 4) HR extraction (overall)
  # ==========================
  extract_hr_ci_p <- function(fit, ref_level = "nonAKD") {
    tidy_res <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
    
    hr_data <- tidy_res %>%
      filter(grepl("^jin_label", term)) %>%
      mutate(
        Group = case_when(
          term == "jin_labelRecovery" ~ "Recovery",
          term == "jin_labelNon-Recovery" ~ "Non-Recovery",
          TRUE ~ NA_character_
        ),
        HR = estimate, CI_low = conf.low, CI_high = conf.high, P_value = p.value
      ) %>%
      select(Group, HR, CI_low, CI_high, P_value)
    
    bind_rows(
      tibble(Group = ref_level, HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
      hr_data
    )
  }
  
  hr_overall <- extract_hr_ci_p(cox_overall) %>% mutate(Subgroup = "Overall")
  
  # ==========================
  # 5) HR by age subgroup from interaction model
  #   ※ age75_fge75 の相互作用を拾う（ASCIIで安定）
  # ==========================
  hr_int_table_age75 <- function(fit){
    cf <- coef(fit); vc <- vcov(fit); nm <- names(cf)
    
    find_name <- function(patterns) {
      hit <- nm[Reduce(`|`, lapply(patterns, function(p) grepl(p, nm)))]
      if (length(hit) == 0) return(NA_character_)
      hit[1]
    }
    
    b_rec  <- "jin_labelRecovery"
    b_nrec <- "jin_labelNon-Recovery"
    if (!b_rec %in% nm)  stop("主効果が見つかりません: ", b_rec)
    if (!b_nrec %in% nm) stop("主効果が見つかりません: ", b_nrec)
    
    # interaction terms for ge75
    int_rec  <- find_name(c(paste0("^", b_rec, ":age75_fge75$"),
                            paste0("^age75_fge75:", b_rec, "$")))
    int_nrec <- find_name(c(paste0("^", b_nrec, ":age75_fge75$"),
                            paste0("^age75_fge75:", b_nrec, "$")))
    
    comp_one <- function(b, int, subgroup){
      if (subgroup == "lt75") {
        est <- cf[b]; v <- vc[b,b]
      } else {
        if (is.na(int)) stop("相互作用の係数が見つかりません（", b, " × ge75）。")
        est <- cf[b] + cf[int]
        v   <- vc[b,b] + vc[int,int] + 2*vc[b,int]
      }
      se <- sqrt(v)
      HR <- exp(est); lo <- exp(est - 1.96*se); hi <- exp(est + 1.96*se)
      z  <- est / se; p_value <- 2 * pnorm(-abs(z))
      c(HR = HR, CI_low = lo, CI_high = hi, P_value = p_value)
    }
    
    out <- rbind(
      c("lt75", "Recovery",     comp_one(b_rec,  int_rec,  "lt75")),
      c("lt75", "Non-Recovery", comp_one(b_nrec, int_nrec, "lt75")),
      c("ge75", "Recovery",     comp_one(b_rec,  int_rec,  "ge75")),
      c("ge75", "Non-Recovery", comp_one(b_nrec, int_nrec, "ge75"))
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
  
  hr_age0 <- hr_int_table_age75(cox_int_age)
  
  hr_age <- bind_rows(
    tibble(Subgroup = "lt75", Group = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
    dplyr::filter(hr_age0, Subgroup == "lt75"),
    tibble(Subgroup = "ge75", Group = "nonAKD", HR = 1, CI_low = 1, CI_high = 1, P_value = NA_real_),
    dplyr::filter(hr_age0, Subgroup == "ge75")
  ) %>%
    distinct()
  
  # ==========================
  # 6) Build forest table data
  #   + RR/RD (descriptive)
  #   + display labels（ここで≧75・Group表記を反映）
  # ==========================
  # ==========================
  # Fix: "≧" glyph fallback
  #  - try "≧" (U+2267)
  #  - if not available, use "≥" (U+2265)
  #  - if that also fails, use ">="
  # ==========================
  # ==========================
  # forest_data_full (your code + fixes)
  # ==========================
  forest_data_full <- bind_rows(
    stats_overall %>% left_join(hr_overall, by = c("Subgroup","Group")),
    stats_age     %>% left_join(hr_age,     by = c("Subgroup","Group"))
  ) %>%
    mutate(
      # Subgroup order
      Subgroup = factor(Subgroup, levels = c("Overall","lt75","ge75")),
      Group    = factor(Group, levels = c("nonAKD","Recovery","Non-Recovery"))
    ) %>%
    arrange(Subgroup, Group) %>%
    group_by(Subgroup) %>%
    mutate(
      risk_ref = risk[Group == "nonAKD"][1],
      RR = if_else(Group == "nonAKD", NA_real_, risk / risk_ref),
      RD = if_else(Group == "nonAKD", NA_real_, risk - risk_ref)
    ) %>%
    ungroup() %>%
    mutate(
      # ---- left subgroup column with indent ----
      Subgroup_disp = case_when(
        Subgroup == "Overall" & Group == "nonAKD" ~ "Overall",
        
        Subgroup == "lt75" & Group == "nonAKD"    ~ "Age subgroup",
        Subgroup == "lt75" & Group == "Recovery"  ~ "  <75",
        
        Subgroup == "ge75" & Group == "nonAKD"    ~ "Age subgroup",
        Subgroup == "ge75" & Group == "Recovery" ~ "  >=75",
        
        TRUE ~ ""
      ),
      
      # ---- display Group labels ----
      Group_disp = factor(
        disp_group(Group),
        levels = c("non-AKD","AKD with recovery","AKD without recovery")
      ),
      
      `Number`     = as.character(Number),
      `Deaths (%)` = paste0(Deaths, " (", Death_pct, ")"),
      `Risk (%)`   = sprintf("%.1f", 100 * risk),
      
      `RR vs non-AKD` = case_when(
        Group == "nonAKD" ~ "Reference",
        is.na(RR) ~ "",
        TRUE ~ sprintf("%.2f", RR)
      ),
      `RD vs non-AKD (pp)` = case_when(
        Group == "nonAKD" ~ "",
        is.na(RD) ~ "",
        TRUE ~ sprintf("%+.1f", 100 * RD)
      ),
      
      `HR (95%CI)` = if_else(
        Group == "nonAKD",
        "Reference",
        sprintf("%.2f (%.2f–%.2f)", HR, CI_low, CI_high)
      ),
      `P-value` = case_when(
        Group == "nonAKD" ~ "",
        is.na(P_value) ~ "",
        P_value < 0.001 ~ "<0.001",
        P_value < 0.01  ~ sprintf("%.3f", P_value),
        TRUE            ~ sprintf("%.2f", P_value)
      )
    )
  
  # ---- forestploter用 空白列（フォレスト描画位置）----
  forest_data_full$` ` <- paste(rep(" ", 20), collapse = " ")
  
  forest_data_plot <- forest_data_full %>%
    select(
      Subgroup = Subgroup_disp,
      Group = Group_disp,
      Number,
      `Deaths (%)`,
      `Risk (%)`,
      `RR vs non-AKD`,
      `RD vs non-AKD (pp)`,
      ` `,
      `HR (95%CI)`,
      `P-value`
    )
  
  # ==========================
  # 7) Forest plot
  # ==========================
  forest_plot_combined <- forestploter::forest(
    data  = forest_data_plot,
    est   = forest_data_full$HR,
    lower = forest_data_full$CI_low,
    upper = forest_data_full$CI_high,
    sizes = 0.6,
    ci_column = 8,
    ref_line  = 1,
    x_trans   = "log",
    xlim      = c(0.5, 10),
    ticks_at  = c(0.5, 1, 2, 4, 8),
    arrow_lab = c("Higher risk", "Lower risk")
  )
  
  # 右揃え（数値列）
  forest_plot_combined <- edit_plot(
    forest_plot_combined,
    col = c(3,4,5,6,7),
    which = "text",
    hjust = unit(1, "npc"),
    x = unit(1, "npc")
  )
  
  # ヘッダー下に横線
  forest_plot_combined <- add_border(
    forest_plot_combined,
    row = 0,
    where = "bottom",
    gp = gpar(lwd = 1)
  )
  
  # 行間
  forest_plot_combined$heights <- rep(unit(7.5, "mm"), nrow(forest_plot_combined))
  
  # ==========================
  # Title
  # ==========================
  title_grob <- grid::textGrob(
    "Supplemental Figure 4. Adjusted Hazard Ratios by Age (Interaction Model)",
    x = grid::unit(0, "npc"),
    hjust = 0,
    gp = grid::gpar(fontface = "bold", cex = 1.1, fontfamily = base_family)
  )
  
  # ==========================
  # 8) Footnote (line-by-line + wrap)
  # ==========================
  legend_lines <- c(
    "Abbreviations: AKD, acute kidney disease; HR, hazard ratio.",
    "RR, relative risk (crude risk ratio) compared with non-AKD within each subgroup.",
    "RD, risk difference (percentage points) compared with non-AKD within each subgroup.",
    "Risk (%), RR, and RD are descriptive estimates based on crude proportions.",
    "HRs are adjusted estimates from Cox proportional hazards models."
  )
  
  # 折り返し（1行が長すぎてはみ出すのを防ぐ）
  # 目安: 95〜110くらいで調整（幅が足りないなら小さく）
  legend_txt <- paste(unlist(strwrap(legend_lines, width = 105)), collapse = "\n")
  
  foot_grob <- grid::textGrob(
    legend_txt,
    x = grid::unit(0, "npc"), hjust = 0, just = "left",
    gp = grid::gpar(cex = 0.72, fontfamily = base_family)  # 少し小さく
  )
  
  # forest + footnote を縦結合（脚注領域を増やす）
  final_grob <- gridExtra::arrangeGrob(
    title_grob,
    forest_plot_combined,
    foot_grob,
    ncol = 1,
    heights = grid::unit.c(
      grid::unit(8, "mm"),                        # タイトル高さ
      grid::unit(1, "npc") - grid::unit(32, "mm"),# 本体
      grid::unit(24, "mm")                        # 脚注
    )
  )
  
  grid::grid.newpage()
  grid::grid.draw(final_grob)
  # ==========================
  # 9) Save (脚注ぶん高さを増やす)
  # ==========================
  w_mm <- 360
  h_mm <- 150   # ← ★ここを増やす（切れるなら 160 まで上げてOK）
  
  # ---- TIFF ----
  tiff(out_tiff, width = w_mm, height = h_mm, units = "mm",
       res = 600, compression = "lzw", family = base_family)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(x = grid::unit(0.002, "npc"), y = grid::unit(0.003, "npc"),
                                    width = grid::unit(0.996, "npc"), height = grid::unit(0.994, "npc"),
                                    just = c("left","bottom")))
  grid::grid.draw(final_grob)
  grid::popViewport()
  dev.off()
  # ---- PDF（Cairoなし）----
  pdf(out_pdf, width = w_mm/25.4, height = h_mm/25.4, family = base_family)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(x = grid::unit(0.002, "npc"), y = grid::unit(0.003, "npc"),
                                    width = grid::unit(0.996, "npc"), height = grid::unit(0.994, "npc"),
                                    just = c("left","bottom")))
  grid::grid.draw(final_grob)
  grid::popViewport()
  dev.off()
}#75歳区切りのHR表
{
  ############################################################
  # Continuous age × AKD interaction (Cox)
  #  - HR(age) curves (Recovery / Non-Recovery vs nonAKD)
  #  - Save as Supplementary Figure 5 (PDF + TIFF 600dpi)
  ############################################################
  
  library(dplyr)
  library(survival)
  library(ggplot2)
  
  # ==========================
  # 1) Analysis dataset (from dat_age)
  # ==========================
  dat <- dat_age %>%
    filter(exclude == "include",
           jin_status %in% c("AKD","nonAKD")) %>%
    distinct(id, .keep_all = TRUE) %>%
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
      time_years = as.numeric(last_follow_death - index_plus_210)/365.25
    ) %>%
    filter(!is.na(jin_label),
           !is.na(time_years), time_years >= 0,
           !is.na(primary_death),
           !is.na(age),
           !is.na(index_cre))
  
  # ==========================
  # 2) Build continuous age centered variable
  # ==========================
  age_mean <- mean(dat$age, na.rm = TRUE)
  dat <- dat %>% mutate(age_c = age - age_mean)
  
  # ==========================
  # 3) Cox models (base vs interaction)
  # ==========================
  cox_base_age_cont <- coxph(
    Surv(time_years, primary_death) ~
      jin_label + age_c + index_cre + arb_acei_use +
      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
    data = dat
  )
  
  cox_int_age_cont <- coxph(
    Surv(time_years, primary_death) ~
      jin_label * age_c + index_cre + arb_acei_use +
      dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
    data = dat
  )
  
  # LRT for interaction
  print(anova(cox_base_age_cont, cox_int_age_cont, test = "LRT"))
  
  # ==========================
  # 4) Extract coef/vcov safely
  # ==========================
  cf <- coef(cox_int_age_cont)
  vc <- vcov(cox_int_age_cont)
  
  # 係数名チェック（念のため）
  print(names(cf)[grepl("^jin_label|age_c", names(cf))])
  
  # ==========================
  # 5) HR(age) calculation function
  # ==========================
  calc_hr_age <- function(age_val, b_main, b_int) {
    age_c_val <- as.numeric(age_val - age_mean)
    
    # interaction model: log(HR) = b_main + b_int * age_c
    est <- as.numeric(cf[b_main] + cf[b_int] * age_c_val)
    
    v1 <- as.numeric(vc[b_main, b_main])
    v2 <- as.numeric(vc[b_int,  b_int])
    v3 <- as.numeric(vc[b_main, b_int])
    
    var <- as.numeric(v1 + (age_c_val^2) * v2 + 2 * age_c_val * v3)
    var <- pmax(var, 0)
    se  <- sqrt(var)
    
    HR  <- exp(est)
    lo  <- exp(est - 1.96 * se)
    hi  <- exp(est + 1.96 * se)
    
    setNames(c(HR, lo, hi), c("HR", "CI_low", "CI_high"))
  }
  
  # 年齢範囲（データに合わせて安全に）
  age_min <- floor(quantile(dat$age, 0.01, na.rm = TRUE))
  age_max <- ceiling(quantile(dat$age, 0.99, na.rm = TRUE))
  age_seq <- seq(age_min, age_max, by = 1)
  
  plot_df <- do.call(rbind, lapply(age_seq, function(a){
    rec  <- calc_hr_age(a, "jin_labelRecovery", "jin_labelRecovery:age_c")
    nonr <- calc_hr_age(a, "jin_labelNon-Recovery", "jin_labelNon-Recovery:age_c")
    
    data.frame(
      age = a,
      HR_rec = rec[["HR"]],
      CI_low_rec = rec[["CI_low"]],
      CI_high_rec = rec[["CI_high"]],
      HR_nonrec = nonr[["HR"]],
      CI_low_nonrec = nonr[["CI_low"]],
      CI_high_nonrec = nonr[["CI_high"]]
    )
  })) %>%
    filter(is.finite(HR_rec), is.finite(HR_nonrec))
  
  stopifnot(nrow(plot_df) > 0)
  
  # y-axis
  ymin <- min(plot_df$CI_low_rec, plot_df$CI_low_nonrec, 1, na.rm = TRUE)
  ymax <- max(plot_df$CI_high_rec, plot_df$CI_high_nonrec, 1, na.rm = TRUE)
  
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  
  # ---- plot_df (wide) -> long へ（凡例を分かりやすくする）----
  plot_long <- plot_df %>%
    transmute(
      age,
      Recovery_HR = HR_rec,
      Recovery_lo = CI_low_rec,
      Recovery_hi = CI_high_rec,
      NonRecovery_HR = HR_nonrec,
      NonRecovery_lo = CI_low_nonrec,
      NonRecovery_hi = CI_high_nonrec
    ) %>%
    pivot_longer(
      cols = -age,
      names_to = c("Group","Metric"),
      names_pattern = "(Recovery|NonRecovery)_(HR|lo|hi)",
      values_to = "value"
    ) %>%
    pivot_wider(names_from = Metric, values_from = value) %>%
    mutate(
      Group = factor(Group, levels = c("Recovery","NonRecovery"),
                     labels = c("AKD with recovery ", "AKD without recovery")),
      linetype = if_else(Group == "AKD with recovery", "solid", "dashed")
    )
  
  # ---- y axis ----
  ymin <- min(plot_long$lo, 1, na.rm = TRUE)
  ymax <- max(plot_long$hi, 1, na.rm = TRUE)
  
  # ---- 色（ここは好きに変更OK）----
  col_rec  <- "#1B9E77"  # recovery
  col_nonr <- "#D95F02"  # non-recovery
  # ---- Figure legend (caption) text ----
  fig_legend <- paste(
    "Solid and dashed lines represent the estimated hazard ratios for AKD with recovery \nand AKD without recovery, respectively.",
    "\nShaded areas indicate 95% confidence intervals.",
    "\nAbbreviations: AKD, acute kidney disease; HR, hazard ratio; CI, confidence interval; \nCox, Cox proportional hazards model."
  )
  
  p_age_int2 <- ggplot(plot_long, aes(x = age, y = HR, color = Group, linetype = linetype)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = Group),
                alpha = 0.15, color = NA, show.legend = FALSE) +
    geom_line(linewidth = 1.2) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    coord_cartesian(ylim = c(ymin, ymax)) +
    scale_color_manual(values = c(col_rec, col_nonr)) +
    scale_fill_manual(values  = c(col_rec, col_nonr)) +
    scale_linetype_identity() +
    labs(
      x = "Age (years)",
      y = "Hazard ratio vs non-AKD",
      title = "Supplementary Figure 5. Age-dependent association \nbetween AKD and all-cause mortality",
      color = "AKD group",
      linetype = NULL,
      caption = fig_legend
    ) +
    guides(
      color = guide_legend(override.aes = list(linewidth = 2.0))
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      
      # ---- caption settings ----
      plot.caption = element_text(
        size = 10.5,
        hjust = 0,          # left align
        lineheight = 1.15,
        margin = margin(t = 8)
      ),
      
      # ---- margins (give room for caption at bottom) ----
      plot.margin = margin(t = 12, r = 8, b = 18, l = 8),
      
      plot.title = element_text(face = "bold", lineheight = 1.05),
      plot.subtitle = element_text(lineheight = 1.05)
    )
  
  print(p_age_int2)
  
  # 1年増えるごとのHR倍率（interaction係数）
  mult_rec  <- exp(cf["jin_labelRecovery:age_c"])
  mult_nonr <- exp(cf["jin_labelNon-Recovery:age_c"])
  cat("\nPer +1 year multiplicative change in HR (vs non-AKD):\n")
  cat("  Recovery     :", mult_rec,  "\n")
  cat("  Non-Recovery :", mult_nonr, "\n")
  
  # ==========================
  # 6) Save
  # ==========================
  outdir <- "X:/R"
  setwd(outdir)
  
  pdf("Supplementary_Figure5_AgeContinuousInteraction.pdf",
      width = 180/25.4, height = 120/25.4)
  print(p_age_int2)
  dev.off()
  
  tiff("Supplementary_Figure5_AgeContinuousInteraction.tiff",
       width = 180, height = 120, units = "mm",
       res = 600, compression = "lzw")
  print(p_age_int2)
  dev.off()
  
  cat("\nSaved files in:", outdir, "\n")
  
}#年齢スロープ


####以下検証用#####
library(dplyr)
library(readr)

# ---- 1) 年齢カットを複数用意（必要なものだけ残してOK）----
cuts <- c(60, 65, 75)

# ---- 2) カットごとに n と死亡数を作る関数 ----
make_n_event_table <- function(data, cut_age){
  data %>%
    mutate(
      age_group = if_else(age < cut_age,
                          paste0("<", cut_age),
                          paste0("≥", cut_age))
    ) %>%
    group_by(age_group, jin_label) %>%
    summarise(
      n = n(),
      deaths = sum(primary_death, na.rm = TRUE),
      death_pct = 100 * deaths / n,
      .groups = "drop"
    ) %>%
    mutate(
      cut = paste0("Age ", "<", cut_age, " vs ≥", cut_age),
      death_pct = sprintf("%.1f", death_pct),
      deaths_fmt = paste0(deaths, " (", death_pct, "%)")
    ) %>%
    select(cut, age_group, jin_label, n, deaths, deaths_fmt)
}

# ---- 3) 実行（複数cutを縦に結合）----
tbl_n_event <- bind_rows(lapply(cuts, function(x) make_n_event_table(dat, x)))

print(tbl_n_event, n = Inf)

# ---- 4) 保存（Supplementary Table用）----
setwd("X:/R")
write_csv(tbl_n_event, "Supplementary_Table_YoungAKD_N_and_Deaths_byAgeCut.csv")
