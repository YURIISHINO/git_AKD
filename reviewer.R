{
  
############################################################
# Sensitivity analysis
#
# 6 groups:
# 1. non-AKD
# 2. AKI only
# 3. AKD with AKI, recovery
# 4. AKD with AKI, non-recovery
# 5. AKD without AKI, recovery
# 6. AKD without AKI, non-recovery
#
# AKI_status:
#   "AKI" = AKI present
#   "nonAKI" / "nd" = AKI absent
#
# Output:
#   Kaplan-Meier / cumulative incidence figure
#   Number at risk
#   Adjusted Cox HR table
############################################################

graphics.off()

# ==========================================================
# Packages
# ==========================================================

pkgs <- c(
  "readr", "dplyr", "survival", "broom",
  "ggplot2", "grid", "gridExtra", "gtable",
  "tibble", "ragg", "officer", "flextable"
)

to_install <- pkgs[
  !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(to_install) > 0) {
  install.packages(to_install, dependencies = TRUE)
}

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
risk_fs <- 10
hr_fs   <- 10

km_line_lwd <- 1.05
km_ci_alpha <- 0.15

num_size <- 3.2
leg_text_size <- 2.8

# ==========================================================
# Paths
# ==========================================================

setwd("X:/R")

in_csv <- "jin1_Eligibile.csv"

out_dir <- file.path(
  "X:/R",
  "word_supp_tables"
)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

out_fig_tif <- file.path(
  out_dir,
  "Figure_AKI_AKD_subgroups_KM.tif"
)

out_fig_pdf <- file.path(
  out_dir,
  "Figure_AKI_AKD_subgroups_KM.pdf"
)

out_tab_tif <- file.path(
  out_dir,
  "Table_AKI_AKD_subgroups_HR.tif"
)

out_tab_pdf <- file.path(
  out_dir,
  "Table_AKI_AKD_subgroups_HR.pdf"
)

out_tab_csv <- file.path(
  out_dir,
  "Table_AKI_AKD_subgroups_HR.csv"
)

out_tab_docx <- file.path(
  out_dir,
  "Table_AKI_AKD_subgroups_HR.docx"
)


# ==========================================================
# Load data
# ==========================================================

jin1_Eligibile <- read_csv(
  in_csv,
  locale = locale(encoding = "SHIFT-JIS"),
  show_col_types = FALSE
)

problems(jin1_Eligibile)

colnames(jin1_Eligibile)


# ==========================================================
# Check original AKD / AKI classification
# ==========================================================

jin1_Eligibile %>%
  group_by(
    AKD_status,
    AKI_status,
    jin_status
  ) %>%
  summarise(
    n = n_distinct(id),
    .groups = "drop"
  ) %>%
  print(n = Inf)


# ==========================================================
# Required columns
# ==========================================================

req_cols <- c(
  "id",
  "exclude",
  "jin_status",
  "AKD_status",
  "AKI_status",
  "index_date",
  "date",
  "150_210recovery",
  "90_150recovery",
  "last_follow_death",
  "index_plus_210",
  "primary_death",
  "age",
  "index_cre",
  "arb",
  "acei",
  "dn1",
  "dn3",
  "dn4",
  "dn5",
  "dn6",
  "dn7",
  "dn8",
  "dn9",
  "dn10",
  "dn12",
  "dn13",
  "dn14",
  "dn15"
)

missing_cols <- setdiff(
  req_cols,
  names(jin1_Eligibile)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}


# ==========================================================
# Build one-row-per-patient dataset
# ==========================================================

dat1 <- jin1_Eligibile %>%
  filter(
    exclude == "include",
    jin_status %in% c(
      "AKD",
      "nonAKD"
    )
  ) %>%
  group_by(id) %>%
  arrange(
    index_date,
    date,
    .by_group = TRUE
  ) %>%
  slice(1) %>%
  ungroup()


# ==========================================================
# Group levels
# ==========================================================

levels_full <- c(
  "non-AKD",
  "AKI only",
  "AKD with AKI, recovery",
  "AKD with AKI, non-recovery",
  "AKD without AKI, recovery",
  "AKD without AKI, non-recovery"
)


# ==========================================================
# Create recovery status + 6 groups
# ==========================================================

dat_km_pre <- dat1 %>%
  mutate(
    
    # ------------------------------------------------------
    # AKD recovery status
    # ------------------------------------------------------
    
    recovery_status = case_when(
      
      `150_210recovery` == 1 ~
        "recovery",
      
      `150_210recovery` == 2 ~
        "non-recovery",
      
      `150_210recovery` == 0 &
        `90_150recovery` == 1 ~
        "recovery",
      
      `150_210recovery` == 0 &
        `90_150recovery` %in% c(0, 2) ~
        "non-recovery",
      
      TRUE ~ NA_character_
    ),
    
    # ------------------------------------------------------
    # AKI simplified
    #
    # AKI = AKI
    # nonAKI / nd = without AKI
    # ------------------------------------------------------
    
    AKI_simple = case_when(
      
      AKI_status == "AKI" ~
        "AKI",
      
      AKI_status %in% c(
        "nonAKI",
        "nd"
      ) ~
        "without AKI",
      
      TRUE ~ NA_character_
    ),
    
    # ------------------------------------------------------
    # Six-group classification
    # ------------------------------------------------------
    
    group = case_when(
      
      # 1. non-AKD
      jin_status == "nonAKD" &
        AKI_status != "AKI" ~
        "non-AKD",
      
      # 2. AKI only
      AKD_status == "nonAKD" &
        AKI_status == "AKI" ~
        "AKI only",
      
      # 3. AKD with AKI + recovery
      AKD_status %in% c(
        "AKD_cre",
        "AKD_egfr"
      ) &
        AKI_status == "AKI" &
        recovery_status == "recovery" ~
        "AKD with AKI, recovery",
      
      # 4. AKD with AKI + non-recovery
      AKD_status %in% c(
        "AKD_cre",
        "AKD_egfr"
      ) &
        AKI_status == "AKI" &
        recovery_status == "non-recovery" ~
        "AKD with AKI, non-recovery",
      
      # 5. AKD without AKI + recovery
      AKD_status %in% c(
        "AKD_cre",
        "AKD_egfr"
      ) &
        AKI_status %in% c(
          "nonAKI",
          "nd"
        ) &
        recovery_status == "recovery" ~
        "AKD without AKI, recovery",
      
      # 6. AKD without AKI + non-recovery
      AKD_status %in% c(
        "AKD_cre",
        "AKD_egfr"
      ) &
        AKI_status %in% c(
          "nonAKI",
          "nd"
        ) &
        recovery_status == "non-recovery" ~
        "AKD without AKI, non-recovery",
      
      TRUE ~ NA_character_
    ),
    
    group = factor(
      group,
      levels = levels_full
    ),
    
    time_years =
      as.numeric(
        last_follow_death -
          index_plus_210
      ) / 365.25,
    
    primary_death =
      as.numeric(primary_death)
  )


# ==========================================================
# Check classification BEFORE KM filtering
# ==========================================================

cat(
  "\n====================================\n",
  "6-group classification before KM filtering\n",
  "====================================\n"
)

dat_km_pre %>%
  group_by(group) %>%
  summarise(
    n = n_distinct(id),
    .groups = "drop"
  ) %>%
  print(n = Inf)


# Cross-check against original variables

dat_km_pre %>%
  group_by(
    group,
    AKD_status,
    AKI_status,
    jin_status
  ) %>%
  summarise(
    n = n_distinct(id),
    .groups = "drop"
  ) %>%
  print(n = Inf)


# ==========================================================
# Check exclusion from landmark KM analysis
# ==========================================================

dat_km_check <- dat_km_pre %>%
  mutate(
    
    reason = case_when(
      
      is.na(group) ~
        "group missing",
      
      is.na(time_years) ~
        "time_years missing",
      
      time_years < 0 ~
        "time_years < 0",
      
      is.na(primary_death) ~
        "primary_death missing",
      
      TRUE ~
        "kept"
    )
  )


dat_km_check %>%
  count(
    group,
    reason
  ) %>%
  print(n = Inf)


# ==========================================================
# Final KM dataset
# ==========================================================

dat_km <- dat_km_pre %>%
  filter(
    !is.na(group),
    !is.na(time_years),
    time_years >= 0,
    !is.na(primary_death)
  )


if (nrow(dat_km) == 0) {
  stop("dat_km has 0 rows after filtering.")
}


# ==========================================================
# Final sample size by group
# ==========================================================

cat(
  "\n====================================\n",
  "Final KM sample size\n",
  "====================================\n"
)

dat_km %>%
  group_by(group) %>%
  summarise(
    n = n_distinct(id),
    deaths = sum(primary_death == 1),
    .groups = "drop"
  ) %>%
  print(n = Inf)


# ==========================================================
# Common axis
# ==========================================================

ticks_show <- seq(
  0,
  8,
  by = 2
)

x_right <- 9.5

x_right_risk <- 8.2


# ==========================================================
# Colors
# ==========================================================

pal <- c(
  
  "non-AKD" =
    "#95A5A6",
  
  "AKI only" =
    "#3498DB",
  
  "AKD with AKI, recovery" =
    "#2ECC71",
  
  "AKD with AKI, non-recovery" =
    "#E74C3C",
  
  "AKD without AKI, recovery" =
    "#16A085",
  
  "AKD without AKI, non-recovery" =
    "#8E44AD"
)


common_left_margin_pt  <- 115
common_right_margin_pt <- 60


# ==========================================================
# KM fit
# ==========================================================

fit <- survfit(
  Surv(
    time_years,
    primary_death
  ) ~ group,
  data = dat_km
)

print(fit)


# ==========================================================
# Create cumulative-incidence dataframe
# ==========================================================

s <- summary(fit)

km_df <- data.frame(
  
  time =
    s$time,
  
  surv =
    s$surv,
  
  lower =
    s$lower,
  
  upper =
    s$upper,
  
  strata =
    s$strata
  
) %>%
  
  mutate(
    
    group = factor(
      sub(
        "^group=",
        "",
        strata
      ),
      levels = levels_full
    ),
    
    cif =
      1 - surv,
    
    cif_l =
      1 - upper,
    
    cif_u =
      1 - lower
  ) %>%
  
  arrange(
    group,
    time
  )


# ==========================================================
# KM panel
#
# y-axis = 0 to 1.0
# ==========================================================

p_km <- ggplot(
  km_df,
  aes(
    x = time,
    y = cif,
    colour = group,
    fill = group
  )
) +
  
  geom_ribbon(
    aes(
      ymin = cif_l,
      ymax = cif_u
    ),
    alpha = km_ci_alpha,
    linewidth = 0,
    show.legend = FALSE
  ) +
  
  geom_step(
    linewidth = km_line_lwd,
    show.legend = FALSE
  ) +
  
  scale_x_continuous(
    breaks = ticks_show,
    limits = c(
      0,
      x_right
    ),
    expand = c(
      0,
      0
    )
  ) +
  
  scale_y_continuous(
    limits = c(
      0,
      1.0
    ),
    breaks = seq(
      0,
      1.0,
      by = 0.2
    ),
    expand = c(
      0,
      0
    )
  ) +
  
  scale_colour_manual(
    values = pal,
    breaks = levels_full
  ) +
  
  scale_fill_manual(
    values = pal,
    breaks = levels_full
  ) +
  
  labs(
    x = NULL,
    y = "Cumulative incidence"
  ) +
  
  theme_classic(
    base_size = base_fs
  ) +
  
  theme(
    
    legend.position =
      "none",
    
    axis.title.x =
      element_blank(),
    
    axis.text.x =
      element_blank(),
    
    axis.ticks.x =
      element_blank(),
    
    plot.margin =
      margin(
        t = 4,
        r = common_right_margin_pt,
        b = 10,
        l = common_left_margin_pt,
        unit = "pt"
      )
  )


# ==========================================================
# Legend panel
#
# 6 groups -> 2 rows ?~ 3 columns
# ==========================================================

leg_df <- tibble(
  
  group =
    factor(
      levels_full,
      levels = levels_full
    ),
  
  x0 =
    c(
      1,
      12,
      23,
      1,
      12,
      23
    ),
  
  y =
    c(
      1.45,
      1.45,
      1.45,
      0.75,
      0.75,
      0.75
    )
  
) %>%
  
  mutate(
    x1 = x0 + 0.8
  )


p_leg <- ggplot(
  leg_df
) +
  
  geom_rect(
    aes(
      xmin = x0,
      xmax = x1,
      ymin = y - 0.15,
      ymax = y + 0.15,
      fill = group
    ),
    alpha = 0.25,
    colour = NA
  ) +
  
  geom_segment(
    aes(
      x = x0,
      xend = x1,
      y = y,
      yend = y,
      colour = group
    ),
    linewidth = 1.1
  ) +
  
  geom_text(
    aes(
      x = x1 + 0.35,
      y = y,
      label = group
    ),
    hjust = 0,
    size = leg_text_size
  ) +
  
  scale_fill_manual(
    values = pal,
    guide = "none"
  ) +
  
  scale_colour_manual(
    values = pal,
    guide = "none"
  ) +
  
  coord_cartesian(
    xlim = c(
      0,
      35
    ),
    ylim = c(
      0.35,
      1.8
    ),
    clip = "off"
  ) +
  
  theme_void(
    base_size = base_fs
  ) +
  
  theme(
    plot.margin =
      margin(
        t = 0,
        r = common_right_margin_pt,
        b = 0,
        l = common_left_margin_pt,
        unit = "pt"
      )
  )


# ==========================================================
# Number at risk
# ==========================================================

y_map <- c(
  
  "non-AKD" =
    0.70,
  
  "AKI only" =
    0.60,
  
  "AKD with AKI, recovery" =
    0.50,
  
  "AKD with AKI, non-recovery" =
    0.40,
  
  "AKD without AKI, recovery" =
    0.30,
  
  "AKD without AKI, non-recovery" =
    0.20
)


sfit <- summary(
  fit,
  times = ticks_show,
  extend = TRUE
)


time0_shift <- 0.22

label_pad_mm <- 10


risk_df <- data.frame(
  
  time =
    sfit$time,
  
  strata =
    sfit$strata,
  
  n_risk =
    sfit$n.risk
  
) %>%
  
  mutate(
    
    group =
      factor(
        sub(
          "^group=",
          "",
          strata
        ),
        levels = levels_full
      ),
    
    y =
      unname(
        y_map[
          as.character(group)
        ]
      ),
    
    time_plot =
      ifelse(
        time == 0,
        time0_shift,
        time
      ),
    
    group_disp =
      case_when(
        
        as.character(group) ==
          "AKD with AKI, recovery" ~
          "AKD with AKI,\nrecovery",
        
        as.character(group) ==
          "AKD with AKI, non-recovery" ~
          "AKD with AKI,\nnon-recovery",
        
        as.character(group) ==
          "AKD without AKI, recovery" ~
          "AKD without AKI,\nrecovery",
        
        as.character(group) ==
          "AKD without AKI, non-recovery" ~
          "AKD without AKI,\nnon-recovery",
        
        TRUE ~
          as.character(group)
      )
  ) %>%
  
  filter(
    !is.na(group),
    !is.na(y)
  )


# ==========================================================
# Risk table plot
# ==========================================================

p_risk <- ggplot() +
  
  geom_text(
    data = risk_df,
    aes(
      x = time_plot,
      y = y,
      label = n_risk
    ),
    size = num_size
  ) +
  
  scale_x_continuous(
    breaks = ticks_show,
    limits = c(
      0,
      x_right_risk
    ),
    expand = c(
      0,
      0
    )
  ) +
  
  scale_y_continuous(
    breaks =
      unname(
        y_map[levels_full]
      ),
    labels =
      rep(
        "",
        length(levels_full)
      ),
    expand = c(
      0,
      0
    )
  ) +
  
  labs(
    x = "Years",
    title = "Number at risk"
  ) +
  
  theme_classic(
    base_size = risk_fs
  ) +
  
  theme(
    
    axis.title.y =
      element_blank(),
    
    axis.text.y =
      element_blank(),
    
    axis.ticks.y =
      element_blank(),
    
    axis.line.y =
      element_blank(),
    
    plot.margin =
      margin(
        t = 10,
        r = common_right_margin_pt,
        b = 2,
        l = common_left_margin_pt,
        unit = "pt"
      ),
    
    plot.title.position =
      "plot",
    
    plot.title =
      element_text(
        margin =
          margin(
            b = 2
          ),
        vjust = 0
      )
  ) +
  
  coord_cartesian(
    ylim = c(
      0.14,
      0.75
    ),
    clip = "off"
  )


# ==========================================================
# Add risk-table group labels
# ==========================================================

for (grp in levels_full) {
  
  yv <-
    unique(
      risk_df$y[
        risk_df$group == grp
      ]
    )[1]
  
  lab <-
    unique(
      risk_df$group_disp[
        risk_df$group == grp
      ]
    )[1]
  
  col <-
    if (
      grp == "non-AKD"
    ) {
      "black"
    } else {
      pal[grp]
    }
  
  p_risk <-
    p_risk +
    
    annotation_custom(
      
      grob =
        textGrob(
          
          lab,
          
          x =
            unit(
              0,
              "npc"
            ) -
            unit(
              label_pad_mm,
              "mm"
            ),
          
          just =
            "right",
          
          gp =
            gpar(
              col = col,
              fontsize = risk_fs
            )
        ),
      
      xmin =
        -Inf,
      
      xmax =
        -Inf,
      
      ymin =
        yv,
      
      ymax =
        yv
    )
}


# ==========================================================
# Bind KM + risk
# ==========================================================

g_km <-
  ggplotGrob(
    p_km
  )

g_risk <-
  ggplotGrob(
    p_risk
  )


g_risk$widths <-
  g_km$widths


km_risk_block <-
  arrangeGrob(
    
    g_km,
    
    g_risk,
    
    ncol = 1,
    
    heights = c(
      2.0,
      1.75
    )
  )


# ==========================================================
# Final figure
# ==========================================================

fig_AKI_AKD <- arrangeGrob(
  
  p_leg,
  
  km_risk_block,
  
  nullGrob(),
  
  ncol = 1,
  
  heights = c(
    0.65,
    5.2,
    0.3
  )
)


# ==========================================================
# Preview
# ==========================================================

grid.newpage()

grid.draw(
  fig_AKI_AKD
)


# ==========================================================
# Save KM figure
# ==========================================================

ragg::agg_tiff(
  
  out_fig_tif,
  
  width = 190,
  
  height = 240,
  
  units = "mm",
  
  res = 600,
  
  compression = "lzw"
)

grid.newpage()

grid.draw(
  fig_AKI_AKD
)

dev.off()


pdf(
  
  out_fig_pdf,
  
  width = 8.2,
  
  height = 9.2
)

grid.newpage()

grid.draw(
  fig_AKI_AKD
)

dev.off()


# ==========================================================
# Cox model
# ==========================================================

dat_cox <- dat_km %>%
  
  mutate(
    
    arb_acei_use =
      if_else(
        
        coalesce(
          arb,
          0
        ) == 1 |
          
          coalesce(
            acei,
            0
          ) == 1,
        
        1L,
        
        0L
      ),
    
    # non-AKD ?𖾎??I??reference?ɂ???
    group =
      relevel(
        group,
        ref = "non-AKD"
      )
  )


# ==========================================================
# Multivariable Cox model
# ==========================================================

fit_main <- coxph(
  
  Surv(
    time_years,
    primary_death
  ) ~
    
    group +
    
    age +
    
    index_cre +
    
    arb_acei_use +
    
    dn1 +
    
    dn3 +
    
    dn4 +
    
    dn5 +
    
    dn6 +
    
    dn7 +
    
    dn8 +
    
    dn9 +
    
    dn10 +
    
    dn12 +
    
    dn13 +
    
    dn14 +
    
    dn15,
  
  data = dat_cox
)


summary(
  fit_main
)


# ==========================================================
# HR table
# ==========================================================

hr_table <- broom::tidy(
  
  fit_main,
  
  exponentiate = TRUE,
  
  conf.int = TRUE
) %>%
  
  filter(
    grepl(
      "^group",
      term
    )
  ) %>%
  
  mutate(
    
    Contrast =
      case_when(
        
        term ==
          "groupAKI only" ~
          "AKI only vs non-AKD",
        
        term ==
          "groupAKD with AKI, recovery" ~
          "AKD with AKI, recovery vs non-AKD",
        
        term ==
          "groupAKD with AKI, non-recovery" ~
          "AKD with AKI, non-recovery vs non-AKD",
        
        term ==
          "groupAKD without AKI, recovery" ~
          "AKD without AKI, recovery vs non-AKD",
        
        term ==
          "groupAKD without AKI, non-recovery" ~
          "AKD without AKI, non-recovery vs non-AKD",
        
        TRUE ~
          term
      ),
    
    HR =
      sprintf(
        "%.2f",
        estimate
      ),
    
    `95% CI` =
      sprintf(
        "%.2f???%.2f",
        conf.low,
        conf.high
      ),
    
    `p-value` =
      case_when(
        
        p.value < 0.001 ~
          "<0.001",
        
        p.value < 0.01 ~
          "<0.01",
        
        TRUE ~
          sprintf(
            "%.2f",
            p.value
          )
      )
  ) %>%
  
  dplyr::select(
    Contrast,
    HR,
    `95% CI`,
    `p-value`
  )


print(
  hr_table
)


# ==========================================================
# Save HR CSV
# ==========================================================

write_csv(
  hr_table,
  out_tab_csv
)


# ==========================================================
# HR table grob
# ==========================================================

hr_grob <- tableGrob(
  
  hr_table,
  
  rows = NULL,
  
  theme =
    ttheme_default(
      
      base_size =
        hr_fs,
      
      core =
        list(
          
          bg_params =
            list(
              col = "black",
              lwd = 0.4
            ),
          
          fg_params =
            list(
              hjust = 0.5,
              x = 0.5
            )
        ),
      
      colhead =
        list(
          
          bg_params =
            list(
              col = "black",
              lwd = 0.6
            ),
          
          fg_params =
            list(
              fontface = "bold",
              hjust = 0.5,
              x = 0.5
            )
        )
    )
)


hr_grob <- gtable_add_grob(
  
  hr_grob,
  
  rectGrob(
    gp =
      gpar(
        fill = NA,
        col = "black",
        lwd = 1
      )
  ),
  
  t = 1,
  
  l = 1,
  
  b = nrow(
    hr_grob
  ),
  
  r = ncol(
    hr_grob
  )
)


# ==========================================================
# Draw HR table
# ==========================================================

draw_hr_table <- function() {
  
  grid.newpage()
  
  x_left <-
    unit(
      8,
      "mm"
    )
  
  y_top <-
    unit(
      145,
      "mm"
    )
  
  table_w <-
    sum(
      hr_grob$widths
    )
  
  table_h <-
    sum(
      hr_grob$heights
    )
  
  pushViewport(
    
    viewport(
      
      x =
        x_left,
      
      y =
        y_top,
      
      width =
        table_w,
      
      height =
        table_h,
      
      just =
        c(
          "left",
          "top"
        )
    )
  )
  
  grid.draw(
    hr_grob
  )
  
  popViewport()
  
  
  grid.text(
    
    paste0(
      "Abbreviations: AKI, acute kidney injury; ",
      "AKD, acute kidney disease; ",
      "HR, hazard ratio; CI, confidence interval."
    ),
    
    x =
      x_left,
    
    y =
      y_top -
      table_h -
      unit(
        6,
        "mm"
      ),
    
    just =
      "left",
    
    gp =
      gpar(
        fontsize = 9
      )
  )
}


# ==========================================================
# Preview HR table
# ==========================================================

draw_hr_table()


# ==========================================================
# Save HR table TIFF
# ==========================================================

ragg::agg_tiff(
  
  out_tab_tif,
  
  width = 200,
  
  height = 160,
  
  units = "mm",
  
  res = 600,
  
  compression = "lzw"
)

draw_hr_table()

dev.off()


# ==========================================================
# Save HR table PDF
# ==========================================================

pdf(
  
  out_tab_pdf,
  
  width = 8.3,
  
  height = 6.2
)

draw_hr_table()

dev.off()


# ==========================================================
# Save HR table Word
# ==========================================================

ft <- flextable(
  hr_table
) %>%
  
  theme_booktabs() %>%
  
  bold(
    part = "header"
  ) %>%
  
  align(
    align = "center",
    part = "all"
  ) %>%
  
  autofit()


doc <- read_docx() %>%
  
  body_add_flextable(
    ft
  ) %>%
  
  body_add_par(
    paste0(
      "Abbreviations: AKI, acute kidney injury; ",
      "AKD, acute kidney disease; ",
      "HR, hazard ratio; CI, confidence interval."
    ),
    style = "Normal"
  )


print(
  doc,
  target = out_tab_docx
)


# ==========================================================
# Final summary
# ==========================================================

cat(
  "\n========================================\n"
)

cat(
  "Analysis completed.\n"
)

cat(
  "KM figure:\n",
  out_fig_tif,
  "\n",
  out_fig_pdf,
  "\n\n"
)

cat(
  "HR table:\n",
  out_tab_tif,
  "\n",
  out_tab_pdf,
  "\n",
  out_tab_csv,
  "\n",
  out_tab_docx,
  "\n"
)

cat(
  "========================================\n"
)
} #AKD???`?ł̊??x????

{
    ############################################################
    # Observed delta eGFR trajectory + annual slope
    # Landmark-aligned version
    # time0 = index_date + 210 days
    # Top: A Observed delta eGFR trajectory (???1 year)
    # Bottom: B Adjusted annual eGFR slope bar plot (???1 year)
    # Save to: X:/R/eGFRslope
    ############################################################
    
    graphics.off()
    
    # =========================
    # Packages
    # =========================
    pkgs <- c(
      "readr","dplyr","tidyr","stringr","purrr",
      "nlme","multcomp","ggplot2","data.table",
      "grid","gridExtra","gtable","ragg"
    )
    to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
    
    library(readr)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(purrr)
    library(nlme)
    library(multcomp)
    library(ggplot2)
    library(data.table)
    library(grid)
    library(gridExtra)
    library(gtable)
    library(ragg)
    
    # =========================
    # 0) Read data
    # =========================
    setwd("X:/R")
    jin1_Eligibile <- read_csv("jin1_Eligibile.csv", locale = locale(encoding = "SHIFT-JIS"))
    
    # =========================
    # 1) Build inclusion data
    #    time0 fixed at index_date + 210
    # =========================
    jin1_inclusion <- jin1_Eligibile %>%
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
        index_plus_210 = index_date + 210
      ) %>%
      filter(!is.na(jin_label))
    
    # baseline eGFR nearest to day 210 within 150-210 days
    time0_egfr_df <- jin1_inclusion %>%
      mutate(
        days_from_index = as.numeric(date - index_date),
        dist_to_210 = abs(days_from_index - 210)
      ) %>%
      filter(days_from_index >= 150, days_from_index <= 210) %>%
      group_by(id) %>%
      arrange(dist_to_210, desc(date), .by_group = TRUE) %>%
      slice(1) %>%
      ungroup() %>%
      transmute(
        id,
        time0_egfr = egfr,
        time0_egfr_date = date
      )
    
    jin1_inclusion <- jin1_inclusion %>%
      mutate(
        time0 = index_plus_210,
        years_from_time0 = as.numeric(date - time0) / 365.25
      ) %>%
      left_join(time0_egfr_df, by = "id")
    
    # =========================
    # 2) Longitudinal data
    # =========================
    akd_time_m <- jin1_inclusion %>%
      filter(years_from_time0 >= 0)
    
    baseline_cov <- jin1_Eligibile %>%
      filter(exclude == "include") %>%
      group_by(id) %>%
      arrange(index_date) %>%
      slice(1) %>%
      ungroup() %>%
      mutate(arb_acei_use = if_else(arb == 1 | acei == 1, 1L, 0L)) %>%
      dplyr::select(
        id, age, sex, arb_acei_use,
        dn1, dn3, dn4, dn5, dn6, dn7, dn8, dn9, dn10, dn12, dn13, dn14, dn15
      )
    
    covars <- c(
      "age","sex","arb_acei_use",
      "dn1","dn3","dn4","dn5","dn6","dn7","dn8","dn9","dn10","dn12","dn13","dn14","dn15"
    )
    
    longdat <- akd_time_m %>%
      dplyr::select(-any_of(covars)) %>%
      left_join(baseline_cov, by = "id") %>%
      mutate(
        age_c        = as.numeric(scale(age, center = TRUE, scale = FALSE)),
        time0_egfr_c = as.numeric(scale(time0_egfr, center = TRUE, scale = FALSE)),
        jin_label    = factor(jin_label, levels = c("nonAKD", "Recovery", "Non-Recovery"))
      ) %>%
      filter(!is.na(time0_egfr))
    
    # =========================
    # 3) Slope-adjusted LME (???1 year only)
    # =========================
    ctrl <- lmeControl(
      maxIter = 1e8, msMaxIter = 1e8,
      opt = "optim", optimMethod = "L-BFGS-B"
    )
    
    fit_slope_adj_1y <- lme(
      egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1 +
        age_c + sex + arb_acei_use +
        dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15 +
        years_from_time0:(age_c + arb_acei_use +
                            dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15),
      random = list(id = pdSymm(~ 1 + years_from_time0)),
      data = longdat %>% filter(years_from_time0 <= 1),
      na.action = na.omit,
      method = "REML",
      control = ctrl
    )
    
    # =========================
    # 4) Extract annual slopes + contrast
    # =========================
    cf <- names(fixef(fit_slope_adj_1y))
    
    v_slope <- function(g){
      vec <- rep(0, length(cf)); names(vec) <- cf
      if ("years_from_time0" %in% cf) vec["years_from_time0"] <- 1
      if (g == "Recovery" && "years_from_time0:jin_labelRecovery" %in% cf)
        vec["years_from_time0:jin_labelRecovery"] <- 1
      if (g == "Non-Recovery" && "years_from_time0:jin_labelNon-Recovery" %in% cf)
        vec["years_from_time0:jin_labelNon-Recovery"] <- 1
      vec
    }
    
    L <- rbind(
      nonAKD         = v_slope("nonAKD"),
      Recovery       = v_slope("Recovery"),
      `Non-Recovery` = v_slope("Non-Recovery")
    )
    
    ci_slope <- suppressMessages(confint(glht(fit_slope_adj_1y, linfct = L)))
    
    df_bar <- tibble(
      group    = rownames(L),
      estimate = ci_slope$confint[, "Estimate"],
      lower    = ci_slope$confint[, "lwr"],
      upper    = ci_slope$confint[, "upr"]
    ) %>%
      mutate(group = factor(group, levels = c("nonAKD","Recovery","Non-Recovery")))
    
    # contrasts vs nonAKD
    b  <- v_slope("nonAKD")
    r  <- v_slope("Recovery")
    nr <- v_slope("Non-Recovery")
    
    K <- rbind(
      `Recovery ?| nonAKD`     = r  - b,
      `Non-Recovery ?| nonAKD` = nr - b
    )
    
    gl_contrast <- glht(fit_slope_adj_1y, linfct = K)
    ci_contrast <- suppressMessages(confint(gl_contrast))
    sm_contrast <- suppressMessages(summary(gl_contrast))
    
    df_contrast <- tibble(
      group      = c("Recovery","Non-Recovery"),
      diff_value = ci_contrast$confint[, "Estimate"],
      diff_lower = ci_contrast$confint[, "lwr"],
      diff_upper = ci_contrast$confint[, "upr"],
      diff_p     = sm_contrast$test$pvalues
    )
    
    sample_sizes <- longdat %>%
      filter(years_from_time0 <= 1) %>%
      distinct(id, jin_label) %>%
      count(jin_label, name = "n") %>%
      transmute(group = as.character(jin_label), n)
    
    df_bar <- df_bar %>%
      left_join(sample_sizes, by = "group") %>%
      left_join(df_contrast, by = "group")
    
    # =========================
    # 5) Observed delta eGFR data
    # =========================
    dt <- as.data.table(
      akd_time_m %>%
        filter(years_from_time0 >= 0, years_from_time0 <= 1) %>%
        dplyr::select(id, years_from_time0, egfr, jin_label)
    )
    
    dt <- dt[!is.na(egfr)]
    
    median_interval <- dt[order(id, years_from_time0),
                          .(d = diff(years_from_time0)), by = id]$d %>%
      median(na.rm = TRUE)
    
    window_width <- median_interval * 0.75
    half_w <- window_width / 2
    
    target_timepoints <- seq(0, 1, by = 0.25)
    ids <- unique(dt$id)
    grid_dt <- CJ(id = ids, target_time = target_timepoints)
    
    setkey(dt, id, years_from_time0)
    
    window_data_obs <- dt[
      grid_dt,
      on = .(id, years_from_time0 = target_time),
      roll = "nearest",
      nomatch = 0L,
      .(id,
        target_time = i.target_time,
        years_from_time0,
        egfr,
        jin_label)
    ]
    
    window_data_obs <- window_data_obs[abs(years_from_time0 - target_time) <= half_w]
    
    base_dt <- window_data_obs[target_time == 0, .(baseline_egfr = egfr[1]), by = id]
    
    wd2 <- merge(window_data_obs, base_dt, by = "id", all.x = FALSE, all.y = FALSE)
    wd2[, delta_egfr := egfr - baseline_egfr]
    
    plot_dat <- as_tibble(wd2) %>%
      group_by(jin_label, target_time) %>%
      summarise(
        mean_delta = mean(delta_egfr, na.rm = TRUE),
        sd         = sd(delta_egfr, na.rm = TRUE),
        n          = sum(!is.na(delta_egfr)),
        se         = sd / sqrt(n),
        lwr        = mean_delta - 1.96 * se,
        upr        = mean_delta + 1.96 * se,
        .groups    = "drop"
      ) %>%
      mutate(
        jin_label = factor(as.character(jin_label),
                           levels = c("nonAKD","Recovery","Non-Recovery"))
      )
    
    # =========================
    # 6) Colors / labels
    # =========================
    values_group <- c(
      nonAKD         = "#95A5A6",
      Recovery       = "#2ECC71",
      `Non-Recovery` = "#E74C3C"
    )
    
    labels_group <- c(
      nonAKD         = "non-AKD",
      Recovery       = "AKD with recovery",
      `Non-Recovery` = "AKD without recovery"
    )
    # ==========================================================
    # 7) Bar plot only
    #    Adjusted annual eGFR slopes (???1 year)
    #    time0 = index_date + 210 days
    # ==========================================================
    
    df_bar <- df_bar %>%
      mutate(
        x = c(1, 2, 3),
        
        diff_label = case_when(
          group == "nonAKD" ~ "Reference",
          TRUE ~ sprintf(
            "Diff: %.2f\n(%.2f, %.2f)",
            diff_value, diff_lower, diff_upper
          )
        ),
        
        p_label = case_when(
          group == "nonAKD" ~ "",
          diff_p < 0.001 ~ "p<0.001",
          diff_p < 0.01  ~ sprintf("p=%.3f", diff_p),
          TRUE           ~ sprintf("p=%.2f", diff_p)
        )
      )
    
    # ----------------------------------------------------------
    # Y-axis range
    # ----------------------------------------------------------
    
    bar_ymin <- floor(min(df_bar$lower, na.rm = TRUE) - 4.5)
    bar_ymax <- max(
      1.6,
      ceiling(max(df_bar$upper, na.rm = TRUE) + 0.8)
    )
    
    # Text positions
    n_y    <- 0.55
    diff_y <- bar_ymin + 1.9
    p_y    <- bar_ymin + 0.5
    
    
    # ==========================================================
    # 8) Adjusted annual eGFR slope figure
    # ==========================================================
    
    p_slope <- ggplot(
      df_bar,
      aes(x = x, y = estimate, fill = group)
    ) +
      
      geom_col(
        width = 0.62,
        color = NA
      ) +
      
      geom_errorbar(
        aes(ymin = lower, ymax = upper),
        width = 0.16,
        linewidth = 0.8
      ) +
      
      # Sample size
      geom_text(
        aes(
          y = n_y,
          label = paste0("n=", n)
        ),
        size = 3.8,
        fontface = "bold",
        color = "#333333"
      ) +
      
      # Difference vs non-AKD
      geom_text(
        aes(
          y = diff_y,
          label = diff_label
        ),
        size = 3.2,
        fontface = "italic",
        color = "#4D4D4D",
        lineheight = 0.95
      ) +
      
      # P-value
      geom_text(
        aes(
          y = p_y,
          label = p_label
        ),
        size = 3.4,
        fontface = "bold",
        color = "#333333"
      ) +
      
      scale_fill_manual(
        values = c(
          nonAKD = "#95A5A6",
          Recovery = "#2ECC71",
          `Non-Recovery` = "#E74C3C"
        ),
        breaks = c(
          "nonAKD",
          "Recovery",
          "Non-Recovery"
        ),
        labels = c(
          nonAKD = "Non-AKD",
          Recovery = "AKD with recovery",
          `Non-Recovery` = "AKD without recovery"
        ),
        name = "Group"
      ) +
      
      scale_x_continuous(
        breaks = NULL,
        limits = c(0.4, 3.6),
        expand = expansion(mult = c(0.02, 0.02))
      ) +
      
      scale_y_continuous(
        limits = c(bar_ymin, bar_ymax),
        expand = expansion(mult = c(0.03, 0.03))
      ) +
      
      labs(
        x = NULL,
        y = "Adjusted annual eGFR slope\n(mL/min/1.73 m2 per year)"
      ) +
      
      theme_bw(base_size = 12) +
      
      theme(
        panel.background = element_rect(
          fill = "white",
          color = NA
        ),
        
        plot.background = element_rect(
          fill = "white",
          color = NA
        ),
        
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        
        panel.border = element_rect(
          color = "black",
          fill = NA,
          linewidth = 0.6
        ),
        
        axis.title = element_text(
          size = 12,
          color = "black"
        ),
        
        axis.text.y = element_text(
          size = 10,
          color = "black"
        ),
        
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        
        legend.position = "bottom",
        legend.direction = "horizontal",
        
        legend.title = element_text(
          size = 11,
          face = "bold"
        ),
        
        legend.text = element_text(
          size = 10
        ),
        
        legend.key.size = unit(
          0.9,
          "lines"
        ),
        
        legend.key = element_rect(
          fill = "white",
          color = NA
        ),
        
        plot.margin = margin(
          t = 4,
          r = 4,
          b = 4,
          l = 4
        )
      ) +
      
      guides(
        fill = guide_legend(
          nrow = 1,
          byrow = TRUE
        )
      )
    
    
    # ==========================================================
    # 9) Display
    # ==========================================================
    
    print(p_slope)
    
    
    # ==========================================================
    # 10) Save
    # ==========================================================
    
    outdir <- "X:/R/eGFRslope"
    
    if (!dir.exists(outdir)) {
      dir.create(
        outdir,
        recursive = TRUE
      )
    }
    
    
    # TIFF
    ragg::agg_tiff(
      filename = file.path(
        outdir,
        "adjusted_annual_eGFR_slope_1y_landmark210.tiff"
      ),
      width = 160,
      height = 120,
      units = "mm",
      res = 600,
      compression = "lzw"
    )
    
    print(p_slope)
    
    dev.off()
    
    
    # PDF
    ggsave(
      filename = file.path(
        outdir,
        "adjusted_annual_eGFR_slope_1y_landmark210.pdf"
      ),
      plot = p_slope,
      width = 160,
      height = 120,
      units = "mm",
      device = cairo_pdf
    )
    
  } #index_date+210?Ɉ??Ԃ?????eGFR?l??time0?ɂ????iindex_date+150~210?Łj

{
  # ==========================================================
  # Sensitivity analysis:
  # Adjusted annual eGFR slopes with 95% CIs
  #
  # time0 =
  #   first available eGFR measurement AFTER
  #   index_date + 210 days
  #
  # Analysis period:
  #   first 1 year from patient-specific time0
  #
  # Figure:
  #   adjusted annual eGFR slopes only
  #
  # Save to:
  #   X:/R/eGFRslope
  # ==========================================================
  
  graphics.off()
  
  
  # ==========================================================
  # 0. Packages
  # ==========================================================
  
  pkgs <- c(
    "readr",
    "dplyr",
    "tidyr",
    "stringr",
    "purrr",
    "nlme",
    "multcomp",
    "ggplot2",
    "ragg",
    "grid"
  )
  
  to_install <- pkgs[
    !vapply(
      pkgs,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]
  
  if (length(to_install) > 0) {
    install.packages(
      to_install,
      dependencies = TRUE
    )
  }
  
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(nlme)
  library(multcomp)
  library(ggplot2)
  library(ragg)
  library(grid)
  
  
  # ==========================================================
  # 1. Read data
  # ==========================================================
  
  setwd("X:/R")
  
  jin1_Eligibile <- read_csv(
    "jin1_Eligibile.csv",
    locale = locale(
      encoding = "SHIFT-JIS"
    )
  )
  
  
  # ==========================================================
  # 2. Define study groups
  # ==========================================================
  
  jin1_inclusion <- jin1_Eligibile %>%
    
    filter(
      exclude == "include",
      jin_status %in% c(
        "AKD",
        "nonAKD"
      )
    ) %>%
    
    mutate(
      
      jin_label = case_when(
        
        jin_status == "nonAKD" ~
          "nonAKD",
        
        jin_status == "AKD" &
          `150_210recovery` == 1 ~
          "Recovery",
        
        jin_status == "AKD" &
          `150_210recovery` == 2 ~
          "Non-Recovery",
        
        jin_status == "AKD" &
          `150_210recovery` == 0 &
          `90_150recovery` == 1 ~
          "Recovery",
        
        jin_status == "AKD" &
          `150_210recovery` == 0 &
          `90_150recovery` %in% c(0, 2) ~
          "Non-Recovery",
        
        TRUE ~ NA_character_
      ),
      
      jin_label = factor(
        jin_label,
        levels = c(
          "nonAKD",
          "Recovery",
          "Non-Recovery"
        )
      ),
      
      # Landmark date
      landmark_date =
        index_date + 210
    ) %>%
    
    filter(
      !is.na(jin_label)
    )
  
  
  # ==========================================================
  # 3. Identify patient-specific time0
  #
  # time0 =
  #   first eGFR measurement strictly AFTER
  #   index_date + 210 days
  #
  # Example:
  #
  # index                 day 210       first eGFR
  # |-----------------------|-------------??
  #                                      ??
  #                                    time0
  #
  # ==========================================================
  
  time0_egfr_df <- jin1_inclusion %>%
    
    mutate(
      
      days_after_landmark =
        as.numeric(
          date - landmark_date
        )
    ) %>%
    
    # strictly AFTER day 210
    filter(
      days_after_landmark > 0,
      !is.na(egfr)
    ) %>%
    
    group_by(id) %>%
    
    # Earliest eGFR after day 210
    arrange(
      date,
      .by_group = TRUE
    ) %>%
    
    slice(1) %>%
    
    ungroup() %>%
    
    transmute(
      
      id,
      
      # Actual measurement date becomes time0
      time0 = date,
      
      # eGFR at actual time0
      time0_egfr = egfr,
      
      # Distance from the 210-day landmark
      days_after_landmark =
        as.numeric(
          date - landmark_date
        )
    )
  
  
  # ==========================================================
  # 4. Attach time0 to longitudinal measurements
  # ==========================================================
  
  jin1_inclusion <- jin1_inclusion %>%
    
    left_join(
      time0_egfr_df,
      by = "id"
    ) %>%
    
    # Exclude patients without any eGFR after day 210
    filter(
      !is.na(time0),
      !is.na(time0_egfr)
    ) %>%
    
    mutate(
      
      # Time measured from each patient's
      # actual first eGFR after day 210
      years_from_time0 =
        as.numeric(
          date - time0
        ) / 365.25
    )
  
  
  # ==========================================================
  # 5. Longitudinal data from time0 onward
  # ==========================================================
  
  akd_time_m <- jin1_inclusion %>%
    
    filter(
      years_from_time0 >= 0
    )
  
  
  # ==========================================================
  # 6. Baseline covariates
  # ==========================================================
  
  baseline_cov <- jin1_Eligibile %>%
    
    filter(
      exclude == "include"
    ) %>%
    
    group_by(id) %>%
    
    arrange(
      index_date,
      .by_group = TRUE
    ) %>%
    
    slice(1) %>%
    
    ungroup() %>%
    
    mutate(
      
      arb_acei_use = if_else(
        coalesce(arb, 0) == 1 |
          coalesce(acei, 0) == 1,
        1L,
        0L
      )
    ) %>%
    
    dplyr::select(
      
      id,
      age,
      sex,
      arb_acei_use,
      
      dn1,
      dn3,
      dn4,
      dn5,
      dn6,
      dn7,
      dn8,
      dn9,
      dn10,
      dn12,
      dn13,
      dn14,
      dn15
    )
  
  
  # ==========================================================
  # 7. Covariate names
  # ==========================================================
  
  covars <- c(
    "age",
    "sex",
    "arb_acei_use",
    "dn1",
    "dn3",
    "dn4",
    "dn5",
    "dn6",
    "dn7",
    "dn8",
    "dn9",
    "dn10",
    "dn12",
    "dn13",
    "dn14",
    "dn15"
  )
  
  
  # ==========================================================
  # 8. Analysis dataset
  # ==========================================================
  
  longdat <- akd_time_m %>%
    
    dplyr::select(
      -any_of(covars)
    ) %>%
    
    left_join(
      baseline_cov,
      by = "id"
    ) %>%
    
    mutate(
      
      age_c =
        as.numeric(
          scale(
            age,
            center = TRUE,
            scale = FALSE
          )
        ),
      
      time0_egfr_c =
        as.numeric(
          scale(
            time0_egfr,
            center = TRUE,
            scale = FALSE
          )
        ),
      
      jin_label =
        factor(
          jin_label,
          levels = c(
            "nonAKD",
            "Recovery",
            "Non-Recovery"
          )
        )
    ) %>%
    
    filter(
      !is.na(time0_egfr)
    )
  
  
  # ==========================================================
  # 9. Check time0 distribution
  #
  # Important diagnostic:
  # How many days after day 210 was actual time0?
  # ==========================================================
  
  time0_summary <- time0_egfr_df %>%
    
    left_join(
      
      jin1_inclusion %>%
        distinct(
          id,
          jin_label
        ),
      
      by = "id"
    ) %>%
    
    distinct(
      id,
      .keep_all = TRUE
    ) %>%
    
    group_by(
      jin_label
    ) %>%
    
    summarise(
      
      n = n(),
      
      median_days =
        median(
          days_after_landmark,
          na.rm = TRUE
        ),
      
      q1_days =
        quantile(
          days_after_landmark,
          0.25,
          na.rm = TRUE
        ),
      
      q3_days =
        quantile(
          days_after_landmark,
          0.75,
          na.rm = TRUE
        ),
      
      min_days =
        min(
          days_after_landmark,
          na.rm = TRUE
        ),
      
      max_days =
        max(
          days_after_landmark,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    )
  
  print(time0_summary)
  
  
  # ==========================================================
  # 10. Number of measurements during first year
  # ==========================================================
  
  measurement_summary <- longdat %>%
    
    filter(
      years_from_time0 >= 0,
      years_from_time0 <= 1
    ) %>%
    
    group_by(
      id,
      jin_label
    ) %>%
    
    summarise(
      
      n_measurements = n(),
      
      .groups = "drop"
    ) %>%
    
    group_by(
      jin_label
    ) %>%
    
    summarise(
      
      n_patients = n(),
      
      median_measurements =
        median(
          n_measurements,
          na.rm = TRUE
        ),
      
      q1_measurements =
        quantile(
          n_measurements,
          0.25,
          na.rm = TRUE
        ),
      
      q3_measurements =
        quantile(
          n_measurements,
          0.75,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    )
  
  print(measurement_summary)
  
  
  # ==========================================================
  # 11. LME control
  # ==========================================================
  
  ctrl <- lmeControl(
    
    maxIter = 1e8,
    msMaxIter = 1e8,
    
    opt = "optim",
    optimMethod = "L-BFGS-B"
  )
  
  
  # ==========================================================
  # 12. Adjusted LME
  #
  # Analysis period:
  #   0 to 1 year after patient-specific time0
  #
  # ==========================================================
  
  fit_slope_adj_1y <- lme(
    
    egfr ~
      years_from_time0 * jin_label +
      time0_egfr_c - 1 +
      
      age_c +
      sex +
      arb_acei_use +
      
      dn1 +
      dn3 +
      dn4 +
      dn5 +
      dn6 +
      dn7 +
      dn8 +
      dn9 +
      dn10 +
      dn12 +
      dn13 +
      dn14 +
      dn15 +
      
      years_from_time0:(
        age_c +
          arb_acei_use +
          dn1 +
          dn3 +
          dn4 +
          dn5 +
          dn6 +
          dn7 +
          dn8 +
          dn9 +
          dn10 +
          dn12 +
          dn13 +
          dn14 +
          dn15
      ),
    
    random = list(
      
      id = pdSymm(
        ~ 1 + years_from_time0
      )
    ),
    
    data = longdat %>%
      filter(
        years_from_time0 >= 0,
        years_from_time0 <= 1
      ),
    
    na.action = na.omit,
    
    method = "REML",
    
    control = ctrl
  )
  
  
  # ==========================================================
  # 13. Show model
  # ==========================================================
  
  print(
    summary(
      fit_slope_adj_1y
    )
  )
  
  
  # ==========================================================
  # 14. Extract group-specific annual slopes
  # ==========================================================
  
  cf <- names(
    fixef(
      fit_slope_adj_1y
    )
  )
  
  
  v_slope <- function(g) {
    
    vec <- rep(
      0,
      length(cf)
    )
    
    names(vec) <- cf
    
    
    # Base slope
    if (
      "years_from_time0" %in% cf
    ) {
      
      vec[
        "years_from_time0"
      ] <- 1
    }
    
    
    # Recovery interaction
    if (
      g == "Recovery" &&
      "years_from_time0:jin_labelRecovery" %in% cf
    ) {
      
      vec[
        "years_from_time0:jin_labelRecovery"
      ] <- 1
    }
    
    
    # Non-Recovery interaction
    if (
      g == "Non-Recovery" &&
      "years_from_time0:jin_labelNon-Recovery" %in% cf
    ) {
      
      vec[
        "years_from_time0:jin_labelNon-Recovery"
      ] <- 1
    }
    
    
    vec
  }
  
  
  # ==========================================================
  # 15. Group-specific slopes
  # ==========================================================
  
  L <- rbind(
    
    nonAKD =
      v_slope(
        "nonAKD"
      ),
    
    Recovery =
      v_slope(
        "Recovery"
      ),
    
    `Non-Recovery` =
      v_slope(
        "Non-Recovery"
      )
  )
  
  
  gl_slope <- glht(
    
    fit_slope_adj_1y,
    
    linfct = L
  )
  
  
  ci_slope <- suppressMessages(
    
    confint(
      gl_slope
    )
  )
  
  
  df_bar <- tibble(
    
    group =
      rownames(L),
    
    estimate =
      ci_slope$confint[
        ,
        "Estimate"
      ],
    
    lower =
      ci_slope$confint[
        ,
        "lwr"
      ],
    
    upper =
      ci_slope$confint[
        ,
        "upr"
      ]
  ) %>%
    
    mutate(
      
      group =
        factor(
          group,
          levels = c(
            "nonAKD",
            "Recovery",
            "Non-Recovery"
          )
        )
    )
  
  
  # ==========================================================
  # 16. Differences vs non-AKD
  # ==========================================================
  
  b <-
    v_slope(
      "nonAKD"
    )
  
  r <-
    v_slope(
      "Recovery"
    )
  
  nr <-
    v_slope(
      "Non-Recovery"
    )
  
  
  K <- rbind(
    
    `Recovery - nonAKD` =
      r - b,
    
    `Non-Recovery - nonAKD` =
      nr - b
  )
  
  
  gl_contrast <- glht(
    
    fit_slope_adj_1y,
    
    linfct = K
  )
  
  
  ci_contrast <- suppressMessages(
    
    confint(
      gl_contrast
    )
  )
  
  
  sm_contrast <- suppressMessages(
    
    summary(
      gl_contrast
    )
  )
  
  
  df_contrast <- tibble(
    
    group = c(
      "Recovery",
      "Non-Recovery"
    ),
    
    diff_value =
      ci_contrast$confint[
        ,
        "Estimate"
      ],
    
    diff_lower =
      ci_contrast$confint[
        ,
        "lwr"
      ],
    
    diff_upper =
      ci_contrast$confint[
        ,
        "upr"
      ],
    
    diff_p =
      sm_contrast$test$pvalues
  )
  
  
  # ==========================================================
  # 17. Sample sizes
  # ==========================================================
  
  sample_sizes <- longdat %>%
    
    filter(
      years_from_time0 >= 0,
      years_from_time0 <= 1
    ) %>%
    
    distinct(
      id,
      jin_label
    ) %>%
    
    count(
      jin_label,
      name = "n"
    ) %>%
    
    transmute(
      
      group =
        as.character(
          jin_label
        ),
      
      n
    )
  
  
  # ==========================================================
  # 18. Merge plotting data
  # ==========================================================
  
  df_bar <- df_bar %>%
    
    mutate(
      group =
        as.character(
          group
        )
    ) %>%
    
    left_join(
      sample_sizes,
      by = "group"
    ) %>%
    
    left_join(
      df_contrast,
      by = "group"
    ) %>%
    
    mutate(
      
      group =
        factor(
          group,
          levels = c(
            "nonAKD",
            "Recovery",
            "Non-Recovery"
          )
        )
    )
  
  
  # ==========================================================
  # 19. Print numerical results
  # ==========================================================
  
  print(
    df_bar
  )
  
  
  # ==========================================================
  # 20. Labels for Figure
  # ==========================================================
  
  df_bar <- df_bar %>%
    
    mutate(
      
      x = c(
        1,
        2,
        3
      ),
      
      diff_label =
        case_when(
          
          group == "nonAKD" ~
            "Reference",
          
          TRUE ~
            sprintf(
              "Diff: %.2f\n(%.2f, %.2f)",
              diff_value,
              diff_lower,
              diff_upper
            )
        ),
      
      p_label =
        case_when(
          
          group == "nonAKD" ~
            "",
          
          diff_p < 0.001 ~
            "p<0.001",
          
          diff_p < 0.01 ~
            sprintf(
              "p=%.3f",
              diff_p
            ),
          
          TRUE ~
            sprintf(
              "p=%.2f",
              diff_p
            )
        )
    )
  
  
  # ==========================================================
  # 21. Y-axis
  # ==========================================================
  
  bar_ymin <- floor(
    
    min(
      df_bar$lower,
      na.rm = TRUE
    ) - 4.5
  )
  
  
  bar_ymax <- max(
    
    1.6,
    
    ceiling(
      max(
        df_bar$upper,
        na.rm = TRUE
      ) + 0.8
    )
  )
  
  
  # Text positions
  n_y <- 0.55
  
  diff_y <-
    bar_ymin + 1.9
  
  p_y <-
    bar_ymin + 0.5
  
  
  # ==========================================================
  # 22. Figure
  #    Adjusted annual eGFR slopes only
  # ==========================================================
  
  p_slope <- ggplot(
    
    df_bar,
    
    aes(
      x = x,
      y = estimate,
      fill = group
    )
  ) +
    
    # --------------------------------------------------------
  # Bars
  # --------------------------------------------------------
  
  geom_col(
    
    width = 0.62,
    
    color = NA
  ) +
    
    
    # --------------------------------------------------------
  # 95% CI
  # --------------------------------------------------------
  
  geom_errorbar(
    
    aes(
      ymin = lower,
      ymax = upper
    ),
    
    width = 0.16,
    
    linewidth = 0.8
  ) +
    
    
    # --------------------------------------------------------
  # N
  # --------------------------------------------------------
  
  geom_text(
    
    aes(
      y = n_y,
      label = paste0(
        "N = ",
        format(
          n,
          big.mark = ",",
          scientific = FALSE
        )
      )
    ),
    
    size = 4.0,
    
    fontface = "bold",
    
    color = "#333333"
  ) +
    
    
    # --------------------------------------------------------
  # Difference vs non-AKD
  # --------------------------------------------------------
  
  geom_text(
    
    aes(
      y = diff_y,
      label = diff_label
    ),
    
    size = 3.4,
    
    fontface = "italic",
    
    color = "#4D4D4D",
    
    lineheight = 0.95
  ) +
    
    
    # --------------------------------------------------------
  # P value
  # --------------------------------------------------------
  
  geom_text(
    
    aes(
      y = p_y,
      label = p_label
    ),
    
    size = 3.5,
    
    fontface = "bold",
    
    color = "#333333"
  ) +
    
    
    # --------------------------------------------------------
  # Colors
  # --------------------------------------------------------
  
  scale_fill_manual(
    
    values = c(
      
      nonAKD =
        "#95A5A6",
      
      Recovery =
        "#2ECC71",
      
      `Non-Recovery` =
        "#E74C3C"
    ),
    
    breaks = c(
      "nonAKD",
      "Recovery",
      "Non-Recovery"
    ),
    
    labels = c(
      
      nonAKD =
        "Non-AKD",
      
      Recovery =
        "AKD with recovery",
      
      `Non-Recovery` =
        "AKD without recovery"
    ),
    
    name = "Group"
  ) +
    
    
    # --------------------------------------------------------
  # X axis
  # --------------------------------------------------------
  
  scale_x_continuous(
    
    breaks = NULL,
    
    limits = c(
      0.4,
      3.6
    ),
    
    expand =
      expansion(
        mult = c(
          0.02,
          0.02
        )
      )
  ) +
    
    
    # --------------------------------------------------------
  # Y axis
  # --------------------------------------------------------
  
  scale_y_continuous(
    
    limits = c(
      bar_ymin,
      bar_ymax
    ),
    
    expand =
      expansion(
        mult = c(
          0.03,
          0.03
        )
      )
  ) +
    
    
    # --------------------------------------------------------
  # Labels
  # --------------------------------------------------------
  
  labs(
    
    x = NULL,
    
    y =
      "Adjusted annual eGFR slope\n(mL/min/1.73 m2 per year)"
  ) +
    
    
    # --------------------------------------------------------
  # Theme
  # --------------------------------------------------------
  
  theme_bw(
    base_size = 12
  ) +
    
    theme(
      
      panel.background =
        element_rect(
          fill = "white",
          color = NA
        ),
      
      plot.background =
        element_rect(
          fill = "white",
          color = NA
        ),
      
      panel.grid.major =
        element_blank(),
      
      panel.grid.minor =
        element_blank(),
      
      panel.border =
        element_rect(
          color = "black",
          fill = NA,
          linewidth = 0.6
        ),
      
      axis.title =
        element_text(
          size = 13,
          color = "black"
        ),
      
      axis.text.y =
        element_text(
          size = 11,
          color = "black"
        ),
      
      axis.text.x =
        element_blank(),
      
      axis.ticks.x =
        element_blank(),
      
      legend.position =
        "bottom",
      
      legend.direction =
        "horizontal",
      
      legend.title =
        element_blank(),
      
      legend.text =
        element_text(
          size = 11
        ),
      
      legend.key.size =
        unit(
          0.9,
          "lines"
        ),
      
      legend.key =
        element_rect(
          fill = "white",
          color = NA
        ),
      
      plot.margin =
        margin(
          t = 5,
          r = 5,
          b = 5,
          l = 5
        )
    ) +
    
    guides(
      
      fill =
        guide_legend(
          nrow = 1,
          byrow = TRUE
        )
    )
  
  
  # ==========================================================
  # 23. Display Figure
  # ==========================================================
  
  print(
    p_slope
  )
  
  
  # ==========================================================
  # 24. Output directory
  # ==========================================================
  
  outdir <-
    "X:/R/eGFRslope"
  
  if (
    !dir.exists(
      outdir
    )
  ) {
    
    dir.create(
      outdir,
      recursive = TRUE
    )
  }
  
  
  # ==========================================================
  # 25. Save TIFF
  # ==========================================================
  
  ragg::agg_tiff(
    
    filename =
      file.path(
        outdir,
        "Figure_eGFR_slope_first_eGFR_after_day210.tiff"
      ),
    
    width = 160,
    
    height = 120,
    
    units = "mm",
    
    res = 600,
    
    compression = "lzw"
  )
  
  print(
    p_slope
  )
  
  dev.off()
  
  
  # ==========================================================
  # 26. Save PDF
  # ==========================================================
  
  ggsave(
    
    filename =
      file.path(
        outdir,
        "Figure_eGFR_slope_first_eGFR_after_day210.pdf"
      ),
    
    plot =
      p_slope,
    
    width = 160,
    
    height = 120,
    
    units = "mm",
    
    device =
      cairo_pdf
  )
  
  
  # ==========================================================
  # 27. Save numerical results
  # ==========================================================
  
  write_csv(
    
    df_bar,
    
    file.path(
      outdir,
      "Results_eGFR_slope_first_eGFR_after_day210.csv"
    )
  )
  
  
  write_csv(
    
    time0_summary,
    
    file.path(
      outdir,
      "Summary_time0_after_day210.csv"
    )
  )
  
  
  write_csv(
    
    measurement_summary,
    
    file.path(
      outdir,
      "Summary_measurements_first_year.csv"
    )
  )
  
  
  # ==========================================================
  # 28. Finished
  # ==========================================================
  
  cat(
    "\n==============================================\n",
    "Sensitivity analysis completed.\n",
    "time0 = first eGFR measurement after day 210.\n",
    "Files saved in:\n",
    outdir,
    "\n==============================================\n"
  )
}#210???ȍ~???ԋ߂?eGFR??time0?ɂ???

{
  # ==========================================================
  # Cohort attrition:
  # Overall cohort -> 210-day landmark KM cohort
  # ==========================================================
  
  library(dplyr)
  library(readr)
  
  # ----------------------------------------------------------
  # 1. Read data
  # ----------------------------------------------------------
  
  jin1_Eligibile <- read_csv(
    "X:/R/jin1_Eligibile.csv",
    locale = locale(encoding = "SHIFT-JIS"),
    show_col_types = FALSE
  )
  
  
  # ----------------------------------------------------------
  # 2. One row per patient
  # ----------------------------------------------------------
  
  dat0 <- jin1_Eligibile %>%
    filter(
      exclude == "include",
      jin_status %in% c("AKD", "nonAKD")
    ) %>%
    group_by(id) %>%
    arrange(index_date, date, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      
      group = case_when(
        
        jin_status == "nonAKD" ~
          "non-AKD",
        
        jin_status == "AKD" &
          `150_210recovery` == 1 ~
          "AKD with recovery",
        
        jin_status == "AKD" &
          `150_210recovery` == 2 ~
          "AKD without recovery",
        
        jin_status == "AKD" &
          `150_210recovery` == 0 &
          `90_150recovery` == 1 ~
          "AKD with recovery",
        
        jin_status == "AKD" &
          `150_210recovery` == 0 &
          `90_150recovery` %in% c(0, 2) ~
          "AKD without recovery",
        
        TRUE ~ NA_character_
      ),
      
      landmark_date = index_plus_210,
      
      time_years =
        as.numeric(
          last_follow_death - index_plus_210
        ) / 365.25,
      
      primary_death =
        as.numeric(primary_death)
    )
  
  
  # ==========================================================
  # 3. Overall cohort numbers
  # ==========================================================
  
  overall_counts <- dat0 %>%
    count(
      group,
      name = "Overall cohort"
    )
  
  print(overall_counts)
  
  
  # ==========================================================
  # 4. Check each possible reason for exclusion
  #    These are NOT mutually exclusive
  # ==========================================================
  
  reason_overlap <- dat0 %>%
    summarise(
      
      total = n(),
      
      group_missing =
        sum(is.na(group)),
      
      last_follow_missing =
        sum(is.na(last_follow_death)),
      
      landmark_missing =
        sum(is.na(index_plus_210)),
      
      time_missing =
        sum(is.na(time_years)),
      
      followup_before_landmark =
        sum(
          !is.na(time_years) &
            time_years < 0
        ),
      
      death_variable_missing =
        sum(is.na(primary_death)),
      
      eligible_for_KM =
        sum(
          !is.na(group) &
            !is.na(time_years) &
            time_years >= 0 &
            !is.na(primary_death)
        )
    )
  
  print(reason_overlap)
  
  
  # ==========================================================
  # 5. Mutually exclusive reason for exclusion
  #
  # Priority:
  # group missing
  # -> follow-up data missing
  # -> ended before landmark
  # -> death status missing
  # -> kept
  # ==========================================================
  
  dat_reason <- dat0 %>%
    mutate(
      
      exclusion_reason = case_when(
        
        is.na(group) ~
          "Recovery group could not be classified",
        
        is.na(last_follow_death) |
          is.na(index_plus_210) |
          is.na(time_years) ~
          "Follow-up/landmark data missing",
        
        time_years < 0 ~
          "Follow-up ended before day-210 landmark",
        
        is.na(primary_death) ~
          "Death status missing",
        
        TRUE ~
          "Included in KM analysis"
      )
    )
  
  
  # ==========================================================
  # 6. Overall attrition table
  # ==========================================================
  
  attrition_total <- dat_reason %>%
    count(
      exclusion_reason,
      name = "N"
    ) %>%
    mutate(
      Percent = round(
        100 * N / sum(N),
        1
      )
    )
  
  print(attrition_total)
  
  
  # ==========================================================
  # 7. Attrition by AKD/recovery group
  # ==========================================================
  
  attrition_group <- dat_reason %>%
    count(
      group,
      exclusion_reason,
      name = "N"
    ) %>%
    arrange(
      group,
      exclusion_reason
    )
  
  print(attrition_group)
  
  
  # ==========================================================
  # 8. Overall vs KM numbers by group
  # ==========================================================
  
  comparison_group <- dat_reason %>%
    group_by(group) %>%
    summarise(
      
      Overall_N =
        n(),
      
      KM_N =
        sum(
          exclusion_reason ==
            "Included in KM analysis"
        ),
      
      Excluded_N =
        Overall_N - KM_N,
      
      Included_percent =
        round(
          100 * KM_N / Overall_N,
          1
        ),
      
      .groups = "drop"
    )
  
  print(comparison_group)
  
  
  # ==========================================================
  # 9. Separate early death from other loss to follow-up
  # ==========================================================
  #
  # This is useful only if primary_death accurately indicates
  # that death occurred before the landmark.
  # ==========================================================
  
  early_landmark_status <- dat0 %>%
    filter(
      !is.na(time_years),
      time_years < 0
    ) %>%
    mutate(
      
      status_before_landmark =
        case_when(
          
          primary_death == 1 ~
            "Death before day 210",
          
          primary_death == 0 ~
            "Follow-up ended before day 210 without recorded death",
          
          TRUE ~
            "Death status missing"
        )
    ) %>%
    count(
      group,
      status_before_landmark,
      name = "N"
    )
  
  print(early_landmark_status)
  
  
  # ==========================================================
  # 10. Final KM cohort
  # ==========================================================
  
  dat_km_check <- dat0 %>%
    filter(
      !is.na(group),
      !is.na(time_years),
      time_years >= 0,
      !is.na(primary_death)
    )
  
  final_KM_counts <- dat_km_check %>%
    count(
      group,
      name = "N_at_landmark"
    )
  
  print(final_KM_counts)
  
  
  # ==========================================================
  # 11. Confirm with survfit N
  # ==========================================================
  
  library(survival)
  
  fit_check <- survfit(
    Surv(
      time_years,
      primary_death
    ) ~ group,
    data = dat_km_check
  )
  
  print(fit_check)
  

  # ==========================================================
  # AKD without recovery:
  # Breakdown by recovery-window data availability
  # and landmark KM eligibility
  # ==========================================================
  
  library(dplyr)
  library(readr)
  
  jin1_Eligibile <- read_csv(
    "X:/R/jin1_Eligibile.csv",
    locale = locale(encoding = "SHIFT-JIS"),
    show_col_types = FALSE
  )
  
  # ----------------------------------------------------------
  # One row per patient
  # ----------------------------------------------------------
  
  dat0 <- jin1_Eligibile %>%
    filter(
      exclude == "include",
      jin_status %in% c("AKD", "nonAKD")
    ) %>%
    group_by(id) %>%
    arrange(index_date, date, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      
      group = case_when(
        jin_status == "nonAKD" ~ "non-AKD",
        
        jin_status == "AKD" &
          `150_210recovery` == 1 ~ "AKD with recovery",
        
        jin_status == "AKD" &
          `150_210recovery` == 2 ~ "AKD without recovery",
        
        jin_status == "AKD" &
          `150_210recovery` == 0 &
          `90_150recovery` == 1 ~ "AKD with recovery",
        
        jin_status == "AKD" &
          `150_210recovery` == 0 &
          `90_150recovery` %in% c(0, 2) ~ "AKD without recovery",
        
        TRUE ~ NA_character_
      ),
      
      time_years =
        as.numeric(last_follow_death - index_plus_210) / 365.25,
      
      km_eligible =
        !is.na(group) &
        !is.na(time_years) &
        time_years >= 0 &
        !is.na(primary_death),
      
      # Recovery assessment data availability
      recovery_data_status = case_when(
        
        # Primary recovery window has usable classification
        `150_210recovery` %in% c(1, 2) ~
          "Data available in day 150-210 window",
        
        # Primary window unavailable, but earlier window available
        `150_210recovery` == 0 &
          `90_150recovery` %in% c(1, 2) ~
          "No day 150-210 data, but day 90-150 data available",
        
        # No usable recovery data in either window
        `150_210recovery` == 0 &
          `90_150recovery` == 0 ~
          "No recovery-window data",
        
        TRUE ~
          "Other/unclear"
      )
    )
  
  
  # ==========================================================
  # 1. AKD without recovery only
  # ==========================================================
  
  akd_nr <- dat0 %>%
    filter(
      group == "AKD without recovery"
    )
  
  
  # ==========================================================
  # 2. Overall breakdown
  # ==========================================================
  
  akd_nr %>%
    count(
      recovery_data_status,
      km_eligible,
      name = "N"
    ) %>%
    arrange(
      recovery_data_status,
      desc(km_eligible)
    ) %>%
    print()
  
  
  # ==========================================================
  # 3. Cleaner table
  # ==========================================================
  
  breakdown_table <- akd_nr %>%
    mutate(
      KM_status = if_else(
        km_eligible,
        "Included in KM",
        "Excluded before landmark"
      )
    ) %>%
    count(
      recovery_data_status,
      KM_status,
      name = "N"
    ) %>%
    tidyr::pivot_wider(
      names_from = KM_status,
      values_from = N,
      values_fill = 0
    ) %>%
    mutate(
      Total = `Included in KM` + `Excluded before landmark`
    )
  
  print(breakdown_table)
  
  
  # ==========================================================
  # 4. Specifically count no recovery-window data
  # ==========================================================
  
  akd_nr %>%
    summarise(
      total_nonrecovery = n(),
      
      no_recovery_data_total =
        sum(recovery_data_status == "No recovery-window data"),
      
      no_recovery_data_in_KM =
        sum(
          recovery_data_status == "No recovery-window data" &
            km_eligible
        ),
      
      no_recovery_data_excluded_before_landmark =
        sum(
          recovery_data_status == "No recovery-window data" &
            !km_eligible
        )
    ) %>%
    print()  
}#?J?v?????}?C???[?Ől???????邱?Ƃɂ??Ă̏W?v

{
  # ==========================================================
  # Follow-up duration summary
  # 1) Mortality: landmark (index_date + 210) -> last follow-up
  # 2) eGFR slope: time0 -> last eGFR measurement
  # ==========================================================
  
  library(dplyr)
  library(readr)
  
  # ----------------------------------------------------------
  # 1. Read data
  # ----------------------------------------------------------
  
  jin1_Eligibile <- read_csv(
    "X:/R/jin1_Eligibile.csv",
    locale = locale(encoding = "SHIFT-JIS"),
    show_col_types = FALSE
  )
  
  
  # ==========================================================
  # 2. Define groups
  # ==========================================================
  
  dat_group <- jin1_Eligibile %>%
    filter(
      exclude == "include",
      jin_status %in% c("AKD", "nonAKD")
    ) %>%
    group_by(id) %>%
    arrange(index_date, date, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      group = case_when(
        jin_status == "nonAKD" ~ "non-AKD",
        
        jin_status == "AKD" &
          `150_210recovery` == 1 ~
          "AKD with recovery",
        
        jin_status == "AKD" &
          `150_210recovery` == 2 ~
          "AKD without recovery",
        
        jin_status == "AKD" &
          `150_210recovery` == 0 &
          `90_150recovery` == 1 ~
          "AKD with recovery",
        
        jin_status == "AKD" &
          `150_210recovery` == 0 &
          `90_150recovery` %in% c(0, 2) ~
          "AKD without recovery",
        
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(
      id,
      group,
      index_date,
      index_plus_210,
      last_follow_death
    )
  
  
  # ==========================================================
  # 3. Mortality follow-up
  #
  # landmark = index_date + 210
  # follow-up = landmark -> last_follow_death
  # ==========================================================
  
  mortality_followup <- dat_group %>%
    mutate(
      followup_days =
        as.numeric(
          last_follow_death - index_plus_210
        ),
      
      followup_years =
        followup_days / 365.25
    ) %>%
    filter(
      !is.na(group),
      !is.na(followup_days),
      followup_days >= 0
    )
  
  
  # ----------------------------------------------------------
  # Overall mortality follow-up
  # ----------------------------------------------------------
  
  mortality_overall <- mortality_followup %>%
    summarise(
      N = n(),
      
      median_days =
        median(followup_days, na.rm = TRUE),
      
      Q1_days =
        quantile(followup_days, 0.25, na.rm = TRUE),
      
      Q3_days =
        quantile(followup_days, 0.75, na.rm = TRUE),
      
      min_days =
        min(followup_days, na.rm = TRUE),
      
      max_days =
        max(followup_days, na.rm = TRUE),
      
      median_years =
        median(followup_years, na.rm = TRUE),
      
      Q1_years =
        quantile(followup_years, 0.25, na.rm = TRUE),
      
      Q3_years =
        quantile(followup_years, 0.75, na.rm = TRUE)
    )
  
  print(mortality_overall)
  
  
  # ----------------------------------------------------------
  # Mortality follow-up by group
  # ----------------------------------------------------------
  
  mortality_by_group <- mortality_followup %>%
    group_by(group) %>%
    summarise(
      N = n(),
      
      median_days =
        median(followup_days, na.rm = TRUE),
      
      Q1_days =
        quantile(followup_days, 0.25, na.rm = TRUE),
      
      Q3_days =
        quantile(followup_days, 0.75, na.rm = TRUE),
      
      min_days =
        min(followup_days, na.rm = TRUE),
      
      max_days =
        max(followup_days, na.rm = TRUE),
      
      median_years =
        median(followup_years, na.rm = TRUE),
      
      Q1_years =
        quantile(followup_years, 0.25, na.rm = TRUE),
      
      Q3_years =
        quantile(followup_years, 0.75, na.rm = TRUE),
      
      .groups = "drop"
    )
  
  print(mortality_by_group)
  
  
  # ==========================================================
  # 4. Define time0 for PRIMARY eGFR slope analysis
  #
  # Priority:
  # 1) date of maximum eGFR within day 90-210
  # 2) if unavailable, latest eGFR within day 0-90
  # 3) if unavailable, index date
  # ==========================================================
  
  egfr_dat <- jin1_Eligibile %>%
    filter(
      exclude == "include",
      jin_status %in% c("AKD", "nonAKD")
    ) %>%
    left_join(
      dat_group %>%
        dplyr::select(id, group),
      by = "id"
    ) %>%
    mutate(
      days_from_index =
        as.numeric(date - index_date)
    )
  
  
  # ----------------------------------------------------------
  # First priority:
  # maximum eGFR during day 90-210
  # ----------------------------------------------------------
  
  time0_90_210 <- egfr_dat %>%
    filter(
      days_from_index >= 90,
      days_from_index <= 210,
      !is.na(egfr)
    ) %>%
    group_by(id) %>%
    arrange(
      desc(egfr),
      date,
      .by_group = TRUE
    ) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      id,
      time0_primary = date,
      time0_source = "Maximum eGFR at day 90-210"
    )
  
  
  # ----------------------------------------------------------
  # Second priority:
  # latest eGFR during day 0-90
  # ----------------------------------------------------------
  
  time0_0_90 <- egfr_dat %>%
    filter(
      days_from_index >= 0,
      days_from_index < 90,
      !is.na(egfr)
    ) %>%
    group_by(id) %>%
    arrange(
      desc(date),
      .by_group = TRUE
    ) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      id,
      time0_second = date
    )
  
  
  # ----------------------------------------------------------
  # Combine hierarchy
  # ----------------------------------------------------------
  
  time0_primary <- dat_group %>%
    dplyr::select(
      id,
      group,
      index_date
    ) %>%
    left_join(
      time0_90_210,
      by = "id"
    ) %>%
    left_join(
      time0_0_90,
      by = "id"
    ) %>%
    mutate(
      time0 = case_when(
        !is.na(time0_primary) ~ time0_primary,
        is.na(time0_primary) & !is.na(time0_second) ~ time0_second,
        TRUE ~ index_date
      ),
      
      time0_source = case_when(
        !is.na(time0_primary) ~
          "Maximum eGFR at day 90-210",
        
        is.na(time0_primary) &
          !is.na(time0_second) ~
          "Latest eGFR at day 0-90",
        
        TRUE ~
          "Index date"
      )
    )
  
  
  # ==========================================================
  # 5. Last eGFR measurement after time0
  # ==========================================================
  
  last_egfr_after_time0 <- egfr_dat %>%
    dplyr::select(
      id,
      date,
      egfr
    ) %>%
    inner_join(
      time0_primary %>%
        dplyr::select(
          id,
          group,
          time0,
          time0_source
        ),
      by = "id"
    ) %>%
    filter(
      !is.na(egfr),
      date >= time0
    ) %>%
    group_by(
      id,
      group,
      time0,
      time0_source
    ) %>%
    summarise(
      last_egfr_date =
        max(date, na.rm = TRUE),
      
      .groups = "drop"
    ) %>%
    mutate(
      followup_days =
        as.numeric(
          last_egfr_date - time0
        ),
      
      followup_years =
        followup_days / 365.25
    )
  
  
  # ==========================================================
  # 6. eGFR follow-up overall
  # ==========================================================
  
  egfr_followup_overall <- last_egfr_after_time0 %>%
    summarise(
      N = n(),
      
      median_days =
        median(followup_days, na.rm = TRUE),
      
      Q1_days =
        quantile(followup_days, 0.25, na.rm = TRUE),
      
      Q3_days =
        quantile(followup_days, 0.75, na.rm = TRUE),
      
      min_days =
        min(followup_days, na.rm = TRUE),
      
      max_days =
        max(followup_days, na.rm = TRUE),
      
      median_years =
        median(followup_years, na.rm = TRUE),
      
      Q1_years =
        quantile(followup_years, 0.25, na.rm = TRUE),
      
      Q3_years =
        quantile(followup_years, 0.75, na.rm = TRUE)
    )
  
  print(egfr_followup_overall)
  
  
  # ==========================================================
  # 7. eGFR follow-up by group
  # ==========================================================
  
  egfr_followup_by_group <- last_egfr_after_time0 %>%
    group_by(group) %>%
    summarise(
      N = n(),
      
      median_days =
        median(followup_days, na.rm = TRUE),
      
      Q1_days =
        quantile(followup_days, 0.25, na.rm = TRUE),
      
      Q3_days =
        quantile(followup_days, 0.75, na.rm = TRUE),
      
      min_days =
        min(followup_days, na.rm = TRUE),
      
      max_days =
        max(followup_days, na.rm = TRUE),
      
      median_years =
        median(followup_years, na.rm = TRUE),
      
      Q1_years =
        quantile(followup_years, 0.25, na.rm = TRUE),
      
      Q3_years =
        quantile(followup_years, 0.75, na.rm = TRUE),
      
      .groups = "drop"
    )
  
  print(egfr_followup_by_group)
  
  
  # ==========================================================
  # 8. Check how time0 was determined
  # ==========================================================
  
  time0_source_summary <- time0_primary %>%
    count(
      group,
      time0_source,
      name = "N"
    )
  
  print(time0_source_summary)
}#?A?E?g?J???ǐՊ??Ԃ̏W?v
