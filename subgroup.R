library(readr)
library(dplyr)
library(survival)
library(broom)
library(ggplot2)

# CSVファイルをtibbleとして読み込む_藤倉用
jin1_Eligibile <- read_csv("/Users/tfuji/Dropbox/臨床研究/石野先生/石野先生_practice/rstudio-export_25.8.15/jin1_Eligibile.csv")

colnames(jin1_Eligibile)
# 1人1行（index_date当日レコードのみ・最初のindex_dateを採用）
jin1_Eligibile_unique_sub <- jin1_Eligibile %>%
  group_by(id) %>%
  arrange(index_date, date, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  # nonAKD と AKDのみ（あなたの方針に合わせて必要ならここは保持）
  filter(jin_status %in% c("nonAKD", "AKD")) %>%
  # 3群ラベル
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
    # 併用フラグ
    arb_acei_use = if_else(coalesce(arb, 0) == 1 | coalesce(acei, 0) == 1, 1L, 0L),
    # 追跡時間（年）：index_date+210日以降に統一
    time_years   = as.numeric(last_follow_death - index_plus_210) / 365.25
  ) %>%
  # 追跡開始以降のみ
  filter(!is.na(jin_label), !is.na(time_years), time_years >= 0)

# --- CKDの2値化（不明"nd"は除外。含めたい場合は0/1へ再符号化を検討） ---
dat <- jin1_Eligibile_unique_sub %>%
  mutate(
    ckd_bin = case_when(
      CKD_status == "CKD"    ~ 1L,
      CKD_status == "nonCKD" ~ 0L,
      TRUE ~ NA_integer_
    )
  ) %>%
  filter(!is.na(ckd_bin), !is.na(jin_label), !is.na(time_years), time_years >= 0)

# サブグループのn・イベント数を把握
dat %>% count(ckd_bin, jin_label, name = "n") %>% print()
dat %>% group_by(ckd_bin, jin_label) %>% summarise(events = sum(primary_death), n = n(), .groups="drop") %>% print()

# --- サブグループ別のCox（CKDはモデルから外す：サブセット内で不変だから） ---
cox_formula_sub <- as.formula(
  Surv(time_years, primary_death) ~
    jin_label + age + index_cre + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15
)

fit_ckd1 <- coxph(cox_formula_sub, data = filter(dat, ckd_bin == 1))  # CKDあり
fit_ckd0 <- coxph(cox_formula_sub, data = filter(dat, ckd_bin == 0))  # CKDなし

# PH仮定チェック（必要に応じて）
# cox.zph(fit_ckd1); cox.zph(fit_ckd0)

# --- 群間HR（Recovery/Non-Recovery vs nonAKD）だけを取り出して整形 ---
tidy_sub <- function(fit, label){
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% c("jin_labelRecovery", "jin_labelNon-Recovery")) %>%
    transmute(
      Subgroup = label,
      Comparison = recode(term,
                          "jin_labelRecovery" = "Recovery vs nonAKD",
                          "jin_labelNon-Recovery" = "Non-Recovery vs nonAKD"),
      HR = estimate,
      CI_low = conf.low,
      CI_high = conf.high,
      p = p.value
    )
}

tbl_ckd <- bind_rows(
  tidy_sub(fit_ckd1, "CKD: yes"),
  tidy_sub(fit_ckd0, "CKD: no")
)

print(tbl_ckd, n = Inf)

# --- サブグループ・フォレスト ---
# 並び順を明示（上：CKDあり → 下：CKDなし）
plot_df <- tbl_ckd %>%
  mutate(
    Comparison = factor(Comparison,
                        levels = c("Recovery vs nonAKD","Non-Recovery vs nonAKD")),
    Subgroup   = factor(Subgroup, levels = c("CKD: yes","CKD: no"))
  )

# CKDありの2本（Recovery/Non-Recovery）を上段、
# CKDなしの2本を下段に配置
ggplot(plot_df,
       aes(x = HR, y = Comparison, xmin = CI_low, xmax = CI_high)) +
  geom_point(size = 3) +
  geom_errorbarh(height = 0.18) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  facet_grid(rows = vars(Subgroup), switch = "y") +  # ← ここで上下に分ける
  labs(
    title = "Death (Subgroup Cox): AKD groups within CKD strata",
    x = "Hazard Ratio (log scale)", y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold") # 見出しを左に
  )


# 交互作用なしのベースモデル
cox_base <- coxph(
  Surv(time_years, primary_death) ~
    jin_label + ckd_bin + age + index_cre + arb_acei_use +
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
  data = dat
)

# 交互作用モデルのcode記載を追記
cox_int <- coxph(
  Surv(time_years, primary_death) ~
    jin_label * ckd_bin + age + index_cre + arb_acei_use + # `*` は主効果と交互作用の両方を含む
    dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15,
  data = dat
)

# 交互作用の尤度比検定（推奨）　#p＝ 0.6739
anova(cox_base, cox_int, test = "LRT")

# --- 最終モデル（cox_base）の結果を整形 ---
# 交互作用が有意でなかったため、このモデルの結果を主たる結果として採用する
final_results <- tidy(cox_base, exponentiate = TRUE, conf.int = TRUE)

print("--- Final Model Results (cox_base) ---")
print(final_results, n = Inf)

# 交互作用モデルから“層ごとのHR”を計算して表に
# 交互作用モデルの係数を確認（Waldでも可）
hr_int_table <- function(fit){
  cf <- coef(fit); vc <- vcov(fit); nm <- names(cf)
  
  # 相互作用の係数名を安全に拾うヘルパー
  find_name <- function(patterns) {
    hit <- nm[Reduce(`|`, lapply(patterns, function(p) grepl(p, nm)))]
    if (length(hit) == 0) return(NA_character_)
    hit[1]
  }
  
  # jin_label の主効果（Recovery/Non-Recovery）
  b_rec  <- "jin_labelRecovery"
  b_nrec <- "jin_labelNon-Recovery"
  if (!b_rec %in% nm)  stop("主効果が見つかりません: ", b_rec)
  if (!b_nrec %in% nm) stop("主効果が見つかりません: ", b_nrec)
  
  # ckd_bin の相互作用（順序が左右どちらでも拾う / factorなら ckdbin1 を許容）
  int_rec  <- find_name(c(paste0("^", b_rec, ":(ckd_bin|ckd_bin1)$"),
                          paste0("^(ckd_bin|ckd_bin1):", b_rec, "$")))
  int_nrec <- find_name(c(paste0("^", b_nrec, ":(ckd_bin|ckd_bin1)$"),
                          paste0("^(ckd_bin|ckd_bin1):", b_nrec, "$")))
  
  # HR計算の小関数（b:主効果名, int:相互作用名 or NA, ckd=0/1）
  comp_one <- function(b, int, ckd){
    if (ckd == 0) {
      est <- cf[b]; v <- vc[b,b]
    } else {
      if (is.na(int)) stop("相互作用の係数が見つかりません（", b, " × CKD）。",
                           "ckd_bin の型や水準、スパースを確認してください。")
      est <- cf[b] + cf[int]
      v   <- vc[b,b] + vc[int,int] + 2*vc[b,int]
    }
    se <- sqrt(v); HR <- exp(est); lo <- exp(est - 1.96*se); hi <- exp(est + 1.96*se)
    c(HR = HR, CI_low = lo, CI_high = hi)
  }
  
  # 組み立て
  out <- rbind(
    c("CKD: no", "Recovery vs nonAKD",      comp_one(b_rec,  int_rec,  0)),
    c("CKD: no", "Non-Recovery vs nonAKD",  comp_one(b_nrec, int_nrec, 0)),
    c("CKD: yes","Recovery vs nonAKD",      comp_one(b_rec,  int_rec,  1)),
    c("CKD: yes","Non-Recovery vs nonAKD",  comp_one(b_nrec, int_nrec, 1))
  )
  
  tibble::tibble(
    Subgroup   = out[,1],
    Comparison = out[,2],
    HR         = as.numeric(out[,3]),
    CI_low     = as.numeric(out[,4]),
    CI_high    = as.numeric(out[,5])
  )
}

# 実行
hr_int <- hr_int_table(cox_int)
print(hr_int)

##fit_ckd1のcoxを走らせるとモデル不安定のアラートが出現するので下記対応がお勧めと、geminiの回答
#報告すべき主要な結果: 交互作用が有意でないため、cox_base（交互作用なしのモデル）の結果を主たる結果として報告します。
#これは「AKDの3群分類が死亡に与える影響は、CKDの有無によって統計的に有意に異なるとは言えない」という検定結果を反映した、最も頑健な結論です。
#hr_int_table の扱い: この結果は「記述的な推定値」または「参考値」として扱います。
#交互作用が有意でない以上、CKDあり/なしでハザード比が異なるという積極的な主張はできません。
#論文の補足資料（Supplement）に載せるか、本文中で「交互作用は有意ではなかったが、参考までに各層でのハザード比を示すと...」と記述する程度に留めるのが一般的です。
#サブグループ別解析（fit_ckd1）の扱い: fit_ckd1で出た警告は、データを分割したことでモデルが不安定になった明確な証拠です。
#したがって、このサブグループ別解析の結果（特にfit_ckd1）は信頼性が低く、主要な結果として採用すべきではありません。
#交互作用モデルのアプローチが、この不安定性を回避するためのより優れた方法であったことが示された、と解釈できます。