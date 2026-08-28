# KM解析とeGFR slope解析のサンプルサイズ不一致 — 解決ガイド

> 対象: eGFR slope landmark210解析の修正担当者
> 作成日: 2026-04-22

---

## 1. 問題の構造: なぜNが合わないのか

KM（カプラン・マイヤー）解析とeGFR slope解析では、解析に必要なデータの性質が根本的に異なる。

```
KM解析（生存解析）
  必要なもの: 患者ごとに「生存/死亡」の1つのイベント情報 + フォローアップ期間
  → eGFR測定値は不要
  → 対象: index_date+210以降に生存確認ができる全患者

eGFR slope解析（縦断解析）
  必要なもの: 患者ごとに複数回のeGFR測定値 + ベースラインeGFR
  → 最低でもベースライン + 1回以上のeGFR測定が必要
  → 対象: eGFRが繰り返し測定されている患者のみ
```

### 現状のN関係

```
KM N  <  オリジナルslope N  >  landmark210 slope N
```

**KM N < オリジナルslope N となる理由**:
オリジナルslope解析は動的time0（例: index_date+120日にmax eGFR）を使う。
ある患者がindex_date+120日にeGFR測定 → +130日にも測定あり → slope解析に寄与。
しかし+200日で通院中断 → index_date+210のランドマーク時点で追跡できない → KMでは除外（`time_years < 0`）。

### N不一致の原因まとめ

| # | 原因 | KM | オリジナルslope | landmark210 slope |
|---|------|-----|----------------|-------------------|
| 1 | ランドマーク要件 | index_date+210以降の生存確認が必要 | 不要（動的time0） | index_date+210以降のeGFR測定が必要 |
| 2 | ベースラインeGFR | 不要 | 必須（90-210日 or 0-90日） | 必須（150-210日のみ、狭い） |
| 3 | 測定回数 | イベント情報1つで十分 | ≥2回のeGFR | ≥2回のeGFR |

- KMから脱落するがslopeに入る患者: day 210前に追跡不能だがday 90-210にeGFR測定あり
- slopeから脱落するがKMに入る患者: day 210以降に生存確認ありだがeGFR測定が不十分
- landmark210はさらにベースラインeGFR取得窓（150-210日）が狭いためN最小

**結論**: N完全一致は原理的に不可能（KM=イベント解析 vs slope=縦断解析で、ランドマーク要件も異なる）。以下の3つの戦略で対処する。

---

## 2. 戦略1: 主解析のN差を透明に報告する（コード変更なし）

### やること
論文のMethods/Supplementに以下の趣旨を記載する。コード変更は不要。

### 記載文言案（英語、Methods向け）

> The eGFR slope analysis included patients who had a valid baseline eGFR value and at least one subsequent eGFR measurement after the baseline timepoint. Sample sizes for the eGFR slope analysis therefore differ from those of the survival analysis, which requires only follow-up status (alive or dead) after the landmark at index date + 210 days.

### 記載文言案（英語、Supplemental Table向け）

> Supplemental Table X shows the number of patients included in each analysis by group. Differences in sample size reflect the distinct data requirements: the Kaplan-Meier analysis requires only event status, whereas the eGFR slope analysis requires repeated eGFR measurements.

---

## 3. 戦略2: Concordant slope感度分析（slope解析をKM患者に限定）

### 概要

現状は **KM患者 ⊆ slope患者**（KMはslopeの部分集合）である。

```
┌─────────────────────────────────────┐
│          slope解析対象者             │
│                                     │
│   ┌─────────────────────┐           │
│   │   KM解析対象者      │           │
│   │  (= slope ∩ KM)     │           │
│   │                     │           │
│   └─────────────────────┘           │
│                                     │  ← この差分: day 210前に追跡不能だが
│                                     │    eGFR測定はある患者（slope only）
└─────────────────────────────────────┘
```

**なぜ KM ⊆ slope か**:
- KMは `last_follow_death >= index_date + 210`（ランドマーク生存）を要求
- slopeはランドマーク要件なし（動的time0で、index_dateまでフォールバック）
- KMはeGFRを不要、slopeのtime0_egfrも階層的フォールバックでほぼ全員取得可能
- → KMに入れる条件はslopeより厳しい → KM患者は全員slope解析にも含まれる

**このため**:
- ~~Concordant KM~~（KMをslope患者に限定）= オリジナルKMと同一 → **不要**
- **Concordant slope**（slopeをKM患者に限定）= 唯一意味のある感度分析

### この感度分析で示すこと
「day 210以降も追跡可能な患者（= KM集団）に限定しても、eGFR slope解析の結果は変わらない」

### 重要ポイント
- **slopeのtime0は動的time0のまま変更なし**（解析手法は変えない、対象集団だけ絞る）
- **選択バイアス**: KM集団に限定 = day 210以降も追跡可能な患者に限定。早期追跡不能（通院中断・転院等）の患者が除外されるため、比較的安定した通院歴の集団に偏る可能性がある。Limitationsで言及が必要。

### 修正箇所A: KM対象IDの保存

**ファイル**: `recvsnon_recvsAKD.R`
**場所**: L148の直後（メインのKM `dat_km` 確認の直後）
**理由**: KM解析に含まれる患者のIDリストを保存し、slope側で使えるようにする

```r
# ========================================
# Save KM cohort IDs for concordant slope analysis
# ========================================
km_ids <- dat_km %>%
  distinct(id, group)

write_csv(km_ids, file.path(getwd(), "km_cohort_ids.csv"))
cat("KM cohort IDs saved to km_cohort_ids.csv\n")
cat("KM N per group:\n")
print(table(km_ids$group, useNA = "ifany"))
```

### 修正箇所B: Concordant slope感度分析の実行

**ファイル**: `eGFRslope_confounding.R`
**場所**: L253の直後（オリジナルslope解析の `sample_sizes` 計算の直後）
**理由**: KM解析対象者に限定してslope解析を再実行し、結果の頑健性を確認する

```r
# ==========================================================
# Sensitivity: Concordant slope
# (restrict slope analysis to KM-eligible patients,
#  i.e., patients who survived past the day-210 landmark)
# ==========================================================
km_ids <- read_csv(file.path(getwd(), "km_cohort_ids.csv"))

longdat_concordant <- longdat %>%
  filter(id %in% km_ids$id)

cat("\n=== Concordant slope (restricted to KM patients) ===\n")
cat("Full slope N:\n")
print(longdat %>% distinct(id, jin_label) %>% count(jin_label))
cat("\nConcordant slope N (= KM N):\n")
print(longdat_concordant %>% distinct(id, jin_label) %>% count(jin_label))

# Re-fit LME on concordant cohort (same model specification)
# time0 remains dynamic (unchanged from original slope analysis)
# Example for ≤1 year:
# fit_concordant_1y <- lme(
#   egfr ~ years_from_time0 * jin_label + time0_egfr_c - 1 +
#     age_c + sex + arb_acei_use +
#     dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15 +
#     years_from_time0:(age_c + arb_acei_use +
#       dn1 + dn3 + dn4 + dn5 + dn6 + dn7 + dn8 + dn9 + dn10 + dn12 + dn13 + dn14 + dn15),
#   random = list(id = pdSymm(~ 1 + years_from_time0)),
#   data = longdat_concordant %>% filter(years_from_time0 <= 1),
#   na.action = na.omit, method = "REML", control = ctrl
# )
```

### 実行順序

1. KM解析を実行 → `km_cohort_ids.csv` を保存（修正箇所A）
2. slope解析を実行（オリジナルのまま）
3. concordant slope感度分析を実行（修正箇所B: KM IDを読み込む）

### Limitations記載文言案（英語）

> In the concordant cohort sensitivity analysis, the eGFR slope analysis was restricted to patients who were also eligible for the Kaplan-Meier analysis (i.e., those with confirmed follow-up beyond the 210-day landmark). This restriction excludes patients who were lost to follow-up before day 210, which may introduce selection bias toward patients with more stable clinical engagement.

---

## 4. 戦略3: landmark210 slope解析のベースラインeGFR取得改善

### 概要
現在のlandmark210では **150-210日の60日間** しかベースラインeGFRを探していない。
これを **階層的フォールバック** に変更し、より多くの患者を取り込む。

### 現行ロジック（`eGFRslope_confounding.R` L2009-2024）

```r
# baseline eGFR nearest to day 210 within 150-210 days
time0_egfr_df <- jin1_inclusion %>%
  mutate(
    days_from_index = as.numeric(date - index_date),
    dist_to_210 = abs(days_from_index - 210)
  ) %>%
  filter(days_from_index >= 150, days_from_index <= 210) %>%   # <-- 150-210日のみ
  group_by(id) %>%
  arrange(dist_to_210, desc(date), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    id,
    time0_egfr = egfr,
    time0_egfr_date = date
  )
```

**問題点**: 150-210日にeGFR測定がない患者は全員脱落。

### 改善ロジック（差し替え案）

L2009-2024を以下に差し替える。time0は `index_plus_210` のまま固定（ランドマーク）。
ベースラインeGFR（共変量として使用）の取得元を広げるだけ。

```r
# ========================================
# baseline eGFR: hierarchical fallback
# Priority 1: 150-210 days (nearest to day 210)
# Priority 2:  90-150 days (nearest to day 150)
# Priority 3:   0- 90 days (latest measurement)
# ========================================

# Priority 1: 150-210 days
p1 <- jin1_inclusion %>%
  mutate(days_from_index = as.numeric(date - index_date)) %>%
  filter(days_from_index >= 150, days_from_index <= 210) %>%
  group_by(id) %>%
  arrange(abs(days_from_index - 210), desc(date), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(id, time0_egfr = egfr, time0_egfr_date = date, egfr_source = "150-210d")

# Priority 2: 90-150 days (for patients NOT found in Priority 1)
p2 <- jin1_inclusion %>%
  filter(!id %in% p1$id) %>%
  mutate(days_from_index = as.numeric(date - index_date)) %>%
  filter(days_from_index >= 90, days_from_index < 150) %>%
  group_by(id) %>%
  arrange(abs(days_from_index - 150), desc(date), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(id, time0_egfr = egfr, time0_egfr_date = date, egfr_source = "90-150d")

# Priority 3: 0-90 days (for patients NOT found in Priority 1 or 2)
p3 <- jin1_inclusion %>%
  filter(!id %in% p1$id, !id %in% p2$id) %>%
  mutate(days_from_index = as.numeric(date - index_date)) %>%
  filter(days_from_index >= 0, days_from_index <= 90) %>%
  group_by(id) %>%
  arrange(desc(days_from_index), desc(date), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(id, time0_egfr = egfr, time0_egfr_date = date, egfr_source = "0-90d")

# Combine (one row per patient, highest priority wins)
time0_egfr_df <- bind_rows(p1, p2, p3)

# Diagnostic: how many patients from each source?
cat("\nBaseline eGFR source distribution:\n")
print(table(time0_egfr_df$egfr_source, useNA = "ifany"))
```

### なぜこの改善が有効か

- **time0（ランドマーク）は変えない**: `index_date + 210` のまま（L2028）
- **変えるのはベースラインeGFRの取得元だけ**: これはLMEモデルの共変量（`time0_egfr_c`）として使うもの。測定日がtime0から多少ずれていても、共変量調整として統計的に問題ない。
- **Nが増える**: 150-210日にeGFRがない患者も、90-150日や0-90日の値で救済される。

### L2064のfilterはそのまま維持

```r
filter(!is.na(time0_egfr))   # L2064: これは維持
```

階層的フォールバックによって `time0_egfr` がNAになる患者が減るため、このfilterで脱落する患者数が自動的に減少する。

### 注意: landmark210セクションは2箇所ある

`eGFRslope_confounding.R` 内に同じlandmark210ロジックが **2箇所** ある:
- **L2009-2024**: 1つ目のlandmark210セクション（≤1年解析用）
- **L2475-2490**: 2つ目のlandmark210セクション（別の解析ブロック）

両方とも `filter(days_from_index >= 150, days_from_index <= 210)` を含む。
**両箇所に同じ階層的フォールバックの改善を適用すること。**

---

## 5. 実装の順序

1. **戦略1**: コード変更なし。論文の記載で対応。
2. **戦略3**: `eGFRslope_confounding.R` のlandmark210ベースラインeGFRロジックを改善 → N増加を確認
3. **戦略2**: slope対象IDを保存 → `recvsnon_recvsAKD.R` でconcordant KMを実行

### 確認ポイント
- [ ] 改善後のlandmark210 slope解析のN数を記録（各群: nonAKD, Recovery, Non-Recovery）
- [ ] KM N と slope N を比較し、共通集団のN数を確認
- [ ] concordant slope の結果がフルコホート slope と大きく乖離しないことを確認
- [ ] concordant KM の結果がフルコホート KM と大きく乖離しないことを確認
- [ ] 各解析のN数を一覧表にまとめる（下記テンプレート参照）

### N数比較テーブル（テンプレート）

| 解析 | nonAKD | Recovery | Non-Recovery | 備考 |
|------|--------|----------|--------------|------|
| KM（フルコホート） | | | | 主解析 |
| オリジナルslope（動的time0） | | | | 主解析（KM N ⊆ slope N） |
| concordant slope（KM患者に限定） | | | | 感度分析（N = KM N と一致するはず） |
| landmark210 slope（改善前） | | | | |
| landmark210 slope（改善後） | | | | 戦略3適用後 |

---

## 6. 補足: なぜ完全一致は不要か

学術論文において、KMとeGFR slopeのNが異なることは一般的であり、査読者も理解している。
重要なのは:

1. **N差の理由が説明されていること**（戦略1）
2. **結果のロバスト性が示されていること**（戦略2: concordant cohort）
3. **不必要にNを減らしていないこと**（戦略3: ウィンドウ拡大）

KDIGOガイドライン関連の論文でも、生存解析と腎機能軌跡解析で異なるNを報告するのは標準的な手法である。
