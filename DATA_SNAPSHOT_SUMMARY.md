# Data snapshot summary (from `seu_with_nk_flags_meta.data.csv`)

Computed with a lightweight Python CSV scan (no pandas dependency).

- Total cells: **158,978**
- Infections (cells): UI 74,685; Lm 28,237; Tg 28,098; Pf 27,958
- Timepoints (cells): 24h 113,028; 48h 45,950
- Donors: 10
- Top cell types by cell count: HBC 40,157; F 24,059; Endo_f 23,505; VCT_fusing 19,584; PV 16,615

NK/cytotoxic flags:
- `NK_like_strict = TRUE`: 12 / 158,978 (~0.0075%)
- `cytotoxic_like = TRUE`: 7,100 / 158,978 (~4.47%)
- `immune_cytotoxic_like = TRUE`: 4,138 / 158,978 (~2.60%)

QC metric quantiles:
- `nCount_RNA`: median 18,160; q01 3,416.77; q99 105,114.84
- `nFeature_RNA`: median 4,356; q01 1,558; q99 8,965
- `percent.mt`: median 0.0393; q01 0.00345; q99 0.1556

Sample-level exploratory correlations (AdhesionScore vs InnateScore):
- Lm_24h: n=3, r=0.097
- Lm_48h: n=2, r=-1.000
- Pf_24h: n=4, r=-0.414
- Pf_48h: n=2, r=-1.000
- Tg_24h: n=4, r=0.859
- Tg_48h: n=4, r=0.424
- UI_24h: n=10, r=0.350
- UI_48h: n=7, r=0.005

Interpretation:
- Correlations with n=2 are not stable/inferentially useful and should be treated as descriptive only.
