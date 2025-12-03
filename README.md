[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=dnafinder/myfriedman&file=myfriedman.m)

# myfriedman

Friedman test for nonparametric two-way ANOVA with optional post-hoc multiple comparisons, implemented in MATLAB.

Compared to MATLAB’s built-in `friedman`, this function:
- uses exact critical values for small block–treatment designs (via tables),
- reports both chi-square and F approximations for larger designs,
- optionally performs post-hoc multiple comparisons between treatments when the global test is significant.

## Syntax

- [stats, mc] = myfriedman(X)
- [stats, mc] = myfriedman(X, 'Alpha', ALPHA, 'Reps', REPS, ...
                                'PostHoc', POSTHOC, 'Display', DISPLAY, ...
                                'Exact', EXACT)

## Inputs

### Required

- X  
  Numeric matrix of size (Blocks × Treatments), real, finite, non-NaN, non-empty.

### Name–Value options

- 'Alpha'  
  Significance level (scalar in (0,1), default 0.05).

- 'Reps'  
  Number of replicate observations per block–treatment cell (positive integer, default 1).  
  If Reps > 1, rows of X are interpreted as stacked replicates for each block.

- 'PostHoc'  
  Logical-like flag to enable/disable post-hoc multiple comparisons (default true).  
  When true and the global Friedman test rejects H0, pairwise comparisons between treatments are computed.

- 'Display'  
  Logical-like flag to control command-window output (default true).  
  If false, no text is printed; only stats and mc are returned.

- 'Exact'  
  Logical-like flag to use exact critical values from tables for small designs (default true).  
  When true, exact critical values are used for specific (blocks, treatments, alpha) combinations; otherwise, the chi-square/F approximations are always used.

For 'PostHoc', 'Display', and 'Exact', the following are accepted as true/false:  
true/false, 1/0, "on"/"off", "yes"/"no", "true"/"false" (case-insensitive).

When using exact critical values, Alpha must be one of: 0.1, 0.05, 0.025, 0.01, 0.005, 0.001.

## Outputs

### STATS structure

- stats.blocks        – number of blocks (effective rows)  
- stats.treatments    – number of treatments (columns)  
- stats.reps          – replicates per cell  
- stats.sumRanks      – K×1 vector of sums of ranks per treatment  
- stats.sumSqRanks    – sum of squared ranks  
- stats.Fr            – Friedman statistic (unscaled Tx)  
- stats.ties          – logical flag indicating presence of ties  
- stats.tieCorrection – correction factor used for ties (if any)  
- stats.alpha         – significance level used  
- stats.method        – 'Exact Friedman', 'Chi-square/F approximation', or similar  
- stats.chi2.value    – chi-square approximation (or NaN if not used)  
- stats.chi2.df       – chi-square degrees of freedom  
- stats.chi2.p        – p-value for chi-square approximation  
- stats.F.value       – F approximation (or NaN if not used)  
- stats.F.df_num      – numerator degrees of freedom  
- stats.F.df_denom    – denominator degrees of freedom  
- stats.F.p           – p-value for F approximation  
- stats.criticalValue – exact critical value (for small designs, otherwise NaN)  
- stats.rejectNull    – true if the treatments do not have identical effects  

### MC structure (post-hoc)

Returned when 'PostHoc' is true and stats.rejectNull is true; otherwise empty.

- mc.method      – 'Exact pairwise (Bioinformatics 2017)' or 'Conover-type LSD'  
- mc.Rdiff       – K×K matrix of absolute differences among rank sums  
- mc.pvalue      – K×K matrix of pairwise p-values (exact method only; empty otherwise)  
- mc.cv          – critical value for differences (Conover-type method only; empty otherwise)  
- mc.significant – logical K×K matrix (lower triangle), true where differences are significant at the chosen alpha  

## Method (breve)

- Each block is ranked across treatments (and replicates if present) using tiedrank, yielding a rank matrix R.  
- The sum of ranks for each treatment is computed, along with the sum of squared ranks.  
- For small designs and suitable (blocks, treatments, alpha), exact critical values are taken from tables: friedman.A/B/C/D in myfriedmantables.mat.  
- For larger designs, a chi-square approximation (with optional tie correction) and a corresponding F approximation are computed.  
- If the F-based p-value is below Alpha, the global null hypothesis of identical treatment effects is rejected.  
- When the global test is significant and 'PostHoc' is true:  
  - If all rank-sum differences are integers, an exact pairwise procedure (Bioinformatics 2017) is used to compute p-values.  
  - Otherwise, a Conover-type nonparametric LSD approach is used, comparing rank-sum differences to a critical value based on the t distribution.

## Example

    x = [115 142  36  91  28;
          28  31   7  21   6;
         220 311 108  51 117;
          82  56  24  46  33;
         256 298 124  46  84;
         294 322 176  54  86;
          98  87  55  84  25];

    [stats, mc] = myfriedman(x);

## Citation

If you use this function in a scientific publication, please cite it as:

Cardillo G. (2009). MYFRIEDMAN: Friedman test for non parametric two way ANalysis Of VAriance. Available on GitHub: https://github.com/dnafinder/myfriedman
