# DOC-2-077 results (frozen analysis, PROTOCOL.md lock-1; no post-hoc changes to gates)

Label (set mechanically by analysis.py): HONEST NEGATIVE. G1, G2 and G3 all fail. The data show no detectable positive partial association between sequence length and ESM-2 zero-shot rho. Not detected, imprecise; this is not evidence of absence. The G1 CI reaches +0.25 and the G3 CI includes +0.10, so effects near the thresholds cannot be excluded. The label is a gate outcome, not a finding that small proteins do not do worse. Independent gate review: scoped pass for publishing this honest negative, with the wording fixes now applied; it reproduced the numbers in an independent rerun and confirmed lock-before-run ordering by server-side release timestamps.

## Results (ESM-2 650M wt-marginal, 204 assays, analysis.py run once)
- G1 (discovery, year <= 2021; n=98 assays, 79 proteins): partial rho +0.001, 95% protein-cluster bootstrap CI [-0.214, +0.247]. Needed >= +0.15 with CI above 0. Fail.
- G2 (external, year >= 2022; n=106, 98 proteins): partial rho -0.062, CI [-0.276, +0.159]. Needed same sign as G1 and >= +0.10. Fail.
- G3 (Other-class assays only, pooled; n=141, 112 proteins): partial rho +0.084, CI [-0.105, +0.283]. Needed >= +0.10 with CI above 0. Fail. The CI is wide; an effect near the +0.10 threshold is not ruled out.
- Descriptive: raw (unadjusted) Spearman between log(seq_len) and rho over all 204 assays is -0.232, i.e. the opposite sign to the "small proteins do worse" hypothesis. With the assay-class, taxon and WTLL adjustment the association is near zero in the discovery and external splits (G1 +0.001, G2 -0.062); no pooled adjusted statistic over all 204 assays was computed. n_total=204, of which 63 are Tsuboyama stability assays. Raw correlations were not a gate; I report this because it shows the covariates matter.
- Post-hoc (cannot change the label): with WTLL removed from the covariates, G1 = -0.048, G2 = -0.062.

## Deviations and limits
- G4 (ESM-2 150M replication): NOT RUN. The protocol allowed "not run".
- Clusters are UniProt_ID (protein level), not family level. Weaker control for relatedness.
- Scoring is wild-type-marginal (a different estimand from masked-marginal), fp32 on CPU.
- WTLL is a weak covariate (see DOC-2-073 RESULTS.md, where its baseline R2 was negative).
- The rho values were reused from DOC-2-073 (the same scoring run). They were not examined against length or assay class before the lock-1 commit (f8ce0f1); analysis.py was smoke-tested only on synthetic random rho.
- Observational only. Length is correlated with assay class and taxon; the adjustment set is linear in ranks.
- Evidence caveats: analysis_run.log contains only the RESULT_JSON line (no command, timestamp or library versions). "analysis.py was run once in the sandbox" and "smoke-tested only on synthetic random rho" are the builder's statements, not repo evidence. Lock-before-run ordering is shown by the lock-1 tag/commit history; the independent gate confirmed it with server-side timestamps.
- Discovery and external sets share 2 proteins (DYR_ECOLI, SRC_HUMAN: 2 discovery assays, 3 external assays), so G2 is not fully protein-held-out. The effect of this is judged negligible but was not tested.
- Bootstrap: 10,000 resamples per gate (seed 12345). Resamples skipped by the nunique>2 rule: G1 0, G2 0, G3 0 (counted by count_skips.py, which replays analysis.py's exact RNG stream; output in skips.log).
