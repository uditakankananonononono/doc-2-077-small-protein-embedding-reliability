# DOC-2-077 "The Small Protein Problem" - frozen protocol (lock-1)

Written 2026-10-09 IST, committed BEFORE any comparison of per-assay rho against sequence length has been opened.

Question: does zero-shot ESM-2 variant-effect accuracy (per-assay Spearman rho between the model score and the measured DMS score) fall for small proteins, beyond the assay-type confound?

## Data (all inputs frozen, from the DOC-2-073 repo, same run)
- per_assay.csv: ESM-2 650M wild-type-marginal scores vs ProteinGym v1 DMS (Zenodo record 15293562), 204 assays, CPU fp32. Columns: DMS_id, n_scored, n_mismatch, rho, wtll, seq_len.
- exposures.tsv: UniProt_ID, taxon, year, seq_len etc. from DOC-2-073 (UniProt release 2026_03).
- md5 of the frozen inputs: INPUT_HASHES.md5.
- Disclosure: these rho values already existed from DOC-2-073. They have NOT been examined against sequence length or assay class before this lock. Their summary appeared only in the DOC-2-073 outputs (RESULTS.md), which do not report rho by length.

## Definitions
- x = log(seq_len). y = per-assay rho. Assays with n_scored < 50 are dropped (as in 073).
- Assay class: Tsuboyama = DMS_id contains "Tsuboyama_2023" (stability); Other = everything else.
- Discovery = year <= 2021. External (held out) = year >= 2022. Same split as 073.
- Partial Spearman: rank(x) and rank(y), each residualized (OLS) on covariates, then Pearson correlation.
- Covariates: Tsuboyama indicator, taxon dummies (Human reference; Eukaryote, Prokaryote, Virus), rank(WTLL).
- CI: protein-cluster bootstrap (cluster = UniProt_ID), 10,000 resamples, seed 12345, 95% percentile.

## Gates
- G1 (discovery): partial rho >= +0.15 and CI lower bound > 0.
- G2 (external, evaluated once): same sign as G1, partial rho >= +0.10, CI lower bound > 0.
- G3 (robustness): within Other-class assays only (discovery + external pooled), covariates taxon dummies + rank(WTLL): partial rho >= +0.10 and CI lower bound > 0.
- G4 (optional second model): ESM-2 150M scores of the same assays reproduce the sign of G1 and G2. If CPU time does not allow, the report says "not run". No partial G4.

## Labels (decided mechanically by analysis.py)
- POSITIVE: G1, G2 and G3 all pass.
- CONFOUNDED: G1 and G2 pass, G3 fails (effect carried by stability assays).
- HONEST NEGATIVE: otherwise (including any failing G1 or G2).
- Sensitivity: a run with WTLL removed from the covariates is reported but labelled POST-HOC; it can never change the label.

## Rules
- No simulated data, no stubs. Every number comes from analysis.py run once on the frozen inputs.
- Deviation budget: one tolerance line. Any amendment is a new dated AMENDMENT-N.md committed before outcomes exist.
- Known limits stated up front: protein-level clustering only (no family clustering); wild-type-marginal scoring; length is correlated with assay type and taxon; observational only.
