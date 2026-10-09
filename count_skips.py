"""Counts bootstrap resamples skipped by the nunique>2 rule, replaying analysis.py's exact RNG stream (partial() consumes no RNG)."""
import numpy as np, pandas as pd
rng = np.random.default_rng(12345)
a = pd.read_csv("per_assay.csv"); e = pd.read_csv("exposures.tsv", sep="\t").drop(columns=["seq_len"])
d = a.merge(e, on="DMS_id"); d = d[d.rho.notna() & (d.n_scored >= 50)].copy()
d["x"] = np.log(d.seq_len); d["tsub"] = d.DMS_id.str.contains("Tsuboyama_2023").astype(float)
disc, ext = d[d.year <= 2021].copy(), d[d.year >= 2022].copy()
def cnt(df, n=10000):
    prots = df.UniProt_ID.unique(); grp = {p: df[df.UniProt_ID == p].x.values for p in prots}; sk = 0
    for _ in range(n):
        s = np.concatenate([grp[p] for p in rng.choice(prots, len(prots))])
        if len(np.unique(s)) <= 2: sk += 1
    return sk
# order must match analysis.py: G1 boot(disc), G2 boot(ext), G3 boot(oth)
out = {"G1": cnt(disc), "G2": cnt(ext), "G3": cnt(d[d.tsub == 0])}
ov = sorted(set(disc.UniProt_ID) & set(ext.UniProt_ID))
out["overlap_proteins"] = ov; out["n_overlap"] = len(ov)
out["assays_in_overlap_disc"] = int(disc.UniProt_ID.isin(ov).sum()); out["assays_in_overlap_ext"] = int(ext.UniProt_ID.isin(ov).sum())
print(out)
