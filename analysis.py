"""DOC-2-077 frozen analysis (PROTOCOL.md lock-1). Run once. Needs per_assay.csv, exposures.tsv."""
import json, numpy as np, pandas as pd
from scipy.stats import rankdata
rng = np.random.default_rng(12345)
a = pd.read_csv("per_assay.csv"); e = pd.read_csv("exposures.tsv", sep="\t").drop(columns=["seq_len"])
d = a.merge(e, on="DMS_id"); d = d[d.rho.notna() & (d.n_scored >= 50)].copy()
d["x"] = np.log(d.seq_len); d["tsub"] = d.DMS_id.str.contains("Tsuboyama_2023").astype(float)
disc, ext = d[d.year <= 2021].copy(), d[d.year >= 2022].copy()
def resid(v, X):
    X1 = np.column_stack([np.ones(len(X)), X]); b = np.linalg.lstsq(X1, v, rcond=None)[0]; return v - X1 @ b
def covs(df, use_tsub=True, use_wtll=True):
    c = [(df.taxon == t).astype(float).values for t in ["Eukaryote", "Prokaryote", "Virus"] if (df.taxon == t).any()]
    if use_tsub and df.tsub.nunique() > 1: c.append(df.tsub.values)
    if use_wtll: c.append(rankdata(df.wtll))
    return np.column_stack(c) if c else np.zeros((len(df), 0))
def partial(df, **k):
    X = covs(df, **k); rx, ry = rankdata(df.x), rankdata(df.rho)
    return np.corrcoef(resid(rx, X), resid(ry, X))[0, 1]
def boot(df, n=10000, **k):
    prots = df.UniProt_ID.unique(); grp = {p: df[df.UniProt_ID == p] for p in prots}; out = []
    for _ in range(n):
        s = pd.concat([grp[p] for p in rng.choice(prots, len(prots))])
        if s.x.nunique() > 2: out.append(partial(s, **k))
    return list(np.percentile(out, [2.5, 97.5]))
res = {}
def gate(name, df, thr, ref_sign=None, **k):
    p = partial(df, **k); ci = boot(df, **k)
    ok = bool(p >= thr and ci[0] > 0 and (ref_sign is None or np.sign(p) == ref_sign))
    res[name] = dict(n=len(df), n_proteins=int(df.UniProt_ID.nunique()), partial_rho=p, ci=ci, pass_=ok); return p
g1 = gate("G1", disc, 0.15)
gate("G2", ext, 0.10, ref_sign=np.sign(g1))
oth = d[d.tsub == 0]; gate("G3", oth, 0.10, use_tsub=False)
res["posthoc_noWTLL"] = dict(G1=partial(disc, use_wtll=False), G2=partial(ext, use_wtll=False), label="POST-HOC, cannot change label")
res["descriptive"] = dict(n_total=len(d), n_tsub=int(d.tsub.sum()), raw_spearman_all=float(np.corrcoef(rankdata(d.x), rankdata(d.rho))[0, 1]))
p = [res[g]["pass_"] for g in ["G1", "G2", "G3"]]
res["LABEL"] = "POSITIVE" if all(p) else ("CONFOUNDED" if p[0] and p[1] else "HONEST NEGATIVE")
print("RESULT_JSON", json.dumps(res, default=float)); open("results.json", "w").write(json.dumps(res, default=float, indent=1))
