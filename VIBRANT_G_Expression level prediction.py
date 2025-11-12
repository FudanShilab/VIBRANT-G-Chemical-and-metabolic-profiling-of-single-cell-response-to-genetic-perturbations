"""
=============================================================================
Copyright (c) 2025 Lixue Shi and collaborators
Licensed under the Apache License, Version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)

Script Name: VIBRANT_G_Expression level prediction.py

Description:
Performs 1D binary LDA (WT vs reference) on single-cell FTIR features to
calibrate distance-to-WB mapping and predict the target WB value (numeric-only).

Corresponding paper:
"VIBRANT-G: Chemical and metabolic profiling of single-cell response
 to genetic perturbations"
Authors: Minqian Wei, Minjie Fu, Tongqi Wang, Yuchen Sun, Mingyu Wang,
         Xinyuan Bi, Wei Hua, Lixue Shi*

Dependencies:
- Python 3.8+
- numpy, pandas, scikit-learn, matplotlib, seaborn

=============================================================================
"""

from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.decomposition import PCA

# =========================
# ===== USER CONFIG =======
# =========================

CSV_PATH      = "/path/to/demo_data.csv"    # <-- change to your file path
GROUP_COLNAME = "group"

# Map raw group names to internal labels (WT / reference / target)
GROUP_NAME_MAP = {
    "WT": "WT",
    "known expression levels": "reference",
    "expression levels to be predicted": "target",
}

# Western blot values (relative to actin)
WT_WB         = 1.0
REFERENCE_WB  = 0.074  # for demo data: # T98G TKT referencce = 0.1460; T98G SLC2A1 referencce = 0.1973; U87 TKT referencce = 0.05226; U87 SLC2A1 referencce = 0.074;

# Target WB will be predicted

# Preprocessing
USE_PCA_BEFORE_LDA   = True
PCA_VARIANCE_TO_KEEP = 0.99
RANDOM_SEED          = 42

# Output
OUT_DIR   = "lda_1d_results"
SAVE_CSV  = True

# =========================
# ====== HELPERS ==========
# =========================

def ensure_outdir(path: Path):
    path.mkdir(parents=True, exist_ok=True)

def to_numpy(df_or_array) -> np.ndarray:
    return (df_or_array.values if isinstance(df_or_array, pd.DataFrame) else np.asarray(df_or_array)).astype(float)

def pooled_variance_1d(values: np.ndarray, labels: np.ndarray) -> float:
    """Compute pooled within-class variance (WT + reference only)."""
    classes = np.unique(labels)
    N, _ = values.shape
    pooled = 0.0
    for c in classes:
        Zc = values[labels == c]
        if Zc.shape[0] > 1:
            pooled += (Zc.shape[0] - 1) * np.var(Zc, ddof=1)
    denom = max(1, N - len(classes))
    return pooled / denom

# =========================
# ========== MAIN =========
# =========================

def main():
    np.random.seed(RANDOM_SEED)
    out_path = Path(OUT_DIR).resolve()
    ensure_outdir(out_path)

    # === Load data ===
    df = pd.read_csv(CSV_PATH)
    if GROUP_COLNAME not in df.columns:
        raise ValueError(f"Group column '{GROUP_COLNAME}' not found in CSV")

    groups_original = df[GROUP_COLNAME].astype(str).tolist()
    merged_map = {"WT": "WT", "reference": "reference", "target": "target", **GROUP_NAME_MAP}
    groups_internal = [merged_map[g] for g in groups_original]

    keep_mask = np.array([g in {"WT", "reference", "target"} for g in groups_internal], bool)
    df_feats = df.drop(columns=[GROUP_COLNAME])
    X_all = to_numpy(df_feats)[keep_mask]
    # 0: WT, 1: reference, 2: target
    map_name_to_id = {"WT": 0, "reference": 1, "target": 2}
    y = np.array([map_name_to_id[g] for g in np.array(groups_internal)[keep_mask]], int)

    # === Preprocess ===
    scaler = StandardScaler()
    Xz = scaler.fit_transform(X_all)
    if USE_PCA_BEFORE_LDA:
        pca = PCA(n_components=None, random_state=RANDOM_SEED)
        Z = pca.fit_transform(Xz)
        cum = np.cumsum(pca.explained_variance_ratio_)
        n_keep = int(np.searchsorted(cum, PCA_VARIANCE_TO_KEEP) + 1)
        X_for_lda = Z[:, :n_keep]
    else:
        X_for_lda = Xz

    # === Binary LDA (WT vs reference) ===
    mask2 = (y == 0) | (y == 1)
    lda = LinearDiscriminantAnalysis(n_components=1)
    lda.fit(X_for_lda[mask2], y[mask2])

    z_all = lda.transform(X_for_lda).ravel()
    z_wt   = float(z_all[y == 0].mean())
    z_ref  = float(z_all[y == 1].mean())
    z_tgt  = float(z_all[y == 2].mean())

    # Pooled variance from WT + reference
    var1d = pooled_variance_1d(z_all[mask2, None], y[mask2])
    s1d   = np.sqrt(var1d) if var1d > 1e-12 else 1.0

    D1_WT_REF = abs(z_wt - z_ref) / s1d
    D1_WT_TGT = abs(z_wt - z_tgt) / s1d

    # Calibration: ΔWB = a * D, where a = (REFERENCE_WB - WT_WB) / D1_WT_REF
    dWB_ref = REFERENCE_WB - WT_WB
    TARGET_WB_pred = np.nan if D1_WT_REF <= 1e-12 else WT_WB + dWB_ref / D1_WT_REF * D1_WT_TGT

    # === Console summary ===
    print("\n== 1D LDA Calibration Results ==")
    print(f"D_1D(WT, reference) = {D1_WT_REF:.6f}")
    print(f"D_1D(WT, target)    = {D1_WT_TGT:.6f}")
    print(f"WT_WB = {WT_WB:.6f}, REFERENCE_WB = {REFERENCE_WB:.6f}")
    print(f"Predicted TARGET_WB = {TARGET_WB_pred:.6f}")

    # === CSV export ===
    if SAVE_CSV:
        summary = pd.DataFrame({
            "metric": [
                "D_1D(WT,reference)",
                "D_1D(WT,target)",
                "WT_WB",
                "REFERENCE_WB",
                "TARGET_WB_pred"
            ],
            "value": [
                D1_WT_REF, D1_WT_TGT,
                WT_WB, REFERENCE_WB, TARGET_WB_pred
            ]
        })
        path = out_path / "summary.csv"
        summary.to_csv(path, index=False)
        print("\nSaved CSV:", path)

if __name__ == "__main__":
    main()
