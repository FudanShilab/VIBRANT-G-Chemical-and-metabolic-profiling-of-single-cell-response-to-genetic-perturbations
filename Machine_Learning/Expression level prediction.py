"""
Script Name: expression_level_prediction.py

Description: Performs 1D binary LDA (WT vs KD1) on single-cell FTIR features to
            calibrate distance-to-WB mapping and predict KD2_WB (numeric-only).

Corresponding paper:
"VIBRANT-G: Chemical and metabolic profiling of single-cell response
 to genetic perturbations"
Authors: Minqian Wei, Minjie Fu, Tongqi Wang, Yuchen Sun, Mingyu Wang,
         Xinyuan Bi, Wei Hua, Lixue Shi*

Dependencies:
- Python 3.8+
- numpy, pandas, scikit-learn, matplotlib, seaborn

License:
For academic use only.
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

CSV_PATH      = "/Users/weiminqian/Downloads/machine_learning/unsupervised model/demodata/data_expression_level_prediction/U87 TKT.csv"    # <-- change to your file path
GROUP_COLNAME = "group"

# Map raw group names to internal labels (WT / KD1 / KD2)
GROUP_NAME_MAP = {
    "WT"   : "WT",
    "KD1" : "KD1",
    "KD2" : "KD2",
}

# Western blot values (relative to actin)
WT_WB  = 1.0
KD1_WB = 0.052  # T98G TKT KD1 = 0.1460; T98G SLC2A1 KD1 = 0.1973; U87 TKT KD1 = 0.05226; U87 SLC2A1 KD1 = 0.074;

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
    """Compute pooled within-class variance (WT + KD1 only)."""
    classes = np.unique(labels)
    N, p = values.shape
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
    merged_map = {"WT":"WT","KD1":"KD1","KD2":"KD2", **GROUP_NAME_MAP}
    groups_internal = [merged_map[g] for g in groups_original]

    keep_mask = np.array([g in {"WT","KD1","KD2"} for g in groups_internal], bool)
    df_feats = df.drop(columns=[GROUP_COLNAME])
    X_all = to_numpy(df_feats)[keep_mask]
    map_name_to_id = {"WT":0,"KD1":1,"KD2":2}
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

    # === Binary LDA (WT vs KD1) ===
    mask2 = (y==0) | (y==1)
    lda = LinearDiscriminantAnalysis(n_components=1)
    lda.fit(X_for_lda[mask2], y[mask2])

    z_all = lda.transform(X_for_lda).ravel()
    z_wt   = float(z_all[y==0].mean())
    z_kd1  = float(z_all[y==1].mean())
    z_kd2  = float(z_all[y==2].mean())

    # Pooled variance from WT + KD1
    var1d = pooled_variance_1d(z_all[mask2,None], y[mask2])
    s1d   = np.sqrt(var1d) if var1d > 1e-12 else 1.0

    D1_WT_KD1 = abs(z_wt - z_kd1) / s1d
    D1_WT_KD2 = abs(z_wt - z_kd2) / s1d

    # Calibration: ΔWB = a * D
    dWB1 = KD1_WB - WT_WB
    KD2_WB_pred = np.nan if D1_WT_KD1<=1e-12 else WT_WB + dWB1/D1_WT_KD1 * D1_WT_KD2

    # === Console summary ===
    print("\n== 1D LDA Calibration Results ==")
    print(f"D_1D(WT, KD1) = {D1_WT_KD1:.6f}")
    print(f"D_1D(WT, KD2) = {D1_WT_KD2:.6f}")
    print(f"WT_WB = {WT_WB:.6f}, KD1_WB = {KD1_WB:.6f}")
    print(f"Predicted KD2_WB = {KD2_WB_pred:.6f}")

    # === CSV export ===
    if SAVE_CSV:
        summary = pd.DataFrame({
            "metric":[
                "D_1D(WT,KD1)",
                "D_1D(WT,KD2)",
                "WT_WB",
                "KD1_WB",
                "KD2_WB_pred"
            ],
            "value":[
                D1_WT_KD1, D1_WT_KD2,
                WT_WB, KD1_WB, KD2_WB_pred
            ]
        })
        path = out_path / "summary.csv"
        summary.to_csv(path, index=False)
        print("\nSaved CSV:", path)

if __name__ == "__main__":
    main()
