"""
=============================================================================
Script Name: VIBRANT_G_phenotype_classification.py
Description:
1) Train/test split on labeled data; StandardScaler → LDA; report accuracyand normalized confusion matrix.
(2) PCA on the full labeled set; re-train LDA in PCA space.
(3) Predict labels for new spectra in the same PCA space.
(4) Compute Mahalanobis distance from each new sample to class centroidsusing pooled within-class covariance.
(5) IsolationForest-based novelty detection.

Data format:
- X_TRAIN_CSV (train_data.csv): rows = samples, cols = features (numeric)
- Y_TRAIN_CSV (train_label.csv): single-column labels aligned to X_TRAIN rows
- X_NEW_CSV   (test_data.csv):   rows = new samples, cols = features (numeric)
Corresponding paper:
"VIBRANT-G: Chemical and metabolic profiling of single-cell response to genetic perturbations"
 Authors: Minqian Wei, Minjie Fu, Tongqi Wang, Yuchen Sun, Mingyu Wang,Dan Ye, Wei Hua, Lixue Shi

Dependencies:
- Python 3.8+
 - numpy, pandas, scikit-learn, matplotlib, seaborn

License:For academic use only.
"""

from __future__ import annotations
import json, os
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg", force=True)      # 保存图用，无需交互
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, confusion_matrix
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import IsolationForest
from scipy.spatial.distance import mahalanobis


# =========================
#           CONFIG
# =========================
X_TRAIN_CSV   = "/path/to/train_data.csv"  
Y_TRAIN_CSV   = "/path/to/train_label.csv"    
X_NEW_CSV     = "/path/to/test_data.csv"   
OUT_DIR       = "results_mid"

TEST_SIZE     = 0.30
SEED          = 109
PCA_NCOMP     = 50

RUN_IFOREST   = True
IF_CONTAM     = 0.1
IF_N_EST      = 1000

# =========================
#        HELPERS
# =========================
def load_matrix_csv(path: str) -> np.ndarray:
    df = pd.read_csv(path)
    if df.shape[0] == 0 or df.shape[1] == 0:
        raise ValueError(f"Empty CSV: {path}")
    arr = df.values
    if not np.isfinite(arr).all():
        raise ValueError(f"Non-finite values in: {path}")
    return arr

def infer_uniform_priors(y: np.ndarray) -> np.ndarray:
    n = len(np.unique(y))
    return np.ones(n, dtype=float) / float(n)

def pooled_within_covariance(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    classes = np.unique(y)
    N, p = X.shape
    C = len(classes)
    S = np.zeros((p, p))
    for c in classes:
        Xc = X[y == c]
        if Xc.shape[0] > 1:
            Sc = np.cov(Xc, rowvar=False, ddof=1)
            S += (Xc.shape[0]-1) * Sc
    S /= max(1, N - C)
    return S

def class_centroids(X: np.ndarray, y: np.ndarray) -> dict[int, np.ndarray]:
    return {int(c): X[y == c].mean(axis=0) for c in np.unique(y)}

def plot_confusion_norm(cm_norm: np.ndarray, classes: np.ndarray, out_png: Path):
    plt.figure(figsize=(8, 6))
    df_cm = pd.DataFrame(cm_norm, index=classes, columns=classes)
    sns.heatmap(df_cm, annot=True, fmt=".2f", cmap="Blues", cbar=True)
    plt.title("Normalized Confusion Matrix (rows sum to 1)")
    plt.xlabel("Predicted"); plt.ylabel("True")
    plt.tight_layout(); plt.savefig(out_png, dpi=300); plt.close()

# =========================
#          MAIN
# =========================
def main():
    np.random.seed(SEED)
    outdir = Path(OUT_DIR); outdir.mkdir(parents=True, exist_ok=True)

    # ---- load data
    X = load_matrix_csv(X_TRAIN_CSV)     # (N, p)
    y = np.ravel(load_matrix_csv(Y_TRAIN_CSV))  # (N,)
    if X.shape[0] != y.shape[0]:
        raise ValueError(f"Row mismatch: X rows={X.shape[0]} vs y rows={y.shape[0]}")
    X_new = load_matrix_csv(X_NEW_CSV)

    # ---- split & scale
    X_tr, X_te, y_tr, y_te = train_test_split(
        X, y, test_size=TEST_SIZE, stratify=y, random_state=SEED
    )
    scaler_sup = StandardScaler().fit(X_tr)
    X_tr_z = scaler_sup.transform(X_tr)
    X_te_z = scaler_sup.transform(X_te)

    # ---- LDA (supervised) on original features
    priors = infer_uniform_priors(y_tr)
    lda_sup = LinearDiscriminantAnalysis(priors=priors, store_covariance=True)
    lda_sup.fit(X_tr_z, y_tr)
    y_pred = lda_sup.predict(X_te_z)
    acc = accuracy_score(y_te, y_pred)

    classes = np.unique(y)
    cm = confusion_matrix(y_te, y_pred, labels=classes)
    cm_norm = cm.astype(float) / np.maximum(cm.sum(axis=1, keepdims=True), 1.0)

    # save metrics
    (outdir / "metrics.json").write_text(
        json.dumps({
            "accuracy": float(acc),
            "classes": [int(c) for c in classes],
            "test_size": TEST_SIZE,
            "seed": SEED,
            "pca_components": PCA_NCOMP
        }, indent=2),
        encoding="utf-8"
    )
    pd.DataFrame(cm_norm, index=classes, columns=classes).to_csv(outdir / "confusion_matrix_norm.csv", index=True)
    plot_confusion_norm(cm_norm, classes, outdir / "confusion_matrix_norm.png")
    print(f"[LDA] accuracy (original space): {acc:.4f}")

    # ---- PCA on full X, then re-train LDA in PCA space
    scaler_full = StandardScaler().fit(X)
    X_z  = scaler_full.transform(X)
    pca  = PCA(n_components=PCA_NCOMP, svd_solver="full", random_state=SEED).fit(X_z)
    Xr   = pca.transform(X_z)

    lda_red = LinearDiscriminantAnalysis(priors=infer_uniform_priors(y), store_covariance=True).fit(Xr, y)

    # ---- transform new data & predict
    X_new_r = pca.transform(scaler_full.transform(X_new))
    new_pred = lda_red.predict(X_new_r)
    pd.DataFrame({"pred_label": new_pred}).to_csv(outdir / "new_predictions.csv", index=False)
    print(f"[LDA] predicted labels for new samples -> {outdir/'new_predictions.csv'}")

    # ---- Mahalanobis distances to class centroids (PCA space)
    cents = class_centroids(Xr, y)
    S_pooled = pooled_within_covariance(Xr, y)
    S_inv = np.linalg.pinv(S_pooled + 1e-8 * np.eye(S_pooled.shape[0]))

    classes_sorted = np.array(sorted(cents.keys()))
    D = np.zeros((X_new_r.shape[0], len(classes_sorted)), dtype=float)
    for i, x in enumerate(X_new_r):
        for j, c in enumerate(classes_sorted):
            mu = cents[c]
            D[i, j] = mahalanobis(x, mu, S_inv)

    pd.DataFrame(D, columns=[f"class_{int(c)}" for c in classes_sorted]) \
        .to_csv(outdir / "mahalanobis_distances.csv", index=False)
    print(f"[MD] saved Mahalanobis distances -> {outdir/'mahalanobis_distances.csv'}")

    # ---- optional: IsolationForest novelty detection per predicted class
    if RUN_IFOREST:
        preds = []
        for i, x in enumerate(X_new_r):
            c_hat = new_pred[i]
            mask = (y == c_hat)
            X_if = Xr if mask.sum() < 10 else Xr[mask]
            clf = IsolationForest(
                contamination=IF_CONTAM,
                n_estimators=IF_N_EST,
                random_state=SEED
            ).fit(X_if)
            preds.append(int(clf.predict(x.reshape(1, -1))[0]))  # 1=inlier, -1=outlier
        pd.DataFrame({"iforest_pred": preds}).to_csv(outdir / "novelty_isolation_forest.csv", index=False)
        print(f"[IForest] saved -> {outdir/'novelty_isolation_forest.csv'}")

    print(f"\nAll done. Outputs in: {outdir.resolve()}")

if __name__ == "__main__":
    main()
