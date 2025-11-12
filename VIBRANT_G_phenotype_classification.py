"""
=============================================================================
Copyright (c) 2025 Lixue Shi and collaborators
Licensed under the Apache License, Version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)

Script Name: VIBRANT_G_phenotype_classification.py
Description:
This script performs phenotype classification and novelty detection of
single-cell spectra from genetic perturbation experiments.
The analysis proceeds in two major parts:
(1) **LDA-based phenotype classification.**
    Labeled spectra from reference genetic perturbations are used to train
    a Linear Discriminant Analysis (LDA) classifier. The model is evaluated
    by stratified train/test splitting, and classification accuracy together
    with the normalized confusion matrix are reported.
(2) **Two-step novelty detection.**
    First, Mahalanobis distances are calculated between spectra of the test
    condition and each reference perturbation class, and the class with the
    minimum distance is identified as the closest phenotypic match.
    Second, an IsolationForest model trained on that closest class is used
    to classify each single-cell spectrum as an inlier (1, phenotypically
    similar) or outlier (-1, phenotypically distinct). Averaging the binary
    predictions across all cells yields a mean novelty score ranging from
    1 (fully known phenotype) to -1 (fully novel phenotype). The outlier
    percentage, defined as the fraction of single cells classified as
    outliers, provides a quantitative measure of phenotypic deviation from
    the reference.

Data format:
- TRAIN_CSV: Combined training data, first column = label, remaining columns = features.
- TEST_CSV:  New data, same feature columns as training data, without label column.

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


import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

from numbers import Real
from scipy.spatial import distance
from scipy.stats import mode
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import IsolationForest
from sklearn.metrics import confusion_matrix, accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.covariance import ledoit_wolf, shrunk_covariance


def main():
    # ========= I/O =========
    # Update these paths to your own files
    input_train = "/path/to/training_data.csv"
    input_test = "/path/to/testing_data.csv"

    # Load training data (first column = label, rest = features)
    train_df = pd.read_csv(input_train)
    label_col = train_df.columns[0]
    y_raw = train_df.iloc[:, 0].astype(str).values
    data_all = train_df.iloc[:, 1:].values

    # Load testing data.
    # NOTE: This assumes the test CSV contains **features only** with the same
    # number/order as the training features. If your test file also has a first
    # label column, change the next line to `pd.read_csv(input_test).iloc[:, 1:].values`.
    newspec = pd.read_csv(input_test).values

    # ========= Label encoding =========
    le = LabelEncoder()
    y_enc = le.fit_transform(y_raw)
    classes = le.classes_
    n_classes = len(classes)
    priors = [1.0 / n_classes] * n_classes
    print(f"[INFO] Number of classes: {n_classes}, class names: {list(classes)}")

    # ========= LDA training & evaluation (original feature space) =========
    X_train, X_test, y_train, y_test = train_test_split(
        data_all, y_enc, test_size=0.3, stratify=y_enc, random_state=109
    )

    lda_clf = LinearDiscriminantAnalysis(priors=priors, store_covariance=True).fit(X_train, y_train)
    lda_pred = lda_clf.predict(X_test)
    lda_pred_accu = accuracy_score(y_test, lda_pred)
    print("Accuracy (LDA, original features):", lda_pred_accu)

    # ========= Confusion matrix (show original class names) =========
    cm = confusion_matrix(y_test, lda_pred, labels=np.arange(n_classes), normalize='true')
    df_cm = pd.DataFrame(cm, index=list(classes), columns=list(classes))
    plt.figure(figsize=(10, 6))
    plt.rcParams.update({'font.size': 8})
    plt.title('Confusion matrix of LDA classifier')
    sns.heatmap(df_cm, annot=True, fmt='.2%', cbar=True)
    plt.xlabel('Predicted label')
    plt.ylabel('True label')
    plt.tight_layout()
    plt.show()

    # ========= PCA reduction (more stable embedding & distances) =========
    n_comp = min(50, data_all.shape[1])
    pca = PCA(n_components=n_comp, svd_solver='full')
    data_mean = np.mean(data_all, axis=0)
    data_all_reduced = pca.fit_transform(data_all - data_mean)
    newspec_reduced = pca.transform(newspec - data_mean)

    # ========= LDA in PCA space (for nearest-class prediction) =========
    lda_clf_novel = LinearDiscriminantAnalysis(priors=priors, store_covariance=True).fit(
        data_all_reduced, y_enc
    )
    idx_matrix = lda_clf_novel.predict(newspec_reduced)
    idx_mode = mode(idx_matrix, keepdims=False).mode
    idx = int(np.asarray(idx_mode).item())
    closest_class_name = classes[idx]
    print("[INFO] Most likely class for test samples (encoded):", idx)
    print("[INFO] Most likely class for test samples (name):", closest_class_name)

    # ========= Covariance & Mahalanobis distance helpers =========
    def empirical_covariance(X, *, assume_centered=False):
        X = np.asarray(X)
        if X.ndim == 1:
            X = np.reshape(X, (1, -1))
        if X.shape[0] == 1:
            warnings.warn("Only one sample available. Consider reshaping your data.")

        if assume_centered:
            cov = np.dot(X.T, X) / X.shape[0]
        else:
            cov = np.cov(X.T, bias=1)
        if cov.ndim == 0:
            cov = np.array([[cov]])
        return cov

    def _cov(X, shrinkage=None, covariance_estimator=None):
        if covariance_estimator is None:
            shrinkage = "empirical" if shrinkage is None else shrinkage
            if isinstance(shrinkage, str):
                if shrinkage == "auto":
                    sc = StandardScaler()
                    Xs = sc.fit_transform(X)
                    s = ledoit_wolf(Xs)[0]
                    s = sc.scale_[:, np.newaxis] * s * sc.scale_[np.newaxis, :]
                elif shrinkage == "empirical":
                    s = empirical_covariance(X)
                else:
                    raise ValueError("Unknown shrinkage specification.")
            elif isinstance(shrinkage, Real):
                s = shrunk_covariance(empirical_covariance(X), shrinkage)
            else:
                raise ValueError("Invalid shrinkage type.")
        else:
            if shrinkage is not None and shrinkage != 0:
                raise ValueError("Cannot set both covariance_estimator and shrinkage.")
            covariance_estimator.fit(X)
            if not hasattr(covariance_estimator, "covariance_"):
                raise ValueError(
                    f"{covariance_estimator.__class__.__name__} lacks covariance_ attribute"
                )
            s = covariance_estimator.covariance_
        return s

    def _class_cov(X, y, shrinkage=None, covariance_estimator=None):
        X = np.asarray(X)
        y = np.asarray(y)
        classes_ = np.unique(y)
        cov = np.zeros((X.shape[1], X.shape[1]))
        priors_vec = np.bincount(y) / float(len(y))
        for c in classes_:
            Xc = X[y == c]
            cov += priors_vec[c] * np.atleast_2d(_cov(Xc, shrinkage, covariance_estimator))
        return cov

    def get_class_stats(X, y, pooled_cov):
        X = np.asarray(X)
        y = np.asarray(y)
        classes_ = np.unique(y)
        inv = np.linalg.pinv(pooled_cov)
        class_means = []
        covinv_list = []
        for c in classes_:
            Xc = X[y == c]
            class_means.append(np.mean(Xc, axis=0))
            covinv_list.append(inv)  # shared pooled covariance inverse (LDA-consistent)
        return np.array(class_means), covinv_list

    def mahalanobis_to_classes(class_means, covinv_list, Xq):
        Xq = np.asarray(Xq)
        K = class_means.shape[0]
        M = Xq.shape[0]
        dist_matrix = np.zeros((M, K))
        for i in range(M):
            for j in range(K):
                dist_matrix[i, j] = distance.mahalanobis(Xq[i, :], class_means[j, :], covinv_list[j])
        return dist_matrix

    # ========= Mahalanobis distances (in PCA space) =========
    pooled_cov = _class_cov(data_all_reduced, y_enc)
    class_means, covinv_list = get_class_stats(data_all_reduced, y_enc, pooled_cov)
    dist_m = mahalanobis_to_classes(class_means, covinv_list, newspec_reduced)
    pd.DataFrame(dist_m, columns=[f"class_{name}" for name in classes]).to_csv("dist_m.csv", index=False)
    print("[INFO] Mahalanobis distance matrix saved to dist_m.csv")

    # ========= IsolationForest novelty detection =========
    mask_closest = (y_enc == idx)
    if np.sum(mask_closest) < 5:
        warnings.warn(
            "Nearest-class sample count is small; IsolationForest results may be unstable."
        )

    iso_clf = IsolationForest(contamination=0.03, random_state=10, n_estimators=1000)
    iso_clf.fit(data_all[mask_closest])

    nov = iso_clf.predict(newspec)  # 1 = inlier, -1 = outlier
    novelty_score = np.mean(nov)
    outlier_pct = np.sum(nov == -1) / len(nov)
    print("IsolationForest predictions (1=inlier, -1=outlier):", nov)
    print("Novelty score (mean of predictions):", novelty_score)
    print("Outlier percentage:", outlier_pct)

    print("[DONE] Pipeline finished.")


if __name__ == "__main__":
    # Optional: choose a non-interactive backend if running headless (CI, servers)
    # matplotlib.use("Agg")
    main()
