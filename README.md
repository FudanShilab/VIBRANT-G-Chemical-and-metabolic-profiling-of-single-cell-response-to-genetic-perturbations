# VIBRANT-G-Chemical-and-metabolic-profiling-of-single-cell-response-to-genetic-perturbations
*Authors: Minqian Wei, Minjie Fu, Tongqi Wang, Yuchen Sun, Mingyu Wang, Dan Ye, Wei Hua, Lixue Shi**  

This package provides MATLAB and Python scripts for preprocessing, feature selection, spectral correction, and phenotype classification described in the paper.

---

## MATLAB Scripts
- **`bc_rubber.m`** – Baseline correction (rubber-band).  
- **`spcnormalize.m`** – Normalization of single-cell spectra.  
- **`Metabolic activity related feature selection.m`** – Feature selection by correlation with metabolic proxies (1616 & 2096 cm⁻¹).  
- **`spectral shift correction.m`** – Penalized Reference Matching (PRM) for segmented spectral alignment.  

## Python Scripts
- **`expression_level_prediction.py`** – Predicts relative expression level (e.g., WB results) from spectral distances.  
- **`VIBRANT_G_phenotype_classification.py`** – Phenotype classification and novelty detection: LDA + PCA, Mahalanobis distances, IsolationForest.  

---

## Data Format
- **Spectral CSV (MATLAB):**  
  Row 1 = wavenumber axis; Col 1 = cell ID; Col 2 = group (optional); Col 3..end = intensities.  
- **Python inputs:**  
  `train_data.csv` (features), `train_label.csv` (labels), `test_data.csv` (new samples).  

---

## Dependencies
- MATLAB R2023b or later.  
- Python 3.8+, with `numpy`, `pandas`, `scikit-learn`, `matplotlib`, `seaborn`.  

---

## Citation
If you use these scripts, please cite:  
**"VIBRANT-G: Chemical and metabolic profiling of single-cell response to genetic perturbations"**

## License
This repository is distributed under the **Apache License, Version 2.0**.  
You may obtain a copy of the License at [http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0).  
Unless required by applicable law or agreed to in writing, this software is provided "AS IS", without warranties or conditions of any kind.


