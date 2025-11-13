**VIBRANT-G: Chemical and metabolic profiling of single-cell response to genetic perturbations**  
*Authors: Minqian Wei, Minjie Fu, Tongqi Wang, Yuchen Sun, Mingyu Wang, Xingyuan Bi, Wei Hua, Lixue Shi** 

This package provides MATLAB and Python scripts for data preprocessing, feature extraction, visualization, and quantitative prediction of gene expression and phenotype classification, as described in the paper.

---

## MATLAB Scripts
- **`Datapreprocessing_BaselineCorrection.m`** – Performs baseline subtraction across spectra.
- **`Datapreprocessing_spectral_shift_correction.m`** – Aligns spectra to correct wavelength/wavenumber shifts.
- **`Datapreprocessing_Metabolic_activity_related_feature_selection.m`** – Selects spectral features correlated with metabolic activity (13C amide I and azido bands).
- **`spcnormalize.m`** – Normalization of single-cell spectra.  
- **`bc_rubber.m`** – Baseline correction using a rubberband algorithm.
- **`drawn3D_EE.m`** – Draws 3D error ellipsoids for confidence visualization.
- **`SingleCell_Metabolic3DScatter.m`** – Generates 3D scatter plots of metabolic activity or feature embeddings. 

## Python Scripts
- **`VIBRANT_G_Expression_level_prediction.py`** – Predicts gene expression levels based on extracted single-cell vibrational features. Uses reference data (known expression) to predict target data (to be predicted expression).  
- **`VIBRANT_G_phenotype_classification.py`** – Classifies metabolic phenotypes from vibrational data using LDA and Mahalanobis distance.
---

## **Data Format**

---

### **1. Datapreprocessing_BaselineCorrection.m**

**Input format:**  
- `.mat` file containing:  
  - `C` — FTIR hyperspectral cube `[Nx, Ny, N_wavenumbers, N_channels]`  
  - `Minfo` — metadata structure (fields: `File.LWN`, `File.UWN`, `File.NOD`, `File.WVS`)  
 
### **2. Datapreprocessing_Metabolic_activity_related_feature_selection.m**

**Input format:**  
- `.csv` file (`processed_data.csv`)  
  - **Row 1:** wavenumber axis (cm⁻¹)  
  - **Column 1:** cell ID  
  - **Column 2:** cell label  
  - **Columns 3..end:** single-cell spectra  

**Output:**  
- `top_features.csv` — top-ranked spectral features correlated with metabolic activity

### **3. Datapreprocessing_spectal_shift_correction.m**

**Input format:**  
- `data.csv` files for both **reference** and **target** datasets  
  - **Row:** single cell (one spectrum per row)  
  - **Column 1:** cell ID  
  - **Columns 2..end:** intensity values at wavenumbers specified by header row  


### **4. VIBRANT_G_Expression_level_prediction.py**

**Input format:**  
- `data_expression_level_prediction_demo_data.csv` — combined dataset containing:  
  - **WT data** (wild-type baseline)  
  - **Reference data** with known expression levels  
  - **Target data** with unknown expression levels to be predicted  

**Output:**  
- predicted expression levels for target data


### **5. VIBRANT_G_phenotype_classification.py**

**Input format:**  
- `phenotype_classification_training_data.csv` — dataset with **known genetic perturbations**  
- `phenotype_classification_testing_data.csv` — dataset with **unknown genetic perturbations**

**Output:**  
- predicted phenotype categories  

---
## **Installation**

### **Software Dependencies**
- **MATLAB:** R2023b or later  
- **Python:** 3.12 or later  
  - Required packages: `numpy`, `pandas`, `matplotlib`, `scikit-learn`, `scipy`, `umap-learn`

### **Hardware Requirements**
- Standard desktop or laptop computer  
- Minimum 8 GB RAM (recommended **16 GB or more** for large datasets)  
- Sufficient storage for hyperspectral data processing  

### **Versioning**
- Scripts have been tested on **MATLAB R2023b** and **Python 3.12**  
- For full compatibility, please use the same software versions  

---

## **Installation Guide**

### **MATLAB Setup**
1. Ensure MATLAB **R2023b or later** is installed.  
2. Place all `.m` scripts in your MATLAB working directory.  
3. Open MATLAB and execute scripts according to the workflow described in the manuscript.  
4. No additional MATLAB toolbox installation is required beyond the standard **Signal Processing Toolbox** and **Statistics and Machine Learning Toolbox**.  
5. Typical setup time: *approximately 5 minutes.*

### **Python Setup**
1. Install **Python 3.12** (or later).  
2. Install all required packages using:  `numpy`, `pandas`, `matplotlib`, `scikit-learn`, `scipy`, `umap-learn`. 
3. Place all .py scripts in the same project directory. 
4. Run analysis scripts from the command line or any Python IDE (e.g., VS Code, PyCharm, Jupyter Notebook).  
5. Typical setup time: *approximately 5 minutes.*
---

## **Citation**

If you use these scripts, please cite:  
> *VIBRANT-G: Chemical and metabolic profiling of single-cell response to genetic perturbations.*  

---

## **License**

This repository is distributed under the **Apache License, Version 2.0**.  
You may obtain a copy of the License at:  
[http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, this software is provided **"AS IS"**,  
without warranties or conditions of any kind, either express or implied.  
See the License for the specific language governing permissions and limitations under the License.

---


