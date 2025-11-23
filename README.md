# Mass Spectrometry Fragmentation Prediction: Ensemble Learning Approach

## Overview

This repository contains supplementary materials for the scientific research **"Improving In-Silico Fragmentation Prediction via Ensemble Machine Learning"** authored by Gregory Kleiner, Alexander Tarnopolsky, and Michael Elgart at the Molecular Discovery Lab, Sheba Medical Center.

The work demonstrates that combining multiple fragmentation prediction tools through ensemble machine learning methods yields substantial performance improvements over individual predictors when analyzing novel chemical compounds.

## Research Summary

**Problem Statement:** Current in-silico mass spectrometry fragmentation tools achieve less than 50% accuracy on novel compounds, with inter-tool agreement below 20%, limiting practical applicability for metabolite identification.

**Solution:** An ensemble learning framework that aggregates predictions from five independent tools (ICEBERG, SCARF, RASSP, CFM-ID, MassFormer) combined with molecular descriptors to predict peak presence/absence across mass spectrum bins.

**Results:**
- Training dataset: 21,023 novel compounds from GNPS repository (submissions after January 2024)
- Optimal model: XGBoost classifier achieving 0.79 accuracy and 0.75 F1 score
- Performance gain: Approximately 40% improvement over best individual predictor
- Feature importance: SHAP analysis identifies precursor mass, structural topology, hydrophobicity, and tool-specific predictions as critical

## Citation

**BibTeX Format:**
```bibtex
@article{kleiner2025msensemble,
  author = {Kleiner, Gregory and Tarnopolsky, Alexander and Elgart, Michael},
  title = {Improving In-Silico Fragmentation Prediction via Ensemble Machine Learning},
  institution = {Molecular Discovery Lab, Metabolic Center, Sheba Medical Center},
  address = {Tel HaShomer, Israel},
  year = {2025},
  howpublished = {\url{https://github.com/MDU-AI-MS-Project/MDU-AI-MS-Project}}
}
```

**Text Format:**
```
Kleiner, G., Tarnopolsky, A., & Elgart, M. (2025).
Improving In-Silico Fragmentation Prediction via Ensemble Machine Learning.
Molecular Discovery Lab, Metabolic Center, Sheba Medical Center, Tel HaShomer, Israel.
Repository: https://github.com/MDU-AI-MS-Project/MDU-AI-MS-Project
```

## License

This work is released under the **MIT License**.

**Copyright Notice:**  
© 2025 Molecular Discovery Lab, Sheba Medical Center

**Terms:**

Permission is granted, without charge, to any individual obtaining a copy of this software and associated documentation (the "Materials"), to utilize the Materials without restriction, including rights to use, copy, modify, merge, publish, distribute, sublicense, and sell copies, subject to the following conditions:

- This copyright notice and permission notice must be included in all copies or substantial portions of the Materials.

**Warranty Disclaimer:**  
THE MATERIALS ARE PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES, OR OTHER LIABILITY ARISING FROM THE MATERIALS OR THEIR USE.

**Third-Party Components:**  
Portions of utility code may derive from publicly available examples in scikit-learn, RDKit, and XGBoost documentation. Such components retain their original licenses. Integration code and ensemble architectures are original works of the authors.

Full license text: https://opensource.org/licenses/MIT

## Disclaimer

**Research Software Notice**

This repository contains experimental software developed for scientific research purposes. Users should be aware of the following:

**Intended Use:**
- Academic research and education
- Method benchmarking and comparison
- Non-clinical computational chemistry applications

**Limitations and Warnings:**
1. **No Warranties:** Software provided without guarantees of correctness, completeness, or suitability for specific purposes
2. **Generalization:** Models trained on specific GNPS dataset; performance may differ on other chemical spaces
3. **External Dependencies:** Predictions require third-party tools (CFM-ID, RASSP, etc.) with separate licensing
4. **Not for Clinical Use:** Results require experimental validation before any diagnostic or therapeutic applications
5. **Liability:** Authors and Sheba Medical Center assume no responsibility for damages resulting from software use
6. **Data Privacy:** Ensure compliance with institutional policies when processing proprietary chemical data

**Recommendations:**
- Validate predictions with experimental MS/MS data when possible
- Consult domain experts for interpretation in critical applications
- Review external tool licenses before commercial deployment
- Perform independent testing on representative data before production use

**Regulatory Compliance:**  
This software has not been evaluated or approved by regulatory agencies (FDA, EMA, etc.). Any clinical, diagnostic, or pharmaceutical applications require appropriate regulatory review and approval processes.

## Contact Information

**Corresponding Author:**  
Michael Elgart, PhD  
Email: michael.elgart@sheba.health.gov.il

**Institution:**  
Molecular Discovery Lab  
Metabolic Center  
Sheba Medical Center  
Tel HaShomer 5265601, Israel

**Contributing Authors:**  
- Gregory Kleiner (Co-first author)
- Alexander Tarnopolsky (Co-first author)

**Repository Issues:**  
For technical questions or bug reports, please open an issue on the GitHub repository.

## Acknowledgments

This research builds upon computational predictions from CFM-ID, RASSP, MassFormer, Iceberg, and SCARF. We gratefully acknowledge the developers of these tools and the GNPS community for maintaining publicly accessible spectral databases.

Data analysis performed using the Python scientific computing ecosystem: NumPy, pandas, scikit-learn, XGBoost, RDKit, matplotlib, and SHAP.

---

# Repository Structure and Technical Documentation

## File Tree Structure

```
MDU-AI-MS-Project/
├── README.md
├── base_methods/
│   ├── cfm-id/
│   │   └── run_cfmid_batch.sh
│   ├── iceberg/
│   │   └── run_model_chunk.sh
│   ├── massformer/
│   │   └── run_model_chunk.sh
│   ├── rassp/
│   │   └── run_model_chunk.sh
│   └── scarf/
│       └── run_model_chunk.sh
├── data/
│   └── GNPSnew_filtered.txt
├── ensemble/
│   ├── SPLIT_DB.sh
│   ├── mslearn_seed_array.sh
│   ├── mslearn_distrib_hypsrch_v2.py
│   ├── mslearn_distrib_valid.py
│   ├── split_orig_data_distrib_v1.py
│   ├── mslearn_input_extended_v2_cleaned_valid.json
│   └── mslearn_input_split[1-9].json
├── feature_engineering/
│   ├── msp_utils.py
│   ├── mdu_batch.py
│   └── table_aggregator.py
├── graphics/
│   ├── binary_plot.py
│   ├── ms_correlation_plot.py
│   ├── ms_correlation_plot.R
│   ├── ms_graphviz.py
│   ├── ms_shap.py
│   └── ms_svg_merge.py
└── models/
    ├── 2107-XGBoost_Paramset449_valid-AllDescAllPred.json
    ├── 2107-XGBoost_Paramset449_valid-AllDescAllPred.pkl
    ├── 2108-XGBoost_Paramset1025_valid-AllDescAllPred.json
    ├── 2108-XGBoost_Paramset1025_valid-AllDescAllPred.pkl
    └── readme
```

## File Descriptions

### Root Directory

- **README.md**: Main documentation file for the repository

### base_methods/

Contains shell scripts for running third-party mass spectrometry prediction tools:

- **base_methods/cfm-id/run_cfmid_batch.sh**: SLURM batch script for CFM-ID predictions using Singularity container
- **base_methods/iceberg/run_model_chunk.sh**: SLURM array job script for ICEBERG predictions (600 chunks)
- **base_methods/massformer/run_model_chunk.sh**: SLURM array job script for MassFormer predictions (600 chunks)
- **base_methods/rassp/run_model_chunk.sh**: SLURM array job script for RASSP predictions (600 chunks)
- **base_methods/scarf/run_model_chunk.sh**: SLURM array job script for SCARF predictions (600 chunks)

### data/

- **GNPSnew_filtered.txt**: Input dataset containing SMILES strings for novel compounds from GNPS (21,147 entries)

### ensemble/

Machine learning ensemble training and validation scripts:

- **SPLIT_DB.sh**: SLURM script for generating train/validation/test splits (splits 2-10)
- **mslearn_seed_array.sh**: SLURM array job for hyperparameter search across multiple configurations
- **mslearn_distrib_hypsrch_v2.py**: Hyperparameter search script for ensemble classifiers
- **mslearn_distrib_valid.py**: Main validation script for training and evaluating ML classifiers
- **split_orig_data_distrib_v1.py**: Data splitting utility (70% train, 15% validation, 15% test)
- **mslearn_input_extended_v2_cleaned_valid.json**: Configuration file for model validation experiments
- **mslearn_input_split[1-9].json**: Individual split configuration files

### feature_engineering/

Data processing and feature extraction modules:

- **msp_utils.py**: Mass spectrum parsing utilities for multiple file formats (MSP, CFM-ID, RASSP, MassFormer, ICEBERG, SCARF, reference JSON)
- **mdu_batch.py**: Batch processing pipeline for aggregating predictions and computing molecular descriptors
- **table_aggregator.py**: Data aggregation class for managing molecular data, canonical SMILES, and synonym handling

### graphics/

Visualization and plotting scripts:

- **binary_plot.py**: Plotly-based binary classification performance visualization
- **ms_correlation_plot.py**: Python/seaborn correlation heatmap generator
- **ms_correlation_plot.R**: R/corrplot correlation heatmap with custom color compression
- **ms_graphviz.py**: GraphViz decision tree visualization
- **ms_shap.py**: SHAP value calculation and visualization for model interpretability
- **ms_svg_merge.py**: SVG figure merging utility for publication-ready graphics

### models/

Pre-trained ensemble model artifacts:

- **2107-XGBoost_Paramset449_valid-AllDescAllPred.json**: Model configuration metadata
- **2107-XGBoost_Paramset449_valid-AllDescAllPred.pkl**: Serialized XGBoost model (primary recommended model)
- **2108-XGBoost_Paramset1025_valid-AllDescAllPred.json**: Alternative model configuration metadata
- **2108-XGBoost_Paramset1025_valid-AllDescAllPred.pkl**: Serialized XGBoost model (alternative parameters)
- **readme**: Model inventory and performance summary

## Third-Party Dependencies

### Python Packages

```
pandas >= 1.5.0
numpy >= 1.23.0
scikit-learn >= 1.2.0
xgboost >= 1.7.0
rdkit >= 2022.09.1
shap >= 0.41.0
matplotlib >= 3.6.0
seaborn >= 0.12.0
plotly >= 5.11.0
graphviz >= 0.20.0
```

### R Packages

```
corrplot >= 0.92
readr >= 2.1.0
```

### System Requirements

- **Singularity**: For CFM-ID containerized execution
- **SLURM**: Workload manager for HPC cluster execution
- **Conda/Miniconda**: Environment management

### External Tools

- **CFM-ID 4.0**: Mass spectrum prediction (containerized)
- **ICEBERG**: Fragment prediction model
- **MassFormer**: Transformer-based MS prediction
- **RASSP**: Neural network-based spectrum prediction
- **SCARF**: Spectral prediction framework

## Installation Instructions

### 1. Clone Repository

```bash
git clone https://github.com/MDU-AI-MS-Project/MDU-AI-MS-Project.git
cd MDU-AI-MS-Project
```

### 2. Create Conda Environment

```bash
conda create -n msensemble python=3.9
conda activate msensemble
```

### 3. Install Python Dependencies

```bash
pip install pandas numpy scikit-learn xgboost
pip install rdkit shap matplotlib seaborn plotly graphviz
```

### 4. Install R and Packages (for correlation plots)

```bash
conda install -c conda-forge r-base r-corrplot r-readr
```

### 5. Configure External Tools

**Note**: External prediction tools (CFM-ID, ICEBERG, MassFormer, RASSP, SCARF) must be installed separately according to their respective documentation.

## Runtime Instructions

### Important: Update Hardcoded Paths

**Before running any scripts**, update the hardcoded file paths to match your working environment:

- `ensemble/split_orig_data_distrib_v1.py`: Update paths in lines 40, 52, 65, 79-98
- `ensemble/mslearn_distrib_valid.py`: Update workdir paths in line 145
- `feature_engineering/mdu_batch.py`: Update rootdir path in line 16
- `graphics/*.py`: Update root variable in each script
- `base_methods/*/run_model_chunk.sh`: Update input/output paths in each SLURM script

### Step 1: Generate Predictions from Base Methods

Run each base method on your compound dataset:

```bash
# Update paths in scripts before execution
sbatch base_methods/cfm-id/run_cfmid_batch.sh
sbatch base_methods/iceberg/run_model_chunk.sh
sbatch base_methods/massformer/run_model_chunk.sh
sbatch base_methods/rassp/run_model_chunk.sh
sbatch base_methods/scarf/run_model_chunk.sh
```

### Step 2: Process Predictions and Extract Features

```bash
python feature_engineering/mdu_batch.py
```

This script:
- Aggregates predictions from all base methods
- Computes molecular descriptors using RDKit
- Processes reference spectra
- Bins spectral data into 2800 bins (0.5 Da resolution, 100-1500 Da range)

### Step 3: Generate Train/Validation/Test Splits

```bash
sbatch ensemble/SPLIT_DB.sh
```

Creates 9 stratified splits with 70% train, 15% validation, 15% test ratios.

### Step 4: Train Ensemble Models

```bash
# Single configuration
python ensemble/mslearn_distrib_valid.py -w /path/to/workdir -i 1 -j mslearn_input_extended_v2_cleaned_valid.json

# Hyperparameter search across all configurations
sbatch ensemble/mslearn_seed_array.sh
```

Optional: Add `--save-model` flag to save trained models to `models` directory.

### Step 5: Generate Visualizations

```bash
# Binary classification performance
python graphics/binary_plot.py

# Correlation analysis
python graphics/ms_correlation_plot.py
Rscript graphics/ms_correlation_plot.R

# SHAP interpretability
python graphics/ms_shap.py

# Merge figures for publication
python graphics/ms_svg_merge.py
```

## Configuration Files

- **mslearn_input_extended_v2_cleaned_valid.json**: Contains experiment configurations with parameters:
  - `id`: Unique experiment identifier
  - `train_path`: Path to training CSV
  - `test_path`: Path to test CSV
  - `classification_mode`: "binary" or "binned"
  - `fitting_set`: Feature set name (AllPred, AllDesc, AllDescAllPred, etc.)
  - `classifier`: ML algorithm name (XGBoost, Random-Forest, Neural-Net, etc.)

## Output Files

- **Training results**: CSV files with accuracy, precision, recall, F1 scores
- **Trained models**: Serialized model files in `models` (if `--save-model` flag used)
- **SHAP values**: CSV files with feature importance rankings
- **Visualizations**: SVG/PNG plots in `graphics` directory

---

**Repository URL:** https://github.com/MDU-AI-MS-Project/MDU-AI-MS-Project  
**Last Updated:** November 2025