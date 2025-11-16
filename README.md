# Mass Spectrometry Fragmentation Prediction: Ensemble Learning Approach

## Overview

This repository contains supplementary materials for the scientific research **"Improving In-Silico Fragmentation Prediction via Ensemble Machine Learning"** authored by Gregory Kleiner, Alexander Tarnopolsky, and Michael Elgart at the Molecular Discovery Lab, Sheba Medical Center.

The work demonstrates that combining multiple fragmentation prediction tools through ensemble machine learning methods yields substantial performance improvements over individual predictors when analyzing novel chemical compounds.

## Research Summary

**Problem Statement:** Current in-silico mass spectrometry fragmentation tools achieve less than 50% accuracy on novel compounds, with inter-tool agreement below 20%, limiting practical applicability for metabolite identification.

**Solution:** An ensemble learning framework that aggregates predictions from five independent tools (ICEBERG, SCARF, RASSP, CFM-ID, MassFormer) combined with molecular descriptors to predict fragment ion presence/absence across mass spectrum bins.

**Results:**
- Training dataset: 21,147 novel compounds from GNPS repository (submissions after January 2024)
- Optimal model: XGBoost classifier achieving 0.70 accuracy and 0.57 F1 score
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

**Repository URL:** https://github.com/MDU-AI-MS-Project/MDU-AI-MS-Project  
**Last Updated:** November 2025  
**Version:** 1.0.0
