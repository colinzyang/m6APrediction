# m6APrediction: A ML Tool for Predicting m6A Modification Sites in RNA

##  Overview

**m6APrediction** is an R package designed for predicting **N6-methyladenosine (m6A)** modification sites in RNA sequences.

m6A is one of the most prevalent internal RNA modifications, playing key roles in gene expression regulation, RNA stability, and disease mechanisms such as cancer. This package leverages a **Random Forest** machine learning model trained on comprehensive sequence and genomic features to identify potential methylation sites with high precision.

##  Installation

You can install the development version of m6APrediction from [GitHub](https://github.com/) with:

```r
# install.packages("devtools")
devtools::install_github("https://github.com/colinzyang/m6APrediction")

Dependencies: randomForest, stats
R Version: R ≥ 3.5.0
```
##  Model Performance
The m6APrediction model has been validated using ROC and PRC analyses, demonstrating strong discriminative power in handling imbalanced genomic data.

AUC (ROC): ~0.93
AUPRC: ~0.87

(Note: Detailed performance plots are available in the project report/vignettes)

##  License
This project is licensed under the GPL-3 License.
