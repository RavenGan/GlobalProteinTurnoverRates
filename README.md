# GlobalProteinTurnoverRates
This repository contains the code used for the calculation of the global protein turnover rates. Explanations for the code are as follows.

- `./main/preprocess.R`: This script preprocesses the data  from the original source. It creates the metadata for samples, filters out missing entries, and did log-transformation.
- `./main/main_analysis.R`: This script contains the main analysis for the calculation of the global protein turnover rates.
- `./fcn/fcn_metric.R`: This script contains the functions used in the main analysis.