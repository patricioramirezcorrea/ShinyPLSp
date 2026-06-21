# ShinyPLSp

### Advanced PLS-SEM Analysis in R Through an Interactive Shiny Interface

ShinyPLSp is an open-source R package that provides an interactive environment for conducting **Partial Least Squares Structural Equation Modeling (PLS-SEM)**. Designed for researchers, academics, and data analysts, the package integrates model specification, estimation, validation, segmentation, and performance analysis within a unified graphical interface.

ShinyPLSp combines the methodological rigor of the **cSEM** ecosystem with advanced segmentation capabilities based on **PATHMOX**, allowing users to perform sophisticated analyses without extensive programming experience.

***

## Key Features

### PLS-SEM Estimation

* Reflective and formative measurement models

* Structural model specification through a graphical interface

* Multiple weighting approaches:

  * PLS-PM

  * GSCA

  * SUMCORR

  * ML

  * OLS

  * 2SLS

* PLSc correction for reflective constructs

* Bootstrap and jackknife resampling procedures

### Measurement Model Assessment

* Indicator loadings

* Reliability assessment

  * Cronbach’s Alpha

  * rho\_A

* Convergent validity

  * Average Variance Extracted (AVE)

* Discriminant validity

  * HTMT

  * Fornell-Larcker criterion

### Structural Model Assessment

* Path coefficients

* Statistical significance testing

* R² values

* Model fit statistics (SRMR)

* Estimated path diagrams

### PATHMOX Segmentation

Advanced heterogeneity detection through:

* Tree-based segmentation

* MICOM-based measurement invariance assessment

* Henseler MGA comparisons

* Variable importance analysis

* Terminal-node SEM evaluation

### Project Management

* Save complete projects as `.rds`

* Restore previous sessions instantly

* Export detailed Excel reports

* Export reproducible R scripts

***

## Installation

Install the development version directly from GitHub:

```r
# install.packages("remotes")
remotes::install_github("patricioramirezcorrea/ShinyPLSp")
```

Launch the application:

```r
library(ShinyPLSp)
run_app()
```

***

## Data Requirements

Input datasets must satisfy the following conditions:

* Excel format (`.xls` or `.xlsx`)

* First row contains variable names

* Subsequent rows contain observations

* SEM indicators must be numeric

* PATHMOX segmentation variables may be numeric or categorical

***

## Recommended Workflow

1. Import data

2. Handle missing values

3. Define constructs

4. Specify structural relations

5. Run PLS-SEM estimation

6. Evaluate measurement quality

7. Assess structural relationships

8. Perform PATHMOX segmentation (optional)

***

## Citation

If you use ShinyPLSp in academic research, please cite:

Ramírez-Correa, P. (2026). _ShinyPLSp: An Interactive Web Application for Advanced PLS-SEM Analysis_ (Version 1.0) \[Computer software].

### Underlying Methodological Engines

Please also cite:

* Rademaker, M., & Schuberth, F. (2020). _cSEM_.

* Lamberti, L., Noci, G., Guo, J., & Zhu, S. (2016). _genpathmox_.

***

## License

ShinyPLSp is released under the **GNU General Public License v3 (GPL-3)**.

You are free to use, modify, and distribute the software provided that derivative works remain under the same license.

***

## Disclaimer

The software is provided _as is_, without warranty of any kind.

Researchers should not rely solely on automated statistical outputs. Results must always be interpreted within appropriate theoretical, methodological, and substantive frameworks. ShinyPLSp facilitates statistical computation; it does not replace scientific judgment.

***

## Author

Patricio Ramírez-Correa

Professor and Researcher

Universidad Católica del Norte

Chile

