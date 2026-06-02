⚖️ Legal Notice and Software License
Open Source Software Agreement
Unlike proprietary commercial software, ShinyPLSp is an open-source project developed for the academic and scientific community. The source code of the software is released under the GNU General Public License (GPL) version 3. This means you are free to use, modify, and distribute the software, provided that any derivative works are also distributed under the same open-source license.
Disclaimer of Warranty and Liability
The Software is provided “as is”, without any warranty of any kind. Under no circumstances shall the authors be held liable for any direct, indirect, incidental, or consequential damages caused by the use of or inability to use the Software.
Important Warning on Statistical Analysis
Users are strongly cautioned against blindly accepting the results provided by ShinyPLSp. You must actively double-check and validate your results against applicable theoretical frameworks. The software computes algorithms based on the data you provide; it does not replace critical thinking.
📌 How to Cite ShinyPLSp
If you use ShinyPLSp to analyze data for a publication, you must properly cite the software and its underlying statistical engines:
Ramírez-Correa, P. (2026). ShinyPLSp: An Interactive Web Application for Advanced PLS-SEM Analysis (Version 1.0) [Computer software].
Engines citations: Users must also cite the cSEM package (Rademaker & Schuberth, 2020) for the PLS-SEM estimations, and the genpathmox package (Lamberti et al., 2016) if the PATHMOX segmentation module is used.
1. General Presentation of the Application
ShinyPLSp is a modular analytical platform designed for researchers and professionals who need to perform Partial Least Squares Structural Equation Modeling (PLS-SEM). The application guides the user through a sequential workflow: loading data, preparing the database, defining the measurement model, defining the structural model, estimating, interpreting results, and performing advanced segmentation and performance analysis.
Purpose: To estimate complex relationships between latent variables (constructs) and observable variables (indicators).
Target User Profile: Researchers and data analysts needing a robust, code-free interface for PLS-SEM.
Main Libraries Used: cSEM (unified estimation engine) and genpathmox (detecting unobserved heterogeneity).
2. Requirements and Preparation
To ensure successful estimation, users must prepare their database following these guidelines:
File Format: The application exclusively accepts Excel files (.xls, .xlsx).
Structure: The first row must contain the names of the variables (indicators). The subsequent rows must be the cases or observations.
Data Types: Indicators intended for the SEM model must be numeric. Segmentation variables for the PATHMOX module can be either categorical (text) or numeric.
3. General View of the Interface
The application is organized using a sidebar layout alongside main top navigation tabs.
Sidebar: Contains critical controls that affect the entire application: project saving/loading, data import, missing data treatment, estimation algorithms, and resampling (Bootstrap) parameters.
Main Panels (Tabs): Organized sequentially from left to right.
Recommendation: Follow this flow strictly from left to right, as later estimations depend on the existence of valid data and constructs defined in earlier tabs.
4. Step-by-Step Manual by Modules
4.1. Project Configuration and Data Import
Objective: Load the database and handle missing values before estimation.
Code to be treated as missing: Enter the numeric value that represents missing data in your survey (e.g., -1 or 999).
Treatment of missing values:
Listwise deletion: Removes any row containing at least one NA (Recommended).
Mean/Median imputation: Replaces the NA with the average or median.
4.2. Definition of Constructs (Measurement Model)
Objective: Build the measurement model linking data columns to theoretical concepts.
Measurement Model Type:
Common factor: The construct causes the indicators (Reflective model).
Composite: The indicators form or create the construct (Formative model).
Restrictions: An indicator assigned to one construct disappears from the "Available" list, preventing cross-loadings.
4.3. Structural Relations (Structural Model)
Objective: Specify the paths (causal links) between the previously created constructs.
Validations & Restrictions: The app automatically blocks a variable from predicting itself. Reciprocal relationships (A -> B and B -> A) and structural cycles (A -> B -> C -> A) are strictly prohibited. The UI will show an error notification ("Relation blocked").
4.4. SEM Estimation and Results
Objective: Execute the statistical estimation. Press the green "▶ Run PLS-SEM" button in the sidebar.
Outputs: Includes tables for Model fit, Construct types, R², Paths, Loadings, Measurement quality, HTMT, and Fornell-Larcker. It also renders the estimated path diagram.
Export: The "Export detailed Excel" button generates a multi-sheet .xlsx file containing the raw data, cleaned data, and all statistical tables.
4.5. WPI Analysis (Weighted Performance Index)
Objective: Evaluates individual subject performance weighted by the PLS-SEM total effects.
Internal Logic: Extracts the total effects associated with a Target Construct, rescales construct scores to a 0-100 range, and computes a weighted performance index.
Outputs: Applies an automated K-Means clustering algorithm to group users into segments and runs a t-SNE algorithm to project multi-dimensional performance into a 2D scatter plot.
4.6. PATHMOX Analysis (Segmentation)
Objective: Discover if the structural model changes significantly depending on certain variables (e.g., Gender, Income).
Advanced Feature (MICOM): Integrates a Measurement Invariance check. If the new child nodes do not interpret the survey questions identically, the split is rejected.
Outputs: Renders the PATHMOX Tree, Terminal path comparisons, Variable importance, MICOM step 2 results, and Henseler MGA p-values.
4.7. Export Reproducible R Script 🆕
Objective: Ensure 100% transparency and reproducibility of your research.
Functionality: This module translates all your visual clicks, data cleaning choices, construct definitions, structural relations, and algorithm settings into a clean, native R script.
Export: You can preview the code live and download the .R file. You can attach this script as supplementary material in your academic papers so other researchers can replicate your findings exactly.
4.8. Save/Load Projects
Save project: Downloads an .rds file containing your raw data, constructs, relations, results, and all sidebar settings.
Load project: Upload a previously saved .rds file to instantly restore your entire session exactly as you left it.
5. Detailed Explanation of Estimation Parameters
Weighting approach: Options include PLS-PM (default), GSCA, SUMCORR, ML, OLS, and 2SLS.
Path estimation: Controls how the structural paths are estimated (OLS or 2SLS).
PLS inner weighting: Path (default) is recommended for structural models.
Default outer mode: Auto is highly recommended as it reads the construct type.
PLSc correction: Activates the Dijkstra-Henseler correction to fix attenuation bias in reflective models.
Resample method & Number: Bootstrap (with replacement), jackknife, or none. 5000+ resamples are standard for final publication.
Handle inadmissibles: Replace (default) repeats the resampling until a valid solution is found.
6. Results Interpretation Guide
Model fit (SRMR): A value < 0.08 generally indicates a good overall fit.
R² (R-squared): > 0.75 is substantial, ~0.50 is moderate, and < 0.25 is weak.
Paths: A P-value < 0.05 means the hypothesis is statistically significant.
Loadings: Ideally should be > 0.708.
Measurement quality: Cronbach's Alpha and rho_A should be > 0.70. AVE must be > 0.50.
HTMT (Discriminant Validity): Values should strictly be < 0.90 (ideally < 0.85).
7. Troubleshooting and Common Errors
"Define constructs first": You attempted to run the model without defining the measurement model.
"Relation blocked (creates cycle)": You tried to draw a path that creates a feedback loop (e.g., A -> B -> A). PLS-SEM requires recursive models.
App crashes during PATHMOX: Often caused by a segmentation variable with zero variance (e.g., 100% of respondents are "Female") or terminal nodes that are too small.
"Processed data still contains NAs": If you chose "None" for missing data treatment, the cSEM algorithm will fail. Choose "Listwise deletion" or "Imputation".
8. Best Practices Recommendations
Recommended Workflow Sequence: Always follow the tabs from left to right. Changing the dataset erases constructs to prevent data mapping errors.
Naming Conventions: Avoid special characters or spaces in construct names (use Brand_Trust instead of Brand Trust!).
PATHMOX Sample Size: It is recommended to have a total sample of at least 150-200 observations before running PATHMOX.
Prototyping vs. Final Run: Set the "Resample method" to 100 resamples for exploratory learning. Switch to 5000+ for final publications.
