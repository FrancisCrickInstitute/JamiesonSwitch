This repository contains analysis code for the project '_Overnight circuit remodelling drives juvenile alloparental care_' by Bradley B. Jamieson, Maxwell X. Chen, Swang Liang, Lina S. H. El Rasheed, Patty Wai, Grace M. K. Chattey, Maria Voronkov, Emma Davis, J. Mark Skehel, I. Lorena Arancibia-Cárcamo, Molly Strom and Johannes Kohl

**ANOVA_vsAll**
Loads a CSV of group columns (e.g., P14, P15), reshapes to long format, and checks normality across groups (Shapiro–Wilk, ≥3 values). For every pair of groups it runs either parametric ANOVA + pairwise t-tests (Holm-corrected) or Kruskal–Wallis + Dunn post-hoc (Holm-corrected), then exports a *_results.csv with tests, p-values, and significance stars.

**ANOVA_vsControl**
Compares a user-specified control column against every other group column after reshaping to long format and assessing normality. Uses OLS ANOVA with a control-vs-group multiple-comparisons adjustment when parametric, or Kruskal–Wallis + Dunn (Holm) when nonparametric, and saves a *_results.csv summary with corrected p-values and stars.

**BrainGlobe_Processing**
Reads a BrainGlobe/Cellfinder summary table of per-structure cell counts and volumes and removes a predefined list of excluded structures (e.g., fibre tracts/ventricles/broad parents). Computes cell density (cells per mm³), aggregates by structure, and exports a cleaned summary to Excel.

**Chi2_Independence**
Loads a contingency table CSV, coerces to numeric, drops missing rows, and runs a chi-square test of independence between two columns (groups). Reports raw and BH-FDR–adjusted p-values (plus stars) and saves results to *_chi2_results.csv.

**Discrimination**
Takes a stimulus-by-response CSV (columns = stimuli), computes one-vs-rest ROC AUC per stimulus with bootstrap mean ± SEM, and builds a pairwise dominance matrix P(AUC_A > AUC_B). It also computes a bootstrapped tuning specificity index per stimulus, tests TSI > 0 with Wilcoxon, applies BH-FDR correction, and exports all summaries to a single results CSV.

**Elbow_points**
Loads and concatenates morphometric CSVs, keeps user-selected features, standardizes them, and fits K-means across K=1–10. Uses the curvature (“elbow”) of inertia to pick an optimal K, plots the elbow curve, and prints the chosen K.

**Fisher's_multiple**
Reads a 2×N In/Out count table and runs Fisher’s exact test on adjacent column pairs (1 vs 2, 3 vs 4, etc.). Applies BH-FDR across all tests, assigns significance stars, and saves a results CSV with raw and corrected p-values.

**Fisher's_vsControl**
Loads an In/Out count table and compares each group column to a specified control column using Fisher’s exact test. Outputs raw p-values, applies BH-FDR correction (as implemented in the script), assigns stars, and saves a *_results.csv summary.

**Gene_ANOVA**
Parses a gene table where columns encode Region and Age (Region_Age), reshapes to long format, and fits a two-way ANOVA (Region * Age). It then performs within-region pairwise age comparisons, applies BH-FDR across all post-hoc tests, and exports ANOVA + post-hoc results to CSV.

**KaplanMeier_multiple**
Loads retrieval-time columns, reshapes to long format, and defines an event as retrieval < 600 (≥600 censored). Runs log-rank tests for adjacent column pairs, applies BH-FDR across comparisons, and saves a *_results.csv with corrected p-values and stars.

**KaplanMeier_vsControl**
Reshapes retrieval-time data and treats Vir. as the control group, with events defined as retrieval < 600. Runs log-rank tests comparing Vir. to each timepoint, BH-FDR–corrects across tests, and exports a *_results.csv summary.

**Linear_regression**
Fits an OLS linear regression between two numeric columns from a CSV and reports slope, intercept, R², and slope p-value. Generates a scatter plot with fitted line and an approximate 90% confidence band, and saves a small results CSV.

**MassSpec_CellChatDB (R)**
Filters a hypothalamus DE protein table (FDR < 0.05, |log2FC| > 1) and uses it to subset CellChat’s mouse ligand–receptor interaction database. It also plots selected genes across region × age groups using expression matrices plus metadata, saving gene-wise PDFs.

**MassSpec_PCA_DEs (R)**
Runs PCA on hypothalamus samples from normalized proteomics (after collapsing duplicate protein names) and saves a PC1–PC2 plot. Then uses limma to compute region- and region-specific P15 vs P14 differential expression contrasts, exports DE tables, and generates volcano plots (including WGCNA module annotations).

**MassSpec_Run_WGCNA (R)**
Builds a signed WGCNA network on filtered proteomics data (bicor; chosen soft-threshold), writes module assignments, and correlates module eigengenes with region×age traits with BH correction. Produces a module–trait heatmap and GO barplot panels per module, exporting PDFs and summary tables.

**MassSpec_Run_fGSEA (R)**
Runs multilevel FGSEA on a ranked hypothalamus-specific DE list using GO gene sets and exports the enrichment table. Produces a filtered, ordered NES barplot (direction-coded by NES sign) and saves it as a PDF.

**Multiple_t-tests**
Loads a numeric CSV and performs comparisons in column pairs (1 vs 2, 3 vs 4, etc.), choosing parametric (paired/unpaired t-test) or nonparametric (Wilcoxon/Mann–Whitney) based on Shapiro–Wilk normality. Applies BH-FDR across paired tests (and uses Tukey/Dunn adjustments for unpaired cases as implemented), then exports a *_results.csv with corrected p-values and stars.

**PCA_Microglia**
Standardizes numeric features from a single CSV, runs PCA (2D) and UMAP (2D), and clusters the scaled data with K-means (k=4). Saves an annotated CSV with embeddings and cluster labels, prints variance/loadings, computes cluster centers and distances, and plots the UMAP colored by cluster.

**PCA_Microglia_Prep**
Scans one or more directories for specific Imaris output CSV patterns, reads each file (skipping header rows), and keeps only selected ID + measurement columns. Renames FilamentID→ID where needed, outer-merges everything into one wide table, and saves a cleaned merged CSV.

**PCA_Spines**
Selects defined spine morphometric features, standardizes them, assigns K-means clusters (k=4), and computes PCA/UMAP embeddings with PCA loadings exported to CSV. It also classifies spines into Stubby/Filopodia/Mushroom/Thin via simple rules, then exports the annotated dataset plus per-filament and cluster×spine-type count summaries.

**PCA_withPermutationTesting**
Concatenates all CSVs in a directory, assigns coarse clusters from the first two ID characters, and runs PCA (2D) followed by UMAP on PCA scores. Quantifies separation by computing cluster-centroid distances in PCA space and uses a feature-wise permutation test to generate null distance distributions and p-values, saving results and plotting UMAP plus observed-vs-null histograms.

**Photometry**
Loads raw 470 nm (Ca-dependent) and 415 nm (isosbestic) traces for two fibres, subtracts 470 background, and constructs time in seconds (assumes 30 fps). Fits a bi-exponential to the isosbestic drift, linearly scales it to the calcium channel, subtracts this baseline to detrend, and exports the processed signal.

**QuPath_CellCounts**
Cleans a QuPath/Cellfinder per-structure summary by removing excluded structures and simplifying region names by collapsing subregions into parent labels. Computes cells/mm³, optionally collapses MPOA subregions into a single MPOA row, sorts by density, and exports a compact Excel summary.
