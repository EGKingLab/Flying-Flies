---
editor_options: 
  markdown: 
    wrap: 72
---

# Readme

The folders are arranged in the following manner:

**Big_data**: contains large genomic data files (not raw genome fastq
files).

**Raw_data**: contains smaller raw data files, mainly phenotypic data.
Processed_data: ..? (delete?)

**Flight_Exp_Evo_1/Code**: contains all the code and scripts used to
process raw data and generate plots and figures.

**Flight_Exp_Evo_1/Code/Phenotypes**: contains all the code and scripts
for phenotypic data

**Flight_Exp_Evo_1/Code/Genomics**: contains all the code and scripts
for genomic data

**Flight_Exp_Evo_1/Plots**: contains all the plots used in the
manuscript

**Flight_Exp_Evo_1/Tables**: contains all the tables used in the
manuscript

## Instructions for generating figures and tables for phenotypic data

### Generating plots for the evolution of flight performance (figure 2)

Run the following QMD file:

`Flight_Exp_Evo_1/Code/Phenotypes/Flight_perf_evol.qmd`

The script reads data stored in the Raw_data folder and generates figure
2 from the manuscript in:

`Flight_Exp_Evo_1/Plots/flight_evol_fig.jpeg`

As well as the supplementary figure X in:

`Flight_Exp_Evo_1/Plots/flight_evol_scatter_all.jpeg`

An RDS file containing the list of all tables is saved in:

`Flight_Exp_Evo_1/Tables/fltevol_tables.RDS`

### Generating plots for the factorial tests of stimuli on flight performance (figure 3)

Run the following QMD file:

`Flight_Exp_Evo_1/Code/Phenotypes/Factorial.qmd`

The script reads data stored in the Raw_data folder and generates figure
3 from the manuscript in:

`Flight_Exp_Evo_1/Plots/fact_flight_fig.jpeg`

An RDS file containing the list of all tables is saved in:

`Flight_Exp_Evo_1/Tables/factflt_tables.RDS`

### Generating plots for the activity measurement assay (figure 4)

Run the following file *first*:
`Flight_Exp_Evo_1/Code/Phenotypes/PreProcess_Activity.R`

Next, run the QMD file: `Flight_Exp_Evo_1/Code/Phenotypes/Activity.qmd`

The script reads data stored in the Raw_data folder and generates figure
4 from the manuscript in:

`Flight_Exp_Evo_1/Plots/act_dist_fig.jpeg`

An RDS file containing the list of all tables is saved in:

`Flight_Exp_Evo_1/Tables/act_tables.RDS`

### Generating plots for wing morphometrics (figure 5)

Render the quarto file `Wing_Landmarking/Landmark_analysis.qmd`, which carries out the analysis and saves figure 5.

### Generating plots for life-history data: dry mass, survivorship and fecundity (figure 6)

#### Dry mass

Run the following QMD file:

`Flight_Exp_Evo_1/Code/Phenotypes/DryMass.qmd`

The script reads data stored in the Raw_data folder and generates an RDS
file for figure 6A:

`Flight_Exp_Evo_1/Plots/ggplot_files/mass_plot.RDS`

An RDS file containing the list of all tables is saved in:

`Flight_Exp_Evo_1/Tables/mass_tables.RDS`

#### Egg-to-adult survivorship

Run the following QMD file:

`Flight_Exp_Evo_1/Code/Phenotypes/Survivorship.qmd`

The script reads data stored in the Raw_data folder and generates an RDS
file for figure 6B:

`Flight_Exp_Evo_1/Plots/ggplot_files/surv_plot.RDS`

An RDS file containing the list of all tables is saved in:

`Flight_Exp_Evo_1/Tables/surv_tables.RDS`

#### Fecundity

Run the following QMD file:

`Flight_Exp_Evo_1/Code/Phenotypes/Fecundity.qmd`

The script reads data stored in the Raw_data folder and generates an RDS
file for figure 6C:

`Flight_Exp_Evo_1/Plots/ggplot_files/fec_plot.RDS`

An RDS file containing the list of all tables is saved in:

`Flight_Exp_Evo_1/Tables/fec_tables.RDS`

The supplementary figure X for fecundity data is saved in:

`Flight_Exp_Evo_1/Plots/fec_supp_fig.jpeg`

#### Compiling life-history data figure 6:

Run the following QMD file:
`Flight_Exp_Evo_1/Code/Phenotypes/Pheno_figures.qmd`

The script compiles the life-history figure from the RDS files of plots
saved in `Flight_Exp_Evo_1/Plots/ggplot_files/`

The figure is saved in:

`Flight_Exp_Evo_1/Plots/LH_fig.jpeg`

### Adjusting P-values and generating tables

Make sure all the table RDS files mentioned above have been saved in
`Flight_Exp_Evo_1/Tables/`

Run the following QMD file:
`Flight_Exp_Evo_1/Code/Phenotypes/Tables_and_PAdjust.qmd`

This script adjusts the P-values and generates the following .docx files
containing compiled tables in `Flight_Exp_Evo_1/Tables/`:
`activity_tables.docx`; `factflight_tables.docx`;
`fecundity_tables.docx`; `flightevol_tables.docx`;
`drymass_tables.docx`; `survivorship_tables.docx`

## Instructions for generating figures and tables for genomic data

### Aligning, cleanup and variant calling, parsing SNP tables

The nextflow code for alignment and variant calling are present in:
`Flight_Exp_Evo_1/Code/Genomics/nf_vc_dsl2.nf`

The following SBatch file was used to run the nextflow code:
`Flight_Exp_Evo_1/Code/Genomics/sh_vc_dsl2.sh`

The generated .vcf files are parsed with the following files to generate
SNP tables: `Flight_Exp_Evo_1/Code/Genomics/SNPtable_parseVCF.pl`

The SNP table is saved as a .txt file in:
`Big_data/FlyingFlies_SNPtable.txt`

### Filtering based on read depth and allele frequency

The following script filters the SNP data, and generates supplementary
figure X: `Flight_Exp_Evo_1/Code/Genomics/Clean_and_Filter_data.qmd`

The filtered data, to be used for all downstream analyses, is saved in:
`Big_data/FF_genomes_filtered.RDS`

Supplementary figure X, containing the upset plot for filtered data and
the read depth distribution is saved in:
`Flight_Exp_Evo_1/Plots/f_upset_depth_fig.jpeg`

### PCA (supplementary figure X)

Run the following QMD file: `Flight_Exp_Evo_1/Code/Genomics/PCA.qmd`

This script runs the PCA on the filtered data
`Big_data/FF_genomes_filtered.RDS`

It saves supplementary figure X, containing the biplot and the screeplot
in: `Flight_Exp_Evo_1/Plots/PCA_fig.jpeg`

### Linear model: SNP differentiation analysis

Run the .txt file:
`Flight_Exp_Evo_1/Code/Genomics/multi_looper_SNP.txt`

which loops over all chromosome arm names and 
runs site-wise analysis using the code in:
`Flight_Exp_Evo_1/Code/Genomics/SNP_multi_loop.R`

This saves the following files:
`Big_data/SNP_diff2L.RDS`; `Big_data/SNP_diff2R.RDS`; `Big_data/SNP_diff3L.RDS`;
`Big_data/SNP_diff3R.RDS`; `Big_data/SNP_diffX.RDS`

Next, run the QMD file:
`Flight_Exp_Evo_1/Code/Genomics/SNP_differentiation.qmd`

This joins the saved data files for the 5 chromosome arms, and performs P-value
FDR correction.
All P-values from linear model are saved in:
`Big_data/lm_adjpvals_SNPs.RDS`
significant sites (P-value < 0.05) are saved in:
`Big_data/lm_targetSNPs.RDS`

A plot for frequency of significant sites is saved as an RDS file
(to be integrated with gene mapping later)
`Flight_Exp_Evo_1/Plots/ggplot_files/sigsites_freq_plt.RDS`

A compiled figure of read depth of significant sites, and a scatter plot of
frequency of significant sites and read depth of significant sites is saved
(supplementary figure X):
`Flight_Exp_Evo_1/Plots/depth_freq_SNPfig.jpeg`

### Gene diversity and Fst calculations using Grenedalf
(Figure 8 and supplementary figures X, Y)

Run the QMD file:
`Flight_Exp_Evo_1/Code/Genomics/grenedalf_genediversity_Fst.qmd`

This loads the data in `Big_data/FF_genomes_filtered.RDS`

and runs code for gene diversity and Fst using Grenedalf
https://github.com/lczech/grenedalf

A formatted input file is saved in:
`Big_data/FF_grenedalf_input.csv`

Gene diversity code saves the files:
`Big_data/grenedalf_diversity_out/diversity.csv`
`Big_data/grenedalf_diversity_out/gdalf_div_log.txt`
`Big_data/grenedalf_diversity_out/gdalf_div_errors.txt`

Fst (Hudson) code saves the files:
`Big_data/grenedalf_fsthudson_out/fst.csv`
`Big_data/grenedalf_fsthudson_out/gdalf_fsthudson_log.txt`
`Big_data/grenedalf_fsthudson_out/gdalf_fsthudson_errors.txt`

Fst (Nei) code saves the files:
`Big_data/grenedalf_fstnei_out/fst.csv`
`Big_data/grenedalf_fstnei_out/gdalf_fstnei_log.txt`
`Big_data/grenedalf_fstnei_out/gdalf_fstnei_errors.txt`

Figure 8 is saved in:
`Flight_Exp_Evo_1/Plots/div_fst_fig.jpeg`

supplementary figure X (gene diveristy distribution):
supplementary figure Y (Fst Hudson CONvsFP):

supplementary figure Z (heatmaps of Fst Hudson and Fst Nei):
`Flight_Exp_Evo_1/Plots/Fst_heatmp_fig.jpeg`

### Haplotype analysis


### Gene Ontology and Enrichment analyses
