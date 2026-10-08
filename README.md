# Readme
The folder structure is arranged in the following manner:
<tree>

**Big_data**: contains large genomic data files (not raw genome fastq files).
**Raw_data**: contains smaller raw data files, mainly phenotypic data
Processed_data: ..? (delete?)
**Flight_Exp_Evo_1/Code**: contains all the code and scripts used to process raw data and generate plots and figures.

**Flight_Exp_Evo_1/Code/Phenotypes**: contains all the code and scripts for phenotypic data

**Flight_Exp_Evo_1/Code/Genomics**: contains all the code and scripts for genomic data

**Flight_Exp_Evo_1/Plots**: contains all the plots used in the manuscript

**Flight_Exp_Evo_1/Tables**: contains all the tables used in the manuscript


## Instructions for generating figures and tables for phenotypic data

### Generating plots for the evolution of flight performance (figure 2)

Run the following QMD file:

Flight_Exp_Evo_1/Code/Phenotypes/Flight_perf_evol.qmd

The script reads data stored in the Raw_data folder and generates figure 2 from the manuscript in:

Flight_Exp_Evo_1/Plots/flight_evol_fig.jpeg

As well as the supplementary figure X in:

Flight_Exp_Evo_1/Plots/flight_evol_scatter_all.jpeg

An RDS file containing the list of all tables is saved in:

Flight_Exp_Evo_1/Tables/fltevol_tables.RDS




### Generating plots for the factorial tests of stimuli on flight performance (figure 3)

Run the following QMD file:

Flight_Exp_Evo_1/Code/Phenotypes/Factorial.qmd

The script reads data stored in the Raw_data folder and generates figure 3 from the manuscript in:

Flight_Exp_Evo_1//Plots/fact_flight_fig.jpeg

An RDS file containing the list of all tables is saved in:

Flight_Exp_Evo_1/Tables/factflt_tables.RDS



