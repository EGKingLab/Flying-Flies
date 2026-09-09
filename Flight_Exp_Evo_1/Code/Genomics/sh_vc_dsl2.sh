#!/bin/bash
#-------------------------------------------------------------------------------
#  SBATCH CONFIG
#-------------------------------------------------------------------------------
## resources
#SBATCH -p general
#SBATCH -N 1
#SBATCH --cpus-per-task 4
#SBATCH --mem=5G
#SBATCH -t 1-00:00
#SBATCH --account=general
#
## labels and outputs
#SBATCH -J cpu_gatk
#SBATCH -o cpu_gatk.out
#
## notifications
##SBATCH --mail-user=svqht@umsystem.edu,may96c@umsystem.edu
##SBATCH --mail-type=END,FAIL
#
#-------------------------------------------------------------------------------
module load miniconda3
module load samtools
module load nextflow
module load minimap2
##module load gatk
module load hisat2
module load r
export PATH=/cluster/pixstor/kingeg-lab/BAC_nf_vc_ex/bin:$PATH

nextflow run nf_vc_dsl2.nf -resume \
	--genome /cluster/pixstor/kingeg-lab/FlyingFlies_vc/dmel-all-chromosome-r6.60.fasta \
	--sample_sheet /cluster/pixstor/kingeg-lab/FlyingFlies_vc/king_vc_samplesheet.csv \
	--bqsr false \
	--cohort fly \
	--cpus_mapping 22 \
	--cpus 6 \
	--col_metrics false \
	--gpu false \
	--run_all true \
	-with-trace logs/trace_cpu.txt \
        -with-report logs/report_cpu.html \
        -with-timeline logs/timeline_cpu.html

