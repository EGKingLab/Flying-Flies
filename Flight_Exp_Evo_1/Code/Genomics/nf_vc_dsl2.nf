#!/usr/bin/env nextflow
nextflow.enable.dsl=2
print("""\
|-------------------------------------------------------------------------------------------------|
|                               MU IRCF VARIANT CALLING PIPELINE                                  |
|-------------------------------------------------------------------------------------------------|
|                                                                                                 |
|             ("`-''-`").___..--''"`-._            ("`-''-`").___..--''"`-._                      |
|              `o_ o  )   `-.  (     ).`-.__.`)     `o_ o  )   `-.  (     ).`-.__.`)              |
|              (_Y_.)'  ._   )  `._ `. ``-..-'      (_Y_.)'  ._   )  `._ `. ``-..-'               |
|              _`--'_..-_/  /--'_.' .'              _`--'_..-_/  /--'_.' .'                       |
|          (il).-''  (li).'` ((!.-'`            (il).-''  (li).'` ((!.-'`                         |
|                                                                                                 |
|  Command-line parameters:                                                                       |
|-------------------------------------------------------------------------------------------------|
|                                                                                                 |
|params.genome  = ["ref/genome.fa"]                                                               |
|params.reads   = ["/cluster/VAST/bac/active-projects/variantCalling_mysk/giab_data/*_end{1,2}.fq"]                                        |
|params.pl      = ["illumina"]                                                                    |
|params.aligner = ["minimap2"]                                                                    |
|-------------------------------------------------------------------------------------------------|

""")

/*
    INPUT PARAMETERS AND READ CHANNEL CREATION:
 */

//LabID,FCID,Lane,Library,SampleID,Platform,R1,R2
Channel
    .fromPath(params.sample_sheet)
    .splitCsv(header: true)
    .map{ row-> tuple(row.LabID, row.SampleID, row.FCID, row.Library, row.Lane, row.Platform, file(row.R1), file(row.R2) )}
    .set { reads_ch }

/*


                                            PIPELINE START


 */

 process split_intervals {
    
    // 1. Create an index (.fai) of the reference genome in case it wasn't provided with samtools faidx.
    // 2. Create a dictionary (.dict) of the reference genome in case it wasn't provided with gatk CreateSequenceDictionary.
    // 3. Use gatk SplitIntervals to split the reference into evenly sized chunks. Defaults to 20 scatter gathers.
    // Note: uses BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW, see GATK documentation for info.
    //       in short, avoids splitting input intervals -- that is, the set of input intervals is split, but individual intervals are left intact.
    // Hong: add if to check if files or folders exist before run the jobs
    
    publishDir "${baseDir}/intervals", mode:'copy'
    
    input:
    path genome
    path make_intervals

    output:
    path "${genome}.fai"
    path "${genome.baseName}.dict"
    path "interval_scatter48/*-scatter.interval_list"

    //TODO
    //May need to play with subdivision mode as BALANCING_WITHOUT_INTERVAL_SUBDIVISION is very unbalanced due to chrUn..so maybe INTERVAL_SUBDIVISION
    //without switching every time we're splitting the last few jos in a process will hang up so only minimally parallelizes 
    //BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW balances bp per interval but still many more (small) intervals in some scatters 
    
    script:
    """
    samtools faidx -o "${genome}.fai" ${params.genome} 

    gatk --java-options "-Xmx10g -Xms1g" \
    CreateSequenceDictionary  \
        -R ${genome} \
        -O ${genome.baseName}.dict

    Rscript ${make_intervals} ${genome.baseName}.dict ${genome.baseName}.interval_list

    gatk --java-options "-Xmx10g -Xms1g" \
        SplitIntervals \
        -R ${genome} \
	-L ${genome.baseName}.interval_list \
	--subdivision-mode INTERVAL_SUBDIVISION \
        --extension -scatter.interval_list \
	--scatter-count 64 \
        -O interval_scatter48
    """
}

process bwa_index {

    publishDir "${baseDir}/reference", mode:"copy"

    input:
    path genome

    output:
    path "${genome.getName()}"
    path "${genome.getName()}.{ann,amb,bwt,pac,sa}"

    script:
    """
    mkdir -p ${baseDir}/reference
    rsync ${genome} ${baseDir}/reference
    bwa index ${genome} ${baseDir}/reference/
    """
}

process map_reads {

    // 1. Create readgroup with unique sample id, library, fcid, platform, and lane info to keep reads well annotated.
    // 2. If the input files are unmapped bams, convert them to fastq, then query sort them for MarkDuplicates.
    // 3. If the input files are fastqs map them to the reference, query-sorted.
    // If you can't populate any fields in the input csv (flowcell id (FCID), library, labid, lane, etc.) you need to change
    // the input tuple here as well as in the channel that sets reads_ch (see preprocessing)

    // TODO: Clean up naming conventions so less specific in expectations of Lane, FCID, etc.

    //echo true
    //publishDir "${baseDir}/aligned_reads", mode:'copy'

    input:
    path genome
    tuple val(labid), val(sample_id), val(fcid), val(library), val(lane), val(platform), path(read1), path(read2)

    output:
    tuple val(sample_id), val("${sample_id}_${fcid}_${library}_${lane}"), path("${sample_id}_${fcid}_${library}_${lane}_aligned.bam")

    script:
    readGroup = "@RG\\tID:${library}.${fcid}.${lane}\\tPU:${library}.${fcid}.${lane}\\tSM:${sample_id}\\tPL:${platform}\\tLB:${library}"
    if( read1.extension == "bam" )
        """
        echo "found a bam..converting to fastq"
        samtools collate --threads ${params.cpus_mapping} -o collated.bam ${read1}
        samtools fastq -F 256 --threads ${params.cpus_mapping} -1 read_R1.fastq.gz -2 read_R2.fastq.gz --reference ${params.genome} -@ ${params.cpus_mapping} -0 /dev/null -s /dev/null -n collated.bam
        minimap2 -t ${params.cpus_mapping} -2 -R \"${readGroup}\" -ax sr ${params.genome} "read_R1.fastq.gz" "read_R2.fastq.gz" > "${sample_id}_${fcid}_${library}_${lane}_aligned.sam"
        samtools sort --threads ${params.cpus_mapping} -m 4G -o "${sample_id}_${fcid}_${library}_${lane}_aligned.bam" -n "${sample_id}_${fcid}_${library}_${lane}_aligned.sam"
        rm "${sample_id}_${fcid}_${library}_${lane}_aligned.sam" collated.bam read_R1.fastq.gz read_R2.fastq.gz
        """
    else
        """
        echo "production minimap2 mode"
        minimap2 -t ${params.cpus_mapping} -2 -R \"${readGroup}\" -ax sr ${params.genome} ${read1} ${read2} > "${sample_id}_${fcid}_${library}_${lane}_aligned.sam"
        samtools sort --threads ${params.cpus_mapping} -m 2G -o "${sample_id}_${fcid}_${library}_${lane}_aligned.bam" -t RG "${sample_id}_${fcid}_${library}_${lane}_aligned.sam"
        rm "${sample_id}_${fcid}_${library}_${lane}_aligned.sam"
        """
}

process mark_duplicates {

   // Mark duplicate reads, but don't remove them from the reads that are used, in query sorted input using MarkDuplicatesSpark
   // Neegatk_nf_dsl2_2.outd to be query sorted for ~2x speedup. The bam should be coordinate sorted after this step, but we funnel everything through a full coordinate sort in the next process anyways. Could inv                                          estigate this and remove next process potentially.
   // Locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA.
   // Works by comparing sequences in the 5 prime positions of both reads and read-pairs in a SAM/BAM file.
   // After duplicate reads are collected, the tool differentiates the primary and duplicate reads using an algorithm that ranks reads by the sums of their base-quality scores.
   // Requires holding all of the readnames in memory while it groups read information, machine configuration and starting sort-order impact tool performance.
   // See: https://gatk.broadinstitute.org/hc/en-us/articles/360037224932-MarkDuplicatesSpark

    //publishDir "${baseDir}/dedup_bam", mode:'copy'

    input:
    tuple val(subject), val(unique_id), path(aligned_bams)

    output:
    tuple val(subject), val(unique_id), path("${unique_id}_dedup.bam")
    path "${unique_id}_dedup_metrics.txt"

    script:
    """
    gatk --java-options "-Xmx120G -Xms120g -Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_samtools=true" MarkDuplicates \
        --INPUT ${aligned_bams} \
        --OUTPUT ${unique_id}_dedup.bam \
        --METRICS_FILE ${unique_id}_dedup_metrics.txt \
        --ASSUME_SORT_ORDER queryname
    """
}

process sort_bams {

    // Sort the bams into coordinate-order coming out of MarkDuplicatesSpark.
    // 1. samtools sort sort the bam
    // 2. samtools index index the bam
    // We don't put this in the previous process so that we're not running anything without great multi-threading that could take a decent chunk of time. (Large performance gains)
    // HA: use --threads ${params.cpus_mapping} -m 4G for samtools sort to match the config file

    input:
    tuple val(subject), val(unique_id), path(dedup_bam)

    output:
    tuple val(subject), val(unique_id), path("${unique_id}_sorted.bam"), path("${unique_id}_sorted.bam.bai")

    script:
    """
    samtools sort --threads ${params.cpus_mapping} -m 2G $dedup_bam > "${unique_id}_sorted.bam"
    samtools index -@ ${params.cpus} "${unique_id}_sorted.bam" "${unique_id}_sorted.bam.bai"
    """
}

process merge_bams {

    //publishDir "${baseDir}/merged_bams", mode:'copy'
    //echo true

    input:
    tuple val(subject), val(unique_id), path(bams), path(bais)

    output:
    tuple val(subject), path("${subject}_merged.bam"), path("${subject}_merged.bam.bai")

    shell:
    """
    echo !subject !bams
    mkdir -p /tmp/$SLURM_JOB_ID/!{subject}
    find . -name '*_sorted.bam' |
        parallel -j8 -N4095 -m --tmpdir /tmp/$SLURM_JOB_ID/!{subject} --files samtools merge -u - |
        parallel --tmpdir /tmp/$SLURM_JOB_ID/!{subject} --xargs samtools merge -@ 8 !{subject}_merged.bam {}";" rm {}
    samtools index -@ !{params.cpus} "!{subject}_merged.bam" "!{subject}_merged.bam.bai"
    rm -rf /tmp/$SLURM_JOB_ID/!{subject}
    """
}

process getMetrics {

    //publishDir "${baseDir}/metrics", mode:'copy'

    //if the user doesn't want (massive) stat files then skip this step
    when: params.col_metrics

    input:
    path genome
    tuple val(subject), val(unique_id), path(sorted_dedup_bam), path(bai)

    output:
    tuple val(subject), \
        path("${subject}_alignment_metrics.txt"), \
        path("${subject}_insert_metrics.txt"), \
        path("${subject}_insert_size_histogram.pdf"), \
        path("${subject}_depth_out.txt") \

    script:
    """
    gatk --java-options "-Xmx10g -Xms10g" \
    CollectAlignmentSummaryMetrics \
        R=${params.genome} \
        I=${sorted_dedup_bam} \
        O=${subject}_alignment_metrics.txt

    gatk --java-options "-Xmx10g -Xms10g" \
    CollectInsertSizeMetrics \
        R=${params.genome} \
        I=${sorted_dedup_bam} \
        O=${subject}_insert_metrics.txt \
        H=${subject}_insert_size_histogram.pdf

    samtools depth -a ${sorted_dedup_bam} > ${subject}_depth_out.txt
    samtools coverage -A ${sorted_dedup_bam} > ${subject}_coverage.txt
    """
}

process haplotypecaller {

    //publishDir "${baseDir}/haplotypecaller", mode:'copy'
    //tag {subject+"-"+file(Interval).baseName}

    input:
    tuple val(subject), path(bam), path(bai), path(Interval)
    path genome
    path fai
    path dict

    output:
    tuple val(subject), path("${subject}.${int_tag}.g.vcf"), path("${subject}.${int_tag}.g.vcf.idx")

    script:
    int_tag=Interval.baseName
    """
    gatk --java-options "-Xmx20g -Xms20g" \
    HaplotypeCaller \
        --gatk-config-file ${baseDir}/config_gatk.properties \
        --native-pair-hmm-threads 10 \
        --sample-ploidy 8 \
        -ERC GVCF \
        -L ${Interval} \
        -R ${genome} \
        -I ${bam} \
        -O ${subject}.${int_tag}.g.vcf
    """
}

process merge_GVCFs {

    publishDir "${baseDir}/gvcfs", mode: 'copy'
    tag { subject }

    input:
    tuple val(subject), path(gvcfs), path(gvcfidxs)
    path genome
    path index

    output:
    path "${subject}.g.vcf"
    path "${subject}.g.vcf.idx"

    script:
    """
    mkdir -p /tmp/$subject
    gatk --java-options "-Xmx20g -Xms20g" \
    SortVcf \
        ${gvcfs.collect { "-I $it " }.join()} \
        --TMP_DIR /tmp/$subject \
        -O ${subject}.g.vcf
    rm -rf /tmp/$subject
    """
}

process parabricks {

    publishDir "${baseDir}/gcvfs", mode: 'copy'
    tag { sample_id }

    input:
    path genome
    tuple val(labid), val(sample_id), val(fcid), val(library), val(lane), val(platform), path(read1), path(read2)

    output:
    path "${sample_id}.g.vcf"
    path "${sample_id}.g.vcf.idx"

    script:
    readGroup = "@RG\\tID:${library}.${fcid}.${lane}\\tPU:${library}.${fcid}.${lane}\\tSM:${sample_id}\\tPL:${platform}\\tLB:${library}"
    """
    singularity exec -B /cluster/VAST/bac/ --nv /cluster/software/containers/parabricks/4.0.1-1/parabricks_4.0.1-1.sif \
      pbrun germline \
      --ref ${baseDir}/reference/${genome.getName()} \
      --in-fq ${read1} ${read2} "${readGroup}" \
      --out-bam ${sample_id}.bam \
      --gvcf \
      --out-variants ${sample_id}.g.vcf \
      --tmp-dir ${sample_id}_temp \
      --num-gpus 2

    rm -rf ${sample_id}_temp
    """
}

process GenomicsDBImport {
    
    //publishDir "${baseDir}/GenomicsDB", mode: 'copy'
    //tag {file(interval).simpleName}

    input:
    path gvcf
    path gvcf_idxs
    each interval

    output:
    tuple val(interval.simpleName), path("${params.cohort}-${interval.simpleName}")

    script:
    """
    mkdir -p /tmp/${params.cohort}-${interval.simpleName}
    gatk GenomicsDBImport --java-options "-Xmx120g -Xms120g -Djava.io.tmpdir=/tmp/${params.cohort}-${interval.simpleName}" \
        --gatk-config-file ${baseDir}/config_gatk.properties \
        --genomicsdb-shared-posixfs-optimizations true \
        --merge-input-intervals true \
        --genomicsdb-workspace-path ${params.cohort}-${interval.simpleName} \
        ${gvcf.collect { " -V $it " }.join()} \
        -L ${interval} \
        --reader-threads 20 \
        --tmp-dir "/tmp/${params.cohort}-${interval.simpleName}"
    rm -rf /tmp/${params.cohort}-${interval.simpleName}
    """
}

process GenotypeGVCFs {

    //publishDir "${baseDir}/GenotypedGVCFs", mode: 'copy', pattern: '*.{vcf,idx}'
    //tag { interval }

    input:
    tuple val(interval), path(workspace)
    path genome
    path faidx
    path dict

    output:
    path "${params.cohort}.${interval}.vcf"
    path "${params.cohort}.${interval}.vcf.idx"

    script:
    """
    mkdir -p /tmp/$SLURM_JOB_ID/GenotypeGVCFs/${interval}
    gatk --java-options "-Xmx40g -Xms40g" \
        GenotypeGVCFs \
        -R ${params.genome} \
        -O ${params.cohort}.${interval}.vcf \
        --genomicsdb-shared-posixfs-optimizations true \
        --tmp-dir /tmp/$SLURM_JOB_ID/GenotypeGVCFs/${interval} \
        -G StandardAnnotation \
	-V gendb://${workspace.baseName}
    rm -rf /tmp/$SLURM_JOB_ID/GenotypeGVCFs/${interval}
    """
}


process gather_vcfs {

    //publishDir "${baseDir}/result_VCF", mode: 'copy'

    input:
    path split_vcf
    path split_vcfidx
    path genome
    path dict

    output:
    path "${params.cohort}-result.vcf.gz"
    path "${params.cohort}-result.vcf.gz.tbi"

    script:
    """
    gatk --java-options "-Xmx40g -Xms40g" \
    SortVcf \
        ${split_vcf.collect { "-I $it " }.sort().join()} \
        --CREATE_INDEX true \
        -O ${params.cohort}-result.vcf.gz
    """
}

process SelectSNPs {

    //publishDir "${baseDir}/result_VCF", mode: 'copy'
    //tag { subject }

    input:
    path vcf
    path vcf_idx

    output:
    path "${params.cohort}-result.SNP.vcf.gz"
    path "${params.cohort}-result.SNP.vcf.gz.tbi"

    script:
    """
    gatk --java-options "-Xmx80G -Xms80g" \
        SelectVariants \
        ${vcf.collect { " -V $it " }.join()} \
        -select-type SNP \
        -O ${params.cohort}-result.SNP.vcf.gz
    """
}

process SelectInDels {
    
    //publishDir "${baseDir}/result_VCF", mode: 'copy'
    //tag { subject }

    input:
    path vcf
    path vcf_idx

    output:
    path "${params.cohort}-result.InDel.vcf.gz"
    path "${params.cohort}-result.InDel.vcf.gz.tbi"

    script:
    """
    gatk --java-options "-Xmx80G -Xms80g" \
        SelectVariants \
        ${vcf.collect { " -V $it " }.join()} \
        -select-type INDEL \
        -O ${params.cohort}-result.InDel.vcf.gz
    """
}

process fltmarkSNP {

    //publishDir "${baseDir}/result_VCF", mode: 'copy'
    //tag { subject }

    input:
    path vcf
    path vcf_idx

    output:
    path "${params.cohort}-result.fltmark.SNP.vcf.gz"
    path "${params.cohort}-result.fltmark.SNP.vcf.gz.tbi"

    script:
    """
    gatk --java-options "-Xmx80G -Xms80g" \
        VariantFiltration \
        ${vcf.collect { " -V $it " }.join()} \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        -O ${params.cohort}-result.fltmark.SNP.vcf.gz
    """
}

process fltmarkInDel {

    publishDir "${baseDir}/result_VCF", mode: 'copy'
    tag { subject }

    input:
    path vcf
    path vcf_idx

    output:
    path "${params.cohort}-result.fltmark.InDel.vcf.gz"
    path "${params.cohort}-result.fltmark.InDel.vcf.gz.tbi"

    script:
    """
    gatk --java-options "-Xmx80G -Xms80g" \
        VariantFiltration \
        ${vcf.collect { " -V $it " }.join()} \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "FS > 200.0" --filter-name "FS200" \
        -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
        -O ${params.cohort}-result.fltmark.InDel.vcf.gz
    """
}

process hdfltSNP {

    publishDir "${baseDir}/result_VCF", mode: 'copy'
    tag { subject }

    input:
    path vcf
    path vcf_idx

    output:
    path "${params.cohort}-result.hdflt.SNP.vcf.gz"
    path "${params.cohort}-result.hdflt.SNP.vcf.gz.tbi"

    script:
    """
    bcftools view -f PASS ${vcf.collect { " $it " }.join()} -O z -o ${params.cohort}-result.hdflt.SNP.vcf.gz
    gatk --java-options "-Xmx80G -Xms80g" IndexFeatureFile -I ${params.cohort}-result.hdflt.SNP.vcf.gz
    """
}

process hdfltInDel {

    publishDir "${baseDir}/result_VCF", mode: 'copy'
    tag { subject }

    input:
    path vcf
    path vcf_idx

    output:
    path "${params.cohort}-result.hdflt.InDel.vcf.gz"
    path "${params.cohort}-result.hdflt.InDel.vcf.gz.tbi"

    script:
    """
    bcftools view -f PASS ${vcf.collect { " $it " }.join()} -O z -o ${params.cohort}-result.hdflt.InDel.vcf.gz
    gatk --java-options "-Xmx80G -Xms80g" IndexFeatureFile -I ${params.cohort}-result.hdflt.InDel.vcf.gz
    """
}

process MergeSNPInDel {

    publishDir "${baseDir}/result_VCF", mode: 'copy'
    tag { subject }

    input:
    path snpvcf
    path snpvcf_idx
    path indelvcf
    path indelvcf_idx

    output:
    path "${params.cohort}-result.hdflt.vcf.gz"
    path "${params.cohort}-result.hdflt.vcf.gz.tbi"

    script:
    """
    gatk --java-options "-Xmx80g -Xms80g" \
        SortVcf \
        ${snpvcf.collect { " -I $it " }.join()} \
        ${indelvcf.collect { " -I $it " }.join()} \
        --CREATE_INDEX true \
        -O ${params.cohort}-result.hdflt.vcf.gz
    """
}


workflow {
	// prepare reference files
        split_intervals(params.genome, params.Rinterval)

	// check for GPU flag
	if( params.gpu ) {

		// call BWA
		bwa_index(params.genome)

		// call parabricks
		parabricks(bwa_index.out[0], reads_ch)
	}
	else {

		// map reads, mark duplicates, sort BAMs, and get metrics
	        map_reads(params.genome, reads_ch) | mark_duplicates
	        mark_duplicates.out[0] | sort_bams
        	getMetrics(params.genome, sort_bams.out[0])

	        // sort all BAMs and merge them
        	sort_bams.out[0] | groupTuple | merge_bams

	        // call haplotypes and sort gVCFs
        	haplotypecaller(merge_bams.out[0].combine(split_intervals.out[2].flatten()), params.genome, split_intervals.out[0], split_intervals.out[1])
	        merge_GVCFs(haplotypecaller.out[0] | groupTuple, params.genome, split_intervals.out[0])
	}

	// check if whole workflow should run
        if( params.run_all ) {

                // DB import
		params.gpu ? GenomicsDBImport(parabricks.out[0].collect(), parabricks.out[1].collect(), split_intervals.out[2]) : GenomicsDBImport(merge_GVCFs.out[0].collect(), merge_GVCFs.out[1].collect(), split_intervals.out[2])

                // genotype and gather VCFs
                GenotypeGVCFs(GenomicsDBImport.out[0], params.genome, split_intervals.out[0], split_intervals.out[1])
                gather_vcfs(GenotypeGVCFs.out[0].collect(), GenotypeGVCFs.out[1].collect(), params.genome, split_intervals.out[1])

                // select and filter SNPs
                SelectSNPs(gather_vcfs.out[0].collect(), gather_vcfs.out[1].collect())
                fltmarkSNP(SelectSNPs.out[0].collect(), SelectSNPs.out[1].collect())
                hdfltSNP(fltmarkSNP.out[0].collect(), fltmarkSNP.out[1].collect())

                // select and filter indels
                SelectInDels(gather_vcfs.out[0].collect(), gather_vcfs.out[1].collect())
                fltmarkInDel(SelectInDels.out[0].collect(), SelectInDels.out[1].collect())
                hdfltInDel(fltmarkInDel.out[0].collect(), fltmarkInDel.out[1].collect())

                // merge SNPs and indels
                MergeSNPInDel(hdfltSNP.out[0].collect(), hdfltSNP.out[1].collect(), hdfltInDel.out[0].collect(), hdfltInDel.out[1].collect())
        }
}

