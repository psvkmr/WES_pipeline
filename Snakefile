# WES pipeline with snakemake

#conda create -n gatk_snake -c bioconda gatk4 fastqc multiqc fastx_toolkit bwa samtools bcftools=1.9 ensembl-vep r-ggplot2 snakemake r-base
#conda activate gatk_snake

# test fastqs from wes dir
#mkdir fastqs
#for fastq in $(ls ../wes/fastqs/); do gunzip -c ../wes/fastqs/$fastq | head -n 50000 | gzip -c > fastqs/$fastq ; done 

import os

#SRA,FRR = glob_wildcards('fastqs/{sra}_{frr}.fastq')
configfile: 'config.yaml'

samples=config['samples']

#base=config['dirs']['base']
#fastqs=config['dirs']['fastqs']
#fastqcs=config['dirs']['fastqcs']
#multiqcs=config['dirs']['multiqcs']
#ref=config['dirs']['ref']
#bams=config['dirs']['bams']
#mpileups=config['dirs']['mpileups']
#gvars=config['dirs']['gvars']
#jvars=config['dirs']['jvars']

# any output that is not used for input into another rule must be included in rule all
rule all:
    input:
        expand('multiqcs/multiqc_{tool}.html', tool=['fastqc']),
        expand('bams/{sample}_fixed_bam_indexstats.txt', sample=samples),
        expand('mpileups/{sample}_mpileup_subset.tab', sample=samples),
        'jvars/variant_call_metrics.txt',
        'jvars/variant_call_metrics_gdb.txt'

rule fastqc:
    input:
        fq='fastqs/{sample}_{read}.fastq.gz'
    output:
        fqzip='fastqcs/{sample}_{read}_fastqc.zip',
        fqhtml='fastqcs/{sample}_{read}_fastqc.html'
    log:
        'logs/fastqc.log'
    shell:
        """
        date > {log}
        fastqc {input} -o fastqcs 2>> {log}
        """

rule multiqc_fastqc:
    input:
        expand('fastqcs/{sample}_{read}_fastqc.zip', sample=samples, read=['1','2'])
    output:
        'multiqcs/multiqc_fastqc/multiqc_report.html'
    log:
        'logs/multiqc_fastqc.log'
    params:
        dir='multiqcs/multiqc_fastqc'
    shell:
        """
        date > {log}
        multiqc {input} -o {params.dir} 2>> {log}
        """

# omitted read group due to format issue: -R {params}
rule bwa_map:
    input:
        fq1='fastqs/{sample}_1.fastq.gz',
        fq2='fastqs/{sample}_2.fastq.gz'
    output:
        'bams/{sample}.sam'
    log:
        'logs/bwa_{sample}.log'
    params:
        id='{sample}'
    shell:
        """
        date > {log}
        bwa mem -M -R "@RG\\tID:{params.id}\\tSM:{params.id}\\tLB:WES\\tPL:Illumina" \
        ../wes/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz {input.fq1} {input.fq2} > {output} 2>> {log}
        """

rule samtools_bam:
    input:
        rules.bwa_map.output
    output:
        'bams/{sample}.bam',
    log:
        'logs/samtools_{sample}.log'
    shell:
        """
        date > {log}
        samtools view -S -b {input} > {output} 2>> {log}
        """

rule samtools_sort:
    input:
        rules.samtools_bam.output
    output:
        'bams/{sample}_sorted.bam'
    log:
        'logs/samtools_{sample}.log'
    shell:
        """
        date >> {log}
        samtools sort -o {output} {input} 2>> {log}
        """

#rule samtools_index:
#    input:
#        'bams/{sample}_sorted.bam'
#    output:
#        'bams/{sample}_sorted.bam.bai'
#    log:
#        'logs/samtools_{sample}.log'
#    shell:
#        """
#        date >> {log}
#        samtools index {input} 2>> {log}
#        """

#rule gatk_add_rgs:
#    input:
#        rules.samtools_sort.output
#    output:
#        'bams/{sample}_readgrouped.bam'
#    log:
#        'logs/gatk_bams_{sample}.log'
#    params:
#        '{sample}'
#    shell:
#        """
#        date >> {log}
#        gatk AddOrReplaceReadGroups -I {input} -O {output} --RGPU 1 --RGID {params} --RGSM {params} --RGLB 'WES' --RGPL 'Illumina' 2>> {log}
#        """

rule gatk_mark_dups:
    input:
        rules.samtools_sort.output
    output:
        bam='bams/{sample}_deduped.bam',
        metrics='bams/{sample}_deduped_metrics.txt'
    log:
        'logs/gatk_bams_{sample}.log'
    shell:
        """
        date > {log}
        gatk MarkDuplicates -I {input} -O {output.bam} -M {output.metrics} 2>> {log}
        """

rule gatk_fix_mates:
    input:
        rules.gatk_mark_dups.output.bam
    output:
        'bams/{sample}_fixed.bam'
    log:
        'logs/gatk_bams_{sample}.log'
    shell:
        """
        date >> {log}
        gatk FixMateInformation -SO coordinate -I {input} -O {output} 2>> {log}
        """

rule gatk_bqsr:
    input:
        rules.gatk_fix_mates.output
    output:
        'bams/{sample}_qual_recal_table.txt'
    log:
        'logs/gatk_bams_{sample}.log'
    shell:
        """
        date >> {log}
        gatk BaseRecalibrator -R ../wes/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa -I {input} \
        --known-sites ../wes/ref/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38_nochr.vcf.gz \ 
        --known-sites ../wes/ref/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38_nochr.vcf.gz \
        --known-sites ../wes/ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138_nochr.vcf.gz \
        -O {output} 2>> {log}
        """

# cannot currently applybqsr as lacking read info, skip to bam checks


rule gatk_bam_index:
    input:
        rules.gatk_fix_mates.output
    output:
        'bams/{sample}_fixed.bam.bai'
    log:
        'logs/gatk_bams_{sample}.log'
    shell:
        """
        date >> {log}
        gatk BuildBamIndex -I {input} -O {output} 2>> {log}
        """

rule gatk_bam_stats:
    input:
        bam=rules.gatk_fix_mates.output,
        bai=rules.gatk_bam_index.output
    output:
        'bams/{sample}_fixed_bam_indexstats.txt'
    log:
        'logs/gatk_bams_{sample}.log'
    shell:
        """
        date >> {log}
        gatk BamIndexStats -I {input.bam} > {output} 2>> {log}
        """

rule samtools_mpileups:
    input:
        rules.gatk_fix_mates.output
    output:
        'mpileups/{sample}_mpileup_subset.tab'
    log:
        'logs/samtools_{sample}.log'
    shell:
        """ 
        date >> {log}
        samtools mpileup -f ../wes/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa -s {input} 2>> {log} | sed -n '500,515p' > {output}
        """

rule gatk_haplotype_caller:
    input:
        bams=rules.gatk_fix_mates.output, 
        bais=rules.gatk_bam_index.output
    output:
        'gvars/{sample}.g.vcf.gz'
    log:
        'logs/gatk_gvars_{sample}.log'
    shell:
        """
        date > {log}
        gatk HaplotypeCaller -ERC GVCF -R ../wes/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        -I {input.bams} -O {output} 2>> {log}
        """

rule gvcf_locate:
    input:
        'sample_ids.txt'
    output:
        paths='gvars/sample_gvcf_paths.txt',
    log:
        'logs/gvcf_locate.log'
    shell:
        """
        date > {log}
        while read line; 
        do 
            readlink -f gvars/$line.g.vcf.gz; 
        done < {input} > {output.paths} 2> {log}
        """

rule gvcf_paths:
    input:
        ids='sample_ids.txt',
        locate=rules.gvcf_locate.output
    output:
        'gvars/sample_gvcf_mapping.txt'
    log: 
        'logs/gvcf_paths.log'
    shell:
        """
        date > {log}
        paste {input.ids} {input.locate} > {output} 2> {log}
        """

rule gvcfs_as_args:
    input:
        rules.gvcf_locate.output
    output:
        'gvars/combine_gvcfs_args_list.txt'
    log:
        'logs/comb_gvcfs_args.txt'
    shell:
        """
        date > {log}
        sed 's/^/--variant /g' {input} > {output} 2>> {log}
        """


#mkdir -p genomicsdb
# get gatk group own interval list
rule gatk_genomicsdb:
    input:
        rules.gvcf_paths.output
    output:
        'gdb_done.txt'
    log:
        'logs/genomicsdb.log'
    params:
        gdb_workspace='genomicsdb'
    shell:
        """
        date > {log}
        gatk GenomicsDBImport -L 1 -L 2 -L 3 -L 4 -L 5 -L 6 -L 7 -L 8 -L 9 -L 10 -L 11 -L 12 -L 13 \
        -L 14 -L 15 -L 16 -L 17 -L 18 -L 19 -L 20 -L 21 -L 22 -L X -L Y \
        --sample-name-map {input} --genomicsdb-workspace-path {params.gdb_workspace} 2>> {log}
        touch gdb_done.txt
        """

# must be more efficient way of listing each gvcf file
rule gatk_combine_vcfs:
    input:
        rules.gvcfs_as_args.output
    output:
        'jvars/cohort.g.vcf.gz'
    log:
        'logs/gatk_combine_vcfs.log'
    shell:
        """
        date > {log}
        gatk CombineGVCFs -R ../wes/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        --arguments_file {input} -O {output} 2>> {log}
        """

rule gatk_genotype_gvcfs:
    input:
        rules.gatk_combine_vcfs.output
    output:
        'jvars/cohort.vcf.gz'
    log:
        'logs/gatk_genotype_gvcfs.log'
    shell:
        """
        date > {log}
        gatk GenotypeGVCFs -R ../wes/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa -V {input} \
        --annotate-with-num-discovered-alleles --create-output-variant-index -O {output} 2>> {log}
        """

rule gatk_genotype_gvcfs_gdb:
    output:
        'jvars/cohort_gdb.vcf.gz'
    log:
        'logs/gatk_genotype_gvcfs_gdb.log'
    shell:
        """
        date > {log}
        gatk GenotypeGVCFs -R ../wes/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa -V gendb://genomicsdb \
        --annotate-with-num-discovered-alleles --create-output-variant-index -O {output} 2>> {log}
        """

rule gatk_calculate_vqsr_snps:
    input:
        rules.gatk_genotype_gvcfs.output
    output:
        recal='jvars/cohort_snps_recal.txt',
        formats='jvars/cohort_snps_recal_format_plots.R',
        tranches='jvars/cohort_snps_recal_tranches_plots.txt'
    log:
        'logs/gatk_calculate_vqsr_snps.log'
    shell:
        """
        date > {log}
        gatk VariantRecalibrator \
        -R ../wes/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        -V {input} \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ../wes/ref/resources_broad_hg38_v0_hapmap_3.3.hg38_nochr.vcf.gz \
        --resource:omini,known=false,training=true,truth=false,prior=12.0 ../wes/ref/resources_broad_hg38_v0_1000G_omni2.5.hg38_nochr.vcf.gz \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 ../wes/ref/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38_nochr.vcf.gz \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ../wes/ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138_nochr.vcf.gz \
        -an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum -mode SNP -tranche 100.0 -tranche 99.0 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
        --rscript-file {output.formats} \
        --tranches-file {output.tranches} \
        -O {output.recal} 2>> {log} 
        """

rule gatk_calculate_vqsr_snps_gdb:
    input:
        rules.gatk_genotype_gvcfs_gdb.output
    output:
        recal='jvars/cohort_snps_recal_gdb.txt',
        formats='jvars/cohort_snps_recal_format_plots_gdb.R',
        tranches='jvars/cohort_snps_recal_tranches_plots_gdb.txt'
    log:
        'logs/gatk_calculate_vqsr_snps_gdb.log'
    shell:
        """
        date > {log}
        gatk VariantRecalibrator \
        -R ../wes/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        -V {input} \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ../wes/ref/resources_broad_hg38_v0_hapmap_3.3.hg38_nochr.vcf.gz \
        --resource:omini,known=false,training=true,truth=false,prior=12.0 ../wes/ref/resources_broad_hg38_v0_1000G_omni2.5.hg38_nochr.vcf.gz \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 ../wes/ref/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38_nochr.vcf.gz \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ../wes/ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138_nochr.vcf.gz \
        -an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum -mode SNP -tranche 100.0 -tranche 99.0 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
        --rscript-file {output.formats} \
        --tranches-file {output.tranches} \
        -O {output.recal} 2>> {log} 
        """

# some annotation caussing issues due to zero variance -an MQRankSum 
rule gatk_calculate_vqsr_indels:
    input:
        rules.gatk_genotype_gvcfs.output
    output:
        recal='jvars/cohort_indels_recal.txt',
        formats='jvars/cohort_indels_recal_format_plots.R',
        tranches='jvars/cohort_indels_recal_tranches_plots.txt'
    log:
        'logs/gatk_calculate_vqsr_indels.log'
    shell:
        """
        date > {log}
        gatk VariantRecalibrator \
        -R ../wes/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        -V {input} \
        --resource:mills,known=false,training=true,truth=true,prior=15.0 ../wes/ref/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38_nochr.vcf.gz \
        -an QD -an MQ -an FS -an SOR -mode INDEL \
        --rscript-file {output.formats} \
        --tranches-file {output.tranches} \
        -O {output.recal} 2>> {log} 
        """

# -an MQRankSum 
rule gatk_calculate_vqsr_indels_gdb:
    input:
        rules.gatk_genotype_gvcfs_gdb.output
    output:
        recal='jvars/cohort_indels_recal_gdb.txt',
        formats='jvars/cohort_indels_recal_format_plots_gdb.R',
        tranches='jvars/cohort_indels_recal_tranches_plots_gdb.txt'
    log:
        'logs/gatk_calculate_vqsr_indels_gdb.log'
    shell:
        """
        date > {log}
        gatk VariantRecalibrator \
        -R ../wes/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        -V {input} \
        --resource:mills,known=false,training=true,truth=true,prior=15.0 ../wes/ref/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38_nochr.vcf.gz \
        -an QD -an MQ -an SOR -an FS -mode INDEL \
        --rscript-file {output.formats} \
        --tranches-file {output.tranches} \
        -O {output.recal} 2>> {log} 
        """

rule gatk_apply_vqsr_snps:
    input:
        vcf=rules.gatk_genotype_gvcfs.output,
        recal=rules.gatk_calculate_vqsr_snps.output.recal,
        tranches=rules.gatk_calculate_vqsr_snps.output.tranches 
    output:
        'jvars/vqsr_snps.vcf.gz'
    log:
        'logs/gatk_apply_vqsr_snps.log'
    shell:
        """
        date > {log}
        gatk ApplyVQSR \
        -R ../wes/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        -V {input.vcf} \
        --truth-sensitivity-filter-level 99.0 \
        --tranches-file {input.tranches} \
        --recal-file {input.recal} \
        -mode SNP \
        -O {output} 2>> {log}
        """

rule gatk_apply_vqsr_snps_gdb:
    input:
        vcf=rules.gatk_genotype_gvcfs_gdb.output,
        recal=rules.gatk_calculate_vqsr_snps_gdb.output.recal,
        tranches=rules.gatk_calculate_vqsr_snps_gdb.output.tranches  
    output:
        'jvars/vqsr_snps_gdb.vcf.gz'
    log:
        'logs/gatk_apply_vqsr_snps_gdb.log'
    shell:
        """
        date > {log}
        gatk ApplyVQSR \
        -R ../wes/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        -V {input.vcf} \
        --truth-sensitivity-filter-level 99.0 \
        --tranches-file {input.tranches} \
        --recal-file {input.recal} \
        -mode SNP \
        -O {output} 2>> {log}
        """

rule gatk_apply_vqsr_indels:
    input:
        vcf=rules.gatk_apply_vqsr_snps.output,
        recal=rules.gatk_calculate_vqsr_indels.output.recal, 
        tranches=rules.gatk_calculate_vqsr_indels.output.tranches
    output:
        'jvars/vqsr_snps_indels.vcf.gz'
    log:
        'logs/gatk_apply_vqsr_indels.log'
    shell:
        """
        date > {log}
        gatk ApplyVQSR \
        -R ../wes/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        -V {input.vcf} \
        --truth-sensitivity-filter-level 99.0 \
        --tranches-file {input.tranches} \
        --recal-file {input.recal} \
        -mode INDEL \
        -O {output} 2>> {log}
        """

rule gatk_apply_vqsr_indels_gdb:
    input:
        vcf=rules.gatk_apply_vqsr_snps_gdb.output,
        recal=rules.gatk_calculate_vqsr_indels_gdb.output.recal, 
        tranches=rules.gatk_calculate_vqsr_indels_gdb.output.tranches
    output:
        'jvars/vqsr_snps_indels_gdb.vcf.gz'
    log:
        'logs/gatk_apply_vqsr_indels_gdb.log'
    shell:
        """
        date > {log}
        gatk ApplyVQSR \
        -R ../wes/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        -V {input.vcf} \
        --truth-sensitivity-filter-level 99.0 \
        --tranches-file {input.tranches} \
        --recal-file {input.recal} \
        -mode INDEL \
        -O {output} 2>> {log}
        """

rule index_vqsr_vcfs:
    input:
        rules.gatk_apply_vqsr_indels.output 
    output:
        'jvars/vqsr_snps_indels.vcf.gz.tbi'
    log:
        'logs/index_vqsr_bams.log'
    shell:
        """
        date > {log}
        tabix {input} -f > {output} 2>> {log}
        """

rule index_vqsr_vcfs_gdb:
    input:
        rules.gatk_apply_vqsr_indels_gdb.output 
    output:
        'jvars/vqsr_snps_indels_gdb.vcf.gz.tbi'
    log:
        'logs/index_vqsr_bams_gdb.log'
    shell:
        """
        date > {log}
        tabix {input} -f > {output} 2>> {log}
        """

rule gatk_genotype_posteriors:
    input:
        vcf=rules.gatk_apply_vqsr_indels.output,
        index=rules.index_vqsr_vcfs.output
    output:
        'jvars/geno_post.vcf.gz'
    log:
        'logs/gatk_geno_post.log'
    shell:
        """
        date > {log}
        gatk CalculateGenotypePosteriors -V {input.vcf} \
        -supporting ../wes/ref/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38_nochr.vcf. \
        -O {output} 2>> {log}
        """

rule gatk_genotype_posteriors_gdb:
    input:
        vcf=rules.gatk_apply_vqsr_indels_gdb.output,
        index=rules.index_vqsr_vcfs_gdb.output
    output:
        'jvars/geno_post_gdb.vcf.gz'
    log:
        'logs/gatk_geno_post_gdb.log'
    shell:
        """
        date > {log}
        gatk CalculateGenotypePosteriors -V {input.vcf} \
        -supporting ../wes/ref/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38_nochr.vcf \
        -O {output} 2>> {log}
        """

#For SNPs:
#QD &lt; 2.0
#MQ &lt; 40.0
#FS &gt; 60.0
#SOR &gt; 3.0
#MQRankSum &lt; -12.5
#ReadPosRankSum &lt; -8.0

#For indels:
#QD &lt; 2.0
#ReadPosRankSum &lt; -20.0
#InbreedingCoeff &lt; -0.8
#FS &gt; 200.0
#SOR &gt; 10.0

rule gatk_filter_vars:
    input:
        rules.gatk_genotype_posteriors.output
    output:
        'jvars/filtered_vars.vcf.gz'
    log:
        'logs/gatk_filt_vars.log'
    shell:
        """
        date > {log}
        gatk VariantFiltration -V {input} \
        --filter-expression "GQ < 20.0" --filter-name "GQ20" \
        --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
        -O {output} 2>> {log}
        """

rule gatk_filter_vars_gdb:
    input:
        rules.gatk_genotype_posteriors_gdb.output
    output:
        'jvars/filtered_vars_gdb.vcf.gz'
    log:
        'logs/gatk_filt_vars_gdb.log'
    shell:
        """
        date > {log}
        gatk VariantFiltration -V {input} \
        --filter-expression "GQ < 20.0" --filter-name "GQ20" \
        --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
        -O {output} 2>> {log}
        """


rule gatk_var_call_metrics:
    input:
        rules.gatk_filter_vars.output 
    output:
        'jvars/variant_call_metrics.txt'
    log:
        'logs/gatk_var_call_metrics.log'
    shell:
        """
        date > {log}
        gatk CollectVariantCallingMetrics -I {input} \
        --DBSNP ../wes/ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
        -SD ../wes/ref/Homo_sapiens.GRCh38.dna.primary_assembly.dict \
        -O {output} 2>> {log}
        """

rule gatk_var_call_metrics_gdb:
    input:
        rules.gatk_filter_vars_gdb.output 
    output:
        'jvars/variant_call_metrics_gdb.txt'
    log:
        'logs/gatk_var_call_metrics_gdb.log'
    shell:
        """
        date > {log}
        gatk CollectVariantCallingMetrics -I {input} \
        --DBSNP ../wes/ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
        -SD ../wes/ref/Homo_sapiens.GRCh38.dna.primary_assembly.dict \
        -O {output} 2>> {log}
        """


rule vep:
    input:
        rules.gatk_filter_vars.output 
    output:
        vcf='annotated/full_anno_vep.vcf.gz',
        stats='annotated/full_anno_vep_stats.html'
    log:
        'logs/vep.log'
    shell:
        """
        date > {log}
        vep -i {input} \
        --format vcf \
        --assembly GRCh38 \
        --fasta ../wes/vep_dir/homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
        --offline \
        --cache \
        --dir_cache ../wes/vep_dir \
        --check_existing \
        --everything \
        --per_gene \
        --pick_order rank \
        --nearest gene \
        --overlaps \
        --gene_phenotype \
        --show_ref_allele \
        --total_length \
        --clin_sig_allele 0 \
        --af_exac \
        --pubmed \
        --plugin CADD,../wes/ref/gnomad.genomes.r3.0.snv.tsv.gz \
        --plugin PostGAP,../wes/ref/postgap_GRCh38.txt.gz,ALL \
        --dir_plugins /home/k2142172/.vep/Plugins \
        --verbose \
        --vcf \
        --stats_file {output.stats} \
        -o {output.vcf} 2>> {log}
        """

rule vep_gdb:
    input:
        rules.gatk_filter_vars_gdb.output 
    output:
        vcf='annotated/full_anno_vep_gdb.vcf.gz',
        stats='annotated/full_anno_vep_stats_gdb.html'
    log:
        'logs/vep_gdb.log'
    shell:
        """
        date > {log}
        vep -i {input} \
        --format vcf \
        --assembly GRCh38 \
        --fasta ../wes/vep_dir/homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
        --offline \
        --cache \
        --dir_cache ../wes/vep_dir \
        --check_existing \
        --everything \
        --per_gene \
        --pick_order rank \
        --nearest gene \
        --overlaps \
        --gene_phenotype \
        --show_ref_allele \
        --total_length \
        --clin_sig_allele 0 \
        --af_exac \
        --pubmed \
        --plugin CADD,../wes/ref/gnomad.genomes.r3.0.snv.tsv.gz \
        --plugin PostGAP,../wes/ref/postgap_GRCh38.txt.gz,ALL \
        --dir_plugins /home/k2142172/.vep/Plugins \
        --verbose \
        --vcf \
        --stats_file {output.stats} \
        -o {output.vcf} 2>> {log}
        """

# validate variants