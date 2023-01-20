# WES pipeline with snakemake

# use bcftools annotate --rename-chrs

#conda create -n gatk_snake -c bioconda gatk4 fastqc multiqc fastx_toolkit bwa samtools bcftools=1.9 ensembl-vep r-ggplot2 snakemake r-base r-data.table r-tidyverse
#conda activate gatk_snake

#mkdir fastqs
#for fastq in $(ls fastqs/); do gunzip -c fastqs/$fastq | head -n 50000 | gzip -c > fastqs/$fastq ; done 

import os

#SRA,FRR = glob_wildcards('fastqs/{sra}_{frr}.fastq')
configfile: 'config.yaml'

samples=config['samples']
hg38=config['hg38']

# any output that is not used for input into another rule must be included in rule all
rule all:
    input:
        expand('multiqcs/multiqc_{tool}/multiqc_report.html', tool=['fastqc','samtools']),
        expand('bams/{sample}_fixed_bam_indexstats.txt', sample=samples),
        expand('mpileups/{sample}_mpileup_subset.tab', sample=samples),
        expand('jvars/cohort_{pipeline}.variant_calling_detail_metrics', pipeline=['cmb','gdb']),
        expand('annotated/impact_filtered_vep_{pipeline}.vcf.processed.tsv', pipeline=['cmb','gdb'])

rule fastqc:
    input:
        'fastqs/{sample}_{read}.fastq.gz'
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

# create index: bwa index Homo_sapiens_assembly38.fasta
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
        {hg38} {input.fq1} {input.fq2} > {output} 2>> {log}
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
        gatk BaseRecalibrator -R {hg38} -I {input} \
        --known-sites ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \ 
        --known-sites ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
        --known-sites ref/Homo_sapiens_assembly38.dbsnp138.vcf \
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
        samtools mpileup -f {hg38} -s {input} 2>> {log} | sed -n '500,515p' > {output}
        """

# for multiqc, idxstats file must have 'idxstat' somewhere in file name
rule samtools_stats:
    input:
        rules.gatk_fix_mates.output 
    output:
        stats='bams/samtools_{sample}_stats.txt',
        flg='bams/samtools_{sample}_flagstats.txt',
        idx='bams/samtools_{sample}_idxstats.txt'
    log:
        'logs/samtools_{sample}.log'
    shell:
        """
        date >> {log}
        samtools stats {input} > {output.stats} 2>> {log}
        samtools flagstat {input} > {output.flg} 2>> {log}
        samtools idxstats {input} > {output.idx} 2>> {log}
        """

rule multiqc_samtools:
    input:
        stats=expand(rules.samtools_stats.output.stats, sample=samples),
        flg=expand(rules.samtools_stats.output.flg, sample=samples),
        idx=expand(rules.samtools_stats.output.idx, sample=samples)
    output:
        'multiqcs/multiqc_samtools/multiqc_report.html'
    log:
        'logs/multiqc_samtools.log'
    params:
        dir='multiqcs/multiqc_samtools'
    shell:
        """
        date > {log}
        multiqc {input.stats} {input.flg} {input.idx} -o {params.dir} 2>> {log}
        """

# need to have fasta dict file created: gatk CreateSequenceDictionary -R {hg38}
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
        gatk --java-options "-Xmx8G -Xms8G" HaplotypeCaller -ERC GVCF -R {hg38} \
        -I {input.bams} -O {output} 2>> {log}
        """

rule gvcf_locate:
    input:
        ids='sample_ids.txt',
        vcfs=expand('gvars/{sample}.g.vcf.gz', sample=samples)
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
        done < {input.ids} > {output.paths} 2> {log}
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
        gatk GenomicsDBImport -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 \
        -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 \
        -L chr22 -L chrX -L chrY \
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
        gatk CombineGVCFs -R {hg38} \
        --arguments_file {input} -O {output} 2>> {log}
        """

rule gatk_genotype_gvcfs:
    input:
        rules.gatk_combine_vcfs.output
    output:
        'jvars/cohort_cmb.vcf.gz'
    log:
        'logs/gatk_genotype_gvcfs.log'
    shell:
        """
        date > {log}
        gatk GenotypeGVCFs -R {hg38} -V {input} \
        --annotate-with-num-discovered-alleles --create-output-variant-index -O {output} 2>> {log}
        """

rule gatk_genotype_gvcfs_gdb:
    input:
        'gdb_done.txt'
    output:
        'jvars/cohort_gdb.vcf.gz'
    log:
        'logs/gatk_genotype_gvcfs_gdb.log'
    shell:
        """
        date > {log}
        gatk GenotypeGVCFs -R {hg38} -V gendb://genomicsdb \
        --annotate-with-num-discovered-alleles --create-output-variant-index -O {output} 2>> {log}
        """

rule gatk_validate_variants:
    input:
        'jvars/cohort_{pipeline}.vcf.gz' 
    output:
        'jvars/validation_done_{pipeline}.txt'
    log:
        'logs/validate_vcfs_{pipeline}.log'
    shell:
        """
        date > {log}
        gatk ValidateVariants -R {hg38} -V {input} \
        --dbsnp ref/Homo_sapiens_assembly38.dbsnp138.vcf
        touch {output} 2>> {log}
        """

# -an MQRankSum -an DP 
# from gatk vqsr page: mapping quality has a very different distribution because it is 
# not a calibrated statistic, so in some cases it can destabilize the model
# Depth of coverage should not be used when working with exome datasets since there 
# is extreme variation in the depth to which targets are captured
rule gatk_calculate_vqsr_snps:
    input:
        'jvars/cohort_{pipeline}.vcf.gz'
    output:
        recal='jvars/cohort_snps_recal_{pipeline}.vcf',
        formats='jvars/cohort_snps_recal_format_plots_{pipeline}.R',
        tranches='jvars/cohort_snps_recal_tranches_plots_{pipeline}.txt'
    log:
        'logs/gatk_calculate_vqsr_snps_{pipeline}.log'
    shell:
        """
        date > {log}
        gatk VariantRecalibrator \
        -R {hg38} \
        -V {input} \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ref/hapmap_3.3.hg38.vcf.gz \
        --resource:omini,known=false,training=true,truth=true,prior=12.0 ref/1000G_omni2.5.hg38.vcf.gz \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ref/Homo_sapiens_assembly38.dbsnp138.vcf \
        -an QD -an FS -an SOR -an ReadPosRankSum -mode SNP -tranche 100.0 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
        --rscript-file {output.formats} \
        --tranches-file {output.tranches} \
        -O {output.recal} 2>> {log} 
        """

# some annotation caussing issues due to zero variance -an MQRankSum -an MQ 
# add https://doi.org/10.1038/nature15393 Axiom (1000G)
rule gatk_calculate_vqsr_indels:
    input:
        'jvars/cohort_{pipeline}.vcf.gz'
    output:
        recal='jvars/cohort_indels_recal_{pipeline}.vcf',
        formats='jvars/cohort_indels_recal_format_plots_{pipeline}.R',
        tranches='jvars/cohort_indels_recal_tranches_plots_{pipeline}.txt'
    log:
        'logs/gatk_calculate_vqsr_indels_{pipeline}.log'
    shell:
        """
        date > {log}
        gatk VariantRecalibrator \
        -R {hg38} \
        -V {input} \
        --resource:mills,known=false,training=true,truth=true,prior=12.0 ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
        --resource:axiom,known=false,training=true,truth=true,prior=10 ref/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ref/Homo_sapiens_assembly38.dbsnp138.vcf \
        -an QD -an FS -an SOR -mode INDEL \
        --rscript-file {output.formats} \
        --tranches-file {output.tranches} \
        -O {output.recal} 2>> {log} 
        """

rule gatk_apply_vqsr_snps:
    input:
        vcf='jvars/cohort_{pipeline}.vcf.gz',
        recal=rules.gatk_calculate_vqsr_snps.output.recal,
        tranches=rules.gatk_calculate_vqsr_snps.output.tranches 
    output:
        'jvars/vqsr_snps_{pipeline}.vcf.gz'
    log:
        'logs/gatk_apply_vqsr_snps_{pipeline}.log'
    shell:
        """
        date > {log}
        gatk ApplyVQSR \
        -R {hg38} \
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
        'jvars/vqsr_snps_indels_{pipeline}.vcf.gz'
    log:
        'logs/gatk_apply_vqsr_indels_{pipeline}.log'
    shell:
        """
        date > {log}
        gatk ApplyVQSR \
        -R {hg38} \
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
        'jvars/vqsr_snps_indels_{pipeline}.vcf.gz.tbi'
    log:
        'logs/index_vqsr_bams_{pipeline}.log'
    shell:
        """
        date > {log}
        tabix {input} -f > {output} 2>> {log}
        """


# have to skip using these because of differences in size of chr15 between vcf and 'truth' set
rule gatk_genotype_posteriors:
    input:
        vcf=rules.gatk_apply_vqsr_indels.output,
        index=rules.index_vqsr_vcfs.output
    output:
        'jvars/geno_post_{pipeline}.vcf.gz'
    log:
        'logs/gatk_geno_post_{pipeline}.log'
    shell:
        """
        date > {log}
        gatk CalculateGenotypePosteriors -V {input.vcf} \
        -supporting ref/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf \
        -O {output} 2>> {log}
        """

# with how VQSR now works, not recommended to hard filter out any variants anymore

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
        vcf=rules.gatk_apply_vqsr_indels.output,
        index=rules.index_vqsr_vcfs.output 
    output:
        'jvars/filtered_vars_{pipeline}.vcf.gz'
    log:
        'logs/gatk_filt_vars_{pipeline}.log'
    shell:
        """
        date > {log}
        gatk VariantFiltration -V {input.vcf} \
        --filter-expression "GQ < 20.0" --filter-name "GQ20" \
        --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
        -O {output} 2>> {log}
        """


rule gatk_var_call_metrics:
    input:
        rules.gatk_filter_vars.output 
    output:
        'jvars/cohort_{pipeline}.variant_calling_detail_metrics'
    log:
        'logs/gatk_var_call_metrics_{pipeline}.log'
    params:
        'jvars/cohort_{pipeline}'
    shell:
        """
        date > {log}
        gatk CollectVariantCallingMetrics -I {input} \
        --DBSNP ref/Homo_sapiens_assembly38.dbsnp138.vcf \
        -SD ref/Homo_sapiens_assembly38.dict \
        -O {params} 2>> {log}
        """

# bcftools renannotate to remove 'chr' from contigs
#for i in {1..22} X Y M;do echo "chr${i} ${i}";done > rename_chrom.txt
#bcftools annotate jvars/cohort_gdb.vcf.gz --rename-chrs rename_chrom.txt -O z -o data_renamed.vcf.gz

# try -tab output too
rule vep:
    input:
        rules.gatk_filter_vars.output 
    output:
        vcf='annotated/full_anno_vep_{pipeline}.vcf',
        stats='annotated/full_anno_vep_stats_{pipeline}.html'
    log:
        'logs/vep_{pipeline}.log'
    shell:
        """
        date > {log}
        vep -i {input} \
        --format vcf \
        --assembly GRCh38 \
        --fasta vep_dir/homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
        --offline \
        --cache \
        --dir_cache vep_dir \
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
        --plugin CADD,ref/gnomad.genomes.r3.0.snv.tsv.gz \
        --plugin PostGAP,ref/postgap_GRCh38.txt.gz,ALL \
        --dir_plugins /home/k2142172/.vep/Plugins \
        --verbose \
        --vcf \
        --stats_file {output.stats} \
        -o {output.vcf} 2>> {log}
        """

rule vep_tab:
    input:
        rules.gatk_filter_vars.output 
    output:
        vcf='annotated/full_anno_vep_{pipeline}.tab',
    log:
        'logs/vep_tab_{pipeline}.log'
    shell:
        """
        date > {log}
        vep -i {input} \
        --format vcf \
        --assembly GRCh38 \
        --fasta vep_dir/homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
        --offline \
        --cache \
        --dir_cache vep_dir \
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
        --plugin CADD,ref/gnomad.genomes.r3.0.snv.tsv.gz \
        --plugin PostGAP,ref/postgap_GRCh38.txt.gz,ALL \
        --dir_plugins /home/k2142172/.vep/Plugins \
        --verbose \
        --tab \
        -o {output.vcf} 2>> {log}
        """

rule index_vep:
    input:
        rules.vep.output.vcf 
    output:
        vcf='annotated/full_anno_vep_{pipeline}.vcf.gz',
        csi='annotated/full_anno_vep_{pipeline}.vcf.gz.csi'
    log:
        'logs/index_vep_{pipeline}.log'
    shell:
        """
        date > {log}
        bgzip {input}
        bcftools index {output.vcf}
        """

rule filter_vep:
    input:
        rules.index_vep.output.vcf 
    output:
        'annotated/impact_filtered_vep_{pipeline}.vcf'
    log:
        'logs/filter_vep_{pipeline}.log'
    shell:
        """
        date > {log}
        bcftools view -f 'PASS' {input} | filter_vep --format vcf \
        --filter '(IMPACT = MODERATE) or (IMPACT = HIGH)' --only_matched -o {output} 2>> {log}
        """
        
rule csqs_and_headers:
    input:
        rules.filter_vep.output
    output:
        csq='annotated/impact_filtered_vep_csq_{pipeline}.txt',
        headers='annotated/impact_filtered_vep_anno_headers_{pipeline}.txt'
    log:
        'logs/csqs_{pipeline}.log'
    shell:
        """
        date > {log}
        bcftools view -h {input} | grep CSQ > {output.csq} 2>> {log}
        bcftools view -h {input} | tail -n 1 > {output.headers} 2>> {log}
        """

rule vcf_to_tsv:
    input:
        vcf=rules.filter_vep.output,
        csq=rules.csqs_and_headers.output.csq,
        headers=rules.csqs_and_headers.output.headers 
    output:
        'annotated/impact_filtered_vep_{pipeline}.vcf.processed.RData',
        'annotated/impact_filtered_vep_{pipeline}.vcf.processed.tsv'
    log:
        'logs/vcf_to_tsv_{pipeline}.log'
    shell:
        """
        date > {log}
        Rscript vcf_to_tsv.R \
        {input.vcf} {input.headers} {input.csq}
        """

# bgzip annotated/full_anno_vep.vcf
# bgzip annotated/full_anno_vep_gdb.vcf
# bcftools index annotated/impact_filtered_vep.vcf.gz
# bcftools index annotated/impact_filtered_vep_gdb.vcf.gz
# bcftools index --stats annotated/impact_filtered_vep.vcf.gz > annotated/impact_filtered_vep_stats.txt
# bcftools index --stats annotated/impact_filtered_vep_gdb.vcf.gz > annotated/impact_filtered_vep_stats_gdb.txt

# validate variants
# gatk SelectVariants --exlcude-filtered # to remove any variants that failed any site-level filter
