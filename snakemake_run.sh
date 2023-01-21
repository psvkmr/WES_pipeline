snakemake -np fastqcs/SRR222689{42..51}_{1,2}_fastqc.zip
snakemake -np multiqcs/multiqc_fastqc/multiqc_report.html 
snakemake -np bams/SRR222689{42..51}.sam 
snakemake -np bams/SRR222689{42..51}_sorted.bam
snakemake -np bams/SRR222689{42..51}_fixed.bam
snakemake -np bams/SRR222689{42..51}_fixed_bam_indexstats.txt
snakemake -np mpileups/SRR222689{42..51}_mpileup_subset.tab
#snakemake -np multiqcs/multiqc_samtools/multiqc_report.html
snakemake -np gvars/SRR222689{42..51}.g.vcf.gz --cores 2
snakemake -np gdb_done.txt --cores 4
snakemake -np jvars/cohort_{cmb,gdb}.g.vcf.gz --cores 4
snakemake -np jvars/cohort_{cmb,gdb}.vcf.gz 
snakemake -np jvars/vqsr_snps_indels_{cmb,gdb}.vcf.gz.tbi
snakemake -np jvars/cohort_{cmb,gdb}.variant_calling_detail_metrics
snakemake -np annotated/full_anno_vep_{cmb,gdb}.vcf
snakemake -np annotated/impact_filtered_vep_{cmb,gdb}.vcf
snakemake -np annotated/impact_filtered_vep_{cmb,gdb}.vcf.processed.tsv