
conda create -n gatk_full -c bioconda gatk4 sra-tools fastqc multiqc fastx_toolkit trim-galore bwa samtools snpeff bcftools qualimap

# conda sra tools v3 causes too many conflicts, download outside of package manager
# download sra-toolkit
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.2/sratoolkit.3.0.2-ubuntu64.tar.gz

# bcftools not working as missing a lib
#sudo apt-get install libopenblas-dev

# example data
# https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP407611&o=acc_s%3Aa
wget https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP407611&o=acc_s%3Aa
# manually download SRR Acc list and Sra Run Table 

# subset samples as not enough space
head -n 10 SRR_Acc_List.txt > subset_SRR_Acc_List.txt
# fetch data from samples list
sratoolkit.3.0.2-ubuntu64/bin/prefetch --option-file subset_SRR_Acc_List.txt

# sras variable
sras=$(ls */*.sra)
# use fastq dump to get paired end fq1 and fq2 per sample
for sra in $sras
do
    echo $sra
    sratoolkit.3.0.2-ubuntu64/bin/fasterq-dump --split-files -O ./fastqs/ $sra
done

mkdir sras
mv SRR2* sras/.
#sras=$(ls sras/*/*.sra)

# run fastqc, multiqc
mkdir fastqcs multiqcs
fastqc fastqs/* -o fastqcs/
multiqc fastqcs/* -o multiqcs/ -n multiqc_fastqc.html


# download ref data if don't have
mkdir ref
wget -v https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget -v https://ftp.ensembl.org/pub/release-108/gff3/homo_sapiens/Homo_sapiens.GRCh38.108.gff3.gz
# samtools mpileup needs bgzip format
gunzip -c Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz | bgzip -c > Homo_sapiens.GRCh38.dna.primary_assembly.fa.bgz
# index ref data
bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
# index gff3 ?
# samtools indexing for bgz not working
#samtools index Homo_sapiens.GRCh38.dna.primary_assembly.fa.bgz

mkdir bams
# align data
for item in $(ls sras)
do
    echo "processing ${item}"
    bwa mem -M -R "@RG\tID:${item}\tSM:${item}\tLB:WES\tPL:Illumina" ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz  fastqs/${item}_1.fastq fastqs/${item}_2.fastq  > bams/${item}.sam
    samtools view -S -b bams/${item}.sam > bams/${item}.bam
    samtools sort -o bams/${item}_sorted.bam bams/${item}.bam
    samtools index -c bams/${item}_sorted.bam
    samtools flagstat bams/${item}_sorted.bam > bams/${item}_sorted_flagstats.txt
done

# zip fastqs to save space
for item in $(ls sras)
do 
    echo "zipping fastqs for ${item}"
    gzip fastqs/${item}_1.fastq
    gzip fastqs/${item}_2.fastq
done

#head -n 100000 fastqs/SRR22268942_1.fastq > tfq1.fastq
#head -n 100000 fastqs/SRR22268942_2.fastq > tfq2.fastq
#bwa mem -M -R "@RG\tID:${item}\tSM:${item}\tLB:WES\tPL:Illumina" ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz  tfq1.fastq tfq2.fastq  > bams/${item}.sam

# -a for secondary alignments from single unpaired reads
# -M mark shorter split hits as secondary
# -I specify stat parameters for insert size distribution

# removed .sam files for space

for item in $(ls sras)
do
    echo "processing bam for ${item}"
    gatk MarkDuplicates -I bams/${item}_sorted.bam -O bams/${item}_duped.bam -M bams/${item}_duped_metrics.txt 
    gatk FixMateInformation -SO coordinate -I bams/${item}_duped.bam -O bams/${item}_fixed.bam
    gatk BuildBamIndex -I bams/${item}_fixed.bam
done

#wget https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
#wget https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
#wget https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
# move to ref/ folder and rename as .gz for two files
mv ref/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf ref/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz
mv ref/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf ref/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
bgzip ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf
bgzip ref/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf

# hacky way to remove 'chr' from CHROM column values, problem is it leaves it in ##contig so they no longer match
bgzip -cd ref/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz | sed 's/chr//g' | bgzip > ref/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38_nochr.vcf.gz
bgzip -cd ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz | sed 's/chr//g' | bgzip > ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138_nochr.vcf.gz
bgzip -cd ref/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz | sed 's/chr//g' | bgzip > ref/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38_nochr.vcf.gz

# index files myself as .tbi
bcftools index -t ref/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38_nochr.vcf.gz
bcftools index -t ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138_nochr.vcf.gz
bcftools index -t ref/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38_nochr.vcf.gz

# gatk may not work with gzipped data so use unzipped
gunzip -c ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa
gatk CreateSequenceDictionary -R ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa
samtools faidx ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa

for item in $(ls sras)
do
    echo "gatk bqsr processing for ${item}"
    gatk BaseRecalibrator -R ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa -I bams/${item}_fixed.bam --known-sites ref/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38_nochr.vcf.gz --known-sites ref/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38_nochr.vcf.gz --known-sites ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138_nochr.vcf.gz -O bams/${item}_qual_recal_table.txt
done 

# need to install ggplot to use
#gatk AnalyzeCovariates --bqsr bams/${item}_qual_recal_table.txt --plots bams/${item}_covariates.pdf

# read group issues: java.lang.IllegalStateException: The covariates table is missing ReadGroup SRR22268942 in RecalTable0
#gatk ApplyBQSR -R ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa -I bams/SRR22268942_fixed.bam --bqsr bams/${item}_qual_recal_table.txt -O bams/${item}_bqsr.bam 
#gatk ValidateSamFile -I bams/${item}_bqsr.bam >> bams/check_bam_log.txt

# for now, using fixed bam as bqsr bam
for item in $(ls sras)
do
    cp bams/${item}_fixed.bam bams/${item}_bqsr.bams
done


# check bam
# https://davetang.org/muse/2015/08/26/samtools-mpileup/
mkdir mpileups
for item in $(ls sras)
do
    samtools mpileup -f ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa -s bams/${item}_bqsr.bam | sed -n '500,515p' > mpileups/${item}_mpileup_subset.tab
done

# get index and stats
for item in $(ls sras)
do
    echo "index stats for ${item}"
    gatk BuildBamIndex -I bams/${item}_bqsr.bam -O bams/${item}_bqsr.bam.bai
    gatk BamIndexStats -I bams/${item}_bqsr.bam > bams/${item}_bqsr_indexstats.txt
done 

# don't have kit details, they used customised panel
#gatk CollectWgsMetrics -R ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa -I bams/${item}_bqsr.bam --INTERVALS ${KIT_1} --INCLUDE_BQ_HISTOGRAM -O bams/${item}_bqsr_WGSmetrics.txt

mkdir variants
# haplotype caller
# add intervals, interval padding 
for item in $(ls sras)
do
    echo "haplotype calling for ${item}"
    gatk HaplotypeCaller -ERC GVCF -R ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa -I bams/${item}_bqsr.bam -O variants/${item}.g.vcf.gz
done

# subset SRR_Acc_List.txt for used IDs
# get paths for each of those files
while read line; 
do 
    readlink -f variants/${line}.g.vcf.gz; 
done < sample_ids.txt > sample_vcf_paths.txt

# vcf mapping file
paste sample_ids.txt sample_vcf_paths.txt > sample_vcf_mapping.txt

# form sample map with 2 cols tab separated: sample name \t path to sample gvcf
gatk GenomicsDBImport -L 1 -L 2 -L 3 -L 4 -L 5 -L 6 -L 7 -L 8 -L 9 -L 10 -L 11 -L 12 -L 13 -L 14 -L 15 -L 16 -L 17 -L 18 -L 19 -L 20 -L 21 -L 22 -L X -L Y --sample-name-map sample_vcf_mapping.txt --genomicsdb-workspace-path genomicsdb
gatk CombineGVCFs -R ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa -O variants/cohort.g.vcf.gz \
--variant /mnt/c/Users/Prasanth/Documents/wes/variants/SRR22268942.g.vcf.gz \
--variant /mnt/c/Users/Prasanth/Documents/wes/variants/SRR22268943.g.vcf.gz \
--variant /mnt/c/Users/Prasanth/Documents/wes/variants/SRR22268944.g.vcf.gz \
--variant /mnt/c/Users/Prasanth/Documents/wes/variants/SRR22268945.g.vcf.gz \
--variant /mnt/c/Users/Prasanth/Documents/wes/variants/SRR22268946.g.vcf.gz \
--variant /mnt/c/Users/Prasanth/Documents/wes/variants/SRR22268947.g.vcf.gz \
--variant /mnt/c/Users/Prasanth/Documents/wes/variants/SRR22268948.g.vcf.gz \
--variant /mnt/c/Users/Prasanth/Documents/wes/variants/SRR22268949.g.vcf.gz \
--variant /mnt/c/Users/Prasanth/Documents/wes/variants/SRR22268950.g.vcf.gz \
--variant /mnt/c/Users/Prasanth/Documents/wes/variants/SRR22268951.g.vcf.gz

# get cohort gvcf 

gatk GenotypeGVCFs -R ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa -V gendb://genomicsdb --annotate-with-num-discovered-alleles --create-output-variant-index -O variants/cohort_gdb.vcf.gz 
gatk GenotypeGVCFs -R ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa -V variants/cohort.g.vcf.gz --annotate-with-num-discovered-alleles --create-output-variant-index -O variants/cohort.vcf.gz 

gatk ValidateVariants -R ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa -V variants/cohort.vcf.gz 
gatk ValidateVariants -R ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa -V variants/cohort_gdb.vcf.gz 

# variant vqsr - assign new confidence score for each vavriant: VQSLOD
# require more ref files
wget https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz
wget https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz
# move to ref
# same hacky edit to remove 'chr'
bgzip -cd resources_broad_hg38_v0_hapmap_3.3.hg38.vcf | sed 's/chr//g' | bgzip > resources_broad_hg38_v0_hapmap_3.3.hg38_nochr.vcf.gz
bgzip -cd resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf | sed 's/chr//g' | bgzip > resources_broad_hg38_v0_1000G_omni2.5.hg38_nochr.vcf.gz
bcftools index -t resources_broad_hg38_v0_hapmap_3.3.hg38_nochr.vcf.gz
bcftools index -t resources_broad_hg38_v0_1000G_omni2.5.hg38_nochr.vcf.gz

# check snps
# rscript requires ggplot2
#conda install r-ggplot2
gatk VariantRecalibrator -R ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa -V variants/cohort.vcf.gz --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ref/resources_broad_hg38_v0_hapmap_3.3.hg38_nochr.vcf.gz --resource:omini,known=false,training=true,truth=false,prior=12.0 ref/resources_broad_hg38_v0_1000G_omni2.5.hg38_nochr.vcf.gz --resource:1000G,known=false,training=true,truth=false,prior=10.0 ref/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38_nochr.vcf.gz --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138_nochr.vcf.gz -an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum -mode SNP -tranche 100.0 -tranche 99.0 -tranche 99.0 -tranche 95.0 -tranche 90.0 --rscript-file variants/output.plots.R --tranches-file variants/output.tranches -O variants/cohort_snps_recal.txt 
gatk VariantRecalibrator -R ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa -V variants/cohort_gdb.vcf.gz --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ref/resources_broad_hg38_v0_hapmap_3.3.hg38_nochr.vcf.gz --resource:omini,known=false,training=true,truth=false,prior=12.0 ref/resources_broad_hg38_v0_1000G_omni2.5.hg38_nochr.vcf.gz --resource:1000G,known=false,training=true,truth=false,prior=10.0 ref/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38_nochr.vcf.gz --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138_nochr.vcf.gz -an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum -mode SNP -tranche 100.0 -tranche 99.0 -tranche 99.0 -tranche 95.0 -tranche 90.0 --rscript-file variants/output_gdb.plots.R --tranches-file variants/output_gdb.tranches -O variants/cohort_gdb_snps_recal.txt 
# check indels
gatk VariantRecalibrator -R ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa -V variants/cohort.vcf.gz --resource:mills,known=false,training=true,truth=true,prior=15.0 ref/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38_nochr.vcf.gz -an QD -an MQ -an MQRankSum -an FS -an SOR -mode INDEL --rscript-file variants/output.plots.indels.R --tranches-file variants/output_indels.tranches -O variants/cohort_indels_recal.txt
gatk VariantRecalibrator -R ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa -V variants/cohort_gdb.vcf.gz --resource:mills,known=false,training=true,truth=true,prior=15.0 ref/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38_nochr.vcf.gz -an QD -an MQ -an MQRankSum -an FS -an SOR -mode INDEL --rscript-file variants/output.plots.indels_gdb.R --tranches-file variants/output_indels_gdb.tranches -O variants/cohort_gdb_indels_recal.txt

# apply vqsr
gatk ApplyVQSR -R ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa -V variants/cohort.vcf.gz --truth-sensitivity-filter-level 99.0 --tranches-file variants/output.tranches --recal-file variants/cohort_snps_recal.txt -mode SNP -O variants/vqsr_snps.vcf.gz
gatk ApplyVQSR -R ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa -V variants/vqsr_snps.vcf.gz --truth-sensitivity-filter-level 99.0 --tranches-file variants/output_indels.tranches --recal-file variants/cohort_indels_recal.txt -mode INDEL -O variants/vqsr_snps_indels.vcf.gz
# and to gdb
gatk ApplyVQSR -R ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa -V variants/cohort_gdb.vcf.gz --truth-sensitivity-filter-level 99.0 --tranches-file variants/output.tranches --recal-file variants/cohort_gdb_snps_recal.txt -mode SNP -O variants/vqsr_gdb_snps.vcf.gz
gatk ApplyVQSR -R ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa -V variants/vqsr_gdb_snps.vcf.gz --truth-sensitivity-filter-level 99.0 --tranches-file variants/output_indels_gdb.tranches --recal-file variants/cohort_gdb_indels_recal.txt -mode INDEL -O variants/vqsr_gdb_snps_indels.vcf.gz

# genotype posteriors - require ped file 

# conda vep install
conda install -c bioconda ensembl-vep
# hs caused downgrade of bcftools from 1.11 > 1.8 - has broken it
# installed bcftools=1.9, caused samtools downgrade 1.12 > 1.6
vep_install -a cf -s homo_sapiens -y GRCh38 -c /mnt/c/Users/Prasanth/Documents/wes/VEP --CONVERT

# dir VEP changed name to vep_dir
# everything adds: --sift b, --polyphen b, --ccds, --hgvs, --symbol, --numbers, --domains, --regulatory, --canonical, --protein, --biotype, --af, --af_1kg, --af_esp, --af_gnomade, --af_gnomadg, --max_af, --pubmed, --uniprot, --mane, --tsl, --appris, --variant_class, --gene_phenotype, --mirna
# add --verbose, stats_file, remove --offline, --gff, --nearest
# no -gff due to index issue --gff ref/Homo_sapiens.GRCh38.108.gff3

# small test file
bcftools filter variants/vqsr_snps_indels.vcf.gz -r 1:30000-800000 -O z -o test.vcf.gz
tabix test.vcf.gz

# calcualte genotype priors
# get 1000 genomes ref file for accuracy
#wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf
#bgzip 1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf
#tabix 1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38_nochr.vcf.gz
sed 's/chr//g' 1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf > 1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38_nochr.vcf
gatk CalculateGenotypePosteriors -V test.vcf.gz -supporting ref/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf.gz -O cgp.vcf.gz

# hard filter on low qual variants
gatk VariantFiltration -V cgp.vcf.gz --filter-expression "GQ < 20.0" --filter-name "GQ20" --filter-expression "QUAL < 30.0" --filter-name "QUAL30" -O vf.vcf.gz

# get call metrics
gatk CollectVariantCallingMetrics -I vf.vcf.gz --DBSNP ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf -SD ref/Homo_sapiens.GRCh38.dna.primary_assembly.dict -O metrics 

# wget https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/108/CADD.pm
vep_install -a p -g CADD
# get ref files for CADD
wget https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.snv.tsv.gz
wget https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.snv.tsv.gz.tbi

# install plugins for GERP and postgap data
vep_install -a p -g Conservation,PostGAP
wget https://storage.googleapis.com/postgap-data/postgap.txt.gz
gunzip postgap.txt.gz
(grep ^"ld_snp_rsID" postgap.txt; grep -v ^"ld_snp_rsID" postgap.txt | sort -k4,4 -k5,5n ) | bgzip > postgap_GRCh38.txt.gz
tabix -s 4 -b 5 -e 5 -c l postgap_GRCh38.txt.gz
# conservation failing due to no BigFile.pm

# full run
vep -i variants/vqsr_snps_indels.vcf.gz --format vcf --offline --assembly GRCh38 --vcf --fasta vep_dir/homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz --cache --dir_cache vep_dir --plugin CADD,ref/gnomad.genomes.r3.0.snv.tsv.gz --plugin PostGAP,ref/postgap_GRCh38.txt.gz,ALL --dir_plugins /home/k2142172/.vep/Plugins --check_existing --everything --per_gene --pick_order rank --verbose --stats_file full_vqsr_snps_indels_vep_stats.html --nearest gene --overlaps --gene_phenotype --show_ref_allele --total_length --clin_sig_allele 0 --af_exac --pubmed -o full_vqsr_snps_indels_vep.vcf.gz
vep -i variants/vqsr_gdb_snps_indels.vcf.gz --format vcf --offline --assembly GRCh38 --vcf --fasta vep_dir/homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz --cache --dir_cache vep_dir --plugin CADD,ref/gnomad.genomes.r3.0.snv.tsv.gz --plugin PostGAP,ref/postgap_GRCh38.txt.gz,ALL --dir_plugins /home/k2142172/.vep/Plugins --check_existing --everything --per_gene --pick_order rank --verbose --stats_file full_vqsr_gdb_snps_indels_vep_stats.html --nearest gene --overlaps --gene_phenotype --show_ref_allele --total_length --clin_sig_allele 0 --af_exac --pubmed -o full_vqsr_gdb_snps_indels_vep.vcf.gz
