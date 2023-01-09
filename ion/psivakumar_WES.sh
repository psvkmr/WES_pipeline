#pipeline

# required: WES fastqs 1 and 2; ref sequence fasta; BWA; GATK; samtools; verifybamid;

# what kit, machine used for sequencing
# default exome intervals list generated with:
# curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz" | gunzip -c | cut -f 3,5,6 | sort -t $'\t' -k1,1 -k2,2n | bedtools merge -i - > exome.bed

args=("$@")
FQ1=${args[0]}
FQ2=${args[1]}
WD=${args[2]}
NAME=${args[3]}
mkdir -p ${WD}

#FQ1=/hades/psivakumar/pipeline/test_files/Swapnil_B_Tatte_R1.fastq.gz
#FQ2=/hades/psivakumar/pipeline/test_files/Swapnil_B_Tatte_R2.fastq.gz
#NAME=/hades/psivakumar/pipeline/test_files/test_Swap
#OUT_BAM1=/hades/psivakumar/pipeline/test_files/test_Swap.bam
#OUT_BAM2=/hades/psivakumar/pipeline/test_files/Swap.bam
BWA=/hades/Software/NGS_Software/bwa/bwa
SAMTOOLS=/hades/Software/NGS_Software/samtools-1.9/samtools
REF=/hades/psivakumar/pipeline/grch38/Homo_sapiens_assembly38.fasta
REF_IND=/hades/psivakumar/pipeline/grch38/Homo_sapiens_assembly38.index
GATK=/hades/Software/NGS_Software/gatk-4.1.4.0/gatk
#GATK=/hades/Software/NGS_Software/gatk-4.0.0.0/gatk-4.0.0.0/gatk
VBI=/hades/Software/NGS_Software/verifyBamID
VEP=/hades/psivakumar/pipeline/runVEP_hades.sh
TABIX=/home/hades/anaconda3/bin/tabix
NOVOALIGN=/hades/Software/NGS_Software/novocraftV3.08.02/novoalign
NOVOSORT=/hades/Software/NGS_Software/novocraftV3.08.02/novosort

INDELS=/hades/psivakumar/pipeline/gatk_hg38_resource_bundle/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
SNPS=/hades/psivakumar/pipeline/gatk_hg38_resource_bundle/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz
DBSNPS=/hades/psivakumar/pipeline/gatk_hg38_resource_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf
KIT_1=/hades/psivakumar/pipeline/SureSelect_Human_All_Exon_V6_UTR_r2/S07604624_Regions.intervallist
KIT_2=/hades/psivakumar/pipeline/SureSelect_Human_All_Exon_V6_UTR_r2/S07604624_Regions.bed
#KIT=/hades/psivakumar/pipeline/liftedover_SureSelect_V6_UTR.bed
#${GATK} BedToIntervalList -I /hades/psivakumar/pipeline/SureSelect_Human_All_Exon_V6_UTR_r2/S07604624_Regions.bed -SD ${REF} -O /hades/psivakumar/pipeline/SureSelect_Human_All_Exon_V6_UTR_r2/S07604624_Regions.intervallist
#--java-options "-Xmx4g"

# bwa -t 320
#${BWA} index ${REF}
#${SAMTOOLS} faidx ${REF}
#${GATK} CreateSequenceDictionary -R ${REF} -O ${REF}.dict
#${GATK} IndexFeatureFile -F ${INDELS}
#${GATK} IndexFeatureFile -F ${SNPS}
#${GATK} IndexFeatureFile -F ${DBSNPS}
#${BWA} mem -t 24 -T 0 -M -R "@RG\tID:1\tPU:Unit\tLB:library\tSM:${NAME}\tPL:illumina" ${REF} ${FQ1} ${FQ2} > ${WD}/${NAME}.aln.sam #add log?
${BWA} mem -t 24 -T 0 -M -H -a -P -R "@RG\tID:1\tPU:Unit\tLB:library\tSM:${NAME}\tPL:illumina" ${REF} ${FQ1} ${FQ2} > ${WD}/${NAME}.aln.sam 2> ${WD}/${NAME}.log
${GATK} SortSam -I ${WD}/${NAME}.aln.sam -SO coordinate -O ${WD}/${NAME}.tmp.sorted.bam 2>> ${WD}/${NAME}.log
${GATK} MarkDuplicates -I ${WD}/${NAME}.tmp.sorted.bam -M ${WD}/${NAME}.markedDups.txt -O ${WD}/${NAME}.tmp.markDups.sorted.bam 2>> ${WD}/${NAME}.log
${GATK} BuildBamIndex -I ${WD}/${NAME}.tmp.markDups.sorted.bam 2>> ${WD}/${NAME}.log
#${SAMTOOLS} view -Shb ${WD}/${NAME}.aln.sam -o ${OUT_BAM1}.unsorted
#${SAMTOOLS} sort ${OUT_BAM1}.unsorted -o ${OUT_BAM1}
#${SAMTOOLS} index ${OUT_BAM1}

# GATK

#${GATK} RealignerTargetCreator -R ${REF} --known ${INDELS} -I ${OUT_BAM1} -o ${WD}/${NAME}.realign.intervals
#${GATK} IndelRealigner -R ${REF} --known ${INDELS} -targetIntervals ${WD}/${NAME}.realign.intervals --noOriginalAlignmentTags -I ${OUT_BAM1} -nWayOut ${WD}/${NAME}.map
${GATK} BaseRecalibrator --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' -R ${REF} -I ${WD}/${NAME}.tmp.markDups.sorted.bam --known-sites ${SNPS} --known-sites ${DBSNPS} --known-sites ${INDELS} -O ${WD}/${NAME}.bqsr.grp 2>> ${WD}/${NAME}.log
#${GATK} BaseRecalibrator -R ${REF} -I ${OUT_BAM1} --BQSR ${WD}/${NAME}.bqsr.grp -O ${WD}/${NAME}.bqsr.grp #--plots ${WD}/${NAME}.afterRecal.pdf
${GATK} AnalyzeCovariates --bqsr ${WD}/${NAME}.bqsr.grp --plots ${WD}/${NAME}.covariates.pdf 2>> ${WD}/${NAME}.log
${GATK} ApplyBQSR  -R ${REF} -I ${WD}/${NAME}.tmp.markDups.sorted.bam --bqsr-recal-file ${WD}/${NAME}.bqsr.grp -O ${WD}/${NAME}.markDups.sorted.bam 2>> ${WD}/${NAME}.log
#${GATK} PrintReads -R ${REF} -I ${OUT_BAM1} --BQSR ${WD}/${NAME}.bqsr.grp -o ${OUT_BAM2} 2>> ${WD}/${NAME}.log
${GATK} BuildBamIndex -I ${WD}/${NAME}.markDups.sorted.bam 2>> ${WD}/${NAME}.log
#${SAMTOOLS} index ${OUT_BAM2}
${VBI} --vcf ${SNPS} --bam ${WD}/${NAME}.markDups.sorted.bam --out ${WD}/${NAME}.markDups.sorted.vbi 2>> ${WD}/${NAME}.log
${GATK} BamIndexStats -I ${WD}/${NAME}.markDups.sorted.bam > ${WD}/${NAME}.markDups.sorted.indexstats.txt 2>> ${WD}/${NAME}.log
${GATK} CollectWgsMetrics -R ${REF} -I ${WD}/${NAME}.markDups.sorted.bam --INTERVALS ${KIT_1} --INCLUDE_BQ_HISTOGRAM -O ${WD}/${NAME}.markDups.sorted.WGSmetrics.txt 2>> ${WD}/${NAME}.log

# add dbsnp file
${GATK} HaplotypeCaller -R ${REF} --intervals ${KIT_2} -I ${WD}/${NAME}.markDups.sorted.bam -ERC GVCF -O ${WD}/${NAME}.g.vcf.gz 2>> ${WD}/${NAME}.log
# {GATK} CombineGVCFs / GenomicsDBImport
${GATK} GenotypeGVCFs -R ${REF} -V ${WD}/${NAME}.g.vcf.gz --annotate-with-num-discovered-alleles --create-output-variant-index -O ${WD}/${NAME}.vcf.gz 2>> ${WD}/${NAME}.log
#${TABIX} -p vcf ${WD}/${NAME}.vcf.gz
${GATK} ValidateVariants -R ${REF} -V ${WD}/${NAME}.vcf.gz 2>> ${WD}/${NAME}.log
#${GATK} VariantAnnotator -R ${REF} -I ${OUT_BAM2} -V ${WD}/${NAME}.g.vcf.gz -L ${WD}/${NAME}.g.vcf.gz -A MQ0 -A SpanningDeletions -o ${WD}/${NAME}.reannotated.g.vcf.gz
#${GATK} CNNScoreVariants -R ${REF} -V ${WD}/${NAME}.g.vcf.gz -O ${WD}/${NAME}.CNNanno.vcf.gz
#${GATK} FilterVariantTranches
# implement VQSR
bash ${VEP} ${WD}/${NAME}.vcf.gz 2>> ${WD}/${NAME}.log

# for control purposes, novoalign method
# -k quality calibration; -H hard clipping of 3' low qual bases; -t maximum alignemtn score acceptable for best alignment; if -k on then write original base qualities as SAM OQ:Z; -hdrhd limit set to check identity between headers in pe reads;
${NOVOALIGN} -c 10 -d ${REF_IND} -f ${FQ1} ${FQ2} -o SAM $"@RG\tID:1\tPU:Unit\tLB:library\tSM:${NAME}\tPL:illumina"  --rOQ --hdrhd 3 -H 15 -k -o Soft -t 320  > ${WD}/${NAME}.novo.sam 2> ${WD}/${NAME}.log
${SAMTOOLS} view -Sb ${WD}/${NAME}.novo.sam > ${WD}/${NAME}.novo.tmp.unsorted.bam 2>> ${WD}/${NAME}.log
# --md mark duplicates; --kt keep sam tags added for duplicate detection; single end read duplicates also can be removed; -c threads; -t tmp dir; -f force coordinate sort even if already done; -i create index
${NOVOSORT} --md --kt --ise -c 10 -t /hades/psivakumar/pipeline/tmp/ -f -i -o ${WD}/${NAME}.novo.tmp.sorted.bam ${WD}/${NAME}.novo.tmp.unsorted.bam 2>> ${WD}/${NAME}.log
${GATK} BaseRecalibrator --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' -R ${REF} -I ${WD}/${NAME}.novo.tmp.sorted.bam --known-sites ${SNPS} --known-sites ${INDELS} --known-sites ${DBSNPS} -O ${WD}/${NAME}.novo.bqsr.grp 2>> ${WD}/${NAME}.log
#${GATK} BaseRecalibrator -R ${REF} -I ${OUT_BAM1} --BQSR ${NAME}.bqsr.grp -O ${NAME}.bqsr.grp #--plots ${NAME}.afterRecal.pdf
${GATK} AnalyzeCovariates --bqsr ${WD}/${NAME}.novo.bqsr.grp --plots ${WD}/${NAME}.novo.covariates.pdf 2>> ${WD}/${NAME}.log
${GATK} ApplyBQSR  -R ${REF} -I ${WD}/${NAME}.novo.tmp.sorted.bam --bqsr-recal-file ${WD}/${NAME}.novo.bqsr.grp -O ${WD}/${NAME}.novo.sorted.bam 2>> ${WD}/${NAME}.log
#${GATK} PrintReads -R ${REF} -I ${OUT_BAM1} --BQSR ${NAME}.bqsr.grp -o ${OUT_BAM2}
${GATK} BuildBamIndex -I ${WD}/${NAME}.novo.sorted.bam 2>> ${WD}/${NAME}.log
#${SAMTOOLS} index ${OUT_BAM2}
${VBI} --vcf ${SNPS} --bam ${WD}/${NAME}.novo.sorted.bam --out ${WD}/${NAME}.novo.sorted.vbi 2>> ${WD}/${NAME}.log
${GATK} BamIndexStats -I ${WD}/${NAME}.novo.sorted.bam > ${WD}/${NAME}.novo.sorted.bam.indexstats 2>> ${WD}/${NAME}.log
${GATK} CollectWgsMetrics -R ${REF} -I ${WD}/${NAME}.novo.sorted.bam --INTERVALS ${KIT_1} --INCLUDE_BQ_HISTOGRAM -O ${WD}/${NAME}.novo.sorted.bam.WGSmetrics 2>> ${WD}/${NAME}.log

# add dbsnp file
${GATK} HaplotypeCaller -R ${REF} --intervals ${KIT_2} -I ${WD}/${NAME}.novo.sorted.bam -ERC GVCF -O ${WD}/${NAME}.novo.g.vcf.gz 2>> ${WD}/${NAME}.log
# {GATK} CombineGVCFs / GenomicsDBImport
${GATK} GenotypeGVCFs -R ${REF} -V ${WD}/${NAME}.novo.g.vcf.gz --annotate-with-num-discovered-alleles --create-output-variant-index -O ${WD}/${NAME}.novo.vcf.gz 2>> ${WD}/${NAME}.log
#${TABIX} -p vcf ${WD}/${NAME}.vcf.gz
${GATK} ValidateVariants -R ${REF} -V ${WD}/${NAME}.novo.vcf.gz 2>> ${WD}/${NAME}.log
bash ${VEP} ${WD}/${NAME}.novo.vcf.gz 2>> ${WD}/${NAME}.log
