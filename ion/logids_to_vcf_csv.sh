#!/bin/bash

args=("$@")
SET=$args{[0]}
WD=${args[1]}/${args[0]}
LOGIDS=${args[2]}

mkdir -p ${WD}
scp ${LOGIDS} ${WD}/${SET}_logIds.txt

CONVERT_TRACKED=/hades/psivakumar/cohorts/tracked_html_to_csv.py
TRACKED=/hades/psivakumar/cohorts/tracked.csv
REF=/hades/psivakumar/pipeline/grch38/Homo_sapiens_assembly38.fasta
REF_IND=/hades/psivakumar/pipeline/grch38/Homo_sapiens_assembly38.index
TMP_SPACE=/hades/psivakumar/pipeline/tmp
KIT_1=/hades/psivakumar/pipeline/SureSelect_Human_All_Exon_V6_UTR_r2/S07604624_Regions.intervallist
KIT_2=/hades/psivakumar/pipeline/SureSelect_Human_All_Exon_V6_UTR_r2/S07604624_Regions.bed
KOIOS=/mnt/qsg-results/inprogress/Koios_release_4_all.vcf.gz
TO_CSV=/hades/psivakumar/pipeline/vcf_to_csv.R

GATK=/hades/Software/NGS_Software/gatk-4.1.4.0/gatk
BCFTOOLS=/hades/Software/NGS_Software/bcftools-1.9/bcftools/bcftools
VEP=/hades/psivakumar/pipeline/runVEP_hades.sh
TABIX=/home/hades/anaconda3/bin/tabix


# update tracked csv from tracked html
python3.5 ${CONVERT_TRACKED}

# check if all samples in tracked
logid_in_tracked(){
  l1=$(grep -F -f ${WD}/${SET}_logIds.txt ${TRACKED} | wc -l)
  l2=$(wc -l ${WD}/${SET}_logIds.txt | cut -d' ' -f1)
  cat ${TRACKED} | grep -F -w -f ${WD}/${SET}_logIds.txt > ${WD}/${SET}_samples_in_tracked.txt
  if [ $l1 != $l2 ]
    then
      echo "some log ids missing from tracked file"
    else
      echo "all log ids in tracked file"
  fi
}

# get gvcf locs
gvcf_in_tracked(){
  cat ${TRACKED} | grep -F -w -f ${WD}/${SET}_logIds.txt | cut -d',' -f2,7 > ${WD}/${SET}_gvcf_locs.csv
  awk '$1=$1' FS="," OFS="\t" ${WD}/${SET}_gvcf_locs.csv > ${WD}/${SET}_gvcf_locs.txt
  #sed -i 's/,/     /g' ${WD}/${SET}_gvcf_locs.txt
  l1=$(grep -F 'g.vcf.gz' ${WD}/${SET}_gvcf_locs.txt | wc -l)
  l2=$(cat ${WD}/${SET}_samples_in_tracked.txt | wc -l)
  if [ $l1 != $l2 ]
    then
      echo "some tracked file log ids have no associated gvcf"
    else
      echo "all tracked file log ids have associated gvcfs"
  fi
}

# merge gvcfs
merge_gvcfs(){
  ${GATK} GenomicsDBImport --genomicsdb-workspace-path ${WD}/gdbworkspace -L ${KIT_2} --sample-name-map ${WD}/${SET}_gvcf_locs.txt --tmp-dir ${TMP_SPACE} --batch-size 0 --create-output-variant-index
}

#--java-options "-Xmx4g -Xms4g"
# alt merge
combine_gvcfs(){
  sed 's/LI.*\/mnt/--variant \/mnt/g' ${WD}/${SET}_gvcf_locs.txt > ${WD}/${SET}_cbgvcf.sh
  sed -i 's/$/ \\/g' ${WD}/${SET}_cbgvcf.sh
  FUN_1="${GATK} CombineGVCFs -R ${REF} --create-output-variant-index"
  sed -i '1 i\'"${FUN_1}"'' ${WD}/${SET}_cbgvcf.sh
  sed -i '1{s/$/ \\/}' ${WD}/${SET}_cbgvcf.sh
  FUN_2="-O ${WD}/${SET}_cbgvcf.vcf.gz"
  sed -i "\$a${FUN_2}" ${WD}/${SET}_cbgvcf.sh
  bash ${WD}/${SET}_cbgvcf.sh
}

# gvcfs to vcf
gvcf_to_vcf(){
  ${GATK} GenotypeGVCFs -R ${REF} -V ${WD}/${SET}_cbgvcf.vcf.gz --allow-old-rms-mapping-quality-annotation-data --annotate-with-num-discovered-alleles --create-output-variant-index -O ${WD}/${SET}.j.vcf.gz
}

# only keep polymorphic variants
get_polymorphics(){
  ${BCFTOOLS} view -c1.0:minor ${WD}/${SET}.j.vcf.gz -O z -o ${WD}/${SET}_poly.vcf.gz
}

# annotate
annotate(){
  bash ${VEP} ${WD}/${SET}.j.vcf.gz
}

# add QC steps

# create list of variants for which to get total het and hom counts
# check tab space matches vcf tab delim, chage to query command
get_variants_list(){
  ${BCFTOOLS} query -f '%CHROM\t%POS\t%REF\t%ALT\n' ${WD}/${SET}.j.vcf.gz > ${WD}/${SET}_variantIds.txt
}

# create info files
info_files(){
  ${BCFTOOLS} view -h ${WD}/${SET}.j.vcf.gz.vep.out | grep CSQ > ${WD}/${SET}.csq
  ${BCFTOOLS} view -h ${WD}/${SET}.j.vcf.gz.vep.out | tail -n 1 > ${WD}/${SET}.headers
}

vcf_to_csv(){
  scp ${WD}/${SET}.j.vcf.gz.vep.out ${WD}/${SET}_vep.j.vcf.gz
  bgzip -d ${WD}/${SET}_vep.j.vcf.gz
  Rscript --vanilla ${TO_CSV} ${WD}/${SET}_vep.j.vcf ${WD}/${SET}.headers ${WD}/${SET}.csq
}

echo "checking for logids in tracked file"
logid_in_tracked 2> ${WD}/${SET}.log
wait
if ! [ -s ${WD}/${SET}_samples_in_tracked.txt ]; then
  echo "no valid logIds" > ${WD}/${SET}.log
  exit
fi
echo "checking for gvcfs in tracked"
gvcf_in_tracked 2>> ${WD}/${SET}.log
wait
echo "combining gvcfs"
combine_gvcfs 2>> ${WD}/${SET}.log
wait
echo "gvcf to vcf"
gvcf_to_vcf 2>> ${WD}/${SET}.log
wait
echo "annotating with VEP"
annotate 2>> ${WD}/${SET}.log
wait
echo "making info files"
info_files 2>> ${WD}/${SET}.log
wait
echo "converting to csv format"
vcf_to_csv 2>> ${WD}/${SET}.log
wait
echo "finished" > ${WD}/${SET}.log
pwd
