#!/bin/bash

args=("$@")
SET=$1
WD=$2/$1
mkdir -p ${WD}
cd ${WD}
VCF_DIR=/mnt/qsg-results/inprogress/Release4
CHR=$(echo {1..22} X)

printf '%s\n' ${args[@]:2} > ${WD}/${SET}_dnaIds.txt

# dna number to logid
dnaid_to_logid(){
cut -d',' -f1,38- /mnt/NG/Sequencing-data/all-sequencing-log/all_Log_sequencing.csv | grep -f ${WD}/${SET}_dnaIds.txt | cut -d',' -f1,2 | sed s/^/LI/g > ${WD}/${SET}_sampleIds.txt
cut -d',' -f1 ${WD}/${SET}_sampleIds.txt > ${WD}/${SET}_logIds.txt
l1=$(wc -l ${WD}/${SET}_dnaIds.txt | cut -d' ' -f1)
l2=$(wc -l ${WD}/${SET}_logIds.txt | cut -d' ' -f1)
if [ $l1 != $l2 ]
  then
    echo "some DNA numbers have no log ids"
  else
    echo "all DNA numbers have log ids"
fi
}

# check if all samples in vcf
logid_in_vcf(){
bcftools query -l $VCF_DIR/Release4.norm.full.vcf.gz > ${WD}/Release4_samples.txt
l1=$(grep -F -f ${WD}/${SET}_logIds.txt Release4_samples.txt | wc -l)
l2=$(wc -l ${WD}/${SET}_logIds.txt | cut -d' ' -f1)
if [ $l1 != $l2 ]
  then
    echo "some log ids missing from vcf"
  else
    echo "all log ids in vcf"
fi
}

# get rare variants only for samples of interest
get_rare_variants(){
for i in $CHR;
do
  bcftools filter -i 'AF < 0.01' ${VCF_DIR}/chr${i}.norm.vcf.gz -O v | bcftools view -S ${WD}/${SET}_logIds.txt -O v -o ${WD}/${SET}_af0.01_chr${i}.vcf &
done
}

# only keep polymorphic variants
get_polymorphics(){
for i in $CHR;
do
  bcftools view -c1.0:minor ${WD}/${SET}_af0.01_chr${i}.vcf -O v -o ${WD}/${SET}_af0.01_poly_chr${i}.vcf &
done
}

# create list of variants for which to get total het and hom counts
# check tab space matches vcf tab delim, chage to query command
get_variants_list(){
for i in $CHR;
do
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ${WD}/${SET}_af0.01_poly_chr${i}.vcf > ${WD}/${SET}_chr${i}_targets.txt &
done
}

# figure out a way of getting variants from jvcf files without duplicating pos
targets_in_all_vcf(){
for i in $CHR;
do
  bcftools view -H -T ${WD}/${SET}_chr${i}_targets.txt ${VCF_DIR}/chr${i}.norm.vcf.gz -O v -o ${WD}/${SET}_targets_in_chr${i}.vcf &
done
}

variant_id_fields(){
for i in $CHR;
do
  cut -f1-5 ${WD}/${SET}_targets_in_chr${i}.vcf > ${WD}/${SET}_targets_in_chr${i}.info &
done
}

# create files per chr of het and hom counts
het_and_hom_counts(){
for i in $CHR;
do
  vcf=${WD}/${SET}_targets_in_chr${i}.vcf
  grep -n -o -F '0/1' $vcf | sort -n | uniq -c | cut -d':' -f1 > ${WD}/${SET}_chr${i}_het_counts.txt &
  grep -n -o -F '1/1' $vcf | sort -n | uniq -c | cut -d':' -f1 > ${WD}/${SET}_chr${i}_hom_counts.txt &
done
}

# reannotate vcfs
reannotate(){
for i in $CHR;
do
  bash /array/psivakumar/other/runVEP_psivakumar.sh ${WD}/${SET}_af0.01_poly_chr${i}.vcf &
done
}

# create info files
info_files(){
bcftools view -h ${WD}/${SET}_af0.01_poly_chr1.vcf.vep.out | grep CSQ > ${WD}/${SET}.csq
bcftools view -h ${WD}/${SET}_af0.01_poly_chr1.vcf.vep.out | tail -n 1 > ${WD}/${SET}.headers
}

# run R script per chr
R_formatting(){
for i in $CHR;
do
Rscript --vanilla /array/psivakumar/other/per_family_vcfs_format_open.R ${WD}/${SET}_af0.01_poly_chr${i}.vcf.vep.out ${WD}/${SET}.headers ${WD}/${SET}.csq ${WD}/${SET}_targets_in_chr${i}.info ${WD}/${SET}_chr${i}_het_counts.txt ${WD}/${SET}_chr${i}_hom_counts.txt ${WD}/${SET}_sampleIds.txt
done
}

# merge files
merge_files(){
rm -f ${WD}/${SET}_all_chr.csv
ls -v ${WD}/${SET}_af0.01_poly_chr*.vcf.vep.out.processed.csv | xargs cat >> ${WD}/${SET}_all_chr.csv
sed -e '1n' -e '/^CHROM/d' ${WD}/${SET}_all_chr.csv > ${WD}/${SET}_processed.csv
}

echo "checking for logids"
dnaid_to_logid
wait
if ! [ -s ${WD}/${SET}_sampleIds.txt ]; then
  echo "no valid logIds"
  exit
fi
echo "checking for logids in full vcf"
logid_in_vcf
wait
echo "filtering for rare variants"
get_rare_variants
wait
echo "filtering for polymorphic variants"
get_polymorphics
wait
echo "getting list of variants"
get_variants_list
wait
echo "filtering full vcf for target variants"
targets_in_all_vcf
wait
echo "getting variant id fields"
variant_id_fields
wait
echo "getting all het and hom counts"
het_and_hom_counts
wait
echo "reannotating"
reannotate
wait
echo "creating info files"
info_files
wait
echo "merging and formatting data"
R_formatting
wait
echo "merging all chromosomes"
merge_files
wait
echo "finished"
pwd
