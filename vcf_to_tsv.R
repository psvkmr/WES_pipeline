#$bcftools view -h  annotated/full_anno_vep.vcf.gz | grep CSQ > ${out_dir}/ofg_vep_anno_csq.txt
#$bcftools view -h ${out_dir}/ofg_chr22.vcf.gz.batch.vep.vcf | tail -n 1 > ${out_dir}/ofg_vep_anno_headers.txt

# create tsv file convert from vcf for each chr with headers and csq as additional args
#for i in `seq 1 22`;
#do
#  sbatch ${base_dir}/scripts/misc_scripts/vcf_to_tsv.sh ${out_dir}/ofg_chr${i}.vcf.gz.batch.vep.vcf ${out_dir}/ofg_vep_anno_headers.txt ${out_dir}/ofg_vep_anno_csq.txt
#done

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  print(length(args))
  stop('select vcf, headers, csq', call.=FALSE)
}

library(data.table)
library(tidyverse)

chr <- args[1]
headers <- args[2]
csq <- args[3]

print('reading files...')
vcf <- read.table(chr, header = F, sep = '\t')
print('finished reading vcf')
vcf.headers <- read.table(headers, header = F, comment.char = '')
print('finished reading headers')
vcf.csq <- fread(csq, sep = '|', header = F)
print('finished reading csq')

vcf.headers <- vcf.headers %>% unlist() %>% gsub('#', '', .)
vcf.csq <- vcf.csq %>%
  separate(V1, c(NA, 'V1'), sep = ': ') %>%
  gsub(pattern = '">', replacement = '') %>%
  unlist()

colnames(vcf) <- vcf.headers

processCsq <- function(vcf, csq) {
  vcf2 <- vcf %>%
    mutate(CSQ = gsub('.*CSQ=', '', INFO)) %>%
    separate_rows(CSQ, sep = ',') %>%
    separate(CSQ, vcf.csq, sep = '\\|')
#    mutate_at(vars(contains('gnomAD'), as.numeric)) %>%
#    mutate('SIFT_score' = as.numeric(str_extract(.$SIFT, '0\\.*[0-9]*')), 
#           'PolyPhen_score' = as.numeric(str_extract(.$PolyPhen, '0\\.*[0-9]*')),
#           'CADD_PHRED' = as.numeric(CADD_PHRED))
  return(vcf2)
}
info.df <- processCsq(vcf, vcf.csq)
print('finished formatting csq')


processInfo <- function(csq_df){
  info <- csq_df %>%
    mutate(ASTATS = gsub(';CSQ.*$', '', INFO)) %>%
    separate(ASTATS, c('VCF_AF', 'AQ', 'AC', 'AN'), sep = ';')
#    mutate_at(vars(c('VCF_AF', 'AQ', 'AC', 'AN'), ~ gsub('^.*=', '', .))) %>%
#    mutate_at(vars(c('VCF_AF', 'AQ', 'AC', 'AN'), as.numeric))
  return(info)
}
info.df <- processInfo(info.df)

print('writing file...')
write.table(info.df, paste0(chr, '.processed.tsv'), sep = '\t', row.names=F, quote=F)
print('done')