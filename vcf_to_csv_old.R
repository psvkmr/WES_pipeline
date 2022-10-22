args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  print(length(args))
  stop("select vcf, headers, csq", call.=FALSE)
}

library(data.table)
library(bigreadr)
library(tidyverse)

chr <- args[1]
headers <- args[2]
csq <- args[3]

cmd.vcf <- paste0("grep -v '^#' ", chr)
vcf <- fread2(cmd=cmd.vcf, header = F, sep = '\t')
vcf.headers <- fread2(headers, header = F)
vcf.csq <- fread2(csq, header = F)
koios.counts <- fread('/hades/psivakumar/cohorts/koios_counts.csv')
omim <- fread('/hades/psivakumar/cohorts/omim/genemap2.txt', sep='\t')
panels <- fread('/hades/psivakumar/cohorts/200521_panelsWithQuality.csv')

vcf.headers <- vcf.headers %>% unlist() %>% gsub("#", "", .)
vcf.csq <- vcf.csq %>%
  separate(V1, c(NA, "V1"), sep = ":") %>%
  gsub(pattern = '">', replacement = "") %>%
  unlist()

colnames(vcf) <- vcf.headers

processInfo <- function(vcf, csq) {
  vcf2 <- vcf %>%
    mutate(CSQ = gsub(".*CSQ=", "", INFO)) %>%
    separate_rows(CSQ, sep = ",") %>%
    separate(CSQ, csq, sep = '\\|') %>%
    mutate_at(c("gnomADg_AF_AFR", "gnomADg_AF_AMR", "gnomADg_AF_ASJ", "gnomADg_AF_EAS", "gnomADg_AF_FIN", "gnomADg_AF_NFE", "gnomADg_AF_OTH"), as.numeric) %>%
    rowwise() %>%
    mutate(MAX_AF_g = max(gnomADg_AF_AFR, gnomADg_AF_AMR, gnomADg_AF_ASJ, gnomADg_AF_EAS, gnomADg_AF_FIN, gnomADg_AF_NFE, gnomADg_AF_OTH)) %>%
    replace_na(list(MAX_AF = 0, MAX_AF_g = 0))
    return(vcf2)
}

info.df <- processInfo(vcf, vcf.csq)

omim.df <- omim %>% select(`Ensembl Gene ID`, Phenotypes) %>% unique() %>% filter(grepl("ENSG",`Ensembl Gene ID`))
omim.df <-  cbind(aggregate(Phenotypes ~ `Ensembl Gene ID`, data = omim, paste, collapse = ",")) %>% filter(grepl("ENSG",`Ensembl Gene ID`))

processed <- left_join(info.df, koios.counts[, 1:6], by = c("CHROM" = "chrom", "POS" = "pos", "REF" = "ref", "ALT" = "alt"))
processed.join <- left_join(processed, omim.df, by = c("Gene"="Ensembl Gene ID")) %>% left_join(panels, by = c('Gene' = 'EnsemblGeneIds', 'SYMBOL' = 'GeneSymbol'))
#info.df <- mutate(info.df, MAX_AF = as.numeric(MAX_AF), MAX_AF_g = as.numeric(MAX_AF_g)) %>% filter((MAX_AF < 0.01 | is.na(MAX_AF)), (MAX_AF_g < 0.01 | is.na(MAX_AF_g)), BIOTYPE == "protein_coding")

#all.sample.df <- all.sample.info %>%
#  rownames_to_column("row_n") %>%
#  mutate_at("row_n", as.numeric) %>%
#  left_join(all.sample.het, by = c("row_n" = "V2")) %>%
#  left_join(all.sample.hom, by = c("row_n" = "V2")) %>%
#  dplyr::select(-row_n) %>%
#  `colnames<-`(c("CHROM", "POS", "ID", "REF", "ALT", "allHETCount", "allHOMCount")) %>%
#  replace_na(list(allHETCount = 0, allHOMCount = 0))

#processed <- left_join(info.df, all.sample.df, by = c("CHROM", "POS", "ID", "REF", "ALT"))

#setnames(processed, new = as.character(vcf.sample.ids$V2), old = vcf.sample.ids$V1)

fwrite(processed.join, paste0(chr, ".processed.csv"))
