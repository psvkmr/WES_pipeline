#This is a script to extract variant annotation from Ensembl VEP version 95. Result is 
#a vcf file (vep.out) whihc includes the most severe variant consequence per gene.

#Shortcut flag to switch on all of the following:
# --sift b, --polyphen b, --ccds, --uniprot, --hgvs, --symbol, 
#--numbers, --domains, --regulatory, --canonical, --protein, 
#--biotype, --uniprot, --tsl, --appris, --gene_phenotype 
#--af, --af_1kg, --af_esp, --af_gnomad, --max_af, --pubmed, 
#--variant_class

#run as qsub -cwd -pe make 10 script.sh VARIANT_FILE ASSEMBLY
#############################################################
#export PERL5LIB=/hades/Software/NGS_Software/VEP_95/ensembl-vep/biodbhts/lib/
#PERL5LIB=${PERL5LIB}:/usr/local/lib/x86_64-linux-gnu/perl/5.22.1/ 
#export PERL5LIB
#PERL5LIB=${PERL5LIB}:/hades/Software/NGS_Software/VEP_95/VEP_plugins/
#export PERL5LIB
#PERL5LIB=${PERL5LIB}:/hades/Software/NGS_Software/VEP_95/ensembl-vep/biodbhts/lib/Bio/DB/HTS/
#export PERL5LIB
#export PERL5LIB=/data/kronos/NGS_Software/VEP_plugins/

echo $PERL5LIB

VARIANT_FILE=$1
ASSEMBLY=38
VEP_DIR="/hades/Software/NGS_Software/VEP_95"

$VEP_DIR/ensembl-vep/vep -i $VARIANT_FILE \
--everything \
--per_gene \
--pick_order rank \
--check_existing \
--offline \
--cache \
--dir_cache $VEP_DIR/ensembl-vep \
--cache_version 95 \
-o $VARIANT_FILE.vep.out \
--assembly GRCh$ASSEMBLY \
--vcf \
--no_stats \
--plugin SpliceRegion \
--fork 10 \
-custom $VEP_DIR/gnomad.genomes.r3.0.sites.vcf.bgz,gnomADg,vcf,exact,0,\
AF_afr,AF_ami,AF_amr,AF_asj,AF_eas,AF_fin,AF_nfe,AF_oth,AF_raw,AF_sas,faf95_adj,faf95_afr,faf95_amr,faf95_eas,faf95_nfe,faf95_sas,faf99_adj,faf99_afr,faf99_amr,faf99_eas,faf99_nfe,faf99_sas \
--fasta $VEP_DIR/Homo_sapiens.GRCh$ASSEMBLY.dna.primary_assembly.fa \
--plugin CADD,/hades/Software/NGS_Reference/plugins/CADD_GRCH38_whole_genome_SNVs.tsv.gz,/hades/Software/NGS_Reference/plugins/CADD_GRCH38_InDels.tsv.gz \
--plugin SpliceAI,snv=/hades/Software/NGS_Reference/plugins/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/hades/Software/NGS_Reference/plugins/spliceai_scores.raw.indel.hg38.vcf.gz \
--plugin MaxEntScan,/hades/Software/NGS_Reference/plugins/fordownload,SWA,NCSS \
--plugin Carol \
--plugin satMutMPRA,file=/hades/Software/NGS_Reference/plugins/satMutMPRA_GRCh38_ALL.gz,cols=ALL
###############################################################

