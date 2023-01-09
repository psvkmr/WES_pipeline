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
#export PERL5LIB=/data/kronos/NGS_Software/VEP_95/ensembl-vep/biodbhts/lib/

VARIANT_FILE=$1
ASSEMBLY="38"
VEP_DIR="/hades/Software/NGS_Software/VEP_95"
#VCF_DIR="/mnt/qsg-results/inprogress/Release4"

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
--plugin MaxEntScan,$VEP_DIR/maxentscan/fordownload \
--plugin CADD,/mnt/qsg-results/inprogress/whole_genome_SNVs.tsv.gz,$VEP_DIR/InDels.tsv.gz\
--plugin Conservation \
--fork 10 \
-custom $VEP_DIR/gnomad.genomes.r2.0.1.sites.GRCh$ASSEMBLY.noVEP.vcf.gz,gnomADg,vcf,exact,0,\
AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH \
--fasta $VEP_DIR/Homo_sapiens.GRCh$ASSEMBLY.dna.primary_assembly.fa

###############################################################
