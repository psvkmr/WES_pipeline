#!/usr/bin/env bash
# ##################################################
# My Generic BASH script template
#
version="1.0.0"               # Sets version variable
#
scriptTemplateVersion="1.5.0" # Version of scriptTemplate.sh that this script is based on
#                               v1.1.0 -  Added 'debug' option
#                               v1.1.1 -  Moved all shared variables to Utils
#                                      -  Added $PASS variable when -p is passed
#                               v1.2.0 -  Added 'checkDependencies' function to ensure needed
#                                         Bash packages are installed prior to execution
#                               v1.3.0 -  Can now pass CLI without an option to $args
#                               v1.4.0 -  checkDependencies now checks gems and mac apps via
#                                         Homebrew cask
#                               v1.5.0 - Now has preferred IFS setting
#                                      - Preset flags now respect true/false
#                                      - Moved 'safeExit' function into template where it should
#                                        have been all along.
#
# HISTORY:
#
# * DATE - v1.0.0  - First Creation
#
# ##################################################
# Provide a variable with the location of this script.
scriptPath="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Source Scripting Utilities
# -----------------------------------
# These shared utilities provide many functions which are needed to provide
# the functionality in this boilerplate. This script will fail if they can
# not be found.
# -----------------------------------
utilsLocation="/home/hades/pipeline/shell-scripts/lib/utils.sh" # Update this path to find the utilities.
if [ -f "${utilsLocation}" ]; then
  source "${utilsLocation}"
else
  echo "Please find the file util.sh and add a reference to it in this script. Exiting."
  exit 1
fi
# trapCleanup Function
# -----------------------------------
# Any actions that should be taken if the script is prematurely
# exited.  Always call this function at the top of your script.
# -----------------------------------
function trapCleanup() {
  echo ""
  # Delete temp files, if any
  if is_dir "${tmpDir}"; then
    rm -r "${tmpDir}"
  fi
  die "Exit trapped. In function: '${FUNCNAME[*]}'"
}
# safeExit
# -----------------------------------
# Non destructive exit for when script exits naturally.
# Usage: Add this function at the end of every script.
# -----------------------------------
function safeExit() {
  # Delete temp files, if any
  if is_dir "${tmpDir}"; then
    rm -r "${tmpDir}"
  fi
  trap - INT TERM EXIT
  exit
}
# Set Flags
# -----------------------------------
# Flags which can be overridden by user input.
# Default values are below
# -----------------------------------
quiet=false
printLog=true
verbose=false
force=false
strict=false
debug=false
args=()
# Set Temp Directory
# -----------------------------------
# Create temp directory with three random numbers and the process ID
# in the name.  This directory is removed automatically at exit.
# -----------------------------------
tmpDir="/tmp/${scriptName}.$RANDOM.$RANDOM.$RANDOM.$$"
(umask 077 && mkdir "${tmpDir}") || {
  die "Could not create temporary directory! Exiting."
}
# Logging
# -----------------------------------
# Log is only used when the '-l' flag is set.
#
# To never save a logfile change variable to '/dev/null'
# Save to Desktop use: $HOME/Desktop/${scriptBasename}.log
# Save to standard user log location use: $HOME/Library/Logs/sync.log
# -----------------------------------
logFile="/hades/dmurphy/pipeline/pipeline.log"
# Check for Dependencies
# -----------------------------------
# Arrays containing package dependencies needed to execute this script.
# The script will fail if dependencies are not installed.  For Mac users,
# most dependencies can be installed automatically using the package
# manager 'Homebrew'.  Mac applications will be installed using
# Homebrew Casks. Ruby and gems via RVM.
# -----------------------------------
homebrewDependencies=()
caskDependencies=()
gemDependencies=()
function mainScript() {
############## Begin Script Here ###################
####################################################
echo ${args[0]}
echo ${args[1]}
echo ${args[2]}
ID=${args[0]}

file1=${args[1]}
file2=${args[2]}
file3=${args[3]}
file4=${args[4]}
file5=${args[5]}
file6=${args[6]}
file7=${args[7]}
file8=${args[8]}
file9=${args[9]}
file10=${args[10]}
file11=${args[11]}
file12=${args[12]}
file13=${args[13]}
file14=${args[14]}
file15=${args[15]}
file16=${args[16]}
file17=${args[17]}
file18=${args[18]}
file19=${args[19]}
file20=${args[20]}
file21=${args[21]}
file22=${args[22]}
file23=${args[23]}
file24=${args[24]}
file25=${args[25]}
file26=${args[26]}
file27=${args[27]}
file28=${args[28]}
file29=${args[29]}
file30=${args[30]}
file31=${args[31]}
file32=${args[32]}

part1=""
part2=""
part3=""
part4=""
part5=""
part6=""
part7=""
part8=""
part9=""
part10=""
part11=""
part12=""
part13=""
part14=""
part15=""
part16=""


part1bam=""
part2bam=""
part3bam=""
part4bam=""
part5bam=""
part6bam=""
part7bam=""
part8bam=""
part9bam=""
part10bam=""
part11bam=""
part12bam=""
part13bam=""
part14bam=""
part15bam=""
part16bam=""



GRCH38="/hades/dmurphy/pipeline/Homo_sapiens_assembly38.fasta"
GRCH38index="/hades/dmurphy/pipeline/Homo_sapiens_assembly38.index"
knownINDELS="/hades/dmurphy/pipeline/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
echo $GRCH38
echo $file1
echo $file2
echo $ID



mkdir /mnt/qsg-results/inprogress/$ID
cd /mnt/qsg-results/inprogress/$ID

#if $# > 3
#for arg in args
#gunzip fastq >> /mnt/qsg-results/inprogress/$ID/fastq$ID.fq


ID=${args[0]}


#set -o xtrace
#set -o verbose

declare -A parts=()
declare -A bams=()

x=1
for ((i=1;i<=32;i=i+2)); 
do 
   if [ ! -z ${args[i]} ]
	then
		parts[$x]='part'$x'_'$ID'.sam'
		bams[$x]='./part'${x}'_'$ID'.bam'
		/hades/Software/NGS_Software/novocraftV3.08.02/nooalign -c 10 -d $GRCH38index -f ${args[i]} ${args[i+1]} -o SAM $"@RG\tID:1\tPU:Unit\tLB:library\tSM:$ID\tPL:illumina"  --rOQ --hdrhd 3 -H 15 -k -o Soft -t 320  > ${parts[$x]}
		samtools view -Sb ${parts[$x]} > ${bams[$x]}
		rm ${parts[$x]}
		x=$((x+1))
	fi
done
echo "complete"

bamlist=""
for ((i=1;i<=${#bams[@]};i=i+1)); 
do 
    bamlist=$bamlist' '${bams[$i]}
done

echo $bamlist

/hades/Software/NGS_Software/novocraftV3.08.02/novosort --md --kt --ise -c 10 -t /hades/pipelinetemp/ -f -i -o Sorted_1_$ID.bam  $part1bam $part2bam $part3bam $part4bam

for ((i=1;i<=${#bams[@]};i=i+1)); 
do 
    bamlist=$bamlist' '${bams[$i]}
done

/hades/Software/NGS_Software/gatk-4.0.0.0/gatk-4.0.0.0/gatk --java-options "-Xmx4g" BaseRecalibrator  -I Sorted_1_$ID.bam  -R $GRCH38  --known-sites /hades/dmurphy/pipeline/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz -O recal_data.table

/hades/Software/NGS_Software/gatk-4.0.0.0/gatk-4.0.0.0/gatk --java-options "-Xmx4g" AnalyzeCovariates -bqsr recal_data.table  -plots AnalyzeCovariates.pdf

/hades/Software/NGS_Software/gatk-4.0.0.0/gatk-4.0.0.0/gatk --java-options "-Xmx4g" ApplyBQSR  -R $GRCH38  -I Sorted_1_$ID.bam  --bqsr-recal-file recal_data.table -O Sorted_$ID.bam



/hades/Software/NGS_Software/gatk-4.0.0.0/gatk-4.0.0.0/gatk --java-options "-Xmx4g" HaplotypeCaller  \
   --intervals /hades/dmurphy/pipeline/liftedover_SureSelect_V6_UTR.bed \
   -R $GRCH38 \
   -I Sorted_$ID.bam \
   -O $ID.g.vcf.gz \
   -ERC GVCF
     
   
/hades/Software/NGS_Software/gatk-4.0.0.0/gatk-4.0.0.0/gatk Mutect2 -L chrM:1-16569 -R  $GRCH38  -I Sorted_$ID.bam   -O MT_$ID.vcf -tumor $ID
/hades/Software/NGS_Software/gatk-4.0.0.0/gatk-4.0.0.0/gatk FilterMutectCalls   -V MT_$ID.vcf    -O Filtered_MT_$ID.vcf  --max-germline-posterior 1000 --max-events-in-region 1000
/hades/Software/NGS_Software/gatk-4.0.0.0/gatk-4.0.0.0/gatk SelectVariants  -R $GRCH38 -V Filtered_MT_$ID.vcf  -O Filtered_Passing_MT_$ID.vcf  --exclude-filtered true

samtools index Sorted_$ID.bam

/hades/Software/NGS_Software/verifyBamID   --vcf /hades/dmurphy/pipeline/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz --bam Sorted_$ID.bam --out VBID$ID


/hades/Software/NGS_Software/gatk-4.0.0.0/gatk-4.0.0.0/gatk --java-options "-Xmx4g"  BamIndexStats  -I Sorted_$ID.bam > indexstats.txt

/hades/Software/NGS_Software/gatk-4.0.0.0/gatk-4.0.0.0/gatk --java-options "-Xmx4g"  CollectWgsMetrics  -I Sorted_$ID.bam -O wgs_metrics.txt -R /hades/dmurphy/pipeline/Homo_sapiens_assembly38.fasta --INTERVALS /hades/dmurphy/pipeline/liftedover_SureSelect_V6_UTR.interval_list


#catching the exit status codes in case this is a re-run and some of the files do not exist. 

t1=$(rm -f  $ID.bam)
t2=$(rm -f   $part1)
t2=$(rm -f  $part2)
t2=$(rm -f  $part3)
t2=$(rm -f  $part4)
t2=$(rm -f  $part1bam)
t2=$(rm -f  $part2bam)
t2=$(rm -f  $part3bam)
t2=$(rm -f  $part4bam)
t3=$(rm -f Sorted_1_$ID.bam)

mv /mnt/qsg-results/inprogress/$ID/ /mnt/qsg-results/pipeline/$ID/

echo "complete"

####################################################
############### End Script Here ####################
}
############## Begin Options and Usage ###################
# Print usage
usage() {
  echo -n "${scriptName} [OPTION]... [FILE]...
This is a script template.  Edit this description to print help to users.
 ${bold}Options:${reset}
  -u, --username    Username for script
  -p, --password    User password
  --force           Skip all user interaction.  Implied 'Yes' to all actions.
  -q, --quiet       Quiet (no output)
  -l, --log         Print log to file
  -s, --strict      Exit script with null variables.  i.e 'set -o nounset'
  -v, --verbose     Output more information. (Items echoed to 'verbose')
  -d, --debug       Runs script in BASH debug mode (set -x)
  -h, --help        Display this help and exit
      --version     Output version information and exit
"
}
# Iterate over options breaking -ab into -a -b when needed and --foo=bar into
# --foo bar
optstring=h
unset options
while (($#)); do
  case $1 in
    # If option is of type -ab
    -[!-]?*)
      # Loop over each character starting with the second
      for ((i=1; i < ${#1}; i++)); do
        c=${1:i:1}
        # Add current char to options
        options+=("-$c")
        # If option takes a required argument, and it's not the last char make
        # the rest of the string its argument
        if [[ $optstring = *"$c:"* && ${1:i+1} ]]; then
          options+=("${1:i+1}")
          break
        fi
      done
      ;;
    # If option is of type --foo=bar
    --?*=*) options+=("${1%%=*}" "${1#*=}") ;;
    # add --endopts for --
    --) options+=(--endopts) ;;
    # Otherwise, nothing special
    *) options+=("$1") ;;
  esac
  shift
done
set -- "${options[@]}"
unset options
# Print help if no arguments were passed.
# Uncomment to force arguments when invoking the script
# [[ $# -eq 0 ]] && set -- "--help"
# Read the options and set stuff
while [[ $1 = -?* ]]; do
  case $1 in
    -h|--help) usage >&2; safeExit ;;
    --version) echo "$(basename $0) ${version}"; safeExit ;;
    -u|--username) shift; username=${1} ;;
    -p|--password) shift; echo "Enter Pass: "; stty -echo; read PASS; stty echo;
      echo ;;
    -v|--verbose) verbose=true ;;
    -l|--log) printLog=true ;;
    -q|--quiet) quiet=true ;;
    -s|--strict) strict=true;;
    -d|--debug) debug=true;;
    --force) force=true ;;
    --endopts) shift; break ;;
    *) die "invalid option: '$1'." ;;
  esac
  shift
done
# Store the remaining part as arguments.
args+=("$@")
############## End Options and Usage ###################
# ############# ############# #############
# ##       TIME TO RUN THE SCRIPT        ##
# ##                                     ##
# ## You shouldn't need to edit anything ##
# ## beneath this line                   ##
# ##                                     ##
# ############# ############# #############
# Trap bad exits with your cleanup function
trap trapCleanup EXIT INT TERM
# Set IFS to preferred implementation
IFS=$'\n\t'
# Exit on error. Append '||true' when you run the script if you expect an error.
set -o errexit
# Run in debug mode, if set
if ${debug}; then set -x ; fi
# Exit on empty variable
if ${strict}; then set -o nounset ; fi
# Bash will remember & return the highest exitcode in a chain of pipes.
# This way you can catch the error in case mysqldump fails in `mysqldump |gzip`, for example.
set -o pipefail
# Invoke the checkDependenices function to test for Bash packages.  Uncomment if needed.
# checkDependencies
# Run your script
mainScript
# Exit cleanlyd
safeExit
