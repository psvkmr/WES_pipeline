SAMPLES=$(ls -d /hades/psivakumar/ICGNMD/24_02_2020/*/ | cut -d'/' -f6)
ARR_S=($SAMPLES)

#${!ARR_S[@]}
for i in `seq 0 1`
do
    bash /hades/psivakumar/pipeline/psivakumar_WES_bwa.sh /hades/psivakumar/ICGNMD/24_02_2020/${ARR_S[i]}/${ARR_S[i]}_R1.fastq.gz /hades/psivakumar/ICGNMD/24_02_2020/${ARR_S[i]}/${ARR_S[i]}_R2.fastq.gz /hades/psivakumar/ICGNMD/24_02_2020/processed_bwa/${ARR_S[i]} ${ARR_S[i]}
done
