source `which env_parallel.bash`
trap "exit" INT

################################################

export PATH_TO_NANOSORT="/Users/heatherwells/Documents/Recombination_Paper_1/final/NanoSort"

################################################

export PARENT1=${1}
export PARENT2=${2}
export INPUT_READS=${3}
export OUT=$(pwd)/${INPUT_READS}_${PARENT1}_${PARENT2}_recombination_test
export PYTHONPATH="${PATH_TO_NANOSORT}/:${OUT}"
export w=200 #sliding window size
export pc=41 #minimum number of required sub-reads
export s=10 #step size

if [ ! -d $OUT ]; then
  mkdir ${OUT}
fi

BLAST_RESULTS=${OUT}/${INPUT_READS}_blast_results.out
echo "BLASTing reads"
seqkit fq2fa ${INPUT_READS} | blastn -db ~/bin/ONT_recomb/data/pairs/${PARENT1}_${PARENT2}/BLAST/${PARENT1}_${PARENT2}.fasta -outfmt "6 qseqid" -out ${BLAST_RESULTS}

export FILTERED_READS=${OUT}/${INPUT_READS}_filtered.fastq.gz
echo "filtering viral reads"
seqkit grep --pattern-file ${BLAST_RESULTS} ${INPUT_READS} | seqkit seq -m 1000 -o ${FILTERED_READS} 

wait

if [ ! -d ${OUT}/reads ]; then
  mkdir ${OUT}/reads
fi

if [ ! -f ${OUT}/overall_summary.tsv ]; then
  echo -e "read_name\trecomb\tnonhom\tsubgen\n" > ${OUT}/overall_summary.tsv
fi

cd $PATH_TO_NANOSORT
python set_config.py -i "${INPUT_READS}" -p1 "${PARENT1}" -p2 "${PARENT2}" -s "${s}" -w "${w}" -pc "${pc}" -wd "${OUT}" -nd ${PATH_TO_NANOSORT}

echo "classifying viral reads"
gzcat ${FILTERED_READS} | env_parallel --env OUT --env PYTHONPATH --env PARENT1 --env PARENT2 --pipe -n 4 -j 32 ./parallel.sh