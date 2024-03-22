read first_line; read second_line; read third_line; read -r fourth_line

name=${first_line}
namearray=(${name})
read_name=${namearray[0]}

FILE=${OUT}/reads/${read_name}

printf "%s\n%s\n%s\n%s" $read_name $second_line $third_line $fourth_line > "${FILE}.fastq"

seqkit sliding -s "${s}" -W "${w}" -o "${FILE}_sliding.fastq" "${FILE}.fastq"

if [ $(stat -f%z $FILE.fastq) -eq 0 ]; then exit 1; fi

minimap2 "data/pairs/${PARENT1}_${PARENT2}/${PARENT1}_${PARENT2}.mmi" "${FILE}_sliding.fastq" -m 20 -w 1 --secondary=no -a | samtools view -F 2048 > "${FILE}_sliding.sam" 

minimap2 "data/pairs/${PARENT1}_${PARENT2}/${PARENT1}_${PARENT2}.mmi" "${FILE}.fastq" -m 20 -w 1 --secondary=no -a  > "${FILE}.sam" 

python3 base_task.py -i $FILE -o ${OUT}

wait

rm ${FILE}_sliding.sam ${FILE}_sliding.fastq ${FILE}.sam ${FILE}.fastq
