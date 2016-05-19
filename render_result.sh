#!/bin/bash
#Check arguments
if [ $# -ne 2  ]
then
        echo "$0 <result_dir> <render_dir>"
        exit
fi

for sample in $1/*/
do
    sample_name=$(basename $sample)
    mkdir $2/$sample_name
    rm -f $sample/tmp_fastq
    #find $sample/ ! -name *.sam -exec cp -t $2/$sample_name/ {} +
    cp $sample/*.txt $sample/*.html $sample/*.fna $sample/*.metagenemark $sample/*.faa   $2/$sample_name/
    cp -R $sample/fetch_result $2/$sample_name/
    cp $sample/spades/scaffolds.fasta $2/$sample_name/${sample_name}_scaffolds.fasta
    cp -R $sample/filter_alientrimmer/*.html $2/$sample_name/
    rm -f $2/$sample_name/*_interest_list.txt
done