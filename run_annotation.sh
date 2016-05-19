#!/bin/bash
PBS_SCRIPT=$HOME/assembly_submission.sh

#Check arguments
if [ $# -ne 5  ]
then
        echo "$0 <fna> <faa> <samplename>  <output_dir> <nb_cpu>"
        exit
fi

gene=$(readlink -e $1)
prot=$(readlink -e $2)
outdir=$(readlink -e $4)
echo """#!/bin/bash
#$ -S /bin/bash
#$ -M amine.ghozlane@pasteur.fr
#$ -m bea
#$ -q hubbioit
##$ -q test
#$ -pe thread $5
#$ -l mem_total=50G
### LIBRARY
source /local/gensoft2/adm/etc/profile.d/modules.sh
module purge
module add blast+/2.2.31 Python/2.7.8 fastqc/0.11.2 bowtie2/2.2.3 AlienTrimmer/0.4.0 SPAdes/3.6.2 hmmer/3.1b1 samtools/1.2 KronaTools/2.4 hmmer/3.1b1 barrnap/0.7
#$ -N "${3}"
/bin/bash $HOME/assembly_illumina/assembly_illumina.sh  -g $gene -p $prot -s $3 -o $outdir -n $5 --tax_annotation --func_annotation &> $outdir/log_assembly_illumina.txt 
""" >"${PBS_SCRIPT}"
PBSID=`qsub ${PBS_SCRIPT}`
#exit
echo "! Soumission PBS :> JOBID = $PBSID"
