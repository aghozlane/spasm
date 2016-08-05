#!/bin/bash
SLURM_SCRIPT=$HOME/assembly_submission.sh

#Check arguments
if [ $# -ne 7  ]
then
        echo "$0 <fna> <faa> <samplename>  <output_dir> <nb_cpu> <email> <queue>"
        exit
fi

SCRIPTPATH=$(dirname "${BASH_SOURCE[0]}")
gene=$(readlink -e $1)
prot=$(readlink -e $2)
mkdir -p $4
outdir=$(readlink -e $4)
echo """#!/bin/bash
#SBATCH --mail-user=$6
#SBATCH --mail-type=ALL
#SBATCH --qos=$7
#SBATCH --cpus-per-task=$5
#SBATCH --mem=50000
#SBATCH --job-name=\"assembly_${3}\"
### LIBRARY
source /local/gensoft2/adm/etc/profile.d/modules.sh
module purge
module add blast+/2.2.31 Python/2.7.8 fastqc/0.11.5 bowtie2/2.2.3 AlienTrimmer/0.4.0 SPAdes/3.9.0 hmmer/3.1b1 samtools/1.2 KronaTools/2.4 hmmer/3.1b1 barrnap/0.7
/bin/bash $SCRIPTPATH/assembly_illumina.sh  -g $gene -p $prot -s $3 -o $outdir -n $5 --tax_annotation --func_annotation &> $outdir/log_assembly_illumina.txt || exit 1
exit 0
""" >"${SLURM_SCRIPT}"
SLURMID=`sbatch ${SLURM_SCRIPT}`
echo "! Soumission SLURM :> JOBID = $SLURMID"

