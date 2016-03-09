#!/bin/bash

#check arguments
if [ $# -ne "2" ]
then
     echo "Usage : $0 <result_dir> <fichier_output>"
    exit
fi

echo "Sample;NbReads;NbFilteredReads;NbScaffolds;N50;MeanScaffoldslength;MedianScaffoldslength;NbGenes;NbFilteredGenes;MeanGenelength;MedianGenelength;NbMarker;NbNCBI;NbEGGNOG;NbCazy" > $2
for dir in $(ls -d $1/*/)
do
#echo "$dir"
#food=$(basename $(dirname $(dirname $dir)))
#baby=$(basename $(dirname $dir))
sample=$(basename $dir)
#echo $sample
NbFilteredReads=`grep "^+" -c ${dir}/filter_2/un-conc-mate_1.fastq|cut -f 2 -d ":"`
NbAlien=`grep "^+" -c ${dir}/filter_alientrimmer/un-conc-mate_1.fastq`
Stats_scaf=`$HOME/get_sequence_length/get_sequence_length.py -i ${dir}/*_scaffolds.fasta | egrep -i "^Mean|^Median|^N50" | cut -d ':' -f2`
Meanscaffoldslength=`echo $Stats_scaf | cut -f1 -d ' '`
Medianscaffoldslength=`echo $Stats_scaf | cut -f2 -d ' '`
N50=`echo $Stats_scaf | cut -f3 -d ' '`

NbScaffolds=`grep "^>" -c ${dir}/*_scaffolds.fasta`
#LonguestScaffold=`grep "^>" ${dir}/*_scaffolds.fasta -m1 |cut -f 2 -d "e"`

#echo "nb contigs: ${nb_contigs}"
Nbgenes=`grep "^>" -c ${dir}/*gene.fna`
Nbfilteredgenes=`grep "^>" -c ${dir}/*gene_60_filtered.fna`
#echo "nb gènes: ${nb_genes}"
# Nbfilteredproteins=`grep "^>" -c ${dir}/*prot_filtered.faa`
#NbCogs=`head -n 1 ${dir}/*prot_cog_stat.txt |cut -c 1-2`
#echo "nb cogs: ${num_cogs}"
NbMarker=$(grep "^>" -c ${dir}/*_gene_marker.fna)
NbNCBI=$(wc -l ${dir}/*_gene_60_blastn_ncbi_genome_cov80.txt|cut -f1 -d " ")
let "NbNCBI=$NbNCBI - 1"
NbEGGNOG=$(wc -l ${dir}/*_prot_eggnog_annotation_filtered.txt|cut -f1 -d " ")
let "NbEGGNOG=$NbEGGNOG - 1"
NbCazy=$(wc -l ${dir}/*_prot_cazy_annotation_filtered.txt|cut -f1 -d " ")
let "NbCazy=$NbCazy - 1"
#echo "nb species: ${num_species}"
#NbGenera=`cut -s -f 2 ${dir}/*gene_100_blastn_annotation_genome_id85_cov90.txt | sort | uniq | wc -l`
#echo "nb genera: ${num_genera}"
#num_n50=`more ${dir}/*n50_stats.txt`
#n50=`grep "^N50 value"  ${dir}/*n50_stats.txt|cut -d ":" -f2  | sed -e 's/^ *//' -e 's/ *$//'`
#echo "nb n50: ${num_n50}"
#Mappedreadstocontigs=`grep "# reads with at least one reported alignment:" ${dir}/log/log_unmapped_reads_*.txt| cut -d '(' -f2| cut -d ')' -f 1`
#Unmappedreadstocontigs=`grep "# reads that failed to align:" ${dir}/log/log_unmapped_reads_*.txt| cut -d '(' -f2| cut -d ')' -f 1`
Stats=`$HOME/get_sequence_length/get_sequence_length.py -i ${dir}/*gene_60_filtered.fna | egrep -i "^Mean|^Median" | cut -d ':' -f2`
Meangeneslength=`echo $Stats | cut -f1 -d ' '`
Mediangeneslength=`echo $Stats | cut -f2 -d ' '`

echo "${sample};$NbFilteredReads;$NbAlien;$NbScaffolds;$N50;$Meanscaffoldslength;$Medianscaffoldslength;$Nbgenes;$Nbfilteredgenes;$Meangeneslength;$Mediangeneslength;$NbMarker;$NbNCBI;$NbEGGNOG;$NbCazy"
#cut -s -f 3 *gene_100_blastn_annotation_genome.txt | sort | uniq -c
#seqs_align=`grep '^AlignedSeqs' ${dir}/*.report | sed 's/\s\+/\t/g' | cut -f3`
#echo "nb seqs alignées: ${seqs_align}"
#bases_align=`grep '^AlignedBases' ${dir}/*.report | sed 's/\s\+/\t/g' | cut -f3`
#echo "nb bases alignées: ${bases_align}"

#output="xstoocky14:/glusterfs/gvrepli/mgps/home/asalvarez/assembly_baby_analysis/$dir/"
#mkdir -p $output

done  >> $2