#!/bin/bash
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html
# ------------------------------------------------------------------
# Authors: Amine Ghozlane (amine.ghozlane@jouy.inra.fr),
# Title:  Assemble proton single end data
# Description : Assembly pipeline for proton reads with annotation
# ------------------------------------------------------------------

function say_parameters {
    # echo blue
    echo -e "\e[34m# $1\e[0m" >&2
}

function say {
    # echo green
    echo -e "\033[1;32m* $1\033[0m" >&2
}

function error {
    # echo red
    echo -e  "\e[31m* $1\e[0m" >&2
}

function check_log {
    # check if log file is not empty
    if [ -s $1 ]
    then
        error "$1 is not empty !"
        exit 1
    fi
}

function check_integer {
    if [[ $1 != [0-9]* ]]
    then
        error "\"$1\" is not an integer value"
        exit 1
    fi
}

function check_file {
    # Check if result is well produced
    if [ ! -f $1 ] && [ ! -s $1 ]
    then
        error "File \"$1\" does not exist or is empty !"
        exit 1
    fi
}

function check_dir {
    # Check if directory doesnt exist
    if [ ! -d $1 ]
    then
        mkdir $1
        if [ ! -d $1 ]
        then
            error "The program cannot create the directory \"$1\""
            exit 1
        fi
    fi
}

display_help() {
    if [ "$1" -eq "0" ]
    then
        echo """Usage : 
case 1 - Assembly: $0 -1 <read.fastq> -2 <read.fastq> -s <sample_name> -o </path/to/result/directory/> --assembly
case 2 - Taxonomic annotation/MTU: $0 -g <gene.fasta> -o </path/to/result/directory/> --tax_annotation
case 3 - Functional annotation: $0 -p <protein.fasta> -o </path/to/result/directory/> --func_annotation
all : $0 -1 <read.fastq> -2 <read.fastq> -s <sample_name> -o </path/to/result/directory/>
        """
        echo "Use '--help' to print detailed descriptions of command line arguments"
    else
        display_parameters
    fi
    exit
}

display_parameters() {
   # Display the parameters of the analysis
   say_parameters "Sample fastq [-1] [-2] :"
   echo $input >&2
   say_parameters "Sample prefix [-r] :"
   echo $SamplePrefix >&2
   say_parameters "Input gene [-g] :"
   echo $input_gene >&2
   say_parameters "Input protein [-p] :"
   echo $input_protein >&2
   say_parameters "Output directory [-o] :"
   echo $resultDir  >&2
   say_parameters "Filtering reads from the databases :"
   for db in ${filterRef[@]}
   do
    echo $db  >&2
   done
   say_parameters "Assembly parameters :"
   if [ "$readTrim" -eq "1" ]
   then
    echo "Trimming reads [-t]" >&2
   fi
   echo "Number of mismatch for the mapping [-m] = $NbMismatchMapping">&2
   echo "Gene length threshold [-c] = $geneLengthThreshold" >&2
   echo "Number of process [-n|--NbProc] = $NbProc" >&2
   say_parameters "Essential genes :"
#    echo "MOTU [--markerCOGRef] = $markerCOGRef" >&2
#    echo "MOTU blast db [--markerCOGRefBlast] = $markerCOGRefBlast" >&2
   say_parameters "Filter gene predited with Prodigal :"
   echo "Use metagenemark [--metagenemark] = $metagenemark"
   echo "Gene type [-f] = $filterGene"
   say_parameters "-Taxonomic annotation-"  
   
   say_parameters "-Functional annotation-"
   say_parameters "Egg/Nog annotation :"
   echo "Egg/Nog sequence database [--eggnogRef] = $eggnogRef" >&2
   echo "Egg/Nog sequence database [--eggnogRefBlast] = $eggnogRefBlast" >&2
   echo "Egg/Nog description [--eggnogdescRef] = $eggnogdescRef" >&2
   echo "Egg/Nog members [--eggnogmemRef] = $eggnogmemRef" >&2
   echo "Egg/Nog functional category description [--eggnogfunccats] = $eggnogfunccats" >&2
   echo "Egg/Nog functional category [--eggnogfunccatRef] = $eggnogfunccatRef" >&2
   say_parameters "Cazy annotation :"
   echo "Cazy HMM db [--cazyRefHmmer] = $cazyRefHmmer" >&2
   echo "Cazy family information [--famInfoRef] = $famInfoRef" >&2
#    say_parameters "Kegg annotation :"
#    echo "Kegg sequence db [--keggRef] = $keggRef" >&2
#    echo "Kegg sequence db [--keggRefBlast] = $keggRefBlast" >&2
#    echo "KO ref [--koRef] = $koRef" >&2
#    echo "KO gene to ko [--koGeneKoRef] = $koGeneKoRef" >&2
#    echo "KO gene to uniprot [--koGeneUniprotRef] = $koGeneUniprotRef" >&2
#    echo "KO to uniprot [--koUniprotRef] = $koUniprotRef" >&2
#    echo "KO gene to geneid [--koGeneGeneidRef] = $koGeneGeneidRef" >&2
#    echo "KO gene to gi [--koGeneGiRef] = $koGeneGiRef" >&2
#    echo "KO gene to pathway [--koGenePathwayRef] = $koGenePathwayRef" >&2
}

function timer()
{
   if [[ $# -eq 0 ]]; then
         echo $(date '+%s')
   else
      local  stime=$1
      etime=$(date '+%s')
      if [[ -z "$stime" ]]; then stime=$etime; fi
      dt=$((etime - stime))
      ds=$((dt % 60))
      dm=$(((dt / 60) % 60))
      dh=$((dt / 3600))
      printf '%d:%02d:%02d' $dh $dm $ds
  fi
}

SCRIPTPATH=$(dirname "${BASH_SOURCE[0]}")

#############
# Databases #
#############
filterRef=("$SCRIPTPATH/databases/homo_sapiens.fna" "$SCRIPTPATH/databases/phi174_seq_NC001422.fa")
alienseq="$SCRIPTPATH/databases/alienTrimmerPF8contaminants.fasta"
eggnogRef="$SCRIPTPATH/databases/eggnog_v3/fasta/eggnog_v3.faa"
eggnogRefBlast="$SCRIPTPATH/databases/eggnog_v3/blast/eggnog_v3.faa"
eggnogdescRef="$SCRIPTPATH/databases/eggnog_v3/info/description/"
eggnogmemRef="$SCRIPTPATH/databases/eggnog_v3/info/members/"
eggnogfunccatRef="$SCRIPTPATH/databases/eggnog_v3/info/funccat/"
eggnogfunccats="$SCRIPTPATH/databases/eggnog_v3/info/eggnogv3.funccats.txt"
cazyRefHmmer="$SCRIPTPATH/databases/cazy/cazy_2013_05_11/hmmer/dbCAN-fam-HMMs.txt"
famInfoRef="$SCRIPTPATH/databases/cazy/cazy_2013_05_11/info/FamInfo.txt"
# keggRef="/services/project/biodatabank/kegg/kegg_2014_04_20/fasta/kegg_genes.faa"
# keggRefBlast="/services/project/biodatabank/kegg/kegg_2014_04_20/blast/kegg_genes.faa"
# koRef="/services/project/biodatabank/kegg/kegg_2014_04_20/info/ko.txt"
# koGeneKoRef="/services/project/biodatabank/kegg/kegg_2014_04_20/info/genes_ko.list"
# koGeneUniprotRef="/services/project/biodatabank/kegg/kegg_2014_04_20/info/genes_uniprot.list"
# koUniprotRef="/services/project/biodatabank/kegg/kegg_2014_04_20/info/ko_uniprot.list"
# koGeneGeneidRef="/services/project/biodatabank/kegg/kegg_2014_04_20/info/genes_ncbi-geneid.list"
# koGeneGiRef="/services/project/biodatabank/kegg/kegg_2014_04_20/info/genes_ncbi-gi.list"
# koGenePathwayRef="/services/project/biodatabank/kegg/kegg_2014_04_20/info/genes_pathway.list"

gitaxidnucl="$SCRIPTPATH/databases/gi_taxid_nucl.dmp"
nodes="/local/databases/release/taxodb/nodes.dmp"
names="/local/databases/release/taxodb/names.dmp"
ntBlast="/local/databases/fasta/nt"
wgsBlast="/local/databases/fasta/genbank_wgsnuc"

minikrakendb="$SCRIPTPATH/databases/MiniKraken_DB/minikraken_20141208"

#######################
# Assembly Parameters #
#######################
NbMismatchMapping=1
NbProc=$(grep -c ^processor /proc/cpuinfo)
scaffoldLengthThreshold=500
geneLengthThreshold=60
# Prodigal filtering
filterGene=11
# Activity
all=1
assembly=0
tax_annotation=0
func_annotation=0
readTrim=1
maxTargetSeqs=10
evalueTaxAnnot="1e-5"
evalueFuncAnnot="1e-5"
numberBestannotation=1
spadesAssembly=0
wgs=0
metagenemark=0
kegg=0

############
# Programs #
############
# Mathieu
FilterFasta="$SCRIPTPATH/FilterFasta/FilterFasta.py"
extractimomi="$SCRIPTPATH/ExtractIMOMIAnnotation/ExtractIMOMIAnnotation.py"
extractMarker="$SCRIPTPATH/ExtractMarker/ExtractMarker.py"
trimReads="$SCRIPTPATH/TrimReads/trimReads.py"
# Amine
filterProdigal="$SCRIPTPATH/FilterProdigal/FilterProdigal.py"
extractEggNog="$SCRIPTPATH/ExtractEggNog/ExtractEggNog.py"
extractCazy="$SCRIPTPATH/ExtractCazy/ExtractCazy.py"
extractProteins="$SCRIPTPATH/ExtractProteins/ExtractProteins.py"
extractKegg="$SCRIPTPATH/ExtractKegg/ExtractKegg.py"
blastAnalyzer="$SCRIPTPATH/BlastAnalyzer/BlastAnalyzer.py"
extractNCBIDB="$SCRIPTPATH/ExtractNCBIDB/ExtractNCBIDB.py"
grabcataloguesequence="$SCRIPTPATH/grab_catalogue_sequence/grab_catalogue_sequence.py"
gettaxonomy="$SCRIPTPATH/get_taxonomy/get_taxonomy"
extractMetagenemark="$SCRIPTPATH/ExtractMetagenemark/ExtractMetagenemark.py"

# Blast
blastp="blastp" #"$SCRIPTPATH/ncbi-blast_2.2.30/blastp"
blastn="blastn" #"$SCRIPTPATH/ncbi-blast_2.2.30/blastn"

# Bowtie
bowtie2="bowtie2" #"$SCRIPTPATH/bowtie2-2.2.3/bowtie2"
bowtie2_build="bowtie2-build" #"$SCRIPTPATH/bowtie2-2.2.3/bowtie2-build"
# Prodigal
prodigal="prodigal" #$(which prodigal) #"$SCRIPTPATH/Prodigal-2.6.1/prodigal"
# Hmmer
hmmscan="hmmscan" #"$SCRIPTPATH/hmmer_3.1b1/hmmscan"
# AlienTrimmer
alientrimmer="AlienTrimmer" #"$SCRIPTPATH/AlienTrimmer_0.4.0/src/AlienTrimmer.jar"
fastqc="fastqc"
# FetchMG
fetchmg="$SCRIPTPATH/fetchMG/fetchMG.pl"
# kraken
kraken="$SCRIPTPATH/kraken_bin/kraken"
krakenreport="$SCRIPTPATH/kraken_bin/kraken-report"
# spades
spades="spades.py" #"$SCRIPTPATH/SPAdes-3.8.2-Linux/bin/spades.py"
# Metagenemark
gmhmmp="$SCRIPTPATH/MetaGeneMark/mgm/gmhmmp"
metagenemarkmod="$SCRIPTPATH/MetaGeneMark/mgm/MetaGeneMark_v1.mod"
# FetchMG
fetchmg="$SCRIPTPATH/fetchMG/fetchMG.pl"
# parallel
#parallel="$SCRIPTPATH/gnu_parallel/bin/parallel"
# Krona
krona="ktImportText" #"$SCRIPTPATH/KronaTools-2.6.1/bin/bin/ktImportText"
# Samtools
samtools="samtools"
# Extract samtools data
extractSamtools="$SCRIPTPATH/ExtractSamtools/ExtractSamtools.py"
# Barrnap
barrnap="barrnap"

########
# Main #
########

# Execute getopt on the arguments passed to this program, identified by the special character $@
PARSED_OPTIONS=$(getopt -n "$0"  -o hs:1:2:o:r:e:n:m:tl:f:p:g:q: --long "help,sample_name:,input1:,input2:,input_protein:,input_gene:,output:,NbProc:,NbMismatchMapping:,noTrim,scaffoldLengthThreshold:,geneLengthThreshold:,filterGene:,SamplePrefix:,metagenemark,spadesAssembly,evalueTaxAnnot:,evalueFuncAnnot:,numberBestannotation:,maxTargetSeqs:,referenceGenome:,wgs,assembly,tax_annotation,func_annotation,imomiRefBlast:,imomiIDtoGenome:,markerCOGRef:,markerCOGRefBlast:,eggnogRef:,eggnogRefBlast:,eggnogdescRef:,eggnogmemRef:,eggnogfunccats:,eggnogfunccatRef:,cazyRefHmmer:,famInfoRef:,kegg,keggRef:,keggRefBlast:,koRef:,koGeneKoRef:,koGeneUniprotRef:,koUniprotRef:,koGeneGeneidRef:,koGeneGiRef:,koGenePathwayRef:"  -- "$@")

#Check arguments
if [ $# -eq 0 ]
then
    display_help 0
fi

#Bad arguments, something has gone wrong with the getopt command.
if [ $? -ne 0 ];
then
    display_help 1
    exit 1
fi

# A little magic, necessary when using getopt.
eval set -- "$PARSED_OPTIONS"

# Get Cmd line arguments depending on options
while true;
do
  case "$1" in
    -h)
        display_help 0
        shift;;
    --help)
        display_help 1
        shift;;
    -s|--sample_name)
        sample_name=$2
        shift 2;;
    -1|--input1) 
        check_file $2
        input1=$2
        shift 2;;
    -2|--input2)
        check_file $2
        input2=$2
        shift 2;;
    -p|--input_protein)
        check_file $2
        input_protein=$2
        shift 2;;
    -g|--input_gene)
        check_file $2
        input_gene=$2
        shift 2;;
    -o|--output)
        resultDir=$2
        logDir=$resultDir/log/
        errorlogDir=$resultDir/error_log/
        check_dir $resultDir
        check_dir $logDir
        check_dir $errorlogDir
        shift 2;;
    -n|--NbProc)
        check_integer $2
        NbProc=$2
        shift 2;;
    -m|--NbMismatchMapping)
        check_integer $2
        NbMismatchMapping=$2
        shift 2;;
    -t|--noTrim)
        readTrim=0
        shift;;
    --scaffoldLengthThreshold)
        check_integer $2
        scaffoldLengthThreshold=$2
        shift 2;;
    -l|--geneLengthThreshold)
        check_integer $2
        geneLengthThreshold=$2
        shift 2;;
    -f|--filterGene)
        check_integer $2
        filterGene=$2
        shift 2;;
    -r|--SamplePrefix)
        SamplePrefix=$2
        shift 2;;
    -e|--insertSize)
        check_integer $2
        insertSize=$2
        shift 2;;
    --metagenemark)
        metagenemark=1
        shift;;
    --spadesAssembly)
        spadesAssembly=1
        shift;;
    --evalueTaxAnnot)
        evalueTaxAnnot=$2
        shift 2;;
    --evalueFuncAnnot)
        evalueFuncAnnot=$2
        shift 2;;
    --numberBestannotation)
        check_integer $2
        numberBestannotation=$2
        shift 2;;
    --maxTargetSeqs)
        check_integer $2
        maxTargetSeqs=$2
        shift 2;;
    --referenceGenome)
        check_file $2
        referenceGenome=$2
        shift 2;;
    --wgs)
        wgs=1
        shift 1;;
    --assembly)
        assembly=1
        shift;;
    --tax_annotation)
        tax_annotation=1
        shift;;
    --func_annotation)
        func_annotation=1
        shift;;
    --imomiRefBlast)
        imomiRefBlast=$2
        shift 2;;
    --imomiIDtoGenome)
        check_file $2
        imomiIDtoGenome=$2
        shift 2;;
    --markerCOGRef)
        check_file $2
        markerCOGRef=$2
        shift 2;;
    --markerCOGRefBlast)
        markerCOGRefBlast=$2
        shift 2;;
    --eggnogRef)
        check_file $2
        eggnogRef=$2
        shift 2;;
    --eggnogRefBlast)
       eggnogRefBlast=$2
       shift 2;;
    --eggnogdescRef)
        check_file $2
        eggnogdescRef=$2
        shift 2;;
    --eggnogmemRef)
        check_file $2
        eggnogmemRef=$2
        shift 2;;
    --eggnogfunccats)
        check_file $2
        eggnogfunccats=$2
        shift 2;;
    --eggnogfunccatRef)
        check_file $2
        eggnogfunccatRef=$2
        shift 2;;
    --cazyRefHmmer)
        check_file $2
        cazyRefHmmer=$2
        shift 2;;
    --famInfoRef)
        check_file $2
        famInfoRef=$2
        shift 2;;
    --kegg)
        kegg=1
        shift;;
    --keggRef)
        check_file $2
        keggRef=$2
        shift 2;;
    --keggRefBlast)
        check_file $2
        keggRefBlast=$2
        shift 2;;
    --koRef)
        check_file $2
        koRef=$2
        shift 2;;
    --koGeneKoRef)
        check_file $2
        koGeneKoRef=$2
        shift 2;;
    --koGeneUniprotRef)
        check_file $2
        koGeneUniprotRef=$2
        shift 2;;
    --koUniprotRef)
        check_file $2
        koUniprotRef=$2
        shift 2;;
    --koGeneGeneidRef)
        check_file $2
        koGeneGeneidRef=$2
        shift 2;;
    --koGeneGiRef)
        check_file $2
        koGeneGiRef=$2
        shift 2;;
    --koGenePathwayRef)
        check_file $2
        koGenePathwayRef=$2
        shift 2;;
    --)
      shift
      break;;
  esac
done

# Detect activity
if [ "$assembly" -eq "1" ] || [ "$tax_annotation" -eq "1" ] || [ "$func_annotation" -eq "1" ]
then
    all=0
fi

if [ "$resultDir" = "" ]
then
    error "Please indicate the output directory."
    exit 1
fi

# Check sample
if [ -f "$input1" ] && [ -f "$input2" ]
then
    if [ "$sample_name" = ""  ]
    then
      error "Please provide a sample name [-s]"
      exit 1
    fi
    filename=$(basename "$input1")
    extension=".${filename##*.}"
    if [ "$extension" != ".fastq" ] && [ "$extension" != ".fq" ]
    then
        error "The input file should be a fastq file."
        display_help
    fi
    SamplePath1=$(dirname $input1)
    SamplePath2=$(dirname $input2)
    if [ "$SamplePath1" != "$SamplePath2" ]
    then
        error "The sample must be in the same directory (for the moment)."
        exit 1
    else
        SamplePath="$SamplePath1"
    fi
    SampleName=$sample_name #"${filename%.*}"
    #SampleName=${filename::-8}
    # Set of gene and protein data
    input_gene=${resultDir}/${SampleName}_gene_${geneLengthThreshold}_filtered.fna
    input_protein=${resultDir}/${SampleName}_prot_filtered.faa
elif [ -f "$input_gene" ]
then
    filename=$(basename "$input_gene")
    SamplePath=$(dirname $input_gene)
    SampleName="${filename%.*}"
    if [ ! -f "$input_protein" ]
    then
        input_protein=${resultDir}/${SampleName}_prot_filtered.faa
    fi
elif [ -f "$input_protein" ]
then
    filename=$(basename "$input_protein")
    SamplePath=$(dirname $input_protein)
    SampleName="${filename%.*}"
else
    error "No input file !"
    exit 1
fi

# display parameters
display_parameters

# Start timer
say "Start analysis"
wall_time=$(timer)


if [ "$assembly" -eq "1" ] || [ "$all" -eq "1" ]
then
    # Filtering the reads
    let "essai=1";
    if [ ! -d "${resultDir}/filter_${#filterRef[@]}" ]
    then
            for db in ${filterRef[@]}
            do
                    let "num=$essai-1";
                    if [ ! -d "${resultDir}/filter_${essai}" ] && [ -d "${resultDir}/filter_${num}" ] || [ "$essai" -ne "1" ] 
                    then
                            say "Filter reads in $db"
                            start_time=$(timer)
                            mkdir ${resultDir}/filter_${essai}
                            # Next mapping
                            $bowtie2  -q -N $NbMismatchMapping -p $NbProc -x $db -1 ${resultDir}/filter_${num}/un-conc-mate.1  -2 ${resultDir}/filter_${num}/un-conc-mate.2 -S /dev/null --un-conc ${resultDir}/filter_${essai}/ -t --very-fast > ${logDir}/log_mapping_${SampleName}_${essai}.txt 2>&1
                            check_file ${resultDir}/filter_${essai}/un-conc-mate.1
                            # Remove old file
                            rm -rf ${resultDir}/filter_${num}
                            say "Elapsed time to filter reads in $db : $(timer $start_time)"
                    elif [ -f "$input1" ] && [ -f "$input2" ]  && [ "$essai" -eq "1" ] && [ ! -d "${resultDir}/filter_${essai}" ]
                    then
                            say "Filter reads in $db"
                            start_time=$(timer)
                            mkdir  ${resultDir}/filter_${essai}
                            # First mapping
                            $bowtie2 -q -N $NbMismatchMapping -p $NbProc -x $db  -1 $input1 -2 $input2 -S /dev/null --un-conc ${resultDir}/filter_${essai} -t  > ${logDir}/log_mapping_${SampleName}_${essai}.txt --very-fast 2>&1
                            check_file ${resultDir}/filter_${essai}/un-conc-mate.1
                            say "Elapsed time to filter reads in $db : $(timer $start_time)"
                    fi
                    let "essai=$essai+1";
            done
            # Change extension to fastq
            mv ${resultDir}/filter_${#filterRef[@]}/un-conc-mate.1 ${resultDir}/filter_${#filterRef[@]}/un-conc-mate_1.fastq
            mv ${resultDir}/filter_${#filterRef[@]}/un-conc-mate.2 ${resultDir}/filter_${#filterRef[@]}/un-conc-mate_2.fastq
    fi

    # Triming
    if [ "$readTrim" -eq "1" ]
    then
        filteredSample=$(readlink -f ${resultDir}/filter_alientrimmer/)
        if [ ! -d ${resultDir}/filter_alientrimmer ]
        then
            say "Triming reads with alientrimmer"
            start_time=$(timer)
            mkdir -p ${resultDir}/filter_alientrimmer
            $alientrimmer -if ${resultDir}/filter_${#filterRef[@]}/un-conc-mate_1.fastq -ir ${resultDir}/filter_${#filterRef[@]}/un-conc-mate_2.fastq -cf 1 -cr 2 -of ${resultDir}/filter_alientrimmer/un-conc-mate_1.fastq -or ${resultDir}/filter_alientrimmer/un-conc-mate_2.fastq -os ${resultDir}/filter_alientrimmer/un-conc-mate_sgl.fastq -c $alienseq  > ${logDir}/log_alientrimmer_${SampleName}.txt  2> ${errorlogDir}/error_log_alientrimmer_${SampleName}.txt
            check_file ${resultDir}/filter_alientrimmer/un-conc-mate_1.fastq
            check_log ${errorlogDir}/error_log_alientrimmer_${SampleName}.txt
        fi
    else
        filteredSample=$(readlink -f ${resultDir}/filter_${#filterRef[@]}/)
    fi

    # Fastqc
    if [ -f "${filteredSample}/un-conc-mate_1.fastq" ] && [ -f "${filteredSample}/un-conc-mate_2.fastq" ] && [ ! -f  "${filteredSample}/un-conc-mate_1_fastqc.html" ] && [ ! -f  "${filteredSample}/un-conc-mate_2_fastqc.html" ]
    then
        say "Quality control with Fastqc"
        start_time=$(timer)
        $fastqc ${filteredSample}/un-conc-mate_1.fastq ${filteredSample}/un-conc-mate_2.fastq --nogroup -q -t $NbProc 2> ${errorlogDir}/error_log_fastqc_${SampleName}.txt
        check_file ${filteredSample}/un-conc-mate_1_fastqc.html
        check_file ${filteredSample}/un-conc-mate_2_fastqc.html
        check_log ${errorlogDir}/error_log_fastqc_${SampleName}.txt
        say "Elapsed time with Fastqc: $(timer $start_time)"
    fi

    # Spades
    contigs=$(readlink -f "${resultDir}/${SampleName}_scaffolds.fasta")
    if [ -f "${resultDir}/filter_alientrimmer/un-conc-mate_1.fastq" ] && [ -f "${resultDir}/filter_alientrimmer/un-conc-mate_2.fastq" ] #&& [ ! -d "${resultDir}/spades/" ] && [ ! -f "$contigs" ]
    then
      say "Assembly insert with spades"
      start_time=$(timer)
      mkdir ${resultDir}/spades/
      $spades --meta -1 ${resultDir}/filter_alientrimmer/un-conc-mate_1.fastq -2 ${resultDir}/filter_alientrimmer/un-conc-mate_2.fastq  -t $NbProc -o ${resultDir}/spades/ > ${logDir}/log_spades_${SampleName}.txt  2> ${errorlogDir}/error_log_spades_${SampleName}.txt
      check_file ${resultDir}/spades/scaffolds.fasta
      ln -s $(readlink -f "${resultDir}/spades/scaffolds.fasta") $contigs
      say "Elapsed time to assembly with spades : $(timer $start_time)"
    fi
    
    # Predict Ribosomal RNA
    if [ -f "$contigs" ] && [ ! -f "${resultDir}/${SampleName}_rrna.txt" ]
    then
        say "Predict Ribosomal RNA"
        start_time=$(timer)
        $barrnap $contigs > ${resultDir}/${SampleName}_rrna.txt 2> ${logDir}/log_barrnap_${SampleName}.txt
        check_file ${resultDir}/${SampleName}_rrna.txt
        say "Elapsed time to predict Ribosomal RNA : $(timer $start_time)"
    fi

    # Predict genes
    if [ -f  "$contigs" ] && [ ! -f "${resultDir}/${SampleName}_gene.fna" ] && [ ! -f "${resultDir}/${SampleName}_prot.faa" ] && [ "$metagenemark" -eq "1" ]
    then
        say "Predict genes using Metagenemark"
        start_time=$(timer)
        cp $SCRIPTPATH/.gm_key $HOME/
        $gmhmmp -a -d  -m $metagenemarkmod $contigs  -o ${resultDir}/${SampleName}.metagenemark 2> ${logDir}/log_metagenemark_${SampleName}.txt
        $extractMetagenemark -i ${resultDir}/${SampleName}.metagenemark -a ${resultDir}/${SampleName}_prot.faa -d ${resultDir}/${SampleName}_gene.fna
        check_file ${resultDir}/${SampleName}_gene.fna
        check_file ${resultDir}/${SampleName}_prot.faa
        say "Elapsed time to predict genes using Metagenemark : $(timer $start_time)"
    elif [ -f  "$contigs" ] && [ ! -f "${resultDir}/${SampleName}_gene.fna" ] && [ ! -f "${resultDir}/${SampleName}_prot.faa" ]
    then
        say "Predict genes using Prodigal"
        start_time=$(timer)
        $prodigal -i $contigs -a ${resultDir}/${SampleName}_prot.faa -d ${resultDir}/${SampleName}_gene.fna -p meta > ${resultDir}/${SampleName}.prodigal 2> ${logDir}/log_prodigal_${SampleName}.txt
        check_file ${resultDir}/${SampleName}_gene.fna
        check_file ${resultDir}/${SampleName}_prot.faa
        say "Elapsed time to predict genes using Prodigal : $(timer $start_time)"
    fi

    # Filter prodigal gene over their length
    if [ -f "${resultDir}/${SampleName}_gene.fna" ] && [ ! -f "${resultDir}/${SampleName}_gene_${geneLengthThreshold}.fna" ]
    then
        say "Remove genes shorter than $geneLengthThreshold nt"
        start_time=$(timer)
        python $FilterFasta -f ${resultDir}/${SampleName}_gene.fna -t $geneLengthThreshold -o ${resultDir}/${SampleName}_gene_${geneLengthThreshold}.fna 2> ${errorlogDir}/error_log_filterfasta_gene_${SampleName}.txt
        check_log ${errorlogDir}/error_log_filterfasta_gene_${SampleName}.txt
        check_file ${resultDir}/${SampleName}_gene_${geneLengthThreshold}.fna
        say "Elapsed time to remove short genes : $(timer $start_time)"
    fi

    # Remove genes that lack both ends
    if [ -f "${resultDir}/${SampleName}_gene_$geneLengthThreshold.fna" ] && [ ! -f "$input_gene" ]
    then
        say "Remove genes that lack both ends"
        start_time=$(timer)
        python $filterProdigal -i ${resultDir}/${SampleName}_gene_$geneLengthThreshold.fna -o $input_gene -c $filterGene 2> ${errorlogDir}/error_log_filterprodial_gene_${SampleName}.txt
        check_log ${errorlogDir}/error_log_filterprodial_gene_${SampleName}.txt
        check_file $input_gene
        # Extract list of element of interest
        grep ">" $input_gene |cut -d" " -f1 |cut -d">" -f2  >  ${resultDir}/${SampleName}_interest_list.txt
        say "Elapsed time to remove genes that lack both ends : $(timer $start_time)"
    fi

    # Filter proteins based on selected genes
    if [ -f "$input_gene" ] && [ -f "${resultDir}/${SampleName}_prot.faa" ] && [ -f "${resultDir}/${SampleName}_interest_list.txt" ] && [ ! -f "$input_protein" ]
    then
        say "Select proteins for which the gene has been selected"
        start_time=$(timer)
        python $extractProteins -i ${resultDir}/${SampleName}_interest_list.txt -d ${resultDir}/${SampleName}_prot.faa -o $input_protein 2> ${errorlogDir}/error_log_ExtractProteins_${SampleName}.txt  > ${logDir}/log_grab_sequence_${SampleName}.txt
        check_log ${errorlogDir}/error_log_ExtractProteins_${SampleName}.txt
        check_file $input_protein
        say "Elapsed time to select proteins for which the gene has been selected : $(timer $start_time)"
    fi
fi

if [ "$tax_annotation" -eq "1" ] || [ "$all" -eq "1" ]
then
    # Read annotation with kraken
    if [ -f "${resultDir}/${SampleName}_alien.fastq" ]  && [ ! -f "${resultDir}/${SampleName}_kraken_annotation.txt" ]
    then
       say "Taxonomic annotation with kraken"
       start_time=$(timer)
       $kraken --db $minikrakendb --threads $NbProc --fastq-input --preload "${resultDir}/${SampleName}_alien.fastq" --output ${resultDir}/${SampleName}_kraken_result.txt 2> ${logDir}/log_kraken_${SampleName}.txt
       check_file ${resultDir}/${SampleName}_kraken_result.txt
       $krakenreport --db $minikrakendb ${resultDir}/${SampleName}_kraken_result.txt > ${resultDir}/${SampleName}_kraken_annotation.txt 2> ${errorlogDir}/error_log_kraken_report_${SampleName}.txt
       #check_log ${errorlogDir}/error_log_kraken_report_${SampleName}.txt
       check_file ${resultDir}/${SampleName}_kraken_annotation.txt
       say "Elapsed time to annotation with kraken :  $(timer $start_time)"
    fi

    # Check for essential genes : MOTU
    if [ -f $input_protein ] && [ ! -f ${resultDir}/${SampleName}_marker_cog.txt ] #&&  [ ! -f ${resultDir}/${SampleName}_prot_cog_essential.txt ]
    then
        say "Check for essential genes : MOTU"
        start_time=$(timer)
#         python $extractMarker -i $input_protein -d $markerCOGRef -db $markerCOGRefBlast -o ${resultDir}/${SampleName}_prot_cog_essential.txt  -n $NbProc -b $(dirname $blastp) -s ${resultDir}/${SampleName}_prot_cog_stat.txt 2> ${errorlogDir}/error_log_extract_marker_${SampleName}.txt
        mkdir -p ${resultDir}/fetch_result/
        $fetchmg -m extraction $input_protein -o ${resultDir}/fetch_result/ -d $input_gene -t $NbProc > ${logDir}/log_extract_marker_${SampleName}.txt #2> ${errorlogDir}/error_log_extract_marker_${SampleName}.txt
        wc -l ${resultDir}/fetch_result/temp/*.txt > ${resultDir}/${SampleName}_marker_cog.txt
        cat ${resultDir}/fetch_result/*.fna > ${resultDir}/${SampleName}_gene_marker.fna
        #check_log ${errorlogDir}/error_log_extract_marker_${SampleName}.txt
        #check_file ${resultDir}/${SampleName}_prot_cog_essential.txt
        say "Elapsed time for MOTU : $(timer $start_time)"
    fi

    # Taxonomic annotation of genes with blastn on nt
    if [ -f "$input_gene" ] && [ ! -f "${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_nt.txt" ]
    then
        say "Taxonomic annotation with blastn on nt"
        start_time=$(timer)
        # -use_index true -db_soft_mask 11
        #cat $input_gene | $parallel -j $NbProc --block 500K --recstart '>' --pipe  "$blastn -query - -db $ntBlast -outfmt 6 -evalue $evalueTaxAnnot -max_target_seqs $maxTargetSeqs -task megablast" > ${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_nt.txt 2> ${errorlogDir}/error_log_blastn_nt_${SampleName}.txt
        $blastn -query $input_gene -db $ntBlast -outfmt "6 qseqid sseqid qlen length mismatch gapopen qstart qend sstart send pident qcovs evalue bitscore" -evalue $evalueTaxAnnot -max_target_seqs $maxTargetSeqs -task megablast -out ${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_nt.txt -num_threads $NbProc 2> ${errorlogDir}/error_log_blastn_nt_${SampleName}.txt
        check_log ${errorlogDir}/error_log_blastn_nt_${SampleName}.txt
        check_file ${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_nt.txt
        say "Elapsed time for taxonomic annotation with blastn on nt : $(timer $start_time)"
    fi


    # Taxonomic annotation of genes with blastn on wgs
    if [ -f "$input_gene" ] && [ ! -f "${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_wgs.txt" ] && [ "$wgs" -eq "1" ]
    then
        say "Taxonomic annotation with blastn on wgs"
        start_time=$(timer)
        #cat $input_gene | $parallel -j $NbProc --block 500K --recstart '>' --pipe  "$blastn -query - -db $wgsBlast -outfmt 6 -evalue $evalueTaxAnnot -max_target_seqs $maxTargetSeqs -task megablast" >  ${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_wgs.txt 
        $blastn -query $input_gene -db $wgsBlast -outfmt "6 qseqid sseqid qlen length mismatch gapopen qstart qend sstart send pident qcovs evalue bitscore" -evalue $evalueTaxAnnot -max_target_seqs $maxTargetSeqs -task megablast -out ${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_wgs.txt -num_threads $NbProc 2> ${errorlogDir}/error_log_blastn_wgs_${SampleName}.txt
        check_log ${errorlogDir}/error_log_blastn_wgs_${SampleName}.txt
        check_file ${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_wgs.txt
        say "Elapsed time for taxonomic annotation with blastn on wgs : $(timer $start_time)"
    fi

     # Extract annotations on ncbi
    if [ -f "${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_nt.txt" ] && [ ! -f "${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_ncbi_genome_cov80.txt" ] && [ -f "$input_gene" ]
    then
        say "Analyze results on ncbi"
        start_time=$(timer)
        # ncbi annotation
        if [ -f ${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_wgs.txt ]
        then
            cat ${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_nt.txt ${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_wgs.txt > ${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_ncbi.txt
            ncbi_annotation="${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_ncbi.txt"
        else
            ncbi_annotation="${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_nt.txt"
        fi
        if [ ! -f "${resultDir}/${SampleName}_taxonomy.txt" ]
        then
            $gettaxonomy -i $ncbi_annotation -t $gitaxidnucl -n $names -d $nodes -o ${resultDir}/${SampleName}_taxonomy.txt 2> ${errorlogDir}/error_log_gettaxonomy_${SampleName}.txt
            check_file ${resultDir}/${SampleName}_taxonomy.txt
            check_log ${errorlogDir}/error_log_gettaxonomy_${SampleName}.txt
        fi
        if [ -f "${resultDir}/${SampleName}_taxonomy.txt" ]
        then
            python $extractNCBIDB -f $ncbi_annotation -g ${resultDir}/${SampleName}_taxonomy.txt -o ${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_ncbi_genome_cov80.txt -nb $numberBestannotation  -fc 80 2> ${errorlogDir}/error_log_extractncbidb_${SampleName}.txt
            check_log ${errorlogDir}/error_log_extractncbidb_${SampleName}.txt
            check_file ${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_ncbi_genome_cov80.txt
        fi
        say "Elapsed time to analyze results on ncbi : $(timer $start_time)"
    fi

     # Visualization of taxonomic annotation
    if [ -f "${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_ncbi_genome_cov80.txt" ] && [ ! -f "${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_ncbi_genome_cov80_krona.html" ]
    then
        say "Visualization of taxonomic annotation on nt with krona"
        start_time=$(timer)
        tail -n +2 ${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_ncbi_genome_cov80.txt |cut -f2-9 > ${resultDir}/tmp_annotation
        while read line
        do
            echo -e "1\t$line"
        done < ${resultDir}/tmp_annotation > ${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_ncbi_genome_cov80_krona.txt
        $krona ${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_ncbi_genome_cov80_krona.txt -o ${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_ncbi_genome_cov80_krona.html -n "Kingdom" -c
        rm -f "${resultDir}/tmp_annotation"
        check_file ${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_ncbi_genome_cov80_krona.html
        say "Elapsed time with krona: $(timer $start_time)"
    fi

    # Map reads against gene set
    if [ -f "$input_gene" ] && [ ! -f "${resultDir}/${SampleName}_vs_gene_${geneLengthThreshold}_filtered.sam" ] && [ -f "${resultDir}/filter_alientrimmer/un-conc-mate_1.fastq" ] && [ -f "${resultDir}/filter_alientrimmer/un-conc-mate_2.fastq" ]
    then
        say "Map reads against gene set with Bowtie"
        start_time=$(timer)
        $bowtie2_build $input_gene $input_gene > ${logDir}/log_bowtie2_build_${SampleName}.txt
        #cat ${resultDir}/filter_alientrimmer/un-conc-mate_1.fastq ${resultDir}/filter_alientrimmer/un-conc-mate_2.fastq  > ${resultDir}/tmp_fastq
        $bowtie2 --fast-local -p $NbProc -x $input_gene -1 ${resultDir}/filter_alientrimmer/un-conc-mate_1.fastq -2 ${resultDir}/filter_alientrimmer/un-conc-mate_2.fastq -S ${resultDir}/${SampleName}_vs_gene_${geneLengthThreshold}_filtered.sam > ${logDir}/log_bowtie2_${SampleName}.txt 2>&1
        check_file ${resultDir}/filter_${essai}/un-conc-mate.1
        #rm -f ${resultDir}/tmp_fastq
        check_file ${resultDir}/${SampleName}_vs_gene_${geneLengthThreshold}_filtered.sam
        say "Elapsed time with bowtie: $(timer $start_time)"
    fi

    # Get abundance
    if [ -f "${resultDir}/${SampleName}_vs_gene_${geneLengthThreshold}_filtered.sam" ] && [ ! -f "${resultDir}/${SampleName}_vs_gene_${geneLengthThreshold}_filtered.bam" ]
    then
        say "Check abundance with samtools"
        start_time=$(timer)
        $samtools view -bS ${resultDir}/${SampleName}_vs_gene_${geneLengthThreshold}_filtered.sam | samtools sort - ${resultDir}/${SampleName}_vs_gene_${geneLengthThreshold}_filtered
        $samtools index ${resultDir}/${SampleName}_vs_gene_${geneLengthThreshold}_filtered.bam
        $samtools idxstats ${resultDir}/${SampleName}_vs_gene_${geneLengthThreshold}_filtered.bam > ${resultDir}/${SampleName}_vs_gene_${geneLengthThreshold}_filtered_count.txt
        check_file ${resultDir}/${SampleName}_vs_gene_${geneLengthThreshold}_filtered_count.txt
        say "Elapsed time with bowtie: $(timer $start_time)"
    fi

    # Associate abundance and taxonomy
    if [ -f "${resultDir}/${SampleName}_vs_gene_${geneLengthThreshold}_filtered_count.txt" ] && [ -f "${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_ncbi_genome_cov80.txt" ] && [ ! -f "${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_ncbi_genome_cov80_abundance.txt" ]
    then
        say "Associate abundance and taxonomy"
        start_time=$(timer)
        python $extractSamtools -i ${resultDir}/${SampleName}_vs_gene_${geneLengthThreshold}_filtered_count.txt -a ${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_ncbi_genome_cov80.txt -k ${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_ncbi_genome_cov80_abundance.txt -o ${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_ncbi_genome_cov80_abundance_detailed.txt
        check_file ${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_ncbi_genome_cov80_abundance.txt
        say "Elapsed time with extractSamtools: $(timer $start_time)"
    fi

    # Visualization of the abundance and annotation
    if [ -f "${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_ncbi_genome_cov80_abundance.txt" ] && [ ! -f "${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_ncbi_genome_cov80_abundance.html" ]
    then
        say "Visualization of abundance/taxonomic annotation on nt with krona"
        start_time=$(timer)
        $krona ${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_ncbi_genome_cov80_abundance.txt -n SuperKingdom -c -o ${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_ncbi_genome_cov80_abundance.html
        check_file ${resultDir}/${SampleName}_gene_${geneLengthThreshold}_blastn_ncbi_genome_cov80_abundance.html
        say "Elapsed time with krona: $(timer $start_time)"
    fi
fi


## Annotation
if [ "$func_annotation" -eq "1" ] || [ "$all" -eq "1" ]
then
    # blastp on egg/nog
    if [ -f "$input_protein" ] && [ ! -f "${resultDir}/${SampleName}_prot_eggnog_annotation.txt" ]
    then
        say "Functional annotation with egg/nog"
        start_time=$(timer)
        #cat $input_protein | $parallel -j $NbProc --block 100K --recstart '>' --pipe  "$blastp -query - -db $eggnogRefBlast -outfmt 6 -evalue $evalueFuncAnnot -max_target_seqs $maxTargetSeqs" > ${resultDir}/${SampleName}_prot_eggnog_annotation.txt 2> ${errorlogDir}/error_log_blastp_eggnog_${SampleName}.txt
        $blastp -query $input_protein -db $eggnogRefBlast -outfmt "6 qseqid sseqid qlen length mismatch gapopen qstart qend sstart send pident qcovs evalue bitscore" -evalue $evalueFuncAnnot -max_target_seqs $maxTargetSeqs -out ${resultDir}/${SampleName}_prot_eggnog_annotation.txt -num_threads $NbProc 2> ${errorlogDir}/error_log_blastp_eggnog_${SampleName}.txt
        check_log ${errorlogDir}/error_log_blastp_eggnog_${SampleName}.txt
        #check_file  ${resultDir}/${SampleName}_prot_eggnog_annotation.txt
        say "Elapsed time to get functional annotation with egg/nog : $(timer $start_time)"
    fi

    # Extract eggnog annotations
    if [ -f  "${resultDir}/${SampleName}_prot_eggnog_annotation.txt" ] && [ ! -f "${resultDir}/${SampleName}_prot_eggnog_annotation_filtered.txt" ] && [ -f "$input_protein" ]
    then
        say "Extract egg/nog result"
        start_time=$(timer)
        # -c -p $(dirname $clustalo)
        python $extractEggNog -b ${resultDir}/${SampleName}_prot_eggnog_annotation.txt -m $eggnogmemRef -n $eggnogdescRef -o ${resultDir}/${SampleName}_prot_eggnog_annotation_filtered.txt -f $eggnogfunccatRef -a $eggnogfunccats -nb $numberBestannotation 2> ${errorlogDir}/error_log_extractegg_${SampleName}.txt
        check_log ${errorlogDir}/error_log_extractegg_${SampleName}.txt
        check_file ${resultDir}/${SampleName}_prot_eggnog_annotation_filtered.txt
        say "Elapsed time to extract egg/nog result : $(timer $start_time)"
    fi

    # hmmer on cazy
    if [ -f "$input_protein" ] && [ ! -f "${resultDir}/${SampleName}_prot_cazy_annotation.txt" ]
    then
        say "Functional annotation with cazy"
        start_time=$(timer)
        $hmmscan -E $evalueFuncAnnot --cpu $NbProc --tblout ${resultDir}/${SampleName}_prot_cazy_annotation.txt $cazyRefHmmer $input_protein > ${logDir}/log_cazy_hmm_${SampleName}.txt 2> ${errorlogDir}/error_log_cazy_hmm_${SampleName}.txt
        check_log ${errorlogDir}/error_log_cazy_hmm_${SampleName}.txt
        #check_file ${resultDir}/${SampleName}_prot_cazy_annotation.txt
        say "Elapsed time to get functional annotation with cazy : $(timer $start_time)"
    fi

    # Extract cazy annotations
    if [ -f "${resultDir}/${SampleName}_prot_cazy_annotation.txt" ] && [ ! -f "${resultDir}/${SampleName}_prot_cazy_annotation_filtered.txt" ]
    then
        say "Extract cazy result"
        start_time=$(timer)
        python $extractCazy -i ${resultDir}/${SampleName}_prot_cazy_annotation.txt -f $famInfoRef -o ${resultDir}/${SampleName}_prot_cazy_annotation_filtered.txt 2> ${errorlogDir}/error_log_extractcazy_${SampleName}.txt
        #check_log ${errorlogDir}/error_log_extractcazy_${SampleName}.txt
        #check_file ${resultDir}/${SampleName}_prot_cazy_annotation_filtered.txt
        say "Elapsed time to extract cazy result : $(timer $start_time)"
    fi

#     # blastp on kegg
#     if [ -f $input_protein ] && [ ! -f ${resultDir}/${SampleName}_prot_kegg_annotation.txt ] && [ "$kegg" -eq "1" ]
#     then
#         say "Functional annotation with kegg"
#         start_time=$(timer)
#         cat $input_protein | $parallel -j $NbProc --block 100K --recstart '>' --pipe  "$blastp -query - -db $keggRefBlast -outfmt 6 -evalue $evalueFuncAnnot -max_target_seqs $maxTargetSeqs" > ${resultDir}/${SampleName}_prot_kegg_annotation.txt 2> ${errorlogDir}/error_log_blastp_kegg_${SampleName}.txt
#         # There is often log in case of change of some amino-acid
#         #check_log ${errorlogDir}/error_log_blastp_kegg_${SampleName}.txt
#         #check_file  ${resultDir}/${SampleName}_prot_kegg_annotation.txt
#         say "Elapsed time to get functional annotation with kegg : $(timer $start_time)"
#     fi
# 
#     # Extract kegg annotation
#     if [  -f ${resultDir}/${SampleName}_prot_kegg_annotation.txt ] && [ ! -f ${resultDir}/${SampleName}_prot_kegg_annotation_filtered.txt ] && [ -f $input_protein ]
#     then
#         say "Extract kegg result"
#         start_time=$(timer)
#         # -p $(dirname $clustalo) -c
#         python $extractKegg -q $input_protein -d $keggRef -b ${resultDir}/${SampleName}_prot_kegg_annotation.txt -k $koRef -g $koGeneKoRef -ku $koUniprotRef  -u $koGeneUniprotRef -gid $koGeneGeneidRef -gi $koGeneGiRef -a $koGenePathwayRef -sa ${SampleName}_prot -r ${resultDir}/ -t $NbProc -nb $numberBestannotation 2> ${errorlogDir}/error_log_extractkegg_${SampleName}.txt > ${resultDir}/log_extractkegg_${SampleName}.txt
#         #check_log ${errorlogDir}/error_log_extractkegg_${SampleName}.txt
#         #check_file ${resultDir}/${SampleName}_prot_kegg_annotation_filtered.txt
#         say "Elapsed time to extract kegg result : $(timer $start_time)"
#     fi
fi
say "Assembly is done. Elapsed time: $(timer $wall_time)"
