#!/bin/bash
set -euo pipefail

echo "Command: $0 $@"

### Functions ###
usage() {
    echo "IFDlong perform isoform-level annotation of long-read RNA-seq data, detect gene fusions, and quantify both fusions and isoforms."
    echo -e "\tDemo working directory is example folder"
    echo "Usage:"
    echo "Options:"
    echo "  -h, --help        Check the usage."
    echo "  -o, --outDir      The directory to save the output."
    echo "  -n, --name        The sample name."
    echo "  -i, --inFile      Input FASTQ file of long read RNAseq data OR BAM file after alignment and indexing."
    echo "  -l, --aligner     The aligner used to generate the BAM file; set it to self_align if the input format is BAM."
    echo "  -g, --ghc         Human (hg38), mouse (mm10) or other self-defined species (the same value as -g in refDataSetup.sh), hg38 by default"
    echo "  -t, --bufferLen   The buffer length for novel isoform identification, 9 by default."
    echo "  -a, --anchorLen   The anthor length for fusion filtering, 10 by default."
    echo "  -c, --ncores      How many cores are assigned to run the pipeline in parallel. Use 4 core by default"

    echo "    Questions or issues? Contact: Silvia (shl96[at].pitt.edu)"
    echo "    Modified date: 08 Aug. 2025"
}


symlink_path () {
    SD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
    echo $SD
}


align () {
    echo "### The process the input ###"
    if [[ -z "$inFile" ]]; then
      echo "Usage: $0 <input_file.fq|input_file.bam>"
      exit 1
    fi
    
    # Get file extension
    ext="${inFile##*.}"
    
    filename=$(basename "$inFile")
    if [[ "$filename" == *.fq.gz ]]; then
        ext="fq.gz"
    elif [[ "$filename" == *.fastq.gz ]]; then
        ext="fastq.gz"
    else
        ext="${inFile##*.}"
    fi
    
    
    if [[ "$ext" == "fastq" || "$ext" == "fq" || "$ext" == "fa" || "$ext" == "fasta" || "$ext" == "fq.gz" || "$ext" == "fastq.gz" ]]; then 
      echo "Detected FASTQ file: $inFile"
      
      echo Begin Alignment by minimap2 $(date '+%Y-%m-%d %H:%M:%S')
      Aligner="MINIMAP2"
      outPath=$mainPath/$Aligner
      BAMfile=$outPath/$sample.bam

      echo $outPath $Aligner

      mkdir -p "$outPath"
    
      $minimap2 -ax splice "$genomeDir" "$inFile" -t "$ncores" | samtools view -Sb | samtools sort -o $BAMfile
      echo Alignment Done $(date '+%Y-%m-%d %H:%M:%S')
      echo samtools index
      $samtools index $BAMfile 
      echo samtools index Done $(date '+%Y-%m-%d %H:%M:%S')
    
    elif [[ "$ext" == "bam" ]]; then
      echo "Detected BAM file: $inFile"
      Aligner="self_align"
      outPath=$mainPath/$Aligner

      mkdir -p "$outPath"

      BAMfile=$inFile
    
    else
      echo "Unsupported file format: $ext"
      exit 1
    fi

}


filter () {
    echo filter out unmapped and multiple alignment for $sample $(date '+%Y-%m-%d %H:%M:%S')
    $samtools view -b -F 4 "$BAMfile" | \
    $samtools view -b -F 256 - > "$outPath/${sample}_mapped_woSecond.bam"
    # rm "$outPath/$sample"_mapped.bam"
    echo Finish Filtering for $sample $(date '+%Y-%m-%d %H:%M:%S')

    echo Generate BED intersect for $sample $(date '+%Y-%m-%d %H:%M:%S')
    $bedtools bamtobed -i "$outPath/${sample}_mapped_woSecond.bam" -split -cigar > $outPath/${sample}"_mapped_woSecond.bed"
    $bedtools intersect -a $outPath/${sample}"_mapped_woSecond.bed" -b "$refFile" -f 0.90 -wao > "$outPath/${sample}_mapped_woSecond_intersectS.bed"
    echo Intersecting Finish for $sample $(date '+%Y-%m-%d %H:%M:%S')

}


blocks () {
    echo Begin EXON-uncovered blocks generating $(date '+%Y-%m-%d %H:%M:%S')
    $Rscript $EXONuncover $mainPath $sample $Aligner $refFile
    echo EXON-uncovered blocks generated in bed file Done!
    
    echo Begin gene range covering
    $bedtools intersect -a $outPath/$sample"_mapped_woSecond_intersectS_EXONuncover.bed" -b $refFiletol  -f 0.90  -wao > $outPath/${sample}"_mapped_woSecond_geneTol500intersectS.bed"
    echo Gene range covering done!
}


anno () {
    echo Begin isoform annotation $(date '+%Y-%m-%d %H:%M:%S')
    $Rscript $report $mainPath $sample $Aligner $buffer $anchorLen $refFile $refAAFile $refPseudoFile $refRootFile $hmmatchFile $ghc $ncores
    echo Isoform annotation done! $(date '+%Y-%m-%d %H:%M:%S')
}


quant () {
    echo Begin isoform quantification $(date '+%Y-%m-%d %H:%M:%S')
    $Rscript $quant $mainPath $sample $Aligner $buffer $anchorLen $refGTFFile $ncores
    echo Isoform quantification done!$(date '+%Y-%m-%d %H:%M:%S')
}





### Initialization ###
# Check if no arguments were provided
if [ $# -eq 0 ]; then
    usage
    exit 0

fi

######## default parameter
#Aligner="self_align"
ghc="hg38"
ncores=1
mainPath="output"
codeBase="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

echo $mainPath 
echo $codeBase

### Argument Parsing ###
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      print_help_and_exit
      ;;

    -o|--outDir)
      mainPath="$2"
      shift 2
      ;;

    -n|--name)
      sample="$2"
      shift 2
      ;;

    -i|--inFile)
      inFile="$2"
      shift 2
      ;;

    -l|--aligner)
      Aligner="$2"
      shift 2
      ;;

    -g|--ghc)
      ghc="$2"
      shift 2
      ;;

    -t|--bufferLen)
      bufferLen="$2"
      shift 2
      ;;

    -a|--anchorLen)
      anchorLen="$2"
      shift 2
      ;;

    -c|--ncores)
      ncores="$2"
      shift 2
      ;;

    *)
      echo "Invalid option: $1"
      usage
      ;;
  esac
done

### Required Argument Validation ###
missing_arg=false

if [ -z "$sample" ] && [ -n "$inFile" ]; then
  sample=$(basename "$inFile")
  sample="${sample%%.*}"
  echo "Sample name not provided. Using '$sample' from input file."
fi

if [ -z "$mainPath" ]; then
  echo "Sample name not provided. Using '$mainPath' as default."
  missing_arg=true
fi

if [ -z "$inFile" ]; then
  echo "Missing required option: -i|--inFile (input file)"
  missing_arg=true
fi

if [ "$missing_arg" = true ]; then
  echo
  print_help_and_exit
fi

echo "Initialization complete."


##### check the dependencies
Rscript=$(which Rscript)
samtools=$(which samtools)
bedtools=$(which bedtools)
minimap2=$(which minimap2)

check_tool() {
    tool_name=$1
    tool_path=$(which "$tool_name" 2>/dev/null)

    if [ -z "$tool_path" ]; then
        echo "Error: $tool_name not found in PATH."
        return 1
    else
        echo "$tool_name found at: $tool_path"
        eval "$tool_name=$tool_path"
    fi
}

# List of required tools
check_tool "Rscript" || exit 1
check_tool "samtools" || exit 1
check_tool "bedtools" || exit 1
check_tool "minimap2" || exit 1



### Main ###
echo "Pipeline Begin $(date '+%Y-%m-%d %H:%M:%S')"


######## ref file path
refFile=$codeBase/refData/$ghc/allexon_NO.bed
refFiletol=$codeBase/refData/$ghc/gene_range_tol500.bed
genomeDir=$codeBase/refData/$ghc/genome.fa


buffer="${bufferLen:-9}"
anchorLen="${anchorLen:-10}"
refAAFile=$codeBase/refData/$ghc/isoformAA.txt
refPseudoFile=$codeBase/refData/$ghc/refData/pseudogenes.rds
refRootFile=$codeBase/refData/$ghc/rootName.txt
refGTFFile=$codeBase/refData/$ghc/genes.gtf
hmmatchFile=$codeBase/refData/$ghc/hg_mm_match.rds

## Path to scripts used by the IFDlong pipeline
echo $codeBase
EXONuncover="${codeBase}/scripts/EXONuncover.R"
report="${codeBase}/scripts/Reportcpp.r"
quant="${codeBase}/scripts/quant.R"


symlink_path
align
filter
blocks
anno
quant

echo "Completed!!!!!"

echo LOG END $(date '+%Y-%m-%d %H:%M:%S')

echo Used memory $(free |grep Mem|awk '{print $3}') 
echo Used memory percentage $(free |grep Mem|awk '{print $3/$2 * 100.0}')


exit 0



