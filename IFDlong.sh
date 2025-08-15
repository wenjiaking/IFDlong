#!/bin/bash
set -uo pipefail

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
sample=""
inFile=""
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
if [ -z "$inFile" ]; then
  echo "Missing required option: -i|--inFile (input file)"
  missing_arg=true
fi

if [ -z "${sample:-}" ]; then
  sample=$(basename "$inFile")
  sample="${sample%%.*}"
  echo "Sample name not provided. Using '$sample' from input file."
fi

# Ensure output directory is set
if [ -z "${mainPath:-}" ]; then
  mainPath="output"
  echo "Output directory not provided. Using default: '$mainPath'."
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

#### check the reference genome and gene 
# URLs for downloading genome and GTF (can be updated)
declare -A GENOME_URL
declare -A GTF_URL

GENOME_URL[hg38]="http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
GTF_URL[hg38]="http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.gtf.gz"

GENOME_URL[mm10]="http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz"
GTF_URL[mm10]="http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.gtf.gz"

FA=$(ls "$codeBase/refData/$ghc/"*.fa | head -n 1)
GTF=$(ls "$codeBase/refData/$ghc/"*.gtf | head -n 1)

# -------------------------
# Download if missing
# -------------------------
if [ ! -f "$FA" ]; then
    echo "Downloading $ghc FASTA..."
    wget -O "$codeBase/refData/$ghc/genome.all.fa.gz" "${GENOME_URL[$ghc]}"
    gunzip -f "$codeBase/refData/$ghc/genome.all.fa.gz"
fi

if [ ! -f "$GTF" ]; then
    echo "Downloading $GENOME GTF..."
    wget -O "$codeBase/refData/$ghc/genes.all.gtf.gz" "${GTF_URL[$ghc]}"
    gunzip -f "$codeBase/refData/$ghc/genes.all.gtf.gz"
fi


#### Keep only main chromosomes
echo "Filtering main chromosomes for $GENOME..."

MAIN_FA="$codeBase/refData/$ghc/genome.fa"
MAIN_GTF="$codeBase/refData/$ghc/genes.gtf"

# Filter FASTA
awk '/^>/ {p=($0 ~ /^>chr([1-9]|1[0-9]|2[0-2]|X|Y)$/)} p' $codeBase/refData/$ghc/genome.all.fa > "$MAIN_FA"
rm $codeBase/refData/$ghc/genome.all.fa

# Filter GTF
awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ {print}' $codeBase/refData/$ghc/genes.all.gtf > "$MAIN_GTF"
rm $codeBase/refData/$ghc/genes.all.gtf

echo "Done! Main chromosomes saved:"
echo "FASTA: $MAIN_FA"
echo "GTF  : $MAIN_GTF"


### Main ###
echo "Pipeline Begin $(date '+%Y-%m-%d %H:%M:%S')"


######## ref file path
refFile=$codeBase/refData/$ghc/allexon_NO.bed
refFiletol=$codeBase/refData/$ghc/gene_range_tol500.bed
genomeDir=$(ls "$codeBase/refData/$ghc/"*.fa | head -n 1)
refGTFFile=$(ls "$codeBase/refData/$ghc/"*.gtf | head -n 1)

buffer="${bufferLen:-9}"
anchorLen="${anchorLen:-10}"
refAAFile=$codeBase/refData/$ghc/isoformAA.txt
refPseudoFile=$codeBase/refData/$ghc/refData/pseudogenes.rds
refRootFile=$codeBase/refData/$ghc/rootName.txt
hmmatchFile=$codeBase/refData/$ghc/hg_mm_match.rds

## Path to scripts used by the IFDlong pipeline
echo $codeBase
EXONuncover="${codeBase}/scripts/EXONuncover.R"
report="${codeBase}/scripts/Reportcpp.r"
speedup="${codeBase}/scripts/speedup.cpp"
quant="${codeBase}/scripts/quant.R"


symlink_path
align
filter
blocks
anno
quant

if [ -f "$outPath/${sample}_mapped_woSecond_intersectS_buffer9bp_Isof_quant.csv" ]; then
    cp "$outPath/${sample}_mapped_woSecond_intersectS_buffer9bp_Isof_quant.csv" "$mainPath/${sample}_Isof_quant.csv"
else
    echo "Warning: Isof quant results not found. If you encounter any issues, feel free to report them on GitHub."
fi

# Copy Fusion quant file if it exists
if [ -f "$outPath/${sample}_mapped_woSecond_intersectS_buffer9bp_Fusion_quant_anchor10bp.csv" ]; then
    cp "$outPath/${sample}_mapped_woSecond_intersectS_buffer9bp_Fusion_quant_anchor10bp.csv" "$mainPath/${sample}_Fusion_quant.csv"
else
    echo "Warning: Fusion quant results not found. If you encounter any issues, feel free to report them on GitHub."
fi




echo "Completed!!!!!"

echo LOG END $(date '+%Y-%m-%d %H:%M:%S')

echo Used memory $(free |grep Mem|awk '{print $3}') 
echo Used memory percentage $(free |grep Mem|awk '{print $3/$2 * 100.0}')


exit 0



