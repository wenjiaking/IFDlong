#!/bin/bash

echo "Command: $0 $@"


### Functions ###
usage() {
    echo "Refrence dataset setting up."
    echo -e "\tDownload the genome.fa and genes.gtf files for hg38 (human) or mm10 (mouse)"
    echo "Usage:"
    echo "Options:"
    echo "  -h, --help        Show this help message"
    echo "  -g, --ghc         Specify genome: hg38 (human) or mm10 (mouse). Default: hg38"
    echo "For questions or issues, contact: Silvia (shl96[at].pitt.edu)"
    echo "Modified: 08 Aug 2025"
}


######## default parameter
codeBase="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
ghc="hg38"  # default genome

### Argument Parsing ###
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      usage
      exit 0
      ;;
    -g|--ghc)
      ghc="$2"
      shift 2
      ;;
    *)
      echo "Invalid option: $1"
      usage
      exit 1
      ;;
  esac
done




#### check the reference genome and gene 
# URLs for downloading genome and GTF (can be updated)
if [ "$ghc" = "hg38" ]; then
    GENOME_URL="http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
    GTF_URL="http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz"
elif [ "$ghc" = "mm10" ]; then
    GENOME_URL="http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz"
    GTF_URL="http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.refGene.gtf.gz"
else
    echo "Unsupported genome: $ghc"
    exit 1
fi

mkdir -p "$codeBase/refData/$ghc"

FA="$codeBase/refData/$ghc/genome.fa"
GTF="$codeBase/refData/$ghc/genes.gtf"

# Download if missing
if [ ! -f "$FA" ]; then
    echo "Downloading $ghc FASTA..."
    wget -O "$codeBase/refData/$ghc/genome.all.fa.gz" "$GENOME_URL"
    gunzip -f "$codeBase/refData/$ghc/genome.all.fa.gz"
fi

if [ ! -f "$GTF" ]; then
    echo "Downloading $ghc GTF..."
    wget -O "$codeBase/refData/$ghc/genes.all.gtf.gz" "$GTF_URL"
    gunzip -f "$codeBase/refData/$ghc/genes.all.gtf.gz"
fi

### Keep only main chromosomes
echo "Filtering main chromosomes for $ghc..."

MAIN_FA="$codeBase/refData/$ghc/genome.fa"
MAIN_GTF="$codeBase/refData/$ghc/genes.gtf"

# Filter FASTA
awk '/^>/ {p=($0 ~ /^>chr([1-9]|1[0-9]|2[0-2]|X|Y)$/)} p' "$codeBase/refData/$ghc/genome.all.fa" > "$MAIN_FA"
rm "$codeBase/refData/$ghc/genome.all.fa"

# Filter GTF
awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ {print}' "$codeBase/refData/$ghc/genes.all.gtf" > "$MAIN_GTF"
rm "$codeBase/refData/$ghc/genes.all.gtf"

echo "Done! Main chromosomes saved:"
echo "FASTA: $MAIN_FA"
echo "GTF  : $MAIN_GTF"

### Check required files
files=(
  "$codeBase/refData/$ghc/genome.fa"
  "$codeBase/refData/$ghc/genes.gtf"
  "$codeBase/refData/$ghc/isoformAA.txt"
  "$codeBase/refData/$ghc/pseudogenes.rds"
  "$codeBase/refData/$ghc/rootName.txt"
  "$codeBase/refData/$ghc/hg_mm_match.rds"
  "$codeBase/refData/$ghc/allexon_NO.bed"
  "$codeBase/refData/$ghc/gene_range_tol500.bed"
)

all_exist=true
for f in "${files[@]}"; do
    if [ ! -f "$f" ]; then
        echo "Missing: $f"
        all_exist=false
    fi
done

if [ "$all_exist" = true ]; then
    echo "All files exist."
else
    echo "Some files are missing. Please download from GitHub."
fi
