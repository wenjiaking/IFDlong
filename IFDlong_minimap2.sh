#!/bin/bash (or #!/usr/bin/env bash)
############################################################
# Help                                                     #
############################################################
Help()
{  
   # Display Help
   echo ""
   echo "Options:"
   echo "  -h, --help        Check the usage."
   echo "  -o, --outDir      The directory to save the output."
   echo "  -n, --name        The sample name."
   echo "  -i, --inFile      Input FASTQ file of long read RNAseq data."
   echo "  -g, --ghc         Human (-g hg38) or mouse data (-g mm10), hg38 by default."
   echo "  -t, --bufferLen   The buffer length for novel isoform identification, 9 by default."
   echo "  -a, --anchorLen   The anthor length for fusion filtering, 10 by default."
   echo "  -c, --ncores      How many cores are assigned to run the pipeline in parallel. Use single core by default"

   exit 1
}

############################################################
# Process the input options. Add options as needed.        #
############################################################

# Default value of argument

ghc="hg38"
Aligner="MINIMAP2"
bufferLen=9
anchorLen=10
ncores=1
codeBase="$(dirname "$0")"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in

    -h|--help)
      Help
      ;;

    -o|--outDir)
      mainPath="$2"
      shift
      shift
      ;;

    -n|--name)
      sample="$2"
      shift
      shift
      ;;

    -i|--inFile)
      inFile="$2"
      shift
      shift
      ;;
      
    -g|--ghc)
      ghc="$2"
      shift
      shift
      ;;

    -t|--bufferLen)
      bufferLen="$2"
      shift
      shift
      ;;

    -a|--anchorLen)
      anchorLen="$2"
      shift
      shift
      ;;

    -c|--ncores)
      ncores="$2"
      shift
      shift
      ;;
      
    *)
      echo "Invalid option: $1"
      Help
      ;;
      
  esac
done

# Check if the required options are provided
if [ -z "$name" ]; then
  echo "Missing required option to specify sample name: -n|--name"
  Help
  
elif [[ -z "$outDir" ]]; then
  echo "Missing required option to set output directory: -o|--outDir"
  Help

elif [[ -z "$inFile" ]]; then
  echo "Missing required option to input FASTQ file: -i|--inFile"
  Help

fi


############################################################
############################################################
# Main program                                             #
############################################################

echo Pipeline Begin $(date '+%Y-%m-%d %H:%M:%S')

source $codeBase/tools.path

# e.g., ghc="hg38"
refFile=$codeBase/refData/$ghc/allexon_NO.bed
refFiletol=$codeBase/refData/$ghc/gene_range_tol500.bed
genomeDir=$codeBase/refData/$ghc/genome.fa 


outPath=$mainPath/$sample/$Aligner
InFile=$inFile
BAMfile=$outPath/$sample.bam

buffer=$bufferLen
anchorLen=$anchorLen
refAAFile=$codeBase/refData/$ghc/isoformAA.txt
refPseudoFile=$codeBase/refData/$ghc/refData/pseudogenes.rds
refRootFile=$codeBase/refData/$ghc/rootName.txt
refGTFFile=$codeBase/refData/$ghc/genes.gtf
hmmatchFile=$codeBase/refData/$ghc/match.rds
#Note that if the ghc is neither hg38 nor mm10, the pseudo genes and family genes may not provided
# only use anchor length to do fusion filtering

if [ $ncores -gt 1 ]; then
  reportFuns=$reportFuns_parallel
  echo "Run the pipeline in parallel with $ncores cores."
fi


mkdir -p $outPath

echo Begin Alignment $(date '+%Y-%m-%d %H:%M:%S')
$minimap2 -ax splice $genomeDir $InFile | samtools view -Sb | samtools sort -o $BAMfile
echo Alignment Done $(date '+%Y-%m-%d %H:%M:%S')
echo samtools index
$samtools index $BAMfile 

echo filter out unmapped and multiple alignment for $sample
$samtools view -b -F 4 $BAMfile >$outPath/$sample"_mapped.bam"
$samtools view -b -F 256 $outPath/$sample"_mapped.bam" >$outPath/$sample"_mapped_woSecond.bam"
rm $outPath/$sample"_mapped.bam"
echo Finish Filtering for $sample

echo Generate BED intersect for $sample
$bedtools bamtobed -i $outPath/$sample"_mapped_woSecond.bam" -split -cigar> $outPath/$sample"_mapped_woSecond.bed"
$bedtools intersect -a $outPath/$sample"_mapped_woSecond.bed" -b $refFile  -f 0.90  -wao >$outPath/$sample"_mapped_woSecond_intersectS.bed"
echo Intersecting Finish for $sample

echo Begin EXON-uncovered blocks generating $(date '+%Y-%m-%d %H:%M:%S')
$Rscript $EXONuncover $mainPath $sample $Aligner $refFile $reportFuns
echo EXON-uncovered blocks generated in bed file Done!

echo Begin gene range covering
$bedtools intersect -a $outPath/$sample"_mapped_woSecond_intersectS_EXONuncover.bed" -b $refFiletol  -f 0.90  -wao >$outPath/$sampleName"_mapped_woSecond_geneTol500intersectS.bed"
echo Gene range covering done!

echo Begin isoform annotation
$Rscript $report $mainPath $sample $Aligner $buffer $anchorLen $refFile $refAAFile $refPseudoFile $refRootFile $hmmatchFile $ghc $reportFuns $ncores
echo Isoform annotation done!$(date '+%Y-%m-%d %H:%M:%S')

echo Begin isoform quantification $(date '+%Y-%m-%d %H:%M:%S')
$Rscript $quant $mainPath $sample $Aligner $buffer $anchorLen $refGTFFile $reportFuns $ncores

echo Isoform quantification done!$(date '+%Y-%m-%d %H:%M:%S')

echo Pipeline End $(date '+%Y-%m-%d %H:%M:%S')

echo Used memory $(free |grep Mem|awk '{print $3}') 
echo Used memory percentage $(free |grep Mem|awk '{print $3/$2 * 100.0}')
