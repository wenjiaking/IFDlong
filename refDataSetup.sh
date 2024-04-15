#!/bin/bash 
############################################################
# Help                                                     #
############################################################
Help()
{  
   # Display Help
   echo ""
   echo "Options:"
   echo "  -h, --help        Check the usage. Please note that reference pseudogene list and root names of gene family are highly suggested to be provided, otherwise fusions consisting of pseudogenes or genes from the same family cannot be filtered out. "
   echo "  -r, --refGTF      The reference genes.gtf file from GRC."
   echo "  -m, --refGenome   The reference genome.fa file from GRC."
   echo "  -p, --pseudoList  The reference list of pseudo genes in RDS format."
   echo "  -f, --roots       The gene family root names in TXT format."
   echo "  -g, --ghc         Speficy the species. The reference database has already built for Human (hg38) and mouse (mm10)."
   echo "  -c, --ncores      How many cores are assigned to generate reference files in parallel. Use single core by default"

   exit 1
}

############################################################
# Process the input options. Add options as needed.        #
############################################################

# Default value of argument

ncores=1
codeBase="$(dirname "$0")"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in

    -h|--help)
      Help
      ;;

    -r|--refGTF)
      refGTFFile="$2"
      shift
      shift
      ;;

    -m|--refGenome)
      refGenome="$2"
      shift
      shift
      ;;

    -p|--pseudoList)
      refpseudo="$2"
      shift
      shift
      ;;

    -f|--roots)
      refroots="$2"
      shift
      shift
      ;;

    -g|--ghc)
      ghc="$2"
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
if [ -z "$refGTFFile" ]; then
  echo "Missing required option to input the reference GTF file from GHC : -r|--refGTF"
  Help
  
elif [[ -z "$ghc" ]]; then
  echo "Missing required option to specify the species (other than human and mouse): -g|--ghc"
  Help

elif [[ -z "$refGenome" ]]; then
  echo "Missing required option to input the reference genome file from GHC: -m|--refGenome"
  Help

fi



############################################################
############################################################
# Main program                                             #
############################################################

echo Begin to build reference $(date '+%Y-%m-%d %H:%M:%S')

source $codeBase/tools.path

refDir=$codeBase/refData/$ghc
mkdir -p $refDir

cp $refGTFFile $refDir
base_name=$(basename ${refGTFFile})
mv $refDir/$base_name $refDir/genes.gtf

cp $refGenome $refDir
base_name=$(basename ${refGenome})
mv $refDir/$base_name $refDir/genome.fa

if [ -z "$refpseudo" ]; then
  echo "Missiong optional reference of pseudo gene list. When run IFDlong, reads containing pseudo genes cannot be fitered out."
else
  cp $refpseudo $refDir
  base_name=$(basename ${refpseudo})
  mv $refDir/$base_name $refDir/pseudogenes.rds
fi

if [ -z "$refroots" ]; then
  echo "Missiong optional reference of gene family root names. When run IFDlong, fusion reads consisting of genes from the same family cannot be filtered out."
else
  cp $refroots $refDir
  base_name=$(basename ${refroots})
  mv $refDir/$base_name $refDir/rootName.txt
fi


refFile=$refDir/allexon_NO.bed
refFiletol=$refDir/gene_range_tol500.bed
refCDS=$refDir/CDS_all.bed 

refCDSfa=$refDir/CDS_all.fa
refIsoformfa=$refDir/isoformAll_ref.fa
refAAFile=$refDir/isoformAA.txt


echo Begin Reference Generating $(date '+%Y-%m-%d %H:%M:%S')
$Rscript $refGen $refGTFFile $refFile $refFiletol $refCDS


$bedtools getfasta -fi $refGenome -bed $refCDS -split -s -name -fo $refCDSfa

echo Begin Animo Acid Reference Predicting $(date '+%Y-%m-%d %H:%M:%S')
$Rscript $refAAGen $refCDSfa $refIsoformfa $refAAFile $ncores

echo Reference has been built $(date '+%Y-%m-%d %H:%M:%S')

echo Used memory $(free |grep Mem|awk '{print $3}') 
echo Used memory percentage $(free |grep Mem|awk '{print $3/$2 * 100.0}')


