# IFDlong #

## About ##

IFDlong is a bioinformatics pipeline that can perform long-read RNA-seq annotation at isoform levels, fusion detection, as well as fusion and isoform quantification.


## Installation ##
### Dependencies ###
R (≥ 4.0.0) along with a compatible version of gcc.  
The following tools are required for running IFDlong: minimap2 (≥ 2.17), bedtools (≥ 2.29), and samtools (≥ 1.9)  
If these tools are not already installed, they will be installed automatically during the install.sh step.
Required R packages include: rlist, parallel, stringr, dplyr, seqRFLP, BiocManager, rtracklayer, and Biostrings.
These packages will also be installed automatically during processing.

### Installation pipeline ###
Methods1: Install maunally  
It will automatically install the dependent tools and put their paths to the `tools.path` file.  
```bash
git clone https://github.com/wenjiaking/IFDlong.git
cd IFDlong
bash install.sh
```

Methods2: Install via conda
```
git clone https://github.com/wenjiaking/IFDlong.git
conda env create -f IFDlong.yaml
```

## Installation Reference Database ##
Our tools include built-in reference datasets for human (hg38) and mouse (mm10). If you would like to use a custom reference for another species, please refer to the optional commands.
Download the genome.fa and genes.gtf files, and place them in the appropriate refDB/species/ folder.  
##### Human
++link
```

```

##### Mouse
```
wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz 
```


## Usage ##
### Quick Star
To verify your installation and run a example.
```bash
bash IFDlong.sh -o output_directory -n sample_name -i input_file -l "self_align" -g "hg38" -t 9 -a 10 -c 1
```

#### Demo
+++

### Run IFDlong Pipeline
Our tool accepts both FASTQ and BAM files as input. If using FASTQ, it will be aligned using Minimap2. If using BAM, please ensure the files are already aligned and indexed.

Running with default settings:
```bash
bash IFDlong.sh -o output_directory -n sample_name -b input_BAM_file -l "self_align" -g "hg38" -t 9 -a 10 -c 1
```

Required options:
```
-h, --help        Check the usage.
-o, --outDir      The directory to save the output
-n, --name        The sample name
-i, --inFile      Input file <sampleID.fq>|<sampleID.bam>
-l, --aligner     The aligner used to generated the bam file. Set to be self_align if missing.
-g, --ghc         Human (hg38), mouse (mm10) or other self-defined species (the same value as -g in refDataSetup.sh), hg38 by default
-t, --bufferLen   The buffer length for novel isoform identification, 9 by default
-a, --anchorLen   The anthor length for fusion filtering, 10 by default
-c, --ncores      How many cores are assigned to run the pipeline in parallel. Use single core by default

```



### Reference Database preparation

1. Download the USCS Genome Research Consortium from [illumina iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html). For human and mouse, unzip and save the genome sequence file "genome.fa" and genome annotation file "genes.gtf" to `refData/hg38` and `refData/mm10` respectively. For other species, create a subfolder under `refData`, e.g., `refData/self_species`, and save the genome sequence and annotation files to the directory.

2. Our tools contains the formatted reference files for human and mouse that have been generated based on GRCh38 (Genome Research Consortium human build 38) and GRCm38 (Genome Reference Consortium Mouse Build 38) respectively under the `refData` folder. So you don't need to run `refDataSetup.sh`. However, for other specises, set up the reference folder by `refDataSetup.sh` as below before runing `IFDlong.sh` or `IFDlong_minimap2.sh`.


The following command can set up the reference database for the species other than human and mouse given reference genome sequence in FASTA file and the reference genome annotation in GTF file downloaded in the first step. You are also recommended to provide the pseudogene list in RDS format and the gene family root names in TXT file (containing two columns: "Family.name" and "Common.root.gene.symbol") in order to filter out the reads annotated to pseudo genes and fusions consisting of genes from the same family. All the files must be saved under the species folder (e.g., self_species) under `refData`.

The reference genome sequence and genome annotation are required, while the reference pseudogene list and root names of gene family are highly suggested to be provided, otherwise reads containing pseudo genes and fusion reads consisting of genes from the same family cannot be filtered out by runing IFDlong pipeline. 


```bash
bash refDataSetup.sh [OPTIONS]
```

Running with default settings:
```bash
bash refDataSetup.sh -r path_to_genome.fa -m path_to_genes.gtf -p path_to_pseudo.rds -f path_to_roots.txt -g self_species -c 1
```

Required options:
```
-h, --help        Check the usage
-r, --refGTF      Full path to the reference genes.gtf file from GRC
-m, --refGenome   Full path to the reference genome.fa file from GRC
-p, --pseudoList  Full path to the reference list of pseudo genes in RDS format
-f, --roots       Full path to the gene family root names in TXT format
-g, --ghc         Speficy the species other than Human (hg38) and mouse (mm10)
-c, --ncores      How many cores are assigned to run in parallel. Use single core by default

```

## Citation ##
The study describing the IFDlong method can be found in: 


