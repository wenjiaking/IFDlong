# IFDlong #

## About ##

IFDlong is a bioinformatics pipeline that can perform long-read RNA-seq annotation at isoform levels, fusion detection, as well as fusion and isoform quantification.


## Dependencies ##

The required tools for IFDlong are: gcc (>= 8.2.0), R (>=4.0.0), minimap2 (>=2.17), bedtools (>=2.29), samtools (>=1.9). 

Please make sure that **gcc (>= 8.2.0)**, **R (>=4.0.0)** and **bedtools (>=2.29)** already exist in your environment by `which` command. For example, check `which bedtools`, and if `bedtools` cannot be found, install [bedtools2](https://github.com/arq5x/bedtools2) following the instrunctions, and then copy the `bedtools` under your downloaded `bin` folder to `./tools/bin`. Below is the reference command to help you manually download and install the `bedtools`:

```bash
cd ./tools
wget https://github.com/arq5x/bedtools2/releases/download/v2.29.0/bedtools-2.29.0.tar.gz
tar -zxvf bedtools-2.29.0.tar.gz; rm bedtools-2.29.0.tar.gz
cd bedtools2
make
cd ..
cp bedtools2/bin/bedtools bin/
```

The other dependent tools can be automatically checked and installed by run `bash install.sh` from the sources below:

- [minimap2](https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2)
- [samtools](http://sourceforge.net/projects/samtools/files/samtools/1.9/samtools-1.9.tar.bz2)

In addition, the R pacakges **rlist, parallel, stringr, dplyr, seqRFLP, BiocManager,rtracklayer, Biostrings** are required and recommended to be installed before running the pipeline, though they can be installed automatically.

The paths to these will be automatically put into your tools.path file after you successfully run `bash install.sh`. 

## Usage ##

### Installing

It will automatically install the dependent tools and put their paths to the `tools.path` file.

```bash
git clone https://github.com/wenjiaking/IFDlong.git
cd IFDlong
bash install.sh
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


### Run IFDlong Pipeline

If your input is BAM file after alignment and indexing, then input the BAM file and run the following bash command:

```bash
bash IFDlong.sh [OPTIONS]
```

Running with default settings:
```bash
bash IFDlong.sh -o output_directory -n sample_name -b input_BAM_file -l "self_align" -g "hg38" -t 9 -a 10 -c 1
```

Required options:
```
-h, --help        Check the usage.
-o, --outDir      The directory to save the output
-n, --name        The sample name
-b, --bamFile     Input BAM file after alignment and indexing
-l, --aligner     The aligner used to generated the bam file. Set to be self_align if missing
-g, --ghc         Human (hg38), mouse (mm10) or other self-defined species (the same value as -g in refDataSetup.sh), hg38 by default
-t, --bufferLen   The buffer length for novel isoform identification, 9 by default
-a, --anchorLen   The anthor length for fusion filtering, 10 by default
-c, --ncores      How many cores are assigned to run the pipeline in parallel. Use single core by default

```

`minimap2` is embeded in the IFDlong tool. So if your input is FASTQ file, then input the FASTQ file and run the following bash command:

```bash
bash IFDlong_minimap2.sh [OPTIONS]
```

Running with default settings:
```bash
bash IFDlong_minimap2.sh -o output_directory -n sample_name -i input_FASTQ_file -g "hg38" -t 9 -a 10 -c 1
```

Required options:
```
-h, --help        Check the usage
-o, --outDir      The directory to save the output
-n, --name        The sample name
-i, --inFile      Input FASTQ file of long read RNAseq data
-g, --ghc         Human (hg38), mouse (mm10) or other self-defined species (the same value as -g in refDataSetup.sh), hg38 by default
-t, --bufferLen   The buffer length for novel isoform identification, 9 by default
-a, --anchorLen   The anthor length for fusion filtering, 10 by default
-c, --ncores      How many cores are assigned to run the pipeline in parallel. Use single core by default
```

