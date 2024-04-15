# IFDlong #

## About ##

IFDlong is a bioinformatics pipeline that can perform long-read RNA-seq annotation at isoform levels, fusion detection, as well as fusion and isoform quantification.


## Dependencies ##

Some of the dependent tools can be automatically downloaded and installed by run `bash install.sh` from the sources below:

- [minimap2](https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2)
- [bedtools2](https://github.com/arq5x/bedtools2/releases/download/v2.29.0/bedtools-2.29.0.tar.gz)
- [samtools](http://sourceforge.net/projects/samtools/files/samtools/1.9/samtools-1.9.tar.bz2)

However, **gcc (>= 8.2.0)** and **R (>=4.0.0)** must be manually installed by users. In addition, the R pacakges **rlist, parallel, stringr, dplyr, seqRFLP, BiocManager,rtracklayer, Biostrings** are required and recommended to be installed before running the pipeline, though they can be installed automatically.

The paths to these will be automatically put into your tools.path file after you successfully run `bash install.sh`. 

## Usage ##

### Installing

It will automatically install the dependent tools and put their paths to the `tools.path` file.

```bash
git clone https://github.com/wenjiaking/IFDlong.git
cd IFDlong
bash install.sh
```

### Run for Human or Mouse Data

Our tools contains all the required reference files for human and mouse that have been generated from GRCh38 (Genome Research Consortium human build 38) and GRCm38 (Genome Reference Consortium Mouse Build 38) respectively under the `refData` folder. So you don't need to prepare the genome references unless you want to refer to other genome assembly.


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
-g, --ghc         Human (-g hg38) or mouse data (-g mm10), hg38 by default
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
-g, --ghc         Human (-g hg38) or mouse data (-g mm10), hg38 by default
-t, --bufferLen   The buffer length for novel isoform identification, 9 by default
-a, --anchorLen   The anthor length for fusion filtering, 10 by default
-c, --ncores      How many cores are assigned to run the pipeline in parallel. Use single core by default
```


### Prepare Reference Files for Species Other than Human and Mouse

If the long-read RNA-seq data is other species, you should prepare the required references before runing `IFDlong.sh` or `IFDlong_minimap2.sh`. The following command can set up the reference database given reference genome sequence in FASTA file (e.g., genome.fa), the reference genome annotation in GTF file (e.g., genes.gtf), the pseudogene list in RDS format, as well as the gene family root names in TXT file (containing two columns: "Family.name" and "Common.root.gene.symbol").

The reference genome sequence and genome annotation are required, while the reference pseudogene list and root names of gene family are highly suggested to be provided, otherwise reads containing pseudo genes and fusion reads consisting of genes from the same family cannot be filtered out by runing IFDlong pipeline. 


```bash
bash refDataSetup.sh [OPTIONS]
```

Running with default settings:
```bash
bash refDataSetup.sh -r genome_GTF_file -m genome_FASTA_file -p pseudo_RDS_file -f roots_TXT_file -g species -c 1
```

Required options:
```
-h, --help        Check the usage
-r, --refGTF      The reference genes.gtf file from GRC
-m, --refGenome   The reference genome.fa file from GRC
-p, --pseudoList  The reference list of pseudo genes in RDS format
-f, --roots       The gene family root names in TXT format
-g, --ghc         Speficy the species other than Human (hg38) and mouse (mm10)
-c, --ncores      How many cores are assigned to run in parallel. Use single core by default

```


