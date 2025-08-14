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
conda env create -n IFDlong -f IFDlong.yaml
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
### Quick Start
To verify your installation and run a example.
```bash
bash IFDlong.sh -o output_directory -n sample_name -i input_file -l "self_align" -g "hg38" -t 9 -a 10 -c 1
```

#### Demo
```
bash IFDlong.sh -o out -n example -i example/demo.fq.gz -g "hg38" -t 9 -a 10 -c 1
```

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
-i, --inFile      Input file <sampleID.fq.gz>|<sampleID.bam>
-l, --aligner     The aligner used to generated the bam file. Set to be self_align if missing.
-g, --ghc         Human (hg38), mouse (mm10) or other self-defined species (the same value as -g in refDataSetup.sh), hg38 by default
-t, --bufferLen   The buffer length for novel isoform identification, 9 by default
-a, --anchorLen   The anthor length for fusion filtering, 10 by default
-c, --ncores      How many cores are assigned to run the pipeline in parallel. Use single core by default

```

## Citation ##
The study describing the IFDlong method can be found in: 


