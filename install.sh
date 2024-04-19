#!/bin/bash

## This script will install the tools required for the IFDlong pipeline.
## It will fetched each tool from the web and placed into the tools/ subdirectory.
## Paths to all installed tools can be found in the file tools.path at the
## end of execution of this script. These paths can be changed if a different
## version of software is required. Note that R (>=4.0.0) must be installed manually
##
## Last Modified: Apr. 2024 by Wenjia Wang

# make the bash scripts executable
chmod +x *.sh

#installation methods
mkdir -p tools/bin 
cd tools 

#a list of which programs need to be installed
commands="samtools minimap2"

function minimap2_install {
   wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2
   #wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17.tar.bz2
   tar -xvf minimap2-2.17_x64-linux.tar.bz2 ; rm minimap2-2.17_x64-linux.tar.bz2
#   make -C minimap2-2.17
   cp minimap2-2.17_x64-linux/minimap2 bin/
}

function bedtools_install {
    wget https://github.com/arq5x/bedtools2/releases/download/v2.29.0/bedtools-2.29.0.tar.gz
    tar -zxvf bedtools-2.29.0.tar.gz; rm bedtools-2.29.0.tar.gz
    cd $PWD/bedtools2
    make
    cd ..
    #cp bedtools2/bin/* bin/
    cp bedtools2/bin/bedtools bin/
}


function samtools_install {
   wget --no-check-certificate http://sourceforge.net/projects/samtools/files/samtools/1.9/samtools-1.9.tar.bz2
   tar -jxvf samtools-1.9.tar.bz2
   rm samtools-1.9.tar.bz2
   make prefix=$PWD install -C samtools-1.9/
}

#function bowtie2_install {
#    wget --no-check-certificate http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.4/bowtie2-2.4.4-linux-x86_64.zip
#    unzip bowtie2-2.4.4-linux-x86_64.zip ; rm bowtie2-2.4.4-linux-x86_64.zip
#    ln -s $PWD/bowtie2-2.4.4-linux-x86_64/bowtie2 $PWD/bin
#    ln -s $PWD/bowtie2-2.4.4-linux-x86_64/bowtie2-build $PWD/bin
#}

#Check if the version of gcc is >= 8.2.0
gcc_version=`gcc -dumpversion`
gcc_check=`echo -e "$gcc_version\n8.2.0" | sort -n | tail -n1`
if [[ $gcc_chek = "8.2.0" ]] 
then 
   echo "Your version of gcc is $gcc_version."
   echo "gcc must be >= 8.2.0 to install IFDlong. Exiting..."
   exit 1
fi

echo "gcc check passed"


#Check if the version of R is >= 4.0.0
#R_version=`R --version`
#R_check=`echo -e "${R_version:10:5}\n4.0.0" | sort -n | tail -n1`
#if [[ $R_chek = "4.0.0" ]] 
#then 
#   echo "Your version of R is ${R_version:10:5}."
#   echo "R must be >= 4.0.0 to install IFDlong. Exiting..."
#   exit 1
#fi
#echo "R check passed"

#finally check that R is install, need to install packages: parallel, rlist, stringr, 
# packages for generating reference files: rtracklayer, Biostrings, dypyr, seqRFLP
R_path=`which Rscript 2>/dev/null`
if [ -z $R_path ] ; then
    echo "R not found!"
    echo "Please go to http://www.r-project.org/ and follow the installation instructions."
    echo "R version >= 4.0.0."
    echo "Require the dependent packages: parallel, rlist, stringr, rtracklayer, Biostrings, dplyr, seqRFLP."
    exit 1
fi
echo "Rscript=\"$R_path\"" > ../tools.path

bedtools_path=`which bedtools 2>/dev/null`
if [ -z $bedtools_path ] ; then

    bedtools_path=`which $PWD/bin/bedtools 2>/dev/null`
    if [ -z $bedtools_path ] ; then
        echo "bedtools not found!"
        echo "Please go to https://github.com/arq5x/bedtools2 to download version 2.29 following the instructions."
        echo "Copy the bedtools under your downloaded bin folder to ./tools/bin."
    fi
    
fi
echo "bedtools=\"$bedtools_path\"" >> ../tools.path


echo "// Path to tools used by the IFDlong pipeline" >> ../tools.path

for c in $commands ; do 
    c_path=`which $c 2>/dev/null`
    if [ -z $c_path ] ; then 
        c_path=`which $PWD/bin/$c 2>/dev/null`
        if [ -z $c_path ] ; then 
            echo "$c not found, fetching it"
            ${c}_install
            c_path=`which $PWD/bin/$c 2>/dev/null`
        fi
    fi
    echo "$c=\"$c_path\"" >> ../tools.path
done

#loop through commands to check they are all installed
echo "Checking that all required tools were installed:"
Final_message="All commands installed successfully!"
source ../tools.path
for t in $commands ; do
    #c_path=`which $PWD/bin/$c 2>/dev/null`
    c_path=`which $c 2>/dev/null`
    if [ -z $c_path ] ; then
	echo -n "WARNING: $t could not be found!!!! " 
	echo "You will need to download and install $t manually, then add its path to tools.path"
	Final_message="WARNING: One or more command did not install successfully. See warning messages above. \
                       You will need to correct this before running IFDlong."
    else 
        echo "$t looks like it has been installed"
    fi
done
echo "**********************************************************"
echo $Final_message

cd ..
mkdir -p scripts
cd scripts

stages="reportFuns reportFuns_parallel EXONuncover report quant refGen refAAGen"
echo "// Path to scripts used by the IFDlong pipeline" >> ../tools.path
for s in $stages ; do 
    s_path=$PWD/$s.R
    echo "$s=\"$s_path\"" >> ../tools.path
done






