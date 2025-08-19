#!/bin/bash


# make the bash scripts executable
chmod +x *.sh

# Define install directory
INSTALL_DIR="$(pwd)/tools/"
mkdir -p "$INSTALL_DIR"


check_and_install() {
    TOOL_NAME=$1
    GITHUB_URL=$2
    EXECUTABLE_NAME=$3

    echo "Checking for $TOOL_NAME"

    if command -v "$EXECUTABLE_NAME" >/dev/null 2>&1; then
        echo "$TOOL_NAME found at: $(command -v $EXECUTABLE_NAME)"
        echo "$TOOL_NAME=\"$(command -v "$TOOL_NAME")\"" >> ./tools.path
    else
        echo "$TOOL_NAME not found. Installing in $INSTALL_DIR"

        cd "$INSTALL_DIR" || exit 1

        case "$TOOL_NAME" in
            samtools)
                # Install samtools 
                wget --no-check-certificate http://sourceforge.net/projects/samtools/files/samtools/1.17/samtools-1.17.tar.bz2
   				tar -jxvf samtools-1.17.tar.bz2
   				cd samtools-1.17 || exit 1
                ./configure --prefix="$INSTALL_DIR" && make && make install
                export PATH="${INSTALL_DIR}/samtools-1.17:$PATH"
                rm samtools-*.tar.bz2
                ;;
            minimap2)
                git clone https://github.com/lh3/minimap2.git
                cd minimap2 || exit 1
                make
                export PATH="${INSTALL_DIR}/minimap2:$PATH"
                ;;
            bedtools)
                git clone https://github.com/arq5x/bedtools2.git
                cd bedtools2 || exit 1
                make
                export PATH="${INSTALL_DIR}/bedtools2/bin:$PATH"
                ;;
        esac

        cd "$INSTALL_DIR" || exit 1

        if command -v "$EXECUTABLE_NAME" >/dev/null 2>&1; then
            echo "$TOOL_NAME successfully installed at: $(command -v $EXECUTABLE_NAME)"
            echo "$TOOL_NAME=\"$(command -v "$TOOL_NAME")\"" >> ./tools.path
        else
            echo "Failed to install $TOOL_NAME"
        fi
    fi
}

# Run the check 
check_and_install "samtools" "https://github.com/samtools/samtools" "samtools"
check_and_install "minimap2" "https://github.com/lh3/minimap2" "minimap2"
check_and_install "bedtools" "https://github.com/arq5x/bedtools2" "bedtools"


#Check if the version of gcc is >= 12.2.0
gcc_version=`gcc -dumpversion`
gcc_check=`echo -e "$gcc_version\n12.2.0" | sort -n | tail -n1`
if [[ $gcc_chek = "12.2.0" ]] 
then 
   echo "Your version of gcc is $gcc_version."
   echo "gcc must be >= 12.2.0 to install IFDlong. Exiting..."
   exit 1
fi

echo "gcc check passed"

R_path=`which Rscript 2>/dev/null`
if [ -z $R_path ] ; then
    echo "R not found!"
    echo "Please go to http://www.r-project.org/ and follow the installation instructions."
    echo "R version >= 4.4.0."
    echo "Require the dependent packages: parallel, rlist, stringr, rtracklayer, Biostrings, dplyr, seqRFLP."
    exit 1
fi
echo "Rscript=\"$R_path\"" >> ./tools.path


