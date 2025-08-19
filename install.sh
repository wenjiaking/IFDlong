#!/bin/bash


# make the bash scripts executable
chmod +x *.sh

# Define install directory
INSTALL_DIR="$(pwd)/tools"
mkdir -p "$INSTALL_DIR"


TOOLS_PATH_FILE="./tools.path"
> "$TOOLS_PATH_FILE" 


check_and_install() {
    TOOL_NAME=$1
    EXECUTABLE_NAME=$2

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
                echo "PATH=\"$INSTALL_DIR/samtools-1.17:\$PATH\"" >> "$TOOLS_PATH_FILE"
                rm "$INSTALL_DIR/samtools-1.17.tar.bz2"
                ;;
            minimap2)
                git clone https://github.com/lh3/minimap2.git
                cd minimap2 || exit 1
                make
                echo "PATH=\"$INSTALL_DIR/minimap2:\$PATH\"" >> "$TOOLS_PATH_FILE"
                ;;
            bedtools)
                git clone https://github.com/arq5x/bedtools2.git
                cd bedtools2 || exit 1
                make
                echo "PATH=\"$INSTALL_DIR/bedtools2/bin:\$PATH\"" >> "$TOOLS_PATH_FILE"
                ;;
        esac

        cd "$INSTALL_DIR" || exit 1

        if command -v "$EXECUTABLE_NAME" >/dev/null 2>&1; then
            echo "$TOOL_NAME successfully installed at: $(command -v $EXECUTABLE_NAME)"
            echo "$EXECUTABLE_NAME=\"$(command -v "$EXECUTABLE_NAME")\"" >> "$TOOLS_PATH_FILE"
        else
            echo "Failed to install $TOOL_NAME"
        fi
    fi
}

# Run the checks
check_and_install "samtools" "samtools"
check_and_install "minimap2" "minimap2"
check_and_install "bedtools" "bedtools"


### check gcc
gcc_version=$(gcc -dumpversion)
required_version="12.2.0"
if [[ $(printf '%s\n%s' "$gcc_version" "$required_version" | sort -V | head -n1) != "$required_version" ]]; then
    echo "Your version of gcc is $gcc_version."
    echo "gcc must be >= $required_version to install IFDlong. Exiting..."
    exit 1
fi
echo "gcc check passed"

### check Rscript
R_path=$(which Rscript 2>/dev/null)
if [ -z "$R_path" ] ; then
    echo "R not found!"
    echo "Please install R version >= 4.4.0 and required packages: parallel, rlist, stringr, rtracklayer, Biostrings, dplyr."
    exit 1
fi
echo "Rscript=\"$R_path\"" >> "$TOOLS_PATH_FILE"

echo "Installation finished. Run 'source ./tools.path' to use installed tools."


