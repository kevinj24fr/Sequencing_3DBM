#!/bin/bash


#####################################################################################################
# Script to prepare your computer to generate count matrices out of fastq files                     #
# Author: Kevin Joseph (kevin.joseph@uniklinik-freiburg.de)                                         #
# Lab: 3DBM, Neurosurgery, University Hospital of Freiburg                                          #
#####################################################################################################

echo "Setting up your computer to process illumina sequencing data"
echo "This script is for MacOS based systems"
# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Check if Homebrew is installed
if ! command_exists brew; then
    echo "Homebrew is not installed. Installing Homebrew..."
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    echo 'export PATH="/opt/homebrew/bin:$PATH"' >> ~/.bash_profile
    source ~/.bash_profile
else
    echo "Homebrew is already installed."
fi

# Install Java for Trimmomatic (required dependency)
if ! command_exists java; then
    echo "Java is not installed. Installing Java..."
    brew install openjdk
    echo 'export PATH="/usr/local/opt/openjdk/bin:$PATH"' >> ~/.bash_profile
    source ~/.bash_profile
else
    echo "Java is already installed."
fi

# Array of prerequisites with alternative installation methods
PREREQUISITES=("fastqc" "trimmomatic" "hisat2" "samtools" "subread")

# Loop through each prerequisite and install if not present
for TOOL in "${PREREQUISITES[@]}"; do
    case $TOOL in
        "trimmomatic")
            if ! command_exists trimmomatic; then
                echo "Trimmomatic is not installed. Installing Trimmomatic..."
                brew install brewsci/bio/trimmomatic || {
                    echo "Failed to install Trimmomatic via Homebrew. Downloading manually..."
                    wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
                    unzip Trimmomatic-0.39.zip -d ~/trimmomatic
                    echo 'export PATH="$HOME/trimmomatic:$PATH"' >> ~/.bash_profile
                    source ~/.bash_profile
                }
            else
                echo "Trimmomatic is already installed."
            fi
            ;;
        "hisat2")
            if ! command_exists hisat2; then
                echo "HISAT2 is not installed. Installing HISAT2..."
                brew install brewsci/bio/hisat2 || {
                    echo "Failed to install HISAT2 via Homebrew. Downloading manually..."
                    wget https://cloud.biohpc.swmed.edu/index.php/s/hisat2-2.2.1/download -O hisat2-2.2.1-linux_x86_64.zip
                    unzip hisat2-2.2.1-linux_x86_64.zip -d ~/hisat2
                    echo 'export PATH="$HOME/hisat2:$PATH"' >> ~/.bash_profile
                    source ~/.bash_profile
                }
            else
                echo "HISAT2 is already installed."
            fi
            ;;
        "subread")
            if ! command_exists featureCounts; then
                echo "Subread (featureCounts) is not installed. Installing Subread..."
                brew install brewsci/bio/subread || {
                    echo "Failed to install Subread via Homebrew. Downloading manually..."
                    wget https://sourceforge.net/projects/subread/files/subread-2.0.3/subread-2.0.3-source.tar.gz/download -O subread-2.0.3.tar.gz
                    tar -xzvf subread-2.0.3.tar.gz
                    cd subread-2.0.3 && make -f Makefile.Linux
                    echo 'export PATH="$HOME/subread-2.0.3/bin:$PATH"' >> ~/.bash_profile
                    source ~/.bash_profile
                }
            else
                echo "Subread (featureCounts) is already installed."
            fi
            ;;
        *)
            if ! command_exists "$TOOL"; then
                echo "$TOOL is not installed. Installing $TOOL..."
                brew install "$TOOL"
            else
                echo "$TOOL is already installed."
            fi
            ;;
    esac
done

echo "All prerequisites are installed."
