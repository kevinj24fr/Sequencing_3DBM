#!/bin/bash


#####################################################################################################
# Script to prepare your computer to generate count matrices out of fastq files, Nanopore Edition   #
# Author: Kevin Joseph (kevin.joseph@uniklinik-freiburg.de)                                         #
# Lab: 3DBM, Neurosurgery, University Hospital of Freiburg                                          #
#####################################################################################################

echo "Setting up your computer to process Nanopore sequencing data"
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

# Install Java (required for some tools)
if ! command_exists java; then
    echo "Java is not installed. Installing Java..."
    brew install openjdk
    echo 'export PATH="/usr/local/opt/openjdk/bin:$PATH"' >> ~/.bash_profile
    source ~/.bash_profile
else
    echo "Java is already installed."
fi

# Array of prerequisites with alternative installation methods
PREREQUISITES=("fastqc" "guppy" "minimap2" "samtools" "featureCounts")

# Loop through each prerequisite and install if not present
for TOOL in "${PREREQUISITES[@]}"; do
    case $TOOL in
        "guppy")
            if ! command_exists guppy_basecaller; then
                echo "Guppy is not installed. Please download and install Guppy from the official Oxford Nanopore website."
            else
                echo "Guppy is already installed."
            fi
            ;;
        "minimap2")
            if ! command_exists minimap2; then
                echo "Minimap2 is not installed. Installing Minimap2..."
                brew install brewsci/bio/minimap2
            else
                echo "Minimap2 is already installed."
            fi
            ;;
        "featureCounts")
            if ! command_exists featureCounts; then
                echo "FeatureCounts is not installed. Installing FeatureCounts..."
                brew install brewsci/bio/subread
            else
                echo "FeatureCounts is already installed."
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
