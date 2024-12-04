# RNA Sequencing analysis pipelines
Welcome to our code repo where you can find counts table generation pipelines for various sequencing platforms for the 3DBM & NE labs at Neurosurgery Freiburg

## Optimization
These scripts have been validated on MacOS based systems

## General Instructions
Each pipeline should be preceeded with its associated **'./check_prerequisites.sh'** file, so that all dependencies are installed as needed for the pipeline to run.

## Pipeline instructions
1. Download the folder of the containing the scripts pertaining to the platform that you have data from.
3. Open terminal and make the files excutable using the command **'chmod +x **.sh'**, with **"**.sh"** being the files that were downloaded.
4. In terminal, run **'./check_prerequisites.sh'** to make sure that all dependencies are installed and available for use.
5. In terminal, run **'./platform_count_matrix.sh path_to_fastq_files path_to_reference_file'**. Ensure you leave a space between each argument

