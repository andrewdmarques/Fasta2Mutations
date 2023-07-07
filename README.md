# fasta2mutations

fasta2mutations is a comprehensive program that provides an analytical pipeline to perform various operations on DNA sequences of SARS-CoV-2 virus, including aligning sequences, detecting mutations, mapping those mutations to the relevant gene products, and generating a mutation table.

## Installation

1. Clone the repository to your local machine. This will create a local copy of the codebase on your system:

   ```shell
   git clone https://github.com/andrewdmarques/Fasta2Mutations.git
   ```

2. Install the required R dependencies. The program relies on the `ape` library, a package for analyses of phylogenetics and evolution in R, which can be installed in R using the following command:

   ```R
   install.packages("ape")
   ```

## Requirements

The following are the prerequisites to run fasta2mutations:

- R version 3.6 or higher installed on your machine.

- The `ape` library installed in your R environment.

## Input

The program requires two key input files in the appropriate format:

1. Reference sequence FASTA file (`nc_045512.2.fasta`): This file should contain the reference DNA sequence. It is used to compare with other sequences to detect variations.

2. Sequence file with variants (`sequences.fasta`): This file should contain the DNA sequences of various SARS-CoV-2 variants. The program will compare these sequences with the reference sequence to identify mutations.

Place these input files in the root directory of the project.

## Output

fasta2mutations generates the following outputs:

- `mutation_table.csv`: This CSV file contains a detailed mutation table that includes information about each identified mutation. The table includes columns for sequence numbers, fasta sequences, file names, sequence names, and mutation details. This table can be further analyzed for insights about the distribution and frequency of mutations.

## Usage

To use fasta2mutations, follow to these steps:

1. Ensure that the input files are in the correct format (FASTA for sequences) and placed in the root directory of the project.

2. Open the R script in an R environment (like RStudio or a Jupyter notebook with an R kernel) and run the script.

The script will then perform the necessary operations, detecting and characterizing mutations in the input sequences relative to the reference sequence. After execution, it will generate a comprehensive mutation table as an output.


