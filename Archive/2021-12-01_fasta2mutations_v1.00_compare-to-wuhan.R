# Load libraries
library(ape)
# library(ggplot2)

# Variables to assign.
file_ref <-   "nc_045512.2.fasta"
file_seqs <-  "2021-12-01_sequence-file-names.txt"


# Make directories if they do not exist.
if(!file.exists("01_Align")){system("mkdir 01_Align")}
if(!file.exists("02_Mutation")){system("mkdir 02_Mutation")}

# Determine the files that should be compared.
files <- read.table(file_seqs)
files$V1 <- as.character(files$V1)

# Make a reference file.
ref <- data.frame(matrix(NA, nrow = length(files$V1), ncol = 3))
colnames(ref) <- c("file_name", "seq_name", "mutations")

# For each sequence, perform an alignment.
for(i in 1:length(files$V1)){
  # i <- 1
  
  # Record the information to the ref file.
  temp_file <- as.character(files$V1[i])
  temp_seq <- read.table(files$V1[i])
  temp_seq <- temp_seq$V1[1]
  temp_seq <- as.character(temp_seq)
  ref$file_name[i] <- files$V1[i]
  ref$seq_name[i] <- gsub(">", "", temp_seq)
 
   
  # Perform a MSA using muscle.
  out_file <- paste0("./01_Align/", gsub(".fasta", "_align.fasta",temp_file))
  command <- paste0("./muscle3.8.31_i86linux64 -profile -in1 ", file_ref, " -in2 ", temp_file, " -out ", out_file)
  
  # Write a bash script that will run the blast.
  filepath01 <- '.'
  p <- c('#!/bin/bash', "export PATH=${PATH}:${HOME}/edirect", command)
  writeLines(p, file.path(filepath01, 'temp.script'))
  system(paste0('chmod 755 ', file.path(filepath01, 'temp.script')))
  comm <- paste0(file.path(filepath01, 'temp.script'))
  
  # Run the bash script.
  system(comm)
  # Open the sequence as a variable.
  seq <- ape::read.dna(out_file, format = "fasta")
  
  # Determine the SNP distances from the MSA.
  align_char <- as.character(seq)
  align_df <- data.frame(align_char)
  l <- length(align_df)
  
  # Transpose the data frame so that it has 1 row for each position.
  align_df <- as.data.frame(t(as.matrix(align_df)))
  
  # Change the row index so that it allows it is the numbers and column names so that they are uniform.
  rownames(align_df) <- seq(1:l)
  col_name <- colnames(align_df)
  colnames(align_df) <- c("seq1", "seq2")
  align_df$seq1 <- as.character(align_df$seq1)
  align_df$seq2 <- as.character(align_df$seq2)
  align_mismatch <- subset(align_df, align_df$seq1 != align_df$seq2)
  
  # Remove the first and last group of lines that might results from sequencing (threshould is 50 continuous mathcing bp).
  cutoff <- 50
  if(length(align_mismatch$seq1) >= 1){
    align_mismatch$filter <- "n"    # Automatically set the filter to remove sequences.
    for(j in 1:length(align_mismatch$filter)){
      if(align_mismatch$seq2[j] == "a" | align_mismatch$seq2[j] == "c" | align_mismatch$seq2[j] == "t" |align_mismatch$seq2[j] == "g"){
        align_mismatch$filter[j] <- "y" # If the nucleotide is not ambiguous, then it will not be set to not be filtered.
      }
    }
  }else{
    align_mismatch <- data.frame(matrix(NA, nrow = 0, ncol = 3))
    colnames(align_mismatch) <- c("seq1", "seq2", "filter")
  }
  if(length(align_mismatch$filter) > 1){
    # Remove the head end of mismatches.
    for(j in 2:length(align_mismatch$filter)){
      line <- j # Perform the function starting with the end first.
      if(as.numeric(rownames(align_mismatch)[line]) - as.numeric(rownames(align_mismatch)[line-1]) >= cutoff){
        break
      }
      align_mismatch$filter[line] <- "n"
      align_mismatch$filter[line - 1] <- "n"
    }
    
    # Remove the tail end of mismatches.
    align_mismatch$filter[length(align_mismatch$filter)] <- "n"
    for(j in 1:length(align_mismatch$filter)){
      line <- length(align_mismatch$filter) - j # Perform the function starting with the end first.
      if(as.numeric(rownames(align_mismatch)[line + 1]) - as.numeric(rownames(align_mismatch)[line]) >= cutoff){
        break
      }
      align_mismatch$filter[line] <- "n"
      align_mismatch$filter[line + 1] <- "n"
    }
  }
  
  # Remove the files that are n's and those from the beginning and end of the sequencing batches.
  align_mismatch2 <- subset(align_mismatch, align_mismatch$filter == "y")
  align_mismatch2 <- subset(align_mismatch2, align_mismatch2$seq1 != "n")
  align_mismatch2 <- subset(align_mismatch2, align_mismatch2$seq2 != "n")
  ref$mutations[i] <- length(align_mismatch2$filter)
  
  # Save the SNP file.
  align_mismatch3 <- align_mismatch2[ , ! names(align_mismatch2) %in% c("filter")]
  colnames(align_mismatch3) <- col_name
  write.csv(align_mismatch3, paste0("./02_Mutation/", gsub(".fasta", "_snp.csv",ref$file_name[i])))
  if (i%%10 == 0) {
    out_text <- paste(as.character(i), " closest fasta sequence processed", "")
    print(out_text, quote = FALSE)
  }
}


# Compile the mutations into one file.
mut <- data.frame(matrix(NA, nrow = 0, ncol = 5))
colnames(mut) <- c("VSP", "POS", "REF", "ALT")

# Condense all of the mutations to one data frame.
for(i in 1:length(ref$file_name)){
  # i <- 3
  temp_mutation <- read.csv(paste0("./02_Mutation/", gsub(".fasta", "_snp.csv",ref$file_name[i])))
  colnames(temp_mutation) <- c("POS", "REF", "ALT")
  temp_mutation$VSP <- ref$seq_name[i]
  mut <- rbind(mut,temp_mutation)
}

write.csv(mut, "nt_mutations.csv")

