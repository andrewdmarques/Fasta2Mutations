# Load libraries
library(ape)
# library(ggplot2)

# Variables to assign.
file_ref <-   "nc_045512.2.fasta"
file_seqs <-  "2021-12-01_sequence-file-names.txt"


# Make directories if they do not exist.
if(!file.exists("00_Results")){system("mkdir 00_Results")}
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
      if(align_mismatch$seq2[j] == "-" | 
         align_mismatch$seq2[j] == "a" | 
         align_mismatch$seq2[j] == "c" | 
         align_mismatch$seq2[j] == "t" |
         align_mismatch$seq2[j] == "g"){
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

# Save the reference file as a csv. 
write.csv(ref, "./00_Results/reference.csv")

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

# Reorder the mutations data frame.
mut <- mut[,c(4,1,2,3)]

write.csv(mut, "./00_Results/nt-mutations.csv")

# Compare the mutations from the samples.
mut1 <- mut

# make a column for the gene/mutation combination.
mut1$mutation <- paste0(mut1$REF, mut1$POS, mut1$ALT)
# for(i in 1:length(mut1$VSP)){
#   if(is.na(mut1$type[i]) | mut1$type[i] == "silent" | grepl("del", mut1$type[i])){
#     mut1$mutation[i] <- paste0(mut1$genes[i], "_", mut1$type[i], "_", as.character(mut1$POS[i]))
#   }
# }

# Determine all of the unique samples.
sam1 <- unique(mut1$VSP)
mut1$sample <- mut1$VSP
sam2 <- unique(mut1$sample)

# Determine all of the unique mutations.
mutation_unique <- unique(mut1$mutation)

# Determine the order to present the mutations.
mutation_positions <- data.frame(matrix(NA, nrow = length(mutation_unique), ncol = 2))
colnames(mutation_positions) <- c("mutation", "POS")
mutation_positions$mutation <- mutation_unique
# Add the genomic positions to each mutation.
for(i in 1:length(mutation_positions$mutation)){
  temp <- subset(mut1, mut1$mutation == mutation_unique[i])
  mutation_positions$POS[i] <- temp$POS[1]
}
mutation_positions <- mutation_positions[order(mutation_positions$POS),]

mutation_sort <- mutation_positions$mutation

# Make a blank data frame that is possible mutations x sample.
mut2 <- data.frame(matrix(NA, nrow = length(sam1), ncol = length(mutation_sort)))
colnames(mut2) <- mutation_sort
rownames(mut2) <- sam1

# Populate the data frame.
for(i in 1:length(mutation_sort)){
  for(j in 1:length(sam1)){
    temp <- subset(mut1, mut1$mutation == mutation_sort[i])
    temp <- subset(temp, temp$VSP == sam1[j])
    if(length(temp$VSP) > 0){
      mut2[j,i] <- "Yes"
    }else{
      mut2[j,i] <- ""
    }
  }
}

# Rename the rows to have both the sample and VSP identifier.
rownames(mut2) <- sam2

mut2 <- data.frame(t(mut2))

# Save the file as a csv. 
write.csv(mut2, "./00_Results/nt-mutations-table.csv")

# For cat samples - Determine if there are more mutations that match USDA or our sequences.
mut3 <- mut2
mut3$match <- "Other Two Only"
for(i in 1:length(mut3$match)){
  if(mut3$VSP3508[i] == 'Yes'){
  if(mut3$VSP3508[i] == mut3$usda.cat.1[i]){
    mut3$match[i] <- "USDA"
  }
  if(mut3$VSP3508[i] == mut3$VSP3454[i]){
    mut3$match[i] <- "Bushman"
  }
  if(mut3$VSP3508[i] == mut3$usda.cat.1[i] && mut3$VSP3508[i] == mut3$VSP3454[i]){
    mut3$match[i] <- "Both"
  }
    if(mut3$match[i] == "Other Two Only"){
      mut3$match[i] <- "only this one"
    }
  }
}

# Determine the number of matches for each group.
table(mut3$match)

# For cat samples - Determine if there are more mutations that match USDA or our sequences.
mut4 <- mut2
mut4$match <- "Other Two Only"
for(i in 1:length(mut3$match)){
  if(mut4$VSP3454[i] == 'Yes'){
  if(mut4$VSP3454[i] == mut4$usda.cat.1[i]){
    mut4$match[i] <- "USDA"
  }
  if(mut4$VSP3454[i] == mut4$VSP3508[i]){
    mut4$match[i] <- "Bushman"
  }
  if(mut4$VSP3454[i] == mut4$usda.cat.1[i] && mut4$VSP3508[i] == mut4$VSP3454[i]){
    mut4$match[i] <- "Both"
  }
    if(mut4$match[i] == "Other Two Only"){
      mut4$match[i] <- "only this one"
    }
  }
}

# Determine the number of matches for each group.
library(ggvenn)

mut5 <- data.frame(matrix(NA, nrow = length(mut1$VSP), ncol = 2))
colnames(mut5) <- c("sample", "mutation")
mut5$sample <- mut1$sample
mut5$mutation <- mut1$mutation

a <- subset(mut5, mut5$sample == "usda-cat-1")
b <- subset(mut5, mut5$sample == "VSP3454")
c <- subset(mut5, mut5$sample == "VSP3508")

x <- list(
  USDA = a$mutation, 
  Bushman_Recent = c$mutation,
  Bushman_Previous = b$mutation
)

ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE
)

