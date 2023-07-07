# Load libraries
library(ape)
library(genbankr)

# Variables to assign.
file_ref <-   "nc_045512.2.fasta"
# file_seqs <-  "2021-12-01_sequence-file-names.txt"
file_original <- "2022-01-10_deer_all.fasta"
# Determine what the coding regions are for the reference.
cds <- data.frame(genbankr::cds(genbankr::readGenBank("Wuhan-Hu-1.gb")))  # cp /home/common/SARS-CoV-2-Philadelphia/data/references/Wuhan-Hu-1.gb ./

# Open the data
fas1 <- read.table(file_original, header = FALSE, sep = "", dec = ".")
names(fas1)[1] <- "fasta"
fas1$fasta <- as.character(fas1$fasta)

# Make a reference file.
ref <- data.frame(matrix(NA, nrow = length(fas1$fasta), ncol = 5))
colnames(ref) <- c("seq_num", "fasta", "file_name", "seq_name", "mutations")

# Condense the fasta sequence to be on one line.
fas2 <- data.frame(matrix(ncol = 1, nrow = length(fas1$fasta)))
names(fas2)[1] <- "fasta"
x <- 0
seq_num <- 1
for_file <- ""
for(i in 1:length(fas1$fasta)){
  if(grepl(">",fas1$fasta[i]) == TRUE){
    x <- x + 1
    fas2$fasta[x] <- paste0(">seq_", as.character(seq_num), "_cat")
    # Record the reference information.
    ref$seq_num[seq_num] <- seq_num
    ref$seq_name[seq_num] <- paste0(">seq_", as.character(seq_num))
    ref$fasta[seq_num] <- fas1$fasta[i]
    ref$file_name[seq_num] <- paste0("./01_Sequence/seq_", as.character(seq_num), ".fasta")
    ref$wgs[seq_num] <- "unknown"
    seq_num <- seq_num + 1
    ref$mutations <- 0
    x <- x + 1
    # Output how many sequences have been adapted.
    if (x%%10 == 0) {
      out_text <- paste(as.character(x/2), " sequences have been curated", "")
      print(out_text, quote = FALSE)
    }
  }else{
    fas2$fasta[x] <- gsub(" ", "", paste(fas2$fasta[x], fas1$fasta[i]))
  }
}

# Remove rows with na.
fas3 <- fas2
fas3 <- na.omit(fas3)
ref <- na.omit(ref)

# Make directories if they do not exist.
if(!file.exists("00_Results")){system("mkdir 00_Results")}
if(!file.exists("01_Sequence")){system("mkdir 01_Sequence")}
if(!file.exists("02_Align")){system("mkdir 02_Align")}
if(!file.exists("03_Mutation")){system("mkdir 03_Mutation")}

# Determine the files that should be compared.
# files <- read.table(file_seqs)
seq_files <- as.character(ref$file_name)

# Make a file for each fasta to be aligned.
temp <- data.frame(matrix(NA, nrow = 2, ncol = 1))
colnames(temp) <- c("fasta")
for(i in 1:length(fas3$fasta)){
  # i <- 1
  if(grepl(">",fas3$fasta[i]) == TRUE){
    if(nchar(fas3$fasta[i + 1]) > 29000){
      ref$wgs[((i+1)/2)] <- "y"
      seq <- sub('.', '', fas3$fasta[i + 1])
      seq <- sub('.', '', seq)
      fasta_name <- ref$fasta[((i+1)/2)]
      file_name <- ref$file_name[((i+1)/2)]
      temp$fasta[1] <- fasta_name
      temp$fasta[2] <- seq
      write.table(temp, file_name, quote = F, col.names=F, row.names = F)
    }
  }
}

# For each sequence, perform an alignment, and determine the nt/aa mutations.
for(i in 17:length(seq_files)){
  # i <- 1
  if(ref$wgs[i] == "y"){
  # Record the information to the ref file.
  temp_file <- as.character(seq_files[i])
  temp_seq <- read.table(seq_files[i])
  temp_seq <- temp_seq$V1[1]
  temp_seq <- as.character(temp_seq)
  ref$file_name[i] <- seq_files[i]
  ref$seq_name[i] <- gsub(">", "", temp_seq)
 
   
  # Perform a MSA using muscle.
  out_file <- paste0("./02_Align/", gsub(".fasta", "_align.fasta",gsub("./01_Sequence/","",temp_file)))
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
  
  # Inserted nucleotides are inserted for the nt position immediately before the position.
  align_df$POS <- as.numeric(rownames(align_df))
  temp_p <- as.numeric(0)
  for(j in 1:length(align_df$POS)){
    if(align_df$seq1[j] == "-"){
      align_df$POS[j] <- temp_p
    }else{
      temp_p <- temp_p + 1
      align_df$POS[j] <- temp_p
    }
  }
  
  align_df_original <- align_df
  
  # Condensed indels so that they are considered one mutation.
  align_df$filter <- "n"
  align_df$filter1 <- "n"
  ins_new <- "y"
  ins_pos <- 0
  del_new <- "y"
  del_pos <- 0
  del_num <- 0
  for(j in 1:length(align_df$POS)){
    
    # If insertion, condense it to the first occurance of insertion. 
    if(align_df$seq1[j] == "-"){
      if(ins_new == "y"){
        align_df$seq2[j] <- paste0("ins_", align_df$seq2[j])
        align_df$filter[j] <- "y"
        ins_pos <- j
        ins_new <- "n"
      }else{
        align_df$seq2[ins_pos] <- paste0(align_df$seq2[ins_pos],align_df$seq2[j])
      }
    }else{
      align_df$filter[j] <- "y"
      ins_new <- "y"
    }
    
    # If deletion, condense it to the first occurrance of deletion.
    if(align_df$seq2[j] == "-"){
      if(del_new == "y"){
        align_df$seq2[j] <- paste0("del_", as.character(del_num))
        align_df$filter1[j] <- "y"
        del_pos <- j
        del_new <- "n"
        del_num <- del_num + 1
      }
      if(del_new == "n"){
        align_df$seq2[del_pos] <- paste0("del_", as.character(del_num))
        del_num <- del_num + 1
      }
    }else{
      align_df$filter1[j] <- "y"
      del_new <- "y"
      del_num <- 0
    }
  }
  
  # Remove the lines that were condensed from being an insertion or deletion.
  align_df <- subset(align_df, align_df$filter == "y")
  align_df <- subset(align_df, align_df$filter1 == "y")
  align_df <- subset(align_df, select = -c(filter, filter1) )
  
  # Select only the nucleotides that do not match.
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
      if(grepl("del",align_mismatch$seq2[j])){
        align_mismatch$filter[j] <- "y"
      }
      if(grepl("ins",align_mismatch$seq2[j])){
        align_mismatch$filter[j] <- "y"
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
  
  # If position 1 is a mismatch, this is because of the deletions.
  if(grepl("del", align_mismatch$seq2[1])){
    align_mismatch$filter[1] <- "n"
  }
  
  # Remove the files that are n's and those from the beginning and end of the sequencing batches.
  align_mismatch2 <- subset(align_mismatch, align_mismatch$filter == "y")
  align_mismatch2 <- subset(align_mismatch2, align_mismatch2$seq1 != "n")
  align_mismatch2 <- subset(align_mismatch2, align_mismatch2$seq2 != "n")
  ref$mutations[i] <- length(align_mismatch2$filter)
  
  # Save the SNP file.
  align_mismatch3 <- align_mismatch2[ , ! names(align_mismatch2) %in% c("filter")]
  colnames(align_mismatch3) <- c(col_name, "POS")
  # write.csv(align_mismatch3, paste0("./03_Mutation/", gsub(".fasta", "_snp.csv",gsub("./01_Sequence/","",temp_file))))
  # if (i%%10 == 0) {
  #   out_text <- paste(as.character(i), " fasta sequence processed", "")
  #   print(out_text, quote = FALSE)
  # }
  
  # Determine the AA Mutation ________________________________________________________
  align_mismatch4 <- align_mismatch3
  align_mismatch4$genes <- ""
  align_mismatch4$protein_product <- ""
  align_mismatch4$type <- ""
  
  # Record the protein product.
  for(j in 1:length(align_mismatch4$POS)){
    if(align_mismatch4$POS[j] >= 266 && align_mismatch4$POS[j] < 805){align_mismatch4$protein_product[j] <- "Nsp1"}
    if(align_mismatch4$POS[j] >= 806 && align_mismatch4$POS[j] < 2719){align_mismatch4$protein_product[j] <- "Nsp2"}
    if(align_mismatch4$POS[j] >= 2720 && align_mismatch4$POS[j] < 8554){align_mismatch4$protein_product[j] <- "Nsp3"}
    if(align_mismatch4$POS[j] >= 8555 && align_mismatch4$POS[j] < 10054){align_mismatch4$protein_product[j] <- "Nsp4"}
    if(align_mismatch4$POS[j] >= 10055 && align_mismatch4$POS[j] < 10972){align_mismatch4$protein_product[j] <- "3CL-PRO"}
    if(align_mismatch4$POS[j] >= 10973 && align_mismatch4$POS[j] < 11842){align_mismatch4$protein_product[j] <- "Nsp6"}
    if(align_mismatch4$POS[j] >= 11843 && align_mismatch4$POS[j] < 12091){align_mismatch4$protein_product[j] <- "Nsp7"}
    if(align_mismatch4$POS[j] >= 12092 && align_mismatch4$POS[j] < 12685){align_mismatch4$protein_product[j] <- "Nsp8"}
    if(align_mismatch4$POS[j] >= 12686 && align_mismatch4$POS[j] < 13024){align_mismatch4$protein_product[j] <- "Nsp9"}
    if(align_mismatch4$POS[j] >= 13025 && align_mismatch4$POS[j] < 13441){align_mismatch4$protein_product[j] <- "Nsp10"}
    if(align_mismatch4$POS[j] >= 13442 && align_mismatch4$POS[j] < 16236){align_mismatch4$protein_product[j] <- "Pol"}
    if(align_mismatch4$POS[j] >= 16237 && align_mismatch4$POS[j] < 18039){align_mismatch4$protein_product[j] <- "Hel"}
    if(align_mismatch4$POS[j] >= 18040 && align_mismatch4$POS[j] < 19620){align_mismatch4$protein_product[j] <- "ExoN"}
    if(align_mismatch4$POS[j] >= 19621 && align_mismatch4$POS[j] < 20658){align_mismatch4$protein_product[j] <- "Nsp15"}
    if(align_mismatch4$POS[j] >= 20659 && align_mismatch4$POS[j] < 21552){align_mismatch4$protein_product[j] <- "Nsp16"}
    if(align_mismatch4$POS[j] >= 21563 && align_mismatch4$POS[j] < 25381){align_mismatch4$protein_product[j] <- "Spike"}
    # if(align_mismatch4$POS[j] >= 21599 && align_mismatch4$POS[j] < 23617){align_mismatch4$protein_product[j] <- "S1"}
    # if(align_mismatch4$POS[j] >= 23618 && align_mismatch4$POS[j] < 25381){align_mismatch4$protein_product[j] <- "S2"}
    if(align_mismatch4$POS[j] >= 25393 && align_mismatch4$POS[j] < 26217){align_mismatch4$protein_product[j] <- "ORF3a"}
    if(align_mismatch4$POS[j] >= 26245 && align_mismatch4$POS[j] < 26469){align_mismatch4$protein_product[j] <- "Envelope"}
    if(align_mismatch4$POS[j] >= 26523 && align_mismatch4$POS[j] < 27188){align_mismatch4$protein_product[j] <- "Membrane"}
    if(align_mismatch4$POS[j] >= 27202 && align_mismatch4$POS[j] < 27384){align_mismatch4$protein_product[j] <- "ORF6"}
    if(align_mismatch4$POS[j] >= 27439 && align_mismatch4$POS[j] < 27756){align_mismatch4$protein_product[j] <- "ORF7a"}
    if(align_mismatch4$POS[j] >= 27756 && align_mismatch4$POS[j] < 27884){align_mismatch4$protein_product[j] <- "ORF7b"}
    if(align_mismatch4$POS[j] >= 27939 && align_mismatch4$POS[j] < 28256){align_mismatch4$protein_product[j] <- "ORF8"}
    if(align_mismatch4$POS[j] >= 28274 && align_mismatch4$POS[j] < 29530){align_mismatch4$protein_product[j] <- "Nucleocapsid"}
    if(align_mismatch4$POS[j] >= 29558 && align_mismatch4$POS[j] < 29671){align_mismatch4$protein_product[j] <- "ORF10"}
  }
 
  # Map the mutations from the sample to the reference, excluding the deletion/insertions.
  seq_sample <- align_df_original
  seq_sample <- seq_sample[ , ! names(seq_sample) %in% c("seq2")]
  # For each mutation that is not a deletion or insertion, map it to the reference sequence.
  for(j in 1:length(align_mismatch2$seq1)){
    if(!grepl("del", align_mismatch2$seq2[j])){
      if(!grepl("ins", align_mismatch2$seq2[j])){
        seq_sample[as.numeric(align_mismatch2$POS[j]),1] <- align_mismatch2$seq2[j]
      }
    }
  }

  # Iterate through each of the coding sequences.
  for(j in 1:length(cds$seqnames)){
    # Based on the known coding regions, slice the mapped sample fasta to be the nt fragment for the coding region.
    p_start <- cds$start[j]
    p_end <- cds$end[j]
    for(k in 1:length(align_mismatch4$POS)){
      if(align_mismatch4$POS[k] >= p_start){
        if(align_mismatch4$POS[k] <= p_end){
          align_mismatch4$genes[k] <- cds$gene[j]
        }
      }
    }
    temp_p <- subset(seq_sample, seq_sample$POS >= p_start)
    temp_p <- subset(temp_p, temp_p$POS <= p_end)
    colnames(temp_p) <- ">seq2"
    temp_p <- temp_p[ c(1) ]
    write.table(temp_p, "temp_nt.fasta", append = FALSE, sep = "", dec = ".",
                row.names = FALSE, col.names = TRUE, quote = FALSE)

    # Translate the sample fasta file.
    temp_p <- ape::read.dna("temp_nt.fasta", format = "fasta")
    temp_p <- ape::trans(temp_p)
    ape::write.FASTA(temp_p, "temp_aa_sample.fasta")

    # Determine what that translated products are for the reference.
    ref_protein <- cds$translation[j]
    
    # Combine the reference and sasmple translation into one data frame.
    ref_protein <- data.frame(ref_protein)
    colnames(ref_protein) <- "seq1"
    temp_p <- read.table("temp_aa_sample.fasta", header = TRUE)
    colnames(temp_p) <- "seq2"
    temp_p <- cbind(temp_p,ref_protein)
    temp_p$seq1 <- as.character(temp_p$seq1)
    temp_p$seq2 <- as.character(temp_p$seq2)
    
    # Iterate through the nucleotide mutations within this coding region.
    for(k in 1:length(align_mismatch4$POS)){
      if(align_mismatch4$POS[k] >= p_start){
        if(align_mismatch4$POS[k] <= p_end){
          # Determine what the aa position is for that mutation.
          aa_pos <- floor((align_mismatch4$POS[k] - p_start)/3) + 1
          # Determien what the ref aa is.
          aa_ref <- substr(temp_p$seq1[1],aa_pos,aa_pos)
          # Determine what the sample aa is.
          aa_sam <- substr(temp_p$seq2[1],aa_pos,aa_pos)
          # Record the mutation.
          if(aa_sam == aa_ref){
            temp_type <- "silent"
          }else{
            temp_type <- paste0(aa_ref, as.character(aa_pos), aa_sam)
          }
          align_mismatch4$type[k] <- temp_type
        }
      }
    }
    
    # Curate the mutation file to show the mutation type for deletion and insertion as correct. 
    temp <- align_mismatch4[,c(2)]
    for(k in 1:length(align_mismatch4$POS)){
      if(align_mismatch4$genes[k] == ""){
        align_mismatch4$genes[k] <- "intergenic"
      }
      # Label the type as deletion or insertion. 
      if(grepl("del", temp[k])){
        align_mismatch4$type[k] <- temp[k]
      }
      if(grepl("ins", temp[k])){
        align_mismatch4$type[k] <- temp[k]
      }
    }
  }

  }
  # Save the SNP file.
  write.csv(align_mismatch4, paste0("./03_Mutation/", gsub(".fasta", "_snp.csv",gsub("./01_Sequence/","",temp_file))))
  if (i%%10 == 0) {
    out_text <- paste(as.character(i), " fasta sequence processed", "")
    print(out_text, quote = FALSE)
  }
  # Remove all of the temp files.
  command <- "rm -r ./temp*"
  # Write a bash script that will run the blast.
  filepath01 <- '.'
  p <- c('#!/bin/bash', "export PATH=${PATH}:${HOME}/edirect", command)
  writeLines(p, file.path(filepath01, 'temp.script'))
  system(paste0('chmod 755 ', file.path(filepath01, 'temp.script')))
  comm <- paste0(file.path(filepath01, 'temp.script'))
  # Run the bash script.
  system(comm)
}

# Save the reference file as a csv. 
write.csv(ref, "./00_Results/reference.csv")

# Compile the mutations into one file.
mut <- data.frame(matrix(NA, nrow = 0, ncol = 7))
colnames(mut) <- c("VSP", "POS", "REF", "ALT", "genes", "protein_product", "type")

# Condense all of the mutations to one data frame.
for(i in 1:length(ref$file_name)){
  # i <- 3
  if(ref$wgs[i] == "y"){
  temp_mutation <- read.csv(paste0("./03_Mutation/", gsub(".fasta", "_snp.csv",gsub("./01_Sequence/","",ref$file_name[i]))))
  # temp_mutation$VSP <- gsub("/.*","",gsub(".*-","",ref$seq_name[i]))
  temp_mutation$VSP <- ref$seq_name[i]
  temp_mutation = subset(temp_mutation, select = -c(X) )
  colnames(temp_mutation) <- c("REF", "ALT", "POS", "genes", "protein_product", "type", "VSP")
  temp_mutation <- temp_mutation[,c(7,3,1,2,4,5,6)]
  mut <- rbind(mut,temp_mutation)
  }
}

write.csv(mut, "./00_Results/mutations.csv", row.names = FALSE)