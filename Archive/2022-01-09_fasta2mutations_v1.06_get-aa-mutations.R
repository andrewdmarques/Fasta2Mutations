# Load libraries
library(ape)
library(genbankr)
# library(ggplot2)

# Variables to assign.
file_ref <-   "nc_045512.2.fasta"
# file_seqs <-  "2021-12-01_sequence-file-names.txt"
file_original <- "2022-01-09_vsp3514-and-vsp3515.fasta"
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
for(i in 1:length(seq_files)){
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
  write.csv(align_mismatch3, paste0("./03_Mutation/", gsub(".fasta", "_snp.csv",gsub("./01_Sequence/","",temp_file))))
  if (i%%10 == 0) {
    out_text <- paste(as.character(i), " fasta sequence processed", "")
    print(out_text, quote = FALSE)
  }
  
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
  # # Save the mapped sample sequence as a fasta file.
  # colnames(seq_sample) <- ">seq2"
  # seq_sample <- seq_sample[ c(1) ]
  # write.table(seq_sample, "temp.fasta", append = FALSE, sep = "", dec = ".",
  #             row.names = FALSE, col.names = TRUE, quote = FALSE)
  # # Translate the sample fasta file.
  # temp <- ape::read.dna("temp.fasta", format = "fasta")
  # temp <- ape::trans(temp)
  # ape::write.dna(temp, "temp_aa_sample.fasta", format = "fasta")
  # 
  # # Translate the reference fasta file.
  # temp <- ape::read.dna(file_ref, format = "fasta")
  # temp <- ape::trans(temp)
  # ape::write.dna(temp, "temp_aa_ref.fasta", format = "fasta")
  # 
  # 
  # cds <- gb@cds
  # seqlevels(cds) <- 'genome'
  # seqnames(cds)  <- 'genome'
  # 
  # # Calculate how the shift left or right caused by deletions and insertions.
  # opt$variantTableMajor$shift <- ifelse(grepl('del', opt$variantTableMajor$ALT), (nchar(opt$variantTableMajor$ALT)-3)*-1, 0)
  # opt$variantTableMajor$shift <- ifelse(grepl('ins', opt$variantTableMajor$ALT), (nchar(opt$variantTableMajor$ALT)-3), opt$variantTableMajor$shift)
  # 
  # 
  # # Remove variant positions flanking indels since they appear to be artifacts.
  # artifacts <- c(opt$variantTableMajor[grep('ins|del', opt$variantTableMajor$ALT),]$POS + abs(opt$variantTableMajor[grep('ins|del', opt$variantTableMajor$ALT),]$shift),
  #                opt$variantTableMajor[grep('ins|del', opt$variantTableMajor$ALT),]$POS -1)
  # 
  # opt$variantTableMajor <- opt$variantTableMajor[! opt$variantTableMajor$POS %in% artifacts,]
  # 
  # opt$variantTableMajor <- bind_rows(lapply(split(opt$variantTableMajor, 1:nrow(opt$variantTableMajor)), function(x){
  # 
  #   # Determine the offset of this position in the concensus sequence because it may not be the same length
  #   # if indels have been applied. Here we sum the indel shifts before this variant call.
  # 
  #   offset <- sum(opt$variantTableMajor[1:grep(x$POS, opt$variantTableMajor$POS),]$shift)
  # 
  #   cds2 <- cds
  #   start(cds2) <- start(cds2) + offset
  #   end(cds2) <- end(cds2) + offset
  # 
  #   v1 <- GRanges(seqnames = 'genome', ranges = IRanges(x$POS, end = x$POS), strand = '+')
  #   o1 <- GenomicRanges::findOverlaps(v1, cds)
  # 
  #   v2 <- GRanges(seqnames = 'genome', ranges = IRanges(x$POS + offset, end = x$POS + offset), strand = '+')
  #   o2 <- GenomicRanges::findOverlaps(v2, cds2)
  # 
  #   if(length(o2) == 0){
  #     x$genes <- 'intergenic'
  # 
  #     if (grepl('ins', as.character(x$ALT))){
  #       x$type <- paste0('ins ', nchar(x$ALT)-3)
  #     } else if (grepl('del', as.character(x$ALT))){
  #       x$type <- paste0('del ', nchar(x$ALT)-3)
  #     } else {
  #       x$type <- ' '
  #     }
  #   } else {
  # 
  #     # Define the gene the variant is within.
  #     hit1 <- cds[subjectHits(o1)]
  #     hit2 <- cds2[subjectHits(o2)]
  # 
  #     x$genes <- paste0(hit2$gene, collapse = ', ')
  # 
  #     # Native gene AA sequence.
  #     orf1  <- as.character(translate(DNAString(substr(as.character(readFasta(opt$refGenomeFasta)@sread), start(hit1), end(hit1))), if.fuzzy.codon = 'solve'))
  # 
  #     # Variant gene AA sequence.
  #     orf2 <- as.character(translate(DNAString(substr(opt$concensusSeq, start(hit2), end(hit2))), if.fuzzy.codon = 'solve'))
  # 
  #     # Determine the offset of this position in the concensus sequence because it may not be the same length
  #     # if indels have been applied. Here we sum the indel shifts before this variant call.
  #     # offset <- sum(opt$variantTableMajor[1:grep(x$POS, opt$variantTableMajor$POS),]$shift)
  # 
  #     #              1   2   3   4   5   6   7   8
  #     # 123 456 789 012 345 678 901 234 567 890 123
  #     # ATG CAT TGA ATG GGC TTA CGA GCT TAA GTA TAG
  #     #             ^             x  21-10 + 2 = 13/3 = 4.3 ~ 4
  #     #                          x   20-10 + 2 = 12/3 = 4.0 = 4
  #     #                         x    19-10 + 2 = 11/3 = 3.6 ~ 4
  #     #                                 x   25-10 + 2 = 17/3 = 5.6 ~ 6
  #     #                                  x  26-10 + 2 = 18/3 = 6.0 = 6
  #     #                                   x 27-10 + 2 = 19/3 = 6.3 ~ 6
  # 
  #     aa <- round(((x$POS - start(hit1)) + 2)/3)
  #     orf_aa <- substr(orf1, aa, aa)
  # 
  #     aa2 <- round((((x$POS + offset) - start(hit2)) + 2)/3)
  #     orf2_aa <- substr(orf2, aa2, aa2)
  # 
  #     maxALTchars <- max(nchar(unlist(strsplit(as.character(x$ALT), ','))))
  # 
  #     if(nchar(as.character(x$REF)) == 1 & nchar(as.character(x$ALT)) > 1 & maxALTchars == 1){
  #       x$type <- paste0(x$POS, '_mixedPop')
  #     } else if (grepl('ins', as.character(x$ALT))){
  #       x$type <- paste0('ins ', nchar(x$ALT)-3)
  #     } else if (grepl('del', as.character(x$ALT))){
  #       x$type <- paste0('del ', nchar(x$ALT)-3)
  #     } else if (orf_aa != orf2_aa){
  #       # Fixes probject of multiple hits.
  #       # x$type <- paste0(orf_aa, aa2, orf2_aa)
  #       x$type <- paste0(orf_aa[1], aa2[1], orf2_aa[1])
  #     } else {
  #       x$type <- 'silent'
  #     }
  #   }
  #   x
  # }))
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
  if(ref$wgs[i] == "y"){
  temp_mutation <- read.csv(paste0("./03_Mutation/", gsub(".fasta", "_snp.csv",gsub("./01_Sequence/","",ref$file_name[i]))))
  temp_mutation <- temp_mutation[,c(4,2,3,1)]
  temp_mutation = subset(temp_mutation, select = -c(X) )
  colnames(temp_mutation) <- c("POS", "REF", "ALT")
  temp_mutation$VSP <- ref$seq_name[i]
  mut <- rbind(mut,temp_mutation)
  }
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

# # For cat samples - Determine if there are more mutations that match USDA or our sequences.
# mut3 <- mut2
# mut3$match <- "Other Two Only"
# for(i in 1:length(mut3$match)){
#   if(mut3$VSP3508[i] == 'Yes'){
#   if(mut3$VSP3508[i] == mut3$usda.cat.1[i]){
#     mut3$match[i] <- "USDA"
#   }
#   if(mut3$VSP3508[i] == mut3$VSP3454[i]){
#     mut3$match[i] <- "Bushman"
#   }
#   if(mut3$VSP3508[i] == mut3$usda.cat.1[i] && mut3$VSP3508[i] == mut3$VSP3454[i]){
#     mut3$match[i] <- "Both"
#   }
#     if(mut3$match[i] == "Other Two Only"){
#       mut3$match[i] <- "only this one"
#     }
#   }
# }
# 
# # Determine the number of matches for each group.
# table(mut3$match)
# 
# # For cat samples - Determine if there are more mutations that match USDA or our sequences.
# mut4 <- mut2
# mut4$match <- "Other Two Only"
# for(i in 1:length(mut3$match)){
#   if(mut4$VSP3454[i] == 'Yes'){
#   if(mut4$VSP3454[i] == mut4$usda.cat.1[i]){
#     mut4$match[i] <- "USDA"
#   }
#   if(mut4$VSP3454[i] == mut4$VSP3508[i]){
#     mut4$match[i] <- "Bushman"
#   }
#   if(mut4$VSP3454[i] == mut4$usda.cat.1[i] && mut4$VSP3508[i] == mut4$VSP3454[i]){
#     mut4$match[i] <- "Both"
#   }
#     if(mut4$match[i] == "Other Two Only"){
#       mut4$match[i] <- "only this one"
#     }
#   }
# }
# 
# # Determine the number of matches for each group.
# library(ggvenn)
# 
# mut5 <- data.frame(matrix(NA, nrow = length(mut1$VSP), ncol = 2))
# colnames(mut5) <- c("sample", "mutation")
# mut5$sample <- mut1$sample
# mut5$mutation <- mut1$mutation
# 
# a <- subset(mut5, mut5$sample == "usda-cat-1")
# b <- subset(mut5, mut5$sample == "VSP3454")
# c <- subset(mut5, mut5$sample == "VSP3508")
# 
# x <- list(
#   USDA = a$mutation, 
#   Bushman_Recent = c$mutation,
#   Bushman_Previous = b$mutation
# )
# 
# ggvenn(
#   x, 
#   fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
#   stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE
# )

