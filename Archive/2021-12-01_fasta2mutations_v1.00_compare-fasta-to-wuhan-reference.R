# # Load libraries
# library(ape)
# library(ggplot2)

# Make the necessary directory architecture.
# if(!file.exists("Queries_Blast")){system("mkdir Queries_Blast")}
if(!file.exists("Results_Blast")){system("mkdir Results_Blast")}
if(!file.exists("Results_Align")){system("mkdir Results_Align")}
if(!file.exists("Results_SNPs")){system("mkdir Results_SNPs")}
if(!file.exists("Output_Sequences")){system("mkdir Output_Sequences")}

# Prepare the fasta files for the blast. ###################################################

# Write the fasta files to individual files.
# Pull the fasta files. ########################################################
# Open the data
file_original <- "2021-11-15_cat_with-artic-v4.fasta"
file_lineage <- "lineage_report.csv"
fas1 <- read.table(file_original, header = FALSE, sep = "", dec = ".")
names(fas1)[1] <- "fasta"
fas1$fasta <- as.character(fas1$fasta)

# Condense the fasta sequence to be on one line.
ref <- data.frame(matrix(NA, nrow = length(fas1$fasta), ncol = 5))
colnames(ref) <- c("seq_num", "seq_name", "fasta", "file_name", "wgs")

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
    ref$seq_name[seq_num] <- paste0(">seq_", as.character(seq_num), "_cat")
    ref$fasta[seq_num] <- fas1$fasta[i]
    ref$file_name[seq_num] <- paste0("seq_", as.character(seq_num), "_cat", ".fasta")
    ref$wgs[seq_num] <- "unknown"
    seq_num <- seq_num + 1
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

# Choose only the sequences that are longer than 29k nucleotides and separate sequences into individual files.
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
      file_name <- paste0("./Queries_Blast/", ref$file_name[((i+1)/2)])
      temp$fasta[1] <- fasta_name
      temp$fasta[2] <- seq
      write.table(temp, file_name, quote = F, col.names=F, row.names = F)

      # Run the command for blast.
      # Prepare the variables needed to perform the blast.
      fasta     <- ref$file_name[((i+1)/2)]
      results   <- gsub(".fasta", "_results.txt", fasta)
      database <- "2021-10-29_Human-WGS-Database/genomic.fna"
      command <- paste0("blastn -db ", database, " -query ./Queries_Blast/", fasta, " -out ./Results_Blast/", results)
      print(command)
      # Write a bash script that will run the blast.
      filepath01 <- '.'
      p <- c('#!/bin/bash', command)
      writeLines(p, file.path(filepath01, 'temp.script'))
      system(paste0('chmod 755 ', file.path(filepath01, 'temp.script')))
      comm <- paste0(file.path(filepath01, 'temp.script'))
      # Run the bash script.
      system(comm)
    }else{
      ref$wgs[((i+1)/2)] <- "n"
    }
  }
}

# Add the top scoring human WGS accession to the ref data frame.
ref$human_accession <- NA
for(i in 1:length(ref$seq_num)){
  if(ref$wgs[i] == "y"){
    file <- paste0("./Results_Blast/",gsub(".fasta", "_results.txt", ref$file_name[i]))
    temp <- read.delim(file)
    temp <- data.frame(temp)
    temp$BLASTN.2.6.0. <- as.character(temp$BLASTN.2.6.0.)
    # Determine which is the first line to be accession numbers.
    accession <- NA
    for(j in 1:length(temp$BLASTN.2.6.0.)){
      word <- gsub(" .*","",temp$BLASTN.2.6.0.[j])
      if(word == 'Sequences'){
        accession <- gsub(" .*","",temp$BLASTN.2.6.0.[j+1])
        break
      }
    }
    ref$human_accession[i] <- accession
    if (i%%10 == 0) {
      out_text <- paste(as.character(i), " sequences accession found", "")
      print(out_text, quote = FALSE)
    }
  }
}

# Pull the closest matching human wgs from ncbi and place it in the Output_Sequences directory.
ref$percent_identity <- NA
ref$distance <- NA
ref$snp_distance <- NA
all_seq <- data.frame(matrix(NA, nrow = 0, ncol = 1))
colnames(all_seq) <- c("fasta")
for(i in 1:length(ref$seq_num)){
  if(ref$wgs[i] == 'y'){
    # Copy all wgs query files from the Queries directory to the Output_Sequences directory.
    file.copy(from = paste0("./Queries_Blast/", ref$file_name[i]),
              to = paste0("./Output_Sequences/", ref$file_name[i]))

    # Download the closest sequence that came from the blast search.
    accession <- ref$human_accession[i]
    file <- paste0("./Output_Sequences/", gsub("cat", "human",ref$file_name[i]))
    command <- paste0("esearch -db nucleotide -query '", accession, "' | efetch -format fasta > ", file)
    # Write a bash script that will run the blast.
    filepath01 <- '.'
    p <- c('#!/bin/bash', "export PATH=${PATH}:${HOME}/edirect", command)
    writeLines(p, file.path(filepath01, 'temp.script'))
    system(paste0('chmod 755 ', file.path(filepath01, 'temp.script')))
    comm <- paste0(file.path(filepath01, 'temp.script'))
    # Run the bash script.
    system(comm)

    # Rename the cat fasta.
    fasta <- read.delim(paste0("./Output_Sequences/", ref$file_name[i]))
    fasta <- data.frame(fasta)
    colnames(fasta) <- paste0(gsub("cat", "", ref$seq_name[i]), gsub(">", "",ref$fasta[i]))
    write.table(fasta, paste0("./Output_Sequences/", ref$file_name[i]), quote = F, col.names=T, row.names = F)
    # Add the cat fasta to the all_seqs file.
    temp1 <- data.frame(matrix(NA, nrow = 1, ncol = 1))
    fasta2 <- fasta
    colnames(fasta2) <- "fasta"
    colnames(temp1) <- c("fasta")
    temp1$fasta[1] <- colnames(fasta)
    temp1 <- rbind(temp1,fasta2)
    all_seq <- rbind(all_seq,temp1)

    # Rename the fasta files to prepared for aligning.
    # Rename the human fasta.
    fasta <- read.delim(paste0("./Output_Sequences/", gsub("cat", "human",ref$file_name[i])))
    fasta <- data.frame(fasta)
    colnames(fasta) <- paste0(gsub("cat", "human", ref$seq_name[i]), "_", accession)
    write.table(fasta, paste0("./Output_Sequences/", gsub("cat", "human",ref$file_name[i])), quote = F, col.names=T, row.names = F)
    # Add the human fasta to the all_seqs file.
    temp2 <- data.frame(matrix(NA, nrow = 1, ncol = 1))
    fasta2 <- fasta
    colnames(fasta2) <- "fasta"
    colnames(temp2) <- c("fasta")
    temp2$fasta[1] <- colnames(fasta)
    temp2 <- rbind(temp2,fasta2)
    all_seq <- rbind(all_seq,temp2)

    # Determine the percent identity.
    file1 <- paste0("./Output_Sequences/", gsub("cat", "human",ref$file_name[i]))
    file2 <- paste0("./Output_Sequences/", ref$file_name[i])
    command <- paste0("blastn -query ", file1, " -subject ", file2, " -task blastn -outfmt '6 pident score' > temp.csv")
    # Write a bash script that will run the blast.
    filepath01 <- '.'
    p <- c('#!/bin/bash', "export PATH=${PATH}:${HOME}/edirect", command)
    writeLines(p, file.path(filepath01, 'temp.script'))
    system(paste0('chmod 755 ', file.path(filepath01, 'temp.script')))
    comm <- paste0(file.path(filepath01, 'temp.script'))
    # Run the bash script.
    system(comm)
    temp <- read.delim('temp.csv', header = F)
    ref$percent_identity[i] <- temp$V1[1]


    # Determine the distance between the two samples.
    # Perform a MSA using muscle.
    file3 <- paste0("./Results_Blast/", gsub(".fasta", "_align.fasta",ref$file_name[i]))
    command <- paste0("./muscle3.8.31_i86linux64 -profile -in1 ", file1, " -in2 ", file2, " -out ", file3)
    # Write a bash script that will run the blast.
    filepath01 <- '.'
    p <- c('#!/bin/bash', "export PATH=${PATH}:${HOME}/edirect", command)
    writeLines(p, file.path(filepath01, 'temp.script'))
    system(paste0('chmod 755 ', file.path(filepath01, 'temp.script')))
    comm <- paste0(file.path(filepath01, 'temp.script'))
    # Run the bash script.
    system(comm)
    # Calculate the pairwise distance between the two samples.
    seq <- ape::read.dna(file3, format = "fasta")
    dis <- ape::dist.dna(seq)
    if(dis[1] < 0.0001){
      ref$distance[i] <- 0
    }else{
      ref$distance[i] <- round(dis[1],4)
    }
    
    
    # Determine the SNP distances from the MSA.
    align_char <- as.character(seq)
    align_df <- data.frame(align_char)
    l <- length(align_df)
    # # colnames(align_df) <- as.character(seq(1:length(align_df)))
    # nt <- align_df[100]
    
    # Transpose the data frame so that it has 1 row for each position.
    align_df <- as.data.frame(t(as.matrix(align_df)))
    
    # Change the row index so that it allows it is the numbers and column names so that they are uniform.
    rownames(align_df) <- seq(1:l)
    col_name <- colnames(align_df)
    colnames(align_df) <- c("seq1", "seq2")
    align_df$seq1 <- as.character(align_df$seq1)
    align_df$seq2 <- as.character(align_df$seq2)
    # align_df <- tail(align_df,-29606) # Remove the last 29k lines
    
    
    # For coding add an artificial mismatch.
    # align_df$seq2[200] <- "0000"
    # Select only the lines that are not aligned.
    align_mismatch <- subset(align_df, align_df$seq1 != align_df$seq2)
    
    # Remove the first and last group of lines that might results from sequencing (threshould is 50 continuous mathcing bp).
    cutoff <- 50
    if(length(align_mismatch$seq1) >= 1){
      align_mismatch$filter <- "y"
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
    ref$snp_distance[i] <- length(align_mismatch2$filter)
    
    # Save the SNP file.
    align_mismatch3 <- align_mismatch2[ , ! names(align_mismatch2) %in% c("filter")]
    colnames(align_mismatch3) <- col_name
    write.csv(align_mismatch3, paste0("./Results_SNPs/", gsub(".fasta", "_snp.csv",ref$file_name[i])))
  }
  if (i%%10 == 0) {
    out_text <- paste(as.character(i), " closest fasta sequence processed", "")
    print(out_text, quote = FALSE)
  }
}
write.table(all_seq, "./Output_Sequences/all_seqs.fasta", quote = F, col.names=F, row.names = F)

# Determine the lineage for the input samples.
pangolin_command <- paste0("pangolin ", file_original)
p <- c('#!/bin/bash', 'source /home/marques/miniconda3/etc/profile.d/conda.sh', 'conda activate pangolin', pangolin_command)
writeLines(p, file.path('./temp.script'))
system(paste0('chmod 755 ', file.path('./temp.script')))

comm <- paste0('./', 'temp.script')
system(comm)

# Add the lineage that was detected for each sample.
lin1 <- read.csv(file_lineage)
ref$lineage <- NA
for(i in 1:length(ref$seq_num)){
  temp <- data.frame((subset(lin1, lin1$taxon == gsub(">","", ref$fasta[i]))))
  ref$lineage[i] <- as.character(temp$lineage[1])
}

# Use the new method for labelling the samples.
# ref <- read.csv("reference.csv")
# ref$distance[1] <- 0
# ref1 <- ref

ref1 <- ref
ref1 <- subset(ref1, !is.na(ref1$distance))
for(i in 1:length(ref1$lineage)){
  if(grepl("AY", ref1$lineage[i]) | ref1$lineage[i] == "B.1.617.2"){ref1$lineage[i] <- "Delta"}
  if(ref1$lineage[i] == "B.1.1.7"){ref1$lineage[i] <- "Alpha"}
  if(ref1$lineage[i] == "B.1.526"){ref1$lineage[i] <- "Iota"}
  if(ref1$lineage[i] == "B.1.621"){ref1$lineage[i] <- "Mu"}
  if(ref1$lineage[i] == "B.1.1.7"){ref1$lineage[i] <- "Gamma"}
}


ref1$lineage_original <- ref1$lineage
ref1$lineage <- "Other"
important_lineage <- data.frame(table(ref1$lineage_original))
important_lineage <- subset(important_lineage, important_lineage$Freq > 2) # Choose only the lineage with 3 or more occurances found.
important_lineage <- as.character(important_lineage$Var1)
# important_lineage <- important_lineage[important_lineage != "B.1.1.298"]
for(i in 1:length(ref1$lineage)){
  for(j in 1:length(important_lineage)){
    if(ref1$lineage_original[i] == important_lineage[j]){
      ref1$lineage[i] <- important_lineage[j]
    }
  }
}

# Perform a one-way anova test 
stat1 <- data.frame(matrix(NA, nrow = length(ref1$lineage), ncol = 2))
colnames(stat1) <- c("lineage", "distance")
stat1$lineage <- ref1$lineage
stat1$distance <- ref1$distance
anova1 <- aov(distance~lineage, stat1)
summary(anova1)
p_value <- summary(anova1)[[1]][["Pr(>F)"]][1]
p_value <- round(p_value,4)
if(p_value <0.0001){
  p_value <- "<0.0001"
}else{
  p_value <- as.character(p_value)
}
TukeyHSD(anova1)

# Perform anova on the reduced groups.
stat3 <- stat1
stat3$lineage_original <- stat3$lineage
stat3$lineage <- "Other"
for(i in 1:length(stat3$lineage)){
  temp <- "B.1"
  if(stat3$lineage_original[i] == temp){stat3$lineage[i] <- temp}
  temp <- "B.1.1.7"
  if(stat3$lineage_original[i] == temp){stat3$lineage[i] <- temp}
  temp <- "AY.29"
  if(stat3$lineage_original[i] == temp){stat3$lineage[i] <- temp}
  temp <- "B.1.617.2"
  if(stat3$lineage_original[i] == temp){stat3$lineage[i] <- temp}
  temp <- "B.1.2"
  if(stat3$lineage_original[i] == temp){stat3$lineage[i] <- temp}
}
anova2 <- aov(distance~lineage, stat3)
summary(anova2)
# p_value <- summary(anova2)[[1]][["Pr(>F)"]][1]

# Make a data frame for graphing using Rick's desired output butin pairwise distance.
# Remove the rows where there was not a nearest human genome found.
data1 <- ref1
data2 <- subset(data1, !is.na(data1$human_accession))

# Make a new column where there are categories for the different distances.
data2$dis_cat <- NA
for(i in 1:length(data2$seq_num)){
  if(data2$distance[i] == 0){
    data2$dis_cat[i] <- "0"}
  if(data2$distance[i] > 0 && data2$distance[i] <= 0.0001 ){
    data2$dis_cat[i] <- "0-0.0001"}
  if(data2$distance[i] > 0.0001 && data2$distance[i] <= 0.0005 ){
    data2$dis_cat[i] <- "0.0001-0.0005"}
  if(data2$distance[i] > 0.0005 && data2$distance[i] <= 0.001 ){
    data2$dis_cat[i] <- "0.0005-0.001"}
  if(data2$distance[i] > 0.001){
    data2$dis_cat[i] <- ">0.001"}
}

# dis <- unique(data2$dis_cat) # Use this for automatic determination of the distance categories.
dis <- c("0", "0-0.0001", "0.0001-0.0005", "0.0005-0.001", ">0.001")

lin <- unique(data2$lineage)
data3 <- data.frame(matrix(NA, nrow = length(lin)*length(dis), ncol = 3))
colnames(data3) <- c("lineage", "distance", "count")
data3$count <- 0
data3$lineage <- rep(lin, length(dis))
data3$distance <- rep(dis, 1, each = length(lin))

# Populate the data frame. 
for(i in 1:length(data3$lineage)){
  temp <- subset(data2, data2$lineage == data3$lineage[i])
  temp <- subset(temp, temp$dis_cat == data3$distance[i])
  data3$count[i] <- length(temp$seq_num)
  
}

# Make the x axis labels as factors.
data3$distance <- factor(data3$distance, levels = unique(data3$distance))

colors <- rainbow(length(important_lineage) + 1)
# Plot the stacked bar graph.
plot <- ggplot(data3, aes(fill=lineage, y=count, x=distance)) + 
  geom_bar(position="stack", stat="identity") + 
  xlab("Distance to Nearest Human")+
  ylab("Number of Infected Animals") +
  ggtitle(label = "Distance Between Feline and Human SARS-CoV-2") + 
  guides(fill=guide_legend(title="Lineage"))+ 
  scale_y_continuous(expand = c(0, 0)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                    panel.background = element_blank(), axis.line = element_line(colour = "black"))
plot
ggsave("distance-plot_v2.pdf", plot)

# Save the reference file.
write.csv(ref, "reference2.csv", row.names = F)



# Make a data frame for graphing using Rick's desired output but in SNP distance.
# Remove the rows where there was not a nearest human genome found.
data1 <- ref1
row.names(data1) <- NULL
data1$snp_distance[67] <- 1 # This is because the NJ cat has 1 snp distance to the human cat that we sequenced.This sequence is not on ncbi yet so it cannot be included in the blast.
data2 <- subset(data1, !is.na(data1$human_accession))

# Make a new column where there are categories for the different distances.
data2$dis_cat <- NA
for(i in 1:length(data2$seq_num)){
  if(data2$snp_distance[i] > -1 && data2$distance[i] <= 1 ){
    data2$dis_cat[i] <- "0-1"}
  if(data2$snp_distance[i] > 1 && data2$distance[i] <= 4 ){
    data2$dis_cat[i] <- "2-4"}
  if(data2$snp_distance[i] > 4 && data2$distance[i] <= 7 ){
    data2$dis_cat[i] <- "5-7"}
  if(data2$snp_distance[i] > 7 && data2$distance[i] <= 10 ){
    data2$dis_cat[i] <- "8-10"}
  if(data2$snp_distance[i] > 10 && data2$distance[i] <= 13 ){
    data2$dis_cat[i] <- "11-13"}
  if(data2$snp_distance[i] > 13 && data2$distance[i] <= 16 ){
    data2$dis_cat[i] <- "14-16"}
  if(data2$snp_distance[i] > 16){
    data2$dis_cat[i] <- ">16"}
}

# dis <- unique(data2$dis_cat) # Use this for automatic determination of the distance categories.
dis <- c("0-1", "2-4", "5-7", "8-10", "11-13", "14-16", ">16")

lin <- unique(data2$lineage)
data3 <- data.frame(matrix(NA, nrow = length(lin)*length(dis), ncol = 3))
colnames(data3) <- c("lineage", "distance", "count")
data3$count <- 0
data3$lineage <- rep(lin, length(dis))
data3$distance <- rep(dis, 1, each = length(lin))

# Populate the data frame. 
for(i in 1:length(data3$lineage)){
  temp <- subset(data2, data2$lineage == data3$lineage[i])
  temp <- subset(temp, temp$dis_cat == data3$distance[i])
  data3$count[i] <- length(temp$seq_num)
  
}

# Make the x axis labels as factors.
data3$distance <- factor(data3$distance, levels = unique(data3$distance))

colors <- rainbow(length(important_lineage) + 1)
# Plot the stacked bar graph.
plot <- ggplot(data3, aes(fill=lineage, y=count, x=distance)) + 
  geom_bar(position="stack", stat="identity") + 
  xlab("SNP Distance to Nearest Human")+
  ylab("Number of Animal-Derived Genomes") +
  ggtitle(label = "Distance Between Cat and Human SARS-CoV-2") + 
  guides(fill=guide_legend(title="Lineage"))+ 
  scale_y_continuous(expand = c(0, 0)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
plot
ggsave("snp-distance-plot_v4.pdf", plot)

# Save the reference file.
write.csv(ref, "reference2.csv", row.names = F)


