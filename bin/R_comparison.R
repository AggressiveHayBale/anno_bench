#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

###########
#Functions#
###########

##Ampir-read faa
read_faa<-function (file = NULL) 
{
  faa_lines <- readLines(file)
  seq_name_index <- grep(">", faa_lines)
  seq_name <- gsub(">", "", faa_lines[seq_name_index])
  seq_aa_start_index <- seq_name_index + 1
  seq_aa_end_index <- c(seq_name_index, length(faa_lines) + 
                          1)[-1] - 1
  seq_aa <- rep(NA, length(seq_name_index))
  for (i in seq_along(seq_name_index)) {
    seq_aa_start <- seq_aa_start_index[i]
    seq_aa_end <- seq_aa_end_index[i]
    seq_aa[i] <- gsub("[[:space:]]", "", paste(faa_lines[seq_aa_start:seq_aa_end], 
                                               collapse = ""))
  }
  data.frame(seq_name, seq_aa, stringsAsFactors = FALSE)
}

###########
#Data load#
###########

library(stringr)
library(readr)
library(dplyr)
library(tidyr)


id <- args[1]

prokka_faa<- args[2]
prokka_gbk<- args[3]

bakta_faa<- args[4]
bakta_gbff<- args[5]

eggnog_faa<- args[6]
eggnog_gff<- args[7]

pgap_faa<- args[8]
pgap_gff<- args[9]

fasta_type <- args[10]


 
########
#Prokka#
########

# File read
annotation <- suppressMessages(read_delim(prokka_gbk,
                                          delim = "\t", escape_double = FALSE,
                                          col_names = FALSE, comment = c("#"), trim_ws = TRUE))
fasta_data <- read_faa(prokka_faa)

#Cleanup
annotation <- invisible(na.omit(annotation))
annotation<- annotation[!(annotation$X3=="region" | annotation$X3=="gene"),]
annotation <- separate(annotation, X9, c(NA, "ID"), remove = FALSE, "ID=")
annotation$ID <- gsub(";.*","",annotation$ID)
annotation$product <- str_extract(annotation$X9,pattern = "product=(.*?;)") 
annotation$product <- gsub("product=", "", annotation$product)
annotation$product <- gsub(";", "", annotation$product)
annotation$product[annotation$product==""] <- NA

fasta_data$seq_name <- sub(" .*", "", fasta_data$seq_name )
# Merging
annotation_merged <- merge(annotation,fasta_data, by.x= "ID", by.y = "seq_name", all=TRUE)
###tRNA and rRNA does not have sequences - only CDS

# Cleanup
colnames(annotation_merged)[colnames(annotation_merged) == 'seq_aa']<- c("sequence")

prokka<- annotation_merged

########
#Bakta #
########

# File read
annotation <- suppressMessages(read_delim(bakta_gbff,
                                          delim = "\t", escape_double = FALSE,
                                          col_names = FALSE, comment = c("#"), trim_ws = TRUE))
fasta_data <- read_faa(prokka_faa)


#Cleanup
annotation <- invisible(na.omit(annotation))
annotation<- annotation[!(annotation$X3=="region" | annotation$X3=="gene"),]
annotation <- separate(annotation, X9, c(NA, "ID"), remove = FALSE, "ID=")
annotation$ID <- gsub(";.*","",annotation$ID)
annotation$product <- str_extract(annotation$X9,pattern = "Name=(.*?;)") 
annotation$product <- gsub("Name=", "", annotation$product)
annotation$product <- gsub(";", "", annotation$product)
annotation$product[annotation$product==""] <- NA

# Merging
fasta_data$seq_name <- sub(" .*", "", fasta_data$seq_name )
annotation_merged <- merge(annotation,fasta_data, by.x= "ID", by.y = "seq_name", all=TRUE)
# Cleanup
colnames(annotation_merged)[colnames(annotation_merged) == 'seq_aa']<- c("sequence")

bakta <- annotation_merged


######
#PGAP#
######

annotation <- suppressMessages(read_delim(pgap_gff,
                                          delim = "\t", escape_double = FALSE,
                                          col_names = FALSE, comment = c("#"), trim_ws = TRUE))
fasta_data <- read_faa(pgap_faa)

annotation<- annotation[!(annotation$X3=="region" | annotation$X3=="gene"),]

annotation <- separate(annotation, X9, c(NA, "ID"), remove = FALSE, "ID=cds-")
annotation$ID <- gsub(";.*","",annotation$ID)

annotation$product <- str_extract(annotation$X9,pattern = "product=(.*?;)") 
annotation$product <- gsub("product=", "", annotation$product)
annotation$product <- gsub(";", "", annotation$product)

# Merging
fasta_data$seq_name <- sub(" .*", "", fasta_data$seq_name )
fasta_data$seq_name <- sapply(str_split(fasta_data$seq_name,"\\|"), function(x) (x[3]))
annotation_merged <- merge(annotation,fasta_data, by.x= "ID", by.y = "seq_name", all=TRUE)
# Cleanup

colnames(annotation_merged)[colnames(annotation_merged) == 'seq_aa']<- c("sequence")

pgap<-annotation_merged

########
#Eggnog#
########

annotation <- suppressMessages(read_delim(eggnog_gff,
                                          delim = "\t", escape_double = FALSE,
                                          col_names = FALSE, comment = c("#"), trim_ws = TRUE))

fasta_data <- read_faa(eggnog_faa)
# Merging
fasta_data$seq_name <- sapply(str_split(fasta_data$seq_name," "), function(x) (x[1]))
annotation_merged <- merge(annotation,fasta_data, by=0, all=TRUE)
# Cleanup
colnames(annotation_merged)[colnames(annotation_merged) == 'seq_aa']<- c("sequence")

#Remove star from end of sequence 
annotation_merged$sequence = substring(annotation_merged$sequence,1, nchar(annotation_merged$sequence)-1)

em_Preferred_name <- str_extract(annotation_merged$X9,pattern = "em_Preferred_name=(.*?;)") 
em_Preferred_name <- gsub("em_Preferred_name=", "", em_Preferred_name)
em_Preferred_name <- gsub(";", "", em_Preferred_name)
em_Preferred_name[em_Preferred_name==""] <- NA
annotation_merged$product<- em_Preferred_name
eggnog<-annotation_merged

##########
#Analysis#
##########

#Overall gene count
prokka_gene_ct <- nrow(prokka)
bakta_gene_ct <-nrow(bakta)
pgap_gene_ct <- nrow(pgap)
eggnog_gene_ct  <- nrow(eggnog)

#Hypothetical proteins
prokka_hyps <- nrow(prokka[prokka$product=="hypothetical protein",])
bakta_hyps <- nrow(bakta[bakta$product=="hypothetical protein",])
pgap_hyps <- nrow(pgap[pgap$product=="hypothetical protein",])
eggnog_hyps <- nrow(pgap[pgap$product==NA,])
#eggnog_hyps<- sum(is.na(eggnog$X8))+nrow(eggnog[eggnog$X8=="-" & eggnog$X9=="-",]) + sum(str_detect(na.omit(eggnog$X8), "unknown"))

#Percent of hyps 
prokka_perc_hyps <- (prokka_hyps/prokka_gene_ct)*100
bakta_perc_hyps <- (bakta_hyps/bakta_gene_ct)*100
pgap_perc_hyps <- (pgap_hyps/pgap_gene_ct)*100
eggnog_perc_hyps<- (eggnog_hyps/eggnog_gene_ct)*100



##Check what could be hyps
#library(wordcloud)
#library(tm)
#corpus <- Corpus(VectorSource(na.omit(eggnog$X8)))
#DTM <- TermDocumentMatrix(corpus)
#mat <- as.matrix(DTM)
#f <- sort(rowSums(mat),decreasing=TRUE)
#dat <- data.frame(word = names(f),freq=f)
#wordcloud(words = dat$word, freq = dat$freq,colors=brewer.pal(8, "Dark2"))

#Type of sequence
prokka_type <-as.list(table(prokka$X3))
bakta_type <- as.list(table(bakta$X3))
pgap_type <- as.list(table(pgap$X3))

prokka_type <-as.list(table(prokka$X3))
bakta_type <- as.list(table(bakta$X3))
pgap_type <- as.list(table(pgap$X3))

prokka_CDS<-ifelse(is.null(prokka_type$CDS),0,prokka_type$CDS )
bakta_CDS<- ifelse(is.null(bakta_type$CDS),0,bakta_type$CDS )
pgap_CDS<- ifelse(is.null(pgap_type$CDS),0,pgap_type$CDS )

prokka_rrna<- ifelse(is.null(prokka_type$rRNA),0,prokka_type$rRNA )
bakta_rrna<-ifelse(is.null(bakta_type$rRNA),0,bakta_type$rRNA )
pgap_rrna<- ifelse(is.null(pgap_type$rRNA),0,pgap_type$rRNA )


prokka_trna<-ifelse(is.null(prokka_type$tRNA),0,pgap_type$tRNA )
bakta_trna<-ifelse(is.null(bakta_type$tRNA),0,bakta_type$tRNA )
pgap_trna<-ifelse(is.null(pgap_type$tRNA),0,pgap_type$tRNA )


prokka_tmrna<-ifelse(is.null(prokka_type$tmRNA),0,pgap_type$tmRNA )
bakta_tmrna<-ifelse(is.null(bakta_type$tmRNA),0,bakta_type$tmRNA )
pgap_tmrna<-ifelse(is.null(pgap_type$tmRNA),0,pgap_type$tmRNA )

#Sequence 

prokka_seq_lenght<-median(nchar(na.omit(prokka$sequence)))
bakta_seq_lenght<-median((nchar(na.omit(bakta$sequence))))
pgap_seq_lenght<-median(nchar(na.omit(pgap$sequence)))
eggnog_seq_lenght<-median((nchar(na.omit(eggnog$sequence))))

#Sequence comparison 

prokka_bakta_intersect<-length(intersect(na.omit(prokka$name),na.omit(bakta$name)))
prokka_eggnog_intersect<-length(intersect(na.omit(prokka$name),na.omit(eggnog$name)))
prokka_pgap_intersect<-length(intersect(na.omit(prokka$name),na.omit(pgap$name)))

bakta_eggnog_intersect<-length(intersect(na.omit(bakta$name),na.omit(eggnog$name)))
bakta_pgap_intersect<-length(intersect(na.omit(bakta$name),na.omit(pgap$name)))

eggnog_bakta_intersect<-length(intersect(na.omit(eggnog$name),na.omit(bakta$name)))
eggnog_pgap_intersect<-length(intersect(na.omit(eggnog$name),na.omit(pgap$name)))

pgap_bakta_intersect<-length(intersect(na.omit(pgap$name),na.omit(bakta$name)))
pgap_eggnog_intersect<-length(intersect(na.omit(pgap$name),na.omit(eggnog$name)))


#Create summary
col_labels_genecount<-c(c("Prokka gene count","Bakta gene count","PGAP gene count","EGGnog gene count"))
col_labels_hyps<-c(c("Prokka hypothetical genes","Bakta hypothetical genes","PGAP hypothetical genes","EGGnog hypothetical genes"))
col_labels_perc_hyps<-c(c("Prokka hypothetical genes(%)","Bakta hypothetical genes(%)","PGAP hypothetical genes(%)","EGGnog hypothetical genes(%)"))
col_labels_CDS<- c(c("Prokka CDS","Bakta CDS","PGAP CDS"))
col_labels_rrna<- c(c("Prokka rRNA","Bakta rRNA","PGAP rRNA"))
col_labels_trna<- c(c("Prokka tRNA","Bakta tRNA","PGAP tRNA"))
col_labels_tmrna<- c(c("Prokka tmRNA","Bakta tmRNA","PGAP tmRNA"))
summary_df<- data.frame(id, fasta_type, prokka_gene_ct,bakta_gene_ct,pgap_gene_ct,eggnog_gene_ct,
                        prokka_hyps,bakta_hyps,pgap_hyps,eggnog_hyps,
                        prokka_perc_hyps,bakta_perc_hyps,pgap_perc_hyps,eggnog_perc_hyps,
                        prokka_CDS,bakta_CDS,pgap_CDS,
                        prokka_rrna,bakta_rrna,pgap_rrna,
                        prokka_trna,bakta_trna,pgap_trna,                        
                        prokka_tmrna,bakta_tmrna,pgap_tmrna
                        )
colnames(summary_df) <- c("id","type",col_labels_genecount,col_labels_hyps,col_labels_perc_hyps,col_labels_CDS,col_labels_rrna,col_labels_trna,col_labels_tmrna)

write_csv(summary_df,"comparison.csv")
