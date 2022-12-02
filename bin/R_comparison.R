#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(gggenomes)
library(ampir)
library(stringr)

prokka_gbk<- args[1]
prokka_faa<- args[2]

bakta_gbff<- args[3]
bakta_faa<- args[4]

pgap_gff<- args[5]
pgap_faa<- args[6]

eggnog_gff<- args[7]
eggnog_faa<- args[8]

########
#Prokka#
########

# File read
annotation <- read_gbk(prokka_gbk)
fasta_data <- read_faa(prokka_faa)
# Merging
fasta_data$seq_name <- sub(" .*", "", fasta_data$seq_name )
annotation_merged <- merge(annotation,fasta_data, by.x= "locus_tag", by.y = "seq_name", all=TRUE)
# Cleanup
annotation_merged<- annotation_merged[!(annotation_merged$type=="region" | annotation_merged$type=="gene"),]
colnames(annotation_merged)[colnames(annotation_merged) == 'seq_aa']<- c("sequence")

prokka<- annotation_merged

########
#Bakta #
########

# File read
annotation <- read_gbk(bakta_gbff)
fasta_data <- read_faa(prokka_faa)
# Merging
fasta_data$seq_name <- sub(" .*", "", fasta_data$seq_name )
annotation_merged <- merge(annotation,fasta_data, by.x= "locus_tag", by.y = "seq_name", all=TRUE)
# Cleanup
annotation_merged<- annotation_merged[!(annotation_merged$type=="region" | annotation_merged$type=="gene"),]
colnames(annotation_merged)[colnames(annotation_merged) == 'seq_aa']<- c("sequence")

bakta <- annotation_merged


######
#PGAP#
######

annotation <- read_gbk(pgap_gff)
fasta_data <- read_faa(pgap_faa)
# Merging
fasta_data$seq_name <- sub(" .*", "", fasta_data$seq_name )
fasta_data$seq_name <- sapply(str_split(fasta_data$seq_name,"\\|"), function(x) (x[3]))
annotation_merged <- merge(annotation,fasta_data, by.x= "locus_tag", by.y = "seq_name", all=TRUE)
# Cleanup
annotation_merged<- annotation_merged[!(annotation_merged$type=="region" | annotation_merged$type=="gene"),]
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
annotation_merged <- merge(annotation,fasta_data, by.x= "X1", by.y = "seq_name", all=TRUE)
# Cleanup
colnames(annotation_merged)[colnames(annotation_merged) == 'seq_aa']<- c("sequence")

#Remove star from end of sequence 
annotation_merged$sequence = substring(annotation_merged$sequence,1, nchar(annotation_merged$sequence)-1)
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
nrow(prokka[prokka$product=="hypothetical protein",])
nrow(bakta[bakta$product=="hypothetical protein",])
nrow(pgap[pgap$product=="hypothetical protein",])
sum(is.na(eggnog$X8))+nrow(eggnog[eggnog$X8=="-" & eggnog$X9=="-",]) + sum(str_detect(na.omit(eggnog$X8), "unknown"))

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
prokka_type <-table(prokka$type)
bakta_type <- table(bakta$type)
pgap_type <- table(pgap$type)

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





