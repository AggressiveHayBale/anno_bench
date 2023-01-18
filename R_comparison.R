
#BiocManager::install("ggtree")
#install.packages("devtools")
#devtools::install_github("thackl/thacklr")
#devtools::install_github("thackl/gggenomes")
library(gggenomes)
#install.packages("ampir")
library(ampir)
library(stringr)
########
#Prokka#
########

setwd("/mnt/15cf8305-a1bc-486c-8f97-fe85e5f4d988/gits/anno_bench/results_old/GCF_0036971652/prokka")
# File read
annotation <- read_gbk("GCF_0036971652_vanilla_prokka.gbk")
fasta_data <- read_faa("GCF_0036971652_vanilla_prokka.faa")
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

setwd("/mnt/15cf8305-a1bc-486c-8f97-fe85e5f4d988/gits/anno_bench/results_old/GCF_0036971652/bakta/")
# File read
annotation <- read_gbk("GCF_0036971652_vanilla_bakta.gbff")
fasta_data <- read_faa("GCF_0036971652_vanilla_bakta.faa")
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

setwd("/mnt/15cf8305-a1bc-486c-8f97-fe85e5f4d988/gits/anno_bench/results_old/GCF_003697165.2/pgap/")

annotation <- read_gbk("annot.gbk")
fasta_data <- read_faa("annot.faa")
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
setwd("/mnt/15cf8305-a1bc-486c-8f97-fe85e5f4d988/gits/anno_bench/results_old/GCF_0036971652/eggnog/")

annotation <- suppressMessages(read_delim("GCF_0036971652_vanilla_eggnog.emapper.annotations",
                                         delim = "\t", escape_double = FALSE,
                                         col_names = FALSE, comment = c("#"), trim_ws = TRUE))

fasta_data <- read_faa("GCF_0036971652_vanilla_eggnog.emapper.genepred.fasta")
# Merging
fasta_data$seq_name <- sapply(str_split(fasta_data$seq_name," "), function(x) (x[1]))
annotation_merged <- merge(annotation,fasta_data, by.x= "X1", by.y = "seq_name", all=TRUE)
# Cleanup
colnames(annotation_merged)[colnames(annotation_merged) == 'seq_aa']<- c("sequence")

#Remove star from end of sequence 
annotation_merged$sequence = substring(annotation_merged$sequence,1, nchar(annotation_merged$sequence)-1)
eggnog<-annotation_merged




nrow(prokka)
nrow(bakta)
nrow(pgap)
nrow(eggnog)

nrow(prokka[prokka$product=="hypothetical protein",])/nrow(prokka)
nrow(bakta[bakta$product=="hypothetical protein",])/nrow(bakta)
nrow(pgap[pgap$product=="hypothetical protein",])/nrow(pgap)
sum(is.na(eggnog$X8))+nrow(eggnog[eggnog$X8=="-" & eggnog$X9=="-",]) + sum(str_detect(na.omit(eggnog$X8), "unknown")) /nrow(eggnog)

##Check what could be hyps
#library(wordcloud)
#library(tm)
#corpus <- Corpus(VectorSource(na.omit(eggnog$X8)))
#DTM <- TermDocumentMatrix(corpus)
#mat <- as.matrix(DTM)
#f <- sort(rowSums(mat),decreasing=TRUE)
#dat <- data.frame(word = names(f),freq=f)
#wordcloud(words = dat$word, freq = dat$freq,colors=brewer.pal(8, "Dark2"))
###
sum(is.na(eggnog$X8))+nrow(eggnog[eggnog$X8=="-" & eggnog$X9=="-",]) + sum(str_detect(na.omit(eggnog$X8), "unknown")) /nrow(eggnog)
#Check
(nrow(pgap[pgap$product=="hypothetical protein",])+nrow(bakta)-nrow(pgap))/nrow(bakta)

table(prokka$type)
table(bakta$type)
table(pgap$type)


#Sequence 

mean(nchar(na.omit(prokka$sequence)))
median(nchar(na.omit(prokka$sequence)))
mean(nchar(na.omit(bakta$sequence)))
median(nchar(na.omit(pgap$sequence)))
mean((nchar(na.omit(pgap$sequence))))
median((nchar(na.omit(pgap$sequence))))
mean((nchar(na.omit(eggnog$sequence))))
median((nchar(na.omit(eggnog$sequence))))

length(intersect(na.omit(prokka$name),na.omit(bakta$name)))
length(setdiff(na.omit(prokka$name),na.omit(bakta$name)))

length(intersect(na.omit(prokka$name),na.omit(pgap$name)))
length(setdiff(na.omit(prokka$name),na.omit(pgap$name)))

length(intersect(na.omit(prokka$sequence),na.omit(bakta$sequence)))
length(setdiff(na.omit(prokka$sequence),na.omit(bakta$sequence)))


length(intersect(na.omit(prokka$sequence),na.omit(pgap$sequence)))
length(setdiff(na.omit(prokka$sequence),na.omit(pgap$sequence)))

length(intersect(na.omit(pgap$sequence),na.omit(eggnog$sequence)))
length(setdiff(na.omit(pgap$sequence),na.omit(eggnog$sequence)))

prokka_gene_name<-sapply(str_split(na.omit(prokka$name),"_"), function(x) (x[1]))
bakta_gene_name<-sapply(str_split(na.omit(bakta$name),"_"), function(x) (x[1]))
eggnog_gene_name<-sapply(str_split(na.omit(eggnog$X9),"_"), function(x) (x[1]))
pgap_gene_name<-sapply(str_split(na.omit(pgap$name),"_"), function(x) (x[1]))
#### 
library(VennDiagram)
venn.diagram(
  x = list(prokka_gene_name, bakta_gene_name, eggnog_gene_name, pgap_gene_name),
  category.names = c("prokka" , "bakta " , "pgap", "eggnog"), filename = '#14_venn_diagramm.png',output=TRUE)




