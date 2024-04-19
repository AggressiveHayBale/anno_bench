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


#id <- args[1]
#
#prokka_faa<- args[2]
#prokka_gbk<- args[3]
#
#bakta_faa<- args[4]
#bakta_gbff<- args[5]
#
#eggnog_faa<- args[6]
#eggnog_gff<- args[7]
#
#pgap_faa<- args[8]
#pgap_gff<- args[9]
#
#fasta_type <- args[10]
id <- "GCA_000005825.2"

prokka_faa<- "/mnt/disk1/gits/anno_bench/results/GCA_000005825.2/prokka/original_GCA_000005825.2_prokka.faa"
prokka_gbk<- "/mnt/disk1/gits/anno_bench/results/GCA_000005825.2/prokka/original_GCA_000005825.2_prokka.gff"

bakta_faa<- "/mnt/disk1/gits/anno_bench/results/GCA_000005825.2/bakta/original_GCA_000005825.2_bakta.faa"
bakta_gbff<- "/mnt/disk1/gits/anno_bench/results/GCA_000005825.2/bakta/original_GCA_000005825.2_bakta.gff3"

eggnog_faa<- "/mnt/disk1/gits/anno_bench/results/GCA_000005825.2/eggnog/original_GCA_000005825.2_eggnog.emapper.genepred.fasta"
eggnog_gff<- "/mnt/disk1/gits/anno_bench/results/GCA_000005825.2/eggnog/original_GCA_000005825.2_eggnog.emapper.decorated.gff"

pgap_faa<- "/mnt/disk1/gits/anno_bench/results/GCA_000005825.2/pgap/original_GCA_000005825.2_pgap.faa"
pgap_gff<- "/mnt/disk1/gits/anno_bench/results/GCA_000005825.2/pgap/original_GCA_000005825.2_pgap.gff"

fasta_type <-  "original"
 
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
annotation$go <- str_count(annotation$X9,"GO")
annotation$cog <- str_count(annotation$X9,"COG[0-9]{1,4}")
annotation$prodigal <- str_count(annotation$X9,"Prodigal:")
annotation$uniprotkb <- str_count(annotation$X9,"UniProtKB:")
annotation$ec <- str_count(annotation$X9,"eC_number=")

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
fasta_data <- read_faa(bakta_faa)


#Cleanup
annotation <- invisible(na.omit(annotation))
annotation<- annotation[!(annotation$X3=="region" | annotation$X3=="gene"),]
annotation <- separate(annotation, X9, c(NA, "ID"), remove = FALSE, "ID=")
annotation$ID <- gsub(";.*","",annotation$ID)
annotation$product <- str_extract(annotation$X9,pattern = "Name=(.*?;)") 
annotation$product <- gsub("Name=", "", annotation$product)
annotation$product <- gsub(";", "", annotation$product)
annotation$product[annotation$product==""] <- NA
annotation$go <- str_count(annotation$X9,"GO:")

annotation$kegg <- str_count(annotation$X9,"KEGG:")
annotation$pfam <- str_count(annotation$X9,"PFAM:")
annotation$rfam <- str_count(annotation$X9,"RFAM:")
annotation$ec <- str_count(annotation$X9,"EC:")
annotation$cog <- str_count(annotation$X9,"COG:COG")
annotation$cog_group <- str_count(annotation$X9, "COG:(?!COG)")
annotation$so <- str_count(annotation$X9,"SO:")
annotation$uniparc <- str_count(annotation$X9,"UniParc:")
annotation$uniref100 <- str_count(annotation$X9,"UniRef:UniRef100_")
annotation$uniref90 <- str_count(annotation$X9,"UniRef:UniRef90_")
annotation$uniref50 <- str_count(annotation$X9,"UniRef:UniRef50_")



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
annotation$go <- str_count(annotation$X9,"GO:")
annotation$refseq <- str_count(annotation$X9,"inference=COORDINATES: similar to AA sequence:RefSeq:")
annotation$go_fun <- str_count(annotation$X9,"go_function=")
annotation$go_pro <- str_count(annotation$X9,"go_process=")
annotation$go_comp <- str_count(annotation$X9,"go_component=")

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

em_CAZy <- str_extract(annotation_merged$X9,pattern = "em_CAZy=(.*?;)")
em_CAZy <- gsub("em_CAZy=", "", em_CAZy)
em_CAZy <- gsub(";", "", em_CAZy)
ct_em_CAZy <- str_count(em_CAZy, pattern = ",") + 1 # Adding 1 because str_count counts separators
ct_em_CAZy[is.na(ct_em_CAZy) | ct_em_CAZy == 1] <- 0
annotation_merged$ct_em_CAZy <- ct_em_CAZy

em_COG_cat <- str_extract(annotation_merged$X9,pattern = "em_COG_cat=(.*?;)")
em_COG_cat <- gsub("em_COG_cat=", "", em_COG_cat)
em_COG_cat <- gsub(";", "", em_COG_cat)
em_COG_cat[em_COG_cat == "None"] <- NA

ct_em_COG_cat <- nchar(gsub("[^,]", "", em_COG_cat)) + 1
annotation_merged$ct_em_COG_cat <- ct_em_COG_cat

em_EC <- str_extract(annotation_merged$X9,pattern = "em_EC=(.*?;)")
em_EC <- gsub("em_EC=", "", em_EC)
em_EC <- gsub(";", "", em_EC)
ct_em_EC <- str_count(em_EC, pattern = ",") + 1 # Adding 1 because str_count counts separators
ct_em_EC[is.na(ct_em_EC) | ct_em_EC == 1] <- 0

annotation_merged$ct_em_EC <- ct_em_EC

em_KEGG_ko <- str_extract(annotation_merged$X9,pattern = "em_KEGG_ko=(.*?;)")
em_KEGG_ko <- gsub("em_KEGG_ko=", "", em_KEGG_ko)
em_KEGG_ko <- gsub(";", "", em_KEGG_ko)
ct_em_KEGG_ko <- str_count(em_KEGG_ko, pattern = ",") + 1 # Adding 1 because str_count counts separators
ct_em_KEGG_ko[is.na(ct_em_KEGG_ko) | ct_em_KEGG_ko == 1] <- 0

annotation_merged$ct_em_KEGG_ko <- ct_em_KEGG_ko

em_BRITE <- str_extract(annotation_merged$X9,pattern = "em_BRITE=(.*?;)")
em_BRITE <- gsub("em_BRITE=", "", em_BRITE)
em_BRITE <- gsub(";", "", em_BRITE)
ct_em_BRITE <- str_count(em_BRITE, pattern = ",") + 1 # Adding 1 because str_count counts separators
ct_em_BRITE[is.na(ct_em_BRITE) | ct_em_BRITE == 1] <- 0

annotation_merged$ct_em_BRITE <- ct_em_BRITE

em_BiGG_Reaction <- str_extract(annotation_merged$X9,pattern = "em_BiGG_Reaction=(.*?;)")
em_BiGG_Reaction <- gsub("em_BiGG_Reaction=", "", em_BiGG_Reaction)
em_BiGG_Reaction <- gsub(";", "", em_BiGG_Reaction)
ct_em_BiGG_Reaction <- str_count(em_BiGG_Reaction, pattern = ",") + 1 # Adding 1 because str_count counts separators
ct_em_BiGG_Reaction[is.na(ct_em_BiGG_Reaction) | ct_em_BiGG_Reaction == 1] <- 0

annotation_merged$ct_em_BiGG_Reaction <- ct_em_BiGG_Reaction

em_KEGG_rclass <- str_extract(annotation_merged$X9,pattern = "em_KEGG_rclass=(.*?;)")
em_KEGG_rclass <- gsub("em_KEGG_rclass=", "", em_KEGG_rclass)
em_KEGG_rclass <- gsub(";", "", em_KEGG_rclass)
ct_em_KEGG_rclass <- str_count(em_KEGG_rclass, pattern = ",") + 1 # Adding 1 because str_count counts separators
ct_em_KEGG_rclass[is.na(ct_em_KEGG_rclass) | ct_em_KEGG_rclass == 1] <- 0

annotation_merged$ct_em_KEGG_rclass <- ct_em_KEGG_rclass

em_KEGG_TC <- str_extract(annotation_merged$X9,pattern = "em_KEGG_TC=(.*?;)")
em_KEGG_TC <- gsub("em_KEGG_TC=", "", em_KEGG_TC)
em_KEGG_TC <- gsub(";", "", em_KEGG_TC)
ct_em_KEGG_TC <- str_count(em_KEGG_TC, pattern = ",") + 1 # Adding 1 because str_count counts separators
ct_em_KEGG_TC[is.na(ct_em_KEGG_TC) | ct_em_KEGG_TC == 1] <- 0

annotation_merged$ct_em_KEGG_TC <- ct_em_KEGG_TC

em_KEGG_Pathway <- str_extract(annotation_merged$X9,pattern = "em_KEGG_Pathway=(.*?;)")
em_KEGG_Pathway <- gsub("em_KEGG_Pathway=", "", em_KEGG_Pathway)
em_KEGG_Pathway <- gsub(";", "", em_KEGG_Pathway)
ct_em_KEGG_Pathway <- str_count(em_KEGG_Pathway, pattern = ",") + 1 # Adding 1 because str_count counts separators
ct_em_KEGG_Pathway[is.na(ct_em_KEGG_Pathway) | ct_em_KEGG_Pathway == 1] <- 0

annotation_merged$ct_em_KEGG_Pathway <- ct_em_KEGG_Pathway

em_KEGG_Reaction <- str_extract(annotation_merged$X9,pattern = "em_KEGG_Reaction=(.*?;)")
em_KEGG_Reaction <- gsub("em_KEGG_Reaction=", "", em_KEGG_Reaction)
em_KEGG_Reaction <- gsub(";", "", em_KEGG_Reaction)
ct_em_KEGG_Reaction <- str_count(em_KEGG_Reaction, pattern = ",") + 1 # Adding 1 because str_count counts separators
ct_em_KEGG_Reaction[is.na(ct_em_KEGG_Reaction) | ct_em_KEGG_Reaction == 1] <- 0

annotation_merged$ct_em_KEGG_Reaction <- ct_em_KEGG_Reaction

em_KEGG_Module <- str_extract(annotation_merged$X9,pattern = "em_KEGG_Module=(.*?;)")
em_KEGG_Module <- gsub("em_KEGG_Module=", "", em_KEGG_Module)
em_KEGG_Module <- gsub(";", "", em_KEGG_Module)
ct_em_KEGG_Module <- str_count(em_KEGG_Module, pattern = ",") + 1 # Adding 1 because str_count counts separators
ct_em_KEGG_Module[is.na(ct_em_KEGG_Module) | ct_em_KEGG_Module == 1] <- 0

annotation_merged$ct_em_KEGG_Module <- ct_em_KEGG_Module

em_PFAMs <- str_extract(annotation_merged$X9,pattern = "em_PFAMs=(.*?;)")
em_PFAMs <- gsub("em_PFAMs=", "", em_PFAMs)
em_PFAMs <- gsub(";", "", em_PFAMs)
ct_em_PFAMs <- str_count(em_PFAMs, pattern = ",") + 1 # Adding 1 because str_count counts separators
ct_em_PFAMs[is.na(ct_em_PFAMs) | ct_em_PFAMs == 1] <- 0

annotation_merged$ct_em_PFAMs <- ct_em_PFAMs

em_desc <- str_extract(annotation_merged$X9,pattern = "em_desc=(.*?;)")
em_desc <- gsub("em_desc=", "", em_desc)
em_desc <- gsub(";", "", em_desc)
em_desc[em_desc==""] <- NA
em_desc[em_desc=="None"] <- NA

ct_Bacteria_prot_un_fun <- str_count(em_desc, pattern = "Bacterial protein of unknown function")
ct_Bacteria_prot_un_fun[is.na(ct_Bacteria_prot_un_fun)==TRUE] <- 0
annotation_merged$ct_Bacteria_prot_un_fun <- ct_Bacteria_prot_un_fun

ct_Prot_un_fun <- str_count(em_desc, pattern = "Protein of unknown function")
ct_Prot_un_fun[is.na(ct_Prot_un_fun)==TRUE] <- 0
annotation_merged$ct_Prot_un_fun <- ct_Prot_un_fun

ct_Domain_un_fun <- str_count(em_desc, pattern = "Domain of unknown function")
ct_Domain_un_fun[is.na(ct_Domain_un_fun)==TRUE] <- 0
annotation_merged$ct_Domain_un_fun <- ct_Domain_un_fun

#Functions
annotation_merged$go <- str_count(annotation_merged$X9,"GO:")

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
eggnog_hyps <- sum(is.na(eggnog$product))


#Percent of hyps 
prokka_perc_hyps <- (prokka_hyps/prokka_gene_ct)*100
bakta_perc_hyps <- (bakta_hyps/bakta_gene_ct)*100
pgap_perc_hyps <- (pgap_hyps/pgap_gene_ct)*100
eggnog_perc_hyps<- (eggnog_hyps/eggnog_gene_ct)*100

#Named
named_prokka_seq_lenght<- mean(nchar(na.omit(prokka[prokka$product!="hypothetical protein",]$sequence)))
named_bakta_seq_lenght<-mean(nchar(na.omit(bakta[bakta$product!="hypothetical protein",]$sequence)))
named_pgap_seq_lenght<-mean(nchar(na.omit(pgap[pgap$product!="hypothetical protein",]$sequence)))
named_eggnog_seq_lenght<-mean(nchar(eggnog$sequence[!is.na(eggnog$product)]))

#Hyp gene lenght
hyp_prokka_seq_lenght<- mean(nchar(na.omit(prokka[prokka$product=="hypothetical protein",]$sequence)))
hyp_bakta_seq_lenght<-mean(nchar(na.omit(bakta[bakta$product=="hypothetical protein",]$sequence)))
hyp_pgap_seq_lenght<-mean(nchar(na.omit(pgap[pgap$product=="hypothetical protein",]$sequence)))
hyp_eggnog_seq_lenght<-mean(nchar(eggnog$sequence[is.na(eggnog$product)]))

#Functional analysis

## GO annotated proteins
prokka_go_sum<- sum(prokka$go!=0)
bakta_go_sum<- sum(bakta$go!=0)
eggnog_go_sum<- sum(eggnog$go!=0)
pgap_go_sum<- sum(pgap$go!=0)

## median GO terms 
prokka_med_go<- prokka %>% filter(go!=0) %>% select(go) %>% summarize(median(go)) %>%  as.integer()
bakta_med_go<-bakta %>% filter(go!=0) %>% select(go) %>% summarize(median(go)) %>%  as.integer()
eggnog_med_go<-eggnog %>% filter(go!=0) %>% select(go) %>% summarize(median(go)) %>%  as.integer()
pgap_med_go<-pgap %>% filter(go!=0) %>% select(go) %>% summarize(median(go)) %>%  as.integer()

## GO annotated proteins
prokka_go_perc<- (sum(prokka$go!=0)/nrow(prokka))*100
bakta_go_perc<- (sum(bakta$go!=0)/nrow(bakta))*100
eggnog_go_perc<- (sum(eggnog$go!=0)/nrow(eggnog))*100
pgap_go_perc<- (sum(pgap$go!=0)/nrow(pgap))*100

##GO in non hyps 
prokka_go_nonhyp_perc<-(sum(prokka$go!=0)/(prokka_gene_ct-prokka_hyps))*100
bakta_go_nonhyp_perc<- (sum(bakta$go!=0)/(bakta_gene_ct- bakta_hyps ))*100
eggnog_go_nonhyp_perc<- (sum(eggnog$go!=0)/(eggnog_gene_ct- eggnog_hyps))*100
pgap_go_nonhyp_perc<- (sum(pgap$go!=0)/(pgap_gene_ct-pgap_hyps))*100



#Prokka
### COG annotated proteins
prokka_cog_sum<- sum(prokka$cog!=0)
prokka_prodigal_sum<- sum(prokka$prodigal!=0)
prokka_uniprotkb_sum<- sum(prokka$uniprotkb!=0)
prokka_ec_sum<- sum(prokka$ec!=0)

## Percentage of annotated proteins
prokka_cog_perc<- (sum(prokka$cog!=0)/nrow(prokka))*100
prokka_prodigal_perc<- (sum(prokka$prodigal!=0)/nrow(prokka))*100
prokka_uniprotkb_perc<- (sum(prokka$uniprotkb!=0)/nrow(prokka))*100
prokka_ec_perc<- (sum(prokka$ec!=0)/nrow(prokka))*100
## Percentage of annotated proteins excluding hypotheticals
prokka_cog_nonhyp_perc<-(sum(prokka$cog!=0)/(prokka_gene_ct-prokka_hyps))*100
prokka_prodigal_nonhyp_perc<- (sum(prokka$prodigal!=0)/(prokka_gene_ct- prokka_hyps ))*100
prokka_uniprotkb_nonhyp_perc<- (sum(prokka$uniprotkb!=0)/(prokka_gene_ct- prokka_hyps))*100
prokka_ec_nonhyp_perc<- (sum(prokka$ec!=0)/(prokka_gene_ct-prokka_hyps))*100
## median COG terms
prokka_med_cog<- prokka %>% filter(cog!=0) %>% select(cog) %>% summarize(median(cog)) %>%  as.integer()
prokka_med_prodigal<-prokka %>% filter(prodigal!=0) %>% select(prodigal) %>% summarize(median(prodigal)) %>%  as.integer()
prokka_med_uniprotkb<-prokka %>% filter(uniprotkb!=0) %>% select(uniprotkb) %>% summarize(median(uniprotkb)) %>%  as.integer()
prokka_med_ec<-prokka %>% filter(ec!=0) %>% select(ec) %>% summarize(median(ec)) %>%  as.integer()

####Bakta
## Sum of annotated proteins for each category
bakta_kegg_sum<- sum(bakta$kegg!=0)
bakta_pfam_sum<- sum(bakta$pfam!=0)
bakta_rfam_sum<- sum(bakta$rfam!=0)
bakta_ec_sum<- sum(bakta$ec!=0)
bakta_cog_sum<- sum(bakta$cog!=0)
bakta_cog_group_sum<- sum(bakta$cog_group!=0)
bakta_so_sum<- sum(bakta$so!=0)
bakta_uniparc_sum<- sum(bakta$uniparc!=0)
bakta_uniref100_sum<- sum(bakta$uniref100!=0)
bakta_uniref90_sum<- sum(bakta$uniref90!=0)
bakta_uniref50_sum<- sum(bakta$uniref50!=0)

## Percentage of annotated proteins
bakta_kegg_perc<- (sum(bakta$kegg!=0)/nrow(bakta))*100
bakta_pfam_perc<- (sum(bakta$pfam!=0)/nrow(bakta))*100
bakta_rfam_perc<- (sum(bakta$rfam!=0)/nrow(bakta))*100
bakta_ec_perc<- (sum(bakta$ec!=0)/nrow(bakta))*100
bakta_cog_perc<- (sum(bakta$cog!=0)/nrow(bakta))*100
bakta_cog_group_perc<- (sum(bakta$cog_group!=0)/nrow(bakta))*100
bakta_so_perc<- (sum(bakta$so!=0)/nrow(bakta))*100
bakta_uniparc_perc<- (sum(bakta$uniparc!=0)/nrow(bakta))*100
bakta_uniref100_perc<- (sum(bakta$uniref100!=0)/nrow(bakta))*100
bakta_uniref90_perc<- (sum(bakta$uniref90!=0)/nrow(bakta))*100
bakta_uniref50_perc<- (sum(bakta$uniref50!=0)/nrow(bakta))*100

## Med
bakta_med_kegg <- bakta %>% filter(kegg!=0) %>% select(kegg) %>% summarize(median(kegg)) %>% as.integer()
bakta_med_pfam <- bakta %>% filter(pfam!=0) %>% select(pfam) %>% summarize(median(pfam)) %>% as.integer()
bakta_med_rfam <- bakta %>% filter(rfam!=0) %>% select(rfam) %>% summarize(median(rfam)) %>% as.integer()
bakta_med_ec <- bakta %>% filter(ec!=0) %>% select(ec) %>% summarize(median(ec)) %>% as.integer()
bakta_med_cog <- bakta %>% filter(cog!=0) %>% select(cog) %>% summarize(median(cog)) %>% as.integer()
bakta_med_cog_group <- bakta %>% filter(cog_group!=0) %>% select(cog_group) %>% summarize(median(cog_group)) %>% as.integer()
bakta_med_so <- bakta %>% filter(so!=0) %>% select(so) %>% summarize(median(so)) %>% as.integer()
bakta_med_uniparc <- bakta %>% filter(uniparc!=0) %>% select(uniparc) %>% summarize(median(uniparc)) %>% as.integer()
bakta_med_uniref100 <- bakta %>% filter(uniref100!=0) %>% select(uniref100) %>% summarize(median(uniref100)) %>% as.integer()
bakta_med_uniref90 <- bakta %>% filter(uniref90!=0) %>% select(uniref90) %>% summarize(median(uniref90)) %>% as.integer()
bakta_med_uniref50 <- bakta %>% filter(uniref50!=0) %>% select(uniref50) %>% summarize(median(uniref50)) %>% as.integer()

## Percentage of annotated proteins excluding hypotheticals
bakta_kegg_nonhyp_perc<-(sum(bakta$kegg!=0)/(bakta_gene_ct-bakta_hyps))*100
bakta_pfam_nonhyp_perc<- (sum(bakta$pfam!=0)/(bakta_gene_ct- bakta_hyps ))*100
bakta_rfam_nonhyp_perc<- (sum(bakta$rfam!=0)/(bakta_gene_ct- bakta_hyps))*100
bakta_ec_nonhyp_perc<- (sum(bakta$ec!=0)/(bakta_gene_ct-bakta_hyps))*100
bakta_cog_nonhyp_perc<- (sum(bakta$cog!=0)/(bakta_gene_ct-bakta_hyps))*100
bakta_cog_group_nonhyp_perc<- (sum(bakta$cog_group!=0)/(bakta_gene_ct-bakta_hyps))*100
bakta_so_nonhyp_perc<- (sum(bakta$so!=0)/(bakta_gene_ct-bakta_hyps))*100
bakta_uniparc_nonhyp_perc<- (sum(bakta$uniparc!=0)/(bakta_gene_ct-bakta_hyps))*100
bakta_uniref100_nonhyp_perc<- (sum(bakta$uniref100!=0)/(bakta_gene_ct-bakta_hyps))*100
bakta_uniref90_nonhyp_perc<- (sum(bakta$uniref90!=0)/(bakta_gene_ct-bakta_hyps))*100
bakta_uniref50_nonhyp_perc<- (sum(bakta$uniref50!=0)/(bakta_gene_ct-bakta_hyps))*100


### EGGNOG
## Sum of annotated proteins for each category
eggnog_ct_em_CAZy_sum <- sum(eggnog$ct_em_CAZy!=0)
eggnog_ct_em_COG_cat_sum <- sum(eggnog$ct_em_COG_cat!=0)
eggnog_ct_em_EC_sum <- sum(eggnog$ct_em_EC!=0)
eggnog_ct_em_KEGG_ko_sum <- sum(eggnog$ct_em_KEGG_ko!=0)
eggnog_ct_em_BRITE_sum <- sum(eggnog$ct_em_BRITE!=0)
eggnog_ct_em_BiGG_Reaction_sum <- sum(eggnog$ct_em_BiGG_Reaction!=0)
eggnog_ct_em_KEGG_rclass_sum <- sum(eggnog$ct_em_KEGG_rclass!=0)
eggnog_ct_em_KEGG_TC_sum <- sum(eggnog$ct_em_KEGG_TC!=0)
eggnog_ct_em_KEGG_Pathway_sum <- sum(eggnog$ct_em_KEGG_Pathway!=0)
eggnog_ct_em_KEGG_Reaction_sum <- sum(eggnog$ct_em_KEGG_Reaction!=0)
eggnog_ct_em_KEGG_Module_sum <- sum(eggnog$ct_em_KEGG_Module!=0)
eggnog_ct_em_PFAMs_sum <- sum(eggnog$ct_em_PFAMs!=0)

## Percentage of proteins annotated in each category
eggnog_ct_em_CAZy_perc <- (sum(eggnog$ct_em_CAZy!=0)/nrow(eggnog))*100
eggnog_ct_em_COG_cat_perc <- (sum(eggnog$ct_em_COG_cat!=0)/nrow(eggnog))*100
eggnog_ct_em_EC_perc <- (sum(eggnog$ct_em_EC!=0)/nrow(eggnog))*100
eggnog_ct_em_KEGG_ko_perc <- (sum(eggnog$ct_em_KEGG_ko!=0)/nrow(eggnog))*100
eggnog_ct_em_BRITE_perc <- (sum(eggnog$ct_em_BRITE!=0)/nrow(eggnog))*100
eggnog_ct_em_BiGG_Reaction_perc <- (sum(eggnog$ct_em_BiGG_Reaction!=0)/nrow(eggnog))*100
eggnog_ct_em_KEGG_rclass_perc <- (sum(eggnog$ct_em_KEGG_rclass!=0)/nrow(eggnog))*100
eggnog_ct_em_KEGG_TC_perc <- (sum(eggnog$ct_em_KEGG_TC!=0)/nrow(eggnog))*100
eggnog_ct_em_KEGG_Pathway_perc <- (sum(eggnog$ct_em_KEGG_Pathway!=0)/nrow(eggnog))*100
eggnog_ct_em_KEGG_Reaction_perc <- (sum(eggnog$ct_em_KEGG_Reaction!=0)/nrow(eggnog))*100
eggnog_ct_em_KEGG_Module_perc <- (sum(eggnog$ct_em_KEGG_Module!=0)/nrow(eggnog))*100
eggnog_ct_em_PFAMs_perc <- (sum(eggnog$ct_em_PFAMs!=0)/nrow(eggnog))*100
## Median terms for each category
eggnog_med_ct_em_CAZy <- eggnog %>% filter(ct_em_CAZy!=0) %>% select(ct_em_CAZy) %>% summarize(median(ct_em_CAZy)) %>% as.integer()
eggnog_med_ct_em_COG_cat <- eggnog %>% filter(ct_em_COG_cat!=0) %>% select(ct_em_COG_cat) %>% summarize(median(ct_em_COG_cat)) %>% as.integer()
eggnog_med_ct_em_EC <- eggnog %>% filter(ct_em_EC!=0) %>% select(ct_em_EC) %>% summarize(median(ct_em_EC)) %>% as.integer()
eggnog_med_ct_em_KEGG_ko <- eggnog %>% filter(ct_em_KEGG_ko!=0) %>% select(ct_em_KEGG_ko) %>% summarize(median(ct_em_KEGG_ko)) %>% as.integer()
eggnog_med_ct_em_BRITE <- eggnog %>% filter(ct_em_BRITE!=0) %>% select(ct_em_BRITE) %>% summarize(median(ct_em_BRITE)) %>% as.integer()
eggnog_med_ct_em_BiGG_Reaction <- eggnog %>% filter(ct_em_BiGG_Reaction!=0) %>% select(ct_em_BiGG_Reaction) %>% summarize(median(ct_em_BiGG_Reaction)) %>% as.integer()
eggnog_med_ct_em_KEGG_rclass <- eggnog %>% filter(ct_em_KEGG_rclass!=0) %>% select(ct_em_KEGG_rclass) %>% summarize(median(ct_em_KEGG_rclass)) %>% as.integer()
eggnog_med_ct_em_KEGG_TC <- eggnog %>% filter(ct_em_KEGG_TC!=0) %>% select(ct_em_KEGG_TC) %>% summarize(median(ct_em_KEGG_TC)) %>% as.integer()
eggnog_med_ct_em_KEGG_Pathway <- eggnog %>% filter(ct_em_KEGG_Pathway!=0) %>% select(ct_em_KEGG_Pathway) %>% summarize(median(ct_em_KEGG_Pathway)) %>% as.integer()
eggnog_med_ct_em_KEGG_Reaction <- eggnog %>% filter(ct_em_KEGG_Reaction!=0) %>% select(ct_em_KEGG_Reaction) %>% summarize(median(ct_em_KEGG_Reaction)) %>% as.integer()
eggnog_med_ct_em_KEGG_Module <- eggnog %>% filter(ct_em_KEGG_Module!=0) %>% select(ct_em_KEGG_Module) %>% summarize(median(ct_em_KEGG_Module)) %>% as.integer()
eggnog_med_ct_em_PFAMs <- eggnog %>% filter(ct_em_PFAMs!=0) %>% select(ct_em_PFAMs) %>% summarize(median(ct_em_PFAMs)) %>% as.integer()

## Percentage of proteins annotated in each category, excluding hypotheticals
eggnog_ct_em_CAZy_nonhyp_perc <- (sum(eggnog$ct_em_CAZy!=0)/(eggnog_gene_ct-eggnog_hyps))*100
eggnog_ct_em_COG_cat_nonhyp_perc <- (sum(eggnog$ct_em_COG_cat!=0)/(eggnog_gene_ct-eggnog_hyps))*100
eggnog_ct_em_EC_nonhyp_perc <- (sum(eggnog$ct_em_EC!=0)/(eggnog_gene_ct-eggnog_hyps))*100
eggnog_ct_em_KEGG_ko_nonhyp_perc <- (sum(eggnog$ct_em_KEGG_ko!=0)/(eggnog_gene_ct-eggnog_hyps))*100
eggnog_ct_em_BRITE_nonhyp_perc <- (sum(eggnog$ct_em_BRITE!=0)/(eggnog_gene_ct-eggnog_hyps))*100
eggnog_ct_em_BiGG_Reaction_nonhyp_perc <- (sum(eggnog$ct_em_BiGG_Reaction!=0)/(eggnog_gene_ct-eggnog_hyps))*100
eggnog_ct_em_KEGG_rclass_nonhyp_perc <- (sum(eggnog$ct_em_KEGG_rclass!=0)/(eggnog_gene_ct-eggnog_hyps))*100
eggnog_ct_em_KEGG_TC_nonhyp_perc <- (sum(eggnog$ct_em_KEGG_TC!=0)/(eggnog_gene_ct-eggnog_hyps))*100
eggnog_ct_em_KEGG_Pathway_nonhyp_perc <- (sum(eggnog$ct_em_KEGG_Pathway!=0)/(eggnog_gene_ct-eggnog_hyps))*100
eggnog_ct_em_KEGG_Reaction_nonhyp_perc <- (sum(eggnog$ct_em_KEGG_Reaction!=0)/(eggnog_gene_ct-eggnog_hyps))*100
eggnog_ct_em_KEGG_Module_nonhyp_perc <- (sum(eggnog$ct_em_KEGG_Module!=0)/(eggnog_gene_ct-eggnog_hyps))*100
eggnog_ct_em_PFAMs_nonhyp_perc <- (sum(eggnog$ct_em_PFAMs!=0)/(eggnog_gene_ct-eggnog_hyps))*100


####PGAP
## Sum of annotated proteins for each category
pgap_refseq_sum<- sum(pgap$refseq!=0)
pgap_go_fun_sum<- sum(pgap$go_fun!=0)
pgap_go_pro_sum<- sum(pgap$go_pro!=0)
pgap_go_comp_sum<- sum(pgap$go_comp!=0)

## Percentage of annotated proteins
pgap_refseq_perc<- (sum(pgap$refseq!=0)/nrow(pgap))*100
pgap_go_fun_perc<- (sum(pgap$go_fun!=0)/nrow(pgap))*100
pgap_go_pro_perc<- (sum(pgap$go_pro!=0)/nrow(pgap))*100
pgap_go_comp_perc<- (sum(pgap$go_comp!=0)/nrow(pgap))*100

## Percentage of annotated proteins excluding hypotheticals
pgap_refseq_nonhyp_perc<-(sum(pgap$refseq!=0)/(pgap_gene_ct-pgap_hyps))*100
pgap_go_fun_nonhyp_perc<- (sum(pgap$go_fun!=0)/(pgap_gene_ct- pgap_hyps ))*100
pgap_go_pro_nonhyp_perc<- (sum(pgap$go_pro!=0)/(pgap_gene_ct- pgap_hyps))*100
pgap_go_comp_nonhyp_perc<- (sum(pgap$go_comp!=0)/(pgap_gene_ct-pgap_hyps))*100

## Median terms for each category
pgap_med_refseq<- pgap %>% filter(refseq!=0) %>% select(refseq) %>% summarize(median(refseq)) %>%  as.integer()
pgap_med_go_fun<-pgap %>% filter(go_fun!=0) %>% select(go_fun) %>% summarize(median(go_fun)) %>%  as.integer()
pgap_med_go_pro<-pgap %>% filter(go_pro!=0) %>% select(go_pro) %>% summarize(median(go_pro)) %>%  as.integer()
pgap_med_go_comp<-pgap %>% filter(go_comp!=0) %>% select(go_comp) %>% summarize(median(go_comp)) %>%  as.integer()



#Create summary

#summary_df <- data.frame(
#  Tool = c(rep("Prokka", 4), rep("Bakta", 11), rep("EggNOG", 12), rep("PGAP", 4)),
#  Category = c("COG", "Prodigal", "UniProtKB", "EC",
#               "KEGG", "PFAM", "RFAM", "EC", "COG", "COG Group", "SO", "UniParc", "UniRef100", "UniRef90", "UniRef50",
#               "CAZy", "COG Categories", "EC", "KEGG ko", "BRITE", "BiGG Reaction", "KEGG rclass", "KEGG TC", "KEGG Pathway", "KEGG Reaction", "KEGG Module", "PFAMs",
#               "RefSeq", "GO Fun", "GO Pro", "GO Comp"),
#  Sum = c(prokka_cog_sum, prokka_prodigal_sum, prokka_uniprotkb_sum, prokka_ec_sum,
#          bakta_kegg_sum, bakta_pfam_sum, bakta_rfam_sum, bakta_ec_sum, bakta_cog_sum, bakta_cog_group_sum, bakta_so_sum, bakta_uniparc_sum, bakta_uniref100_sum, bakta_uniref90_sum, bakta_uniref50_sum,
#          eggnog_ct_em_CAZy_sum, eggnog_ct_em_COG_cat_sum, eggnog_ct_em_EC_sum, eggnog_ct_em_KEGG_ko_sum, eggnog_ct_em_BRITE_sum, eggnog_ct_em_BiGG_Reaction_sum, eggnog_ct_em_KEGG_rclass_sum, eggnog_ct_em_KEGG_TC_sum, eggnog_ct_em_KEGG_Pathway_sum, eggnog_ct_em_KEGG_Reaction_sum, eggnog_ct_em_KEGG_Module_sum, eggnog_ct_em_PFAMs_sum,
#          pgap_refseq_sum, pgap_go_fun_sum, pgap_go_pro_sum, pgap_go_comp_sum),
#  Percentage = c(prokka_cog_perc, prokka_prodigal_perc, prokka_uniprotkb_perc, prokka_ec_perc,
#                 bakta_kegg_perc, bakta_pfam_perc, bakta_rfam_perc, bakta_ec_perc, bakta_cog_perc, bakta_cog_group_perc, bakta_so_perc, bakta_uniparc_perc, bakta_uniref100_perc, bakta_uniref90_perc, bakta_uniref50_perc,
#                 eggnog_ct_em_CAZy_perc, eggnog_ct_em_COG_cat_perc, eggnog_ct_em_EC_perc, eggnog_ct_em_KEGG_ko_perc, eggnog_ct_em_BRITE_perc, eggnog_ct_em_BiGG_Reaction_perc, eggnog_ct_em_KEGG_rclass_perc, eggnog_ct_em_KEGG_TC_perc, eggnog_ct_em_KEGG_Pathway_perc, eggnog_ct_em_KEGG_Reaction_perc, eggnog_ct_em_KEGG_Module_perc, eggnog_ct_em_PFAMs_perc,
#                 pgap_refseq_perc, pgap_go_fun_perc, pgap_go_pro_perc, pgap_go_comp_perc),
#  PercentageExclHypotheticals = c(prokka_cog_nonhyp_perc, prokka_prodigal_nonhyp_perc, prokka_uniprotkb_nonhyp_perc, prokka_ec_nonhyp_perc,
#                 bakta_kegg_nonhyp_perc, bakta_pfam_nonhyp_perc, bakta_rfam_nonhyp_perc, bakta_ec_nonhyp_perc, bakta_cog_nonhyp_perc, bakta_cog_group_nonhyp_perc, bakta_so_nonhyp_perc, bakta_uniparc_nonhyp_perc, bakta_uniref100_nonhyp_perc, bakta_uniref90_nonhyp_perc, bakta_uniref50_nonhyp_perc,
#                 eggnog_ct_em_CAZy_nonhyp_perc, eggnog_ct_em_COG_cat_nonhyp_perc, eggnog_ct_em_EC_nonhyp_perc, eggnog_ct_em_KEGG_ko_nonhyp_perc, eggnog_ct_em_BRITE_nonhyp_perc, eggnog_ct_em_BiGG_Reaction_nonhyp_perc, eggnog_ct_em_KEGG_rclass_nonhyp_perc, eggnog_ct_em_KEGG_TC_nonhyp_perc, eggnog_ct_em_KEGG_Pathway_nonhyp_perc, eggnog_ct_em_KEGG_Reaction_nonhyp_perc, eggnog_ct_em_KEGG_Module_nonhyp_perc, eggnog_ct_em_PFAMs_nonhyp_perc,
#                 pgap_refseq_nonhyp_perc, pgap_go_fun_nonhyp_perc, pgap_go_pro_nonhyp_perc, pgap_go_comp_nonhyp_perc),
#  MedianTerms = c(prokka_med_cog, prokka_med_prodigal, prokka_med_uniprotkb, prokka_med_ec,
#                  bakta_med_kegg, bakta_med_pfam, bakta_med_rfam, bakta_med_ec, bakta_med_cog, bakta_med_cog_group, bakta_med_so, bakta_med_uniparc, bakta_med_uniref100, bakta_med_uniref90, bakta_med_uniref50,
#                  eggnog_med_ct_em_CAZy, eggnog_med_ct_em_COG_cat, eggnog_med_ct_em_EC, eggnog_med_ct_em_KEGG_ko, eggnog_med_ct_em_BRITE, eggnog_med_ct_em_BiGG_Reaction, eggnog_med_ct_em_KEGG_rclass, eggnog_med_ct_em_KEGG_TC, eggnog_med_ct_em_KEGG_Pathway, eggnog_med_ct_em_KEGG_Reaction, eggnog_med_ct_em_KEGG_Module, eggnog_med_ct_em_PFAMs,
#                  pgap_med_refseq, pgap_med_go_fun, pgap_med_go_pro, pgap_med_go_comp)
#)
summary_df <- data.frame(
  Prokka_COG_Sum = prokka_cog_sum,
  Prokka_Prodigal_Sum = prokka_prodigal_sum,
  Prokka_UniprotKB_Sum = prokka_uniprotkb_sum,
  Prokka_EC_Sum = prokka_ec_sum,
  Bakta_KEGG_Sum = bakta_kegg_sum,
  Bakta_PFAM_Sum = bakta_pfam_sum,
  Bakta_RFAM_Sum = bakta_rfam_sum,
  Bakta_EC_Sum = bakta_ec_sum,
  Bakta_COG_Sum = bakta_cog_sum,
  Bakta_COG_Group_Sum = bakta_cog_group_sum,
  Bakta_SO_Sum = bakta_so_sum,
  Bakta_UniParc_Sum = bakta_uniparc_sum,
  Bakta_UniRef100_Sum = bakta_uniref100_sum,
  Bakta_UniRef90_Sum = bakta_uniref90_sum,
  Bakta_UniRef50_Sum = bakta_uniref50_sum,
  EggNOG_gene_ct = eggnog_gene_ct,
  EggNOG_hyps = eggnog_hyps,
  EggNOG_perc_hyps = eggnog_perc_hyps,
  named_EggNOG_seq_lenght = named_eggnog_seq_lenght,
  hyp_EggNOG_seq_lenght = hyp_eggnog_seq_lenght,
  EggNOG_CAZy_Sum = eggnog_ct_em_CAZy_sum,
  EggNOG_COG_Cat_Sum = eggnog_ct_em_COG_cat_sum,
  EggNOG_EC_Sum = eggnog_ct_em_EC_sum,
  EggNOG_KEGG_ko_Sum = eggnog_ct_em_KEGG_ko_sum,
  EggNOG_BRITE_Sum = eggnog_ct_em_BRITE_sum,
  EggNOG_BiGG_Reaction_Sum = eggnog_ct_em_BiGG_Reaction_sum,
  EggNOG_KEGG_rclass_Sum = eggnog_ct_em_KEGG_rclass_sum,
  EggNOG_KEGG_TC_Sum = eggnog_ct_em_KEGG_TC_sum,
  EggNOG_KEGG_Pathway_Sum = eggnog_ct_em_KEGG_Pathway_sum,
  EggNOG_KEGG_Reaction_Sum = eggnog_ct_em_KEGG_Reaction_sum,
  EggNOG_KEGG_Module_Sum = eggnog_ct_em_KEGG_Module_sum,
  EggNOG_PFAMs_Sum = eggnog_ct_em_PFAMs_sum,
  PGAP_RefSeq_Sum = pgap_refseq_sum,
  PGAP_GO_Fun_Sum = pgap_go_fun_sum,
  PGAP_GO_Pro_Sum = pgap_go_pro_sum,
  PGAP_GO_Comp_Sum = pgap_go_comp_sum,
  Prokka_COG_Perc = prokka_cog_perc,
  Prokka_Prodigal_Perc = prokka_prodigal_perc,
  Prokka_UniprotKB_Perc = prokka_uniprotkb_perc,
  Prokka_EC_Perc = prokka_ec_perc,
  Bakta_KEGG_Perc = bakta_kegg_perc,
  Bakta_PFAM_Perc = bakta_pfam_perc,
  Bakta_RFAM_Perc = bakta_rfam_perc,
  Bakta_EC2_Perc = bakta_ec_perc,
  Bakta_COG_Perc = bakta_cog_perc,
  Bakta_COG_Group_Perc = bakta_cog_group_perc,
  Bakta_SO_Perc = bakta_so_perc,
  Bakta_UniParc_Perc = bakta_uniparc_perc,
  Bakta_UniRef100_Perc = bakta_uniref100_perc,
  Bakta_UniRef90_Perc = bakta_uniref90_perc,
  Bakta_UniRef50_Perc = bakta_uniref50_perc,
  EggNOG_CAZy_Perc = eggnog_ct_em_CAZy_perc,
  EggNOG_COG_Cat_Perc = eggnog_ct_em_COG_cat_perc,
  EggNOG_EC_Perc = eggnog_ct_em_EC_perc,
  EggNOG_KEGG_ko_Perc = eggnog_ct_em_KEGG_ko_perc,
  EggNOG_BRITE_Perc = eggnog_ct_em_BRITE_perc,
  EggNOG_BiGG_Reaction_Perc = eggnog_ct_em_BiGG_Reaction_perc,
  EggNOG_KEGG_rclass_Perc = eggnog_ct_em_KEGG_rclass_perc,
  EggNOG_KEGG_TC_Perc = eggnog_ct_em_KEGG_TC_perc,
  EggNOG_KEGG_Pathway_Perc = eggnog_ct_em_KEGG_Pathway_perc,
  EggNOG_KEGG_Reaction_Perc = eggnog_ct_em_KEGG_Reaction_perc,
  EggNOG_KEGG_Module_Perc = eggnog_ct_em_KEGG_Module_perc,
  EggNOG_PFAMs_Perc = eggnog_ct_em_PFAMs_perc,
  PGAP_RefSeq_Perc = pgap_refseq_perc,
  PGAP_GO_Fun_Perc = pgap_go_fun_perc,
  PGAP_GO_Pro_Perc = pgap_go_pro_perc,
  PGAP_GO_Comp_Perc = pgap_go_comp_perc,
  Prokka_COG_NonHyp_Perc = prokka_cog_nonhyp_perc,
  Prokka_Prodigal_NonHyp_Perc = prokka_prodigal_nonhyp_perc,
  Prokka_UniprotKB_NonHyp_Perc = prokka_uniprotkb_nonhyp_perc,
  Prokka_EC_NonHyp_Perc = prokka_ec_nonhyp_perc,
  Bakta_KEGG_NonHyp_Perc = bakta_kegg_nonhyp_perc,
  Bakta_PFAM_NonHyp_Perc = bakta_pfam_nonhyp_perc,
  Bakta_RFAM_NonHyp_Perc = bakta_rfam_nonhyp_perc,
  Bakta_EC2_NonHyp_Perc = bakta_ec_nonhyp_perc,
  Bakta_COG_NonHyp_Perc = bakta_cog_nonhyp_perc,
  Bakta_COG_Group_NonHyp_Perc = bakta_cog_group_nonhyp_perc,
  Bakta_SO_NonHyp_Perc = bakta_so_nonhyp_perc,
  Bakta_UniParc_NonHyp_Perc = bakta_uniparc_nonhyp_perc,
  Bakta_UniRef100_NonHyp_Perc = bakta_uniref100_nonhyp_perc,
  Bakta_UniRef90_NonHyp_Perc = bakta_uniref90_nonhyp_perc,
  Bakta_UniRef50_NonHyp_Perc = bakta_uniref50_nonhyp_perc,
  EggNOG_CAZy_NonHyp_Perc = eggnog_ct_em_CAZy_nonhyp_perc,
  EggNOG_COG_Cat_NonHyp_Perc = eggnog_ct_em_COG_cat_nonhyp_perc,
  EggNOG_EC_NonHyp_Perc = eggnog_ct_em_EC_nonhyp_perc,
  EggNOG_KEGG_ko_NonHyp_Perc = eggnog_ct_em_KEGG_ko_nonhyp_perc,
  EggNOG_BRITE_NonHyp_Perc = eggnog_ct_em_BRITE_nonhyp_perc,
  EggNOG_BiGG_Reaction_NonHyp_Perc = eggnog_ct_em_BiGG_Reaction_nonhyp_perc,
  EggNOG_KEGG_rclass_NonHyp_Perc = eggnog_ct_em_KEGG_rclass_nonhyp_perc,
  EggNOG_KEGG_TC_NonHyp_Perc = eggnog_ct_em_KEGG_TC_nonhyp_perc,
  EggNOG_KEGG_Pathway_NonHyp_Perc = eggnog_ct_em_KEGG_Pathway_nonhyp_perc,
  EggNOG_KEGG_Reaction_NonHyp_Perc = eggnog_ct_em_KEGG_Reaction_nonhyp_perc,
  EggNOG_KEGG_Module_NonHyp_Perc = eggnog_ct_em_KEGG_Module_nonhyp_perc,
  EggNOG_PFAMs_NonHyp_Perc = eggnog_ct_em_PFAMs_nonhyp_perc,
  PGAP_RefSeq_NonHyp_Perc = pgap_refseq_nonhyp_perc,
  PGAP_GO_Fun_NonHyp_Perc = pgap_go_fun_nonhyp_perc,
  PGAP_GO_Pro_NonHyp_Perc = pgap_go_pro_nonhyp_perc,
  PGAP_GO_Comp_NonHyp_Perc = pgap_go_comp_nonhyp_perc,
  Prokka_Med_COG = prokka_med_cog,
  Prokka_Med_Prodigal = prokka_med_prodigal,
  Prokka_Med_UniprotKB = prokka_med_uniprotkb,
  Prokka_Med_EC = prokka_med_ec,
  Bakta_Med_KEGG = bakta_med_kegg,
  Bakta_Med_PFAM = bakta_med_pfam,
  Bakta_Med_RFAM = bakta_med_rfam,
  Bakta_Med_EC2 = bakta_med_ec,
  Bakta_Med_COG = bakta_med_cog,
  Bakta_Med_COG_Group = bakta_med_cog_group,
  Bakta_Med_SO = bakta_med_so,
  Bakta_Med_UniParc = bakta_med_uniparc,
  Bakta_Med_UniRef100 = bakta_med_uniref100,
  Bakta_Med_UniRef90 = bakta_med_uniref90,
  Bakta_Med_UniRef50 = bakta_med_uniref50,
  EggNOG_Bacteria_prot_un_fun=eggnog$ct_Bacteria_prot_un_fun,
  EggNOG_Prot_un_fun=eggnog$ct_Prot_un_fun,
  EggNOG_Domain_un_fun=eggnog$ct_Domain_un_fun,
  EggNOG_Med_CAZy = eggnog_med_ct_em_CAZy,
  EggNOG_Med_COG_Cat = eggnog_med_ct_em_COG_cat,
  EggNOG_Med_EC = eggnog_med_ct_em_EC,
  EggNOG_Med_KEGG_ko = eggnog_med_ct_em_KEGG_ko,
  EggNOG_Med_BRITE = eggnog_med_ct_em_BRITE,
  EggNOG_Med_BiGG_Reaction = eggnog_med_ct_em_BiGG_Reaction,
  EggNOG_Med_KEGG_rclass = eggnog_med_ct_em_KEGG_rclass,
  EggNOG_Med_KEGG_TC = eggnog_med_ct_em_KEGG_TC,
  EggNOG_Med_KEGG_Pathway = eggnog_med_ct_em_KEGG_Pathway,
  EggNOG_Med_KEGG_Reaction = eggnog_med_ct_em_KEGG_Reaction,
  EggNOG_Med_KEGG_Module = eggnog_med_ct_em_KEGG_Module,
  EggNOG_Med_PFAMs = eggnog_med_ct_em_PFAMs,
  PGAP_Med_RefSeq = pgap_med_refseq,
  PGAP_Med_GO_Fun = pgap_med_go_fun,
  PGAP_Med_GO_Pro = pgap_med_go_pro,
  PGAP_Med_GO_Comp = pgap_med_go_comp
)

write_csv(summary_df,"comparison.csv")
