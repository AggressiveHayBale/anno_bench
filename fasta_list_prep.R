library(readr)
library(stringr)
library(dplyr)
bac120_metadata_r207 <- read.delim("/mnt/15cf8305-a1bc-486c-8f97-fe85e5f4d988/gits/genome_download/bac120_metadata/bac120_metadata_r207.tsv", header=TRUE)
setwd("/mnt/15cf8305-a1bc-486c-8f97-fe85e5f4d988/gits/genome_download/bac120_metadata/")


#Data check

table(bac120_metadata_r207$gtdb_representative)
quantile(bac120_metadata_r207$checkm_contamination)
quantile(bac120_metadata_r207$checkm_completeness)

sum(bac120_metadata_r207$checkm_contamination ==0)
sum(bac120_metadata_r207$checkm_contamination<=1 )
sum(bac120_metadata_r207$checkm_completeness >=99.5&bac120_metadata_r207$checkm_contamination ==0)

bac120_metadata_r207$ncbi_names_proc <- word(bac120_metadata_r207$ncbi_organism_name, 1, 2) 
bac120_metadata_r207$gtdb_representative
#################
# Count species #
#################
# species <- bac120_metadata_r207 %>%  filter(gtdb_representative=="t") %>% select(ncbi_names_proc) 
# species<- as.data.frame(table(species))

#species <- bac120_metadata_r207 %>%  filter(gtdb_representative=="t") %>%  select(ncbi_organism_name) 
#species <- bac120_metadata_r207 %>%  filter(gtdb_representative=="t") %>% select(ncbi_organism_name) 

######################
#1 sample per species#
######################
bac120_parsed <- bac120_metadata_r207
bac120_parsed$name <- str_extract(bac120_parsed$ncbi_taxonomy, "s__[a-zA-Z\\s]+")
bac120_parsed$name <- gsub("s__", "",bac120_parsed$name)

##1 sample per spec gtdb rep
bac120_parsed_gtdb_rep <- bac120_parsed 
bac120_parsed_gtdb_rep <- filter(bac120_parsed_gtdb_rep,gtdb_representative=="t") 
one_sample_per_name_rep<- distinct(bac120_parsed_gtdb_rep,name,.keep_all = TRUE)

#output check
#check <- as.data.frame(table(one_sample_per_name$name))
#setdiff(one_sample_per_name_rep$name,bac120_parsed$name)


###########
#FILE PREP#
###########

### To download
for_csv <- to_download[,c("ncbi_genbank_assembly_accession")]
write_csv(data.frame(for_csv),"download_rest.csv")


## For workflow 
for_workflow <- to_download[,c("ncbi_genbank_assembly_accession", "name")]
for_workflow$Path <- paste0(rep("/mnt/15cf8305-a1bc-486c-8f97-fe85e5f4d988/gits/genome_download/ncbi_dataset/",nrow(for_workflow)),for_workflow$ncbi_genbank_assembly_accession)

for_workflow$Noise<-rep("0.005",nrow(for_workflow))
for_workflow$Noise2<-rep("0.01",nrow(for_workflow))
for_workflow$Noise3<-rep("0.02",nrow(for_workflow))
colnames(for_workflow)<- c("Accession", "Species","Path","Noise","Noise2","Noise3")

write_csv(for_workflow,"list_workflow_ecoli_1batch.csv")



# File check
setwd("/mnt/15cf8305-a1bc-486c-8f97-fe85e5f4d988/gits/genome_download/e_coli_first_batch/ncbi_dataset")
file_check<-list.files(pattern="GCA*")






# Just E. coli: 
e_coli <-bac120_metadata_r207
e_coli$name <- str_extract(e_coli$ncbi_taxonomy, "s__[a-zA-Z\\s]+")

e_coli$name <- gsub("s__", "",e_coli$name)

e_coli <- e_coli %>%  filter(name=="Escherichia coli")


# 24415/2 = 12207/2 = 6104 
to_download <- e_coli[6104:24415 ,]






#Species count
#count<-as.data.frame(table(bac120_parsed$name))
# count_with_rep<- bac120_parsed %>% filter(gtdb_representative=="t") 
# count_with_rep<- as.data.frame(table(count_with_rep$name))

# first_15<- bac120_metadata_r207[1:15,]
# for_csv <- first_15[,c("ncbi_genbank_assembly_accession","ncbi_organism_name")]
# write_delim(data.frame(for_csv[,1]),"list.tsv")

# for_workflow <- for_csv
# for_csv$ncbi_organism_name <- word(for_csv$ncbi_organism_name, 1, 2)
# colnames(for_workflow)<- c("Accession", "Species")





species_parsed <- str_extract(bac120_metadata_r207$ncbi_taxonomy, "s__[a-zA-Z\\s]+")
sum(is.na(species_parsed))


# Just E. coli: 
e_coli <-bac120_metadata_r207
e_coli$name <- str_extract(e_coli$ncbi_taxonomy, "s__[a-zA-Z\\s]+")

e_coli$name <- gsub("s__", "",e_coli$name)

e_coli <- e_coli %>%  filter(name=="Escherichia coli")


# 24415/2 = 12207/2 = 6104 
to_download <- e_coli[6104:24415 ,]



# # Enterobacter
# e_coli <- e_coli1 %>%  filter(name=="Escherichia coli")
# entero <- e_coli1 %>%  filter( name=="Enterobacter cloacae") 
# kleb <- e_coli1 %>%  filter( name=="Klebsiella pneumoniae") 
# citro <- e_coli1 %>%  filter( name=="Citrobacter freundii")
# pseudo <- e_coli1 %>%  filter( name=="Pseudomonas aeruginosa")
# acineto <- e_coli1 %>%  filter( name=="Acinetobacter baumannii")






### To download
for_csv <- to_download[,c("ncbi_genbank_assembly_accession")]
write_csv(data.frame(for_csv),"download_rest.csv")


## For workflow 
for_workflow <- to_download[,c("ncbi_genbank_assembly_accession", "name")]
for_workflow$Path <- paste0(rep("/mnt/15cf8305-a1bc-486c-8f97-fe85e5f4d988/gits/genome_download/ncbi_dataset/",nrow(for_workflow)),for_workflow$ncbi_genbank_assembly_accession)

for_workflow$Noise<-rep("0.005",nrow(for_workflow))
for_workflow$Noise2<-rep("0.01",nrow(for_workflow))
for_workflow$Noise3<-rep("0.02",nrow(for_workflow))
colnames(for_workflow)<- c("Accession", "Species","Path","Noise","Noise2","Noise3")

write_csv(for_workflow,"list_workflow_ecoli_1batch.csv")



# File check
setwd("/mnt/15cf8305-a1bc-486c-8f97-fe85e5f4d988/gits/genome_download/e_coli_first_batch/ncbi_dataset")
file_check<-list.files(pattern="GCA*")
