#title: "normalising_and_aggregating_nitrifier_genes"
#author: "Kalinka Sand Knudsen"
#update: "2024-07-30"


#Loading packages
.libPaths(c("path_to_R_package_library", .libPaths()))
library(vroom)
library(tidyverse)
library(readxl)

#Set WD
setwd("path_to_working_directory")


#Loading bacterial amoA graftM output
amoA_bac <- vroom("bac_amoA_combined_count_table_e10.txt", delim = "\t") %>%
rename(OTU = ConsensusLineage)  %>%
mutate(OTU = if_else(OTU == "Root", "Root_bac", OTU))

#Linking sample IDs
linking<-vroom("2023-10-11_multi_files.csv", delim=",",  col_select = c(1,2), col_names = c("seq1","seq2" ))%>%
  mutate(seq_com = paste0(seq1,'_v_', seq2))

# Loop through the combining_info data frame
for (i in 1:nrow(linking)) {
  col1 <- linking$seq1[i]
  col2 <- linking$seq2[i]
  new_name <- linking$seq_com[i]
  
  # Use dplyr to create a new column with the combined values
  amoA_bac <- amoA_bac %>%
    mutate(!!new_name := rowSums(select(., col1, col2)))
  
  # Remove the original columns if needed
  amoA_bac <- amoA_bac %>%
    select(-col1, -col2)
}




# Importing the microbial read fraction:
frac <- vroom("output_otu_profile_microbial_frac_stdout_merged.txt", delim = "\t")%>%
  select(sample, microbial_fraction)%>%
  mutate(microbial_fraction = as.numeric(str_replace(microbial_fraction, "%", "")))%>%
  mutate(microbial_fraction=microbial_fraction/100)

#Importing and merging metadata
habitat <-read_excel("2024-02-13_mfd_db.xlsx")
seq_meta <- vroom("2023-10-11_samples_minimal_metadata_collapsed.csv", delim = ",") %>%
  mutate(flat_name=gsub(".fastq.gz","", flat_name))%>%
  rename(SeqId = flat_name) %>%
  relocate(SeqId)
com_meta<-merge(seq_meta, habitat, by="fieldsample_barcode")%>%
  relocate(SeqId)

com_meta_frac<-merge(com_meta, frac, by.x = "SeqId", by.y = "sample")

# Creating a dataframe to scale the read count
bac_scale<-com_meta_frac%>%
  mutate(per_million_scale=(after_total_reads/2)*microbial_fraction/1000000)%>% ###divide by 2 for only forward reads
  mutate(kb_length=(248*3)/1000)%>% #length of bacterial amoA = 248
  select(SeqId, per_million_scale, kb_length)

# Transpose the amoA_bac dataframe to have samples as rows and OTUs as columns
amoA_bac_transposed <- amoA_bac %>%
 select(-OTU) %>%
  t() %>%
  as.data.frame()

amoA_bac_transposed<-amoA_bac_transposed%>%
  mutate(SeqId=rownames(amoA_bac_transposed))%>%
  relocate(SeqId)
row.names(amoA_bac_transposed)=NULL
colnames(amoA_bac_transposed)[2:ncol(amoA_bac_transposed)] <- amoA_bac$OTU
  
# Merge the transposed amoA_bac dataframe with the bac_scale dataframe based on SeqId
merged_df_bac <- merge(bac_scale, amoA_bac_transposed, by = "SeqId", all.x = T)%>%
  replace(is.na(.), 0)

# Divide each column in the merged dataframe by the corresponding value in per_million_scale
scaled_otu_table_bac <- merged_df_bac %>%
  mutate(across(starts_with("Root"), ~./(per_million_scale*kb_length))) %>%
  select(-per_million_scale, -kb_length)

#saveRDS(scaled_otu_table_bac, file = "scaled_amoA_otu_table_bac.rds")


#Loading archaeal amoA graftM output
otutable_arc <- vroom("arch_amoA_combined_count_table_e10.txt", delim = "\t") %>%
rename(OTU = ConsensusLineage)  %>%
mutate(OTU = if_else(OTU == "Root", "Root_arc", OTU))

arc_scale<-com_meta_frac%>%
  mutate(per_million_scale=(after_total_reads/2)*microbial_fraction/1000000)%>%
  mutate(kb_length=(216*3)/1000)%>% #length of archaeal amoA = 216
  select(SeqId, per_million_scale, kb_length)

otutable_arc_transposed <- otutable_arc %>%
 select(-OTU) %>%
  t() %>%
  as.data.frame()

otutable_arc_transposed<-otutable_arc_transposed%>%
  mutate(SeqId=rownames(otutable_arc_transposed))%>%
  relocate(SeqId)
row.names(otutable_arc_transposed)=NULL
colnames(otutable_arc_transposed)[2:ncol(otutable_arc_transposed)] <- otutable_arc$OTU
  
# Merge the transposed otutable_arc dataframe with the arc_scale dataframe based on SeqId
merged_df_arc <- merge(arc_scale, otutable_arc_transposed, by = "SeqId", all.x = T)%>%
  replace(is.na(.), 0)

# Divide each column in the merged dataframe by the corresponding value in per_million_scale
scaled_otu_table_arc <- merged_df_arc %>%
  mutate(across(starts_with("Root"), ~./(per_million_scale*kb_length))) %>%
  select(-per_million_scale, -kb_length)

#saveRDS(scaled_otu_table_arc, file = "scaled_amoA_otu_table_arc.rds")




#Loading bacterial pNXR graftM output
otutable_pNXR <- vroom("pNXR_combined_count_table_e10.txt", delim = "\t") %>%
rename(OTU = ConsensusLineage)  %>%
mutate(OTU = if_else(OTU == "Root", "Root_pNXR", OTU))

pNXR_scale<-com_meta_frac%>% 
  mutate(per_million_scale=(after_total_reads/2)*microbial_fraction/1000000)%>%
  mutate(kb_length=(1137*3)/1000)%>% #length of pNXR search HMM = 1137
  select(SeqId, per_million_scale, kb_length)

otutable_pNXR_transposed <- otutable_pNXR %>%
 select(-OTU) %>%
  t() %>%
  as.data.frame()

otutable_pNXR_transposed<-otutable_pNXR_transposed%>%
  mutate(SeqId=rownames(otutable_pNXR_transposed))%>%
  relocate(SeqId)
row.names(otutable_pNXR_transposed)=NULL
colnames(otutable_pNXR_transposed)[2:ncol(otutable_pNXR_transposed)] <- otutable_pNXR$OTU
  
# Merge the transposed otutable_pNXR dataframe with the pNXR_scale dataframe based on SeqId
merged_df_pNXR <- merge(pNXR_scale, otutable_pNXR_transposed, by = "SeqId", all.x = T)%>%
  replace(is.na(.), 0)

# Divide each column in the merged dataframe by the corresponding value in per_million_scale
scaled_otu_table_pNXR <- merged_df_pNXR %>%
  mutate(across(starts_with("Root"), ~./(per_million_scale*kb_length))) %>%
  select(-per_million_scale, -kb_length)

#saveRDS(scaled_otu_table_pNXR, file = "scaled_pNXR_otu_table.rds")





###### And the cNXR #####
otutable_cNXR <- vroom("cNXR_combined_count_table_e10.txt", delim = "\t") %>%
rename(OTU = ConsensusLineage)  %>%
mutate(OTU = if_else(OTU == "Root", "Root_cNXR", OTU))

cNXR_scale<-com_meta_frac%>% 
  mutate(per_million_scale=(after_total_reads/2)*microbial_fraction/1000000)%>%
  mutate(kb_length=(1364*3)/1000)%>% #length of cNXR search HMM = 1364
  select(SeqId, per_million_scale, kb_length)

otutable_cNXR_transposed <- otutable_cNXR %>%
 select(-OTU) %>%
  t() %>%
  as.data.frame()

otutable_cNXR_transposed<-otutable_cNXR_transposed%>%
  mutate(SeqId=rownames(otutable_cNXR_transposed))%>%
  relocate(SeqId)
row.names(otutable_cNXR_transposed)=NULL
colnames(otutable_cNXR_transposed)[2:ncol(otutable_cNXR_transposed)] <- otutable_cNXR$OTU
  
# Merge the transposed otutable dataframe with the scale dataframe based on SeqId
merged_df_cNXR <- merge(cNXR_scale, otutable_cNXR_transposed, by = "SeqId", all.x = T)%>%
  replace(is.na(.), 0)

# Divide each column in the merged dataframe by the corresponding value in per_million_scale
scaled_otu_table_cNXR <- merged_df_cNXR %>%
  mutate(across(starts_with("Root"), ~./(per_million_scale*kb_length))) %>%
  select(-per_million_scale, -kb_length)

#saveRDS(scaled_otu_table_cNXR, file = "scaled_cNXR_otu_table.rds")




### Merging the dataframes

habitat <-read_excel("2024-02-13_mfd_db.xlsx")
seq_meta <- vroom("2023-10-11_samples_minimal_metadata_collapsed.csv", delim = ",") %>%
  mutate(flat_name=gsub(".fastq.gz","", flat_name))%>%
  rename(SeqId = flat_name) %>%
  relocate(SeqId)
com_meta<-merge(seq_meta, habitat, by="fieldsample_barcode")%>%
  relocate(SeqId)

rm(habitat, seq_meta)
otu_meta_bac<-merge(scaled_otu_table_bac, com_meta, by="SeqId")
otu_merged<-merge(otu_meta_bac, scaled_otu_table_arc, by="SeqId")
otu_merged<-merge(otu_merged, scaled_otu_table_pNXR, by="SeqId")
otu_merged<-merge(otu_merged, scaled_otu_table_cNXR, by="SeqId")


#### Cleaning up the taxonomy 

otu_filtered <- otu_merged %>%
  select(
    matches("^Root; amoA;"),
    matches("^Root; o_Nitrososphaerales"),
    matches("^Root; Root_narG; Thiocapsa"),
    matches("^Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade"),
    matches("^Root; cNarG_Nxr; Methylomirabilis_acidimicrobia"),
    matches("^Root; Nitrotoga"),
    matches("^Root; periplasmic_umbrella; Nitrohelix_Nitrospina_Nitrospiraceae_cluster"),
    matches("^Root; periplasmic_umbrella; Nitrospira_umbrella"),
    matches("^Root; periplasmic_umbrella; NS_4_NS_12"),
    matches("^Root; periplasmic_umbrella; UBA8639"),
    matches("Seq|mfd")
  )

OTU_filtered_long<-pivot_longer(otu_filtered, starts_with("Root;"), values_to = "RPKM", names_to = "Tax")%>%
  mutate(Tax_short=trimws(str_extract(Tax, "[^;]*$")))%>%
  mutate(type=if_else(grepl("periplasmic_umbrella|Nitrotoga", Tax), "pNXR", "Bacterial amoA"))%>%
  mutate(type=if_else(grepl("Root; cNarG_Nxr|Thiocapsa", Tax), "cNXR", type))%>%
  mutate(type=if_else(grepl("o_Nitrososphaerales", Tax),"Archaeal amoA", type))


#### Aggregating groups, and removing groups that are found in nearly no samples for easier interpretability
tax_remove=c("Root; o_Nitrososphaerales; f_Nitrosopumilaceae; Nitrosotalea_TA20_cluster", "Root; o_Nitrososphaerales; f_Nitrososphaeraceae; TH5893_Nitrososphaera_cluster; JAJPHM01_Nitrososphaera_cluster", "Root; periplasmic_umbrella; Nitrohelix_Nitrospina_Nitrospiraceae_cluster", "Root; periplasmic_umbrella; Nitrohelix_Nitrospina_Nitrospiraceae_cluster; f_Nitrospinaceae_unknown_cluster", "Root; periplasmic_umbrella; Nitrohelix_Nitrospina_Nitrospiraceae_cluster; Nitrohelix_Nitrospina", "Root; periplasmic_umbrella; Nitrohelix_Nitrospina_Nitrospiraceae_cluster; Nitrohelix_Nitrospina; Nitrospina", "Root; periplasmic_umbrella; Nitrospira_umbrella; VGXG01_SCVQ01_clade" )

otu_filtered2 <- otu_filtered %>%
  mutate("Root; amoA; Gammaproteobacteria_amoA; Nitrosomonas" = rowSums(select(., matches("Nitrosomonas_clade_")), na.rm = TRUE)) %>%
  mutate("Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; Nitrobacter-like_NxrA; BOG-931" = rowSums(select(., matches("BOG-931_")), na.rm = TRUE)) %>%
  mutate("Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; Nitrobacter-like_NxrA; Bradyrhizobium_Nitrobacter; Bradyrhizobium" = rowSums(select(., matches("Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; Nitrobacter-like_NxrA; Bradyrhizobium_Nitrobacter; Brady")), na.rm = TRUE)) %>%
  mutate("Root; o_Nitrososphaerales; Other Nitrosopumilaceae" = rowSums(select(., matches("g_CSP1|g_DRGT01_1|g_JACEMX01|g_Nitrosarchaerum|g_Nitrosopelagicus|g_Nitrosopumilus|g_Nitrosotenuis|g_UBA8516|g_TA20|Root; o_Nitrososphaerales; f_Nitrosopumilaceae_cluster")), na.rm = TRUE)) %>%
  mutate("Root; o_Nitrososphaerales; Other Nitrosopumilaceae" = rowSums(select(.,c("Root; o_Nitrososphaerales; Other Nitrosopumilaceae", "Root; o_Nitrososphaerales; f_Nitrosopumilaceae")), na.rm = TRUE)) %>%
  mutate("Root; o_Nitrososphaerales; f_Nitrososphaeraceae; TH5893_Nitrososphaera_cluster; g_TH5893" = rowSums(select(., matches("Root; o_Nitrososphaerales; f_Nitrososphaeraceae; TH5893_Nitrososphaera_cluster; g_TH5893_")), na.rm = TRUE)) %>%
  mutate("Root; o_Nitrososphaerales; f_Nitrososphaeraceae; TH5896" = rowSums(select(., matches("Root; o_Nitrososphaerales; f_Nitrososphaeraceae; g_TH5896_")), na.rm = TRUE)) %>%
  mutate("Root; o_Nitrososphaerales; f_Nitrososphaeraceae; g_TA21" = rowSums(select(., matches("g_TA21_")), na.rm = TRUE)) %>%
  mutate("Root; periplasmic_umbrella; Nitrospira_umbrella; Nitrospira_C" = rowSums(select(., matches("Nitrospira_C1|Nitrospira_C2")), na.rm = TRUE)) %>%  
  mutate("Root; periplasmic_umbrella; Nitrospira_umbrella; Nitrospira_A" = rowSums(select(., matches("Nitrospira_A")), na.rm = TRUE)) %>%
  mutate("Root; periplasmic_umbrella; Nitrospira_umbrella; Nitrospira_F" = rowSums(select(., matches("Nitrospira_F")), na.rm = TRUE)) %>%
  mutate("Root; periplasmic_umbrella; Nitrospira_umbrella; Palsa-1315" = rowSums(select(., matches("Palsa-1315")), na.rm = TRUE)) %>%  
  select(-matches("Nitrosomonas_clade_"), -matches("BOG-931_"), -matches("Bradyrhizobium_1"), -matches("Bradyrhizomium_2"),
         -matches("g_CSP1|g_DRGT01_1|g_JACEMX01|g_Nitrosarchaerum|g_Nitrosopelagicus|g_Nitrosopumilus|g_Nitrosotenuis|g_UBA8516|g_TA20"), 
         -matches("Root; o_Nitrososphaerales; f_Nitrososphaeraceae; TH5893_Nitrososphaera_cluster; g_TH5893_"), 
         -matches("Root; o_Nitrososphaerales; f_Nitrososphaeraceae; g_TH5896_"), 
         -matches("g_TA21_"), 
         -matches("Nitrospira_C1|Nitrospira_C2"), 
         -matches("Nitrospira_A1|Nitrospira_A2"), 
         -matches("Nitrospira_F1|Nitrospira_F2"), 
         -matches("Palsa-1315_"),
         -matches("Root; o_Nitrososphaerales; f_Nitrosopumilaceae_cluster"))%>%
  select(!"Root; o_Nitrososphaerales; f_Nitrosopumilaceae")
  

### Testing aggregation
test1<-otu_filtered2%>%
  filter(SeqId %in% c("LIB-MJ044-A1_03", "LIB-MJ166-B1_01", "LIB-MJ176-C4_03", "LIB-MJ172-B5_04", "LIB-MJ260-D1_01"))%>%
  select(SeqId, "Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; Nitrobacter-like_NxrA; Bradyrhizobium_Nitrobacter; Bradyrhizobium")

test2<-otu_filtered%>%
  filter(SeqId %in% c("LIB-MJ044-A1_03", "LIB-MJ166-B1_01", "LIB-MJ176-C4_03", "LIB-MJ172-B5_04", "LIB-MJ260-D1_01"))%>%
  select(., matches("Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; Nitrobacter-like_NxrA; Bradyrhizobium_Nitrobacter; Brady"))


OTU_filtered_long<-pivot_longer(otu_filtered2, starts_with("Root;"), values_to = "RPKM", names_to = "Tax")%>%
  filter(!Tax %in% tax_remove)%>%
  mutate(Tax_short=trimws(str_extract(Tax, "[^;]*$")))%>%
  mutate(type=if_else(grepl("periplasmic_umbrella|Nitrotoga", Tax), "pNxrA", "Bacterial amoA"))%>%
  mutate(type=if_else(grepl("Root; cNarG_Nxr|Thiocapsa", Tax), "cNxrA", type))%>%
  mutate(type=if_else(grepl("o_Nitrososphaerales", Tax),"Archaeal amoA", type))


#### To use .rds object in other scripts
saveRDS(OTU_filtered_long, "OTU_filtered_long.rds")




