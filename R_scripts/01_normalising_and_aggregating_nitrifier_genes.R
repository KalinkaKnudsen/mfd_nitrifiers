#title: "01_normalising_and_aggregating_nitrifier_genes"
#author: "Kalinka Sand Knudsen"
#update: "2025-04-23"


#Loading packages
.libPaths(c("path_to_R_package_library", .libPaths()))
library(vroom)
library(tidyverse)
library(readxl)

#Set WD
setwd("path_to_working_directory")


#Loading bacterial amoA graftM output
amoA_bac <- vroom("./data/bac_amoA_combined_count_table_e10.txt", delim = "\t") %>%
rename(OTU = ConsensusLineage)  %>%
mutate(OTU = if_else(OTU == "Root", "Root_bac", OTU))

#Linking sample IDs
linking<-vroom("./metadata/2023-10-11_multi_files.csv", delim=",",  col_select = c(1,2), col_names = c("seq1","seq2" ))%>%
  mutate(seq_com = paste0(seq1,'_v_', seq2))

# Loop through the combining_info data frame
for (i in 1:nrow(linking)) {
  col1 <- linking$seq1[i]
  col2 <- linking$seq2[i]
  new_name <- linking$seq_com[i]
  
  # Use dplyr to create a new column with the combined values
  amoA_bac <- amoA_bac %>%
    mutate(!!new_name := rowSums(select(., col1, col2)))
  
  # Remove the original columns 
  amoA_bac <- amoA_bac %>%
    select(-col1, -col2)
}




habitat <-read_excel("./metadata/2025-02-19_mfd_db.xlsx")
seq_meta <- vroom("./metadata/2023-10-11_samples_minimal_metadata_collapsed.csv", delim = ",") %>%
  mutate(flat_name=gsub(".fastq.gz","", flat_name))%>%
  rename(SeqId = flat_name) %>%
  relocate(SeqId)
com_meta<-merge(seq_meta, habitat, by="fieldsample_barcode")%>%
  relocate(SeqId)

bac_scale<-com_meta%>%
  mutate(per_million_scale=(after_total_reads/2)/1000000)%>% ###divide by 2 for only forward reads
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
  
# Merge the transposed amoA_bac dataframe with the arc_scale dataframe based on SeqId
merged_df_bac <- merge(bac_scale, amoA_bac_transposed, by = "SeqId", all.x = T)%>%
  replace(is.na(.), 0)

# Divide each column in the merged dataframe by the corresponding value in per_million_scale
scaled_otu_table_bac <- merged_df_bac %>%
  mutate(across(starts_with("Root"), ~./(per_million_scale*kb_length))) %>%
  select(-per_million_scale, -kb_length)

saveRDS(scaled_otu_table_bac, file = "output/scaled_amoA_otu_table_bac.rds")





#### Now to the archaeal 


otutable_arc <- vroom("data/arch_amoA_combined_count_table_e10.txt", delim = "\t") %>%
rename(OTU = ConsensusLineage)  %>%
mutate(OTU = if_else(OTU == "Root", "Root_arc", OTU))

arc_scale<-com_meta%>%
  mutate(per_million_scale=(after_total_reads/2)/1000000)%>%
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

saveRDS(scaled_otu_table_arc, file = "output/scaled_amoA_otu_table_arc.rds")




###### And the p-nxr:

otutable_pNXR <- vroom("./data/pNXR_combined_count_table_e10.txt", delim = "\t") %>%
rename(OTU = ConsensusLineage)  %>%
mutate(OTU = if_else(OTU == "Root", "Root_pNXR", OTU))

pNXR_scale<-com_meta%>% 
  mutate(per_million_scale=(after_total_reads/2)/1000000)%>%
  mutate(kb_length=(1137*3)/1000)%>% #length of nxr search HMM = 1137
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
  
# Merge the transposed otutable_arc dataframe with the arc_scale dataframe based on SeqId
merged_df_pNXR <- merge(pNXR_scale, otutable_pNXR_transposed, by = "SeqId", all.x = T)%>%
  replace(is.na(.), 0)

# Divide each column in the merged dataframe by the corresponding value in per_million_scale
scaled_otu_table_pNXR <- merged_df_pNXR %>%
  mutate(across(starts_with("Root"), ~./(per_million_scale*kb_length))) %>%
  select(-per_million_scale, -kb_length)

saveRDS(scaled_otu_table_pNXR, file = "output/scaled_pNXR_otu_table.rds")





###### And the cNXR #####

otutable_cNXR <- vroom("./data/cNXR_combined_count_table_e10.txt", delim = "\t") %>%
rename(OTU = ConsensusLineage)  %>%
mutate(OTU = if_else(OTU == "Root", "Root_cNXR", OTU))

cNXR_scale<-com_meta%>% 
  mutate(per_million_scale=(after_total_reads/2)/1000000)%>%
  mutate(kb_length=(1364*3)/1000)%>% #length of nxr search HMM = 1364
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
  
# Merge the transposed otutable_arc dataframe with the arc_scale dataframe based on SeqId
merged_df_cNXR <- merge(cNXR_scale, otutable_cNXR_transposed, by = "SeqId", all.x = T)%>%
  replace(is.na(.), 0)

# Divide each column in the merged dataframe by the corresponding value in per_million_scale
scaled_otu_table_cNXR <- merged_df_cNXR %>%
  mutate(across(starts_with("Root"), ~./(per_million_scale*kb_length))) %>%
  select(-per_million_scale, -kb_length)

saveRDS(scaled_otu_table_cNXR, file = "output/scaled_cNXR_otu_table.rds")




### Merging the dataframes

habitat <-read_excel("./metadata/2025-02-19_mfd_db.xlsx")
seq_meta <- vroom("./metadata/2023-10-11_samples_minimal_metadata_collapsed.csv", delim = ",") %>%
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

#saveRDS(otu_merged, file = "output/merged_otu_table_25_02_21.rds")

#otu_merged<-readRDS("output/merged_otu_table_25_02_21.rds")
habitat <-read_excel("./metadata/2025-02-19_mfd_db.xlsx")
seq_meta <- vroom("/metadata/2023-10-11_samples_minimal_metadata_collapsed.csv", delim = ",") %>%
  mutate(flat_name=gsub(".fastq.gz","", flat_name))%>%
  rename(SeqId = flat_name) %>%
  relocate(SeqId)
com_meta<-merge(seq_meta, habitat, by="fieldsample_barcode")%>%
  relocate(SeqId)



########## Now filtering on samples with at least 500000 reads after processing #####
png("filtering_of_reads.png", width=3000, height=1600, res=350)
jpeg("filtering_of_reads.jpeg", width=3000, height=1600, res=350)
hist(com_meta$after_total_reads, breaks=1500, xlim=c(0,50000000), 
     xlab="Total reads after filtering", main="Histogram of reads per sample")
abline(v=500000, col="red3", lwd=2)
# Add a label for the red line
text(500000, 150, "cutoff: 500 000 reads", col="grey20", pos=4, cex=0.8)
dev.off()

 
 too_few_reads<-com_meta%>%filter(after_total_reads<500000)%>%pull(SeqId)

 
otu_merged<-otu_merged%>%filter(!SeqId %in% too_few_reads)
  

otu_filtered <- otu_merged %>%
  select(
    matches("^Root; amoA;"),
    matches("^Root; Nitrosococcus"),
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




### First removing all that is NA in mfd_hab1. It is 77 samples natural soil
OTU_filtered_long%>%select(mfd_sampletype, mfd_areatype)%>%distinct()

OTU_filtered_long_plot<-OTU_filtered_long%>%
  filter(!is.na(mfd_hab1))%>%
  mutate(mfd_areatype=gsub("Subterranean", "S", mfd_areatype),
         mfd_areatype=gsub("Natural", "N", mfd_areatype),
         mfd_areatype=gsub("Urban","U", mfd_areatype),
         mfd_areatype=gsub("Agriculture \\(reclaimed lowland\\)", "A(r)", mfd_areatype),
         mfd_areatype=gsub("Agriculture","A", mfd_areatype))


hab1_order<-OTU_filtered_long_plot %>%
  arrange(mfd_hab1)%>%
  select(SeqId)%>%distinct()%>%
  pull(SeqId)

heat <- OTU_filtered_long_plot %>%
  mutate(SeqId = factor(SeqId, levels = hab1_order, ordered = TRUE)) %>%
 # filter(mfd_hab1 %in% c("Freshwater", "Grassland formations", "Forests", "Greenspaces", "Fields", "Bogs, mires and fens", "Temperate heath and scrub", "Sclerophyllous scrub", "Coastal", "Dunes"))%>%
 # mutate(mfd_hab1 = if_else(grepl("scrub", mfd_hab1), "Scrub", mfd_hab1))%>%
  ggplot(aes(y = Tax, x = SeqId, fill = RPKM)) +
  geom_tile(width=10.0, linewidth=0.0) +
  scale_fill_gradientn(
        name = "RPKM",  # The label
        colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  
        trans = "sqrt",  # Same as in viridis
        limits = c(0, 5),  # Set the limits manually
        na.value = "black"  # Handling values larger than the limit
      ) +
  facet_nested(type ~ mfd_sampletype+mfd_areatype, scales = "free", space = "free", 
               nest_line = element_line(color="grey20", linewidth = 0.4), resect=unit(2, "pt"),
               , switch = "y", strip = strip_nested(clip = "off")) +
  labs(x = "", y = "") +
 # theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        #  axis.text.y = element_blank(),
        axis.text.y = element_text(size = 7, margin = margin(r = -2)),
        legend.key.size = unit(0.45, "cm"),
        legend.title = element_text(size=7),
        legend.text = element_text(size=7),
        legend.position = "right",
        strip.clip = "off",
        strip.text.x = element_text(size = 7, face="bold"),  # Size for x-axis facet labels
        strip.text.y = element_text(size = 7, face="bold"),
        panel.border = element_rect(colour="black", fill=NA, linewidth = 0.1) ,
        strip.background = element_blank(),
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.1, "lines", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))



#heat


habitat_colors <- c(
  "Freshwater" = "cyan3",
  "Grassland formations" = "#33a02c",
  "Rocky habitats and caves" = "#e31a1c",
  "Greenspaces" = "lightgreen",
  "Fields" = "#b2df8a",
  "Dunes" = "#b15928",
  "Roadside" = "#a6cee3",
  "Bogs, mires and fens" = "darkgreen",
  "Sclerophyllous scrub" = "purple",
  "Temperate heath and scrub" = "gold3",
  "Coastal" = "#C2B280",
  "Saltwater" = "darkblue",
  "Forests" = "#8dd3c7",
  "Urban" = "darkred",
  "Wastewater" = "purple4",
  "Drinking water" = "#fb8072",
  "Other" = "#ff7f00",
  "Sandfilter" = "#fdb462",
  "Biogas" = "#b3de69"
)


bar_plot <- OTU_filtered_long_plot %>%
  mutate(SeqId = factor(SeqId, levels = hab1_order, ordered = TRUE)) %>%
  ggplot(aes(x = SeqId, fill = mfd_hab1)) +
  geom_tile(aes(y = 1), width=10.0, linewidth=0.0) +
  facet_grid(. ~ mfd_sampletype+mfd_areatype, scales = "free", space = "free") +
  scale_fill_manual(values=habitat_colors)+
  theme(        legend.key.size = unit(0.3, "cm"),
                legend.text = element_text(size=7),
                legend.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "bottom",
    panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2) ,
    panel.spacing = unit(0.1, "lines", data = NULL)
  )+   scale_y_continuous(expand = c(0,0)) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))



combined_plot <- heat / bar_plot + plot_layout(heights = c(20, 0.35)) + theme(panel.background = element_blank(), plot.background = element_blank())



ggsave("output/heatmap_uncollapsed_all_25_02_21.jpeg",
       combined_plot,
       height = 11,
       width = 22, dpi=600)




#ggsave("heatmap_uncollapsed.png", heat, height=14, width=27, dpi=600)


##################### Now filtering based on the image above ######
# 
# ##Groups with more than one clade, where only a single clade is significant:
# g_TH5893_2
# g_TA21_3
# g_TH5896_2
# g_TA20
# Nitrosotalea_TA20_cluster
# g_UBA8516
# g_Nitrosotenuis_3
# g_Nitrosopumilus_2
# g_Nitrosopelagicus_2
# g_JACEMX01
# g_DRGT01_1
# g_CSP1
# 
# #### AOB ####
# Nitrosomonas (only needs to be aggregated)
# Nitrosococcus (no hits)
# 
# #### cNXR ####
# BOG-931_2 and BOG-931_3 collapse to one
# Methylomirabilis_acidimicrobia 
# 
# ### pNXR ####
# VGXG01_SCVQ01_clade
# Nitrospira_C2
# Nitrohelix_Nitrospina
# Nitrohelix_Nitrospina_Nitrospiraceae_cluster
# f_Nitrospinaceae_unknown_cluster
# Nitrospina


tax_remove<-c("Nitrosococcus", "g_TH5893_2","g_TA21_3","g_TH5896_2","g_TA20","Nitrosotalea_TA20_cluster","g_UBA8516","g_Nitrosotenuis_3","g_Nitrosopumilus_2","g_Nitrosopelagicus_2","g_JACEMX01","g_DRGT01_1","g_CSP1",
              "Methylomirabilis_acidimicrobia", "VGXG01_SCVQ01_clade", "Nitrospira_C2", "Nitrohelix_Nitrospina", "Nitrohelix_Nitrospina_Nitrospiraceae_cluster", 
              "f_Nitrospinaceae_unknown_cluster", "Nitrospina", "JAJPHM01_Nitrososphaera_cluster")



#otu_filtered<-readRDS("output/OTU_filtered_long_24_03_15.rds")

#### Now to collapsing some groups.

otu_filtered2 <- otu_filtered %>%
  mutate("Root; amoA; Gammaproteobacteria_amoA; Nitrosomonas" = rowSums(select(., matches("Nitrosomonas_clade_")), na.rm = TRUE)) %>%
  mutate("Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; Nitrobacter-like_NxrA; BOG-931" = rowSums(select(., matches("BOG-931_")), na.rm = TRUE)) %>%
  mutate("Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; Nitrobacter-like_NxrA; Bradyrhizobium_Nitrobacter; Bradyrhizobium" = rowSums(select(., matches("Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; Nitrobacter-like_NxrA; Bradyrhizobium_Nitrobacter; Brady")), na.rm = TRUE)) %>%
  mutate("Root; o_Nitrososphaerales; f_Nitrososphaeraceae; TA21_Nitrosopolaris_cluster; g_TA21" = rowSums(select(., matches("Root; o_Nitrososphaerales; f_Nitrososphaeraceae; TA21_Nitrosopolaris_cluster; g_TA21")), na.rm = TRUE)) %>%
  mutate("Root; periplasmic_umbrella; Nitrospira_umbrella; Nitrospira_clade_1; Nitrospira_A" = rowSums(select(., matches("Nitrospira_A")), na.rm = TRUE)) %>%
  mutate("Root; periplasmic_umbrella; Nitrospira_umbrella; Nitrospira_clade_2; Nitrospira_D_lenta" = rowSums(select(., matches("Nitrospira_D")), na.rm = TRUE)) %>%
  mutate("Root; periplasmic_umbrella; Nitrospira_umbrella; Nitrospira_clade_2; Nitrospira_D_inopinata" = rowSums(select(., matches("Nitrospira_F")), na.rm = TRUE)) %>%
  mutate("Root; periplasmic_umbrella; Nitrospira_umbrella; Nitrospira_clade_2; Palsa-1315" = rowSums(select(., matches("Palsa-1315")), na.rm = TRUE)) %>%
  mutate("Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; Nitrobacter-like_NxrA; Pseu_JAFA_VAZQ_clade" = rowSums(select(., c("Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; Nitrobacter-like_NxrA; Pseudolabrys_JAFA_VAZQ_umbrella; JAFAXD01_VAZQ01_Pseudolabrys_cluster", "Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; Nitrobacter-like_NxrA; Pseudolabrys_JAFA_VAZQ_umbrella"))))%>%
  select(-matches("Nitrosomonas_clade_"), -matches("BOG-931_"), -matches("Bradyrhizobium_1"), -matches("Bradyrhizomium_2"),
         -matches("Root; o_Nitrososphaerales; f_Nitrososphaeraceae; TA21_Nitrosopolaris_cluster; g_TA21_"), 
         -matches("Nitrospira_A1|Nitrospira_A2"), 
         -matches("Nitrospira_F1|Nitrospira_F2"), 
         -matches("Palsa-1315_"))%>%
  select(-c("Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; Nitrobacter-like_NxrA; Pseudolabrys_JAFA_VAZQ_umbrella; JAFAXD01_VAZQ01_Pseudolabrys_cluster",
         "Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; Nitrobacter-like_NxrA; Pseudolabrys_JAFA_VAZQ_umbrella", "Root; periplasmic_umbrella; Nitrospira_umbrella; Nitrospira_clade_2; Nitrospira_D"))
  


# test1<-otu_filtered2%>%
#   filter(SeqId %in% c("LIB-MJ044-A1_03", "LIB-MJ166-B1_01", "LIB-MJ176-C4_03", "LIB-MJ172-B5_04", "LIB-MJ260-D1_01"))%>%
#   select(SeqId, "Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; Nitrobacter-like_NxrA; Bradyrhizobium_Nitrobacter; Bradyrhizobium")
# 
# test2<-otu_filtered%>%
#   filter(SeqId %in% c("LIB-MJ044-A1_03", "LIB-MJ166-B1_01", "LIB-MJ176-C4_03", "LIB-MJ172-B5_04", "LIB-MJ260-D1_01"))%>%
#   select(., matches("Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; Nitrobacter-like_NxrA; Bradyrhizobium_Nitrobacter; Brady"))
# 

OTU_filtered_long<-pivot_longer(otu_filtered2, starts_with("Root;"), values_to = "RPKM", names_to = "Tax")%>%
  mutate(Tax_short=trimws(str_extract(Tax, "[^;]*$")))%>%
  filter(!Tax_short %in% tax_remove)%>%
  mutate(type=if_else(grepl("periplasmic_umbrella|Nitrotoga", Tax), "pNxrA", "Bacterial amoA"))%>%
  mutate(type=if_else(grepl("Root; cNarG_Nxr|Thiocapsa", Tax), "cNxrA", type))%>%
  mutate(type=if_else(grepl("o_Nitrososphaerales", Tax),"Archaeal amoA", type))

saveRDS(OTU_filtered_long, "output/OTU_filtered_long_25_02_21.rds")




heat <- OTU_filtered_long %>%
  mutate(SeqId = factor(SeqId, levels = hab1_order, ordered = TRUE)) %>%
  # filter(mfd_hab1 %in% c("Freshwater", "Grassland formations", "Forests", "Greenspaces", "Fields", "Bogs, mires and fens", "Temperate heath and scrub", "Sclerophyllous scrub", "Coastal", "Dunes"))%>%
  # mutate(mfd_hab1 = if_else(grepl("scrub", mfd_hab1), "Scrub", mfd_hab1))%>%
  ggplot(aes(y = Tax, x = SeqId, fill = RPKM)) +
  geom_tile(width=10.0, linewidth=0.0) +
  scale_fill_gradientn(
    name = "RPKM",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  
    trans = "sqrt",  # Same as in viridis
    limits = c(0, 5),  # Set the limits manually
    na.value = "black"  # Handling values larger than the limit
  ) +
  facet_nested(type ~ mfd_sampletype+mfd_areatype, scales = "free", space = "free", nest_line = element_line(color="grey20", linewidth = 0.4), resect=unit(2, "pt")) +
  labs(x = "", y = "") +
  # theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        #  axis.text.y = element_blank(),
        axis.text.y = element_text(size = 7),
        legend.key.size = unit(0.45, "cm"),
        legend.title = element_text(size=8),
        legend.text = element_text(size=7),
        legend.position = "right",
        strip.clip = "off",
        strip.text.x = element_text(size = 6.5, face="bold"),  # Size for x-axis facet labels
        strip.text.y = element_text(size = 7, face="bold"),
        panel.border = element_rect(colour="black", fill=NA, linewidth = 0.1) ,
        strip.background = element_blank(),
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.1, "lines", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))
