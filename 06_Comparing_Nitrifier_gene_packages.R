#title: "06_Comparing_Nitrifier_gene_packages"
#author: "Kalinka Sand Knudsen"
#update: "2024-07-31"


#Loading packages
.libPaths(c("path_to_R_package_library", .libPaths()))
library(ggplot2)
library(vroom)
library(tidyverse)
library(readxl)
library(ggh4x)


#Set WD
setwd("path_to_working_directory")

## Importing dataframes
cNXR_old<-vroom("cNXR_old_combined_count_table_e10.txt", delim="\t")%>%
  rename(OTU = ConsensusLineage)  %>%
  filter(grepl("Root; cNarG_Nxr; Nitrobacter_Nitrococcus", OTU)) %>%
  pivot_longer(starts_with("LIB"), names_to = "SeqId", values_to="Read_count")%>%
  mutate(Package="- MFD")%>%
  mutate(Group="Nitrobacter-\nlike nxrA")

cNXR_new<-vroom("cNXR_combined_count_table_e10.txt", delim="\t")%>%
  rename(OTU = ConsensusLineage)  %>%
  filter(grepl("Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade", OTU))%>%
 pivot_longer(starts_with("LIB"), names_to = "SeqId", values_to="Read_count")%>%
  mutate(Package="+ MFD")%>%
  mutate(Group="Nitrobacter-\nlike nxrA")

Arc_old<-vroom("arch_AmoA_old_combined_count_table_e10.txt", delim="\t")%>%
  rename(OTU = ConsensusLineage)  %>%
  filter(grepl("Root; Nitrososphaerales; Nitrososphaeraceae", OTU))%>%
  pivot_longer(starts_with("LIB"), names_to = "SeqId", values_to="Read_count")%>%
  mutate(Package="- MFD")%>%
  mutate(Group="Nitrososphaeraceae-\namoA")

Arc_new<-vroom("arch_AmoA_combined_count_table_e10.txt", delim="\t")%>%
  rename(OTU = ConsensusLineage)  %>%
  filter(grepl("Root; o_Nitrososphaerales; f_Nitrososphaeraceae", OTU))%>%
  pivot_longer(starts_with("LIB"), names_to = "SeqId", values_to="Read_count")%>%
  mutate(Package="+ MFD")%>%
  mutate(Group="Nitrososphaeraceae-\namoA")


#### combining into one
comb<-rbind(cNXR_old, cNXR_new)%>%
  rbind(Arc_old)%>%
  rbind(Arc_new)

##### Importing and merging metadata
linkage <- vroom("2023-10-11_samples_minimal_metadata_collapsed.csv", delim = ",") %>%
  mutate(flat_name=gsub(".fastq.gz","", flat_name))%>%
  rename(SeqId = flat_name) %>%
  relocate(SeqId)%>%
  select(SeqId, fieldsample_barcode)

metadata.sub <- readxl::read_excel('2024-02-13_mfd_db.xlsx') %>%
  left_join(linkage)%>%
  relocate(SeqId)%>%
  select(SeqId:mfd_hab3) %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", "))

d<-comb%>%
  left_join(metadata.sub)%>%
mutate(across(mfd_hab1, ~str_replace(., "Sclerophyllous scrub", "Temperate heath and scrub")),
       across(mfd_hab2, ~str_replace(., "Scrub", "Sclerophyllous scrub"))) %>%
  mutate(mfd_hab2=gsub("\\s*\\(non-habitat type\\)\\s*", "", mfd_hab2))%>%
  mutate(complex_long = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", ")) %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = "\n")) %>%
  group_by(complex) %>%
  mutate(complex_size = n()) %>%
  mutate(label = str_c(complex, ", n = ", complex_size, sep = "")) %>%
  ungroup()%>%
  filter(!complex_long %in% c("Water, Subterranean, Freshwater",
                         "Water, Urban, Sandfilter",
                         "Soil, Urban, Other",
                         "Sediment, Urban, Other",
                         "Soil, Urban, Roadside",
                         "Soil, Subterranean, Urban",
                         "Sediment, Urban, Saltwater",
                         "Sediment, Subterranean, Saltwater",
                         "Soil, Natural, Rocky habitats and caves",
                         "Water, Urban, Drinking water",
                         "Sediment, Natural, Saltwater",
                         "Water, Natural, Saltwater",
                         "Water, Urban, Biogas",
                         "Soil, Natural, Coastal",
                         "Soil, Natural, Temperate heath and scrub",
                         "Soil, Natural, Dunes"))%>%
  filter(!is.na(complex))%>%
  mutate(mfd_hab2 = if_else(complex == "Soil\nNatural\nForests", mfd_hab3, mfd_hab2))%>%
  mutate(mfd_hab2 = if_else(grepl("Beech", mfd_hab2), "Beech", mfd_hab2)) %>%
  mutate(mfd_hab2 = if_else(grepl("Birch", mfd_hab2), "Birch", mfd_hab2)) %>%
  mutate(mfd_hab2 = if_else(grepl("Oak", mfd_hab2), "Oak", mfd_hab2)) %>%
  mutate(mfd_hab2 = if_else(complex == "Sediment\nUrban\nFreshwater", paste0(mfd_hab3), paste0(mfd_hab2)))%>%
  filter(!mfd_hab2 %in% c("Enclosed water, Dried", "Birch", "Pine", "Mire (non-habitat type)"))%>%
  mutate(mfd_hab2 = if_else(is.na(mfd_hab2) & complex == "Soil\nNatural\nForests", paste0(mfd_hab1, " - no MFDO2"), mfd_hab2)) %>%
  mutate(mfd_hab2 = if_else(mfd_hab2=="NA" & complex == "Soil\nNatural\nForests", paste0(mfd_hab1, " - no MFDO2"), mfd_hab2)) %>%
  filter(!is.na(mfd_hab2))%>%
  filter(!mfd_hab2=="NA")%>%
  filter(!mfd_hab2=="Mire")%>%
  mutate(Read_count=if_else(is.na(Read_count), 0, Read_count))%>%
  mutate(Tax_curated=gsub("Root; cNarG_Nxr; Nitrobacter_Nitrococcus; ","", OTU))%>%
  mutate(Tax_curated=gsub("Root; cNarG_Nxr; Nitrobacter_Nitrococcus","Nitrobacter_Nitrococcus_umbrella", Tax_curated))%>%
  mutate(Tax_curated=gsub("Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; ","", Tax_curated))%>%
  mutate(Tax_curated=gsub("Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade","Nitrobacter_Nitrococcus_umbrella", Tax_curated))%>%
  mutate(Tax_curated=gsub("Root; Nitrososphaerales; ","", Tax_curated))%>%
  mutate(Tax_curated=gsub("Root; o_Nitrososphaerales; ","", Tax_curated))


#### Converting OTU to wide format for clustering
all_groups_wide <- d %>%
    pivot_wider(names_from = OTU, values_from = Read_count, id_cols = c(SeqId)) %>%
    as.data.frame(.) %>%
    filter(rowSums(select(., -1)) > 0)
  
  OTU_sums <- all_groups_wide[,-1]
  row.names(OTU_sums) <- all_groups_wide$SeqId
  
  bray_curtis_dist <- vegan::vegdist(vegan::decostand(OTU_sums, method = "hellinger"))
  hclust_ward <- hclust(bray_curtis_dist, method = "ward.D2")
  ward_dendrogram <- as.dendrogram(hclust_ward)
  ward_order <- order.dendrogram(ward_dendrogram)
  level<- hclust_ward$labels[order.dendrogram(ward_dendrogram)]
  

  
### Plotting a heatmap
heat <- d %>%
    mutate(SeqId = factor(SeqId, levels = level, ordered = TRUE)) %>%
    mutate(complex=gsub("Soil\nNatural\nBogs, mires and fens", "Soil\nNatural\nBogs and\nfens", complex))%>%
    mutate(complex=gsub("Soil\nNatural\nGrassland formations", "Soil\nNatural\nGrassland\nformations", complex))%>%
    filter(!complex=="Water\nUrban\nWastewater")%>%
    #arrange(SeqId) %>% 
    ggplot(aes(y = Tax_curated, x = SeqId, fill = Read_count)) +
      geom_tile(width=10.0, linewidth=0.0) +
    scale_fill_viridis_c(name = "Read count\n", trans = "sqrt", limits = c(0, 100), na.value = "#F8E622")+
   # scale_fill_gradient(low = "white", high = "darkred",  na.value = "darkred", name = "RPKM", trans = "sqrt", limits = c(0.00, 3)) +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          #axis.text.y = element_blank(),
          axis.text.y = element_text(size = 35),
          strip.placement = "outside",
          legend.position = "right",
          strip.text = element_text(size = 40, face = "bold"),
          strip.background = element_rect(fill = "grey90", color = "white"),
          legend.key.size = unit(2, 'cm'),
          legend.title = element_text(size=35, hjust=0.9, vjust=0.9),
          #legend.position=c(0.85, 0.3),
          legend.text = element_text(size=35))+
      facet_nested(Group+Package ~ complex, scales = "free", space = "free") +
   scale_y_discrete(expand = c(0,0))
  

heat
  

ggsave("./HM_comparing_packs.png",
       heat,
        height = 10,
  width = 25)

ggsave("./HM_comparing_packs.svg",
       heat,
        height = 10,
  width = 25)


ggsave("./HM_comparing_packs.svg",
       heat,
        height = 10,
  width = 25)

tiff(file = 'Supplemenatry_package_update.tiff',
    height = 1900,
    width = 5500) 
heat
dev.off()


pdf(file = 'Supplemenatry_package_update.pdf',
    height = 30,
    width = 90) 
heat
dev.off()










