#title: "06_Comparing_Nitrifier_gene_packages"
#author: "Kalinka Sand Knudsen"
#update: "2025-04-24"


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
cNXR_old<-vroom("data/cNXR_old_combined_count_table_e10.txt", delim="\t")%>%
  rename(OTU = ConsensusLineage)  %>%
  filter(grepl("Root; cNarG_Nxr; Nitrobacter_Nitrococcus", OTU)) %>%
  pivot_longer(starts_with("LIB"), names_to = "SeqId", values_to="Read_count")%>%
  mutate(Package="- MFD")%>%
  mutate(Group="Nitrobacter-\nlike nxrA")

cNXR_new<-vroom("data/cNXR_combined_count_table_e10.txt", delim="\t")%>%
  rename(OTU = ConsensusLineage)  %>%
  filter(grepl("Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade", OTU))%>%
 pivot_longer(starts_with("LIB"), names_to = "SeqId", values_to="Read_count")%>%
  mutate(Package="+ MFD")%>%
  mutate(Group="Nitrobacter-\nlike nxrA")

Arc_old<-vroom("data/arch_AmoA_old_combined_count_table_e10.txt", delim="\t")%>%
  rename(OTU = ConsensusLineage)  %>%
  filter(grepl("Root; Nitrososphaerales; Nitrososphaeraceae", OTU))%>%
  pivot_longer(starts_with("LIB"), names_to = "SeqId", values_to="Read_count")%>%
  mutate(Package="- MFD")%>%
  mutate(Group="Nitrososphaeraceae-\namoA")

Arc_new<-vroom("data/arch_AmoA_combined_count_table_e10.txt", delim="\t")%>%
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
linkage <- vroom("metadata/2023-10-11_samples_minimal_metadata_collapsed.csv", delim = ",") %>%
  mutate(flat_name=gsub(".fastq.gz","", flat_name))%>%
  rename(SeqId = flat_name) %>%
  relocate(SeqId)%>%
  select(SeqId, fieldsample_barcode, after_total_reads)

metadata.sub <- readxl::read_excel('metadata/2025-02-19_mfd_db.xlsx') %>%
  left_join(linkage)%>%
  relocate(SeqId)%>%
  select(SeqId:mfd_hab3, after_total_reads) %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", "))%>%
    filter(after_total_reads>500000)

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
  mutate(mfd_hab2 = if_else(is.na(mfd_hab2) & complex == "Soil\nNatural\nForests", paste0(mfd_hab1, "\nno MFDO2"), mfd_hab2)) %>%
  mutate(mfd_hab2 = if_else(mfd_hab2=="NA" & complex == "Soil\nNatural\nForests", paste0(mfd_hab1, "\nno MFDO2"), mfd_hab2)) %>%
  filter(!is.na(mfd_hab2))%>%
  filter(!mfd_hab2=="NA")%>%
  filter(!mfd_hab2 %in% c("Mire", "Spruce", "Willow"))%>%
  mutate(mfd_hab2=gsub("\\s*\\(exotic\\)", "", mfd_hab2))%>%
  mutate(mfd_hab2=gsub("Standing freshwater", "Standing\nfreshwater", mfd_hab2))%>%
  mutate(mfd_hab2=gsub("Running freshwater", "Running\nfreshwater", mfd_hab2))%>%
  mutate(mfd_hab2=gsub("Activated sludge", "Activated\nsludge", mfd_hab2))%>%
  mutate(Read_count=if_else(is.na(Read_count), 0, Read_count))%>%
  mutate(Tax_curated=gsub("Root; cNarG_Nxr; Nitrobacter_Nitrococcus; ","", OTU))%>%
  mutate(Tax_curated=gsub("Root; cNarG_Nxr; Nitrobacter_Nitrococcus","Nitrobacter_Nitrococcus_umbrella", Tax_curated))%>%
  mutate(Tax_curated=gsub("Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; ","", Tax_curated))%>%
  mutate(Tax_curated=gsub("Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade","Nitrobacter_Nitrococcus_umbrella", Tax_curated))%>%
  mutate(Tax_curated=gsub("Root; Nitrososphaerales; ","", Tax_curated))%>%
  mutate(Tax_curated=gsub("Root; o_Nitrososphaerales; ","", Tax_curated))%>%
  filter(after_total_reads>500000)

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
  

## plotting
 
heat <- d %>%
    mutate(SeqId = factor(SeqId, levels = level, ordered = TRUE)) %>%
    mutate(complex=gsub("Soil\nNatural\nBogs, mires and fens", "Soil\nNatural\nBogs and\nfens", complex))%>%
    mutate(complex=gsub("Soil\nNatural\nGrassland formations", "Soil\nNatural\nGrassland\nformations", complex))%>%
    filter(!complex=="Water\nUrban\nWastewater")%>%
    #arrange(SeqId) %>% 
    ggplot(aes(y = Tax_curated, x = SeqId, fill = Read_count)) +
    geom_tile() +
    scale_fill_gradientn(
    name = "Read\ncount",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 100),  # Set the limits manually
    na.value = "black"  # Handling values larger than the limit
  ) +
    labs(x = "", y = "") +
    theme_minimal() +
   theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        #  axis.text.y = element_blank(),
        axis.text.y = element_text(size = 7),
        legend.key.size = unit(0.38, "cm"),
        legend.title = element_text(size=7.5),
        legend.text = element_text(size=7),
        strip.clip = "off",
        strip.text.x = element_text(size = 7, face="bold"),  # Size for x-axis facet labels
        strip.text.y = element_text(size = 7, face="bold"),
        panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2) ,
        strip.background = element_blank(),
      #  text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.15, "lines", data = NULL))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+
      facet_nested(Group+Package ~ complex, scales = "free", space = "free") +
   scale_y_discrete(expand = c(0,0))
  

heat
  

ggsave("output/Supplemenatry_package_update.png",
       heat,
        height = 8,
  width = 20, dpi=600)

ggsave("output/Supplemenatry_package_update.svg",
       heat,
        height = 8,
  width = 20)


tiff(file = 'output/Supplemenatry_package_update.tiff',
    height = 8,
  width = 20, res=600, units="in")
heat
dev.off()


pdf(file = 'output/Supplemenatry_package_update.pdf',
    height = 8,
    width = 20) 
heat
dev.off()





