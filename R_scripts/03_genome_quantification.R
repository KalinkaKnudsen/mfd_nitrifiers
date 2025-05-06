#title: "03_genome_quantification"
#author: "Kalinka Sand Knudsen"
#update: "2025-04-24"

#Loading packages
.libPaths(c("path_to_R_package_library", .libPaths()))
library(ggplot2)
library(vroom)
library(tidyverse)
library(readxl)
library(ggstats)
library(ggpubr)
library(ggh4x)
library(patchwork)


#Set WD
setwd("path_to_working_directory")

df<-vroom("data/MFD_SRnodrep_tax_relative_abundance.tsv", delim = "\t")%>%
  filter(grepl("t__", clade_name))%>%
  pivot_longer(cols=starts_with("/projects"), names_to="SeqId", values_to = "Taxonomic_abundance", names_prefix="/projects/microflora_danica/data/sequencing/data_flat/trimmed/")%>%
  filter(Taxonomic_abundance>0)



sp<-df%>%
  mutate(SeqId=gsub("_R1.fastq.gz", "", SeqId))%>%
  separate(clade_name, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"), sep = "\\|")%>%
  mutate(genome=gsub("t__", "",
                     gsub(".fa", "", Strain)))


drep<-vroom("data/mags_shallow_all.tsv", delim = "\t")%>%
  filter(bin %in% sp$genome)%>%
  rename(genome=bin)%>%
  select(genome, primary_cluster, secondary_cluster)



s2<-sp%>%
  left_join(drep)%>%
  group_by(secondary_cluster, SeqId)%>%
  summarise(drep_abundance = sum(Taxonomic_abundance))%>%ungroup()



linkage <- vroom("metadata/2023-10-11_samples_minimal_metadata_collapsed.csv", delim = ",") %>%
  mutate(flat_name=gsub(".fastq.gz","", flat_name))%>%
  rename(SeqId = flat_name) %>%
  relocate(SeqId)%>%
  filter(after_total_reads>1000)%>%
  select(SeqId, fieldsample_barcode)
#
metadata.sub <- readxl::read_excel('metadata/2025-02-19_mfd_db.xlsx') %>%
  left_join(linkage)%>%
  filter(SeqId %in% linkage$SeqId)%>%
  relocate(SeqId)%>%
  select(SeqId:mfd_hab3, cell.10km, cell.1km) %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", "))

merge<-sp%>%
  left_join(drep)%>%
  left_join(s2)%>%
  left_join(metadata.sub)


levels_hab2<-readRDS("output/levels_hab2_25_02_21.rds")


top_clades <- merge %>%
  filter(SeqId %in% levels_hab2)%>%
  filter(mfd_hab1 %in% c("Fields", "Greenspaces","Bogs, mires and fens", "Grassland formations", "Forests", "Freshwater")) %>%
  group_by(secondary_cluster) %>%
  summarize(total_abundance = sum(Taxonomic_abundance, na.rm = TRUE)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 50) %>%
  pull(secondary_cluster)

Genera_save<-merge%>%filter(secondary_cluster %in% top_clades)%>%select(Genus)%>%distinct%>%pull(Genus)

saveRDS(top_clades, "output/most_abundant_clades.rds")
saveRDS(Genera_save, "output/most_abundant_Genera.rds")

# 
# top_clades<-readRDS("output/most_abundant_clades.rds")
# Genera_save<-readRDS("output/most_abundant_Genera.rds")
pattern <- paste(Genera_save, collapse = "|")

# Use the combined pattern in grepl

df3<-vroom("data/MFD_SRnodrep_tax_relative_abundance.tsv", delim = "\t")%>%
  filter(grepl("t__", clade_name))%>%
  filter(grepl(pattern, clade_name))

drep<-vroom("data/mags_shallow_all.tsv", delim = "\t")%>%
  filter(secondary_cluster %in% top_clades)%>%
  rename(genome=bin)%>%
  select(genome, primary_cluster, secondary_cluster)


df4 <- df3 %>%
  filter(str_detect(clade_name, paste(drep$genome, collapse = "|")))

df4<-df4%>%
  pivot_longer(cols=starts_with("/projects"), names_to="SeqId", values_to = "Taxonomic_abundance", names_prefix="/projects/microflora_danica/data/sequencing/data_flat/trimmed/")


sp<-df4%>%
  mutate(SeqId=gsub("_R1.fastq.gz", "", SeqId))%>%
  separate(clade_name, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"), sep = "\\|")%>%
  mutate(genome=gsub("t__", "",
                     gsub(".fa", "", Strain)))


drep<-vroom("data/mags_shallow_all.tsv", delim = "\t")%>%
  filter(bin %in% sp$genome)%>%
  rename(genome=bin)%>%
  select(genome, primary_cluster, secondary_cluster)



s2<-sp%>%
  left_join(drep)%>%
  group_by(secondary_cluster, SeqId)%>%
  summarise(drep_abundance = sum(Taxonomic_abundance))%>%ungroup()



linkage <- vroom("metadata/2023-10-11_samples_minimal_metadata_collapsed.csv", delim = ",") %>%
  mutate(flat_name=gsub(".fastq.gz","", flat_name))%>%
  rename(SeqId = flat_name) %>%
  relocate(SeqId)%>%
  filter(after_total_reads>500000)%>%
  select(SeqId, fieldsample_barcode)
#
metadata.sub <- readxl::read_excel('metadata/2025-02-19_mfd_db.xlsx') %>%
  left_join(linkage)%>%
  filter(SeqId %in% linkage$SeqId)%>%
  relocate(SeqId)%>%
  select(SeqId:mfd_hab3, cell.10km, cell.1km) %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", "))

merge<-sp%>%
  left_join(drep)%>%
  left_join(s2)%>%
  left_join(metadata.sub)

saveRDS(merge, "output/sylph_abund_24_12_15.rds")

levels_hab2<-readRDS("output/levels_hab2_25_02_21.rds")

sylph_heat<-merge%>%
  filter(SeqId %in% levels_hab2)%>%
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
  mutate(mfd_hab2=gsub("\\s*\\(exotic\\)", "", mfd_hab2))%>%
  mutate(mfd_hab2=gsub("Standing freshwater", "Standing\nfreshwater", mfd_hab2))%>%
  mutate(mfd_hab2=gsub("Running freshwater", "Running\nfreshwater", mfd_hab2))%>%
  mutate(mfd_hab2=gsub("Activated sludge", "Activated\nsludge", mfd_hab2))%>%
  mutate(hab1_label=gsub("Bogs, mires and fens", "Bogs, mires\nand fens", mfd_hab1),
         hab1_label=gsub("Grassland formations", "Grassland\nformations", hab1_label),
         hab1_label=gsub("Temperate heath and scrub","Heath and\nscrub", hab1_label),
         hab1_label=gsub("Greenspaces","Green-\nspaces", hab1_label),
         hab1_label=gsub("Coastal","Coast", hab1_label),
         hab1_label=gsub("Freshwater", "Freshwater\nsediment", hab1_label))%>%
  filter(!mfd_hab2 %in% c("Spruce", "Willow"))



saveRDS(sylph_heat, "output/AOA_SR_24_12_15.rds")

################### Plotting #########################

sylph_heat<-readRDS("output/AOA_SR_24_12_15.rds")
levels_hab2<-readRDS("output/levels_hab2_25_02_21.rds")

unique(sylph_heat$secondary_cluster)



names <- sylph_heat %>%
  filter(secondary_cluster %in% top_clades)%>%
  select(genome, Species, secondary_cluster) %>% # Select relevant columns
  distinct() %>% # Remove duplicates
  group_by(Species, secondary_cluster) %>% # Group by Species and secondary_cluster
  summarise(genome_count = n_distinct(genome)) # Count unique genomes in each group


species_counts <- sylph_heat %>%
  select(Species, secondary_cluster) %>%
  group_by(secondary_cluster, Species) %>%
  summarise(entry_count = n(), .groups = "drop")

# Step 2: Find the Species with the highest count for each secondary_cluster
top_species_per_cluster <- species_counts %>%
  group_by(secondary_cluster) %>%
  slice_max(entry_count, n = 1, with_ties = FALSE) %>%
  ungroup()%>%
  mutate(Species_label = Species)%>%
  select(secondary_cluster, Species_label)





t<-sylph_heat%>%left_join(top_species_per_cluster)%>%  mutate(sp_label=paste0(Genus, ";", Species_label, " - ", secondary_cluster))


heat <- t %>%
  filter(secondary_cluster %in% top_clades)%>%
  # mutate(sp_label=paste0(Genus, "_", Species, " - ", secondary_cluster))   %>%
  filter(mfd_hab1 %in% c("Fields", "Greenspaces","Bogs, mires and fens", "Grassland formations", "Forests", "Freshwater"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  ggplot(aes(y = sp_label, x = SeqId, fill = `drep_abundance`)) +
  geom_tile() +
  #scale_fill_viridis_c(name = "Taxonomic\nAbundance\n(%)", trans = "sqrt", limits = c(0.00, 2), na.value = "#F8E622")+
  scale_fill_gradientn(
    name = "Taxonomic\nAbundance\n(%)",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  
    trans = "sqrt",  # Same as in viridis
    limits = c(0, 20),  # Set the limits manually
    na.value = "black"  # Handling values larger than 20
  ) +
  facet_nested(Phylum+Family~ factor(hab1_label,levels=c("Fields", "Green-\nspaces", "Grassland\nformations","Forests", "Bogs, mires\nand fens", "Freshwater\nsediment")), scales = "free", space = "free") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.y = element_blank(),
        axis.text.y = element_text(size = 7),
        legend.position = "right",
        legend.margin=margin(-1,0,0,0),
        legend.box.margin=margin(-1,0,0,0),
        legend.title = element_text(size=7),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.3, "cm"),
        strip.text.x = element_text(size = 8, face="bold"),
        strip.text.y = element_text(angle=0, size=8),
        strip.background = element_blank(),
        strip.clip = "off",
        #  strip.text.y = element_text(size = 8.2, angle = 0),  # Size and orientation for y-axis facet labels
        #   strip.background = element_rect(fill = "grey90", color = "grey90"),
        panel.border = element_rect(colour="black", fill=NA) ,
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.03, "cm", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))

ggsave("output/most_abund5.png", heat, width=20, height=11, dpi=300)



p_combine<-readRDS("output/palette_mfd_hab2.rds")
names(p_combine)[names(p_combine) == "Activated sludge"] <- "Activated\nsludge"
names(p_combine)[names(p_combine) == "Non-native trees (exotic)"] <- "Non-native trees"
names(p_combine)[names(p_combine) == "Standing freshwater"] <- "Standing\nfreshwater"
names(p_combine)[names(p_combine) == "Running freshwater"] <- "Running\nfreshwater"
#names(p_combine)[names(p_combine) == "Forests - no MFDO2"] <- "Forests\nno MFDO2"

bar_plot <- t %>%
  filter(mfd_hab1 %in% c("Fields", "Greenspaces","Bogs, mires and fens", "Grassland formations", "Forests"))%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  ggplot(aes(x = SeqId, fill = mfd_hab2)) +
  geom_tile(aes(y = 1)) +
  facet_grid(. ~ factor(hab1_label,levels=c("Fields", "Green-\nspaces", "Grassland\nformations","Forests", "Bogs, mires\nand fens")), scales = "free", space = "free") +
  scale_fill_manual(values = p_combine) + # Specify your colors
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "none",
    panel.spacing = unit(0.1, "lines"),
    panel.border = element_rect(colour="black", fill=NA),
    legend.text = element_text(size=7),
    legend.title = element_blank(),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing = unit(0.1, "cm"),
  )+   scale_y_continuous(expand = c(0,0)) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))



# Step 3: Combine the heatmap plot and the bar plot

combined_plot <- heat / bar_plot + plot_layout(heights = c(12, 0.3)) + theme(panel.background = element_blank(), plot.background = element_blank())
#combined_plot

ggsave("output/Sylph_AOA_SR_24_12_13.png", combined_plot, width=15, height=7, dpi=200)





############################################################## Making this into an average by habitat ###################################################


df<-vroom("data/MFD_SRnodrep_tax_relative_abundance.tsv", delim = "\t")%>%
  filter(grepl("t__", clade_name))%>%
  pivot_longer(cols=starts_with("/projects"), names_to="SeqId", values_to = "Taxonomic_abundance", names_prefix="/projects/microflora_danica/data/sequencing/data_flat/trimmed/")%>%
  filter(Taxonomic_abundance>0)



sp<-df%>%
  mutate(SeqId=gsub("_R1.fastq.gz", "", SeqId))%>%
  separate(clade_name, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"), sep = "\\|")%>%
  mutate(genome=gsub("t__", "",
                     gsub(".fa", "", Strain)))


drep<-vroom("data/mags_shallow_all.tsv", delim = "\t")%>%
  filter(bin %in% sp$genome)%>%
  rename(genome=bin)%>%
  select(genome, primary_cluster, secondary_cluster)



s2<-sp%>%
  left_join(drep)%>%
  group_by(secondary_cluster, SeqId)%>%
  summarise(drep_abundance = sum(Taxonomic_abundance))%>%ungroup()



linkage <- vroom("metadata/2023-10-11_samples_minimal_metadata_collapsed.csv", delim = ",") %>%
  mutate(flat_name=gsub(".fastq.gz","", flat_name))%>%
  rename(SeqId = flat_name) %>%
  relocate(SeqId)%>%
  filter(after_total_reads>1000)%>%
  select(SeqId, fieldsample_barcode)
#
metadata.sub <- readxl::read_excel('metadata/2025-02-19_mfd_db.xlsx') %>%
  left_join(linkage)%>%
  filter(SeqId %in% linkage$SeqId)%>%
  relocate(SeqId)%>%
  select(SeqId:mfd_hab3, cell.10km, cell.1km) %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", "))

merge<-sp%>%
  left_join(drep)%>%
  left_join(s2)%>%
  left_join(metadata.sub)


levels_hab2<-readRDS("output/levels_hab2_25_02_21.rds")


top_clades <- merge %>%
  filter(SeqId %in% levels_hab2) %>%
  filter(mfd_hab1 %in% c("Fields", "Greenspaces", "Bogs, mires and fens", "Grassland formations", "Forests", "Freshwater")) %>%
  group_by(mfd_hab1) %>%
  mutate(sample_count = n_distinct(SeqId)) %>%  # Count unique samples in each habitat
  group_by(secondary_cluster, mfd_hab1) %>%
  summarize(
    total_abundance = sum(Taxonomic_abundance, na.rm = TRUE),
    sample_count = first(sample_count)  # Ensure consistent sample count for each habitat
  ) %>%
  mutate(normalized_abundance = total_abundance / sample_count) %>%  # Normalize by sample count
  ungroup() %>%
  group_by(secondary_cluster) %>%
  summarize(total_normalized_abundance = sum(normalized_abundance, na.rm = TRUE)) %>%  # Sum normalized abundances across habitats
  arrange(desc(total_normalized_abundance)) %>%
  slice_head(n = 50) %>%
  pull(secondary_cluster)


Genera_save<-merge%>%filter(secondary_cluster %in% top_clades)%>%select(Genus)%>%distinct%>%pull(Genus)


top_clades_old<-readRDS("output/most_abundant_clades.rds")
Genera_save_old<-readRDS("output/most_abundant_Genera.rds")


setdiff(Genera_save, Genera_save_old)
setdiff(top_clades, top_clades_old)




pattern <- paste(Genera_save, collapse = "|")

# Use the combined pattern in grepl

df3<-vroom("data/MFD_SRnodrep_tax_relative_abundance.tsv", delim = "\t")%>%
  filter(grepl("t__", clade_name))%>%
  filter(grepl(pattern, clade_name))

drep<-vroom("data/mags_shallow_all.tsv", delim = "\t")%>%
  filter(secondary_cluster %in% top_clades)%>%
  rename(genome=bin)%>%
  select(genome, primary_cluster, secondary_cluster)


df4 <- df3 %>%
  filter(str_detect(clade_name, paste(drep$genome, collapse = "|")))
#
df4<-df4%>%
  pivot_longer(cols=starts_with("/projects"), names_to="SeqId", values_to = "Taxonomic_abundance", names_prefix="/projects/microflora_danica/data/sequencing/data_flat/trimmed/")


sp<-df4%>%
  mutate(SeqId=gsub("_R1.fastq.gz", "", SeqId))%>%
  separate(clade_name, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"), sep = "\\|")%>%
  mutate(genome=gsub("t__", "",
                     gsub(".fa", "", Strain)))


drep<-vroom("data/mags_shallow_all.tsv", delim = "\t")%>%
  filter(bin %in% sp$genome)%>%
  rename(genome=bin)%>%
  select(genome, primary_cluster, secondary_cluster)



s2<-sp%>%
  left_join(drep)%>%
  group_by(secondary_cluster, SeqId)%>%
  summarise(drep_abundance = sum(Taxonomic_abundance))%>%ungroup()



linkage <- vroom("metadata/2023-10-11_samples_minimal_metadata_collapsed.csv", delim = ",") %>%
  mutate(flat_name=gsub(".fastq.gz","", flat_name))%>%
  rename(SeqId = flat_name) %>%
  relocate(SeqId)%>%
  filter(after_total_reads>1000)%>%
  select(SeqId, fieldsample_barcode)
#
metadata.sub <- readxl::read_excel('metadata/2025-02-19_mfd_db.xlsx') %>%
  left_join(linkage)%>%
  filter(SeqId %in% linkage$SeqId)%>%
  relocate(SeqId)%>%
  select(SeqId:mfd_hab3, cell.10km, cell.1km) %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", "))

merge<-sp%>%
  left_join(drep)%>%
  left_join(s2)%>%
  left_join(metadata.sub)




#saveRDS(merge, "output/sylph_abund_24_12_16_v2.rds")

############################################# Data can be loaded from here ########################

merge<-readRDS("output/sylph_abund_24_12_16_v2.rds")

levels_hab2<-readRDS("output/levels_hab2_25_02_21.rds")

top_clades <- merge %>%
  filter(SeqId %in% levels_hab2) %>%
  filter(mfd_hab1 %in% c("Fields", "Greenspaces", "Bogs, mires and fens", "Grassland formations", "Forests", "Freshwater")) %>%
  group_by(mfd_hab1) %>%
  mutate(sample_count = n_distinct(SeqId)) %>%  # Count unique samples in each habitat
  group_by(secondary_cluster, mfd_hab1) %>%
  summarize(
    total_abundance = sum(Taxonomic_abundance, na.rm = TRUE),
    sample_count = first(sample_count)  # Ensure consistent sample count for each habitat
  ) %>%
  mutate(normalized_abundance = total_abundance / sample_count) %>%  # Normalize by sample count
  ungroup() %>%
  group_by(secondary_cluster) %>%
  summarize(total_normalized_abundance = sum(normalized_abundance, na.rm = TRUE)) %>%  # Sum normalized abundances across habitats
  arrange(desc(total_normalized_abundance)) %>%
  slice_head(n = 50) %>%
  arrange(total_normalized_abundance)%>%
  pull(secondary_cluster)


top_clades_df <- merge %>%
  filter(SeqId %in% levels_hab2) %>%
  filter(mfd_hab1 %in% c("Fields", "Greenspaces", "Bogs, mires and fens", "Grassland formations", "Forests", "Freshwater")) %>%
  group_by(mfd_hab1) %>%
  mutate(sample_count = n_distinct(SeqId)) %>%  # Count unique samples in each habitat
  group_by(secondary_cluster, mfd_hab1) %>%
  summarize(
    total_abundance = sum(Taxonomic_abundance, na.rm = TRUE),
    sample_count = first(sample_count)  # Ensure consistent sample count for each habitat
  ) %>%
  mutate(normalized_abundance = total_abundance / sample_count) %>%  # Normalize by sample count
  ungroup() %>%
  group_by(secondary_cluster) %>%
  summarize(total_normalized_abundance = sum(normalized_abundance, na.rm = TRUE),
            total_abundance=first(total_abundance)) %>%  # Sum normalized abundances across habitats
  arrange(desc(total_normalized_abundance)) %>%
  slice_head(n = 50)%>%
  arrange(total_normalized_abundance)


top_clades<-top_clades_df%>%pull(secondary_cluster)
top_clades_total<-top_clades_df%>%arrange(total_abundance)%>%pull(secondary_cluster)


##### This here is for picking the cluster names in case they have been called differently by GTDB-tk
names <- merge %>%
  filter(secondary_cluster %in% top_clades)%>%
  select(genome, Species, secondary_cluster) %>% # Select relevant columns
  distinct() %>% # Remove duplicates
  group_by(Species, secondary_cluster) %>% # Group by Species and secondary_cluster
  summarise(genome_count = n_distinct(genome)) # Count unique genomes in each group


species_counts <- merge %>%
  select(Species, secondary_cluster) %>%
  group_by(secondary_cluster, Species) %>%
  summarise(entry_count = n(), .groups = "drop")

# Find the Species with the highest count for each secondary_cluster
top_species_per_cluster <- species_counts %>%
  group_by(secondary_cluster) %>%
  slice_max(entry_count, n = 1, with_ties = FALSE) %>%
  ungroup()%>%
  mutate(Species_label = Species)%>%
  select(secondary_cluster, Species_label)


t<-merge%>%left_join(top_species_per_cluster)%>%  mutate(sp_label=paste0(Genus, ";", Species_label, " - ", secondary_cluster))%>%
  mutate(sp_label_long=paste0(Phylum, ";", Family, ";", sp_label ))%>%
  select(starts_with("mfd"), SeqId, starts_with("sp_label"), drep_abundance, secondary_cluster)%>%
  distinct()

#order<-t%>%select(sp_label_long)%>%distinct()%>%pull(sp_label_long)


t<-t%>%
  filter(secondary_cluster %in% top_clades)%>%
  filter(SeqId %in% levels_hab2)%>%
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
  mutate(mfd_hab2=gsub("\\s*\\(exotic\\)", "", mfd_hab2))%>%
  mutate(mfd_hab2=gsub("Standing freshwater", "Standing\nfreshwater", mfd_hab2))%>%
  mutate(mfd_hab2=gsub("Running freshwater", "Running\nfreshwater", mfd_hab2))%>%
  mutate(mfd_hab2=gsub("Activated sludge", "Activated\nsludge", mfd_hab2))%>%
  mutate(hab1_label=gsub("Bogs, mires and fens", "Bogs, mires\nand fens", mfd_hab1),
         hab1_label=gsub("Grassland formations", "Grassland\nformations", hab1_label),
         hab1_label=gsub("Temperate heath and scrub","Heath and\nscrub", hab1_label),
         hab1_label=gsub("Greenspaces","Green-\nspaces", hab1_label),
         hab1_label=gsub("Coastal","Coast", hab1_label),
         hab1_label=gsub("Freshwater", "Freshwater\nsediment", hab1_label))%>%
  filter(!mfd_hab2 %in% c("Spruce", "Willow"))


saveRDS(t, "output/sylph_abund_24_12_23.rds")

t<-t%>%mutate(secondary_cluster=factor(secondary_cluster, levels=top_clades, ordered=T))

t <- t[order(t$secondary_cluster),]
t$sp_label_long <- factor(t$sp_label_long, levels = unique(t$sp_label_long))

t_sum<-t%>%
  group_by(sp_label_long)%>%
  summarize(sum_abund=sum(drep_abundance))%>%
  arrange((sum_abund))


### Removing unknowns with double labelling

heat <- t %>%
  filter(!sp_label %in% c("g__;s__ - 4552_1", "g__;s__ - 2307_1", "g__;s__ - 2910_1",
                          "g__;s__ - 2910_3", "g__;s__ - 1684_2"))%>%
  # mutate(sp_label_long=paste0(Phylum, ";", Family, ";", sp_label ))   %>%
  mutate(sp_label_long=factor(sp_label_long, levels = t_sum$sp_label_long, ordered = TRUE))%>%
  filter(mfd_hab1 %in% c("Fields", "Greenspaces","Bogs, mires and fens", "Grassland formations", "Forests", "Freshwater"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  ggplot(aes(y = sp_label_long, x = SeqId, fill = `drep_abundance`)) +
  geom_tile() +
  #scale_fill_viridis_c(name = "Taxonomic\nAbundance\n(%)", trans = "sqrt", limits = c(0.00, 2), na.value = "#F8E622")+
  scale_fill_gradientn(
    name = "Taxonomic\nAbundance\n(%)",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  
    trans = "sqrt",  # Same as in viridis
    limits = c(0, 30),  # Set the limits manually
    na.value = "black"  # Handling values larger than the limit
  ) +
  facet_nested(.~ factor(hab1_label,levels=c("Fields", "Green-\nspaces", "Grassland\nformations","Forests", "Bogs, mires\nand fens", "Freshwater\nsediment")), scales = "free", space = "free") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.y = element_blank(),
        axis.text.y = element_text(size = 7),
        legend.position = "right",
        legend.margin=margin(-1,0,0,0),
        legend.box.margin=margin(-1,0,0,0),
        legend.title = element_text(size=7),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.3, "cm"),
        strip.text.x = element_text(size = 8, face="bold"),
        strip.text.y = element_text(angle=0, size=8),
        strip.background = element_blank(),
        strip.clip = "off",
        #  strip.text.y = element_text(size = 8.2, angle = 0),  # Size and orientation for y-axis facet labels
        #   strip.background = element_rect(fill = "grey90", color = "grey90"),
        panel.border = element_rect(colour="black", fill=NA) ,
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.03, "cm", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))

ggsave("output/most_abund_avg_24_12_23.png", heat, width=18, height=8, dpi=600)


