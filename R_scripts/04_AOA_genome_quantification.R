#title: "04_AOA_genome_quantification"
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


# Importing the dataframe with the taxonomic relative abundance of taxa
df<-vroom("data/MFD_SRnodrep_tax_relative_abundance.tsv", delim = "\t")%>%
  filter(grepl("f__Nitrososphaeraceae|f__Nitrosopumilaceae", clade_name))%>%
  filter(grepl("t__", clade_name))%>%
  pivot_longer(cols=starts_with("/projects"), names_to="SeqId", values_to = "Taxonomic_abundance", names_prefix="/projects/microflora_danica/data/sequencing/data_flat/trimmed/")


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
metadata.sub <- readxl::read_excel('2025-02-19_mfd_db.xlsx') %>%
  left_join(linkage)%>%
  filter(SeqId %in% linkage$SeqId)%>%
  relocate(SeqId)%>%
  select(SeqId:mfd_hab3, cell.10km, cell.1km) %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", "))

merge<-sp%>%
  left_join(drep)%>%
  left_join(s2)%>%
  left_join(metadata.sub)


#sp<-readRDS("output/sylph_AOA_All_species.rds")


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



#saveRDS(sylph_heat, "output/AOA_SR_24_12_13.rds")
#sylph_heat<-readRDS("output/AOA_SR_24_12_13.rds")

names <- sylph_heat %>%
  select(genome, Species, secondary_cluster) %>% # Select relevant columns
  distinct() %>% # Remove duplicates
  group_by(Species, secondary_cluster) %>% # Group by Species and secondary_cluster
  summarise(genome_count = n_distinct(genome)) # Count unique genomes in each group


species_counts <- sylph_heat %>%
  select(Species, secondary_cluster) %>%
  group_by(secondary_cluster, Species) %>%
  summarise(entry_count = n(), .groups = "drop")

#Find the Species with the highest count for each secondary_cluster
top_species_per_cluster <- species_counts %>%
  group_by(secondary_cluster) %>%
  slice_max(entry_count, n = 1, with_ties = FALSE) %>%
  ungroup()%>%
  mutate(Species_label = Species)%>%
  select(secondary_cluster, Species_label)





t<-sylph_heat%>%left_join(top_species_per_cluster)%>%  mutate(sp_label=paste0(Species_label, " - ", secondary_cluster))

Genus_order<-t%>%select(Family, Genus)%>%distinct()%>%arrange(Family)%>%pull(Genus)


heat <- t %>%
  # filter(bin %in% Nitro_tree$tip.label)%>%
  filter(mfd_hab1 %in% c("Fields", "Greenspaces","Bogs, mires and fens", "Grassland formations", "Forests"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(Genus=factor(Genus, levels = Genus_order, ordered = TRUE)) %>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  ggplot(aes(y = sp_label, x = SeqId, fill = `drep_abundance`)) +
  geom_tile() +
  #scale_fill_viridis_c(name = "Taxonomic\nAbundance\n(%)", trans = "sqrt", limits = c(0.00, 2), na.value = "#F8E622")+
  scale_fill_gradientn(
    name = "Taxonomic\nAbundance\n(%)",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  
    trans = "sqrt",  # Same as in viridis
    limits = c(0, 20),  # Set the limits manually
    na.value = "black"  # Handling values larger than the limit
  ) +
  facet_nested( #factor(Genus, levels=c("g__Nitrosotalea", "g__Nitrosotenuis", "g__Nitrosarchaeum", "g__JACQFM01", "g__VHBM01", "g__Nitrosopolaris", "g__Nitrososphaera", "g__TA-21", "g__TH1177", "g__TH5893", "g__TH5896", "g__JAFAQB01", "g__JARBAU01"))
    Genus           
    ~ factor(hab1_label,levels=c("Fields", "Green-\nspaces", "Grassland\nformations","Forests", "Bogs, mires\nand fens")), scales = "free", space = "free") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.y = element_blank(),
        axis.text.y = element_text(size = 8),
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



############ Investigating the abundances of different species clusters ########
more_than_zero<-t%>%filter(mfd_hab1 %in% c("Fields", "Greenspaces","Bogs, mires and fens", "Grassland formations", "Forests"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  filter(secondary_cluster=="81_1")%>%
  filter(drep_abundance > 0) %>%
  select(SeqId, mfd_hab1, mfd_hab2)%>%
  distinct()%>%
  group_by(mfd_hab1) %>%
  summarise(count = n())

less_than_zero<-t%>%filter(mfd_hab1 %in% c("Fields", "Greenspaces","Bogs, mires and fens", "Grassland formations", "Forests"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  filter(secondary_cluster=="81_1")%>%
  filter(drep_abundance == 0) %>%
  select(SeqId, mfd_hab1, mfd_hab2)%>%
  distinct()%>%
  group_by(mfd_hab2) %>%
  summarise(count = n())


p_combine<-readRDS("output/palette_mfd_hab2.rds")
names(p_combine)[names(p_combine) == "Activated sludge"] <- "Activated\nsludge"
names(p_combine)[names(p_combine) == "Non-native trees (exotic)"] <- "Non-native trees"
names(p_combine)[names(p_combine) == "Standing freshwater"] <- "Standing\nfreshwater"
names(p_combine)[names(p_combine) == "Running freshwater"] <- "Running\nfreshwater"
names(p_combine)[names(p_combine) == "Forests - no MFDO2"] <- "Forests\nno MFDO2"

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

ggsave("Sylph_AOA_SR.png", combined_plot, width=15, height=7, dpi=600)




### Creating df with median and max

max_med<-t%>%group_by(Genus, sp_label, hab1_label)%>%summarise(median_abund=median(drep_abundance),
                                                     max_abund=max(drep_abundance))




#################### Making a model ######################


#sylph_heat<-readRDS("output/AOA_SR_24_12_13.rds")

set.seed(5)

names <- sylph_heat %>%
  select(genome, Species, secondary_cluster) %>% # Select relevant columns
  distinct() %>% # Remove duplicates
  group_by(Species, secondary_cluster) %>% # Group by Species and secondary_cluster
  summarise(genome_count = n_distinct(genome)) # Count unique genomes in each group


species_counts <- sylph_heat %>%
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


t<-sylph_heat%>%left_join(top_species_per_cluster)%>%  mutate(sp_label=paste0(Species_label, " - ", secondary_cluster))


TA_sp_98_1<-t%>%
  filter(secondary_cluster=="98_1")%>%
  select(!c("genome", "Strain", "Taxonomic_abundance", "Species"))%>%
  distinct()


#################### Overview plot ############
TA_sp_98_1<-TA_sp_98_1%>%
  filter(!is.na(mfd_hab2))%>%mutate(mfd_hab3=if_else(mfd_hab1=="Greenspaces", mfd_hab2, mfd_hab3))%>%
  mutate(mfd_hab3=gsub("Oak forest  Old acidophilous woods with Q. robur on sandy plains", "Oak forest old acidophilous woods", mfd_hab3))


### Filtering with at least 10 samples per mfd_hab3 ####
samp_10k <- TA_sp_98_1 %>% 
  filter(!is.na(cell.10km))%>%
  group_by(mfd_hab3, cell.10km)%>%
  slice_sample(n=1)%>%
  select(SeqId, mfd_hab3, cell.10km)

keep<-samp_10k %>% 
  group_by(mfd_hab3)%>%
  filter(!is.na(mfd_hab3))%>%
  summarise(Samples = n())%>%
  filter(Samples>9)


thinned_10k <- TA_sp_98_1 %>% 
  filter(mfd_hab3 %in% keep$mfd_hab3)%>%
  filter(SeqId %in% samp_10k$SeqId)%>%
  group_by(mfd_hab3)%>%
 # slice_sample(n=200)%>%
  arrange(mfd_hab1)%>%
  arrange(mfd_hab2)%>%
  arrange(mfd_hab3)%>%
  distinct()

medians <- thinned_10k %>%
  group_by(mfd_hab3) %>%
  summarise(label = paste0("N=", n(),"\nM=",round(median(drep_abundance), 1)))

levels_mfd_hab3<-thinned_10k %>%
  group_by(mfd_hab3) %>%
  summarise(M=mean(drep_abundance),
            n=n())%>%
  arrange(M)



p_combine<-readRDS("output/palette_mfd_hab2.rds")
names(p_combine)[names(p_combine) == "Activated sludge"] <- "Activated\nsludge"
names(p_combine)[names(p_combine) == "Non-native trees (exotic)"] <- "Non-native trees"
names(p_combine)[names(p_combine) == "Standing freshwater"] <- "Standing\nfreshwater"
names(p_combine)[names(p_combine) == "Running freshwater"] <- "Running\nfreshwater"
names(p_combine)[names(p_combine) == "Forests - no MFDO2"] <- "Forests\nno MFDO2"


thinned_10k$mfd_hab2 <- factor(thinned_10k$mfd_hab2, levels = c(
  "Poales, Cereal", "Mixed crops", "Malvids", "Asterids", "Superasterids", "Fabids", 
  "Poales, grass", "Fallow", "Parks", "Other", 
  "Semi-natural tall-herb humid meadows", "Semi-natural dry grasslands", "Natural grasslands", 
  "Alluvial woodland", "Bog woodland", "Non-native trees", "Beech", "Oak", "Coniferous forest", 
  "Deciduous trees", "Calcareous fens", "Sphagnum acid bogs", 
  "Standing\nfreshwater", "Enclosed water", 
  "Rainwater basin, City", "Rainwater basin, Roadside", "Rainwater basin, Dried"
))


# Arrange the data frame based on the factor levels


thinned_10k<-thinned_10k%>%
  left_join(levels_mfd_hab3)%>%
  mutate(label=paste0(mfd_hab3, ", n= ", n))


fields <- thinned_10k%>%
  filter(mfd_hab1 %in% c("Fields"))%>%
  arrange(M)%>%
  arrange(mfd_hab2)


##### Statistics #####

library(rstatix)

plot_test <- thinned_10k%>%
  filter(mfd_hab1 %in% c("Fields"))%>%
  mutate(label = factor(label, levels = unique(fields$label), ordered = TRUE)) %>%
  ggplot(., aes(x = label, y = drep_abundance)) + 
  geom_jitter(aes(fill=mfd_hab2), size = 1.3, alpha = 0.8, color = "black", pch = 21,
              position = position_jitter(width = 0.2, height = 0.03) # Adjust jitter spread if needed
  ) +
  geom_violin(aes(fill=mfd_hab2),  # Add median quantile line
              color = "black",     # Median line color
              size = 0.6,
              alpha=0.4
              #fill="transparent"# Median line thickness
  ) +
  scale_fill_manual(values = p_combine)+
  stat_summary(fun = mean, geom = "crossbar", width = 0.8, size = 0.6, color = "red3") +
  # stat_compare_means(method = "wilcox.test", aes(label=..p.adj..),
  #                    ref.group = ".all.", method.args = list(alternative = "less"))+ ###Adjust here for one sided "less" 
  ggpubr::geom_pwc(method = "wilcox_test",
                   label = 'p.adj.signif',
                   ref.group = "all",
                   hide.ns = T,
                   #angle = 90, vjust = 0.75,
                   p.adjust.method = "fdr",
                   group.by = "x.var",
                   p.adjust.by = "panel",
                   method.args = list(alternative = "greater"),
                  # method.args = list(alternative = "greater"),
                   remove.bracket = T)+
  
   geom_hline(yintercept = mean(fields%>%pull(drep_abundance)), linetype = 2)+ # Add horizontal line at base mean
  scale_y_sqrt() +
  #  facet_grid(.~mfd_hab2, scales="free_x", space = "free")+
  theme_minimal()+
  theme(legend.position="bottom", panel.grid=element_blank(), axis.text.x = element_text(angle=45, hjust = 1))+
  xlab("MFD ontology level 3")+
  ylab("Taxonomic abundance\ns__TA-21 sp022545895")

#label = "p.signif", #aes(label=..p.adj..),
#p.adjust.methods = "fdr",

thinned_10k%>%
  ungroup()%>%
  filter(mfd_hab1 %in% c("Fields"))%>%
  #group_by(mfd_hab3)%>%
  rstatix::wilcox_test(drep_abundance~label,alternative = "greater", p.adjust.method = "fdr", ref.group = "all", detailed=T)%>%print(n = 24)


# ggsave("Scripts_review_1/mfd_hab3_gradient_fields_two_sided.png", plot_test, height=6, width=15, dpi=400)
ggsave("mfd_hab3_fields_one_sided_less_fdr.png", plot_test, height=5, width=10, dpi=400)
ggsave("mfd_hab3_fields_one_sided_less_fdr.svg", plot_test, height=5, width=10, dpi=400)




######################### Other habitats #########################


non_fields <- thinned_10k%>%
  filter(mfd_hab1 %in% c("Forests", "Greenspaces","Grassland formations"))%>%
  arrange(M)%>%
  arrange(mfd_hab1)


plot_test2 <- thinned_10k%>%
  filter(mfd_hab1 %in% c("Forests", "Greenspaces","Grassland formations"))%>%
  mutate(label = factor(label, levels = unique(non_fields$label), ordered = TRUE)) %>%
  ggplot(., aes(x = label, y = drep_abundance)) + 
  geom_jitter(aes(fill=mfd_hab1), size = 1.3, alpha = 0.8, color = "black", pch = 21,
              position = position_jitter(width = 0.2, height = 0.03) # Adjust jitter spread if needed
  ) +
  geom_violin(aes(fill=mfd_hab1),  # Add median quantile line
              color = "black",     # Median line color
              size = 0.6,
              alpha=0.4
              #fill="transparent"# Median line thickness
  ) +
  #  scale_fill_manual(values = p_combine)+
  stat_summary(fun = mean, geom = "crossbar", width = 0.8, size = 0.6, color = "red3") +
  # stat_compare_means(label = "..p.format..", method = "wilcox.test",
  #                    ref.group = ".all.", method.args = list(alternative = "greater"))+ ###Adjust here for one sided "less"

  ggpubr::geom_pwc(method = "wilcox_test",
                   label = 'p.adj.signif',
                   ref.group = ".all.",
                   hide.ns = T,
                   #angle = 90, vjust = 0.75,
                   p.adjust.method = "fdr",
                   group.by = "x.var",
                   p.adjust.by = "panel",
                   method.args = list(alternative = "less"),
                   # method.args = list(alternative = "greater"),
                   remove.bracket = T)+
   geom_hline(yintercept = mean(non_fields%>%pull(drep_abundance)), linetype = 2)+ 
  scale_y_sqrt() +
  #  facet_grid(.~mfd_hab2, scales="free_x", space = "free")+
  theme_minimal()+
  theme(legend.position="bottom", panel.grid=element_blank(), 
        axis.text.x = element_text(angle=45, hjust = 1))+
  xlab("MFD ontology level 3")+
  ylab("Taxonomic abundance\ns__TA-21 sp022545895")


plot_test2

thinned_10k%>%
  ungroup()%>%
  filter(mfd_hab1 %in% c("Forests", "Greenspaces","Grassland formations"))%>%
  #group_by(mfd_hab3)%>%
  rstatix::wilcox_test(drep_abundance~label,alternative = "less", p.adjust.method = "fdr", ref.group = "all", detailed=T)%>%print(n = 24)


#

ggsave("mfd_hab3_natural_one_sided_greater_fdr.png", plot_test2, height=9, width=12, dpi=400)
ggsave("mfd_hab3_natural_one_sided_greater_fdr.svg", plot_test2, height=9, width=12, dpi=400)


ggsave("mfd_hab3_natural_one_sided_greater_col1.png", plot_test2, height=5, width=9, dpi=400)
ggsave("mfd_hab3_natural_one_sided_greater_col1.svg", plot_test2, height=5, width=9, dpi=400)


#




combined<- plot_test + plot_test2 + plot_layout(widths = c(1.2, 1))
#combined

ggsave("mfd_hab3_combined_fdr.png", combined, height=7.5, width=15, dpi=400)
ggsave("mfd_hab3_combined_fdr.svg", combined, height=7.5, width=15, dpi=400)

#





##############################################################################
##############################################################################
##############################################################################




################### GlobDB + SR_MAGs MFD comparison #########################


 globAOA<-vroom("data/SR_AOA_glob_tax_relative_abundance.tsv", delim = "\t")%>%
   filter(grepl("f__Nitrososphaeraceae", clade_name))%>%
   filter(grepl("t__", clade_name))%>%
   pivot_longer(cols=starts_with("/projects"), names_to="SeqId", values_to = "Taxonomic_abundance", names_prefix="/projects/microflora_danica/data/sequencing/data_flat/trimmed/")
 
 globAOA<-globAOA%>%
   mutate(SeqId=gsub("_R1.fastq.gz", "", SeqId))%>%
   separate(clade_name, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"), sep = "\\|")


saveRDS(globAOA, "output/sylph_SR_AOA_glob.rds")

globAOA<-readRDS("output/sylph_SR_AOA_glob.rds")
globAOA<-globAOA%>%
  mutate(genome=gsub(".fa", "", gsub("t__", "", Strain)))

drep<-vroom("data/AOA_drep.csv", delim = ",")%>%
  select(genome, secondary_cluster)%>%
  mutate(genome=gsub(".fa", "", genome))%>%
  filter(genome %in% globAOA$genome)

globAOA2<-globAOA%>%
  left_join(drep)%>%
  group_by(secondary_cluster, SeqId)%>%
  summarise(drep_abundance = sum(Taxonomic_abundance))%>%ungroup()

merge<-drep%>%
  left_join(globAOA)%>%
  left_join(globAOA2)


linkage <- vroom("metadata/2023-10-11_samples_minimal_metadata_collapsed.csv", delim = ",") %>%
  mutate(flat_name=gsub(".fastq.gz","", flat_name))%>%
  rename(SeqId = flat_name) %>%
  relocate(SeqId)%>%
  filter(after_total_reads>1000)%>%
  select(SeqId, fieldsample_barcode)

metadata.sub <- readxl::read_excel('2025-02-19_mfd_db.xlsx') %>%
  left_join(linkage)%>%
  filter(SeqId %in% linkage$SeqId)%>%
  relocate(SeqId)%>%
  select(SeqId:mfd_hab3, cell.10km, cell.1km) %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", "))

sp<-merge%>%left_join(metadata.sub)%>%
  mutate(sp_label=paste0(Genus, " - ", secondary_cluster))


#saveRDS(sp, "output/sylph_AOA_All_species_24_12_20.rds")

#####

#sp<-readRDS("output/sylph_AOA_All_species_24_12_20.rds")

levels_hab2<-readRDS("output/levels_hab2_25_02_21.rds")

sylph_heat<-sp%>%
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


names <- sylph_heat %>%
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



t<-sylph_heat%>%left_join(top_species_per_cluster)%>%  mutate(sp_label=paste0(Species_label, " - ", secondary_cluster))%>%
  select(!c("genome", "Strain", "Taxonomic_abundance"))%>%
  distinct()%>%
  mutate(method="Glob+AOA")


saveRDS(t, "output/Glob_MFD_AOA.rds")



t<-readRDS("output/Glob_MFD_AOA.rds")
levels_hab2<-readRDS("output/levels_hab2_25_02_21.rds")



Genus_order<-t%>%select(Family, Genus)%>%distinct()%>%arrange(Family)%>%pull(Genus)

heat <- t %>%
  # filter(bin %in% Nitro_tree$tip.label)%>%
  filter(mfd_hab1 %in% c("Fields", "Greenspaces","Bogs, mires and fens", "Grassland formations", "Forests"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(Genus=factor(Genus, levels = Genus_order, ordered = TRUE)) %>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  ggplot(aes(y = sp_label, x = SeqId, fill = `drep_abundance`)) +
  geom_tile() +
  #scale_fill_viridis_c(name = "Taxonomic\nAbundance\n(%)", trans = "sqrt", limits = c(0.00, 2), na.value = "#F8E622")+
  scale_fill_gradientn(
    name = "Taxonomic\nAbundance\n(%)",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  
    trans = "sqrt",  # Same as in viridis
    limits = c(0, 100),  # Set the limits manually
    na.value = "black"  # Handling values larger than the limit
  ) +
  facet_nested( #factor(Genus, levels=c("g__Nitrosotalea", "g__Nitrosotenuis", "g__Nitrosarchaeum", "g__JACQFM01", "g__VHBM01", "g__Nitrosopolaris", "g__Nitrososphaera", "g__TA-21", "g__TH1177", "g__TH5893", "g__TH5896", "g__JAFAQB01", "g__JARBAU01"))
    Genus           
    ~ factor(hab1_label,levels=c("Fields", "Green-\nspaces", "Grassland\nformations","Forests", "Bogs, mires\nand fens")), scales = "free", space = "free") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.y = element_blank(),
        axis.text.y = element_text(size = 8),
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


p_combine<-readRDS("output/palette_mfd_hab2.rds")
names(p_combine)[names(p_combine) == "Activated sludge"] <- "Activated\nsludge"
names(p_combine)[names(p_combine) == "Non-native trees (exotic)"] <- "Non-native trees"
names(p_combine)[names(p_combine) == "Standing freshwater"] <- "Standing\nfreshwater"
names(p_combine)[names(p_combine) == "Running freshwater"] <- "Running\nfreshwater"
#names(p_combine)[names(p_combine) == "Forests - no MFDO2"] <- "Forests\nno MFDO2"

setdiff(names(p_combine), unique(t$mfd_hab2))

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

combined_plot <- heat / bar_plot + plot_layout(heights = c(12, 0.1)) + theme(panel.background = element_blank(), plot.background = element_blank())
#combined_plot

ggsave("output/Sylph_Glob_AOA_SR_24_12_20.png", combined_plot, width=15, height=18, dpi=400)




###############################################################################################
####################### Doing the same now for MFD MAGs exclusively ########################


sylph_heat<-readRDS("output/AOA_SR_24_12_13.rds")

sylph_heat<-sylph_heat%>%
  select(!drep_abundance)%>%
  select(!c("secondary_cluster", "primary_cluster"))

drep<-vroom("data/AOA_drep.csv", delim = ",")%>%
  select(genome, secondary_cluster)%>%
  mutate(genome=gsub(".fa", "", genome))%>%
  filter(genome %in% sylph_heat$genome)

sylph_heat2<-sylph_heat%>%
  left_join(drep)%>%
  group_by(secondary_cluster, SeqId)%>%
  summarise(drep_abundance = sum(Taxonomic_abundance))%>%ungroup()

merge<-drep%>%
  left_join(sylph_heat)%>%left_join(sylph_heat2)



names <- merge %>%
  select(genome, Species, secondary_cluster) %>% # Select relevant columns
  distinct() %>% # Remove duplicates
  group_by(Species, secondary_cluster) %>% # Group by Species and secondary_cluster
  summarise(genome_count = n_distinct(genome)) # Count unique genomes in each group


species_counts <- merge %>%
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


t<-readRDS("output/Glob_MFD_AOA.rds")

# t<-merge%>%left_join(top_species_per_cluster)%>%  mutate(sp_label=paste0(Species_label, " - ", secondary_cluster))%>%
#   select(!c("genome", "Strain", "Taxonomic_abundance", "Species"))%>%
#   distinct()

df <- t %>%
  group_by(SeqId) %>%
  mutate(drep_abundance = (drep_abundance / sum(drep_abundance)) * 100) %>%
  ungroup()%>%
  mutate(drep_abundance=as.numeric(gsub("NaN", 0, drep_abundance)))%>%
  mutate(method="MFD_only")

saveRDS(df, "output/MFD_only_AOA.rds")

MFD_only<-readRDS("output/MFD_only_AOA.rds")

t2<-t%>%select(names(MFD_only))

combined<-rbind(MFD_only, t2)%>%
  filter(SeqId %in% t2$SeqId)


Genus_order<-combined%>%select(Family, Genus)%>%distinct()%>%arrange(Family)%>%pull(Genus)
levels_hab2<-readRDS("output/levels_hab2_25_02_21.rds")

remove<-combined%>% filter(mfd_hab1 %in% c("Fields", "Greenspaces","Bogs, mires and fens", "Grassland formations", "Forests"))%>% 
  group_by(secondary_cluster)%>%summarize(row_sum=sum(drep_abundance))%>%ungroup()%>%
  filter(row_sum<30)

combined <- combined %>%
  mutate(method=gsub("MFD_only", "MFD only", method))%>%
  mutate(method=gsub("Glob\\+AOA", "Glob and MFD", method))

##########################################################
#saveRDS(combined, "output/comparison_25_01_23.rds")

combined<-readRDS("output/comparison_25_01_23.rds")
#########################################################

Genus_order<-combined%>%select(Family, Genus)%>%distinct()%>%arrange(Family)%>%pull(Genus)
levels_hab2<-readRDS("output/levels_hab2_25_02_21.rds")

remove<-combined%>% filter(mfd_hab1 %in% c("Fields", "Greenspaces","Bogs, mires and fens", "Grassland formations", "Forests"))%>% 
  group_by(secondary_cluster)%>%summarize(row_sum=sum(drep_abundance))%>%ungroup()%>%
  filter(row_sum<30)





####### Make the MFD only rows that are missing into NA ########

genomes1<-combined%>%
  filter(method=="MFD only")%>%
  select(sp_label)%>%
  distinct()

genomes2<-combined%>%
  filter(method=="Glob and MFD")%>%
  select(sp_label)%>%
  distinct()

genomes_glob<- setdiff(genomes2$sp_label, genomes1$sp_label)

add<-combined%>%
  filter(sp_label %in% genomes_glob)%>%
  mutate(drep_abundance=NA)%>%
  mutate(method="MFD only")%>%
  distinct()


# Merge the combinations with the ad dataframe, keeping all rows from ad
combined_add <- rbind(add, combined)


heat <- combined_add %>%
  filter(!Genus=="g__NA")%>%
  # filter(bin %in% Nitro_tree$tip.label)%>%
  filter(!secondary_cluster %in% remove$secondary_cluster)%>%
  filter(mfd_hab1 %in% c("Fields", "Greenspaces","Bogs, mires and fens", "Grassland formations", "Forests"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(Genus=factor(Genus, levels = Genus_order, ordered = TRUE)) %>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  ggplot(aes(y = sp_label, x = SeqId, fill = `drep_abundance`)) +
  geom_tile() +
  #scale_fill_viridis_c(name = "Taxonomic\nAbundance\n(%)", trans = "sqrt", limits = c(0.00, 2), na.value = "#F8E622")+
  scale_fill_gradientn(
    name = "Taxonomic\nAbundance\n(%)",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  
    trans = "sqrt",  # Same as in viridis
    limits = c(0, 100),  # Set the limits manually
    na.value = "grey85"  # Handling NA values
  ) +
  facet_nested( Genus ~ method+factor(hab1_label,
    levels=c("Fields", "Green-\nspaces", "Grassland\nformations","Forests", "Bogs, mires\nand fens")),
    scales = "free", space = "free", strip = strip_nested(clip = "off"), 
    nest_line = element_line(color="grey20", linewidth = 0.6), resect=unit(4, "pt")) +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.y = element_blank(),
        axis.text.y = element_text(size = 9),
        legend.position = "right",
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        strip.text.x = element_text(size = 9, face="bold"),
        strip.text.y = element_text(angle=0, size=9),
        strip.background = element_blank(),
        #strip.clip = "off",
        #  strip.text.y = element_text(size = 8.2, angle = 0),  # Size and orientation for y-axis facet labels
        #   strip.background = element_rect(fill = "grey90", color = "grey90"),
        panel.border = element_rect(colour="black", fill=NA) ,
      #  text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.06, "cm", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))


ggsave("output/Sylph_AOA_SR_GLOB_combined_25_01_24.png", 
       heat, width=18, height=13, dpi=400)

ggsave("output/Sylph_AOA_SR_GLOB_combined_25_01_24.svg",
       heat,
       width=18, height=13, dpi=400)


tiff(file = 'output/Sylph_AOA_SR_GLOB_combined_25_01_24.tiff',
     height = 13,
     width = 18, res=400, units="in")
heat
dev.off()


pdf(file = './Scripts_review_1/Sylph_testing/Sylph_AOA_SR_GLOB_combined_25_01_24.pdf',
    height = 13,
    width = 18) 
heat
dev.off()




p_combine<-readRDS("output/palette_mfd_hab2.rds")
names(p_combine)[names(p_combine) == "Activated sludge"] <- "Activated\nsludge"
names(p_combine)[names(p_combine) == "Non-native trees (exotic)"] <- "Non-native trees"
names(p_combine)[names(p_combine) == "Standing freshwater"] <- "Standing\nfreshwater"
names(p_combine)[names(p_combine) == "Running freshwater"] <- "Running\nfreshwater"
#names(p_combine)[names(p_combine) == "Forests - no MFDO2"] <- "Forests\nno MFDO2"

setdiff(names(p_combine), unique(t$mfd_hab2))

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

combined_plot <- heat / bar_plot + plot_layout(heights = c(12, 0.1)) + theme(panel.background = element_blank(), plot.background = element_blank())
#combined_plot

ggsave("output/Sylph_AOA_SR_24_12_20.png", combined_plot, width=15, height=12, dpi=100)

