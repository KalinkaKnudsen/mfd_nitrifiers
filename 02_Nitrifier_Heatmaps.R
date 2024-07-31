#title: "02_Nitrifier_Heatmaps"
#author: "Kalinka Sand Knudsen"
#update: "2024-07-30"


#Loading packages
.libPaths(c("path_to_R_package_library", .libPaths()))
library(ggplot2)
library(vroom)
library(tidyverse)
library(readxl)

#Set WD
setwd("path_to_working_directory")

#Load dataframe generated with "normalising_and_aggregating_nitrifier_genes"
OTU_filtered_long<-readRDS("OTU_filtered_long.rds")

#Ordering the tax:
tax<-OTU_filtered_long%>%
  select(Tax)%>%

### Curating the taxonomy
tax_curated<-tax%>%
  mutate(Tax_curated = if_else(Tax=="Root; amoA; Gammaproteobacteria_amoA", "Nitrosomonadaceae_cluster", Tax),
         Tax_curated = if_else(Tax=="Root; amoA; Gammaproteobacteria_amoA; Nitrosospira", "Nitrosomonadaceae; Nitrosospira", Tax_curated),
         Tax_curated = if_else(Tax=="Root; amoA; Nitrospira; Nitrospira_clade_A", "Nitrospira_amoA; Nitrospira_clade_A", Tax_curated),
         Tax_curated = if_else(Tax=="Root; amoA; Nitrospira; Nitrospira_clade_B", "Nitrospira_amoA; Nitrospira_clade_B", Tax_curated),
         Tax_curated = if_else(Tax=="Root; amoA; Nitrospira", "Nitrospira_amoA_cluster", Tax_curated),
         Tax_curated = if_else(Tax=="Root; o_Nitrososphaerales; f_Nitrosopumilaceae; Nitrosotalea_TA20_cluster; g_Nitrosotalea", "Nitrosopumilaceae; Nitrosotalea", Tax_curated),
         Tax_curated = if_else(Tax=="Root; o_Nitrososphaerales; f_Nitrososphaeraceae", "Nitrosospheraceae_cluster", Tax_curated),
         Tax_curated = if_else(Tax=="Root; o_Nitrososphaerales; f_Nitrososphaeraceae; TH5893_Nitrososphaera_cluster", "Nitrosospheraceae; TH5893_Nitrososphaera_cluster", Tax_curated),
         Tax_curated = if_else(Tax=="Root; o_Nitrososphaerales; Other Nitrosopumilaceae", "Other Nitrosopumilaceae", Tax_curated),
         Tax_curated = if_else(Tax=="Root; o_Nitrososphaerales; f_Nitrososphaeraceae; TA21_Nitrosopolaris_cluster; g_Nitrosopolaris", "Nitrosospheraceae; Nitrosopolaris", Tax_curated),
         Tax_curated = if_else(Tax=="Root; o_Nitrososphaerales; f_Nitrososphaeraceae; TA21_Nitrosopolaris_cluster", "Nitrosospheraceae; TA21_Nitrosopolaris_cluster", Tax_curated),
         Tax_curated = if_else(Tax=="Root; o_Nitrososphaerales; f_Nitrososphaeraceae; g_Nitrosocosmicus", "Nitrosospheraceae; Nitrosocosmicus", Tax_curated),
         Tax_curated = if_else(Tax=="Root; o_Nitrososphaerales", "Nitrosospherales_cluster", Tax_curated),
         Tax_curated = if_else(Tax=="Root; o_Nitrososphaerales; f_Nitrososphaeraceae; TA21_Nitrosopolaris_cluster; g_JAFAQB01", "Nitrosospheraceae; JAFAQB01", Tax_curated),
         Tax_curated = if_else(Tax=="Root; o_Nitrososphaerales; f_Nitrososphaeraceae; TH5893_Nitrososphaera_cluster; JAJPHM01_Nitrososphaera_cluster; g_Nitrososphaera", "Nitrosospheraceae; Nitrososphaera", Tax_curated),
         Tax_curated = if_else(Tax=="Root; o_Nitrososphaerales; f_Nitrososphaeraceae; g_TH1177", "Nitrosospheraceae; TH1177", Tax_curated),
         Tax_curated = if_else(Tax=="Root; o_Nitrososphaerales; f_Nitrosopumilaceae; Nitrosotalea_TA20_cluster; g_TA20", "Nitrosopumilaceae; TA20", Tax_curated),
         Tax_curated = if_else(Tax=="Root; Root_narG; Thiocapsa", "Thiocapsa", Tax_curated),
         Tax_curated = if_else(Tax=="Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; Nitrobacter-like_NxrA", "Nitrobacter_like_nxrA_cluster", Tax_curated),
         Tax_curated = if_else(Tax=="Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; Nitrobacter-like_NxrA; Pseudolabrys_JAFA_VAZQ_umbrella", "Nitrobacter_like_nxrA; Put. NOB cluster", Tax_curated),
         Tax_curated = if_else(Tax=="Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; Nitrobacter-like_NxrA; Pseudolabrys_JAFA_VAZQ_umbrella; JAFAXD01_VAZQ01_Pseudolabrys_cluster", "Nitrobacter_like_nxrA; JAF. VAZ. Pseudo. cluster", Tax_curated),
         Tax_curated = if_else(Tax=="Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; Nitrobacter-like_NxrA; Bradyrhizobium_Nitrobacter", "Nitrobacter_like_nxrA; Bradyrh._Nitrobacter_cluster", Tax_curated),
         Tax_curated = if_else(Tax=="Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; Nitrobacter-like_NxrA; Bradyrhizobium_Nitrobacter; Nitrobacter", "Nitrobacter_like_nxrA; Nitrobacter", Tax_curated),
         Tax_curated = if_else(Tax=="Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; Nitrobacter-like_NxrA; Pseudolabrys_JAFA_VAZQ_umbrella; JAFAXD01_VAZQ01_Pseudolabrys_cluster; Pseudolabrys", "Nitrobacter_like_nxrA; Pseudolabrys", Tax_curated),
         Tax_curated = if_else(Tax=="Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; Nitrobacter-like_NxrA; Pseudolabrys_JAFA_VAZQ_umbrella; JAFAXD01_VAZQ01_Pseudolabrys_cluster; VAZQ01", "Nitrobacter_like_nxrA; VAZQ01", Tax_curated),
         Tax_curated = if_else(Tax=="Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade", "Nitrococcus_Nitrobacter_clade", Tax_curated),
         Tax_curated = if_else(Tax=="Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; Nitrococcus", "Nitrococcus", Tax_curated),
         Tax_curated = if_else(Tax=="Root; cNarG_Nxr; Methylomirabilis_acidimicrobia", "Methylomirabilis_Acidimicrobiia_cluster", Tax_curated),
         Tax_curated = if_else(Tax=="Root; Nitrotoga", "Nitrotoga", Tax_curated),
         Tax_curated = if_else(Tax=="Root; periplasmic_umbrella; Nitrospira_umbrella", "Nitrospirales_nxrA_cluster", Tax_curated),
         Tax_curated = if_else(Tax=="Root; periplasmic_umbrella; Nitrospira_umbrella; Nitrospira_clade_1", "Nitrospira_clade_1", Tax_curated),
         Tax_curated = if_else(Tax=="Root; periplasmic_umbrella; Nitrospira_umbrella; Nitrospira_clade_2", "Nitrospira_clade_2", Tax_curated),
         Tax_curated = if_else(Tax=="Root; periplasmic_umbrella; Nitrospira_umbrella; Nitrospira_clade_2; Nitrospira_D", "Nitrospirales; Nitrospira_D", Tax_curated),
         Tax_curated = if_else(Tax=="Root; periplasmic_umbrella; NS_4_NS_12", "NS-4_NS-12", Tax_curated),
         Tax_curated = if_else(Tax=="Root; periplasmic_umbrella; UBA8639", "UBA8639", Tax_curated),
         Tax_curated = if_else(Tax=="Root; amoA; Gammaproteobacteria_amoA; Nitrosomonas", "Nitrosomonadaceae; Nitrosomonas*", Tax_curated),
         Tax_curated = if_else(Tax=="Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; Nitrobacter-like_NxrA; BOG-931", "Nitrobacter_like_nxrA; BOG-931*", Tax_curated),
         Tax_curated = if_else(Tax=="Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; Nitrobacter-like_NxrA; Bradyrhizobium_Nitrobacter; Bradyrhizobium", "Nitrobacter_like_nxrA; Bradyrhizobium*", Tax_curated),
         Tax_curated = if_else(Tax=="Root; o_Nitrososphaerales; f_Nitrosopumilaceae_clade", "Nitrosopumilaceae_clade", Tax_curated),
         Tax_curated = if_else(Tax=="Root; o_Nitrososphaerales; f_Nitrososphaeraceae; TH5893_Nitrososphaera_cluster; g_TH5893", "Nitrosospheraceae; TH5893*", Tax_curated),
         Tax_curated = if_else(Tax=="Root; o_Nitrososphaerales; f_Nitrososphaeraceae; TH5896", "Nitrosospheraceae; TH5896*", Tax_curated),
         Tax_curated = if_else(Tax=="Root; o_Nitrososphaerales; f_Nitrososphaeraceae; g_TA21", "Nitrosospheraceae; TA21*", Tax_curated),
         Tax_curated = if_else(Tax=="Root; periplasmic_umbrella; Nitrospira_umbrella; Nitrospira_C", "Nitrospirales; Nitrospira_C*", Tax_curated),
         Tax_curated = if_else(Tax=="Root; periplasmic_umbrella; Nitrospira_umbrella; Nitrospira_A", "Nitrospirales; Nitrospira_A*", Tax_curated),
         Tax_curated = if_else(Tax=="Root; periplasmic_umbrella; Nitrospira_umbrella; Nitrospira_F", "Nitrospirales; Nitrospira_F*", Tax_curated),
         Tax_curated = if_else(Tax=="Root; periplasmic_umbrella; Nitrospira_umbrella; Palsa-1315", "Nitrospirales; Palsa_1315*", Tax_curated)
         )


#Merging the curated tax onto the dataframe:
OTU_filtered_long<-merge(OTU_filtered_long, tax_curated, by= "Tax")

#Loading a file with ordering of taxonomy as is present in the gene-phylognetic trees
tax2 <- readxl::read_excel("tax_order.xlsx")

# Filtering of metadata and generation of labels
OTU_filtered_long<-OTU_filtered_long %>%
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
  filter(!mfd_hab2=="Mire")%>%
  filter(!is.na(Tax_curated))%>%
  mutate(mfd_hab2=gsub("\\s*\\(exotic\\)", "", mfd_hab2))%>%
  mutate(mfd_hab2=gsub("Standing freshwater", "Standing\nfreshwater", mfd_hab2))%>%
  mutate(mfd_hab2=gsub("Running freshwater", "Running\nfreshwater", mfd_hab2))%>%
  mutate(mfd_hab2=gsub("Activated sludge", "Activated\nsludge", mfd_hab2))

# Making a wide dataframe to be used in clustering
 OTU_wide<-OTU_filtered_long%>%
  pivot_wider(names_from = Tax_curated, values_from = RPKM, id_cols = c(SeqId))%>%
  as.data.frame(.)%>%
  filter(rowSums(select(., -1)) > 0)
OTU_sums<-OTU_wide[,-1]
row.names(OTU_sums)<-OTU_wide$SeqId
saveRDS(file="OTU_wide.rds", OTU_wide) ## To use in PCoA
groups<-unique(OTU_filtered_long$complex_long)


# Clustering of MFD_habo2 individually
groups2<-unique(OTU_filtered_long$mfd_hab2)


for (i in seq_along(groups2)) {
  # Use paste0 to dynamically create variable names gr_1, gr_2, etc.
  gr_name <- paste0("gr_", i)
  
  # Use the filter condition to select the current group
  current_group <- OTU_filtered_long %>%
    filter(mfd_hab2 %in% groups2[i]) %>%
    pivot_wider(names_from = Tax_curated, values_from = RPKM, id_cols = c(SeqId)) %>%
    as.data.frame(.) %>%
    filter(rowSums(select(., -1)) > 0)
  
  # Assign the result to the dynamically named variable
  assign(gr_name, current_group)
  

  OTU_sums <- get(gr_name)[,-1]
  row.names(OTU_sums) <- get(gr_name)$SeqId
  
  # Calculating Bray Curtis distance to use in clustering
  bray_curtis_dist <- vegan::vegdist(vegan::decostand(OTU_sums, method = "hellinger"))
  hclust_ward <- hclust(bray_curtis_dist, method = "ward.D2")
  ward_dendrogram <- as.dendrogram(hclust_ward)
  ward_order <- order.dendrogram(ward_dendrogram)
  assign(paste0("levels_", gr_name), hclust_ward$labels[order.dendrogram(ward_dendrogram)])

}

## creating a vector of the habitat variables
levels_hab2 <- c(levels_gr_1, levels_gr_2, levels_gr_3, levels_gr_4, levels_gr_5,
  levels_gr_6, levels_gr_7, levels_gr_8, levels_gr_9, levels_gr_10,
  levels_gr_11, levels_gr_12, levels_gr_13, levels_gr_14, levels_gr_15,
  levels_gr_16, levels_gr_17, levels_gr_18, levels_gr_19, levels_gr_20,
  levels_gr_21, levels_gr_22, levels_gr_23, levels_gr_24, levels_gr_25,
  levels_gr_26, levels_gr_27, levels_gr_28, levels_gr_29, levels_gr_30,
  levels_gr_31, levels_gr_32, levels_gr_33, levels_gr_34, levels_gr_35)


#Plotting the heatmaps
options("aplot_guides" = "keep")

 OTU_wide<-OTU_filtered_long%>%
  pivot_wider(names_from = Tax_curated, values_from = RPKM, id_cols = c(SeqId))%>%
  as.data.frame(.)%>%
  filter(rowSums(select(., -1)) > 0)
OTU_sums<-OTU_wide[,-1]
row.names(OTU_sums)<-OTU_wide$SeqId

groups<-unique(OTU_filtered_long$complex_long)

#Loading color palette
p_combine<-readRDS("palette_mfd_hab2.rds") 
names(p_combine)[names(p_combine) == "Activated sludge"] <- "Activated\nsludge"
names(p_combine)[names(p_combine) == "Non-native trees (exotic)"] <- "Non-native trees"
names(p_combine)[names(p_combine) == "Standing freshwater"] <- "Standing\nfreshwater"
names(p_combine)[names(p_combine) == "Running freshwater"] <- "Running\nfreshwater"
names(p_combine)[names(p_combine) == "Forests - no MFDO2"] <- "Forests\nno MFDO2"


#Plotting and saving one by one
for (i in seq_along(groups)) {
  # Use paste0 to dynamically create variable names gr_1, gr_2, etc.
  gr_name <- paste0("gr_", i)
  
  # Use the filter condition to select the current group
  current_group <- OTU_filtered_long %>%
    filter(complex_long %in% groups[i]) %>%
    pivot_wider(names_from = Tax_curated, values_from = RPKM, id_cols = c(SeqId)) %>%
    as.data.frame(.) %>%
    filter(rowSums(select(., -1)) > 0)
  
  # Assign the result to the dynamically named variable
  assign(gr_name, current_group)
  
  # Continue with the rest of your code using the dynamic variable name
  OTU_sums <- get(gr_name)[,-1]
  row.names(OTU_sums) <- get(gr_name)$SeqId
  
  # Continue with the analysis for the current group
  
  level <- get(paste0("levels_", gr_name))
  
  
hab2_sort<-OTU_filtered_long%>%
      mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
    mutate(Tax_curated=factor(Tax_curated, levels = rev(tax2$Tax_order), ordered=TRUE))%>%
  arrange(SeqId) 
  
 
heat <- OTU_filtered_long %>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(Tax_curated=factor(Tax_curated, levels = rev(tax2$Tax_order), ordered=TRUE))%>%
  filter(complex_long %in% groups[i]) %>%
    ggplot(aes(y = Tax_curated, x = SeqId, fill = RPKM)) +
    geom_tile(width=10.0, linewidth=0.0) +
    scale_fill_viridis_c(name = "RPKM", trans = "sqrt", limits = c(0.00, 10), na.value = "#F8E622", guide = "none",)+
   # scale_fill_gradient(low = "white", high = "darkred",  na.value = "darkred", name = "RPKM", trans = "sqrt", limits = c(0.00, 3), guide = "none") +
    labs(x = "", y = "") +
    facet_grid(type ~ complex, scales = "free", space = "free") +
    #theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          #axis.text.y = element_text(size = 8.5),
          legend.position = "right",
          strip.text = element_text(size = 10, face = "bold"),
          strip.background = element_rect(fill = "grey90", color = "white"),
          strip.background.y = element_blank(),
          strip.text.y.right = element_blank(),
          text = element_text(family = "Arial"))+
   scale_y_discrete(expand = c(0,0))
  
  
 tile <- OTU_filtered_long %>%
  filter(complex_long %in% groups[i]) %>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>% 
  mutate(mfd_hab2 = factor(mfd_hab2, levels = unique(hab2_sort$mfd_hab2))) %>%
  ggplot(., aes(x = SeqId, fill = mfd_hab2, y = 1)) +
  geom_tile()+
  theme(legend.position = "bottom",
        legend.text = element_text(size=9), legend.key.size = unit(0.5, "cm"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        text = element_text(family = "Arial")) +
  scale_fill_manual(values = p_combine)+
   scale_y_continuous(expand = c(0,0)) + 
    guides(fill = guide_legend(nrow = 4, title.theme = element_blank())) 
  
assign(paste0("HM_plot_", gsub(", ", "_", groups[i])),
       heat %>% aplot::insert_bottom(tile, height = 0.02))
}



#saving the plots

ggsave("Heatmaps/HM_Sediment_Urban_Freshwater.png",
       HM_plot_Sediment_Urban_Freshwater,
         height = 10,
  width = 1.9)


ggsave("Heatmaps/HM_Soil_Agriculture_Fields.png",
       HM_plot_Soil_Agriculture_Fields,
         height = 10,
  width = 6.8)


ggsave("Heatmaps/HM_Soil_Natural_Bogs.png",
       `HM_plot_Soil_Natural_Bogs_mires and fens`,
         height = 10,
  width = 1.8)


ggsave("Heatmaps/HM_Soil_Natural_Forests.png",
       `HM_plot_Soil_Natural_Forests`,
         height = 10,
  width = 3.5)



ggsave("Heatmaps/HM_Soil_Natural_Grassland.png",
       `HM_plot_Soil_Natural_Grassland formations`,
         height = 10,
  width = 3.7)


ggsave("Heatmaps/HM_Soil_Urban_Greenspaces.png",
       HM_plot_Soil_Urban_Greenspaces,
         height = 10,
  width = 2)



# Making the outermost heatmaps seperately


heat <- OTU_filtered_long %>%
    filter(complex_long %in% groups[7]) %>%
    mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
    mutate(Tax_curated=factor(Tax_curated, levels = rev(tax2$Tax_order), ordered=TRUE))%>%
    ggplot(aes(y = Tax_curated, x = SeqId, fill = RPKM)) +
    geom_tile(width=10.0, linewidth=0.0) +
    scale_fill_viridis_c(name = "RPKM", trans = "sqrt", limits = c(0.00, 10), na.value = "#F8E622", guide = "none",)+
   # scale_fill_gradient(low = "white", high = "darkred",  na.value = "darkred", name = "RPKM", trans = "sqrt", limits = c(0.00, 3), guide = "none") +
    labs(x = "", y = "") +
    facet_grid(type ~ complex, scales = "free", space = "free") +
    #theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          #axis.text.y = element_blank(),
          axis.text.y = element_text(size = 9),
          legend.position = "right",
          strip.text = element_text(size = 10, face = "bold"),
          strip.background = element_rect(fill = "grey90", color = "white"),
          strip.background.y = element_blank(),
          strip.text.y.right = element_blank(),
          text = element_text(family = "Arial"),
          plot.background = element_rect(fill = "transparent"))+
   scale_y_discrete(expand = c(0,0))
  

ggsave("Heatmaps/HM_Sediment_Natural_Freshwater2.png",
       heat,
        height = 8.55,
  width = 3)
  
 tile <- OTU_filtered_long %>%
  filter(complex_long %in% groups[7]) %>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>% 
  mutate(mfd_hab2 = factor(mfd_hab2, levels = unique(hab2_sort$mfd_hab2))) %>%
  ggplot(., aes(x = SeqId, fill = mfd_hab2, y = 1)) +
  geom_tile()+
  theme(legend.position = "bottom",       legend.text = element_text(size=9), legend.key.size = unit(0.5, "cm"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent")) +
  scale_fill_manual(values = p_combine)+
  scale_y_continuous(expand = c(0,0)) + 
  guides(fill = guide_legend(nrow = 4, title.theme = element_blank())) 
  
assign(paste0("HM_plot_", gsub(", ", "_", groups[7])),
       heat %>% aplot::insert_bottom(tile, height = 0.02))

ggsave("Heatmaps/HM_Sediment_Natural_Freshwater.png",
       HM_plot_Sediment_Natural_Freshwater,
        height = 10,
  width = 4.5)

ggsave("Heatmaps/HM_Sediment_Natural_Freshwater.svg",
       HM_plot_Sediment_Natural_Freshwater,
        height = 10,
  width = 4.5)


### Now to wastewater=gr8
heat <- OTU_filtered_long %>%
    filter(complex_long %in% groups[8]) %>%
    mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
    mutate(Tax_curated=factor(Tax_curated, levels = rev(tax2$Tax_order), ordered=TRUE))%>%
    mutate(mfd_hab2 = factor(mfd_hab2, levels = unique(hab2_sort$mfd_hab2))) %>%
    mutate(complex=gsub("Wastewater", "WWTP", complex))%>%
    ggplot(aes(y = Tax_curated, x = SeqId, fill = RPKM)) +
    geom_tile(width=10.0, linewidth=0.0) +
    scale_fill_viridis_c(name = "RPKM", trans = "sqrt", limits = c(0.00, 10), na.value = "#F8E622")+
    labs(x = "", y = "") +
    facet_grid(type ~ complex, scales = "free", space = "free") +
    #theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          #axis.text.y = element_text(size = 8.5),
          axis.text.y = element_blank(),
          legend.position = "right",
          strip.text = element_text(size = 10, face = "bold"),
          strip.background = element_rect(fill = "grey90", color = "white"),
          text = element_text(family = "Arial"))+
          #strip.background.y = element_blank(),
          #strip.text.y.right = element_blank())+
   scale_y_discrete(expand = c(0,0))
  
 tile <- OTU_filtered_long %>%
  filter(complex_long %in% groups[8]) %>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>% 
  mutate(mfd_hab2 = factor(mfd_hab2, levels = unique(hab2_sort$mfd_hab2))) %>%
  ggplot(., aes(x = SeqId, fill = mfd_hab2, y = 1)) +
  geom_tile()+
  theme(legend.position = "bottom",        legend.text = element_text(size=9), legend.key.size = unit(0.5, "cm"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        text = element_text(family = "Arial")) +
  scale_fill_manual(values = p_combine)+
   scale_y_continuous(expand = c(0,0)) + 
    guides(fill = guide_legend(nrow = 4, title.theme = element_blank())) 
  
assign(paste0("HM_plot_", gsub(", ", "_", groups[8])),
       heat %>% aplot::insert_bottom(tile, height = 0.02))



ggsave("Heatmaps/HM_Water_Urban_Wastewater.png",
       HM_plot_Water_Urban_Wastewater,
         height = 10,
  width = 1.9)

ggsave("Heatmaps/HM_Water_Urban_Wastewater.svg",
       HM_plot_Water_Urban_Wastewater,
         height = 10,
  width = 1.9)




######### Adding SingleM for nitrifiers #####


SingleM<-readRDS("SingleM_nitrifiers.RData")


 samples1<-SingleM%>%
   select(flat_name)%>%
   distinct()

samples2<-OTU_filtered_long%>%
   select(SeqId)%>%
   distinct()
 
### When filtering for singleM, some samples have rel. abundance of zero, and are not included. Thus, those are added now
abundance_zero_in_SingleM<- anti_join(samples2, samples1, by=join_by(SeqId==flat_name))
  
add<-OTU_filtered_long%>%
  filter(SeqId %in% abundance_zero_in_SingleM$SeqId)%>%
  select(!starts_with("Tax"))%>%
  select(!c("RPKM", "type"))%>%
  rename(flat_name=SeqId)%>%
  mutate(total_abundance=0)

distinct_combinations <- unique(SingleM[c("label_tax", "group")])

# Merge the combinations with the ad dataframe, keeping all rows from ad
combined_add <- merge(add, distinct_combinations, all.x = TRUE)%>%distinct()


SingleM<-SingleM%>%
  rbind(combined_add) %>%
  mutate(group_2=if_else(grepl("g__Nitrosospira|g__Nitrosomonas", label_tax), "AOB", group))%>%
  mutate(group_2=if_else(grepl("Ammonia\noxidizers", group_2), "AOA", group_2))%>%
  mutate(group_2=if_else(group_2=="Ammonia &\nnitrite\noxidizers", "CMX", group_2),
         group_2=if_else(group_2=="Nitrite\noxidizers", "NOB", group_2),
mutate(mfd_hab2 = if_else(mfd_hab2=="Forests - no MFD_02", paste0(mfd_hab1, "\nno MFDO2"), mfd_hab2))%>%  
  mutate(mfd_hab2=gsub("\\s*\\(exotic\\)", "", mfd_hab2))%>%
  mutate(mfd_hab2=gsub("Standing freshwater", "Standing\nfreshwater", mfd_hab2))%>%
  mutate(mfd_hab2=gsub("Running freshwater", "Running\nfreshwater", mfd_hab2))%>%
  mutate(mfd_hab2=gsub("Activated sludge", "Activated\nsludge", mfd_hab2))  %>%
  mutate(label_tax=if_else(label_tax=="d__Archaea ; o__Nitrososphaerales ; NA ; NA", paste0(label_tax), gsub("^d__.+?; o__.+?; ", "", label_tax)))%>%
  mutate(label_tax=gsub("d__Archaea ; o__Nitrososphaerales ; NA ; NA", "o__Nitrososphaerales ; NA ; NA", label_tax))%>%
  distinct()


#Looping through all the habitats:

for (i in seq_along(groups)) {
  # Use paste0 to dynamically create variable names gr_1, gr_2, etc.
  gr_name <- paste0("gr_", i)

 
heat <- SingleM %>%
  mutate(flat_name = factor(flat_name, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(group_2=factor(group_2, levels = c("AOA", "AOB", "NOB", "CMX"), ordered=TRUE))%>%
  filter(complex_long %in% groups[i]) %>%
    ggplot(aes(y = label_tax, x = flat_name, fill = total_abundance)) +
    geom_tile(width=10.0, linewidth=0.0) +
    scale_fill_viridis_c(name = "RPKM", trans = "sqrt", limits = c(0.00, 12), na.value = "#F8E622", guide = "none",)+
   # scale_fill_gradient(low = "white", high = "darkred",  na.value = "darkred", name = "RPKM", trans = "sqrt", limits = c(0.00, 3), guide = "none") +
    labs(x = "", y = "") +
    facet_grid(group_2 ~ complex, scales = "free", space = "free") +
    #theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          #axis.text.y = element_text(size = 8.5),
          legend.position = "right",
          strip.text = element_text(size = 9.3, face = "bold"),
          strip.background = element_rect(fill = "grey90", color = "white"),
          strip.background.y = element_blank(),
          strip.text.y.right = element_blank(),
          text = element_text(family = "Arial"))+
   scale_y_discrete(expand = c(0,0))
  
  
 tile <- SingleM %>%
  filter(complex_long %in% groups[i]) %>%
  mutate(flat_name = factor(flat_name, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(mfd_hab2 = factor(mfd_hab2, levels = unique(hab2_sort$mfd_hab2))) %>%
  ggplot(., aes(x = flat_name, fill = mfd_hab2, y = 1)) +
  geom_tile()+
  theme(legend.position = "bottom",   legend.text = element_text(size=9), legend.key.size = unit(0.5, "cm"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        text = element_text(family = "Arial")) +
  scale_fill_manual(values = p_combine)+
   scale_y_continuous(expand = c(0,0)) + 
    guides(fill = guide_legend(nrow = 4, title.theme = element_blank())) 
  
assign(paste0("SM_", gsub(", ", "_", groups[i])),
       heat %>% aplot::insert_bottom(tile, height = 0.04))
}







#saving the files




ggsave("SM_Heatmaps/Sediment_Urban_Freshwater.png",
       SM_Sediment_Urban_Freshwater,
         height = 6,
  width = 2)


ggsave("SM_Heatmaps/Soil_Agriculture_Fields.png",
       SM_Soil_Agriculture_Fields,
         height = 6,
  width = 7)


ggsave("SM_Heatmaps/Soil_Natural_Bogs.png",
       `SM_Soil_Natural_Bogs_mires and fens`,
         height = 6,
  width = 1.8)


ggsave("SM_Heatmaps/Soil_Natural_Forests.png",
       `SM_Soil_Natural_Forests`,
         height = 6,
  width = 3.5)



ggsave("SM_Heatmaps/Soil_Natural_Grassland.png",
       `SM_Soil_Natural_Grassland formations`,
         height = 6,
  width = 3.7)


ggsave("SM_Heatmaps/Soil_Urban_Greenspaces.png",
       SM_Soil_Urban_Greenspaces,
         height = 6,
  width = 2)



### And the outermost plots 
heat <- SingleM %>%
    filter(complex_long %in% groups[7]) %>%
  mutate(flat_name = factor(flat_name, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(group_2=factor(group_2, levels = c("AOA", "AOB", "NOB", "CMX"), ordered=TRUE))%>%
    ggplot(aes(y = label_tax, x = flat_name, fill = total_abundance)) +
    geom_tile(width=10.0, linewidth=0.0) +
    scale_fill_viridis_c(name = "RPKM", trans = "sqrt", limits = c(0.00, 12), na.value = "#F8E622", guide = "none",)+
   # scale_fill_gradient(low = "white", high = "darkred",  na.value = "darkred", name = "RPKM", trans = "sqrt", limits = c(0.00, 3), guide = "none") +
    labs(x = "", y = "") +
    facet_grid(group_2 ~ complex, scales = "free", space = "free") +
    #theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          #axis.text.y = element_blank(),
          axis.text.y = element_text(size = 9),
          legend.position = "right",
          strip.text = element_text(size = 9.3, face = "bold"),
          strip.background = element_rect(fill = "grey90", color = "white"),
          strip.background.y = element_blank(),
          strip.text.y.right = element_blank(),
          text = element_text(family = "Arial"),
          plot.background = element_rect(fill = "transparent"))+
   scale_y_discrete(expand = c(0,0))
  
  
 tile <- SingleM %>%
  filter(complex_long %in% groups[7]) %>%
  mutate(flat_name = factor(flat_name, levels = levels_hab2, ordered = TRUE)) %>% 
  mutate(mfd_hab2 = factor(mfd_hab2, levels = unique(hab2_sort$mfd_hab2))) %>%
  ggplot(., aes(x = flat_name, fill = mfd_hab2, y = 1)) +
  geom_tile()+
  theme(legend.position = "bottom", legend.text = element_text(size=9), legend.key.size = unit(0.5, "cm"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent")) +
  scale_fill_manual(values = p_combine)+
  scale_y_continuous(expand = c(0,0)) + 
  guides(fill = guide_legend(nrow = 4, title.theme = element_blank())) 
  
assign(paste0("SM_", gsub(", ", "_", groups[7])),
       heat %>% aplot::insert_bottom(tile, height = 0.04))

ggsave("SM_Heatmaps/Sediment_Natural_Freshwater.png",
       SM_Sediment_Natural_Freshwater,
        height = 6,
  width = 4.1)



## wastewater=gr8
 
 heat <- SingleM %>%
    filter(complex_long %in% groups[8]) %>%
  mutate(flat_name = factor(flat_name, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(group_2=factor(group_2, levels = c("AOA", "AOB", "NOB", "CMX"), ordered=TRUE))%>%
       mutate(complex=gsub("Wastewater", "WWTP", complex))%>%
    ggplot(aes(y = label_tax, x = flat_name, fill = total_abundance)) +
    geom_tile(width=10.0, linewidth=0.0) +
    scale_fill_viridis_c(name = "Rel.\nAbun\n[%]", trans = "sqrt", limits = c(0.00, 12), na.value = "#F8E622")+
    labs(x = "", y = "") +
    facet_grid(group_2 ~ complex, scales = "free", space = "free") +
    ##theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "right",
          legend.text = element_text(size=9),
          legend.title = element_text(size=10),
          strip.text = element_text(size = 10, face = "bold"),
          strip.background = element_rect(fill = "grey90", color = "white"),
        #  strip.background.y = element_blank(),
        #  strip.text.y.right = element_blank(),
          text = element_text(family = "Arial"),
          plot.background = element_rect(fill = "transparent"))+
   scale_y_discrete(expand = c(0,0))
 
 

  
 tile <- SingleM %>%
  filter(complex_long %in% groups[8]) %>%
  mutate(flat_name = factor(flat_name, levels = levels_hab2, ordered = TRUE)) %>% 
  mutate(mfd_hab2 = factor(mfd_hab2, levels = unique(hab2_sort$mfd_hab2))) %>%
  ggplot(., aes(x = flat_name, fill = mfd_hab2, y = 1)) +
  geom_tile()+
  theme(legend.position = "bottom", legend.text = element_text(size=9), legend.key.size = unit(0.5, "cm"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent")) +
  scale_fill_manual(values = p_combine)+
  scale_y_continuous(expand = c(0,0)) + 
  guides(fill = guide_legend(nrow = 4, title.theme = element_blank())) 
  
assign(paste0("SM_", gsub(", ", "_", groups[8])),
       heat %>% aplot::insert_bottom(tile, height = 0.04))

ggsave("SM_Heatmaps/Water_Urban_Wastewater.png",
       SM_Water_Urban_Wastewater,
        height = 6,
  width = 1.84)

ggsave("SM_Heatmaps/Water_Urban_Wastewater.svg",
       SM_Water_Urban_Wastewater,
        height = 6,
  width = 1.84)




#### To generate the colour palette, that is also present as rds object.
p_combine <- structure(list(`Poales, Cereal` = "#cdad00", 
                            `Mixed crops` = "#32849f", 
                            `Poales, grass` = "#7301A8FF", 
                            Malvids = "#33ccff", 
                            Superasterids = "#F8E622", 
                            Fabids = "#ff9999", 
                            Asparagales = "#CC4678FF", 
                            `Running freshwater` = "#cdad00", 
                            `Standing freshwater` = "#32849f", 
                            `Deciduous trees` = "#0D0887FF" , 
                            `Forests - no MFDO2` = "#32849f", 
                            `Non-native trees (exotic)` = "#cd2626",
                            `Coniferous forest` = "#cc33ff", 
                             Beech = "#33ccff", 
                            `Alluvial woodland` = "#F8E622", 
                             Oak = "#ff9999", 
                            `Bog woodland` = "#7301A8FF",
                             Willow = "#cdad00", 
                            Spruce = "#000000",
                            `Mixed crops` = "#32849f", 
                            Fallow = "#0D0887FF", 
                            Asterids = "#cd2626",
                            `Running freshwater` = "#cdad00", 
                            `Standing freshwater` = "#32849f", 
                            `Rainwater basin, City` = "#cdad00",
                            `Rainwater basin, Dried` = "#32849f", 
                            `Enclosed water` = "#7301A8FF", 
                            `Rainwater basin, Roadside` = "#cd2626", 
                            `Semi-natural dry grasslands` = "#32849f", 
                            `Semi-natural tall-herb humid meadows` = "#cdad00", 
                            `Natural grasslands` = "#7301A8FF", 
                            Parks = "#cdad00",
                            Other = "#32849f", 
                            `Calcareous fens` = "#cdad00", 
                            `Sphagnum acid bogs` = "#32849f", 
                            `Wet thicket` = "#7301A8FF", 
                            `Activated sludge` = "#cdad00",
                            Influent = "#32849f"), class = "character", 
                             row.names = c("Poales, Cereal", 
                                     "Mixed crops", 
                                     "Poales, grass", 
                                     "Malvids", 
                                     "Superasterids", 
                                     "Fabids", 
                                     "Asparagales",
                                     "Running freshwater", 
                                     "Standing freshwater", 
                                     "Deciduous trees", 
                                     "Forests - no MFDO2",
                                     "Non-native trees (exotic)", 
                                     "Coniferous forest",
                                     "Beech", 
                                     "Alluvial woodland", 
                                     "Oak", 
                                     "Bog woodland", 
                                     "Willow", 
                                     "Spruce",
                                     "Mixed crops", 
                                     "Fallow", 
                                     "Asterids", 
                                     "Running freshwater", 
                                     "Standing freshwater", 
                                     "Rainwater basin, City",
                                     "Rainwater basin, Dried", 
                                     "Enclosed water", 
                                     "Rainwater basin, Roadside", 
                                     "Semi-natural dry grasslands",
                                     "Semi-natural tall-herb humid meadows",
                                     "Natural grasslands", 
                                     "Parks",
                                     "Other", 
                                     "Calcareous fens", 
                                     "Sphagnum acid bogs", 
                                     "Wet thicket", 
                                     "Activated sludge", 
                                     "Influent"))


saveRDS(p_combine, "palette_mfd_hab2.rds")

