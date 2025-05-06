#title: "02_Nitrifier_Heatmaps"
#author: "Kalinka Sand Knudsen"
#update: "2025-04-23"


#Loading packages
.libPaths(c("path_to_R_package_library", .libPaths()))
library(ggplot2)
library(vroom)
library(tidyverse)
library(readxl)

#Set WD
setwd("path_to_working_directory")

#Load dataframe generated with "normalising_and_aggregating_nitrifier_genes"
OTU_filtered_long<-readRDS("output/OTU_filtered_long_25_02_21.rds")

#Ordering the tax:
tax_order<-OTU_filtered_long%>%
  select(Tax, Tax_short, type)%>%
  distinct()%>%
  arrange(Tax)%>%
  arrange(type)%>%
  mutate(Tax_curated=gsub("o_", "", Tax_short),
         Tax_curated=gsub("f_", "", Tax_curated),
         Tax_curated=gsub("g_", "", Tax_curated))%>%
  mutate(Tax_curated = if_else(Tax_short %in% c("Nitrosomonas", "BOG-931", "Bradyrhizobium", "g_TA21", "Nitrospira_A", "Palsa-1315"), paste0(Tax_curated, "*"), Tax_curated))%>%
  mutate(Tax_curated=gsub("Nitrosopumilus_1", "Nitrosopumilus", Tax_curated),
         Tax_curated=gsub("TH5893_1", "TH5893", Tax_curated),
         Tax_curated=gsub("TH5896_1", "TH5896", Tax_curated),
         Tax_curated=gsub("Nitrospira_C1", "Nitrospira_C", Tax_curated),
         Tax_curated=gsub("cluster", "clade", Tax_curated),
         Tax_curated=gsub("Gammaproteobacteria_amoA", "Gammaproteobacteria", Tax_curated))

Tax_order_curated<-factor(tax_order$Tax_curated,
  levels=c("Nitrososphaerales", "Nitrosopumilaceae", "Nitrosotalea", "Nitrosarchaerum", 
  "Nitrosopumilus", "Nitrososphaeraceae",  "Nitrosocosmicus", "TH1177", "TH5896",  "TA21_Nitrosopolaris_clade", 
  "JAFAQB01", "Nitrosopolaris", "TA21*", "TH5893_Nitrososphaera_clade", 
  "Nitrososphaera", "TH5893", 
  "Gammaproteobacteria", "Nitrosomonas*", "Nitrosospira", "Nitrospira", 
  "Nitrospira_clade_A", "Nitrospira_clade_B", "Thiocapsa", 
  "Nitrococcus_Nitrobacter_clade", "Nitrococcus", "Nitrobacter-like_NxrA", "BOG-931*", 
  "Bradyrhizobium_Nitrobacter", "Bradyrhizobium*", "Nitrobacter", 
  "Pseu_JAFA_VAZQ_clade",
  "Pseudolabrys", "VAZQ01", "Nitrotoga", "NS_4_NS_12", "UBA8639",
  "Nitrospira_umbrella", "Nitrospira_clade_1", "Nitrospira_A*", "Nitrospira_C", 
  "Nitrospira_clade_2", "Palsa-1315*", "Nitrospira_D_inopinata", "Nitrospira_D_lenta"), ordered=TRUE
)


# Merging the curated tax onto the dataframe:
OTU_filtered_long<-left_join(OTU_filtered_long, tax_order)


# Filtering of metadata and generation of labels
OTU_filtered_long<-OTU_filtered_long %>%
  mutate(mfd_areatype=gsub("Agriculture \\(reclaimed lowland\\)", "Agriculture", mfd_areatype))  %>%
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
groups<-unique(OTU_filtered_long$complex_long)


# Clustering of MFD_habo2 individually

OTU_filtered_long$mfd_hab2 <- factor(OTU_filtered_long$mfd_hab2, levels = c(
  "Poales, Cereal", "Mixed crops", "Asparagales", "Malvids", "Asterids", "Superasterids", "Fabids", 
  "Poales, grass", "Fallow", "Parks", "Other", 
  "Semi-natural tall-herb humid meadows", "Semi-natural dry grasslands", "Natural grasslands", 
  "Alluvial woodland", "Bog woodland", "Non-native trees", "Beech", "Oak", "Coniferous forest", 
  "Deciduous trees", "Forests\nno MFDO2", "Wet thicket", "Calcareous fens", "Sphagnum acid bogs", 
  "Running\nfreshwater", "Standing\nfreshwater", "Enclosed water", 
  "Rainwater basin, City", "Rainwater basin, Roadside", "Rainwater basin, Dried", 
  "Activated\nsludge", "Influent"
))

# Arrange the data frame based on the factor levels
OTU_filtered_long <- OTU_filtered_long %>%
  arrange(mfd_hab2)




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
  
  # Continue with the rest of your code using the dynamic variable name
  OTU_sums <- get(gr_name)[,-1]
  row.names(OTU_sums) <- get(gr_name)$SeqId
  
  bray_curtis_dist <- vegan::vegdist(vegan::decostand(OTU_sums, method = "hellinger"))
  hclust_ward <- hclust(bray_curtis_dist, method = "ward.D2")
  ward_dendrogram <- as.dendrogram(hclust_ward)
  ward_order <- order.dendrogram(ward_dendrogram)
  assign(paste0("levels_", gr_name), hclust_ward$labels[order.dendrogram(ward_dendrogram)])
  
  
 # OTU_filtered_long<-OTU_filtered_long%>%
  #   mutate(mfd_hab2=factor(mfd_hab2, levels = groups2[i]), ordered = TRUE)
  
}

levels_hab2 <- c(levels_gr_1, levels_gr_2, levels_gr_3, levels_gr_4, levels_gr_5,
  levels_gr_6, levels_gr_7, levels_gr_8, levels_gr_9, levels_gr_10,
  levels_gr_11, levels_gr_12, levels_gr_13, levels_gr_14, levels_gr_15,
  levels_gr_16, levels_gr_17, levels_gr_18, levels_gr_19, levels_gr_20,
  levels_gr_21, levels_gr_22, levels_gr_23, levels_gr_24, levels_gr_25,
  levels_gr_26, levels_gr_27, levels_gr_28, levels_gr_29, levels_gr_30,
  levels_gr_31, levels_gr_32, levels_gr_33)

saveRDS(levels_hab2,"output/levels_hab2_25_02_21.rds")

#Plotting the heatmaps

 OTU_wide<-OTU_filtered_long%>%
  pivot_wider(names_from = Tax_curated, values_from = RPKM, id_cols = c(SeqId))%>%
  as.data.frame(.)%>%
  filter(rowSums(select(., -1)) > 0)
OTU_sums<-OTU_wide[,-1]
row.names(OTU_sums)<-OTU_wide$SeqId

groups<-unique(OTU_filtered_long$complex_long)


############# Generating colour palette ##################
unique(OTU_filtered_long$mfd_hab2)
p_combine <- structure(list(`Poales, Cereal` = "#cdad00", 
                            `Mixed crops` = "#32849f", 
                            `Poales, grass` = "#7301A8FF", 
                            Malvids = "#33ccff", 
                            Superasterids = "#F8E622", 
                            Fabids = "#ff9999", 
                            Fallow = "#CC4678FF", 
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
                            Asparagales = "#0D0887FF", 
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


saveRDS(p_combine, "output/palette_mfd_hab2.rds")





p_combine<-readRDS("output/palette_mfd_hab2.rds")
names(p_combine)[names(p_combine) == "Activated sludge"] <- "Activated\nsludge"
names(p_combine)[names(p_combine) == "Non-native trees (exotic)"] <- "Non-native trees"
names(p_combine)[names(p_combine) == "Standing freshwater"] <- "Standing\nfreshwater"
names(p_combine)[names(p_combine) == "Running freshwater"] <- "Running\nfreshwater"
names(p_combine)[names(p_combine) == "Forests - no MFDO2"] <- "Forests\nno MFDO2"


intersect(unique(OTU_filtered_long$mfd_hab2), names(p_combine))
setdiff(unique(OTU_filtered_long$mfd_hab2), names(p_combine))



OTU_filtered_long<-OTU_filtered_long%>%mutate(type=gsub("amoA", "*amoA*", type))%>%
                                      mutate(type=gsub("cNxrA", "*cnxrA*", type),
                                             type=gsub("pNxrA", "*pnxrA*", type))          

levels(OTU_filtered_long$type)<-c("Bacterial *amoA*",
                          "Archaeal *amoA*",
                          "*cnxrA*", "*pnxrA*")


heat <- OTU_filtered_long %>%
   mutate(Tax_curated=factor(Tax_curated, levels = rev(levels(Tax_order_curated)), ordered=TRUE))%>%
   filter(mfd_hab1 %in% c("Fields", "Bogs, mires and fens", "Forests", "Grassland formations", "Greenspaces"))%>%
     mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(type=factor(type, levels=c("Archaeal *amoA*", "Bacterial *amoA*","*cnxrA*","*pnxrA*")))%>%
  #mutate(complex=gsub("Wastewater", "WWTP", complex))%>%
    ggplot(aes(y = Tax_curated, x = SeqId, fill = RPKM)) +
    geom_tile() +
   scale_fill_gradientn(
    name = "RPKM",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 8),  # Set the limits manually
    na.value = "black"  # Handling values larger than the limit
  ) +
    #scale_fill_viridis_c(name = "RPKM", trans = "sqrt", limits = c(0.00, 10), na.value = "#F8E622")+
  #labs(x = "", y = "") +
  facet_grid(type ~ factor(complex, levels=c("Soil\nAgriculture\nFields", "Soil\nUrban\nGreenspaces", "Soil\nNatural\nGrassland formations","Soil\nNatural\nForests", "Soil\nNatural\nBogs, mires and fens")), scales = "free", space = "free", switch = "y") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        #  axis.text.y = element_blank(),
        axis.text.y = element_text(size = 7, margin = margin(r = -2)),
        legend.key.size = unit(0.38, "cm"),
        legend.title = element_text(size=7.5),
        legend.text = element_text(size=7),
        strip.clip = "off",
        strip.text.x = element_text(size = 7, face="bold"),  # Size for x-axis facet labels
        strip.text.y = element_text(size = 6, face="bold"),
        strip.text.y.left = ggtext::element_markdown(),
        panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2) ,
        strip.background = element_blank(),
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.15, "lines", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))
heat

bar_plot <- OTU_filtered_long %>%
     filter(mfd_hab1 %in% c("Fields", "Bogs, mires and fens", "Forests", "Grassland formations", "Greenspaces"))%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  #mutate(mfd_hab2 = factor(mfd_hab2, levels = unique(hab2_sort$mfd_hab2))) %>%
  ggplot(aes(x = SeqId, fill = mfd_hab2)) +
  geom_tile(aes(y = 1)) +
  facet_grid(. ~ factor(complex, levels=c("Soil\nAgriculture\nFields", "Soil\nUrban\nGreenspaces", "Soil\nNatural\nGrassland formations","Soil\nNatural\nForests", "Soil\nNatural\nBogs, mires and fens")), scales = "free", space = "free") +
  scale_fill_manual(values = p_combine) + # Specify your colors
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "none",
    panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2) ,
    panel.spacing = unit(0.15, "lines", data = NULL)
  )+   scale_y_continuous(expand = c(0,0)) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))



# Step 3: Combine the heatmap plot and the bar plot

combined_plot <- heat / bar_plot + plot_layout(heights = c(20, 0.4)) + theme(panel.background = element_blank(), plot.background = element_blank())
# 
# ggsave("output/Nitrifiers_review_24_12_10.png",
#        combined_plot,
#          height = 6,
#   width = 11, dpi=600)
#combined_plot







############# Adding in Single M ##############################

SingleM<-vroom("./data/SingleM_nitrifiers.csv", delim=",")

samples1<-SingleM%>%
  # filter(mfd_hab2=="Influent")%>%
  select(flat_name)%>%
  distinct()

samples2<-OTU_filtered_long%>%
  # filter(mfd_hab2=="Influent")%>%
  select(SeqId)%>%
  distinct()

### WHen filtering for singleM, some samples have rel. abundance of zero, and are not included. Thus, we add those now
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

meta<-OTU_filtered_long%>% select(!starts_with("Tax"))%>%
  select(!c("RPKM", "type"))%>%
  rename(flat_name=SeqId)

SingleM<-SingleM%>%
  left_join(meta)%>%
  rbind(combined_add) %>%
  mutate(group_2=if_else(grepl("g__Nitrosospira|g__Nitrosomonas", label_tax), "AOB", group))%>%
  mutate(group_2=if_else(grepl("Ammonia\noxidizers", group_2), "AOA", group_2))%>%
  mutate(group_2=if_else(group_2=="Ammonia &\nnitrite\noxidizers", "CMX", group_2),
         group_2=if_else(group_2=="Nitrite\noxidizers", "NOB", group_2),
         group_2=if_else(grepl("g__Nitrospira_D", label_tax), "CMX", group_2))%>%
  mutate(mfd_hab2 = if_else(mfd_hab2=="Forests - no MFD_02", paste0(mfd_hab1, "\nno MFDO2"), mfd_hab2))%>%  
  mutate(mfd_hab2=gsub("\\s*\\(exotic\\)", "", mfd_hab2))%>%
  mutate(mfd_hab2=gsub("Standing freshwater", "Standing\nfreshwater", mfd_hab2))%>%
  mutate(mfd_hab2=gsub("Running freshwater", "Running\nfreshwater", mfd_hab2))%>%
  mutate(mfd_hab2=gsub("Activated sludge", "Activated\nsludge", mfd_hab2))  %>%
  mutate(label_tax=if_else(label_tax=="d__Archaea ; o__Nitrososphaerales ; NA ; NA", paste0(label_tax), gsub("^d__.+?; o__.+?; ", "", label_tax)))%>%
  mutate(label_tax=gsub("d__Archaea ; o__Nitrososphaerales ; NA ; NA", "o__Nitrososphaerales ; NA ; NA", label_tax))%>%
  distinct()



SingleM$label_tax_short <- sapply(
  strsplit(as.character(SingleM$label_tax), " ; "), 
  function(x) tail(x[x != "NA"], 1)
)


Genera_keep<-SingleM%>%select(label_tax)%>%filter(grepl("g__JAFAQB01|g__JARBAU01|g__TH5896|g__TH5893|g__Nitrosomonas|g__Nitrosospira|g__Nitrosotalea|f__Nitrososphaeraceae ; NA|f__Nitrosopumilaceae ; NA|g__JAFAQB01|g__Nitrosocosmicus|g__Nitrosopolaris|g__Nitrososphaera|g__TA-21|g__TH1177|g__Nitrospira_C|g__Nitrospira_A|g__Nitrospira_D|g__Palsa-1315|g__Nitrobacter|g__Nitrococcus|o__Nitrososphaerales ; NA ; NA", label_tax))%>%
  distinct()%>%pull(label_tax)


SingleM_order<-SingleM%>%
  select(label_tax, label_tax_short, group_2)%>%
  filter(label_tax %in% Genera_keep)%>%
  distinct()%>%
  arrange(label_tax_short)%>%
  arrange(group_2)



SingleM_order_curated<-factor(SingleM_order$label_tax_short,
                              levels = c("o__Nitrososphaerales", "f__Nitrosopumilaceae", "g__Nitrosotalea",
                                         "f__Nitrososphaeraceae", "g__Nitrosopolaris", "g__Nitrosocosmicus",
                                         "g__Nitrososphaera", "g__TA-21", "g__TH5893", 
                                         "g__TH1177", "g__TH5896", "g__JAFAQB01", "g__JARBAU01",
                                         "g__Nitrosomonas", "g__Nitrosospira", "g__Nitrospira_A",
                                         "g__Nitrospira_C", "g__Nitrobacter", "g__Palsa-1315",
                                         "g__Nitrospira_D"), ordered=TRUE)






heat_SingleM <- SingleM %>%
  mutate(flat_name = factor(flat_name, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(group_2=factor(group_2, levels = c("AOA", "AOB", "NOB", "CMX"), ordered=TRUE))%>%
  filter(label_tax %in% Genera_keep)%>%
  mutate(label_tax_short=factor(label_tax_short, levels = rev(levels(SingleM_order_curated)), ordered=TRUE))%>%
  filter(mfd_hab1 %in% c("Fields", "Bogs, mires and fens", "Forests", "Grassland formations", "Greenspaces"))%>%
  #mutate(complex=gsub("Wastewater", "WWTP", complex))%>%
  ggplot(aes(y = label_tax_short, x = flat_name, fill = total_abundance)) +
  geom_tile() +
  scale_fill_gradientn(
    name = "Rel.\nAbund\n[%]",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 5),  # Set the limits manually
    na.value = "black"  # Handling values larger than the limit
  ) +
  #scale_fill_viridis_c(name = "RPKM", trans = "sqrt", limits = c(0.00, 10), na.value = "#F8E622")+
  #labs(x = "", y = "") +
  facet_grid(group_2 ~ factor(complex, levels=c("Soil\nAgriculture\nFields", "Soil\nUrban\nGreenspaces", "Soil\nNatural\nGrassland formations","Soil\nNatural\nForests", "Soil\nNatural\nBogs, mires and fens")), scales = "free", space = "free", switch = "y") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        #  axis.text.y = element_blank(),
        axis.text.y = element_text(size = 7, margin = margin(r = -2)),
        legend.key.size = unit(0.38, "cm"),
        legend.title = element_text(size=7.5),
        legend.text = element_text(size=7),
        legend.position = "right",
        strip.clip = "off",
        strip.text.x = element_blank(),
        #strip.text.x = element_text(size = 7, face="bold"),  # Size for x-axis facet labels
        strip.text.y = element_text(size = 6, face="bold"),
        panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2),
        strip.background = element_blank(),
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.15, "lines", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  theme(plot.margin = unit(c(0.3,0,0,0), "cm"))


combined_plot <- heat / heat_SingleM / bar_plot + plot_layout(heights = c(19, 8.8, 0.5))
#combined_plot



ggsave("output/Nitrifiers_SingleM_review_25_03_20.png",
       combined_plot,
         height = 7,
  width = 11, dpi=600)





################################ Getting all the legends from the plot above #####################################


p_combine<-readRDS("output/palette_mfd_hab2.rds")
names(p_combine)[names(p_combine) == "Activated sludge"] <- "Activated\nsludge"
names(p_combine)[names(p_combine) == "Non-native trees (exotic)"] <- "Non-native trees"
names(p_combine)[names(p_combine) == "Standing freshwater"] <- "Standing\nfreshwater"
names(p_combine)[names(p_combine) == "Running freshwater"] <- "Running\nfreshwater"



OTU_filtered_long<-OTU_filtered_long%>%
 mutate(mfd_hab2 = gsub("Forests\nno MFDO2", "Forests - no MFDO2", mfd_hab2))%>%
  arrange(mfd_hab2)


OTU_filtered_long$mfd_hab2 <- factor(OTU_filtered_long$mfd_hab2, levels = c(
  "Poales, Cereal", "Mixed crops", "Asparagales", "Malvids", "Asterids", "Superasterids", "Fabids", 
  "Poales, grass", "Fallow", "Parks", "Other", 
  "Semi-natural tall-herb humid meadows", "Semi-natural dry grasslands", "Natural grasslands", 
  "Alluvial woodland", "Bog woodland", "Non-native trees", "Beech", "Oak", "Coniferous forest", 
  "Deciduous trees", "Forests - no MFDO2", "Wet thicket", "Calcareous fens", "Sphagnum acid bogs", 
  "Running\nfreshwater", "Standing\nfreshwater", "Enclosed water", 
  "Rainwater basin, City", "Rainwater basin, Roadside", "Rainwater basin, Dried", 
  "Activated\nsludge", "Influent"
))

bar_plot <- OTU_filtered_long %>%
     filter(mfd_hab1 %in% c("Fields", "Bogs, mires and fens", "Forests", "Grassland formations", "Greenspaces"))%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  ggplot(aes(x = SeqId, fill = mfd_hab2)) +
  geom_tile(aes(y = 1)) +
  facet_grid(. ~ factor(complex, levels=c("Soil\nAgriculture\nFields", "Soil\nUrban\nGreenspaces", "Soil\nNatural\nGrassland formations","Soil\nNatural\nForests", "Soil\nNatural\nBogs, mires and fens")), scales = "free", space = "free") +
  scale_fill_manual(values = p_combine) + # Specify your colors
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(0.45, "cm"),
    panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2) ,
    panel.spacing = unit(0.15, "lines", data = NULL)
  )+   scale_y_continuous(expand = c(0,0)) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))



# Getting the legend
legend <- ggpubr::get_legend(bar_plot)
leg<-ggpubr::as_ggplot(legend)
ggsave("output/legend.svg", leg)

