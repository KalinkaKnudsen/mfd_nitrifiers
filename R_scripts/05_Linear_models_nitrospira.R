#title: "05_Linear_models_nitrospira.R"
#author: "Kalinka Sand Knudsen"
#update: "2025-04-23"


#Loading packages
.libPaths(c("path_to_R_package_library", .libPaths()))
library(ggplot2)
library(cowplot)
library(vroom)
library(tidyverse)
library(readxl)
library(patchwork)

#Set WD
setwd("path_to_working_directory")

OTU_filtered_long<-readRDS("output/OTU_filtered_long_25_02_21.rds")

extreme_samples<-OTU_filtered_long%>%
  filter(RPKM>100)

#Ordering the tax:

OTU_filtered_long<-OTU_filtered_long%>%mutate(Tax_curated=gsub("o_", "", Tax_short),
            Tax_curated=gsub("f_", "", Tax_curated),
            Tax_curated=gsub("g_", "", Tax_curated))%>%
  mutate(Tax_curated = if_else(Tax_short %in% c("Nitrosomonas", "BOG-931", "Bradyrhizobium", "g_TA21", "Nitrospira_A", "Palsa-1315"), paste0(Tax_curated, "*"), Tax_curated))%>%
  mutate(Tax_curated=gsub("Nitrosopumilus_1", "Nitrosopumilus", Tax_curated),
         Tax_curated=gsub("TH5893_1", "TH5893", Tax_curated),
         Tax_curated=gsub("TH5896_1", "TH5896", Tax_curated),
         Tax_curated=gsub("Nitrospira_C1", "Nitrospira_C", Tax_curated),
         Tax_curated=gsub("cluster", "clade", Tax_curated))




#Merging the curated tax onto the dataframe:

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
  filter(!mfd_hab2 %in% c("Mire", "Spruce", "Willow"))%>%
  filter(!is.na(Tax_curated))%>%
  mutate(mfd_hab2=gsub("\\s*\\(exotic\\)", "", mfd_hab2))%>%
  mutate(mfd_hab2=gsub("Standing freshwater", "Standing\nfreshwater", mfd_hab2))%>%
  mutate(mfd_hab2=gsub("Running freshwater", "Running\nfreshwater", mfd_hab2))%>%
  mutate(mfd_hab2=gsub("Activated sludge", "Activated\nsludge", mfd_hab2))


OTU_wide<-OTU_filtered_long%>%
  pivot_wider(names_from = Tax_curated, values_from = RPKM, id_cols = c(SeqId))%>%
  as.data.frame(.)%>%
  filter(rowSums(select(., -1)) > 0)
OTU_sums<-OTU_wide[,-1]
row.names(OTU_sums)<-OTU_wide$SeqId

#groups<-unique(OTU_filtered_long$complex_long)



# Correlation between Nitrospira_clade_B and Palsa-1315



#colors <- c("gold3",  "#33a08c","darkred", "darkgreen", "#1a0060","#697e47","#ff7f00",  "#5fa3ca", "#fb9a99", "darkblue", "orange", "#6a3d9a","#1f78b4")
p_combine<-readRDS("output/palette_mfd_hab2.rds")
names(p_combine)[names(p_combine) == "Activated sludge"] <- "Activated\nsludge"
names(p_combine)[names(p_combine) == "Non-native trees (exotic)"] <- "Non-native trees"
names(p_combine)[names(p_combine) == "Standing freshwater"] <- "Standing\nfreshwater"
names(p_combine)[names(p_combine) == "Running freshwater"] <- "Running\nfreshwater"
names(p_combine)[names(p_combine) == "Forests - no MFDO2"] <- "Forests\nno MFDO2"

#### Combination of all #####
Nitrospira<-OTU_filtered_long%>%
  filter(Tax_curated %in% c("Nitrospira_clade_B", "Palsa-1315*"))%>%
  select(SeqId, Tax_curated, RPKM, mfd_hab1, mfd_hab2)%>%
  pivot_wider(names_from = Tax_curated, values_from = RPKM)


Nitrospira_BOG<-Nitrospira%>%
  filter(mfd_hab1=="Bogs, mires and fens")


lm_cladeB_P <- lm(`Nitrospira_clade_B`~`Palsa-1315*`, data = Nitrospira_BOG)
summary(lm_cladeB_P)

num_samples <- n_distinct(Nitrospira_BOG$SeqId)

P_BOGS <-Nitrospira_BOG %>%
  filter(mfd_hab1=="Bogs, mires and fens") %>%
  ggplot(., aes(
    y = `Nitrospira_clade_B`,
    x= `Palsa-1315*`)) +
  geom_point(aes(fill = mfd_hab2), size = 2, alpha = 0.6, color = "black", pch = 21, stroke=.8)+
  scale_fill_manual(values=p_combine, name = paste0("Bogs, Mires and Fens,\nn = ",num_samples))+
  geom_smooth(method = "lm", se = FALSE, col = "darkred", lty = 2) +
  theme_classic() +
  labs(
    y = "Nitrospira clade B amoA [RPKM]",
    x = "Nitrospiraceae Palsa-1315 nxrA [RPKM]",
    subtitle = paste0("r-squared: ", round(summary(lm_cladeB_P)$r.squared, 3), "\np-value: ",formatC(summary(lm_cladeB_P)$coefficients[2,4], format = "e", digits = 2), "\ncoefficient: ", round(summary(lm_cladeB_P)$coefficients[2,1], 2))
  )+
  coord_cartesian(xlim = c(0, max(Nitrospira_BOG$`Palsa-1315*`)), ylim = c(0,max(Nitrospira_BOG$`Palsa-1315*`)))+
  theme(legend.position = c(0.6, 1),
        legend.background = element_rect(fill = "transparent"))+
  guides(fill = guide_legend(ncol = 2))

P_BOGS



Nitrospira_Agri<-Nitrospira%>%
  filter(mfd_hab1=="Fields")


boxplot(Nitrospira_Agri$`Nitrospira_clade_B`)


num_samples <- n_distinct(Nitrospira_Agri$SeqId)

lm_cladeB_P <- lm(`Nitrospira_clade_B`~`Palsa-1315*`, data = Nitrospira_Agri)
summary(lm_cladeB_P)


P_Agri <-Nitrospira_Agri %>%
  ggplot(., aes(
    y = `Nitrospira_clade_B`,
    x= `Palsa-1315*`)) +
  geom_point(aes(fill = mfd_hab2), size = 2, alpha = 0.6, color = "black", pch = 21, stroke=.8)+
  scale_fill_manual(values=p_combine, name = paste0("Agriculture,\nn = ",num_samples))+
 # geom_point(fill="white", alpha=0.9, pch=1, stroke=1) +
  geom_smooth(method = "lm", se = FALSE, col = "darkred", lty = 2) +
  theme_classic() +
  labs(
    y = "Nitrospira clade B amoA [RPKM]",
    x = "Nitrospiraceae Palsa-1315 nxrA [RPKM]",
    subtitle = paste0("r-squared: ", round(summary(lm_cladeB_P)$r.squared, 3), "\np-value: ",formatC(summary(lm_cladeB_P)$coefficients[2,4], format = "e", digits = 2), "\ncoefficient: ", round(summary(lm_cladeB_P)$coefficients[2,1], 2))
  )+
  #coord_cartesian(xlim = c(0, 3), ylim = c(0, 3))+
  theme(legend.position = c(0.6, 0.9),
        legend.background = element_rect(fill = "transparent"))+
  guides(fill = guide_legend(ncol = 2)) +
  coord_cartesian(xlim = c(0, max(Nitrospira_Agri$`Palsa-1315*`)), ylim = c(0,max(Nitrospira_Agri$`Palsa-1315*`)))




Nitrospira_Forest<-Nitrospira%>%
  filter(mfd_hab1=="Forests")

max(Nitrospira_Forest$`Palsa-1315*`)

num_samples <- n_distinct(Nitrospira_Forest$SeqId)

lm_cladeB_P <- lm(`Nitrospira_clade_B`~`Palsa-1315*`, data = Nitrospira_Forest)
summary(lm_cladeB_P)

P_Forest <-Nitrospira_Forest %>%
  ggplot(., aes(
    y = `Nitrospira_clade_B`,
    x= `Palsa-1315*`)) +
  #geom_point(fill="white", alpha=0.9, pch=1, stroke=1) +
  geom_point(aes(fill = mfd_hab2), size = 2, alpha = 0.6, color = "black", pch = 21, stroke=.8)+
  scale_fill_manual(values=p_combine, name = paste0("Forests,\nn = ",num_samples))+
  geom_smooth(method = "lm", se = FALSE, col = "darkred", lty = 2) +
  theme_classic() +
  labs(
    y = "Nitrospira clade B amoA [RPKM]",
    x = "Nitrospiraceae Palsa-1315 nxrA [RPKM]",
    subtitle = paste0("r-squared: ", round(summary(lm_cladeB_P)$r.squared, 3), "\np-value: ",formatC(summary(lm_cladeB_P)$coefficients[2,4], format = "e", digits = 2), "\ncoefficient: ", round(summary(lm_cladeB_P)$coefficients[2,1], 2))
  )+
  coord_cartesian(xlim = c(0, max(Nitrospira_Forest$`Palsa-1315*`)), ylim = c(0,max(Nitrospira_Forest$`Palsa-1315*`)))+
  theme(legend.position = c(0.6, 0.9),
        legend.background = element_rect(fill = "transparent")
        )+
  guides(fill = guide_legend(ncol = 2))





Nitrospira_GS<-Nitrospira%>%
  filter(mfd_hab1=="Greenspaces")

num_samples <- n_distinct(Nitrospira_GS$SeqId)

lm_cladeB_P <- lm(`Nitrospira_clade_B`~`Palsa-1315*`, data = Nitrospira_GS)
summary(lm_cladeB_P)

P_GS <-Nitrospira_GS %>%
  ggplot(., aes(
    y = `Nitrospira_clade_B`,
    x= `Palsa-1315*`)) +
  geom_point(aes(fill = mfd_hab2), size = 2, alpha = 0.6, color = "black", pch = 21, stroke=.8)+
  scale_fill_manual(values=p_combine, name =  paste0("Greenspaces,\nn = ",num_samples))+
  geom_smooth(method = "lm", se = FALSE, col = "darkred", lty = 2) +
  theme_classic() +
  labs(
    y = "Nitrospira clade B amoA [RPKM]",
    x = "Nitrospiraceae Palsa-1315 nxrA [RPKM]",
    subtitle = paste0("r-squared: ", round(summary(lm_cladeB_P)$r.squared, 3), "\np-value: ",formatC(summary(lm_cladeB_P)$coefficients[2,4], format = "e", digits = 2), "\ncoefficient: ", round(summary(lm_cladeB_P)$coefficients[2,1], 2))
  )+
  theme(legend.position = c(0.6, 0.9),
        legend.background = element_rect(fill = "transparent")) +
  guides(fill = guide_legend(ncol = 2))+
  coord_cartesian(xlim = c(0, max(Nitrospira_GS$`Palsa-1315*`)), ylim = c(0,max(Nitrospira_GS$`Palsa-1315*`)))




Nitrospira_Grass<-Nitrospira%>%
  filter(mfd_hab1=="Grassland formations") 


num_samples <- n_distinct(Nitrospira_Grass$SeqId)

lm_cladeB_P <- lm(`Nitrospira_clade_B`~`Palsa-1315*`, data = Nitrospira_Grass)
summary(lm_cladeB_P)

P_Grass <-Nitrospira_Grass %>%
  ggplot(., aes(
    y = `Nitrospira_clade_B`,
    x= `Palsa-1315*`)) +
  geom_point(aes(fill = mfd_hab2), size = 2, alpha = 0.6, color = "black", pch = 21, stroke=.8)+
  scale_fill_manual(values=p_combine, name = paste0("Grassland formations,\nn = ",num_samples))+
  geom_smooth(method = "lm", se = FALSE, col = "darkred", lty = 2) +
  theme_classic() +
  labs(
    y = "Nitrospira clade B amoA [RPKM]",
    x = "Nitrospiraceae Palsa-1315 nxrA [RPKM]",
    subtitle = paste0("r-squared: ", round(summary(lm_cladeB_P)$r.squared, 3), "\np-value: ",formatC(summary(lm_cladeB_P)$coefficients[2,4], format = "e", digits = 2), "\ncoefficient: ", round(summary(lm_cladeB_P)$coefficients[2,1], 2))
  )+
  theme(legend.position = c(0.5, 0.8),
        legend.background = element_rect(fill = "transparent")) +
  guides(fill = guide_legend(ncol = 2))+
  coord_cartesian(xlim = c(0, max(Nitrospira_Grass$`Palsa-1315*`)), ylim = c(0,max(Nitrospira_Grass$`Palsa-1315*`)))




Nitrospira_Sedi<-Nitrospira%>%
  filter(mfd_hab2 %in% c("Running\nfreshwater", "Standing\nfreshwater")) 


num_samples <- n_distinct(Nitrospira_Sedi$SeqId)


lm_cladeB_P <- lm(`Nitrospira_clade_B`~`Palsa-1315*`, data = Nitrospira_Sedi)
summary(lm_cladeB_P)

P_Sedi<-Nitrospira_Sedi %>%
  ggplot(., aes(
    y = `Nitrospira_clade_B`,
    x= `Palsa-1315*`)) +
  geom_point(aes(fill = mfd_hab2), size = 2, alpha = 0.6, color = "black", pch = 21, stroke=.8)+
  scale_fill_manual(values=p_combine, name = paste0("Natural Sediment,\nn = ",num_samples))+
  geom_smooth(method = "lm", se = FALSE, col = "darkred", lty = 2) +
  theme_classic() +
  labs(
    y = "Nitrospira clade B amoA [RPKM]",
    x = "Nitrospiraceae Palsa-1315 nxrA [RPKM]",
    subtitle = paste0("r-squared: ", round(summary(lm_cladeB_P)$r.squared, 3), "\np-value: ",formatC(summary(lm_cladeB_P)$coefficients[2,4], format = "e", digits = 2), "\ncoefficient: ", round(summary(lm_cladeB_P)$coefficients[2,1], 2))
  )+
  theme(legend.position = c(0.6, 0.9),
        legend.background = element_rect(fill = "transparent")) +
  guides(fill = guide_legend(ncol = 2))+
  coord_cartesian(xlim = c(0, max(Nitrospira_Sedi$`Palsa-1315*`)), ylim = c(0,max(Nitrospira_Sedi$`Palsa-1315*`)))



P_Agri
P_BOGS
P_Forest
P_GS
P_Grass
P_Sedi

p <- plot_grid( P_Agri, P_BOGS, P_Forest, P_GS, P_Grass, P_Sedi,
                nrow = 2, ncol = 3,
                align = "hv",
                rel_widths = c(1, 1),  # Adjust the relative widths of the columns
                scale = 1  # Adjust the scale to add space between the plots
)

p

ggsave("output/Correlations_Nitrospira.png",
       p,
       height = 11,
       width = 19,
       dpi=400)



ggsave("output/Correlations_Nitrospira.svg",
       p,
       height = 11,
       width = 19)
