#title: "05_linear_models_Nitrospira"
#author: "Kalinka Sand Knudsen"
#update: "2024-07-31"


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
OTU_filtered_long<-readRDS("OTU_filtered_long.rds") # from normalising_and_aggregating_nitrifier_genes


#Ordering the tax:

tax<-OTU_filtered_long%>%
  select(Tax)%>%
  distinct()


## Curating the gene-taxonomy
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
         Tax_curated = if_else(Tax=="Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; Nitrobacter-like_NxrA; Bradyrhizobium_Nitrobacter", "Nitrobacter_like_nxrA; Bradyrhizobium_Nitrobacter_cluster", Tax_curated),
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

tax2 <- readxl::read_excel("tax_order.xlsx")

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
  mutate(mfd_hab2 = if_else(is.na(mfd_hab2) & complex == "Soil\nNatural\nForests", paste0(mfd_hab1, " - no MFDO2"), mfd_hab2)) %>%
  mutate(mfd_hab2 = if_else(mfd_hab2=="NA" & complex == "Soil\nNatural\nForests", paste0(mfd_hab1, " - no MFDO2"), mfd_hab2)) %>%
  filter(!is.na(mfd_hab2))%>%
  filter(!mfd_hab2=="NA")%>%
  filter(!mfd_hab2=="Mire")%>%
  filter(!is.na(Tax_curated))


OTU_wide<-OTU_filtered_long%>%
  pivot_wider(names_from = Tax_curated, values_from = RPKM, id_cols = c(SeqId))%>%
  as.data.frame(.)%>%
  filter(rowSums(select(., -1)) > 0)
OTU_sums<-OTU_wide[,-1]
row.names(OTU_sums)<-OTU_wide$SeqId

#Trying to do correlation between Nitrospira_clade_B and Palsa-1315
p_combine<-readRDS("palette_mfd_hab2.rds")

Nitrospira<-OTU_filtered_long%>%
  filter(Tax_curated %in% c("Nitrospira_amoA; Nitrospira_clade_B", "Nitrospirales; Palsa_1315*"))%>%
  select(SeqId, Tax_curated, RPKM, mfd_hab1, mfd_hab2)%>%
  pivot_wider(names_from = Tax_curated, values_from = RPKM)




#### Subsetting to selected habitats for investigation 
Nitrospira_BOG<-Nitrospira%>%
  filter(mfd_hab1=="Bogs, mires and fens")

#### Removing anything more than 3 sigma away ####
x<-sd(Nitrospira_BOG$`Nitrospira_amoA; Nitrospira_clade_B`)
y<-sd(Nitrospira_BOG$`Nitrospirales; Palsa_1315*`)

Nitrospira_BOG<-Nitrospira_BOG%>%
  filter(`Nitrospira_amoA; Nitrospira_clade_B`< mean(`Nitrospira_amoA; Nitrospira_clade_B`) +3*x)%>%
  filter(`Nitrospirales; Palsa_1315*`< mean(`Nitrospirales; Palsa_1315*`)+3*y)


### Making linear model
lm_cladeB_P <- lm(`Nitrospira_amoA; Nitrospira_clade_B`~`Nitrospirales; Palsa_1315*`, data = Nitrospira_BOG)
summary(lm_cladeB_P)

num_samples <- n_distinct(Nitrospira_BOG$SeqId)

P_BOGS <-Nitrospira_BOG %>%
  filter(mfd_hab1=="Bogs, mires and fens") %>%
  ggplot(., aes(
    y = `Nitrospira_amoA; Nitrospira_clade_B`,
    x= `Nitrospirales; Palsa_1315*`, color = mfd_hab2)) +
  scale_color_manual(values=p_combine, name = paste0("Bogs, Mires and Fens,\nn = ",num_samples))+
  geom_point(fill="white", alpha=0.9, pch=1, stroke=1) +
  geom_smooth(method = "lm", se = FALSE, col = "darkred", lty = 2) +
  theme_classic() +
  labs(
    title = "Nitrospira clade B amoA vs Palsa-1315 nxrA group",
    y = "Nitrospira clade B amoA [RPKM]",
    x = "Nitrospiraceae Palsa 1315 [RPKM]",
    subtitle = paste0("r-squared: ", round(summary(lm_cladeB_P)$r.squared, 3), "\np-value: ",formatC(summary(lm_cladeB_P)$coefficients[2,4], format = "e", digits = 2), "\ncoefficient: ", round(summary(lm_cladeB_P)$coefficients[2,1], 2))
  )+
  coord_cartesian(xlim = c(0, 6), ylim = c(0, 6))+
  theme(legend.position = c(0.6, 0.9),
        legend.background = element_rect(fill = "transparent"))+
  guides(color = guide_legend(ncol = 1))

P_BOGS


#### Doing the same for field samples
Nitrospira_Agri<-Nitrospira%>%
  filter(mfd_hab1=="Fields")

boxplot(Nitrospira_Agri$`Nitrospira_amoA; Nitrospira_clade_B`)

x<-sd(Nitrospira_Agri$`Nitrospira_amoA; Nitrospira_clade_B`)
y<-sd(Nitrospira_Agri$`Nitrospirales; Palsa_1315*`)

Nitrospira_Agri<-Nitrospira_Agri%>%
  filter(`Nitrospira_amoA; Nitrospira_clade_B`< mean(`Nitrospira_amoA; Nitrospira_clade_B`) +3*x)%>%
  filter(`Nitrospirales; Palsa_1315*`< mean(`Nitrospirales; Palsa_1315*`)+3*y)

num_samples <- n_distinct(Nitrospira_Agri$SeqId)

lm_cladeB_P <- lm(`Nitrospira_amoA; Nitrospira_clade_B`~`Nitrospirales; Palsa_1315*`, data = Nitrospira_Agri)
summary(lm_cladeB_P)

num_samples <- n_distinct(lm_cladeB_P$model)

P_Agri <-Nitrospira_Agri %>%
  ggplot(., aes(
    y = `Nitrospira_amoA; Nitrospira_clade_B`,
    x= `Nitrospirales; Palsa_1315*`, color = mfd_hab2)) +
  scale_color_manual(values=p_combine, name = paste0("Agriculture,\nn = ",num_samples))+
  geom_point(fill="white", alpha=0.9, pch=1, stroke=1) +
  geom_smooth(method = "lm", se = FALSE, col = "darkred", lty = 2) +
  theme_classic() +
  labs(
    title = "Nitrospira clade B amoA vs Palsa-1315 nxrA group",
    y = "Nitrospira clade B amoA [RPKM]",
    x = "Nitrospiraceae Palsa 1315 [RPKM]",
    subtitle = paste0("r-squared: ", round(summary(lm_cladeB_P)$r.squared, 3), "\np-value: ",formatC(summary(lm_cladeB_P)$coefficients[2,4], format = "e", digits = 2), "\ncoefficient: ", round(summary(lm_cladeB_P)$coefficients[2,1], 2))
  )+
  coord_cartesian(xlim = c(0, 3), ylim = c(0, 3))+
  theme(legend.position = c(0.6, 0.9),
        legend.background = element_rect(fill = "transparent"))+
  guides(color = guide_legend(ncol = 2))



#### Doing the same for Forest samples
Nitrospira_Forest<-Nitrospira%>%
  filter(mfd_hab1=="Forests")

x<-sd(Nitrospira_Forest$`Nitrospira_amoA; Nitrospira_clade_B`)
y<-sd(Nitrospira_Forest$`Nitrospirales; Palsa_1315*`)

Nitrospira_Forest<-Nitrospira_Forest%>%
  filter(`Nitrospira_amoA; Nitrospira_clade_B`< mean(`Nitrospira_amoA; Nitrospira_clade_B`) +3*x)%>%
  filter(`Nitrospirales; Palsa_1315*`< mean(`Nitrospirales; Palsa_1315*`)+3*y)

num_samples <- n_distinct(Nitrospira_Forest$SeqId)

lm_cladeB_P <- lm(`Nitrospira_amoA; Nitrospira_clade_B`~`Nitrospirales; Palsa_1315*`, data = Nitrospira_Forest)
summary(lm_cladeB_P)

P_Forest <-Nitrospira_Forest %>%
  ggplot(., aes(
    y = `Nitrospira_amoA; Nitrospira_clade_B`,
    x= `Nitrospirales; Palsa_1315*`, color = mfd_hab2)) +
  scale_color_manual(values=p_combine, name = paste0("Forests,\nn = ",num_samples))+
  geom_point(fill="white", alpha=0.9, pch=1, stroke=1) +
  geom_smooth(method = "lm", se = FALSE, col = "darkred", lty = 2) +
  theme_classic() +
  labs(
    title = "Nitrospira clade B amoA vs Palsa-1315 nxrA group",
    y = "Nitrospira clade B amoA [RPKM]",
    x = "Nitrospiraceae Palsa 1315 [RPKM]",
    subtitle = paste0("r-squared: ", round(summary(lm_cladeB_P)$r.squared, 3), "\np-value: ",formatC(summary(lm_cladeB_P)$coefficients[2,4], format = "e", digits = 2), "\ncoefficient: ", round(summary(lm_cladeB_P)$coefficients[2,1], 2))
  )+
  coord_cartesian(xlim = c(0, 5), ylim = c(0, 5))+
  theme(legend.position = c(0.7, 0.9),
        legend.background = element_rect(fill = "transparent"))+
  guides(color = guide_legend(ncol = 2))


#### Doing the same for Greenspaces samples
Nitrospira_GS<-Nitrospira%>%
  filter(mfd_hab1=="Greenspaces")

x<-sd(Nitrospira_GS$`Nitrospira_amoA; Nitrospira_clade_B`)
y<-sd(Nitrospira_GS$`Nitrospirales; Palsa_1315*`)


Nitrospira_GS<-Nitrospira_GS%>%
  filter(`Nitrospira_amoA; Nitrospira_clade_B`< mean(`Nitrospira_amoA; Nitrospira_clade_B`) +3*x)%>%
  filter(`Nitrospirales; Palsa_1315*`< mean(`Nitrospirales; Palsa_1315*`)+3*y)

num_samples <- n_distinct(Nitrospira_GS$SeqId)

lm_cladeB_P <- lm(`Nitrospira_amoA; Nitrospira_clade_B`~`Nitrospirales; Palsa_1315*`, data = Nitrospira_GS)
summary(lm_cladeB_P)

P_GS <-Nitrospira_GS %>%
  ggplot(., aes(
    y = `Nitrospira_amoA; Nitrospira_clade_B`,
    x= `Nitrospirales; Palsa_1315*`, color = mfd_hab2)) +
  scale_color_manual(values=p_combine, name = paste0("Greenspaces,\nn = ",num_samples))+
  geom_point(fill="white", alpha=0.9, pch=1, stroke=1) +
  geom_smooth(method = "lm", se = FALSE, col = "darkred", lty = 2) +
  theme_classic() +
  labs(
    title = "Nitrospira clade B amoA vs Palsa-1315 nxrA group",
    y = "Nitrospira clade B amoA [RPKM]",
    x = "Nitrospiraceae Palsa 1315 [RPKM]",
    subtitle = paste0("r-squared: ", round(summary(lm_cladeB_P)$r.squared, 3), "\np-value: ",formatC(summary(lm_cladeB_P)$coefficients[2,4], format = "e", digits = 2), "\ncoefficient: ", round(summary(lm_cladeB_P)$coefficients[2,1], 2))
  )+
  coord_cartesian(xlim = c(0, 4), ylim = c(0, 4))+
  theme(legend.position = c(0.6, 0.9),
        legend.background = element_rect(fill = "transparent"))




#### Doing the same for Grassland samples
Nitrospira_Grass<-Nitrospira%>%
  filter(mfd_hab1=="Grassland formations") #%>%   filter(`Nitrospirales; Palsa_1315*`<30)


x<-sd(Nitrospira_Grass$`Nitrospira_amoA; Nitrospira_clade_B`)
y<-sd(Nitrospira_Grass$`Nitrospirales; Palsa_1315*`)


Nitrospira_Grass<-Nitrospira_Grass%>%
  filter(`Nitrospira_amoA; Nitrospira_clade_B`< mean(`Nitrospira_amoA; Nitrospira_clade_B`) +3*x)%>%
  filter(`Nitrospirales; Palsa_1315*`< mean(`Nitrospirales; Palsa_1315*`)+3*y)

num_samples <- n_distinct(Nitrospira_Grass$SeqId)

lm_cladeB_P <- lm(`Nitrospira_amoA; Nitrospira_clade_B`~`Nitrospirales; Palsa_1315*`, data = Nitrospira_Grass)
summary(lm_cladeB_P)

P_Grass <-Nitrospira_Grass %>%
  ggplot(., aes(
    y = `Nitrospira_amoA; Nitrospira_clade_B`,
    x= `Nitrospirales; Palsa_1315*`, color = mfd_hab2)) +
  scale_color_manual(values=p_combine, name = paste0("Grassland formations,\nn = ",num_samples))+
  geom_point(fill="white", alpha=0.9, pch=1, stroke=1) +
  geom_smooth(method = "lm", se = FALSE, col = "darkred", lty = 2) +
  theme_classic() +
  labs(
    title = "Nitrospira clade B amoA vs Palsa-1315 nxrA group",
    y = "Nitrospira clade B amoA [RPKM]",
    x = "Nitrospiraceae Palsa 1315 [RPKM]",
    subtitle = paste0("r-squared: ", round(summary(lm_cladeB_P)$r.squared, 3), "\np-value: ",formatC(summary(lm_cladeB_P)$coefficients[2,4], format = "e", digits = 2), "\ncoefficient: ", round(summary(lm_cladeB_P)$coefficients[2,1], 2))
  )+
  coord_cartesian(xlim = c(0, 6), ylim = c(0, 6))+
  theme(legend.position = c(0.6, 0.9),
        legend.background = element_rect(fill = "transparent"))


#### Doing the same for Sediment samples
Nitrospira_Sedi<-Nitrospira%>%
  filter(mfd_hab2 %in% c("Running freshwater", "Standing freshwater")) #%>%
# filter(`Nitrospirales; Palsa_1315*`<30)



x<-sd(Nitrospira_Sedi$`Nitrospira_amoA; Nitrospira_clade_B`)
y<-sd(Nitrospira_Sedi$`Nitrospirales; Palsa_1315*`)


Nitrospira_Sedi<-Nitrospira_Sedi%>%
  filter(`Nitrospira_amoA; Nitrospira_clade_B`< mean(`Nitrospira_amoA; Nitrospira_clade_B`) +3*x)%>%
  filter(`Nitrospirales; Palsa_1315*`< mean(`Nitrospirales; Palsa_1315*`)+3*y)


num_samples <- n_distinct(Nitrospira_Sedi$SeqId)


lm_cladeB_P <- lm(`Nitrospira_amoA; Nitrospira_clade_B`~`Nitrospirales; Palsa_1315*`, data = Nitrospira_Sedi)
summary(lm_cladeB_P)

P_Sedi<-Nitrospira_Sedi %>%
  ggplot(., aes(
    y = `Nitrospira_amoA; Nitrospira_clade_B`,
    x= `Nitrospirales; Palsa_1315*`, color = mfd_hab2)) +
  scale_color_manual(values=p_combine, name = paste0("Natural Sediment,\nn = ",num_samples))+
  geom_point(fill="white", alpha=0.9, pch=1, stroke=1) +
  geom_smooth(method = "lm", se = FALSE, col = "darkred", lty = 2) +
  theme_classic() +
  labs(
    title = "Nitrospira clade B amoA vs Palsa-1315 nxrA group",
    y = "Nitrospira clade B amoA [RPKM]",
    x = "Nitrospiraceae Palsa 1315 [RPKM]",
    subtitle = paste0("r-squared: ", round(summary(lm_cladeB_P)$r.squared, 3), "\np-value: ",formatC(summary(lm_cladeB_P)$coefficients[2,4], format = "e", digits = 2), "\ncoefficient: ", round(summary(lm_cladeB_P)$coefficients[2,1], 2))
  )+
  coord_cartesian(xlim = c(0, 5), ylim = c(0, 5))+
  theme(legend.position = c(0.6, 0.9),
        legend.background = element_rect(fill = "transparent"))


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
                scale = 0.98  # Adjust the scale to add space between the plots
)

p

ggsave("Correlations_3sd.png",
       p,
       height = 11,
       width = 19)



ggsave("Correlations_3sd.svg",
       p,
       height = 11,
       width = 19)




