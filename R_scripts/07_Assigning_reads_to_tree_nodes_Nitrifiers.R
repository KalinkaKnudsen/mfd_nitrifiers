#title: "07_Assigning_reads_to_tree_nodes_Nitrifiers"
#author: "Kalinka Sand Knudsen"
#update: "2025-04-24"


#Loading packages
.libPaths(c("path_to_R_package_library", .libPaths()))
library(ggridges)
library(tidyr)
library(stringr)
library(purrr)
library(dplyr)
library(vroom)
library(ggplot2)
library(svglite)

#Set WD
setwd("path_to_working_directory")

#Setting label size for plotting globally
label_max<-37

#Importing the allocation of reads to each tree node
cNXR<-vroom("data/cNXR_positions_of_hits.tsv", delim = "\t")

#Investigating the Nitrobacter clade
Nitrobacter_new <- cNXR %>%
  filter(grepl("Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade", Tax)) %>%
  select(Sequence, Tax) %>%
  group_by(Tax) %>%
  summarise(No_of_seqs = n()) %>%
  mutate(fraction_of_group = No_of_seqs / sum(No_of_seqs))%>%
  mutate(percentage_of_group=100*fraction_of_group)%>%
  mutate(Label=label_max*sqrt(fraction_of_group))
sum(Nitrobacter_new$No_of_seqs)


####### Turning this into a coverage plot also #######
Nitrobacter_cov <- cNXR %>%
  filter(grepl("Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade", Tax)) 

gene_positions <- Nitrobacter_cov %>%
  summarise(start = min(hmm_start),
            end = max(hmm_end)) %>%
  expand(Position = seq(start, end))

reads_by_position <- gene_positions %>%
  left_join(Nitrobacter_cov %>% 
              mutate(Position = map2(hmm_start, hmm_end, seq)) %>%
              unnest(Position) %>%
              group_by(Position, Tax) %>%
              count()) %>%
  replace_na(list(n = 0)) %>%
  rename(Count=n)

all_combinations <- expand.grid(
  Position = 1:1364, ##Based on length of alignment HMM
  Tax = unique(reads_by_position$Tax)
)

# Merge with the original data to fill in missing values
result_data <- merge(all_combinations, reads_by_position, by = c("Position", "Tax"), all.x = TRUE)%>%
  filter(!is.na(Tax))%>%
  mutate(Count=if_else(is.na(Count), 0.1, Count))


Nitrobacter_new<-arrange(Nitrobacter_new, No_of_seqs)%>%
  mutate(Tax=gsub("Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; ", "", Tax),
         Tax=gsub("Root; cNarG_Nxr; ", "", Tax),
         Tax=gsub("Nitrobacter-like_NxrA; ", "", Tax),
         Tax=gsub("Pseudolabrys_JAFA_VAZQ_umbrella; JAFAXD01_VAZQ01_Pseudolabrys_cluster; ","",Tax),
         Tax=gsub("Pseudolabrys_JAFA_VAZQ_umbrella; ","",Tax),
         Tax=gsub("JAFAXD01_VAZQ01_Pseudolabrys_cluster","JAF. VAZ. Pseudo. cluster",Tax),
         Tax=gsub("Pseudolabrys_JAFA_VAZQ_umbrella","Put. NOB cluster",Tax),
         Tax=gsub("Bradyrhizobium_Nitrobacter; ","",Tax))

Nitrobacter_cov_plot<-result_data %>%
  mutate(Tax=gsub("Root; cNarG_Nxr; Nitrococcus_Nitrobacter_clade; ", "", Tax),
         Tax=gsub("Root; cNarG_Nxr; ", "", Tax),
         Tax=gsub("Nitrobacter-like_NxrA; ", "", Tax),
         Tax=gsub("Pseudolabrys_JAFA_VAZQ_umbrella; JAFAXD01_VAZQ01_Pseudolabrys_cluster; ","",Tax),
         Tax=gsub("Pseudolabrys_JAFA_VAZQ_umbrella; ","",Tax),
         Tax=gsub("JAFAXD01_VAZQ01_Pseudolabrys_cluster","JAF. VAZ. Pseudo. cluster",Tax),
         Tax=gsub("Pseudolabrys_JAFA_VAZQ_umbrella","Put. NOB cluster",Tax),
         Tax=gsub("Bradyrhizobium_Nitrobacter; ","",Tax))%>%
  filter(!is.na(Tax))%>%
  filter(!grepl("Nitrococcus|BOG-931_2", Tax))%>%
  mutate(Tax=factor(Tax, levels = Nitrobacter_new$Tax, ordered=TRUE))%>%
ggplot(., aes(x = Position, y = Tax, fill = Tax, height=Count)) +
  geom_density_ridges(scale = 10, rel_min_height = 0, stat="identity", linewidth=.5) +
  scale_x_continuous(n.breaks = 7,expand = c(0, 0)) +
  scale_fill_viridis_d(alpha=0.3, guide="none", begin=0, end=1) +
  labs(x="Cytoplasmic NxrA/NarG HMM-alignment position") +
  theme_ridges(font_size = 11, grid = FALSE, center_axis_labels = TRUE) + 
 # theme_bw() +
  theme(axis.title.y = element_blank(), 
        axis.title.x=element_text(size=12, color="grey10"),
        axis.text.x = element_text(size=11, color="grey40"),
        axis.text.y = element_text(size=12, color="grey10"),
        text = element_text(family = "Arial"))

ggsave("output/Nitrobacter_coverage.svg",
       Nitrobacter_cov_plot,
       height = 3.5,
       width = 8)



arch_AmoA<-vroom("data/arch_AmoA_positions_of_hits.tsv", delim = "\t")


### Investigating Nitrososphaeraceae
Nitrososphaeraceae_new <- arch_AmoA %>%
  filter(grepl("Root; o_Nitrososphaerales; f_Nitrososphaeraceae", Tax)) %>%
  select(Sequence, Tax) %>%
  group_by(Tax) %>%
  summarise(No_of_seqs = n()) %>%
  mutate(fraction_of_group = No_of_seqs / sum(No_of_seqs))%>%
  mutate(percentage_of_group=100*fraction_of_group)%>%
  mutate(Label=label_max*sqrt(fraction_of_group))
#sum(Nitrososphaeraceae_new$No_of_seqs)


## Turning this into a coverage plot
Nitrososphaeraceae_cov <- arch_AmoA %>%
  filter(grepl("Root; o_Nitrososphaerales; f_Nitrososphaeraceae", Tax)) %>%
  filter(!Tax %in% c("Root; o_Nitrososphaerales; f_Nitrososphaeraceae; TH5893_Nitrososphaera_cluster; JAJPHM01_Nitrososphaera_cluster"))%>%  
  mutate(Tax=gsub("Root; o_Nitrososphaerales; f_Nitrososphaeraceae; ", "", Tax),
         Tax=gsub("TH5893_Nitrososphaera_cluster; ", "", Tax),
         Tax=gsub("JAJPHM01_Nitrososphaera_cluster; ", "", Tax),
         Tax=gsub("Root; o_Nitrososphaerales; f_Nitrososphaeraceae", "Nitrososphaeraceae", Tax),
         Tax=gsub("TA21_Nitrosopolaris_cluster; ", "", Tax),
         Tax=gsub("g_","", Tax))

gene_positions <- Nitrososphaeraceae_cov %>%
  summarise(start = min(hmm_start),
            end = max(hmm_end)) %>%
  expand(Position = seq(start, end))

reads_by_position <- gene_positions %>%
  left_join(Nitrososphaeraceae_cov %>% 
              mutate(Position = map2(hmm_start, hmm_end, seq)) %>%
              unnest(Position) %>%
              group_by(Position, Tax) %>%
              count()) %>%
  replace_na(list(n = 0)) %>%
  rename(Count=n)

all_combinations <- expand.grid(
  Position = 1:216,
  Tax = unique(reads_by_position$Tax)
)

# Merge with the original data to fill in missing values
result_data <- merge(all_combinations, reads_by_position, by = c("Position", "Tax"), all.x = TRUE)%>%
  filter(!is.na(Tax))%>%
  mutate(Count=if_else(is.na(Count), 0.1, Count))


Nitrososphaeraceae_new<-arrange(Nitrososphaeraceae_new, No_of_seqs)%>%
  filter(grepl("Root; o_Nitrososphaerales; f_Nitrososphaeraceae", Tax)) %>%
  filter(!Tax %in% c("Root; o_Nitrososphaerales; f_Nitrososphaeraceae; TH5893_Nitrososphaera_cluster; JAJPHM01_Nitrososphaera_cluster"))%>%  
  mutate(Tax=gsub("Root; o_Nitrososphaerales; f_Nitrososphaeraceae; ", "", Tax),
         Tax=gsub("TH5893_Nitrososphaera_cluster; ", "", Tax),
         Tax=gsub("JAJPHM01_Nitrososphaera_cluster; ", "", Tax),
         Tax=gsub("Root; o_Nitrososphaerales; f_Nitrososphaeraceae", "Nitrososphaeraceae", Tax),
         Tax=gsub("TA21_Nitrosopolaris_cluster; ", "", Tax),
         Tax=gsub("g_","", Tax))


Nitrososphaeraceae_cov_plot<-result_data %>%
  mutate(Tax=gsub("Root; o_Nitrososphaerales; f_Nitrososphaeraceae; ", "", Tax),
         Tax=gsub("TH5893_Nitrososphaera_cluster; ", "", Tax),
         Tax=gsub("JAJPHM01_Nitrososphaera_cluster; ", "", Tax),
         Tax=gsub("Root; o_Nitrososphaerales; f_Nitrososphaeraceae", "Nitrososphaeraceae", Tax),
         Tax=gsub("TA21_Nitrosopolaris_cluster; ", "", Tax),
         Tax=gsub("g_","", Tax))%>%
  filter(!is.na(Tax))%>%
  filter(!grepl("TA21_3|TH5893_2|TH5896_2", Tax))%>%
  mutate(Tax=factor(Tax, levels = Nitrososphaeraceae_new$Tax, ordered=TRUE))%>%
  #mutate(Tax=gsub("_"," ", Tax))%>%
  ggplot(., aes(x = Position, y = Tax, fill = Tax, height=Count)) +
  geom_density_ridges(scale = 6, rel_min_height = 0, stat="identity", linewidth=.5) +
  scale_x_continuous(n.breaks = 10,expand = c(0, 0)) +
  scale_fill_viridis_d(alpha=0.3, guide="none", begin=1, end=0.4) +
  labs(x="Archaeal AmoA HMM-alignment position") +
  theme_ridges(font_size = 10, grid = FALSE, center_axis_labels = TRUE) + 
  # theme_bw() +
  theme(axis.title.y = element_blank(), 
        axis.title.x=element_text(size=12, color="grey10"),
        axis.text.x = element_text(size=12, color="grey40"),
        axis.text.y = element_text(size=10, color="grey10"),
        text = element_text(family = "Arial"))


Nitrososphaeraceae_cov_plot


ggsave("output/Nitrososphaeraceae_coverage.svg",
       Nitrososphaeraceae_cov_plot,
       height = 2.8,
       width = 7)






