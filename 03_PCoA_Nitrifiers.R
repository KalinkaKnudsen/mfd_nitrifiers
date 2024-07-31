#title: "03_PCoA_Nitrifiers"
#author: "Kalinka Sand Knudsen"
#update: "2024-07-30"


#Loading packages
.libPaths(c("path_to_R_package_library", .libPaths()))
library(tidyverse)
library(vegan)
library(ape)
library(patchwork)
library(vroom)

#Set WD
setwd("path_to_working_directory")


# Load data
p_combine<-readRDS("palette_mfd_hab2.rds") #Color palette


OTU<-readRDS("OTU_wide.rds") ###Generated in previous script 02_Nitrifier_Heatmaps

### Transforming the data before hellinger transformation
OTU[, -1] <- apply(OTU[, -1], 2, function(x) as.numeric(as.character(x)))

transposed_OTU <- as.data.frame(t(OTU))
colnames(transposed_OTU) <- transposed_OTU[1, ]
transposed_OTU <- transposed_OTU[-1, ]
transposed_OTU <- transposed_OTU%>%
  mutate(OTU=rownames(transposed_OTU))%>%
  relocate(OTU)
rownames(transposed_OTU) <- NULL

transposed_OTU[, -1] <- apply(transposed_OTU[, -1], 2, function(x) as.numeric(as.character(x)))


data2<-transposed_OTU%>%
  select(where(is.numeric))
x2<-transposed_OTU[, which(colSums(data2) != 0)]
data <- transposed_OTU %>%
  select(intersect(names(.), names(x2)), OTU)

dist <- data %>%
  select(where(is.numeric), OTU) %>%
  column_to_rownames(var = "OTU") %>%
  t() %>%
  decostand(., method = "hellinger") %>%
  parallelDist::parDist(., method = "bray", threads = 4) %>% #### Calculating Bray-Curtis distance based on hellinger transformed data
  as.matrix() %>%
  data.frame()%>%rename_with(~ gsub(".", "-", .x, fixed = TRUE))

### PCoA
PCOA <- pcoa(dist)

# plot the eigenvalues and interpret
barplot(PCOA$values$Relative_eig[1:10])

sum(as.vector(PCOA$value$Relative_eig)[1])
sum(as.vector(PCOA$value$Relative_eig)[2])

sum(as.vector(PCOA$value$Relative_eig)[1:2])



##### Linking the samples ####

linkage <- vroom("2023-10-11_samples_minimal_metadata_collapsed.csv", delim = ",") %>%
  mutate(flat_name=gsub(".fastq.gz","", flat_name))%>%
  rename(SeqId = flat_name) %>%
  relocate(SeqId)%>%
  select(SeqId, fieldsample_barcode)


## Loading metadata ####
metadata.sub <- readxl::read_excel('2024-02-13_mfd_db.xlsx') %>%
  left_join(linkage)%>%
  filter(SeqId %in% linkage$SeqId)%>%
  relocate(SeqId)%>%
  filter(SeqId %in% colnames(dist)) %>%
  select(SeqId:mfd_hab3) %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", "))



#### Obtaining the PCoA scores and filtering the data (same procedure as the heatmaps)
PCOA.scores <- PCOA$vectors %>%
  as.data.frame() %>%
  select(1:10) %>%
  rename_with(., ~str_replace(., "Axis.", "PCO")) %>%
  cbind(SeqId = colnames(dist)) %>%
  relocate(SeqId, .before = "PCO1") %>%
  filter(SeqId %in% metadata.sub$SeqId)%>%
  left_join(metadata.sub) %>%
   mutate(across(mfd_hab1, ~str_replace(., "Sclerophyllous scrub", "Temperate heath and scrub")),
         across(mfd_hab2, ~str_replace(., "Scrub", "Sclerophyllous scrub"))) %>%
  mutate(mfd_hab2=gsub("\\s*\\(non-habitat type\\)\\s*", "", mfd_hab2))%>%
  mutate(complex_long = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", ")) %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = "\n")) %>%
  group_by(complex) %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", ")) %>%
  group_by(complex) %>%
  mutate(complex_size = n()) %>%
  mutate(label = str_c(complex, "\nn = ", complex_size, sep = "")) %>%
  ungroup()%>%
  filter(!complex %in% c("Water, Subterranean, Freshwater",
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
  mutate(mfd_hab2 = if_else(complex == "Soil, Natural, Forests", mfd_hab3, mfd_hab2))%>%
  mutate(mfd_hab2 = if_else(grepl("Beech", mfd_hab2), "Beech", mfd_hab2)) %>%
  mutate(mfd_hab2 = if_else(grepl("Birch", mfd_hab2), "Birch", mfd_hab2)) %>%
  mutate(mfd_hab2 = if_else(grepl("Oak", mfd_hab2), "Oak", mfd_hab2)) %>%
  mutate(mfd_hab2 = if_else(complex == "Sediment, Urban, Freshwater", paste0(mfd_hab3), paste0(mfd_hab2)))%>%
  filter(!mfd_hab2 %in% c("Enclosed water, Dried", "Birch", "Pine", "Mire (non-habitat type)"))%>%
  mutate(mfd_hab2 = if_else(is.na(mfd_hab2) & mfd_hab2 == "Forests", paste0(mfd_hab1, " - no MFDO2"), mfd_hab2)) %>%
  mutate(mfd_hab2 = if_else(mfd_hab2=="NA" & mfd_hab2 == "Forests", paste0(mfd_hab1, " - no MFDO2"), mfd_hab2)) %>%
  filter(!is.na(mfd_hab2))%>%
  filter(!mfd_hab2=="NA")%>%
  filter(!mfd_hab2=="Mire")



### Plotting the PCoA
p.pcoa.1v2 <- PCOA.scores %>%
  ggplot() + 
  geom_point(
    aes(x = PCO1, y = PCO2, fill = complex), 
    size = 4, alpha = 1, color = "black", pch = 21) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(colour = "black")) +
  xlab("PCO1 - 23%") +
  ylab("PCO2 - 19%")

#p.pcoa.1v2

#png(file = 'output/plot_PCOA.png',
 #   width = 1900,
  #  height = 1200)
#p.pcoa.1v2
#dev.off()


sum(as.vector(PCOA$value$Relative_eig)[1])
sum(as.vector(PCOA$value$Relative_eig)[2])


#### Plotting the PCoA and faceting by habitats

p.pcoa.1v2.all <- ggplot() +
  geom_point(data = PCOA.scores[-27], aes(x = PCO1, y = PCO2), 
             size = 2, alpha = 0.5, color = "black", fill = "white", pch = 21) +
  geom_point(data = PCOA.scores, aes(x = PCO1, y = PCO2, fill = mfd_hab2), 
             size = 2, alpha = 1, color = "black", pch = 21) +
  theme_minimal() +
  #scale_x_continuous(limits = c(-0.6, 0.3)) +
  #scale_y_continuous(limits = c(-0.3, 0.6)) +
  scale_fill_manual(values = p_combine) +
  guides(fill = guide_legend(override.aes = list(shape = 22, size = 5, alpha = 1))) +
  theme(legend.position = "none",
        aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 16),
        plot.margin = margin(0,0,0,0),
        text = element_text(family = "Arial")) +
  xlab(paste0("PCO1 - ", round((sum(as.vector(PCOA$value$Relative_eig)[1])*100),0), "%")) +
  ylab(paste0("PCO2 - ", round((sum(as.vector(PCOA$value$Relative_eig)[2])*100),0), "%")) +
  facet_wrap(~label, ncol = 4) +
 # ggtitle("Nitrifier amoA and nxr", subtitle="All samples") +  # Add an overall header
  theme(axis.title.x = element_text(size = 14),  # Adjust x-axis label size
        axis.title.y = element_text(size = 14))  # Adjust y-axis label size


ggsave("PCoA.png", limitsize = FALSE,
       p.pcoa.1v2.all,
       height = 10,
       width = 20) 

