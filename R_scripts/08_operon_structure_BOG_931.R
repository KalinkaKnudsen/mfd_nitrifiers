#title: "08_operon_structure_BOG_931",
#author: "Kalinka Sand Knudsen"
#update: "2025-05-05"

#Loading packages
.libPaths(c("path_to_R_package_library", .libPaths()))

library(readxl)
library(tidyverse)
library(vroom)
library(gggenes)
library(svglite)

setwd("path_to_working_directory")


KEGG <- vroom("data/annotations_nitrifiers_in_paper.tsv", delim="\t", col_select=c("fasta", "ko_id", "kegg_genes_id", "kegg_hit", "pfam_hits", "scaffold", "gene_position", "start_position", "end_position", "strandedness")) %>%
  separate(kegg_hit, into = c("kegg_genes_id_2", "gene"), sep = '\\)  ', remove = T)%>%
  select(!kegg_genes_id_2)


KOs <- c("K01601","K01602", "K03775")

d <- KEGG%>%
  filter(grepl("nitrate|molybdopterin-dependent oxidoreductase|putative oxidoreductase molybdopterin-binding subunit|cytochrome c|cytochrome C|FKBP|transposase|ribulose|cbb|integrase|nitrite transporter", gene) |
         grepl("K01601|K01602", ko_id)) # %>% filter(gene =="nitrate oxidoreductase alpha subunit")


length<-vroom("data/Nitrobacter_like_nxr_lengths_sort.tsv", delim="\t", col_names = c("contig","length"))
tax<-vroom("data/Nitrobacter_like_nxr_mags_combined.tsv", delim="\t", col_names = c("fasta", "contig","gene","quality","GTDB_tax"))%>%
  left_join(.,length)%>%
  filter(quality=="HQ")%>%
  separate(contig, into = c("contig_text","scaffold","start_position", "frame", "gene_no"), sep = '_', remove=FALSE)%>%
  mutate(scaffold=paste0("contig_", scaffold))%>%
  mutate(start_position=as.numeric(start_position), length=as.numeric(length))%>%
  filter(length>200)%>%
  mutate(end_position=start_position+3*length)%>%
  select(-c(contig_text, contig, frame, length, quality, gene_no))%>%
  mutate(gene=gsub("Root; Cytoplasmic_nxr_narG; ", "", gene))%>%
  mutate(strandedness="1")

d_filter<-d%>%
  filter(scaffold %in% unique(tax$scaffold))%>% #### Modify here to filter on contigs 
  select(fasta, gene, scaffold, start_position, end_position, strandedness) %>%
  merge(tax[, c("fasta", "GTDB_tax")], by = "fasta")%>%
  distinct()

d2<-rbind(d_filter, tax)%>%
  mutate(gene=if_else(grepl("ABC", gene), paste0("put. ABC-nitrate transp."), paste0(gene)))%>%
  mutate(gene=gsub("respiratory nitrate reductase", "nar", gene))%>%
  mutate(gene=gsub("cbbL; Ribulose bisphosphate carboxylase large chain \\(RuBisCO large subunit\\)", "cbbL", gene)) %>%
  mutate(gene=gsub("ribulose bisphosphate carboxylase small subunit", "cbbS", gene))%>%
  mutate(gene=gsub("ribulose-bisphosphate carboxylase large subunit", "cbbL", gene))%>%
  mutate(gene=gsub("form I ribulose bisphosphate carboxylase large subunit", "cbbL", gene))%>%
  mutate(gene=gsub("ribulose 1,5-bisphosphate carboxylase", "cbb", gene))%>%
  mutate(gene=gsub("Ribulose-bisphosphate carboxylase", "cbb", gene))%>%
  mutate(gene=gsub("Ribulose bisphosphate carboxylase large chain", "cbbL", gene))%>%
  mutate(gene=gsub("cbbS; Ribulose bisphosphate carboxylase small chain \\(RuBisCO small subunit\\)", "cbbS", gene))%>%
  mutate(gene=gsub("rbcL; cbbL", "cbbL", gene))%>%
  mutate(gene_short=if_else(grepl("Nitrobacter", gene), paste0("nxrA"), gene))%>%
  mutate(gene_short=if_else(grepl("nar beta", gene_short), paste0("nxrB/narH"), gene_short))%>%
  mutate(gene_short=if_else(grepl("NarJ", gene_short), paste0("narJ"), gene_short))%>%
  mutate(gene_short=if_else(grepl("nar gamma", gene_short), paste0("narI"), gene_short))%>%
  mutate(gene_short=if_else(grepl("alpha subunit", gene_short), paste0("nar a"), gene_short))%>%
  mutate(gene_short=if_else(grepl("cytochrome c, class I", gene_short), paste0("cyt. c class I"), gene_short))%>%
  mutate(gene_short=if_else(grepl("NarK/NasA", gene_short), paste0("NarK/NasA"), gene_short))%>%
  mutate(gene_short=if_else(grepl("FKBP", gene_short), paste0("nxrX"), gene_short))%>%
  mutate(gene_short=if_else(grepl("transposase", gene_short), paste0("transposase"), gene_short))%>%
  mutate(gene_short=if_else(grepl("integrase", gene_short), paste0("integrase"), gene_short))%>%
  mutate(gene_short=if_else(grepl("cbbF", gene_short), paste0("cbbF"), gene_short))%>%
  mutate(gene_short=if_else(grepl("cbbX", gene_short), paste0("cbbX"), gene_short))%>%
  mutate(gene_short=if_else(grepl("aa3-type", gene_short), paste0("cyt. c subunit IV"), gene_short))%>%
  mutate(gene_short=if_else(grepl("CtaG", gene_short), paste0("cyt. c CtaG"), gene_short))%>%
  mutate(gene_short=if_else(grepl("rpe", gene_short), paste0("rpe"), gene_short))%>%
  mutate(gene_short=if_else(grepl("cytochrome C oxidase subunit II|cytochrome c oxidase subunit II", gene_short), paste0("cyt. c subunit II"), gene_short))%>%
  mutate(gene_short=if_else(grepl("cytochrome c oxidase, subunit III", gene_short), paste0("cyt. c subunit III"), gene_short))%>%
  mutate(gene_short=if_else(grepl("cytochrome c oxidase subunit 3", gene_short), paste0("cyt. c subunit III"), gene_short))%>%
  mutate(gene_short=if_else(grepl("cytochrome c family protein", gene_short), paste0("cyt. c"), gene_short))%>%
  mutate(gene_short=if_else(grepl("ctaD", gene_short), paste0("cyt. c subunit I"), gene_short))%>%
  mutate(gene_short=if_else(grepl("Formate/nitrite transporter", gene_short), paste0("formate/nitrite transporter"), gene_short))



custom_colors <-c("#679b60","#5e4fa2", "turquoise4", "darkred","#d97512","#c46ca1" , "#b3943c","darkblue","red","yellow" , "purple", "coral")


p_BOG<- d2 %>%
  mutate(direction=if_else(strandedness=="1", T, F))%>%
  filter(!gene_short=="nxrA")%>% 
  mutate(gene_short=gsub("nar a", "nxrA", gene_short))%>%
  filter(if_else(fasta == "MFD09603.bin.c.1", 
                 start_position > 1680000 & end_position < 1725000, 
                 start_position > 0))%>%
  filter(if_else(fasta == "MFD00991.bin.1.40", 
                 start_position > 820062 & end_position < 890062,
                 start_position > 0))%>% 
  filter(if_else(fasta == "MFD01872.bin.2.154", 
                 start_position < 272000 & end_position > 210000, 
                 start_position > 0))%>% 
  filter(if_else(fasta == "MFD03399.bin.2.43", 
                 start_position > 300002 , 
                 start_position > 0))%>%
  filter(if_else(fasta == "MFD05705.bin.1.139", 
                 start_position < 65000 & start_position > 30000, 
                 start_position > 0))%>%
 mutate(start = (start_position + end_position)/2)  %>%
  ggplot(., aes(xmin = start_position, xmax = end_position, y = fasta, fill = gene_short)) +
  geom_gene_arrow(aes(forward=direction)) +
  facet_wrap(~ fasta, scales = "free", ncol=1) +
  theme_genes()+
  theme(axis.title = element_blank(),
       # axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right") +
  #scale_fill_manual(values = custom_colors) +
  geom_text(aes(x=start_position, label = gene_short), 
            nudge_y = 0.3, angle = 8,
            hjust= 0,
            size= 3, family="Arial")+
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size=10, color="grey40"),
        axis.text.y = element_text(size=10, color="grey10"))+
  labs(title="HQ BOG-931 MAGs")

p_BOG


###################################### From this plot, it is evident that MFD05705.bin.1.139 should be flipped. I will do this below ##########################

d3<-d2 %>%
  filter(if_else(fasta == "MFD09603.bin.c.1", 
                 start_position > 1680000 & end_position < 1725000, 
                 start_position > 0))%>%
  filter(if_else(fasta == "MFD00991.bin.1.40", 
                 start_position > 820062 & end_position < 890062,
                 start_position > 0))%>% 
  filter(if_else(fasta == "MFD01872.bin.2.154", 
                 start_position < 272000 & end_position > 210000, 
                 start_position > 0))%>% 
  filter(if_else(fasta == "MFD03399.bin.2.43", 
                 start_position > 300002 , 
                 start_position > 0))%>%
  filter(if_else(fasta == "MFD05705.bin.1.139", 
                 start_position < 65000 & start_position > 30000, 
                 start_position > 0))

max_nxrA_value <- max(d3$end_position[d3$fasta=="MFD05705.bin.1.139"])
  
BOG_genes_core_operon<-d3 %>%
  filter(grepl("BOG-931", GTDB_tax))%>%
  mutate(length=end_position-start_position)%>%
  mutate(start_position2=if_else(fasta=="MFD05705.bin.1.139", -start_position+max_nxrA_value, start_position))%>%
  mutate(end_position2=if_else(fasta=="MFD05705.bin.1.139", start_position2-length, end_position))%>%
  filter(!gene_short=="nxrA")%>%
  mutate(gene_short=gsub("nar a", "nxrA", gene_short))



BOG_genes_core_operon<-d3 %>%
  filter(grepl("BOG-931", GTDB_tax))%>%
  mutate(length=end_position-start_position)%>%
  mutate(end_position2=if_else(fasta=="MFD05705.bin.1.139", -start_position+max_nxrA_value, end_position))%>%
  mutate(start_position2=if_else(fasta=="MFD05705.bin.1.139", end_position2-length, start_position))%>%
  mutate(strandedness=as.numeric(strandedness))%>%
  mutate(strandedness=if_else(fasta=="MFD05705.bin.1.139", -strandedness, strandedness))%>%
  filter(!gene_short=="nxrA")%>%
  mutate(gene_short=gsub("nar a", "nxrA", gene_short))


p_BOG2<- BOG_genes_core_operon %>%
  filter(grepl("BOG-931", GTDB_tax))%>%
  mutate(direction=if_else(strandedness=="1", T, F))%>%
  filter(if_else(fasta == "MFD09603.bin.c.1", 
                 start_position > 1680000 & end_position < 1725000, 
                 start_position > 0))%>%
  filter(if_else(fasta == "MFD00991.bin.1.40", 
                 start_position > 820062 & end_position < 890062,
                 start_position > 0))%>% 
  filter(if_else(fasta == "MFD01872.bin.2.154", 
                 start_position < 272000 & end_position > 210000, 
                 start_position > 0))%>% 
  filter(if_else(fasta == "MFD03399.bin.2.43", 
                 start_position > 300002 , 
                 start_position > 0))%>%
  filter(if_else(fasta == "MFD05705.bin.1.139", 
                 start_position < 65000 & start_position > 30000, 
                 start_position > 0))%>%
  mutate(start = (start_position2 + end_position2)/2)  %>%
  ggplot(., aes(xmin = start_position2, xmax = end_position2, y = fasta, fill = gene_short)) +
  geom_gene_arrow(aes(forward=direction)) +
  facet_wrap(~ fasta, scales = "free", ncol=1) +
  theme_genes()+
  theme(axis.title = element_blank(),
        # axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right") +
  #scale_fill_manual(values = custom_colors) +
  geom_text(aes(x=start_position2, label = gene_short), 
            nudge_y = 0.3, angle = 8,
            hjust= 0,
            size= 3, family="Arial")+
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size=10, color="grey40"),
        axis.text.y = element_text(size=10, color="grey10"))+
  labs(title="HQ BOG-931 MAGs")


p_BOG2

####### Cool, now it is looking like it is supposed to #####
############## Next up, I need to scale it all to get the nxrA aligned ##############

min_nxrA_value <- min(BOG_genes_core_operon$start_position2[BOG_genes_core_operon$gene_short == "nxrA"])
max_nxrA_value <- max(BOG_genes_core_operon$start_position2[BOG_genes_core_operon$gene_short == "nxrA"])

### I want to add the difference between the nxr_gene and the other genes to each id:
BOG_genes_core_operon_2<-BOG_genes_core_operon%>%
  group_by(fasta)%>%
  filter(gene_short == "nxrA")%>%
  mutate(nxr_diff=max_nxrA_value-start_position2)

BOG_genes_core_operon <- left_join(BOG_genes_core_operon, BOG_genes_core_operon_2 %>% select(fasta, nxr_diff), by = "fasta")%>%
  mutate(start_scaled=start_position2+nxr_diff)%>%
  mutate(end_scaled=start_scaled+(end_position2-start_position2))



#### Now I need to generate colours for the unique genes ####
unique(BOG_genes_core_operon$gene_short)


BOG_genes_core_operon_col <- BOG_genes_core_operon %>%
  mutate(domain=if_else(gene_short=="nxrA", paste0("RE|", start_scaled,"|", end_scaled, "|#71A966|", gene_short), gene_short))%>%
  mutate(domain=if_else(gene_short=="nxrB/narH", paste0("RE|", start_scaled,"|", end_scaled, "|#8DBC80|", gene_short), domain))%>%
  mutate(domain=if_else(gene_short=="nxrX", paste0("RE|", start_scaled,"|", end_scaled, "|#B8D9A9|", gene_short), domain))%>%
  mutate(domain=if_else(gene_short=="narJ", paste0("RE|", start_scaled,"|", end_scaled, "|#859900|", gene_short), domain))%>%
  mutate(domain=if_else(gene_short=="narI", paste0("RE|", start_scaled,"|", end_scaled, "|#5D9D52|", gene_short), domain))%>%
  mutate(domain=if_else(gene_short=="NarK/NasA", paste0("RE|", start_scaled,"|", end_scaled, "|#74AC69|", gene_short), domain))%>%
  mutate(domain=if_else(gene_short=="cbbS", paste0("RE|", start_scaled,"|", end_scaled, "|#65A6C7|", gene_short), domain))%>%
  mutate(domain=if_else(gene_short=="cbbL", paste0("RE|", start_scaled,"|", end_scaled, "|#115379|", gene_short), domain))%>%
  mutate(domain=if_else(gene_short=="cbbF", paste0("RE|", start_scaled,"|", end_scaled, "|#8BC3DD|", gene_short), domain))%>%
  mutate(domain=if_else(gene_short=="cbbX", paste0("RE|", start_scaled,"|", end_scaled, "|#ADD0E1|", gene_short), domain))%>%
  mutate(domain=if_else(gene_short=="rpe", paste0("RE|", start_scaled,"|", end_scaled, "|#7EB0CB|", gene_short), domain))%>%
  mutate(domain=if_else(gene_short=="cyt. c class I", paste0("RE|", start_scaled,"|", end_scaled, "|#af90d6|", gene_short), domain))%>%
  mutate(domain=if_else(gene_short=="cyt. c subunit II", paste0("RE|", start_scaled,"|", end_scaled, "|#FDD474|", gene_short), domain))%>%
  mutate(domain=if_else(gene_short=="cyt. c subunit III", paste0("RE|", start_scaled,"|", end_scaled, "|#F28278|", gene_short), domain))%>%
  mutate(domain=if_else(gene_short=="cyt. c subunit IV", paste0("RE|", start_scaled,"|", end_scaled, "|#F59F93|", gene_short), domain))%>%
  mutate(domain=if_else(gene_short=="cyt. c CtaG", paste0("RE|", start_scaled,"|", end_scaled, "|#FFD17A|", gene_short), domain))%>%
  mutate(domain=if_else(gene_short=="transposase", paste0("RE|", start_scaled,"|", end_scaled, "|#B395C1|", gene_short), domain))%>%
  mutate(domain=if_else(gene_short=="integrase", paste0("RE|", start_scaled,"|", end_scaled, "|#4C2367|", gene_short), domain))%>%
  mutate(domain=if_else(gene_short=="formate/nitrite transporter", paste0("RE|", start_scaled,"|", end_scaled, "|#12d1be|", gene_short), domain))


genes_palette <- structure(list(
  `nxrA` = "#70A866FF",
  `nxrB/narH` = "#97C08CFF",
  `nxrX` = "#11451EFF",
  `narJ` = "#859900",
  `narI` = "#B8D9A9",
  `NarK/NasA` = "#40BFA2",
  `cbbS` = "#65A6C7",
  `cbbL` = "#115379",
  `cbbF` = "#8BC3DD",
  `cbbX` = "#ADD0E1",
  `rpe` = "#7EB0CB",
  `cyt. c class I` = "#FFD17A",
  `cyt. c subunit I` =  "#FC5200FF", 
  `cyt. c subunit II` = "#FEA03CFF",
  `cyt. c subunit III` = "#F28278",
  `cyt. c subunit IV` = "#F59F93",
  `cyt. c CtaG` = "#FC5200FF",
  `cyt. c` ="coral",
  `transposase` = "#B395C1",
  `integrase` = "#4C2367",
  `formate/nitrite transporter` = "#12d1be"
), class = "character", row.names = c(
  "nxrA", "nxrB/narH", "nxrX", "narJ", "narI", "NarK/NasA", "cbbS",
  "cbbL", "cbbF", "cbbX", "rpe", "cyt. c class I","cyt. c subunit I", "cyt. c subunit II",
  "cyt. c subunit III", "cyt. c subunit IV", "cyt. c CtaG","cyt. c","transposase",
  "integrase", "formate/nitrite transporter"
))



p_BOG2<- BOG_genes_core_operon %>%
  filter(grepl("BOG-931", GTDB_tax))%>%
  mutate(direction=if_else(strandedness=="1", T, F))%>%
  filter(if_else(fasta == "MFD09603.bin.c.1", 
                 start_position > 1680000 & end_position < 1725000, 
                 start_position > 0))%>%
  filter(if_else(fasta == "MFD00991.bin.1.40", 
                 start_position > 830062 & end_position < 856062, 
                 start_position > 0))%>% 
  filter(if_else(fasta == "MFD01872.bin.2.154", 
                 start_position < 252000 & end_position > 210000, 
                 start_position > 0))%>% 
  filter(if_else(fasta == "MFD03399.bin.2.43", 
                 start_position > 300002 , 
                 start_position > 0))%>%
  filter(if_else(fasta == "MFD05705.bin.1.139", 
                 start_position < 65000 & start_position > 30000, 
                 start_position > 0))%>%
  mutate(fasta = factor(fasta, levels = c("MFD09603.bin.c.1", "MFD05705.bin.1.139", "MFD00991.bin.1.40", "MFD01872.bin.2.154"), ordered = TRUE)) %>%
  mutate(start = (start_scaled + end_scaled)/2)  %>%
  ggplot(., aes(xmin = start_scaled, xmax = end_scaled, y = fasta, fill = gene_short)) +
  geom_gene_arrow(aes(forward=direction)) +
  facet_wrap(~ fasta, scales = "free_y", ncol=1) +
  theme_genes()+
  theme(axis.title = element_blank(),
        # axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = genes_palette) +
  geom_text(aes(x=start_scaled, label = gene_short), 
            nudge_y = 0.3, angle = 8,
            hjust= 0.0,
            size= 3, family="Arial")+
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size=10, color="grey40"),
        axis.text.y = element_text(size=10, color="grey10"))+
  labs(title="HQ BOG-931 MAGs")

p_BOG2

#### Now to saving this ####

ggsave("output/BOG_HQ.png", p_BOG2, width=20, height=4)
ggsave("output/BOG_HQ.svg", p_BOG2, width=20, height=4)


library(paletteer) 
paletteer_dynamic("cartography::orange.pal", 10)

##### Now I want to investigate the other regions of the genomes for anything interesting #######



MFD00991.bin.1.40<-d2 %>%
  filter(fasta=="MFD00991.bin.1.40")%>%
  mutate(start_2=start_position-min(start_position))%>%
  mutate(end_2=start_2+(end_position-start_position))


###### Nothing more nitrate-related here


MFD01872.bin.2.154<-d2 %>%
  filter(fasta=="MFD01872.bin.2.154")%>%
  mutate(start_2=start_position-min(start_position))%>%
  mutate(end_2=start_2+(end_position-start_position))
#### Want to expand to 295750


p_MFD01872.bin.2.154<- MFD01872.bin.2.154 %>%
  filter(if_else(fasta == "MFD01872.bin.2.154", 
                 start_position < 300000 & end_position > 210000, 
                 start_position > 0))%>% 
  mutate(start = (start_2 + end_2)/2)  %>%
  ggplot(., aes(xmin = start_2, xmax = end_2, y = fasta, fill = gene_short)) +
  geom_gene_arrow() +
  facet_wrap(~ fasta, scales = "free_y", ncol=1) +
  theme_genes()+
  theme(axis.title = element_blank(),
        # axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
  #scale_fill_manual(values = genes_palette) +
  geom_text(aes(x=start_2, label = gene_short), 
            nudge_y = 0.3, angle = 8,
            hjust= 0.0,
            size= 3, family="Arial")+
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size=10, color="grey40"),
        axis.text.y = element_text(size=10, color="grey10"))+
  labs(title="HQ BOG-931 MAGs")

p_MFD01872.bin.2.154




MFD05705.bin.1.139<-d2 %>%
  filter(fasta=="MFD05705.bin.1.139")%>%
  mutate(start_2=start_position-min(start_position))%>%
  mutate(end_2=start_2+(end_position-start_position))

### Only more transposases and other stuff


MFD09603.bin.c.1<-d2 %>%
  filter(fasta=="MFD09603.bin.c.1")%>%
  mutate(start_2=start_position-min(start_position))%>%
  mutate(end_2=start_2+(end_position-start_position))

#### Okay, here from 229652 to 351414


p_MFD09603.bin.c.1<- MFD09603.bin.c.1 %>%
  filter(start_position > 239652)%>% 
  filter(end_position < 246662 )%>%
  mutate(start = (start_2 + end_2)/2)  %>%
  ggplot(., aes(xmin = start_2, xmax = end_2, y = fasta, fill = gene_short)) +
  geom_gene_arrow() +
  facet_wrap(~ fasta, scales = "free_y", ncol=1) +
  theme_genes()+
  theme(axis.title = element_blank(),
        # axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = genes_palette) +
  geom_text(aes(x=start_2, label = gene_short), 
            nudge_y = 0.2, angle = 8,
            hjust= 0.0,
            size= 3, family="Arial")+
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size=10, color="grey40"),
        axis.text.y = element_text(size=10, color="grey10"))+
  labs(title="HQ BOG-931 MAGs")

p_MFD09603.bin.c.1

ggsave("output/p_MFD09603.bin.c.1.png", p_MFD09603.bin.c.1, width=4.4, height=2)
ggsave("output/p_MFD09603.bin.c.1.svg", p_MFD09603.bin.c.1, width=4.4, height=2)

p_BOG2_full_length<- BOG_genes_core_operon %>%
  filter(grepl("BOG-931", GTDB_tax))%>%
  filter(if_else(fasta == "MFD00991.bin.1.40", 
                 start_position > 830062, 
                 start_position > 0))%>% 
  mutate(direction=if_else(strandedness=="1", T, F))%>%
  mutate(fasta = factor(fasta, levels = c("MFD09603.bin.c.1", "MFD05705.bin.1.139", "MFD00991.bin.1.40", "MFD01872.bin.2.154"), ordered = TRUE)) %>%
  mutate(start = (start_scaled + end_scaled)/2)  %>%
  ggplot(., aes(xmin = start_scaled, xmax = end_scaled, y = fasta, fill = gene_short)) +
  geom_gene_arrow(aes(forward=direction)) +
  facet_wrap(~ fasta, scales = "free_y", ncol=1) +
  theme_genes()+
  theme(axis.title = element_blank(),
        # axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = genes_palette) +
  geom_text(aes(x=start_scaled, label = gene_short), 
            nudge_y = 0.3, angle = 8,
            hjust= 0.0,
            size= 3, family="Arial")+
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size=10, color="grey40"),
        axis.text.y = element_text(size=10, color="grey10"))+
  labs(title="HQ BOG-931 MAGs")

ggsave("output/p_BOG2_full_length.png", p_BOG2_full_length, width=30, height=4)
ggsave("output/p_BOG2_full_length.svg", p_BOG2_full_length, width=30, height=4)





#################### Searching in the SR-HQ MAGs #####################

SR_KEGG <- vroom("data/annotations_nitrifiers_in_paper.tsv", delim="\t", col_select=c("fasta", "ko_id", "kegg_genes_id", "kegg_hit", "pfam_hits", "scaffold", "gene_position", "start_position", "end_position", "strandedness")) %>%
 separate(kegg_hit, into = c("kegg_genes_id_2", "gene"), sep = '\\)  ', remove = T)%>%
  select(!kegg_genes_id_2)

unique(SR_KEGG$fasta)

KOs <- c("K01601","K01602", "K03775")

SR_HQ <- SR_KEGG%>%
  filter(grepl("LIB-MJ138-E2_03.20|LIB-MJ138-E2_03.2|LIB-MJ342-B3_02.17|LIB-MJ050-E10_02.3|LIB-MJ176-C4_03.14", fasta))%>%
  filter(grepl("nitrate|molybdopterin-dependent oxidoreductase|putative oxidoreductase molybdopterin-binding subunit|cytochrome c|cytochrome C|FKBP|transposase|ribulose|cbb|integrase|nitrite transporter", gene) |
           grepl("K01601|K01602|K03775", ko_id)) # %>% filter(gene =="nitrate oxidoreductase alpha subunit")


unique(SR_HQ$fasta)



cbbL_SR<-SR_HQ%>%
  filter(ko_id %in% c("K01601","K01602"))




##### Investigating operon structure for the SR-MAGs ####


SR_KEGG <- vroom("data/annotations_nitrifiers_in_paper.tsv", delim="\t", col_select=c("fasta", "ko_id", "kegg_genes_id", "kegg_hit", "pfam_hits", "scaffold", "gene_position", "start_position", "end_position", "strandedness")) %>%
separate(kegg_hit, into = c("kegg_genes_id_2", "gene"), sep = '\\)  ', remove = T)%>%
  select(!kegg_genes_id_2)

unique(SR_KEGG$fasta)

KOs <- c("K01601","K01602", "K03775")

SR_HQ <- SR_KEGG%>%
  filter(grepl("LIB-MJ138-E2_03.20|LIB-MJ138-E2_03.2|LIB-MJ342-B3_02.17|LIB-MJ050-E10_02.3|LIB-MJ176-C4_03.14", fasta))%>%
  filter(grepl("respiratory nitrate reductase|FKBP|cbb|cytochrome c, class I", gene) |
           grepl("K01601|K01602|K03775|K00370", ko_id)) # %>% filter(gene =="nitrate oxidoreductase alpha subunit")


tax<-vroom("data/nxr_mags_combined.tsv", delim="\t", col_names = c("fasta", "contig","gene","quality","GTDB_tax"))%>%
  separate(contig, into = c("contig_text","scaffold","start_position", "frame", "gene_no"), sep = '_', remove=FALSE)%>%
  mutate(scaffold=paste0("contig_", scaffold))%>%
  mutate(gene=gsub("Root; Cytoplasmic_nxr_narG; ", "", gene))%>%
  mutate(strandedness="1")

d_filter<-SR_HQ%>%
  select(fasta, gene, scaffold, start_position, end_position, strandedness) %>%
  merge(tax[, c("fasta", "GTDB_tax")], by = "fasta")%>%
  distinct()%>%
  filter()

d2<-d_filter%>%
  mutate(gene=if_else(grepl("ABC", gene), paste0("put. ABC-nitrate transp."), paste0(gene)))%>%
  mutate(gene=gsub("respiratory nitrate reductase", "nar", gene))%>%
  mutate(gene=gsub("cbbL; Ribulose bisphosphate carboxylase large chain \\(RuBisCO large subunit\\)", "cbbL", gene)) %>%
  mutate(gene=gsub("ribulose bisphosphate carboxylase small subunit", "cbbS", gene))%>%
  mutate(gene=gsub("ribulose-bisphosphate carboxylase large subunit", "cbbL", gene))%>%
  mutate(gene=gsub("ribulose-bisphosphate carboxylase large subunit", "cbbL", gene))%>%
  mutate(gene=gsub("form I ribulose bisphosphate carboxylase large subunit", "cbbL", gene))%>%
  mutate(gene=gsub("ribulose 1,5-bisphosphate carboxylase", "cbb", gene))%>%
  mutate(gene=gsub("Ribulose-bisphosphate carboxylase", "cbbL", gene))%>%
  mutate(gene=gsub("Ribulose bisphosphate carboxylase large chain", "cbbL", gene))%>%
  mutate(gene=gsub("cbbS; Ribulose bisphosphate carboxylase small chain \\(RuBisCO small subunit\\)", "cbbS", gene))%>%
  mutate(gene=gsub("rbcL; cbbL", "cbbL", gene))%>%
  mutate(gene_short=if_else(grepl("Nitrobacter", gene), paste0("nxrA"), gene))%>%
  mutate(gene_short=if_else(grepl("nar beta", gene_short), paste0("nxrB/narH"), gene_short))%>%
  mutate(gene_short=if_else(grepl("NarJ", gene_short), paste0("narJ"), gene_short))%>%
  mutate(gene_short=if_else(grepl("nar gamma", gene_short), paste0("narI"), gene_short))%>%
  mutate(gene_short=if_else(grepl("alpha subunit", gene_short), paste0("nar a"), gene_short))%>%
  mutate(gene_short=if_else(grepl("cytochrome c, class I", gene_short), paste0("cyt. c class I"), gene_short))%>%
  mutate(gene_short=if_else(grepl("NarK/NasA", gene_short), paste0("NarK/NasA"), gene_short))%>%
  mutate(gene_short=if_else(grepl("FKBP", gene_short), paste0("nxrX"), gene_short))%>%
  mutate(gene_short=if_else(grepl("transposase", gene_short), paste0("transposase"), gene_short))%>%
  mutate(gene_short=if_else(grepl("integrase", gene_short), paste0("integrase"), gene_short))%>%
  mutate(gene_short=if_else(grepl("cbbF", gene_short), paste0("cbbF"), gene_short))%>%
  mutate(gene_short=if_else(grepl("cbbX", gene_short), paste0("cbbX"), gene_short))%>%
  mutate(gene_short=if_else(grepl("cbbR", gene_short), paste0("cbbR"), gene_short))%>%
  mutate(gene_short=if_else(grepl("aa3-type", gene_short), paste0("cyt. c subunit IV"), gene_short))%>%
  mutate(gene_short=if_else(grepl("CtaG", gene_short), paste0("cyt. c CtaG"), gene_short))%>%
  mutate(gene_short=if_else(grepl("rpe", gene_short), paste0("rpe"), gene_short))%>%
  mutate(gene_short=if_else(grepl("cytochrome C oxidase subunit II|cytochrome c oxidase subunit II", gene_short), paste0("cyt. c subunit II"), gene_short))%>%
  mutate(gene_short=if_else(grepl("cytochrome c oxidase, subunit III", gene_short), paste0("cyt. c subunit III"), gene_short))%>%
  mutate(gene_short=if_else(grepl("cytochrome c oxidase subunit 3", gene_short), paste0("cyt. c subunit III"), gene_short))%>%
  mutate(gene_short=if_else(grepl("cytochrome c family protein", gene_short), paste0("cyt. c"), gene_short))%>%
  mutate(gene_short=if_else(grepl("ctaD", gene_short), paste0("cyt. c subunit I"), gene_short))%>%
  mutate(gene_short=if_else(grepl("Formate/nitrite transporter", gene_short), paste0("formate/nitrite transporter"), gene_short))







genes_palette <- structure(list(
  `nxrA` = "#70A866FF",
  `nxrB/narH` = "#97C08CFF",
  `nxrX` = "#11451EFF",
  `narJ` = "#859900",
  `narI` = "#B8D9A9",
  `NarK/NasA` = "#40BFA2",
  `cbbS` = "#65A6C7",
  `cbbL` = "#115379",
  `cbbF` = "#8BC3DD",
  `cbbX` = "#ADD0E1",
  `cbbR`="steelblue",
  `rpe` = "#7EB0CB",
  `cyt. c class I` = "#FFD17A",
  `cyt. c subunit I` =  "#FC5200FF", 
  `cyt. c subunit II` = "#FEA03CFF",
  `cyt. c subunit III` = "#F28278",
  `cyt. c subunit IV` = "#F59F93",
  `cyt. c CtaG` = "#FC5200FF",
  `cyt. c` ="coral",
  `transposase` = "#B395C1",
  `integrase` = "#4C2367",
  `formate/nitrite transporter` = "#12d1be"
), class = "character", row.names = c(
  "nxrA", "nxrB/narH", "nxrX", "narJ", "narI", "NarK/NasA", "cbbS",
  "cbbL", "cbbF", "cbbX","cbbR", "rpe", "cyt. c class I","cyt. c subunit I", "cyt. c subunit II",
  "cyt. c subunit III", "cyt. c subunit IV", "cyt. c CtaG","cyt. c","transposase",
  "integrase", "formate/nitrite transporter"
))


unique(d2$scaffold)

p_BOG<- d2 %>%
  filter(grepl("BOG-931", GTDB_tax))%>%
  filter(!grepl("cbb3-type",gene_short))%>%
  filter(!grepl("197778|23284|20177|245492|196551|124944", scaffold))%>%
  mutate(contig=paste0(fasta," | ",gsub("k127","contig", scaffold)))%>%
  mutate(direction=if_else(strandedness=="1", T, F))%>%
  mutate(gene_short=gsub("nar a", "nxrA", gene_short))%>%
  mutate(start = (start_position + end_position)/2)  %>%
  ggplot(., aes(xmin = start_position, xmax = end_position, y = contig, fill = gene_short)) +
  geom_gene_arrow(aes(forward=direction)) +
  facet_wrap(~ fasta, scales = "free_y", ncol=1) +
  theme_genes()+
  theme(axis.title = element_blank(),
        # axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right") +
  #scale_fill_manual(values = custom_colors) +
  scale_fill_manual(values = genes_palette) +
  geom_text(aes(x=start_position, label = gene_short), 
            nudge_y = 0.3, angle = 8,
            hjust= 0,
            size= 3, family="Arial")+
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size=10, color="grey40"),
        axis.text.y = element_text(size=10, color="grey10"),
        legend.position = "none")+
  labs(title="SR BOG-931 MAGs")

p_BOG

ggsave("output/SR_operons_BOG-931.png", p_BOG, width=12, height=10)
ggsave("output/SR_operons_BOG-931.svg", p_BOG, width=12, height=10)







####### Formate/Nitrite transporters ######

KEGG <- vroom("data/annotations_nitrifiers_in_paper.tsv", delim="\t", col_select=c("fasta", "ko_id", "kegg_genes_id", "kegg_hit", "pfam_hits", "scaffold", "gene_position", "start_position", "end_position", "strandedness")) %>%
  separate(kegg_hit, into = c("kegg_genes_id_2", "gene"), sep = '\\)  ', remove = T)%>%
  select(!kegg_genes_id_2)


d <- KEGG%>%
  filter(grepl("Formate/nitrite transporter [PF01226.21]", pfam_hits) |
         grepl("/nitrite transporter", gene)) 


tax<-vroom("data/Nitrobacter_like_nxr_mags_combined.tsv", delim="\t", col_names = c("fasta", "contig","gene","quality","GTDB_tax"))%>%
  filter(quality=="HQ")%>%
  separate(contig, into = c("contig_text","scaffold","start_position", "frame", "gene_no"), sep = '_', remove=FALSE)%>%
  select(-c(contig_text, contig, frame, quality, gene_no))%>%
  filter(grepl("BOG", GTDB_tax))

d_filter<-d%>%
#  filter(scaffold %in% unique(tax$scaffold))%>% #### Modify here to filter on contigs 
  select(fasta, gene, scaffold, start_position, end_position, strandedness, gene_position) %>%
  merge(tax[, c("fasta", "GTDB_tax")], by = "fasta")%>%
  distinct()%>%
  mutate(length=end_position-start_position)


export<-d_filter%>%
  mutate(sequence=paste0(fasta,"_",scaffold,"_",gene_position))%>%
  select(fasta, sequence)

write_excel_csv(export, "output/searchfile.txt", col_names = F, quote = "none")



#### For the SR:
SR_KEGG <- vroom("data/annotations_nitrifiers_in_paper.tsv", delim="\t", col_select=c("fasta", "ko_id", "kegg_genes_id", "kegg_hit", "pfam_hits", "scaffold", "gene_position", "start_position", "end_position", "strandedness")) %>%
 separate(kegg_hit, into = c("kegg_genes_id_2", "gene"), sep = '\\)  ', remove = T)%>%
  select(!kegg_genes_id_2)

unique(SR_KEGG$fasta)

KOs <- c("K01601","K01602", "K03775")

SR_HQ <- SR_KEGG%>%
  filter(grepl("LIB-MJ138-E2_03.20|LIB-MJ138-E2_03.2|LIB-MJ342-B3_02.17|LIB-MJ050-E10_02.3|LIB-MJ176-C4_03.14", fasta))%>%
  filter(grepl("Formate/nitrite transporter [PF01226.21]", pfam_hits) |
           grepl("/nitrite transporter", gene)) %>%
  mutate(length=end_position-start_position)  %>%
  mutate(sequence=paste0(fasta,"_",scaffold,"_",gene_position))%>%
  select(fasta, sequence)

write_excel_csv(SR_HQ, "output/SR_MAGs_searchfile.txt", col_names = F, quote = "none")


#### Just obtaining the genomes that have the transporters ####
unique(SR_HQ$fasta) ### That is all nxr-containing genomes
unique(export$fasta) ### That is all nxr-containing genomes
