#!/usr/bin/Rscript

# Get working directory from command line arguments
args <- commandArgs(trailingOnly = TRUE)
working_directory <- args[1]

# Set working directory
setwd(working_directory)

.libPaths(c("/home/bio.aau.dk/vj52ou/software/R_packages.v.4.3.2", .libPaths()))

#library(tidyr)
library(stringr)
library(purrr)
library(dplyr)
library(vroom)
library(ggplot2)

#Importing the position file 
hits<-vroom("position_files/position_of_hits.tsv", delim = "\t")
hits <- hits %>%
  mutate(Taxonomy = str_extract(Tax, "[^;]*$"))

library(tidyr)
# Step 1: Create a data frame with every position of the gene
gene_positions <- hits %>%
  summarise(start = min(hmm_start),
            end = max(hmm_end)) %>%
  expand(Position = seq(start, end))


# Step 2: Join the gene_positions data frame with the hits data frame and group by "Taxonomy"
reads_by_position <- gene_positions %>%
  left_join(hits %>% 
              mutate(Position = map2(hmm_start, hmm_end, seq)) %>%
              unnest(Position) %>%
              group_by(Position, Taxonomy) %>%
              count()) %>%
  replace_na(list(n = 0)) %>%
  rename(Count=n)



# Step 4: Visualize the distribution of reads along the gene by Taxonomy
cov_tax<-ggplot(reads_by_position, aes(x = Position, y = Count)) + 
  geom_line(linewidth=1)+
  facet_wrap(~Taxonomy, scales = "free_y")+  
  ylab("Total coverage within each taxonomic group")+
  theme(panel.background = element_rect(fill = "transparent", colour = NA))+
  labs(title="E value cutoff e-10")+
  theme(plot.title = element_text(hjust=0.5, face="bold"))

ggsave("cov_e10.png",
       cov_tax,
       height = 10,
       width = 17)



#### I also want to plot parent level and sub 1;

hits<-vroom("position_files/position_of_hits.tsv", delim = "\t")
hits <- hits %>%
separate(Tax, into = c("Root", "Parent"), sep = '; ', remove=FALSE) %>%
mutate(Parent = if_else(is.na(Parent), paste0("Root"), Parent))

# Step 1: Create a data frame with every position of the gene
gene_positions <- hits %>%
  summarise(start = min(hmm_start),
            end = max(hmm_end)) %>%
  expand(Position = seq(start, end))


# Step 2: Join the gene_positions data frame with the hits data frame and group by "Taxonomy"
reads_by_position <- gene_positions %>%
  left_join(hits %>% 
              mutate(Position = map2(hmm_start, hmm_end, seq)) %>%
              unnest(Position) %>%
              group_by(Position, Parent) %>%
              count()) %>%
  replace_na(list(n = 0)) %>%
  rename(Count=n)



# Step 4: Visualize the distribution of reads along the gene by Taxonomy
cov_tax<-ggplot(reads_by_position, aes(x = Position, y = Count)) + 
  geom_line(linewidth=1)+
  facet_wrap(~Parent, scales = "free_y")+  
  ylab("Total coverage within each taxonomic group")+
  theme(panel.background = element_rect(fill = "transparent", colour = NA))+
  labs(title="Coverage parent level")+
  theme(plot.title = element_text(hjust=0.5, face="bold"))

ggsave("cov_parent.png",
       cov_tax,
       height = 10,
       width = 15)


cov_tax_com<-ggplot(reads_by_position, aes(x = Position, y = Count, color =Parent)) + 
  geom_line(linewidth=1.3)+
  #facet_wrap(~Taxonomy, scales = "free_y")+  
  #scale_color_manual(values=c("#4c5e30", "#7e2c2c", "#275c8d","#c26129"))+
  ylab("Total coverage within each taxonomic group")+
  theme(panel.background = element_rect(fill = "transparent", colour = NA))+
  labs(title="Coverage parent level", subtitle = "Combined")+
  theme(plot.title = element_text(hjust=0, face="bold"))


ggsave("cov_parent_com.png",
       cov_tax_com,
       height = 4,
       width = 12)



###### Sub 1 ######
#Importing the position file 
hits<-vroom("position_files/position_of_hits.tsv", delim = "\t")
hits <- hits %>%
separate(Tax, into = c("Root", "Parent", "Sub1"), sep = '; ', remove=FALSE) %>%
mutate(Parent = if_else(is.na(Parent), paste0("Root"), Parent)) %>%
mutate(Sub1 = if_else(is.na(Sub1), paste0(Parent), Sub1))

# Step 1: Create a data frame with every position of the gene
gene_positions <- hits %>%
  summarise(start = min(hmm_start),
            end = max(hmm_end)) %>%
  expand(Position = seq(start, end))


# Step 2: Join the gene_positions data frame with the hits data frame and group by "Taxonomy"
reads_by_position <- gene_positions %>%
  left_join(hits %>% 
              mutate(Position = map2(hmm_start, hmm_end, seq)) %>%
              unnest(Position) %>%
              group_by(Position, Parent, Sub1) %>%
              count()) %>%
  replace_na(list(n = 0)) %>%
  rename(Count=n)



# Step 4: Visualize the distribution of reads along the gene by Taxonomy
cov_tax_com<-ggplot(reads_by_position, aes(x = Position, y = Count, color = Sub1)) + 
  geom_line(linewidth=1.3)+
  #facet_wrap(~Taxonomy, scales = "free_y")+  
  #scale_color_manual(values=c("#4c5e30", "#7e2c2c", "#275c8d","#c26129"))+
  ylab("Total coverage within each taxonomic group")+
  theme(panel.background = element_rect(fill = "transparent", colour = NA))+
  labs(title="Coverage profile", subtitle = "Taxonomic depth 1")+
  theme(plot.title = element_text(hjust=0, face="bold"))

ggsave("cov_sub1.png",
       cov_tax_com,
       height = 4,
       width = 12)


###### And sub 2


#Importing the position file 
hits<-vroom("position_files/position_of_hits.tsv", delim = "\t")
hits <- hits %>%
separate(Tax, into = c("Root", "Parent", "Sub1", "Sub2"), sep = '; ', remove=FALSE) %>%
mutate(Parent = if_else(is.na(Parent), paste0("Root"), Parent)) %>%
mutate(Sub1 = if_else(is.na(Sub1), paste0(Parent), Sub1))%>%
mutate(Sub2 = if_else(is.na(Sub2), paste0(Sub1), Sub2))

# Step 1: Create a data frame with every position of the gene
gene_positions <- hits %>%
  summarise(start = min(hmm_start),
            end = max(hmm_end)) %>%
  expand(Position = seq(start, end))


# Step 2: Join the gene_positions data frame with the hits data frame and group by "Taxonomy"
reads_by_position <- gene_positions %>%
  left_join(hits %>% 
              mutate(Position = map2(hmm_start, hmm_end, seq)) %>%
              unnest(Position) %>%
              group_by(Position, Sub2) %>%
              count()) %>%
  replace_na(list(n = 0)) %>%
  rename(Count=n)



# Step 4: Visualize the distribution of reads along the gene by Taxonomy
cov_tax_com<-ggplot(reads_by_position, aes(x = Position, y = Count, color = Sub2)) + 
  geom_line(linewidth=1.3)+
  #facet_wrap(~Taxonomy, scales = "free_y")+  
  #scale_color_manual(values=c("#4c5e30", "#7e2c2c", "#275c8d","#c26129"))+
  ylab("Total coverage within each taxonomic group")+
  theme(panel.background = element_rect(fill = "transparent", colour = NA))+
  labs(title="Coverage profile", subtitle = "Taxonomic depth 2")+
  theme(plot.title = element_text(hjust=0, face="bold"))

ggsave("cov_sub2.png",
       cov_tax_com,
       height = 4,
       width = 12)



###### And sub 3


#Importing the position file 
hits<-vroom("position_files/position_of_hits.tsv", delim = "\t")
hits <- hits %>%
separate(Tax, into = c("Root", "Parent", "Sub1", "Sub2", "Sub3"), sep = '; ', remove=FALSE) %>%
mutate(Parent = if_else(is.na(Parent), paste0("Root"), Parent)) %>%
mutate(Sub1 = if_else(is.na(Sub1), paste0(Parent), Sub1))%>%
mutate(Sub2 = if_else(is.na(Sub2), paste0(Sub1), Sub2))%>%
mutate(Sub3 = if_else(is.na(Sub3), paste0(Sub2), Sub3))

# Step 1: Create a data frame with every position of the gene
gene_positions <- hits %>%
  summarise(start = min(hmm_start),
            end = max(hmm_end)) %>%
  expand(Position = seq(start, end))


# Step 2: Join the gene_positions data frame with the hits data frame and group by "Taxonomy"
reads_by_position <- gene_positions %>%
  left_join(hits %>% 
              mutate(Position = map2(hmm_start, hmm_end, seq)) %>%
              unnest(Position) %>%
              group_by(Position, Sub3) %>%
              count()) %>%
  replace_na(list(n = 0)) %>%
  rename(Count=n)



# Step 4: Visualize the distribution of reads along the gene by Taxonomy
cov_tax_com<-ggplot(reads_by_position, aes(x = Position, y = Count, color = Sub3)) + 
  geom_line(linewidth=1.3)+
  #facet_wrap(~Taxonomy, scales = "free_y")+  
  #scale_color_manual(values=c("#4c5e30", "#7e2c2c", "#275c8d","#c26129"))+
  ylab("Total coverage within each taxonomic group")+
  theme(panel.background = element_rect(fill = "transparent", colour = NA))+
  labs(title="Coverage profile", subtitle = "Taxonomic depth 3")+
  theme(plot.title = element_text(hjust=0, face="bold"))

ggsave("cov_sub3.png",
       cov_tax_com,
       height = 4,
       width = 12)



###### And sub 4


#Importing the position file 
hits<-vroom("position_files/position_of_hits.tsv", delim = "\t")
hits <- hits %>%
separate(Tax, into = c("Root", "Parent", "Sub1", "Sub2", "Sub3", "Sub4"), sep = '; ', remove=FALSE) %>%
mutate(Parent = if_else(is.na(Parent), paste0("Root"), Parent)) %>%
mutate(Sub1 = if_else(is.na(Sub1), paste0(Parent), Sub1))%>%
mutate(Sub2 = if_else(is.na(Sub2), paste0(Sub1), Sub2))%>%
mutate(Sub3 = if_else(is.na(Sub3), paste0(Sub2), Sub3))%>%
mutate(Sub4 = if_else(is.na(Sub4), paste0(Sub3), Sub4))

# Step 1: Create a data frame with every position of the gene
gene_positions <- hits %>%
  summarise(start = min(hmm_start),
            end = max(hmm_end)) %>%
  expand(Position = seq(start, end))


# Step 2: Join the gene_positions data frame with the hits data frame and group by "Taxonomy"
reads_by_position <- gene_positions %>%
  left_join(hits %>% 
              mutate(Position = map2(hmm_start, hmm_end, seq)) %>%
              unnest(Position) %>%
              group_by(Position, Sub4) %>%
              count()) %>%
  replace_na(list(n = 0)) %>%
  rename(Count=n)



# Step 4: Visualize the distribution of reads along the gene by Taxonomy
cov_tax_com<-ggplot(reads_by_position, aes(x = Position, y = Count, color = Sub4)) + 
  geom_line(linewidth=1.3)+
  #facet_wrap(~Taxonomy, scales = "free_y")+  
  #scale_color_manual(values=c("#4c5e30", "#7e2c2c", "#275c8d","#c26129"))+
  ylab("Total coverage within each taxonomic group")+
  theme(panel.background = element_rect(fill = "transparent", colour = NA))+
  labs(title="Coverage profile", subtitle = "Taxonomic depth 4")+
  theme(plot.title = element_text(hjust=0, face="bold"))

ggsave("cov_sub4.png",
       cov_tax_com,
       height = 4,
       width = 12)