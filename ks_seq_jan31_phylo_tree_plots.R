library(ggtree)
library(phytools)
library(tidytree)
library(seqinr)
library(data.table)
library(ggplot2)
library(treeio)
library(RColorBrewer)

seq_metadata = data.table(read.csv("/Users/tbrockfi/Desktop/projects/covid_data_analysis/lhdseq/ks_seq/phylo_tree/genbank_data.csv"))
BA1_samples = seq_metadata[grepl("BA.1",seq_metadata$Pangolin),]
BA2_samples = seq_metadata[grepl("BA.2",seq_metadata$Pangolin),]

#for each MA and USA get 10 BA.1 or desc, and 50 BA.2 or desc
#ma_seqs = c(BA1_samples[USA=="MA",][sample(.N,10)]$Accession, BA2_samples[USA=="MA",][sample(.N,50)]$Accession)
#usa_seqs = c(BA1_samples[USA!="MA",][sample(.N,10)]$Accession, BA2_samples[USA!="MA",][sample(.N,50)]$Accession)

tree=read.tree("/Users/tbrockfi/Desktop/projects/covid_data_analysis/lhdseq/ks_seq/phylo_tree/rooted_dec5_phylo_tree.nwk")
study_tips = grepl("-",tree$tip.label)
ma_tips = tree$tip.label %in% paste("'",seq_metadata[USA=="MA",]$Accession,".1'",sep="")

meta_data = data.frame(ID = tree$tip.label, group = rep("USA",length(tree$tip.label))) 
meta_data[study_tips,]$group = paste("HH",substr(tree$tip.label[study_tips],1,3),sep="")
meta_data[ma_tips,]$group = "MA"
merged_tree_data = full_join(as_tibble(tree), meta_data, by = c("label" = "ID"))

annotated_tree = as.treedata(merged_tree_data)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

col_pal = c("lightgrey","grey68",col_vector[1:11],col_vector[12:22]) #12 rows per col

group_list = c("USA", "MA",sort(unique(paste("HH",substr(tree$tip.label[study_tips],1,3),sep="")))[1:11], 
               sort(unique(paste("HH",substr(tree$tip.label[study_tips],1,3),sep="")))[12:22])


ggtree(annotated_tree,color="black",linewidth=0.3) + coord_flip() + scale_x_reverse() +
  geom_tippoint(aes(color=group),size=1.4) + 
  guides(color = guide_legend(override.aes = list(size = 5),ncol=6),reverse=FALSE) + 
  scale_color_manual(values=col_pal, limits = group_list) +
  theme(legend.title=element_blank(),legend.position = c(0.8,0.72)) +
  geom_treescale(x=-0.0004,y=60,offset=1,width=(3/29903),fontsize = 0,linesize=1.1) + #
  guides(fill = guide_legend(reverse = TRUE)) +
  annotate(geom="text", x=0.00035, y=68, label="5 substitutions") + 
  geom_segment(x=0,xend=2,y=17,yend=18, lwd=0.2) + 
  annotate(geom="text", x=0.0001, y=10, label="BA.1 and \ndescendants") + 
  annotate(geom="text", x=0.0001, y=37, label="BA.2 and \ndescendants") 
