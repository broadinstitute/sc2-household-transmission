#new master script for KS analysis using updated genomes
  #read in metadata, find longest genome per person, 
  #apply usual quality filters, find distance between hh pairs, and check iSNVs
library(data.table)
library(seqinr)
library(readxl)
library(phylotools) 
library(stringr) 
library(tidyr) 
library(purrr) #for max mutually read pair
library(dplyr) #for max mutually read pair
library(lubridate) #for correcting collection dates
library(tibble) 

#create depth table
setwd("/Users/tbrockfi/Desktop/projects/covid_data_analysis/lhdseq/ks_seq/all_seqs_corrected_sep5/all_seqs/all_depths")
dp_files = dir(pattern = ".*\\.txt")
dp_list = lapply(dp_files, function(x) read.delim(x, header=FALSE, comment.char="#"))
dp_sample_names = unlist(lapply(dp_files,function(x){strsplit(x,"\\.")[[1]][1]}))
merged_dps = data.table()
for(i in 1:length(dp_list)){
  merged_dps[,dp_sample_names[i]] = dp_list[[i]]$V3
}

setwd("/Users/tbrockfi/Desktop/projects/covid_data_analysis/lhdseq/ks_seq/all_seqs_corrected_sep5")
master_metadata = data.table(read_excel("/Users/tbrockfi/Desktop/projects/covid_data_analysis/lhdseq/ks_seq/ks_seq_dec5_clean/ks_seq_dec5_master_metadata.xlsx"))
colnames(master_metadata)[1] = "master_id"
colnames(master_metadata)[7] = "ks_id"
master_metadata$collection_date = master_metadata$date+days(master_metadata$day - 1)
master_metadata = master_metadata[order(household,pt,day)]
all_fastas = data.table(read.fasta("all_seqs/ks_nextclade_alignment.fasta"))
all_fastas = all_fastas[!duplicated(all_fastas)]
all_fastas$seq.name = unlist(lapply(all_fastas$seq.name, function(x){unlist(strsplit(x, "\\."))[1]}))
seq_and_metadata = master_metadata[all_fastas, on=.("broad_id"=seq.name)][order(household,pt,day)]
most_complete_sample_by_pt = master_metadata[,.SD[,.(household,pt,master_id,length_unamb)][which.max(length_unamb)],.(household,pt)]
most_complete_sample_by_pt_broad_id = seq_and_metadata[master_id %in% most_complete_sample_by_pt$master_id,]
filtered_most_complete_sample_by_pt_broad_id = most_complete_sample_by_pt_broad_id[length_unamb>=20000]


#making positional quality filters 
dp_cutoff = 20
problematic_sites = data.table(read.table("~/Desktop/projects/covid_data_analysis/lhdseq/problematic_sites_sarsCov2.vcf", quote="\""))
problematic_sites_list = problematic_sites[V7=="mask"]$V2
include_table = merged_dps
include_table = apply(include_table,2,function(x){x>dp_cutoff})
include_table[problematic_sites_list,] = FALSE

masked_most_complete_sample = filtered_most_complete_sample_by_pt_broad_id
#filter for depth and problematic sites
for(i in 1:dim(filtered_most_complete_sample_by_pt_broad_id)[1]){
  seq_name = filtered_most_complete_sample_by_pt_broad_id[i,broad_id]
  char_seq = s2c(filtered_most_complete_sample_by_pt_broad_id[i,seq.text])
  include_list = include_table[,seq_name]
  masked_char_seq = char_seq
  masked_char_seq[!include_list] = "-"
  masked_most_complete_sample[i,]$seq.text = c2s(masked_char_seq)
}


for(hh in unique(masked_most_complete_sample$household)){
  all_hh_seqs_strings = masked_most_complete_sample[household==hh,seq.text]
  all_hh_seqs_char = lapply(all_hh_seqs_strings,s2c)
  names(all_hh_seqs_char) = masked_most_complete_sample[household==hh,broad_id]
  seqs = as.DNAbin(all_hh_seqs_char)
  print(paste("Results for household:",hh,sep=" "))
  print(dist.dna(seqs,model="N"))
}


##### iSNV analysis #####
setwd("/Users/tbrockfi/Desktop/projects/covid_data_analysis/lhdseq/ks_seq/all_seqs_corrected_sep5/all_seqs/all_vcfs")

all_drop_na = function(x){
  if(sum(is.na(x))==length(x)){
    return(FALSE)}
  else {
    x = x[!is.na(x)]
    return(all(x))}}

#reading and cleaning files: 
vcf_files = dir(pattern = ".*\\.vcf")
vcf_list = lapply(vcf_files, function(x) read.delim(x, header=FALSE, comment.char="#")[,c(2,4,5,6,8)])
vcf_list = lapply(vcf_list, function(x) 
  cbind(x[,1:4],data.frame(do.call("rbind", strsplit(as.character(x$V8), ";", fixed = TRUE)))))

vcf_list = lapply(vcf_list, function(y){
  y$X1 = as.numeric(unlist(lapply(strsplit(y$X1,"="), function(x) x[2])))
  return(y)})
vcf_list = lapply(vcf_list, function(y){
  y$X2 = as.numeric(unlist(lapply(strsplit(y$X2,"="), function(x) x[2])))
  return(y)})
vcf_list = lapply(vcf_list, function(y){
  y$X3 = as.numeric(unlist(lapply(strsplit(y$X3,"="), function(x) x[2])))
  return(y)})
vcf_list = lapply(vcf_list, function(x) 
  cbind(x[,1:7],data.frame(do.call("rbind", strsplit(as.character(x$X4), ",", fixed = TRUE)))))
vcf_list = lapply(vcf_list, function(x) cbind(x,paste(x[,2], ">", x[,3], sep="")))
vcf_list = lapply(vcf_list, setNames, c('position','ref_base','alt_base',
                                        'qual','depth','freq','sb',
                                        'ref_fwd_reads','ref_rev_reads','alt_fwd_reads','alt_rev_reads','mut'))
vcf_list = lapply(vcf_list, function(y){
  y$ref_fwd_reads = as.numeric(unlist(lapply(strsplit(y$ref_fwd_reads,"="), function(x) x[2])))
  return(y)})
vcf_list = lapply(vcf_list, function(x) data.table(x))

vcf_list = lapply(vcf_list,function(y){y$strand_bias = apply(y,1,function(x){fisher.test(matrix(c(as.numeric(x[8]),
                                                                                                  as.numeric(x[9]),
                                                                                                  as.numeric(x[10]),
                                                                                                  as.numeric(x[11])),2,2))$p.value > 0.05}); return(y)})

#filtering parameters:
minor_cutoff = 0.5 #exclude any position with ALL samples above this frequency (as consensus)
min_freq = 0.03 #exclude any position with ALL sampled below this frequency
min_depth = 100 #exclude any position that ANY sample has fewer total reads than this
#apply min depth, min freq, and sb filters
vcf_list = lapply(vcf_list,function(x){y = x[x$freq>min_freq]; return(y)})
vcf_list = lapply(vcf_list,function(x){y = x[x$depth>min_depth]; return(y)})
vcf_list = lapply(vcf_list,function(x){y = x[x$strand_bias==TRUE]; return(y)})

#create merged data tables for depth and af (can extend to automatically look for SB and fwd/rev read minimums)
long_depth_data_merged = data.table()
for(i in 1:length(vcf_list)){
  name = strsplit(vcf_files[i],".",2)[[1]][1] 
  current_vcf_long = cbind(vcf_list[[i]][,.(position,depth,mut)], name)
  long_depth_data_merged = rbind(long_depth_data_merged,current_vcf_long)
}
wide_depth_data_merged = spread(long_depth_data_merged,name,depth)
wide_depth_data_merged = as.data.table(wide_depth_data_merged)

long_freq_data_merged = data.table()
for(i in 1:length(vcf_list)){
  name = strsplit(vcf_files[i],".",2)[[1]][1] 
  current_vcf_long = cbind(vcf_list[[i]][,.(position,freq,mut)], name)
  long_freq_data_merged = rbind(long_freq_data_merged,current_vcf_long)
}
wide_freq_data_merged = spread(long_freq_data_merged,name,freq)
wide_freq_data_merged = as.data.table(wide_freq_data_merged)


hh3_ids = masked_most_complete_sample[household==3,broad_id]
hh7_ids = masked_most_complete_sample[household==7,broad_id]
hh9_ids = masked_most_complete_sample[household==9,broad_id]
hh11_ids = masked_most_complete_sample[household==11,broad_id]
hh13_ids = masked_most_complete_sample[household==13,broad_id]
hh15_ids = masked_most_complete_sample[household==15,broad_id]
hh18_ids = masked_most_complete_sample[household==18,broad_id]
hh19_ids = masked_most_complete_sample[household==19,broad_id]
hh21_ids = masked_most_complete_sample[household==21,broad_id]
hh22_ids = masked_most_complete_sample[household==22,broad_id]
hh23_ids = masked_most_complete_sample[household==23,broad_id]
hh25_ids = masked_most_complete_sample[household==25,broad_id]
hh27_ids = masked_most_complete_sample[household==27,broad_id]
hh28_ids = masked_most_complete_sample[household==28,broad_id]
hh36_ids = masked_most_complete_sample[household==36,broad_id]

hh3_ids_isnv = c(c("position","mut"),hh3_ids,paste(hh3_ids,"_depth",sep=""))
hh7_ids_isnv = c(c("position","mut"),hh7_ids,paste(hh7_ids,"_depth",sep=""))
hh9_ids_isnv = c(c("position","mut"),hh9_ids,paste(hh9_ids,"_depth",sep=""))
hh11_ids_isnv = c(c("position","mut"),hh11_ids,paste(hh11_ids,"_depth",sep=""))
hh13_ids_isnv = c(c("position","mut"),hh13_ids,paste(hh13_ids,"_depth",sep=""))
hh15_ids_isnv = c(c("position","mut"),hh15_ids,paste(hh15_ids,"_depth",sep=""))
hh18_ids_isnv = c(c("position","mut"),hh18_ids,paste(hh18_ids,"_depth",sep=""))
hh19_ids_isnv = c(c("position","mut"),hh19_ids,paste(hh19_ids,"_depth",sep=""))
hh21_ids_isnv = c(c("position","mut"),hh21_ids,paste(hh21_ids,"_depth",sep=""))
hh22_ids_isnv = c(c("position","mut"),hh22_ids,paste(hh22_ids,"_depth",sep=""))
hh23_ids_isnv = c(c("position","mut"),hh23_ids,paste(hh23_ids,"_depth",sep=""))
hh25_ids_isnv = c(c("position","mut"),hh25_ids,paste(hh25_ids,"_depth",sep=""))
hh27_ids_isnv = c(c("position","mut"),hh27_ids,paste(hh27_ids,"_depth",sep=""))
hh28_ids_isnv = c(c("position","mut"),hh28_ids,paste(hh28_ids,"_depth",sep=""))
hh36_ids_isnv = c(c("position","mut"),hh36_ids,paste(hh36_ids,"_depth",sep=""))

colnames(wide_depth_data_merged) = c(c("position","mut"),paste(colnames(wide_depth_data_merged)[-c(1,2)],"_depth",sep=""))
wide_isnv_data = wide_freq_data_merged[wide_depth_data_merged,on=.(position,mut)]

hh3_isnv_freq = na.omit(wide_isnv_data[,..hh3_ids_isnv])
hh7_isnv_freq = na.omit(wide_isnv_data[,..hh7_ids_isnv])
hh9_isnv_freq = na.omit(wide_isnv_data[,..hh9_ids_isnv])
hh11_isnv_freq = na.omit(wide_isnv_data[,..hh11_ids_isnv])
hh13_isnv_freq = na.omit(wide_isnv_data[,..hh13_ids_isnv])
hh15_isnv_freq = na.omit(wide_isnv_data[,..hh15_ids_isnv])
hh18_isnv_freq = na.omit(wide_isnv_data[,..hh18_ids_isnv])
hh19_isnv_freq = na.omit(wide_isnv_data[,..hh19_ids_isnv])
hh21_isnv_freq = na.omit(wide_isnv_data[,..hh21_ids_isnv])
hh22_isnv_freq = na.omit(wide_isnv_data[,..hh22_ids_isnv])
hh23_isnv_freq = na.omit(wide_isnv_data[,..hh23_ids_isnv])
hh25_isnv_freq = na.omit(wide_isnv_data[,..hh25_ids_isnv])
hh27_isnv_freq = na.omit(wide_isnv_data[,..hh27_ids_isnv])
hh28_isnv_freq = na.omit(wide_isnv_data[,..hh28_ids_isnv])
hh36_isnv_freq = na.omit(wide_isnv_data[,..hh36_ids_isnv])


#hh3_isnv_table = hh3_isnv_freq[which(hh3_isnv_freq[,3]<0.5 | hh3_isnv_freq[,4]<0.5),]
#hh7_isnv_table = hh7_isnv_freq[which(hh7_isnv_freq[,3]<0.5 | hh7_isnv_freq[,4]<0.5),]
hh9_isnv_table = hh9_isnv_freq[which(hh9_isnv_freq[,3]<0.5 | hh9_isnv_freq[,4]<0.5),]
#hh11_isnv_table = hh11_isnv_freq[which(hh11_isnv_freq[,3]<0.5 | hh11_isnv_freq[,4]<0.5),]
#hh13_isnv_table = hh13_isnv_freq[which(hh13_isnv_freq[,3]<0.5 | hh13_isnv_freq[,4]<0.5),]
#hh15_isnv_table = hh15_isnv_freq[which(hh15_isnv_freq[,3]<0.5 | hh15_isnv_freq[,4]<0.5 | hh15_isnv_freq[,5]<0.5),]
hh18_isnv_table = hh18_isnv_freq[which(hh18_isnv_freq[,3]<0.5 | hh18_isnv_freq[,4]<0.5),]
#hh19_isnv_table = hh19_isnv_freq[which(hh19_isnv_freq[,3]<0.5 | hh19_isnv_freq[,4]<0.5),]
#hh21_isnv_table = hh21_isnv_freq[which(hh21_isnv_freq[,3]<0.5 | hh21_isnv_freq[,4]<0.5 | hh21_isnv_freq[,5]<0.5 | hh21_isnv_freq[,6]<0.5),]
#hh22_isnv_table = hh22_isnv_freq[which(hh22_isnv_freq[,3]<0.5 | hh22_isnv_freq[,4]<0.5 | hh22_isnv_freq[,5]<0.5),]
#hh23_isnv_table = hh23_isnv_freq[which(hh23_isnv_freq[,3]<0.5 | hh23_isnv_freq[,4]<0.5),]
hh25_isnv_table = hh25_isnv_freq[which(hh25_isnv_freq[,3]<0.5 | hh25_isnv_freq[,4]<0.5),]
#hh27_isnv_table = hh27_isnv_freq[which(hh27_isnv_freq[,3]<0.5 | hh27_isnv_freq[,4]<0.5),]
#hh28_isnv_table = hh28_isnv_freq[which(hh28_isnv_freq[,3]<0.5 | hh28_isnv_freq[,4]<0.5),]
#hh36_isnv_table = hh36_isnv_freq[which(hh36_isnv_freq[,3]<0.5 | hh36_isnv_freq[,4]<0.5 | hh36_isnv_freq[,5]<0.5),]

