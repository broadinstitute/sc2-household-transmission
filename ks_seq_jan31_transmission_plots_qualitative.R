library(data.table)
library(readxl)
source("/Users/tbrockfi/Desktop/projects/covid_data_analysis/lhdseq/ks_seq/all_seqs_corrected_sep5/transmission_plots/transmission_plot_master_script.R")
setwd("/Users/tbrockfi/Desktop/projects/covid_data_analysis/lhdseq/ks_seq/ks_seq_dec5_clean/hh_plots")
all_hh_data = data.table(read_excel("/Users/tbrockfi/Desktop/projects/covid_data_analysis/lhdseq/ks_seq/ks_seq_dec5_clean/transmission_plot_data_added_earliest_date_dec5.xlsx")) #adjusting to earliest test date
all_hh_data$hh = as.numeric(unlist(lapply(strsplit(all_hh_data$`PTID:`,"-"), function(x){x[1]})))
all_hh_data$id = gsub(":","",all_hh_data$`PTID:`)
formatted_hh_data = all_hh_data[,c(22,2:21)]

hh3 = data.frame(formatted_hh_data[hh==3][,-21]) 
hh7 = data.frame(formatted_hh_data[hh==7][,-21]) 
hh9 = data.frame(formatted_hh_data[hh==9][,-21]) 
hh11 = data.frame(formatted_hh_data[hh==11][,-21]) 
hh13 = data.frame(formatted_hh_data[hh==13][,-21]) 
hh15 = data.frame(formatted_hh_data[hh==15][,-21]) 
hh18 = data.frame(formatted_hh_data[hh==18][,-21]) 
hh19 = data.frame(formatted_hh_data[hh==19][,-21]) 
hh21 = data.frame(formatted_hh_data[hh==21][,-21]) 
hh22 = data.frame(formatted_hh_data[hh==22][,-21]) 
hh23 = data.frame(formatted_hh_data[hh==23][,-21]) 
hh25 = data.frame(formatted_hh_data[hh==25][,-21]) 
hh27 = data.frame(formatted_hh_data[hh==27][,-21]) 
hh28 = data.frame(formatted_hh_data[hh==28][,-21]) 
hh36 = data.frame(formatted_hh_data[hh==36][,-21]) 
all_color_map = data.frame(type=c("BA.5.5","BA.2.13","BA.2.40.1","BA.2.3","BA.2.48","BA.2.12.1","BA.2","BA.2.23.1","positive","negative"),
                           col=c("#C8D4EB","thistle","#C8D4EB","#C8D4EB","#C8D4EB","#C8D4EB","#C8D4EB","#C8D4EB","#F0F0FA","white")) #old positive color = "#E0EBF6"
all_leg_order = all_color_map$type
hh27_color_map = data.frame(type=c("BA.2.12.1","BA.2","positive","negative"),
                            col=c("#C8D4EB","thistle","#F0F0FA","white"))
hh27_leg_order = hh27_color_map$type

hh9_color_map = data.frame(type=c("BA.2.12.1","BA.2.10","positive","negative"),
                            col=c("#C8D4EB","thistle","#F0F0FA","white"))
hh9_leg_order = hh9_color_map$type


### plot 3
hh3_transmission_data = data.frame(from="003-01",to="003-02",start_day=1,end_day=5)
hh3_plot = longitudinal_transmission_plot(hh3, main_lab = "Household 3",
                                          color_map = all_color_map, legend_order = all_leg_order, transmission_data = hh3_transmission_data) 

### plot 7 
hh7_transmission_data = data.frame(from="007-01",to="007-04",start_day=1,end_day=5)
hh7_order = c("007-01","007-04","007-02","007-03")
hh7_plot = longitudinal_transmission_plot(hh7, main_lab = "Household 7", custom_order = hh7_order,
                                          color_map = all_color_map, legend_order = all_leg_order, transmission_data = hh7_transmission_data) 

### plot 9 
hh9_transmission_data = data.frame(from="009-02",to="009-01",start_day=3,end_day=1)
hh9_plot = longitudinal_transmission_plot_reverse(hh9, main_lab = "Household 9",
                                          color_map = hh9_color_map, legend_order = hh9_leg_order, transmission_data = hh9_transmission_data)

### plot 11
hh11_plot = longitudinal_transmission_plot(hh11, main_lab = "Household 11",
                                           color_map = all_color_map, legend_order = all_leg_order)

### plot 13 
hh13_transmission_data = data.frame(from="013-02",to="013-01",start_day=3,end_day=1)
hh13_plot = longitudinal_transmission_plot_reverse(hh13, main_lab = "Household 13",
                                           color_map = all_color_map, legend_order = all_leg_order, transmission_data = hh13_transmission_data) 

### plot 15
hh15_transmission_data = data.frame(from=c("015-02","015-02","015-04"),
                                    to=c("015-03","015-04","015-03"),
                                    start_day=c(4,4,1),
                                    end_day=c(7,1,7))
hh15_order = c("015-04","015-02","015-03","015-01")
hh15_plot = longitudinal_transmission_plot_reverse(hh15, main_lab = "Household 15", custom_order = hh15_order,
                                           color_map = all_color_map, legend_order = all_leg_order, transmission_data = hh15_transmission_data) 

### plot 18
hh18_transmission_data = data.frame(from="018-03",to="018-01",start_day=2,end_day=1)
hh18_order = c("018-01","018-03","018-02","018-04")
hh18_plot = longitudinal_transmission_plot_reverse(hh18, main_lab = "Household 18", custom_order = hh18_order,
                                           color_map = all_color_map, legend_order = all_leg_order, transmission_data = hh18_transmission_data) 
### plot 19
hh19_order = c("019-02","019-01")
hh19_transmission_data = data.frame(from="019-01",
                                    to="019-02",start_day=3,end_day=1)
hh19_plot = longitudinal_transmission_plot_reverse(hh19, main_lab = "Household 19", custom_order = hh19_order,
                                           color_map = all_color_map, legend_order = all_leg_order, transmission_data = hh19_transmission_data) 

### plot 21 
hh21_transmission_data = data.frame(from=c("021-04","021-02","021-04"),
                                    to=c("021-02","021-03","021-03"),
                                    start_day=c(3,6,3),end_day=c(6,8,8))
hh21_order = c("021-01","021-04","021-02","021-03")
hh21_plot = longitudinal_transmission_plot(hh21, main_lab = "Household 21", custom_order = hh21_order,
                                           color_map = all_color_map, legend_order = all_leg_order, transmission_data = hh21_transmission_data) 

### plot 22 
hh22_transmission_data = data.frame(from=c("022-04"),
                                    to=c("022-05"),
                                    start_day=c(2),end_day=c(7))
hh22_order = c("022-01","022-04","022-05","022-02","022-03")
hh22_plot = longitudinal_transmission_plot(hh22, main_lab = "Household 22", custom_order = hh22_order,
                                           color_map = all_color_map, legend_order = all_leg_order, transmission_data = hh22_transmission_data)

### plot 23
hh23_transmission_data = data.frame(from="023-01",to="023-04",start_day=1,end_day=7)
hh23_order = c("023-01","023-04","023-02","023-03")
hh23_plot = longitudinal_transmission_plot(hh23, main_lab = "Household 23", custom_order = hh23_order,
                                           color_map = all_color_map, legend_order = all_leg_order, transmission_data = hh23_transmission_data)

### plot 25
hh25_plot = longitudinal_transmission_plot(hh25, main_lab = "Household 25", 
                                           color_map = all_color_map, legend_order = all_leg_order)

### plot 27
hh27_order = c("027-01","027-03","027-02")
hh27_plot = longitudinal_transmission_plot(hh27, main_lab = "Household 27", custom_order = hh27_order, 
                                           color_map = hh27_color_map, legend_order = hh27_leg_order)

### plot 28
hh28_transmission_data = data.frame(from="028-01",to="028-02",start_day=1,end_day=6)
hh28_plot = longitudinal_transmission_plot(hh28, main_lab = "Household 28",
                                           color_map = all_color_map, legend_order = all_leg_order, transmission_data = hh28_transmission_data)

### plot36 #to fix
hh36_transmission_data = data.frame(from=c("036-04","036-02","036-02"),
                                    to=c("036-01","036-01","036-04"),
                                    start_day=c(3,3,3),
                                    end_day=c(1,1,3))
hh36_order = c("036-01","036-02","036-04","036-03")
hh36_plot = longitudinal_transmission_plot_reverse(hh36, main_lab = "Household 36", custom_order = hh36_order,
                                           color_map = all_color_map, legend_order = all_leg_order, transmission_data = hh36_transmission_data)



hh3_plot
hh7_plot
hh9_plot
hh11_plot
hh13_plot
hh15_plot
hh18_plot
hh19_plot
hh21_plot
hh22_plot
hh23_plot
hh25_plot
hh27_plot
hh28_plot
hh36_plot 





ggsave(
  filename = paste(deparse(substitute(hh3_plot)),".pdf",sep=""),
  plot = hh3_plot,
  device = "pdf",
  path = "/Users/tbrockfi/Desktop/projects/covid_data_analysis/lhdseq/ks_seq/ks_seq_dec5_clean/hh_plots",
  width = 14,
  height = 4,
  units = "in")
ggsave(
  filename = paste(deparse(substitute(hh7_plot)),".pdf",sep=""),
  plot = hh7_plot,
  device = "pdf",
  path = "/Users/tbrockfi/Desktop/projects/covid_data_analysis/lhdseq/ks_seq/ks_seq_dec5_clean/hh_plots",
  width = 14,
  height = 4,
  units = "in")
ggsave(
  filename = paste(deparse(substitute(hh9_plot)),".pdf",sep=""),
  plot = hh9_plot,
  device = "pdf",
  path = "/Users/tbrockfi/Desktop/projects/covid_data_analysis/lhdseq/ks_seq/ks_seq_dec5_clean/hh_plots",
  width = 14,
  height = 4,
  units = "in")
ggsave(
  filename = paste(deparse(substitute(hh11_plot)),".pdf",sep=""),
  plot = hh11_plot,
  device = "pdf",
  path = "/Users/tbrockfi/Desktop/projects/covid_data_analysis/lhdseq/ks_seq/ks_seq_dec5_clean/hh_plots",
  width = 14,
  height = 4,
  units = "in")
ggsave(
  filename = paste(deparse(substitute(hh13_plot)),".pdf",sep=""),
  plot = hh13_plot,
  device = "pdf",
  path = "/Users/tbrockfi/Desktop/projects/covid_data_analysis/lhdseq/ks_seq/ks_seq_dec5_clean/hh_plots",
  width = 14,
  height = 4,
  units = "in")
ggsave(
  filename = paste(deparse(substitute(hh15_plot)),".pdf",sep=""),
  plot = hh15_plot,
  device = "pdf",
  path = "/Users/tbrockfi/Desktop/projects/covid_data_analysis/lhdseq/ks_seq/ks_seq_dec5_clean/hh_plots",
  width = 14,
  height = 4,
  units = "in")
ggsave(
  filename = paste(deparse(substitute(hh18_plot)),".pdf",sep=""),
  plot = hh18_plot,
  device = "pdf",
  path = "/Users/tbrockfi/Desktop/projects/covid_data_analysis/lhdseq/ks_seq/ks_seq_dec5_clean/hh_plots",
  width = 14,
  height = 4,
  units = "in")
ggsave(
  filename = paste(deparse(substitute(hh19_plot)),".pdf",sep=""),
  plot = hh19_plot,
  device = "pdf",
  path = "/Users/tbrockfi/Desktop/projects/covid_data_analysis/lhdseq/ks_seq/ks_seq_dec5_clean/hh_plots",
  width = 14,
  height = 4,
  units = "in")
ggsave(
  filename = paste(deparse(substitute(hh21_plot)),".pdf",sep=""),
  plot = hh21_plot,
  device = "pdf",
  path = "/Users/tbrockfi/Desktop/projects/covid_data_analysis/lhdseq/ks_seq/ks_seq_dec5_clean/hh_plots",
  width = 14,
  height = 4,
  units = "in")
ggsave(
  filename = paste(deparse(substitute(hh22_plot)),".pdf",sep=""),
  plot = hh22_plot,
  device = "pdf",
  path = "/Users/tbrockfi/Desktop/projects/covid_data_analysis/lhdseq/ks_seq/ks_seq_dec5_clean/hh_plots",
  width = 14,
  height = 4,
  units = "in")
ggsave(
  filename = paste(deparse(substitute(hh23_plot)),".pdf",sep=""),
  plot = hh23_plot,
  device = "pdf",
  path = "/Users/tbrockfi/Desktop/projects/covid_data_analysis/lhdseq/ks_seq/ks_seq_dec5_clean/hh_plots",
  width = 14,
  height = 4,
  units = "in")
ggsave(
  filename = paste(deparse(substitute(hh25_plot)),".pdf",sep=""),
  plot = hh25_plot,
  device = "pdf",
  path = "/Users/tbrockfi/Desktop/projects/covid_data_analysis/lhdseq/ks_seq/ks_seq_dec5_clean/hh_plots",
  width = 14,
  height = 4,
  units = "in")
ggsave(
  filename = paste(deparse(substitute(hh27_plot)),".pdf",sep=""),
  plot = hh27_plot,
  device = "pdf",
  path = "/Users/tbrockfi/Desktop/projects/covid_data_analysis/lhdseq/ks_seq/ks_seq_dec5_clean/hh_plots",
  width = 14,
  height = 4,
  units = "in")
ggsave(
  filename = paste(deparse(substitute(hh28_plot)),".pdf",sep=""),
  plot = hh28_plot,
  device = "pdf",
  path = "/Users/tbrockfi/Desktop/projects/covid_data_analysis/lhdseq/ks_seq/ks_seq_dec5_clean/hh_plots",
  width = 14,
  height = 4,
  units = "in")
ggsave(
  filename = paste(deparse(substitute(hh36_plot)),".pdf",sep=""),
  plot = hh36_plot,
  device = "pdf",
  path = "/Users/tbrockfi/Desktop/projects/covid_data_analysis/lhdseq/ks_seq/ks_seq_dec5_clean/hh_plots",
  width = 14,
  height = 4,
  units = "in")
