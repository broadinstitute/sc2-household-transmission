#a function to plot longitudinal testing and sequencing data
#acecpts a dataframe with: case ID (one per row) and a column for each day of testing, result being either: NA, negative, positive, or the name of the lineage
library(ggplot2)
#library(ggpubr)
create_patient_entry = function(pt,start_day,end_day,type,r = 0.2){
  if(start_day==end_day){
    #custom
    offset = 0.48-r
    coords = data.frame(x=c(start_day-offset+r*cos(seq(0.5*pi,1.5*pi,length.out=100)),
                            seq(start_day,(end_day),length.out=100),
                            end_day+offset+r*cos(seq(1.5*pi,2.5*pi,length.out=100)),
                            seq(start_day,(end_day),length.out=100)),
                        y=c(pt+r*sin(seq(0.5*pi,1.5*pi,length.out=100)),
                            rep(pt-r,100),
                            pt+r*sin(seq(1.5*pi,2.5*pi,length.out=100)),
                            rep(pt+r,100)),
                        group=type)
  }
  else{
    offset = 0.48-r
    coords = data.frame(x=c(start_day-offset+r*cos(seq(0.5*pi,1.5*pi,length.out=100)),
                            seq(start_day-offset,(end_day-start_day+offset),length.out=100),
                            end_day+offset+r*cos(seq(1.5*pi,2.5*pi,length.out=100)),
                            seq(start_day-offset,(end_day-start_day+offset),length.out=100)),
                        y=c(pt+r*sin(seq(0.5*pi,1.5*pi,length.out=100)),
                            rep(pt-r,100),
                            pt+r*sin(seq(1.5*pi,2.5*pi,length.out=100)),
                            rep(pt+r,100)),
                        group=type)
  }
  return(coords)
}

create_patient_data = function(pt_testing_data,na_offset = 0){
  #using letters to order factors, only allows 25 break points!
  #also using it for offset indexing, dates longer than 25 will break!
  pt_id = as.numeric(unname(unlist(pt_testing_data[1])))
  pt_data = pt_testing_data[-1]
  break_points = c(1,which(pt_data[-1]!=pt_data[-length(pt_data)])+1)
  pt_plot_data = create_patient_entry(pt_id,break_points[1],length(pt_data),paste(pt_id,letters[na_offset],letters[1],pt_data[break_points[1]],sep="-"))
  for(i in break_points[-1]){
    pt_plot_data = rbind(pt_plot_data,create_patient_entry(pt_id,i,length(pt_data),paste(pt_id,letters[na_offset],letters[i],pt_data[i],sep="-")))
  }
  pt_plot_data$x = pt_plot_data$x+na_offset
  return(pt_plot_data)
}

#other issues: no automatic trimming at beginning
#now transmission arrows defintiely more broken
#also want method to suppress legend

longitudinal_transmission_plot = function(testing_data, transmission_data = NA, custom_order = NA, legend_order = NA,
                                          day_names = NA, main_lab = NA,color_map = NA, text_annotations = NA){
  #testing_data=hh9 #REMOVE
  #day_names=NA #REMOVE
  #transmission_data = hh9_transmission_data #remove
  first_day_all_negative = min(dim(testing_data)[2],
                               dim(testing_data)[2]-unname(which.max(rev(!apply(testing_data,2,function(x){all(x=="negative",na.rm = TRUE)}))))+2)
  #testing_data = testing_data[,1:first_day_all_negative]
  
  first_day_data = unname(which.min(apply(testing_data[,-1],2, function(x){all(is.na(x))})))+1
  
  testing_data = testing_data[,c(1,first_day_data:first_day_all_negative)]
  if(!any(is.na(custom_order))){
    testing_data = testing_data[match(custom_order,testing_data$id),] #nede to check if NA
  }
  
  case_list = rev(testing_data[,1])
  testing_data[,1] = rev(seq(1,length(case_list))) #maintain drawing order
  case_names = paste("Pt",case_list,sep=" ")
  day_sequence = colnames(testing_data)[-1]
  testing_data = testing_data[order(testing_data$id,decreasing = FALSE),] #reverse plotting order
  if(is.na(day_names)){
    day_names = paste("Day",1:length(day_sequence),sep=" ")
  }
  trimmed_testing_data = testing_data[,-1]
  output_plot = ggplot(data.frame()) + geom_blank() + theme_bw() +
    scale_x_continuous(breaks=seq(1,length(day_sequence),1), limits=c(0.5,length(day_sequence)+0.5),labels=day_names) +
    scale_y_continuous(breaks=seq(1,length(case_list),1), limits=c(0.5,length(case_list)+2.5),labels=case_names) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid = element_line(color = "grey",
                                    linewidth = 0.75,
                                    linetype = 2),
          panel.border = element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y = element_text(margin = margin(r = -15)),
          axis.text.x = element_text(margin = margin(t = -7))) + 
    coord_fixed(ratio = 1, xlim=c(0.5,length(day_sequence)+0.5), ylim=c(0.5,length(case_list)+1)) 
  if(any(is.na(testing_data))){ #should check that no rows are all NA
    element_types = c()
    pt_list = list()
    for(i in 1:dim(testing_data)[1]){
      na_breaks = which(is.na(testing_data[i,]))
      if(length(na_breaks)==0){
        pt_list = c(pt_list,list(create_patient_data(testing_data[i,])))
        element_types = c(unlist(lapply(unlist(lapply(pt_list,function(x){unique(x[,3])})), function(x){strsplit(x,"-")[[1]][4]})),element_types)
      }
      else{
        exclude_na_portions_table = data.frame(interval_start=c(2,na_breaks+1),
                                               interval_end=c(na_breaks-1,dim(testing_data)[2]),
                                               na_offset=c(0,na_breaks-1))
        exclude_na_portions_table = exclude_na_portions_table[!exclude_na_portions_table$interval_start >
                                                                exclude_na_portions_table$interval_end,]
        pt_list = c(pt_list, apply(exclude_na_portions_table,1,
                                   function(x){create_patient_data(testing_data[i,c(1,x[1]:x[2])],x[3])}))
        element_types = c(unlist(lapply(unlist(lapply(pt_list,function(x){unique(x[,3])})), function(x){strsplit(x,"-")[[1]][4]})),element_types)
      }
    }
    #pt_list  = rev(pt_list)
    num_elements = sum(unlist(lapply(pt_list,function(x){length(unique(x[,3]))})))
    outline_vec = rep("darkgrey",num_elements)
  }

  else{ #no NA in data, easy case
    pt_list = apply(testing_data,1,create_patient_data)
    #pt_list = rev(pt_list)
    num_elements = sum(unlist(lapply(pt_list,function(x){length(unique(x[,3]))})))
    outline_vec = rep("darkgrey",num_elements)
    element_types = unlist(lapply(unlist(lapply(pt_list,function(x){unique(x[,3])})), function(x){strsplit(x,"-")[[1]][4]}))
  }
  for(pt in pt_list){
    output_plot = output_plot + geom_polygon(data=pt, aes(x=x,y=y,col=factor(group),fill=factor(group)))
  }
  element_types = unname(element_types)

  if(any(is.na(color_map))){
    color_map=data.frame(type=unique(element_types),col=palette(rainbow(length(unique(element_types)))))
  }
  
  #would like to allow for legend ordering as well
  
  legend_text = sort(unique(element_types))
  if(!any(is.na(legend_order))){
    legend_text = legend_text[order(match(legend_text,legend_order))]
  }
  element_types = color_map$col[match(element_types,color_map$type)]
  
  legend_fills = color_map$col[match(legend_text,color_map$type)] 
  legend_cols = color_map$col[match(legend_text,color_map$type)] 
  legend_cols[which(legend_cols=="white")]="darkgrey"
  
  legend_ploygon = data.frame(x=c(1,length(day_sequence),length(day_sequence),1),
                              y=c(length(case_list)+0.4,
                                  length(case_list)+0.4,
                                  length(case_list)+2.5,
                                  length(case_list)+2.5),
                              group="legend")
  legend_start_pos = length(day_sequence)/2 - length(legend_text)/2 + 0.5
  legend_pos_list = seq(legend_start_pos,legend_start_pos+length(legend_text),length.out=length(legend_text))
  output_plot = output_plot + scale_color_manual(values=outline_vec) + scale_fill_manual(values=element_types) +
    geom_vline(xintercept = seq(1.5,length(day_sequence)-0.5,1), linetype="dashed",col="darkgrey") +
    geom_hline(yintercept = 0.5, linetype="dashed",col="grey") +
    theme(axis.title.y = element_blank(),axis.title.x = element_blank(),legend.position = "none") +
    theme(text=element_text(size=14)) + geom_polygon(data=legend_ploygon, aes(x=x,y=y),inherit.aes=FALSE,fill="white") 
  for(i in 1:length(legend_pos_list)){
    output_plot = output_plot + geom_polygon(data=create_patient_entry(length(case_list)+0.7,legend_pos_list[i],legend_pos_list[i],"test",0.15),
                                             aes(x=x,y=y),
                                             inherit.aes=FALSE,
                                             fill=legend_fills[i],
                                             color=legend_cols[i])
  } 
  output_plot = output_plot + annotate("text", x = legend_pos_list, y = length(case_list)+0.7, label = legend_text)
  if(!is.na(main_lab)){
    output_plot = output_plot + annotate("text", x = 1, 
                                         y = length(case_list)+0.7, label = main_lab,
                                         size=6)
  }
  if(!all(is.na(transmission_data))){ #I think if I convert the path into two segments, I might be able to pass the whole df at once
    transmission_data = data.table(transmission_data)
    transmission_data$reverse = transmission_data$start_day > transmission_data$end_day
    transmission_data[reverse==TRUE, c("from","to","start_day", "end_day") := .(to,from,end_day,start_day)]
    
    transmission_data$from_offset = 0 
    transmission_data$to_offset = 0
    transmission_data$end_day_offset = 0
    master_offset = 0.15
    
    transmission_data$single_day = FALSE
    transmission_data[transmission_data$start_day==transmission_data$end_day,]$single_day = TRUE
    transmission_data[single_day==FALSE,]$end_day = transmission_data[single_day==FALSE,]$end_day-1
    
    from_dupes = names(table(transmission_data$from)[table(transmission_data$from)>1])
    to_dupes = names(table(transmission_data$to)[table(transmission_data$to)>1])
    end_day_dupes = names(table(transmission_data$end_day)[table(transmission_data$end_day)>1])
    
    for(fd in from_dupes){
      num_dupes = unname(table(transmission_data$from)[fd])
      offset_seq = seq(-master_offset/num_dupes,master_offset/num_dupes,length.out = num_dupes)
      transmission_data[transmission_data$from==fd,]$from_offset = offset_seq
    }
    
    for(td in to_dupes){
      num_dupes = unname(table(transmission_data$to)[td])
      offset_seq = seq(-master_offset/num_dupes,master_offset/num_dupes,length.out = num_dupes)
      transmission_data[transmission_data$to==td,]$to_offset = offset_seq
    }
    
    for(ed in end_day_dupes){
      num_dupes = unname(table(transmission_data$end_day)[ed])
      offset_seq = seq(-master_offset/num_dupes,master_offset/num_dupes,length.out = num_dupes)
      transmission_data[transmission_data$end_day==ed,]$end_day_offset = offset_seq
    }
    
    single_day_transmissions = transmission_data[single_day==TRUE,]
    one_day_transmissions = transmission_data[transmission_data$start_day==transmission_data$end_day & single_day==FALSE,]
    reverse_transmissions = transmission_data[transmission_data$start_day!=transmission_data$end_day & reverse==TRUE,]
    transmission_data = transmission_data[single_day==FALSE & reverse==FALSE & transmission_data$start_day!=transmission_data$end_day,]
    
    output_plot = output_plot + 
      geom_segment(data = one_day_transmissions[reverse==FALSE],
                   aes(x=end_day+end_day_offset,
                       xend = end_day+end_day_offset,
                       y=match(from,case_list)+from_offset, #returns incorrect length when same case in here multiple times
                       yend = match(to,case_list)+to_offset), linewidth=0.8, lineend = "round") +
      geom_segment(data = one_day_transmissions[reverse==FALSE],
                   aes(x=end_day+end_day_offset,
                       xend = end_day+0.5,
                       y=match(to,case_list)+to_offset,
                       yend = match(to,case_list)+to_offset), linewidth=0.8,
                   arrow = arrow(length = unit(0.1, "inches")),lineend = "round", linejoin = "bevel")
    
    output_plot = output_plot + 
      geom_segment(data = one_day_transmissions[reverse==TRUE],
                   aes(xend=end_day+end_day_offset,
                       x = end_day+end_day_offset,
                       yend=match(from,case_list)+from_offset, #returns incorrect length when same case in here multiple times
                       y = match(to,case_list)+to_offset), linewidth=0.8,
                   arrow = arrow(length = unit(0.1, "inches")),lineend = "round", linejoin = "bevel") + 
      geom_segment(data = one_day_transmissions[reverse==TRUE],
                   aes(x=end_day+end_day_offset,
                       xend = end_day+0.5,
                       y=match(to,case_list)+to_offset,
                       yend = match(to,case_list)+to_offset), linewidth=0.8)
    
    output_plot = output_plot + 
      geom_segment(data = single_day_transmissions,
                   aes(x=end_day+end_day_offset*1.1,
                       xend = end_day+end_day_offset*1.1,
                       y=match(from,case_list), #returns incorrect length when same case in here multiple times
                       yend = match(to,case_list)), linewidth=0.8,
                   arrow = arrow(length = unit(0.1, "inches")),lineend = "round", linejoin = "bevel")
    
    output_plot = output_plot + 
      geom_segment(data = transmission_data, 
                   aes(x = start_day+0.25,
                       xend = end_day+end_day_offset, 
                       y = match(from,case_list)+from_offset,
                       yend = match(from,case_list)+from_offset), linewidth=0.8, lineend = "round") + 
      geom_segment(data = transmission_data,
                   aes(x=end_day+end_day_offset,
                       xend = end_day+end_day_offset,
                       y=match(from,case_list)+from_offset, #returns incorrect length when same case in here multiple times
                       yend = match(to,case_list)+to_offset), linewidth=0.8, lineend = "round") +
      geom_segment(data = transmission_data,
                   aes(x=end_day+end_day_offset,
                       xend = end_day+0.5,
                       y=match(to,case_list)+to_offset,
                       yend = match(to,case_list)+to_offset), linewidth=0.8,
                   arrow = arrow(length = unit(0.1, "inches")),lineend = "round", linejoin = "bevel")
    
    output_plot = output_plot + 
      geom_segment(data = reverse_transmissions, 
                   aes(x = end_day+end_day_offset, 
                       xend = start_day+0.25,
                       y = match(from,case_list)+from_offset,
                       yend = match(from,case_list)+from_offset), linewidth=0.8, lineend = "round",
                   arrow = arrow(length = unit(0.1, "inches"))) + 
      geom_segment(data = reverse_transmissions,
                   aes(x=end_day+end_day_offset,
                       xend = end_day+end_day_offset,
                       y=match(from,case_list)+from_offset, #returns incorrect length when same case in here multiple times
                       yend = match(to,case_list)+to_offset), linewidth=0.8, lineend = "round") +
      geom_segment(data = reverse_transmissions,
                   aes(x=end_day+end_day_offset,
                       xend = end_day+0.5,
                       y=match(to,case_list)+to_offset,
                       yend = match(to,case_list)+to_offset), linewidth=0.8,
                   lineend = "round", linejoin = "bevel")
  }
  
  if(!all(is.na(text_annotations))){
    output_plot = output_plot + geom_label(aes(x = text_annotations$x, y = text_annotations$y, 
                                               label = text_annotations$text), fill = "white",size=3)
    
    #annotate("text", x = text_annotations$x, 
    #                                     y = text_annotations$y, label = text_annotations$text,
    #                                     size=3)
  }
  return(output_plot)
}

