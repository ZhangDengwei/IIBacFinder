if (!require("ggplot2", quietly = TRUE)){
  install.packages("ggplot2", repos = "http://cran.us.r-project.org")
  require("ggplot2", quietly = TRUE)
}
if (!require("gggenes", quietly = TRUE)){
  install.packages("gggenes", repos = "http://cran.us.r-project.org")
  require("gggenes", quietly = TRUE)
}
if (!require("dplyr", quietly = TRUE)){
  install.packages("dplyr", repos = "http://cran.us.r-project.org")
  require("dplyr", quietly = TRUE)
}
if (!require("cowplot", quietly = TRUE)){
  install.packages("cowplot", repos = "http://cran.us.r-project.org")
  require("cowplot", quietly = TRUE)
}
if (!require("ggiraphExtra", quietly = TRUE)){
  install.packages("ggiraphExtra", repos = "http://cran.us.r-project.org")
  require("ggiraphExtra", quietly = TRUE)
}
if (!require("lattice", quietly = TRUE)){
  install.packages("lattice", repos = "http://cran.us.r-project.org")
  require("lattice", quietly = TRUE)
}
if (!require("ggrepel", quietly = TRUE)){
  install.packages("ggrepel", repos = "http://cran.us.r-project.org")
  require("ggrepel", quietly = TRUE)
}
if (!require("gridExtra", quietly = TRUE)){
  install.packages("gridExtra", repos = "http://cran.us.r-project.org")
  require("gridExtra", quietly = TRUE)
}


f_plot <- function(in_file, out_file){
  df_region <- read.delim2(in_file)
  df_region$seq_domain_start <- df_region$Start + df_region$Domain_start - 1
  df_region$seq_domain_end <- df_region$Start + df_region$Domain_end- 1
  
  df_region$middle <- (df_region$seq_domain_start + df_region$seq_domain_end) / 2
  
  #------------- Specify different color
  color_value <- vector()
  
  for (x in df_region$CDs_Class){
    if (x == "precursor"){
      color_value <- append(color_value, c("precursor"="#FF0000"))
    }else if(x == "immunity"){
      color_value <- append(color_value, c("immunity"="#FFA500"))
    }else if(x == "regulator"){
      color_value <- append(color_value, c("regulator"="#008000"))
    }else if(x == "peptidase/transporter"){
      color_value <- append(color_value, c("peptidase/transporter"="#0000FF"))
    }else if(x == "peptidase"){
      color_value <- append(color_value, c("peptidase"="#00CCCC"))
    }else if(x == "transporter"){
      color_value <- append(color_value, c("transporter"="#FF00FF"))
    }else if(x == "others"){
      color_value <- append(color_value, c("others"="#FFFFFF"))
    }
  }
  
  if (sum(!is.na(df_region$Pfam_domain))>0){
    color_value <- append(color_value, c("domain"="#000000"))
  }
  
  #------------- Plot gene cluster
  p1 <- ggplot(df_region, 
               aes(xmin = Start, xmax = End, x = middle, y = Region, 
                   fill = CDs_Class, forward=Orientation)) +
    facet_wrap(~ Region, scales = "free", ncol = 1) +
    geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"))+
    # geom_text(vjust = -0.7, size = 4) +
    ggrepel::geom_text_repel(df_region, 
                             mapping=aes(x = middle, y = Region, label = Domain_index),
                             inherit.aes = F, nudge_y = 1, lineheight = 0.1)+
    geom_subgene_arrow(
      data = subset(df_region, !is.na(df_region$Domain_index)),
      aes(xsubmin = seq_domain_start, xsubmax = seq_domain_end, fill="domain")
    )+
    scale_fill_brewer(palette = "Set3") +
    #theme_genes() +
    theme(legend.position = "right")+
    labs(x=paste("Contig: ", unique(df_region$Contig), sep=""), 
         y="")+
    theme(axis.text.y = element_blank(),
          axis.line.x = element_line(colour = "black"))+
    theme(panel.grid.major.y = element_line(colour = "black"))+
    scale_fill_manual(values=color_value)
  
  #------------- Plot description
  df_region_plot <- df_region %>% select(c("CDs", "Domain_index", "Partial_index", "Domain_description", "CDs_Class"))
  tt <- gridExtra::ttheme_default(colhead=list(fg_params = list(parse=TRUE)),
                                  base_size = 10,
                                  padding = unit(c(2, 4), "mm"))
  p2 <- gridExtra::tableGrob(df_region_plot, rows=NULL, theme = tt)
  
  #------------- Combine two figures
  q <- cowplot::plot_grid(p1, p2, ncol = 1, rel_heights=c(2,14), align = "v", axis = "tb")
  
  svg(file = out_file, height = 16,width = 12)
  print(q)
  dev.off()
}



args <- commandArgs(T)

f_plot(args[1], args[2])

