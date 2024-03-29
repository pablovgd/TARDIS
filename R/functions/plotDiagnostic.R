### Plotting function for QC's ###
plotDiagnostic <- function(compound_info,output_directory,rt_list,int_list,x_list,y_list,batchnr,sample_names) {
  # Define colors for rt,int pairs and x,y pairs
  
  #Select only the  QC's for plotting 
  rt_list <- rt_list[grep("QC",sample_names)]
  int_list <- int_list[grep("QC",sample_names)]
  x_list <- x_list[grep("QC",sample_names)]
  y_list <- y_list[grep("QC",sample_names)]
  
  
  #Select 5 relevant QC's from beginning middle and end of the batch
  
  if(length(rt_list) > 5){
    QCs <- unique(round(seq(from = 1,to = length(rt_list),length.out =5)))
    rt_list <- rt_list[QCs]
    int_list <- int_list[QCs]
    x_list <- x_list[QCs]
    y_list <- y_list[QCs]
    sample_names <- sample_names[grep("QC",sample_names)]
    sample_names <- sample_names[QCs]
  } else{
    sample_names <- sample_names[grep("QC",sample_names)]
  }

  
  c25 <- palette36.colors(length(rt_list))
  
  plot_title <- paste("Component:", compound_info$ID)
  sub <- compound_info$NAME
  plot_file <-
    paste("Component_", compound_info$ID, ".png", sep = "")
  dir.create(paste0(output_directory, "QCbatch_", batchnr))
  png(filename = file.path(paste0(output_directory, "QCbatch_", batchnr),
                           plot_file))
  par(mar=c(5.1, 4.1, 4.1, 10.1), xpd= TRUE)
  plot(
    NULL,
    xlim = range(unlist(rt_list), unlist(x_list)),
    ylim = c(0,max(unlist(int_list), unlist(y_list))),
    type = "n",
    main = plot_title,
    sub = sub,
    xlab = "rt",
    ylab = "int"
  )
  
  rt_int_colors <-
    c25[1:length(rt_list)] # You can add more colors as needed
  # You can add more colors as needed
  # Plot lines and points using mapply with colors
  mapply(function(rt,
                  int,
                  x,
                  y,
                  rt_int_color)
    #x_y_color)
  {
    lines(rt, int, col = rt_int_color)      # Line plot for rt and int
    points(rt, int, col = rt_int_color)     # Points for rt and int
    a = x[1] #left integration border
    b = tail(x,1) #right integration border
    index_a <- which.min(abs(rt - a))
    index_b <- which.min(abs(rt - b))
    polygon(c(rt[index_a], rt[index_a:index_b], rt[index_b]), c(0, int[index_a:index_b], 0), col = adjustcolor(rt_int_color,alpha.f = 0.3), border = NA)
    
  },
  rt_list,
  int_list,
  x_list,
  y_list,
  rt_int_color = rt_int_colors)
  # x_y_color = x_y_colors)
  
  legend(
    "topright",
    inset=c(-0.45, 0),
    legend = unlist(sample_names),
    col = rt_int_colors,
    lty = c(1, 1, 1, 1),
    pch = c(1, 16, 1, 16),
    cex = 0.7,
    title = "Sample",
    title.cex = 0.8
    
  )
  
  dev.off()
}


