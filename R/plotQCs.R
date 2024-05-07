### Plotting function for QC's ###
plotQCs <- function(compound_info,output_directory,rt_list,int_list,x_list,y_list,batchnr,sample_names) {
  # Define colors for rt,int pairs and x,y pairs
  c25 <- palette36.colors(length(rt_list))

  plot_title <- paste("Component:", compound_info$ID)
  sub <- compound_info$NAME
  plot_file <-
    paste("Component_", compound_info$ID, ".png", sep = "")
  if(dir.exists(paste0(output_directory, "QCbatch_", batchnr)) == FALSE){
    dir.create(paste0(output_directory, "QCbatch_", batchnr))
     }
  png(filename = file.path(paste0(output_directory, "QCbatch_", batchnr),
                           plot_file))
  par(mar=c(5.1, 4.1, 4.1, 10.1), xpd= TRUE)
  plot(
    NULL,
    xlim = range(unlist(rt_list), unlist(x_list)),
    ylim = range(unlist(int_list), unlist(y_list)),
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
    abline(v = x[1], xpd = FALSE,col = rt_int_color)
    abline(v = tail(x, 1), xpd = FALSE ,col = rt_int_color)
    points(rt, int, col = rt_int_color)     # Points for rt and int

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


