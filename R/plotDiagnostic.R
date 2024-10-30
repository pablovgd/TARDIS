#' Diagnostic plots
#'
#' Function to make plots of QC samples to screen peak quality and retention
#' time shifts.
#'
#' @param compound_info
#' @param output_directory
#' @param rt_list
#' @param int_list
#' @param x_list
#' @param y_list
#' @param batchnr
#' @param sample_names
#'
#' @import ggplot2
#' @import RColorBrewer


plotDiagnostic <-
  function(compound_info,
           output_directory,
           rt_list,
           int_list,
           x_list,
           y_list,
           batchnr,
           sample_names) {
    # Define colors for rt,int pairs and x,y pairs

    if (length(grep("QC", sample_names)) != 0) {
      #Select only the  QC's for plotting
      rt_list <- rt_list[grep("QC", sample_names)]
      int_list <- int_list[grep("QC", sample_names)]
      x_list <- x_list[grep("QC", sample_names)]
      y_list <- y_list[grep("QC", sample_names)]


      #Select 5 relevant QC's from beginning middle and end of the batch

      if (length(rt_list) > 5) {
        QCs <-
          unique(round(seq(
            from = 1,
            to = length(rt_list),
            length.out = 5
          )))
        rt_list <- rt_list[QCs]
        int_list <- int_list[QCs]
        x_list <- x_list[QCs]
        y_list <- y_list[QCs]
        sample_names <- sample_names[grep("QC", sample_names)]
        sample_names <- sample_names[QCs]
      } else{
        sample_names <- sample_names[grep("QC", sample_names)]
      }


      c25 <- RColorBrewer::brewer.pal(n = length(rt_list), "Set1")

      plot_title <- paste("Component:", compound_info$ID)
      sub <- compound_info$NAME
      plot_file <-
        paste("Component_", compound_info$ID, ".png", sep = "")
      if(dir.exists(paste0(output_directory, "Diagnostic_QCs_Batch_",
                           batchnr)) == FALSE){

        dir.create(paste0(output_directory, "Diagnostic_QCs_Batch_", batchnr))
       }
      png(filename = file.path(paste0(output_directory, "Diagnostic_QCs_Batch_",
                                      batchnr),
                               plot_file))
      par(mar = c(5.1, 4.1, 4.1, 10.1), xpd = TRUE)
      plot(
        NULL,
        xlim = range(unlist(rt_list) / 60, unlist(x_list) / 60),
        ylim = c(0, max(unlist(int_list), unlist(y_list))),
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
        lines(rt / 60, int, col = rt_int_color)      # Line plot for rt and int
        points(rt / 60, int, col = rt_int_color)     # Points for rt and int
        a = x[1] #left integration border
        b = tail(x, 1) #right integration border
        index_a <- which.min(abs(rt - a))
        index_b <- which.min(abs(rt - b))
        polygon(
          c(rt[index_a] / 60, rt[index_a:index_b] / 60, rt[index_b] / 60),
          c(0, int[index_a:index_b], 0),
          col = adjustcolor(rt_int_color, alpha.f = 0.3),
          border = NA
        )

      },
      rt_list,
      int_list,
      x_list,
      y_list,
      rt_int_color = rt_int_colors)
      # x_y_color = x_y_colors)

      legend(
        "topright",
        inset = c(-0.45, 0),
        legend = unlist(sample_names),
        col = rt_int_colors,
        lty = c(1, 1, 1, 1),
        pch = c(1, 16, 1, 16),
        cex = 0.7,
        title = "QC",
        title.cex = 0.8

      )

      dev.off()

    } else {

      chunk_size <- 5
      num_chunks <- ceiling(length(rt_list) / chunk_size)
      for(i in 1:num_chunks){

        start_index <- (i - 1) * chunk_size + 1
        end_index <- min(i * chunk_size, length(rt_list))

        rt_list_chunk  <- rt_list[start_index:end_index]
        int_list_chunk  <- int_list[start_index:end_index]
        x_list_chunk  <- x_list[start_index:end_index]
        y_list_chunk  <- y_list[start_index:end_index]

        c25 <- RColorBrewer::brewer.pal(n = length(rt_list_chunk), "Set1")


        plot_title <- paste("Component:", compound_info$ID)
        sub <- compound_info$NAME
        plot_file <-
          paste("Component_", compound_info$ID, ".png", sep = "")
        if(dir.exists(paste0(output_directory, "evalbatch_", batchnr,i)) ==
           FALSE){
          dir.create(paste0(output_directory, "evalbatch_", batchnr,i))
        }
        png(filename = file.path(paste0(output_directory, "evalbatch_", batchnr,
                                        i),
                                 plot_file))
        par(mar = c(5.1, 4.1, 4.1, 10.1), xpd = TRUE)
        plot(
          NULL,
          xlim = range(unlist(rt_list_chunk) / 60, unlist(x_list_chunk) / 60),
          ylim = c(0, max(unlist(int_list_chunk), unlist(y_list_chunk))),
          type = "n",
          main = plot_title,
          sub = sub,
          xlab = "rt",
          ylab = "int"
        )

        rt_int_colors <-
          c25[1:length(rt_list_chunk)] # You can add more colors as needed
        # You can add more colors as needed
        # Plot lines and points using mapply with colors
        mapply(function(rt,
                        int,
                        x,
                        y,
                        rt_int_color)
          #x_y_color)
        {
          lines(rt / 60, int, col = rt_int_color)      # Line plot for rt and int
          points(rt / 60, int, col = rt_int_color)     # Points for rt and int
          a = x[1] #left integration border
          b = tail(x, 1) #right integration border
          index_a <- which.min(abs(rt - a))
          index_b <- which.min(abs(rt - b))
          polygon(
            c(rt[index_a] / 60, rt[index_a:index_b] / 60, rt[index_b] / 60),
            c(0, int[index_a:index_b], 0),
            col = adjustcolor(rt_int_color, alpha.f = 0.3),
            border = NA
          )

        },
        rt_list_chunk,
        int_list_chunk,
        x_list_chunk,
        y_list_chunk,
        rt_int_color = rt_int_colors)
        # x_y_color = x_y_colors)

        legend(
          "topright",
          inset = c(-0.45, 0),
          legend = unlist(sample_names[start_index:end_index]),
          col = rt_int_colors,
          lty = c(1, 1, 1, 1),
          pch = c(1, 16, 1, 16),
          cex = 0.7,
          title = "QC",
          title.cex = 0.8

        )

        dev.off()

      }

    }

  }
