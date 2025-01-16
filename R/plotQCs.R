#' Function to make plots of QC samples
#'
#' @param compound_info info for target compounds like name, id, ...
#' @param output_directory `character(1)` directory to save plots.
#' @param rt_list `list` list with retention times of complete window.
#' @param int_list `list` list with intensities in complete window.
#' @param x_list `list` list with retention times of integrated peak.
#' @param y_list `list` list with intensities of integrated peak.
#' @param batchnr `numeric(1)` batch number.
#' @param sample_names names of samples.
#'
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#'
#' @author Pablo Vangeenderhuysen
#' @export

plotQCs <- function(compound_info, output_directory, rt_list, int_list, x_list,
                    y_list, batchnr, sample_names) {
  c25 <- brewer.pal(n = 5, "Set1")
  # Create directory if it doesn't exist
  batch_dir <- paste0(output_directory, "QCbatch_", batchnr)
  if (!dir.exists(batch_dir)) {
    dir.create(batch_dir)
  }
  # Split the data into groups of 5
  n <- length(rt_list)
  groups <- ceiling(n / 5)
  for (i in seq_len(groups)) {
    start_idx <- (i - 1) * 5 + 1
    end_idx <- min(i * 5, n)
    rt_subset <- rt_list[start_idx:end_idx]
    int_subset <- int_list[start_idx:end_idx]
    x_subset <- x_list[start_idx:end_idx]
    y_subset <- y_list[start_idx:end_idx]
    sample_names_subset <- sample_names[start_idx:end_idx]
    # Create file name for the plot
    plot_file <- paste("Component_", compound_info$ID, "_part_", i,
                       ".png", sep = "")
    # Define the output file path
    png(filename = file.path(batch_dir, plot_file))
    # Set up the plot with defined limits and labels
    plot(
      NULL,
      xlim = range(c(unlist(rt_subset) / 60, unlist(x_subset) / 60)),
      ylim = range(c(unlist(int_subset), unlist(y_subset))),
      type = "n",
      main = paste("Component:", compound_info$ID, "- Part", i),
      sub = compound_info$NAME,
      xlab = "Retention Time (minutes)",
      ylab = "Intensity"
    )
    # Plot lines and points using mapply with colors
    mapply(function(rt, int,x,y, rt_int_color) {
      lines(rt / 60, int, col = rt_int_color) # Line plot for rt and int
      a = x[1] #left integration border
      b = tail(x, 1) #right integration border
      index_a <- which.min(abs(rt - a))
      index_b <- which.min(abs(rt - b))
      polygon(
        c(rt[index_a]/60, rt[index_a:index_b]/60, rt[index_b]/60),
        c(0, int[index_a:index_b], 0),
        col = adjustcolor(rt_int_color, alpha.f = 0.3),
        border = NA)
    }, rt_subset, int_subset, x_subset,y_subset,
    rt_int_color = c25[1:length(rt_subset)])
    # Add a legend to the plot
    legend(
      "topright",
      legend = sample_names_subset,
      col = c25[1:length(rt_subset)],
      lty = 1,
      pch = 1,
      cex = 0.8,
      title = "QC"
    )
    # Save the plot to file
    dev.off()
  }
}
