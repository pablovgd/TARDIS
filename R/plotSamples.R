plotSamples <- function(compound_info, output_directory,rt_list,int_list,x_list,y_list,batchnr,sample_names){


c25 <- c(
  "dodgerblue2",
  "#E31A1C",
  # red
  "green4",
  "#6A3D9A",
  # purple
  "#FF7F00",
  # orange
  "black",
  "gold1",
  "skyblue2",
  "#FB9A99",
  # lt pink
  "palegreen2",
  "#CAB2D6",
  # lt purple
  "#FDBF6F",
  # lt orange
  "gray70",
  "khaki2",
  "maroon",
  "orchid1",
  "deeppink1",
  "blue1",
  "steelblue4",
  "darkturquoise",
  "green1",
  "yellow4",
  "yellow3",
  "darkorange4",
  "brown"
)
plot_title <- paste("Component:", compound_info$ID)
sub <- compound_info$NAME
plot_file <-
  paste("Component_", compound_info$ID, ".png", sep = "")
if(dir.exists(paste0(output_directory, "Samplebatch_", batchnr)) == FALSE){
  dir.create(paste0(output_directory, "Samplebatch_", batchnr))
}
png(filename = file.path(
  paste0(output_directory, "Samplebatch_", batchnr),
  plot_file
))
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
  lines(rt,
        int)#col = rt_int_color)      # Line plot for rt and int
  #abline(v = x[1], col = rt_int_color)
  #abline(v = tail(x, 1), col = rt_int_color)
  #points(rt, int, col = rt_int_color)     # Points for rt and int

},
rt_list,
int_list,
x_list,
y_list,
rt_int_color = rt_int_colors)
# x_y_color = x_y_colors)

# legend(
#   "topright",
#   legend = unlist(sample_names),
#   col = c(rt_int_colors[1], rt_int_colors[2]),
#   lty = c(1, 1, 1, 1),
#   pch = c(1, 16, 1, 16),
#   cex = 0.8,
#   title = "Sample"
# )

dev.off()

}
