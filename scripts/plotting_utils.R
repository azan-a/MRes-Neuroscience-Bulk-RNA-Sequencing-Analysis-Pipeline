# Save a plot to a TIFF file with given dimensions and resolution.
save_tiff <- function(
  file_path,
  width = 150,
  height = 84,
  units = "mm",
  res = 300,
  compression = "lzw",
  plot_code
) {
  tiff(
    file_path,
    width = width,
    height = height,
    units = units,
    res = res,
    compression = compression
  )

  plot_code()
  dev.off()
}

# Removes grid in ggplot2 plots.
# Dependencies: ggplot2.
remove_grid <- function(x = TRUE, y = TRUE) {
  p <- ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  if (x) {
    p <- p + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
  }

  if (y) {
    p <- p + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
  }

  p
}
