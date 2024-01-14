#' Function for saving png of biplot
#' @keywords internal

save_png <- function(code, width = 400, height = 400) {
  path <- tempfile(fileext = ".png")
  set.seed(5)
  png(path, width = width, height = height)
  on.exit(dev.off())
  code
  
  path
}