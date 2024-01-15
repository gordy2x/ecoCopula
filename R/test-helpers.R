#' Function for saving png of biplot
#' @keywords internal

save_png <- function(code, width = 400, height = 400) {
  path <- tempfile(fileext = ".png")
  set.seed(5)
  grDevices::png(path, width = width, height = height)
  on.exit(grDevices::dev.off())
  code
  
  path
}

#' Output for tests 
load("data/fixtures.RData")