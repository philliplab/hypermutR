#' Restart R
#'
#' Useful for developing with the vim-R plugin
#' From https://stackoverflow.com/questions/6313079/quit-and-restart-a-clean-r-session-from-within-r

restart_r <- function(status = 0, debug = TRUE) {
  if (FALSE){
    restart_r()
    library(devtools); load_all()
  }
  if (debug) message("restart_r(): Customizing .Last() to relaunch R ...")
  assign(".Last", function() {
    args <- commandArgs()
    system2(args[1], args = args[-1])
  }, envir = globalenv())   
  if (debug) message("restart_r(): Quitting current R session and starting a new one ...")
  quit(save = "no", status = status, runLast = TRUE)
}
