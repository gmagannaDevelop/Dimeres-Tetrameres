
load.if.exists <- function(x){
  #' Load a data file if it exists
  if (file.exists(x)) {
    load(x, envir = .GlobalEnv)
  }
}

load.from.file <- function(var.name, file.path) {
  #' Load a variable directly to global environment
  objs <- load(file.path)
  if (var.name %in% objs) {
    assign(var.name, eval(as.symbol(var.name)), envir = .GlobalEnv)
  }

}

load.from.file2 <- function(var.name, file.path){
  #' explicitly return value of variable stored in filepath
  objs <- load(file.path)
  if (var.name %in% objs){
    return(eval(as.symbol(var.name)))
  } else {
    return(NULL)
  }
}

unhere <- function(filepath){
  #' Simplify an absolute path to parent/filename
  path.vec <- stringr::str_split(filepath, pattern = "/")[[1]]
  file.name <- path.vec[length(path.vec)]
  file.parent <- path.vec[(length(path.vec)-1)]
  glue::glue("{file.parent}/{file.name}")
}


