
source("renv/activate.R")

# automatically bind virtualenv to reticulate : 
if (!require(reticulate)){
  cat("package `reticulate` not installed, Python will be unavailable\n")
} else {
  .POETRY_GET_ENV_CMD <-"poetry env info -p"
  if (0 == system(.POETRY_GET_ENV_CMD)){
    .py_env <- system(.POETRY_GET_ENV_CMD, intern=TRUE)
    cat(glue::glue("Setting RETICULATE_PYTHON to {.py_env}\n"))
    Sys.setenv(RETICULATE_PYTHON=glue::glue("{.py_env}/bin/python"))
    reticulate::use_virtualenv(virtualenv = .py_env, required = TRUE)
  } else {
    cat("\n No poetry virtual envs found!\n")
    cat(" Please specify manually your python installation\n")
    cat(" in order to correctly use reticulate\n")
    cat(" Sys.setenv(RETICULATE_PYTHON='/path/to/python/virtualenv/')\n")
  }
}

