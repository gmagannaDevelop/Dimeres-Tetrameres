
library(purrr)

message("This module depends on the following packages :\n\n",
        "\t\"purrr\"\n", 
        "\t\"tibble\"\n",
        "\t\"dplyr\"\n")

describe <- function(df){
  #' Create printable description table
  # Verifier si toutes les colonnes sont des nombres :
  if ( !all(apply(df, 2, is.numeric)) ){
    message("All columns must be numeric, aborting")
    return(NULL)
  } else {
    df.summary <- as.data.frame(apply(df, 2, summary))
    df.rows <- rownames(df.summary)
    tibble.summary <- df.summary %>% tibble::tibble() %>% 
      dplyr::mutate(Stat = df.rows) %>%
        dplyr::select(Stat, everything())
    return(tibble.summary)
  }
}

columnwise.t.test <- function(df){
  if ( !all(apply(df, 2, is.numeric)) ){
    message("All columns must be numeric, aborting")
    return(NULL)
  } else {
    cols <- colnames(df)
    t.confint <- function(x){ as.numeric(x$conf.int) }
    t.pval <- function(x){ x$p.value }
    t.estimate <- function(x){ x$estimate }
    
    tests.ls <- apply(df, 2, t.test)
    resume.tib <- t(as.data.frame(sapply(tests.ls, t.confint, simplify = F)))
    resume.tib <- cbind(
      resume.tib, t(as.data.frame(sapply(tests.ls, t.pval, simplify = F))),
      t(as.data.frame(sapply(tests.ls, t.estimate, simplify = F)))
    ) 
    colnames(resume.tib) <- c("2.5%", "92.5%", "p-value", "mean")
    resume.tib
  }
  
}


table.lm <- function(lm.summary){
  #' Create printable lm summary
  lm.table <- lm.summary$coefficients %>% as.data.frame()
  lm.table.rows <- rownames(lm.table)
  lm.table %>% tibble::tibble() %>% 
    dplyr::mutate(Term = lm.table.rows) %>%
      dplyr::select(Term, everything())
}

confint.table <- function(model, alpha=0.05){
  #' Create printable confint
  original.conf <- confint(model, level = 1 - alpha) %>% as.data.frame()
  conf.rows <- rownames(original.conf)
  original.conf %>% tibble::tibble %>% 
    dplyr::mutate(Term = conf.rows) %>%
      dplyr::select(Term, everything())
}


