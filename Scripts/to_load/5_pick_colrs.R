# Pick custom colours based in levels
pick_colrs <- function(data, n = 5, end = "Others") {
  nlev <- nlevels(data)
  
  t <- tribble(
    ~palette, ~rev,
    "Mako", "frnt", 
    "Reds", "back",
    "Emrld", "back", 
    "YlOrBr", "back",
    "Magenta", "frnt",
    "Blues 2", "frnt",
    "Earth", "back"
  )
  
  n = n
  mycols <- c()
  col <- 1
  
  while (nlev > 0) {
    range = n + 1
    
    frnt <- c(2:range)
    back <- c(1:(range - 1))
    
    ext <- grDevices::hcl.colors(
      palette = paste0(t$palette[col]), 
      n = range)[eval(parse(text = t$rev[col]))]
    
    mycols <- append(mycols, ext)
    
    nlev <- nlev - n
    col <- col + 1
    if (nlev > 0) {
      if (nlev < 5 | col == nrow(t)) {
        n <- nlev
      }
    }
  }
  
  if (!is.na(end)) {
    mycols <- mycols %>%
      purrr::set_names(levels(data)) %>%
      base::replace(end, scico::scico(1, 
                                      palette = "grayC", 
                                      begin = 0.65))
  } else {
    mycols <- mycols %>%
      purrr::set_names(levels(data))
  }
}
