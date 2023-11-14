# 7. PERMANOVA helpers

# cleaner
adonis_clean <- function(x) {
  x <- x %>%
    tibble::rownames_to_column(var = "Var") %>%
    dplyr::rename_with(~ paste("Pr_F_stat"), starts_with("Pr")) %>%
    dplyr::rename_with(~ paste("F_statistic"), starts_with("F")) %>%
    dplyr::mutate(
      Sigif = case_when(
        (Pr_F_stat > 0 & Pr_F_stat <= 0.001) ~ "***",
        (Pr_F_stat > 0.001 & Pr_F_stat <= 0.01) ~ "**",
        (Pr_F_stat > 0.01 & Pr_F_stat <= 0.05) ~ "*",
        (Pr_F_stat > 0.05 & Pr_F_stat <= 0.1) ~ "."
      ),
      Sigif = tidyr::replace_na(.data$Sigif, replace = ""),
      Pr_F_stat = as.numeric(.data$Pr_F_stat)
    ) %>%
    dplyr::filter(Pr_F_stat > 0) %>%
    dplyr::arrange(Pr_F_stat)
}

# runner
PRMNV_run <- function(
    df, 
    int_vat, 
    var_nm, 
    mthd = "bray",
    permt = 999, 
    strt = NULL, 
    strt_nm = NULL
    ) {
  set.seed(26)
  
  # display the variable being tested
  print(var_nm)
  
  # remove NAs from input data
  data_in <- df %>%
    dplyr::filter(rownames(.) %in% int_vat$Sample) %>% 
    as.matrix()
  
  if (!is.null(strt)) {
  strt <- strt %>% 
    dplyr::filter(Sample %in% int_vat$Sample)
  }
  
  # create formula
  frm <- as.formula(paste0("data_in ~ ", var_nm))
  
  # run adonis
  adonis_hold <- vegan::adonis2(
    formula = frm,
    data = int_vat,
    method = mthd,
    strata = strt[[strt_nm]],
    permutations = permt,
    parallel = 5
  )

  adonis_hold
}


# factorise
fctrs <- function(
    df, 
    var_nm, 
    excl
    ) {
  if (!var_nm %in% excl$Var) {
    var <- sym(var_nm)
    
    out <- df %>% 
      dplyr::mutate(!!var := factor(!!var, levels = unique(sort(!!var))))
  } else {
    out <- df
  }
  
  out
}
