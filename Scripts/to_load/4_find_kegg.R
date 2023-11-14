# Extract KEGGs from EC numbers ----
find_kegg <- function(x) {
  kegg <- KEGGREST::keggGet(x)
  check <- kegg[[1]]$NAME

  if (stringr::str_detect(
    check[1],
    pattern = "Transferred"
  )) {
    print("replacing faulty EC")

    kegg_name_check <- str_extract(
      check,
      pattern = regex("\\d.*")
    )

    to_replace <- as.character(x)

    smpl_df <- smpl_df %>%
      dplyr::mutate(
        EC_number = stringr::str_replace_all(
          smpl_df$EC_number,
          pattern = to_replace,
          replacement = kegg_name_check
        )
      )
    
    assign("smpl_df", smpl_df, envir = .GlobalEnv)
    x <- kegg_name_check
    kegg <- KEGGREST::keggGet(x)
  }
  pb_assign$tick()
  kegg
}
