# HEMI - data preparation

# Session parameters and packages ----
pkgs_to_load <- list.files(
  here::here(
    "Scripts",
    "to_load"
  ),
  pattern = "^1_|^6_",
  full.names = T
)

pkgs_to_load <- purrr::map(
  pkgs_to_load,
  source
)

## define major directory paths ----
path_to_anlys <- here::here(
  "Data",
  "Analysis"
)

remake_trees <- FALSE


# Load data ----
ASV_runs <- sort(list.files(
  path = path_to_anlys,
  pattern = "_seqtab",
  full.names = TRUE
)) %>%
  as.list() %>% 
  purrr::set_names(~ basename(.) %>% 
                     stringr::str_extract("MiSeq\\d*")) %>% 
  purrr::map(readRDS) %>%
  purrr::map(~ .x %>%
    data.frame() %>%
    tibble::rownames_to_column(., var = "Sample_ID") %>%
    dplyr::select(Sample_ID)) 


## ASVs ----
ASV_tbl_HEMI <- sort(list.files(
  path = path_to_anlys,
  pattern = "_seqtab",
  full.names = TRUE
)) %>%
  purrr::map(readRDS) %>%
  purrr::map(data.frame) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(across(
    .cols = everything(),
    .fns = ~ tidyr::replace_na(.x, 0)
  )) %>%
  as.matrix()


## taxonomy ----
taxa_tbl_HEMI <- sort(list.files(
  path = path_to_anlys,
  pattern = "_taxonomy",
  full.names = TRUE
)) %>%
  purrr::map(readRDS) %>%
  purrr::map(data.frame) %>%
  purrr::map(~ tibble::rownames_to_column(
    .x,
    var = "seqs"
  )) %>%
  dplyr::bind_rows() %>%
  dplyr::group_by(seqs, species) %>%
  dplyr::distinct(seqs, .keep_all = TRUE) %>%
  dplyr::arrange(seqs) %>%
  dplyr::mutate(
    species = paste(species, " ", sep = ""),
    species = ifelse(str_ends(species, ".1"),
      "uncultured bacteria",
      species
    ),
    species = str_extract(species, "^(?:[^ ]+ ){2}"),
    species = str_replace_all(species, " ", "."),
    species = str_replace(species, "\\.", " "),
    species = str_remove_all(species, "\\."),
    across(
      .cols = everything(),
      .fns = ~ stringr::str_replace(
        .x,
        "^Eubacteriales Family",
        "unclassified Eubacteriales"
      )
    ),
    across(
      .cols = everything(),
      .fns = ~ stringr::str_remove(
        .x,
        " Family"
      )
    )
  ) %>%
  tibble::column_to_rownames(var = "seqs")

if (nrow(taxa_tbl_HEMI) == ncol(ASV_tbl_HEMI)) {
  print("Good to go")
} else {
  stop("ASVs and taxa table do not match")
}

cln <- c("I100",
         "I101",
         "I102",
         "I103",
         "I104",
         "I106",
         "I107",
         "I108",
         "I109",
         "I110",
         "I111",
         "I113",
         "I137",
         "I138",
         "I139",
         "I140")

# metatab
meta_tbl <- data.frame(ASV_tbl_HEMI) %>%
  dplyr::mutate(
    Richness = vegan::specnumber(ASV_tbl_HEMI),
    Diversity = vegan::diversity(ASV_tbl_HEMI, index = "simpson"),
    ReadCounts = rowSums(ASV_tbl_HEMI)
  ) %>%
  dplyr::select(ReadCounts, Richness, Diversity) %>%
  tibble::rownames_to_column(var = "Sample_IDs") %>%
  dplyr::mutate(
    ID = str_extract(
      Sample_IDs,
      pattern = "(?<=x)(.+)|(?<=_)(.+)(?=_)"
    ),
    sub_grp = ifelse(!str_detect(ID, "Blank|Empty|Control|PCR"),
      str_extract(ID, pattern = "\\D$"),
      "Ctrls"
    ),
    ID = ifelse(!str_detect(ID, "Blank|Empty|Control|PCR"),
      str_remove(ID, "\\D$"),
      ID
    ),
    Sample_type = case_when(
      (str_ends(sub_grp, "C")) ~ "C_DHS",
      (str_ends(sub_grp, "H")) ~ "H_LHL",
      (str_ends(sub_grp, "J")) ~ "J_FrET",
      (str_ends(sub_grp, "L")) ~ "L_FoET",
      (str_ends(sub_grp, "N|O")) ~ "N_PU8",
      (str_ends(sub_grp, "A")) ~ "A?",
      (str_ends(sub_grp, "E")) ~ "E?",
      TRUE ~ "Blnk_ctrl"
    ),
    sub_grp = ifelse(str_detect(ID, "Blank|Empty|Control|PCR"),
      "Ctrls",
      sub_grp
    ),
    run = ifelse(Sample_IDs %in% ASV_runs$MiSeq42$Sample_ID,
      "MiSeq42",
      "MiSeq47"
    ), 
    clinic = ifelse(ID %in% cln, 
                 "OYS", 
                 "HUS")
  ) %>%
  tibble::column_to_rownames(var = "Sample_IDs") %>%
  dplyr::arrange(nchar(ID), ID) %>%
  dplyr::filter(!is.na(sub_grp))


## metadata ----
metadata <- readxl::read_xlsx(
  here::here(
    "Data",
    "HEMI_meta_n76.xlsx"
  )
) %>%
  dplyr::select(where(
    ~ !all(is.na(.x))
  )) %>% 
  dplyr::mutate(
    age_strat = case_when((age < 30) ~ 1,
                          (age >= 30 & age < 35) ~ 2,
                          (age >= 35) ~ 3),
    bmi = round(bmi, digits = 1),
    bmi_strat = case_when((bmi < 24.9) ~ 1, 
                          (bmi >= 24.9 & bmi < 29.9) ~ 2,
                          (bmi >= 29.9) ~ 3), 
    g_strat = ifelse(g == 0, 
                     0, 
                     1)) %>% 
  dplyr::select(ID, sort(names(.)))

meta_tbl_HEMI <- meta_tbl %>%
  tibble::rownames_to_column(var = "Sample_ID") %>%
  dplyr::left_join(metadata, by = "ID") %>%
  dplyr::select(where(~ !all(is.na(.x)))) %>%
  dplyr::mutate(
    across(
      .cols = where(is.character),
      .fns = ~ tidyr::replace_na(
        data = .x,
        replace = "unknown"
      )
    )
  ) %>% 
  tibble::column_to_rownames(var = "Sample_ID") 


view(summarytools::dfSummary(meta_tbl_HEMI))
view(summarytools::descr(meta_tbl_HEMI))

# Phyloseq objects ----
## compile ----
physeq_HEMI <- ASV_tbl_HEMI %>%
  phyloseq::phyloseq(
    phyloseq::otu_table(ASV_tbl_HEMI, taxa_are_rows = F),
    phyloseq::tax_table(as.matrix(taxa_tbl_HEMI)),
    phyloseq::sample_data(meta_tbl_HEMI)
  ) %>%
  phyloseq::tax_glom(taxrank = "species") 


phy_test <- check_reads(physeq_HEMI)

phy_test <- try({
  physeq_HEMI %>%
    phyloseq::`taxa_names<-`(., .@tax_table@.Data[, 7])
})



## make phylogenetic trees ----
if (remake_trees) {
  trees <- physeq_HEMI %>%
    mk_trees()

  saveRDS(
    trees,
    file = here::here(
      path_to_anlys,
      "3_HEMI_phy_trees.rds"
    )
  )
} else {
  trees <- readRDS(
    file = here::here(
      path_to_anlys,
      "3_HEMI_phy_trees.rds"
    )
  )
}


## merge with phylogenetic objects ----
physeq_wtrees <- physeq_HEMI %>%
  phyloseq::merge_phyloseq(
    trees$tree
  )

phy_test <- check_reads(physeq_wtrees)


## Save phyloseq objects ----
saveRDS(
  physeq_wtrees,
  file = here::here(
    path_to_anlys,
    "3_HEMI_physeq.rds"
  )
)
