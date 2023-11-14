# 5. HEMI - PERMANOVA

# Session parameters and packages ----
pkgs_to_load <- list.files(
  here::here(
    "Scripts",
    "to_load"
  ),
  pattern = "^1_|^2_|^7_",
  full.names = T
)

pkgs_to_load <- purrr::map(
  pkgs_to_load,
  source
)

## define major directory paths ----
pth_anlys <- here::here(
  "Data",
  "Analysis"
)

pth_rslts <- here::here(
  "Results",
  "stats"
)

if (!dir.exists(here::here(pth_rslts))) {
  dir.create(here::here(pth_rslts), recursive = T)
}

pth_figrs <- here::here(
  "Results", 
  "HEMI_plots"
)

if(!dir.exists(here::here(pth_figrs))) {
  dir.create(here::here(pth_figrs), recursive = T)
}

# Load data ----
# (both J and N samples, read-count cutoff 500)
physeq_HEMI <- readr::read_rds(
  here::here(
    pth_anlys,
    "3_HEMI_physeq.rds"
  )
) %>%
  phyloseq::subset_samples(
    sub_grp %in% c("J", "N")
  ) %>%
  phyloseq::subset_samples(
    ReadCounts >= 500
  ) %>% 
  phyloseq::filter_taxa(
    function(x) sum(x) > 0,
    prune = TRUE
  )


# Convert raw counts to relative abundance
physeq_HEMI_rel <- phyloseq::transform_sample_counts(
  physeq_HEMI,
  function(OTU) OTU / sum(OTU)
) %>%
  phyloseq::transform_sample_counts(
    function(OTU) ifelse(is.na(OTU), 0, OTU)
  )


##### Summary of metadata variables ##-----
 view(summarytools::dfSummary(
   data.frame(physeq_HEMI_rel@sam_data),
   plain.ascii = FALSE,
   style = "grid"))


## Define plot theme
plot_theme(text_size = 8)

dt <- phyloseq::psmelt(physeq_HEMI_rel) %>% 
  dplyr::rowwise() %>% 
  dplyr::filter(
    any(sub_grp == "N" & prgcy == 1) | 
      any(sub_grp == "J" & !is.na(prgcy)
    )
  )

tally <- dt %>% 
  dplyr::group_by(ID, sub_grp) %>% 
  dplyr::distinct(ID, .keep_all = TRUE) %>% 
  dplyr::group_by(sub_grp) %>% 
  dplyr::tally()


data_step <- dt %>% 
  dplyr::select(
    OTU, 
    Sample, 
    Abundance
  ) %>% 
  tidyr::pivot_wider(
    id_cols = "Sample", 
    names_from = OTU, 
    values_from = Abundance
  ) %>% 
  dplyr::arrange(Sample) %>% 
  tibble::column_to_rownames("Sample")


# subset metadata
meta_step <- dt %>% 
  dplyr::group_by(ID, sub_grp) %>% 
  dplyr::distinct(ID, .keep_all = TRUE) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-c(OTU, superkingdom:species, ID)) %>% 
  dplyr::arrange(Sample) %>% 
  tibble::column_to_rownames("Sample")

if (
  all(rownames(data_step) == rownames(meta_step))
) {
  print("Sample order is correct")
} else {
  stop("Sample order incorrect")
}

sub <- meta_step %>% 
  dplyr::mutate(across(
    .cols = everything(), 
    .fns = ~ factor(.x, levels = unique(.x))
  )) %>% 
  dplyr::summarise(across(
    .cols = everything(), 
    .fns = ~ nlevels(.x)
  )) %>% 
  t() %>% 
  data.frame() %>% 
  purrr::set_names("levels") %>% 
  tibble::rownames_to_column(var = "Var") %>% 
  dplyr::filter(levels > 10)


meta_list <- meta_step %>%
  dplyr::mutate(ReadCounts = as.character(ReadCounts)) %>%
  dplyr::select(-c(sub_grp, Sample_type)) %>%
  as.list() %>%
  purrr::imap(.f = ~ .x %>% 
                data.frame() %>%
      purrr::set_names(.y) %>%
      dplyr::mutate(Sample = rownames(meta_step)) %>%
      dplyr::filter(across(.cols = everything(), 
                           .fns = ~ !is.na(.x))
      )
  ) %>%
  purrr::imap(
    .f = ~ fctrs(
      df = .x, 
      var_nm = .y, 
      excl = sub
    )
  ) %>% 
  purrr::map_at(
    .at = "ReadCounts",
    .f = function(x) {
      x <- x %>%
        dplyr::mutate(ReadCounts = as.numeric(ReadCounts))
    }
  )


adonis_res <- meta_list %>%
  purrr::imap(.f = ~ PRMNV_run(
    df = data_step,
    int_vat = .x,
    var_nm = .y, 
    permt = 99999, 
    mthd = "bray"
  )) %>%
  dplyr::bind_rows()

# clean output
adonis_scrub <- adonis_clean(adonis_res)

# write excel file
writexl::write_xlsx(adonis_scrub,
  path = here::here(
    pth_rslts,
    "HEMI_PERMANOVA.xlsx"
  )
)

# save PERMANOVA results
saveRDS(adonis_scrub,
  file = here::here(
    "Data",
    "MiSeq42_HEMI_PERMANOVA.rds"
  )
)

clev_plotdata <- adonis_scrub %>%
  dplyr::arrange(Pr_F_stat) %>%
  dplyr::select(Var, Pr_F_stat) %>%
  dplyr::mutate(Var = factor(Var, levels = unique(Var)))

plt_title <- meta_step %>% 
  dplyr::group_by(sub_grp) %>% 
  dplyr::tally() %>% 
  tidyr::unite(col = "comb", sep = ": ") %>% 
  as.matrix() %>% 
  paste(., collapse = ", ")

p <- ggplot2::ggplot(
  clev_plotdata,
  aes(
    x = Pr_F_stat,
    y = Var
  )
) +
  ggplot2::geom_segment(aes(
    x = 0,
    xend = Pr_F_stat,
    y = Var,
    yend = Var
  ),
  colour = "grey80"
  ) +
  ggplot2::geom_jitter(aes(colour = "red"),
    size = 4,
    height = 0.1,
    width = 0.01,
    show.legend = F
  ) +
  ggplot2::geom_vline(
    xintercept = 0.05,
    colour = viridis::rocket(1, begin = 0.4),
    linetype = "dashed",
    size = 0.2
  ) +
  ggplot2::scale_x_log10(limits = c(
    min(na.omit(clev_plotdata$Pr_F_stat)) * 0.5,
    max(clev_plotdata$Pr_F_stat) * 1.5
  )) +
  viridis::scale_colour_viridis(
    option = "rocket",
    begin = 0.2,
    end = 0.8,
    discrete = T,
    alpha = 0.7,
    name = ""
  ) +
  ggplot2::ggtitle(paste0(
    "HEMI - ", 
    plt_title
    )) +
  ggplot2::theme(plot.title = element_text(
    hjust = 0.5,
    colour = "grey40"
  )) +
  ggplot2::xlab("p-values") +
  ggplot2::ylab("Variables")

ggsave(p,
  path = pth_figrs,
  filename = "HEMI_pvals_PERMANOVA.pdf",
  device = "pdf",
  height = 30,
  width = 20,
  units = "cm",
  dpi = 320,
  limitsize = F
)
