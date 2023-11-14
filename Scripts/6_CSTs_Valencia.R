# 9. HEMI - valencia CSTs

# Load packages ----
pkgs_to_load <- list.files(
  here::here(
    "Scripts",
    "to_load"
  ),
  pattern = "^1_",
  full.names = T
)

pkgs_to_load <- purrr::map(
  pkgs_to_load,
  source
)

devtools::install_github("EdwinTh/dutchmasters")


# Define paths ----
pth_anlys <- here::here(
  "Data",
  "analysis"
)

pth_figrs <- here::here(
  "Results",
  "HEMI_plots"
)

if (!dir.exists(here::here(pth_figrs))) {
  dir.create(here::here(pth_figrs), recursive = T)
}


# Load data ----
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
  # remove taxa that are not present in any of the 'J|N' samples
  phyloseq::filter_taxa(
    function(x) sum(x) > 0,
    prune = TRUE
  )


# Raw data
Valencia_input <- physeq_HEMI@otu_table %>% 
  t() %>%
  data.frame() %>%
  dplyr::mutate(sp = physeq_HEMI@tax_table@.Data[,7]) %>%
  dplyr::mutate(g = physeq_HEMI@tax_table@.Data[,6]) %>%
  dplyr::mutate(across(where(is.factor), 
                       as.character)) %>%
  dplyr::select(c(sp, g), everything()) %>%
  dplyr::mutate(g = ifelse(
    stringr::str_detect(
      sp, 
      "uncultured"), 
    "uncultured", 
    g)
  ) %>%
  dplyr::mutate(g = ifelse(
    is.na(g), 
    sp, 
    g)) %>%
  dplyr::mutate(sp = stringr::str_replace(
    sp, 
    " ", 
    "_")) %>%
  dplyr::mutate(sp = ifelse(
    stringr::str_starts(
      g, 
      "Lactobacillus|Gardnerella|Prevotella|Fannyhessea|Sneathia"), 
    sp, 
    paste0("g_", g))) %>%  
  dplyr::select(-g) %>%  
  dplyr::group_by(sp) %>%
  dplyr::summarise(across(everything(), 
                          ~ sum(.))) %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "sp") %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "sampleID") %>%
  dplyr::mutate(read_count = rowSums(
    dplyr::select(., -sampleID))) %>%
  dplyr::select(sampleID, 
                read_count, 
                everything())

write_csv(Valencia_input,
          file = here::here("Data",
                            "HEMI_Valencia_input.csv"))

system2(command = "python3", 
        args = paste("~/Documents/Valencia/Valencia.py", 
                     "-r ~/Documents/Valencia/CST_centroids_012920.csv", 
                     "-i Data/HEMI_Valencia_input.csv", 
                     "-o Data/HEMI_Val_CSTs", 
                     "-p Data/HEMI_Val_plot"), 
        wait = TRUE)


Val_tab <- read.csv(
  here::here("Data",
             "HEMI_Val_CSTs.csv")) %>%
  dplyr::select(sampleID, CST, subCST) %>%
  dplyr::filter(!str_detect(sampleID, "Blank")) %>%
  dplyr::arrange(CST) 


#write_excel_csv2(path = "Valencia_CSTs.csv")
writexl::write_xlsx(Valencia, 
                    path = here::here(pth_figrs,
                                      "HEMI_Valencia_CSTs.xlsx"))

Valencia <- Val_tab %>%
  dplyr::mutate(
    ID = str_extract(
      sampleID,
      pattern = "(?<=x)(.+)|(?<=_)(.+)(?=_)"
    ),
    sub_grp = ifelse(!str_detect(ID, "Blank|Empty|Control|PCR"),
                     str_extract(ID, pattern = "\\D$"),
                     "Ctrls"
    )) %>% 
  dplyr::select(-sampleID) %>% 
  dplyr::group_by(CST) %>% 
  dplyr::arrange(subCST, sub_grp) %>% 
  dplyr::ungroup() %>% 
  tibble::column_to_rownames(var = "ID")


physeq_HEMI_rel <- physeq_HEMI %>% 
  phyloseq::transform_sample_counts(
    function(OTU) OTU/sum(OTU)) %>%
  phyloseq::transform_sample_counts(
    function(OTU) ifelse(is.na(OTU), 0, OTU)
  ) %>%
  phyloseq::`taxa_names<-`(., .@tax_table@.Data[, 7])

ASV_table <- physeq_HEMI_rel %>%
  phyloseq::psmelt() %>% 
  dplyr::rename("taxa" = "OTU") %>% 
  tidyr::unite("ID", ID:sub_grp, sep = "") %>% 
  dplyr::mutate(
    taxa = if_else(Abundance < 0.02, "Others", taxa)) %>%
  dplyr::group_by(taxa, ID) %>% 
  dplyr::summarise(Abundance = sum(Abundance)) %>% 
  dplyr::group_by(taxa) %>% 
  dplyr::add_tally(Abundance > 0) %>% 
  dplyr::filter(n > 2) %>% 
  tidyr::pivot_wider(id_cols = ID, 
                     names_from = taxa, 
                     values_from = Abundance) %>%
  tibble::column_to_rownames(var = "ID") %>%
  dplyr::mutate(
    across(
      .cols = everything(), 
      .fns = ~ tidyr::replace_na(.x, 0))) %>% 
  t() %>% 
  data.frame() %>% 
  dplyr::select(rownames(Valencia)) %>% 
  as.matrix()

# Specify color
## CSTs
CST_colours <- grDevices::hcl.colors(
  length(unique(Valencia$CST)), 
  palette = "TealRose")
names(CST_colours) <- unique(Valencia$CST)

## subCSTs
subCST_colours <- dutchmasters::dutchmasters_pal(
  palette = "milkmaid")(length(unique(Valencia$subCST)))
names(subCST_colours) <- unique(Valencia$subCST)

##sub-groups
sub_grp_clrs <- grDevices::hcl.colors(
  2, 
  palette = "LightGrays"
)
names(sub_grp_clrs) <- unique(Valencia$sub_grp)

## annotation col
annot_colours <- list(
  CST = CST_colours, 
  subCST = subCST_colours, 
  sub_grp = sub_grp_clrs)

hm <- pheatmap::pheatmap(ASV_table, 
                        color = rev(grDevices::hcl.colors(
                          100, 
                          alpha = 1, 
                          rev = TRUE, 
                          palette = "Rocket"
                        )),
                        annotation_col = Valencia,
                        annotation_colors = annot_colours,
                        clustering_method = "ward.D2", 
                        cutree_cols = 5,
                        fontsize = 14, 
                        treeheight_row = 200,
                        treeheight_col = 200,
                        angle_col = 90,
                        border_color = "grey80")
dev.off()

# modify x-axis labels
hm$gtable$grobs[[4]]$gp = grid::gpar(
  col = "grey40", 
  fontsize = 10
  )

# modify y-axis labels
hm$gtable$grobs[[5]]$gp = grid::gpar(
  fontface = "italic", 
  col = "grey40"
  )

# modify annotation col labels
hm$gtable$grobs[[7]]$gp = grid::gpar(
  col = "grey40", 
  fontsize = 10, 
  fontface = "plain"
  )

ggplot2::ggsave(
  hm, 
  path = here::here(
    pth_figrs),
  filename = "HEMI_heatmap_CSTs.pdf",
  device = "pdf", 
  dpi = 1200, 
  width = 25, 
  height = 10, 
  units = "in")

ggplot2::ggsave(
  hm, 
  path = here::here(
    pth_figrs),
  filename = "HEMI_heatmap_CSTs.eps",
  device = cairo_ps, 
  dpi = 1200, 
  width = 25, 
  height = 10, 
  units = "in")

ids <- data.frame(physeq_HEMI@sam_data) %>% 
  dplyr::group_by(ID) %>% 
  dplyr::mutate(pairs = n_distinct(sub_grp)) %>% 
  dplyr::filter(pairs == 2) %>% 
  dplyr::filter(prgcy == 1)

tab <- Valencia %>% 
  tibble::rownames_to_column(var = "ID") %>% 
  dplyr::mutate(
    sub_grp = stringr::str_extract(
      ID, "\\D$"
    ), 
    ID = stringr::str_remove(
      ID, "\\D$"
    )
  ) %>% 
  dplyr::filter(ID %in% ids$ID) %>% 
  dplyr::arrange(ID, sub_grp)

writexl::write_xlsx(tab, 
                    path = here::here(pth_figrs, 
                                      "HEMI_CSTs_JvsN.xlsx"))
