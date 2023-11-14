# 4. HEMI - Manuscript plots

# Load packages ----
pkgs_to_load <- list.files(
  here::here(
    "Scripts",
    "to_load"
  ),
  pattern = "^1_|^2_|^5_|^8_|^circ",
  full.names = T
)

pkgs_to_load <- purrr::map(
  pkgs_to_load,
  source
)


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
phy_HEMI <- readr::read_rds(
  here::here(
    pth_anlys,
    "3_HEMI_physeq.rds"
  )
) %>%
  phyloseq::subset_samples(
    !is.na(prgcy)
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


phy_HEMI_rel <- phyloseq::transform_sample_counts(
  phy_HEMI, function(OTU) OTU / sum(OTU)
) %>%
  phyloseq::transform_sample_counts(
    function(OTU) ifelse(is.na(OTU), 0, OTU)
  ) %>%
  phyloseq::`taxa_names<-`(., .@tax_table@.Data[, 7])


HEMI_scrs <- list.files(
  here::here(pth_anlys),
  pattern = "^1_.*scr|^2_.*scr",
  full.names = TRUE
) %>%
  purrr::map(
    ~ readRDS(
      .x
    )
  ) %>%
  dplyr::bind_rows() %>%
  dplyr::arrange(species) %>%
  dplyr::group_by(.) %>%
  dplyr::distinct(species, .keep_all = TRUE) %>%
  dplyr::group_by(species) %>%
  dplyr::rename(taxa = species) %>%
  dplyr::mutate(
    taxa = stringr::str_replace(
      taxa, "\\.", ""
    )
  ) %>%
  dplyr::ungroup()


# data prep
plt_dt <- phy_HEMI_rel %>%
  phyloseq::psmelt() %>%
  dplyr::arrange(Sample) %>%
  dplyr::rename("taxa" = "OTU") %>%
  dplyr::filter(!is.na(prgcy)) %>%
  dplyr::rowwise() %>% 
  dplyr::filter(
    any(sub_grp == "N" & prgcy == 1) | sub_grp == "J"
  )

mk_tree <- FALSE


# Illustrations - J ----

# 46 J (0), 29 J (1), 22 N (1)
tally <- plt_dt %>% 
  dplyr::group_by(sub_grp) %>% 
  dplyr::distinct(Sample, .keep_all = TRUE) %>% 
  dplyr::select(Sample, sub_grp, prgcy) %>% 
  
  dplyr::arrange(sub_grp, prgcy) %>% 
  dplyr::group_by(sub_grp, prgcy) %>% 
  dplyr::tally()

plt_dt <- plt_dt %>% 
  dplyr::mutate(
    prgcy = ifelse(
      prgcy == 1,
      "pregnancy",
      "non-pregnancy"
    ),
    birth = ifelse(
      birth == 1, 
      "live birth", 
      "no live birth"))


## Figure 1 A----
### polar plot for clinical pregnancy ----
#### clean ----
plt_dt_plr_prep <- plt_dt %>%
  dplyr::filter(sub_grp == "J") %>% 
  dplyr::left_join(
    HEMI_scrs,
    by = "taxa"
  ) %>%
  dplyr::filter(
    !is.na(qual)
  ) %>%
  dplyr::mutate(
    Abundance = tidyr::replace_na(
      Abundance, 0
    )
  ) %>%
  dplyr::group_by(taxa) %>%
  dplyr::mutate(
    taxa = stringr::str_replace_all(
      taxa, "^c .|\\.\\.$", ""
    ),
    freq = sum(Abundance > 0.50),
  ) %>%
  dplyr::ungroup()


# factor levels for taxa
plt_dt_plr_tx <- plt_dt_plr_prep %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(
    taxa = if_else(
      Abundance < 0.02, "Others", taxa
    ),
    qual = ifelse(
      taxa == "Others", NA, qual
    ), 
    tax_ord = ifelse(stringr::str_detect(
      taxa, "Lactobacillus"), 1, 0)
  ) %>%
  dplyr::group_by(ID) %>%
  dplyr::arrange(
    desc(prgcy),
    desc(tax_ord),
    desc(freq),
    taxa
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    taxa = factor(
      taxa,
      levels = unique(taxa)
    ),
    taxa = forcats::fct_relevel(taxa, "Others", after = Inf)
  )


abd_test <- plt_dt_plr_tx %>% 
  dplyr::group_by(Sample) %>% 
  dplyr::summarise(Abundance = sum(Abundance)) %>%
  dplyr::mutate(Abundance = signif(Abundance, digits = 5))

if (all(abd_test$Abundance == 1)) {
  print("Good to go")
} else {
  stop("Check again")
}


# create variable for sample levels
plt_dt_plr_lev <- plt_dt_plr_tx %>%
  dplyr::mutate(
    high = ifelse(Abundance >= 0.50,
                  sum(Abundance >= 0.50) / as.integer(taxa),
                  0
    ),
    mid = ifelse(Abundance >= 0.30 & Abundance < 0.50,
                 sum(Abundance >= 0.30) / as.integer(taxa),
                 0
    ),
    low = ifelse(Abundance >= 0.10 & Abundance < 0.30,
                 sum(Abundance >= 0.10) / as.integer(taxa),
                 0
    )
  )

# clean
plt_dt_plr_cln <- plt_dt_plr_lev %>%
  dplyr::arrange(
    desc(high),
    desc(mid),
    desc(low),
    taxa,
    desc(Abundance)
  ) %>%
  dplyr::mutate(
    ID = factor(ID, levels = unique(ID)),
    Sample_type = factor(Sample_type, levels = unique(Sample_type)),
    prgcy = factor(prgcy, levels = unique(prgcy))
  ) %>%
  as.data.frame()


#### modify  ----
labls <- plt_dt_plr_cln %>%
  dplyr::distinct(taxa, .keep_all = T) %>%
  dplyr::arrange(taxa) %>%
  dplyr::select(taxa, qual) %>%
  dplyr::mutate(
    taxa = as.character(taxa),
    taxa = stringr::str_replace_all(
      taxa, " ", "\\."
    )
  ) %>%
  purrr::transpose() %>%
  purrr::modify_depth(
    .depth = 1,
    ~ purrr::map2(
      .x = .x[[1]],
      .y = .x[[2]],
      ~ paste0("italic(", .x, "^", .y, ")")
    )
  ) %>%
  unlist() %>%
  purrr::map(~ str2expression(.x))


# parameters
clrs <- pick_colrs(plt_dt_plr_cln$taxa)
pal2 <- PNWColors::pnw_palette(name = "Sunset2", n = 5) %>%
  .[c(1, 5)]


#### plot ----
polar <- circ_bar(
  df = plt_dt_plr_cln,
  x = "prgcy", 
  y = "Abundance", 
  id = "ID",
  taxa = "taxa",
  labels = labls,
  clr_fill = clrs,
  breaks = 5, 
  lgd_rows = nlevels(plt_dt_plr_cln$taxa), 
  lgd_pos = "left"
)



## Figure 1 B ----
### stacked and summarised plot ----
sample_tally <- plt_dt %>%
  dplyr::filter(sub_grp == "J") %>% 
  dplyr::distinct(ID, .keep_all = T) %>%
  dplyr::group_by(prgcy) %>%
  dplyr::count()


plt_dt_stck <- plt_dt_plr_prep %>%
  dplyr::group_by(taxa, prgcy) %>%
  dplyr::summarise(Abundance = mean(Abundance)) %>%
  dplyr::ungroup() %>%
  #dplyr::mutate(
  #  taxa = ifelse(Abundance < 0.01, 
  #                "Others",
  #                taxa)
  #) %>% 
  dplyr::mutate(taxa = str_remove(taxa, "X.")) %>%
  dplyr::mutate(taxa = str_replace(taxa, "\\.", " "), 
                tax_ord = ifelse(stringr::str_detect(
      taxa, "Lactobacillus"), 1, 0)
  ) %>%
  dplyr::arrange(
    #desc(tax_ord),
    desc(Abundance)) %>% 
  dplyr::mutate(taxa = factor(
    #taxa, levels = unique(taxa))
    taxa, levels(plt_dt_plr_tx$taxa)
  )) %>% 
  dplyr::arrange(taxa)


# define species colours
mycols <- clrs %>% 
  #pick_colrs(plt_dt_stck$taxa, n = 4) %>%
  data.frame() %>%
  purrr::set_names("colours") %>%
  tibble::rownames_to_column(var = "taxa") %>%
  dplyr::mutate(lightness = shades::lightness(colours)) %>%
  dplyr::mutate(shades = if_else(lightness > 65, "light", "dark"))


if (all(mycols$colours == clrs)) {
  print("identical colour palette - good to go")
} else {
  stop("colour palettes do not match")
}


plt_dt_stck_wclrs <- plt_dt_stck %>%
  dplyr::inner_join(mycols) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    taxa = factor(taxa, levels = unique(taxa))) %>%
  dplyr::arrange(taxa) %>%
  dplyr::mutate(shades = factor(shades, levels = c("dark", "light"))) %>%
  data.frame()


plot_theme(text_size = 10)

stack_bar_summ <- ggplot2::ggplot(
  plt_dt_stck_wclrs, aes(
    x = prgcy,
    y = Abundance
  )
) +
  
  # add stacked bar
  ggplot2::geom_bar(aes(fill = taxa),
                    stat = "identity",
                    position = position_stack(reverse = T), 
                    show.legend = FALSE
  ) +
  
  # add manual colour to bars
  ggplot2::scale_fill_manual(
    values = unique(plt_dt_stck_wclrs$colours), 
    name = "Species"
  ) +
  ggplot2::geom_text(
    aes(
      label = if_else(Abundance > 0.01,
                      as.character(round(Abundance * 100, digits = 2)), 
                      ""
      ),
      size = Abundance,
      colour = shades,
      group = taxa,
      angle = if_else(Abundance > 0.10, 0, 270)
    ),
    position = position_stack(0.5, reverse = T),
    show.legend = F
  ) +
  
  # add scale colour to text
  ggplot2::scale_colour_manual(
    values = c("grey85", "grey35"), name = "Species") +
  
  # add sample tally
  ggplot2::geom_text(
    data = sample_tally,
    aes(
      x = prgcy,
      y = 1.05,
      label = paste("n = ", n, sep = "")
    ),
    colour = "grey40"
  ) +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 270, colour = "grey40")
  ) +
  ggplot2::xlab("") +
  ggplot2::ylab("Relative Abundance (%)") +
  ggplot2::coord_flip() +
  ggplot2::guides(fill = guide_legend(nrow = 5), direction = "vertical")




### grouptest ----
# (Only J samples, read-count cutoff 500)
phy_J <- phy_HEMI %>%
  phyloseq::subset_samples(
    sub_grp == "J"
  ) %>%
  phyloseq::`sample_data<-`(
    .@sam_data %>% 
      data.frame() %>% 
      dplyr::rename("ReadCount" = ReadCounts)) %>% 
  
  # remove taxa that are not present in any of the 'J' samples
  phyloseq::filter_taxa(
    function(x) sum(x) > 0,
    prune = TRUE
  ) %>%
  phyloseq::`taxa_names<-`(., .@tax_table@.Data[, 7])


phy_J_rel <- phyloseq::transform_sample_counts(
  phy_J, 
  function(OTU) OTU / sum(OTU)
) %>%
  phyloseq::transform_sample_counts(
    function(OTU) ifelse(is.na(OTU), 0, OTU)
  ) 


plot_data <- phy_J_rel %>% 
  phyloseq::psmelt() %>% 
  dplyr::rename("Taxa" = OTU) %>% 
  dplyr::mutate(Taxa = stringr::str_replace(Taxa, " ", "\\."))

plot_theme()


# J = 75, 0 = 46, 1 = 29
title <- plot_data %>% 
  dplyr::distinct(ID, .keep_all = TRUE) %>% 
  dplyr::group_by(sub_grp, prgcy) %>% 
  dplyr::tally() %>%
  dplyr::mutate(
    prgcy = ifelse(
      prgcy == 1,
      "pregnancy",
      "non-pregnancy")
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-sub_grp) %>% 
  tidyr::unite(col = "comb", sep = ": ") %>% 
  as.matrix() %>% 
  paste(., collapse = ", ")

# set conditions
ref <- sym("prgcy")

cnfr <- list(NULL, 
             c("age_strat", "g_strat", "p"))

names <- list("None", "PGA")

step <- 1

gt_dt <- data.frame()

# Run tests and create plots
while (step <= length(cnfr)) {
  
  gt_res <- try({GroupTest_mod(
    species.table = phy_J,
    meta = phy_J,
    group = paste(ref), 
    group_name = paste(ref), 
    xl_titl = names[[step]],
    confounders = paste(cnfr[[step]]),
    dir_for_res = here::here("Results", "mare"), 
    min.prevalence = 0.1, 
    min.abundance = 0.1,
    outlier.cutoff = 3, 
    p.cutoff = 0.05, 
    keep.result = T, 
    nonzero = F,
    pdf = F, 
    show_quartz = F)
  }) 
  
  if (class(gt_res) == "try-error") {
    gt_res <- c()
  }
  
  gt_res_clean <- gt_res %>%
    dplyr::select(taxon, ends_with("_FDR")) %>%
    tidyr::pivot_longer(cols = -taxon, 
                        names_to = "comp", 
                        values_to = "pvals") %>%
    dplyr::group_by(taxon) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(variable = str_extract(
      comp, 
      pattern = paste(
        levels(
          as.factor(
            data.frame(phy_J@sam_data)[,paste(ref)])), 
        collapse = "|"))) %>%
    dplyr::rename(!!ref := variable) %>%
    dplyr::select(-comp) %>%
    base::as.data.frame() %>%
    {
      p <- .
      if (!is.na(as.numeric(p[,paste(ref)][1])))
        p <- p %>%
          dplyr::mutate(!!ref := as.numeric(!!ref))
      else 
        .
    } %>%
    dplyr::rename("Taxa" = taxon) %>% 
    dplyr::mutate(test = names[[step]])
  
  gt_dt <- gt_dt %>% 
    dplyr::bind_rows(gt_res_clean)
  
  step <- step + 1 
}



## Clean and rearrange data
p2 <- plot_data %>%
  dplyr::select(Taxa, Abundance, !!ref) %>%
  dplyr::filter(Taxa %in% gt_dt$Taxa) %>%
  dplyr::full_join(y = gt_dt) %>%
  dplyr::mutate(
    Taxa = str_replace(Taxa, "\\.", " "), 
    sig = case_when((pvals > 0 & pvals <= 0.001) ~ "***",
                    (pvals > 0.001 & pvals <= 0.01) ~ "**",
                    (pvals > 0.01 & pvals <= 0.05) ~ "*",
                    (pvals > 0.05) ~ "")
  ) %>%
  dplyr::mutate(prgcy = ifelse(
    prgcy == 1, 
    "pregnancy",
    "non-pregnancy"
  )) %>% 
  dplyr::mutate(!!ref := as.factor(!!ref))


## significance labels
p2_lab <- p2 %>%
  dplyr::filter(!is.na(sig)) %>% 
  dplyr::group_by(Taxa, !!ref, test) %>% 
  dplyr::distinct(Taxa, .keep_all = TRUE) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(
    pvals = signif(pvals, digits = 2),
    pvals = ifelse(pvals < 0.001, "<0.001", pvals), 
    pvals = ifelse(test == "None", 
                   paste0("<sup>1</sup>", pvals, sig), 
                   paste0("<sup>2</sup>", pvals, sig))) %>% 
  dplyr::group_by(Taxa, Abundance, prgcy) %>% 
  dplyr::summarise(pvals = paste0(pvals, collapse = "<br>"))


## assign colours for plot
pal <- c("#946666", "#66948e")


#### viobox ----
gt_viobox <- ggplot2::ggplot(
  data = p2 %>%
    dplyr::filter(Abundance > 0), 
  aes(x = !!ref,
      y = Abundance,
      fill = !!ref)) +
  
  ## Add jitter plot
  ggplot2::geom_jitter(aes(colour = !!ref),
                       show.legend = F, 
                       size = 5,
                       alpha = 0.8) +
  
  ## half violin on the left
  gghalves::geom_half_violin(
    mapping = aes(colour = !!ref),
    trim = F, 
    side = "l",
    draw_quantiles = c(0.25, 0.5, 0.75), 
    scale = "area",
    #colour = "grey70",
    size = .5,
    alpha = 0.7) +
  
  ## half box on the right
  gghalves::geom_half_boxplot(
    mapping = aes(colour = !!ref),
    side = "r", 
    nudge = .05,
    size = .5, 
    #colour = "grey70",
    outlier.shape = NA, 
    alpha = 0.7,
    show.legend = F) +
  
  ## facet by Taxa
  ggplot2::facet_grid(. ~ Taxa, scales = "free_y") +
  ggtext::geom_richtext(
    data = p2_lab, 
    aes(x = !!ref,
        y = 20,
        label = pvals),
    colour = "grey60",
    size = 6,
    alpha = 0.85,
    inherit.aes = F, 
    show.legend = FALSE
  ) +
  
  ## add colours
  ggplot2::scale_fill_manual(values = pal, name = "") +
  ggplot2::scale_colour_manual(values = pal, name = "")  +
  
  ## log10 transform y-axis
  ggplot2::scale_y_log10() +
  
  ## add point for median
  ggplot2::stat_summary(
    fun = "mean", 
    geom = "point", 
    colour = "#c1e0dd", 
    size = 5, 
    show.legend = FALSE) +
  
  ## adjust labels and theme
  ggplot2::xlab("") +
  ggplot2::ylab("Relative abundance (log10)") +
  ggplot2::theme(
    axis.text.x = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(1, "cm"), 
    plot.title = element_text(colour = "grey50")) +
  ggplot2::guides(fill = guide_legend(ncol = 2))


gt_viobox

### compile and save infographic ----
comb_plot <- cowplot::plot_grid(
  polar,
  stack_bar_summ, 
  gt_viobox, 
  nrow = 3, 
  rel_heights = c(1, 0.3, 0.6), 
  labels = "auto", 
  label_colour = "grey30", 
  label_size = 24
) 


ggplot2::ggsave(
  path = here::here(pth_figrs),
  plot = comb_plot,
  filename = paste0(
    "Fig1.pdf"), 
  device = "pdf", 
  dpi = 1200,
  width = 35,
  height = 45,
  units = "cm", 
  limitsize = F)


ggplot2::ggsave(
  path = here::here(pth_figrs),
  plot = comb_plot,
  filename = paste0(
    "Fig1.eps"), 
  device = cairo_ps, 
  dpi = 1200,
  width = 35,
  height = 45,
  units = "cm",
  limitsize = F)




# Illustrations - birth ----
## Figure 2 A----
### polar plot for live birth ----
#### clean ----
plt_dt_brt_prep <- plt_dt %>%
  dplyr::filter(!is.na(birth), 
                sub_grp == "J") %>% 
  dplyr::left_join(
    HEMI_scrs,
    by = "taxa"
  ) %>%
  dplyr::filter(
    !is.na(qual)
  ) %>%
  dplyr::mutate(
    Abundance = tidyr::replace_na(
      Abundance, 0
    )
  ) %>%
  dplyr::group_by(taxa) %>%
  dplyr::mutate(
    taxa = stringr::str_replace_all(
      taxa, "^c .|\\.\\.$", ""
    ),
    freq = sum(Abundance > 0.50),
  ) %>%
  dplyr::ungroup() 


# 50 (0) & 25 (1)
tally <- plt_dt_brt_prep %>% 
  dplyr::group_by(sub_grp, birth) %>%
  dplyr::distinct(Sample, .keep_all = TRUE) %>% 
  dplyr::count()

# factor levels for taxa
plt_dt_brt_tx <- plt_dt_brt_prep %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(
    taxa = if_else(
      Abundance < 0.02, "Others", taxa
    ),
    qual = ifelse(
      taxa == "Others", NA, qual
    ), 
    tax_ord = ifelse(stringr::str_detect(
      taxa, "Lactobacillus"), 1, 0)
  ) %>%
  dplyr::group_by(ID) %>%
  dplyr::arrange(
    desc(birth),
    desc(tax_ord),
    desc(freq),
    taxa
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    taxa = factor(
      taxa,
      levels = unique(taxa)
    ),
    taxa = forcats::fct_relevel(taxa, "Others", after = Inf), 
    taxa = forcats::fct_relevel(taxa, "Lactobacillus crispatus", after = 0L)
  )


t <- plt_dt_brt_tx %>% 
  dplyr::group_by(Sample) %>% 
  dplyr::summarise(Abundance = sum(Abundance)) %>%
  dplyr::mutate(Abundance = signif(Abundance, digits = 5))

if (all(t$Abundance == 1)) {
  print("Good to go")
} else {
  stop("Check again")
}


# create variable for sample levels
plt_dt_brt_lev <- plt_dt_brt_tx %>%
  dplyr::mutate(
    high = ifelse(Abundance >= 0.50,
                  sum(Abundance >= 0.50) / as.integer(taxa),
                  0
    ),
    mid = ifelse(Abundance >= 0.30 & Abundance < 0.50,
                 sum(Abundance >= 0.30) / as.integer(taxa),
                 0
    ),
    low = ifelse(Abundance >= 0.10 & Abundance < 0.30,
                 sum(Abundance >= 0.10) / as.integer(taxa),
                 0
    )
  )

# clean
plt_dt_brt_cln <- plt_dt_brt_lev %>%
  dplyr::arrange(
    desc(high),
    desc(mid),
    desc(low),
    taxa,
    desc(Abundance)
  ) %>%
  dplyr::mutate(
    ID = factor(ID, levels = unique(ID)),
    Sample_type = factor(Sample_type, levels = unique(Sample_type)),
    birth = factor(birth, levels = unique(birth))
  ) %>%
  as.data.frame()


#### modify ----
labls <- plt_dt_brt_cln %>%
  dplyr::distinct(taxa, .keep_all = T) %>%
  dplyr::arrange(taxa) %>%
  dplyr::select(taxa, qual) %>%
  dplyr::mutate(
    taxa = as.character(taxa),
    taxa = stringr::str_replace_all(
      taxa, " ", "\\."
    )
  ) %>%
  purrr::transpose() %>%
  purrr::modify_depth(
    .depth = 1,
    ~ purrr::map2(
      .x = .x[[1]],
      .y = .x[[2]],
      ~ paste0("italic(", .x, "^", .y, ")")
    )
  ) %>%
  unlist() %>%
  purrr::map(~ str2expression(.x))



# parameters
clrs_brt <- pick_colrs(plt_dt_brt_cln$taxa)
pal2 <- PNWColors::pnw_palette(name = "Sunset2", n = 5) %>%
  .[c(1, 5)]


#### plot ----
polar_brt <- circ_bar(
  df = plt_dt_brt_cln,
  x = "birth", 
  y = "Abundance", 
  id = "ID",
  taxa = "taxa",
  labels = labls,
  clr_fill = clrs_brt,
  breaks = 5, 
  lgd_rows = nlevels(plt_dt_brt_cln$taxa), 
  lgd_pos = "left"
)



## Figure 2 B----
### stacked and summarised plot ----
sample_tally <- plt_dt %>%
  dplyr::filter(!is.na(birth), 
                sub_grp == "J") %>% 
  dplyr::distinct(ID, .keep_all = T) %>%
  dplyr::group_by(birth) %>%
  dplyr::count()


plt_dt_brt_stck <- plt_dt_brt_prep %>%
  dplyr::group_by(taxa, birth) %>%
  dplyr::summarise(Abundance = mean(Abundance)) %>%
  dplyr::ungroup() %>%
  #dplyr::mutate(
  #  taxa = ifelse(Abundance < 0.01, 
  #                "Others",
  #                taxa)
  #) %>% 
  dplyr::mutate(taxa = str_remove(taxa, "X.")) %>%
  dplyr::mutate(taxa = str_replace(taxa, "\\.", " ")) %>%
  dplyr::arrange(desc(Abundance)) %>% 
  dplyr::mutate(taxa = factor(
    taxa, levels = levels(plt_dt_brt_tx$taxa))
  ) %>% 
  dplyr::arrange(taxa)


# define species colours
mycols_brt <- clrs_brt %>% 
  data.frame() %>%
  purrr::set_names("colours") %>%
  tibble::rownames_to_column(var = "taxa") %>%
  dplyr::mutate(lightness = shades::lightness(colours)) %>%
  dplyr::mutate(shades = if_else(lightness > 65, "light", "dark"))


if (all(mycols_brt$colours == clrs_brt)) {
  print("identical colour palette - good to go")
} else {
  stop("colour palettes do not match")
}


plt_dt_brt_stck_wclrs <- plt_dt_brt_stck %>%
  dplyr::inner_join(mycols_brt) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    taxa = factor(taxa, levels = levels(plt_dt_brt_tx$taxa))) %>%
  dplyr::arrange(taxa) %>%
  dplyr::mutate(shades = factor(shades, levels = c("dark", "light"))) %>%
  data.frame()


plot_theme(text_size = 10)

stack_bar_brt_summ <- ggplot2::ggplot(
  plt_dt_brt_stck_wclrs, aes(
    x = birth,
    y = Abundance
  )
) +
  
  # add stacked bar
  ggplot2::geom_bar(aes(fill = taxa),
                    stat = "identity",
                    position = position_stack(reverse = T), 
                    show.legend = FALSE
  ) +
  
  # add manual colour to bars
  ggplot2::scale_fill_manual(
    values = unique(plt_dt_brt_stck_wclrs$colours), 
    name = "Species"
  ) +
  ggplot2::geom_text(
    aes(
      label = if_else(Abundance > 0.01,
                      as.character(round(Abundance * 100, digits = 2)), 
                      ""
      ),
      size = Abundance,
      colour = shades,
      group = taxa,
      angle = if_else(Abundance > 0.10, 0, 270)
    ),
    position = position_stack(0.5, reverse = T),
    show.legend = F
  ) +
  
  # add scale colour to text
  ggplot2::scale_colour_manual(
    values = c("grey85", "grey35"), name = "Species") +
  
  # add sample tally
  ggplot2::geom_text(
    data = sample_tally,
    aes(
      x = birth,
      y = 1.05,
      label = paste("n = ", n, sep = "")
    ),
    colour = "grey40"
  ) +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 270, colour = "grey40")
  ) +
  ggplot2::xlab("") +
  ggplot2::ylab("Relative Abundance (%)") +
  ggplot2::coord_flip() +
  ggplot2::guides(fill = guide_legend(nrow = 5), direction = "vertical")




### grouptest ----
# (Only J samples, read-count cutoff 500)
phy_J_birth <- phy_HEMI %>% 
  phyloseq::subset_samples(
    sub_grp == "J"
  ) %>% 
  phyloseq::`sample_data<-`(
    .@sam_data %>% 
      data.frame() %>% 
      dplyr::rename("ReadCount" = ReadCounts)) %>% 
  
  # remove taxa that are not present in any of the 'J' samples
  phyloseq::filter_taxa(
    function(x) sum(x) > 0,
    prune = TRUE
  ) %>%
  phyloseq::`taxa_names<-`(., .@tax_table@.Data[, 7])


phy_J_birth_rel <- phyloseq::transform_sample_counts(
  phy_J_birth, 
  function(OTU) OTU / sum(OTU)
) %>%
  phyloseq::transform_sample_counts(
    function(OTU) ifelse(is.na(OTU), 0, OTU)
  ) 


plot_data <- phy_J_birth_rel %>% 
  phyloseq::psmelt() %>% 
  dplyr::rename("Taxa" = OTU) %>% 
  dplyr::mutate(Taxa = stringr::str_replace(Taxa, " ", "\\."))

plot_theme()


# J = 75, 0 = 50, 1 = 25
title <- plot_data %>% 
  dplyr::distinct(ID, .keep_all = TRUE) %>% 
  dplyr::group_by(sub_grp, birth) %>% 
  dplyr::tally() %>%
  dplyr::mutate(
    birth = ifelse(
      birth == 1,
      "live birth",
      "no live birth")
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-sub_grp) %>% 
  tidyr::unite(col = "comb", sep = ": ") %>% 
  as.matrix() %>% 
  paste(., collapse = ", ")

# set conditions
ref <- sym("birth")

cnfr <- list(NULL, 
             c("age_strat", "g_strat", "p"))

names <- list("None", "PGA")

step <- 1

gt_dt_birth <- data.frame()

# Run tests and create plots
while (step <= length(cnfr)) {
  
  gt_res <- try({GroupTest_mod(
    species.table = phy_J_birth,
    meta = phy_J_birth,
    group = paste(ref), 
    group_name = paste(ref), 
    xl_titl = names[[step]],
    confounders = paste(cnfr[[step]]),
    dir_for_res = here::here("Results", "mare"), 
    min.prevalence = 0.1, 
    min.abundance = 0.1,
    outlier.cutoff = 3, 
    p.cutoff = 0.05, 
    keep.result = T, 
    nonzero = F,
    pdf = F, 
    show_quartz = F)
  }) 
  
  if (class(gt_res) == "try-error") {
    gt_res <- c()
  }
  
  gt_res_clean <- gt_res %>%
    dplyr::select(taxon, ends_with("_FDR")) %>%
    tidyr::pivot_longer(cols = -taxon, 
                        names_to = "comp", 
                        values_to = "pvals") %>%
    dplyr::group_by(taxon) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(variable = str_extract(
      comp, 
      pattern = paste(
        levels(
          as.factor(
            data.frame(phy_J_birth@sam_data)[,paste(ref)])), 
        collapse = "|"))) %>%
    dplyr::rename(!!ref := variable) %>%
    dplyr::select(-comp) %>%
    base::as.data.frame() %>%
    {
      p <- .
      if (!is.na(as.numeric(p[,paste(ref)][1])))
        p <- p %>%
          dplyr::mutate(!!ref := as.numeric(!!ref))
      else 
        .
    } %>%
    dplyr::rename("Taxa" = taxon) %>% 
    dplyr::mutate(test = names[[step]])
  
  gt_dt_birth <- gt_dt_birth %>% 
    dplyr::bind_rows(gt_res_clean)
  
  step <- step + 1 
}



## Clean and rearrange data
p2 <- plot_data %>%
  dplyr::select(Taxa, Abundance, !!ref) %>%
  dplyr::filter(Taxa %in% gt_dt_birth$Taxa) %>%
  dplyr::full_join(y = gt_dt_birth) %>%
  dplyr::mutate(
    Taxa = str_replace(Taxa, "\\.", " "), 
    sig = case_when((pvals > 0 & pvals <= 0.001) ~ "***",
                    (pvals > 0.001 & pvals <= 0.01) ~ "**",
                    (pvals > 0.01 & pvals <= 0.05) ~ "*",
                    (pvals > 0.05) ~ "")
  ) %>%
  dplyr::mutate(birth = ifelse(
    birth == 1, 
    "live birth",
    "no live birth"
  )) %>% 
  dplyr::mutate(!!ref := as.factor(!!ref))


## significance labels
p2_lab <- p2 %>%
  dplyr::filter(!is.na(sig)) %>% 
  dplyr::group_by(Taxa, !!ref, test) %>% 
  dplyr::distinct(Taxa, .keep_all = TRUE) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(
    pvals = signif(pvals, digits = 2),
    pvals = ifelse(pvals < 0.001, "<0.001", pvals), 
    pvals = ifelse(test == "None", 
                   paste0("<sup>1</sup>", pvals, sig), 
                   paste0("<sup>2</sup>", pvals, sig))) %>% 
  dplyr::group_by(Taxa, Abundance, birth) %>% 
  dplyr::summarise(pvals = paste0(pvals, collapse = "<br>"))


## assign colours for plot
pal <- c("#946666", "#66948e")


## construct viobox
gt_viobox <- ggplot2::ggplot(
  data = p2 %>%
    dplyr::filter(Abundance > 0), 
  aes(x = !!ref,
      y = Abundance,
      fill = !!ref)) +
  
  ## Add jitter plot
  ggplot2::geom_jitter(aes(colour = !!ref),
                       show.legend = F, 
                       size = 5,
                       alpha = 0.8) +
  
  ## half violin on the left
  gghalves::geom_half_violin(
    mapping = aes(colour = !!ref),
    trim = F, 
    side = "l",
    draw_quantiles = c(0.25, 0.5, 0.75), 
    scale = "area",
    #colour = "grey70",
    size = .5,
    alpha = 0.7) +
  
  ## half box on the right
  gghalves::geom_half_boxplot(
    mapping = aes(colour = !!ref),
    side = "r", 
    nudge = .05,
    size = .5, 
    #colour = "grey70",
    outlier.shape = NA, 
    alpha = 0.7,
    show.legend = F) +
  
  ## facet by Taxa
  ggplot2::facet_grid(. ~ Taxa, scales = "free_y") +
  ggtext::geom_richtext(
    data = p2_lab, 
    aes(x = !!ref,
        y = 20,
        label = pvals),
    colour = "grey60",
    size = 6,
    alpha = 0.85,
    inherit.aes = F, 
    show.legend = FALSE
  ) +
  
  ## add colours
  ggplot2::scale_fill_manual(values = pal, name = "") +
  ggplot2::scale_colour_manual(values = pal, name = "")  +
  
  ## log10 transform y-axis
  ggplot2::scale_y_log10() +
  
  ## add point for median
  ggplot2::stat_summary(
    fun = "mean", 
    geom = "point", 
    colour = "#c1e0dd", 
    size = 5, 
    show.legend = FALSE) +
  
  ## adjust labels and theme
  ggplot2::xlab("") +
  ggplot2::ylab("Relative abundance (log10)") +
  ggplot2::theme(
    axis.text.x = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(1, "cm"), 
    plot.title = element_text(colour = "grey50")) +
  ggplot2::guides(fill = guide_legend(ncol = 2))


gt_viobox


### compile and save infographic ----
comb_plot <- cowplot::plot_grid(
  polar_brt,
  stack_bar_brt_summ, 
  gt_viobox,
  nrow = 3, 
  rel_heights = c(1, 0.3, 0.6), 
  labels = "auto", 
  label_colour = "grey30", 
  label_size = 24
)


## compile and save infographic 
ggplot2::ggsave(
  path = here::here(pth_figrs),
  plot = comb_plot,
  filename = paste0(
    "Fig2.pdf"), 
  device = "pdf", 
  dpi = 1200,
  width = 35,
  height = 45,
  units = "cm", 
  limitsize = F)


ggplot2::ggsave(
  path = here::here(pth_figrs),
  plot = comb_plot,
  filename = paste0(
    "Fig2.eps"), 
  device = cairo_ps, 
  dpi = 1200,
  width = 35,
  height = 45,
  units = "cm",
  limitsize = F)




# Illustrations - J vs N ----
## Figure 3 A ----
### stacked bar plot ----
plt_dt_JN_prep <- plt_dt %>%
  dplyr::filter(prgcy == "pregnancy") %>% 
  dplyr::mutate(
    sub_grp = ifelse(sub_grp == "J", 
                     paste0("Fresh embryo\n", "transfer"), 
                     paste0("8th week\n",  "of gestation"))
  ) %>%
  dplyr::left_join(
    HEMI_scrs,
    by = "taxa"
  ) %>%
  dplyr::filter(
    !is.na(qual)
  ) %>%
  dplyr::mutate(
    Abundance = tidyr::replace_na(
      Abundance, 0
    )
  ) %>%
  dplyr::group_by(taxa) %>%
  dplyr::mutate(
    taxa = stringr::str_replace_all(
      taxa, "^c .|\\.\\.$", ""
    ),
    freq = sum(Abundance > 0.50),
  ) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(ID) %>% 
  dplyr::mutate(pairs = n_distinct(sub_grp)) %>% 
  dplyr::filter(pairs == 2)


# 21 (J) & 21 (N)
tally <- plt_dt_JN_prep %>% 
  dplyr::group_by(sub_grp, prgcy) %>%
  dplyr::distinct(Sample, .keep_all = TRUE) %>% 
  dplyr::count()


# factor levels for taxa
plt_dt_JN_tx <- plt_dt_JN_prep %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(
    taxa = if_else(
      Abundance < 0.02, "Others", taxa
    ),
    qual = ifelse(
      taxa == "Others", NA, qual
    ),
    tax_ord = ifelse(stringr::str_detect(
      taxa, "Lactobacillus"), 1, 0)
  ) %>%
  dplyr::group_by(ID) %>%
  dplyr::arrange(
    desc(prgcy),
    desc(tax_ord),
    desc(freq),
    taxa
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    taxa = factor(
      taxa,
      levels = unique(taxa)
    ),
    taxa = forcats::fct_relevel(taxa, "Others", after = Inf), 
  ) %>% 
  dplyr::arrange(taxa)


# create variable for sample levels
plt_dt_JN_lev <- plt_dt_JN_tx %>%
  dplyr::mutate(
    high = ifelse(Abundance >= 0.50,
                  sum(Abundance >= 0.50) / as.integer(taxa),
                  0
    ),
    mid = ifelse(Abundance >= 0.30 & Abundance < 0.51,
                 sum(Abundance >= 0.30) / as.integer(taxa),
                 0
    ),
    low = ifelse(Abundance >= 0.10 & Abundance < 0.30,
                 sum(Abundance >= 0.10) / as.integer(taxa),
                 0
    )
  )

# clean
plt_dt_JN_cln <- plt_dt_JN_lev %>%
  dplyr::arrange(
    desc(high),
    desc(mid),
    desc(low),
    taxa,
    desc(Abundance)  
    ) %>%
  dplyr::mutate(
    ID = factor(ID, levels = unique(ID)),
    sub_grp = factor(sub_grp, levels = unique(sub_grp)),
    prgcy = factor(prgcy, levels = unique(prgcy))
  ) %>%
  as.data.frame()

mycols_JN <- pick_colrs(plt_dt_JN_tx$taxa)

plot_theme(text_size = 14)

stack_comp_JN <- ggplot2::ggplot(
  plt_dt_JN_cln, 
  aes(
    x = ID, 
    y = Abundance, 
    fill = taxa
  )
) +
  ggplot2::facet_grid(sub_grp ~ birth, 
                      scales = "free", 
                      space = "free") +
  ggplot2::geom_bar(stat = "identity", 
                    position = position_stack(reverse = TRUE)) +
  ggplot2::scale_fill_manual(
    values = mycols_JN, 
    guide = guide_legend(
      ncol = 3)
  ) + 
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 270)
  ) +
  ggplot2::xlab("") +
  ggplot2::ylab("Relative Abundance")



## Figure 3 B ----
### stacked and summarised plot ----
sample_tally_JN <- plt_dt %>%
  dplyr::filter(prgcy == "pregnancy") %>% 
  dplyr::mutate(
    sub_grp = ifelse(sub_grp == "J", 
                     paste0("Fresh embryo\n", "transfer"), 
                     paste0("8th week\n",  "of gestation"))
  ) %>% 
  dplyr::group_by(ID) %>% 
  dplyr::mutate(pairs = n_distinct(sub_grp)) %>% 
  dplyr::filter(pairs == 2) %>% 
  dplyr::distinct(sub_grp, .keep_all = TRUE) %>% 
  dplyr::group_by(sub_grp) %>%
  dplyr::count()

plt_dt_stck_JN <- plt_dt_JN_prep %>%
  dplyr::group_by(taxa, sub_grp) %>%
  dplyr::summarise(Abundance = mean(Abundance)) %>%
  dplyr::mutate(
    taxa = str_remove(taxa, "X."), 
    taxa = str_replace(taxa, "\\.", " ") 
  ) %>%
  dplyr::ungroup() %>%
  #dplyr::mutate(taxa = ifelse(Abundance < 0.005, 
  #                            "Others", 
  #                            taxa)) %>% 
  dplyr::arrange(desc(Abundance)) %>%
  dplyr::mutate(
    taxa = factor(taxa, levels = levels(plt_dt_JN_tx$taxa)), 
    taxa = forcats::fct_relevel(taxa, "Others", after = Inf)
  )


# define species colours
mycols_JN_summ <- mycols_JN %>% 
  data.frame() %>%
  purrr::set_names("colours") %>%
  tibble::rownames_to_column(var = "taxa") %>%
  dplyr::mutate(lightness = shades::lightness(colours)) %>%
  dplyr::mutate(shades = if_else(lightness > 65, "light", "dark"))

if (all(mycols_JN_summ$colours == mycols_JN)) {
  print("identical colour palette - good to go")
} else {
  stop("colour palettes do not match")
}


plt_dt_stck_JN_wclrs <- plt_dt_stck_JN %>%
  dplyr::inner_join(mycols_JN_summ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(Abundance)) %>%
  dplyr::mutate(
    taxa = factor(taxa, levels = levels(plt_dt_JN_tx$taxa)),
    taxa = forcats::fct_relevel(taxa, "Others", after = Inf)
  ) %>%
  dplyr::mutate(sub_grp = factor(
    sub_grp,
    levels = unique(sort(sub_grp))
  )) %>%
  dplyr::mutate(shades = factor(shades, levels = c("dark", "light"))) %>%
  data.frame()


stack_bar_summ_JN <- ggplot2::ggplot(
  plt_dt_stck_JN_wclrs, aes(
    x = sub_grp,
    y = Abundance
  )
) +
  
  # add stacked bar
  ggplot2::geom_bar(aes(fill = taxa),
                    stat = "identity",
                    position = position_stack(reverse = T),
                    show.legend = FALSE
  ) +
  
  # add manual colour to bars
  ggplot2::scale_fill_manual(
    values = mycols_JN_summ$colours,
    name = "Species"
  ) +
  ggplot2::geom_text(
    aes(
      label = if_else(Abundance > 0.005,
                      as.character(round(Abundance * 100, digits = 2)), 
                      ""
      ),
      size = Abundance,
      colour = shades,
      group = taxa,
      angle = if_else(Abundance > 0.10, 0, 270)
    ),
    position = position_stack(0.5, reverse = T),
    show.legend = F
  ) +
  
  # add scale colour to text
  ggplot2::scale_colour_manual(
    values = c("grey85", "grey35"), 
    name = "Species") +
  
  # add sample tally
  ggplot2::geom_text(
    data = sample_tally_JN,
    aes(
      x = sub_grp,
      y = 1.05,
      label = paste("n = ", n, sep = "")
    ),
    colour = "grey40"
  ) +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 270, colour = "grey40")
  ) +
  ggplot2::xlab("") +
  ggplot2::ylab("Relative Abundance (%)") +
  ggplot2::coord_flip() +
  ggplot2::guides(fill = guide_legend(nrow = 5))



### grouptest ----
# (Only J samples, read-count cutoff 500)
phy_JN <- phy_HEMI %>%
  phyloseq::subset_samples(
    sub_grp %in% c("J", "N")
  ) %>%
  phyloseq::subset_samples(
    prgcy == 1
  ) %>% 
  phyloseq::`sample_data<-`(
    .@sam_data %>% 
      data.frame() %>% 
      dplyr::rename("ReadCount" = ReadCounts)) %>% 
  
  # remove taxa that are not present in any of the samples
  phyloseq::filter_taxa(
    function(x) sum(x) > 0,
    prune = TRUE
  ) %>%
  phyloseq::`taxa_names<-`(., .@tax_table@.Data[, 7])


ids <- data.frame(phy_JN@sam_data) %>% 
  dplyr::group_by(ID) %>% 
  dplyr::mutate(pairs = n_distinct(sub_grp)) %>% 
  dplyr::filter(pairs == 2)


phy_JN <- phy_JN %>% 
  phyloseq::subset_samples(
    ID %in% ids$ID
  )

test <- phy_JN@sam_data %>% 
  data.frame() %>% 
  dplyr::group_by(Sample_type) %>% 
  dplyr::tally()


if (all(test[2] == 21)) {
  print("Good to go")
} else {
  stop("Different sample pairs - check again")
}


phy_JN_rel <- phyloseq::transform_sample_counts(
  phy_JN, 
  function(OTU) OTU / sum(OTU)
) %>%
  phyloseq::transform_sample_counts(
    function(OTU) ifelse(is.na(OTU), 0, OTU)
  ) 



plot_data <- phy_JN_rel %>% 
  phyloseq::psmelt() %>% 
  dplyr::rename("Taxa" = OTU) %>% 
  dplyr::mutate(Taxa = stringr::str_replace(Taxa, " ", "\\."))

plot_theme()


# set conditions
ref <- sym("sub_grp")

cnfr <- list(NULL, 
             c("age_strat", "g_strat", "p"))

names <- list("None", "PGA")

step <- 1

gt_dt_JN <- data.frame()

# Run tests and create plots
while (step <= length(cnfr)) {
  
  gt_res_JN <- try({GroupTest_mod(
    species.table = phy_JN,
    meta = phy_JN,
    group = paste(ref), 
    group_name = paste(ref), 
    xl_titl = names[[step]],
    confounders = paste(cnfr[[step]]),
    dir_for_res = here::here("Results", "mare"), 
    min.prevalence = 0.1, 
    min.abundance = 0.1,
    outlier.cutoff = 3, 
    p.cutoff = 0.05, 
    keep.result = T, 
    nonzero = F,
    pdf = F, 
    show_quartz = F)
  }) 
  
  if (class(gt_res_JN) == "try-error") {
    gt_res_JN <- c()
  }
  
  gt_res_clean_JN <- gt_res_JN %>%
    dplyr::select(taxon, ends_with("_FDR")) %>%
    tidyr::pivot_longer(cols = -taxon, 
                        names_to = "comp", 
                        values_to = "pvals") %>%
    dplyr::group_by(taxon) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(variable = str_extract(
      comp, 
      pattern = paste(
        levels(
          as.factor(
            data.frame(phy_JN@sam_data)[,paste(ref)])), 
        collapse = "|"))) %>%
    dplyr::rename(!!ref := variable) %>%
    dplyr::select(-comp) %>%
    base::as.data.frame() %>%
    {
      p <- .
      if (!is.na(as.numeric(p[,paste(ref)][1])))
        p <- p %>%
          dplyr::mutate(!!ref := as.numeric(!!ref))
      else 
        .
    } %>%
    dplyr::rename("Taxa" = taxon) %>% 
    dplyr::mutate(test = names[[step]])
  
  gt_dt_JN <- gt_dt_JN %>% 
    dplyr::bind_rows(gt_res_clean_JN)
  
  step <- step + 1 
}

## Clean and rearrange data
p2_JN <- plot_data %>% 
  dplyr::select(Taxa, Abundance, !!ref) %>%
  dplyr::filter(Taxa %in% gt_dt_JN$Taxa) %>%
  dplyr::full_join(y = gt_dt_JN, relationship = "many-to-many") %>%
  dplyr::mutate(
    Taxa = str_replace(Taxa, "\\.", " "), 
    sig = case_when((pvals <= 0.001) ~ "***",
                    (pvals > 0.001 & pvals <= 0.01) ~ "**",
                    (pvals > 0.01 & pvals <= 0.05) ~ "*",
                    (pvals > 0.05) ~ "")
  ) %>%
  dplyr::mutate(sub_grp = ifelse(sub_grp == "J", 
                                 paste0("Fresh embryo\n", "transfer"), 
                                 paste0("8th week\n",  "of gestation") 
  )) %>% 
  dplyr::mutate(!!ref := as.factor(!!ref))


## significance labels
p2_lab_JN <- p2_JN %>%
  dplyr::filter(!is.na(sig)) %>% 
  dplyr::group_by(Taxa, !!ref, test) %>% 
  dplyr::distinct(Taxa, .keep_all = TRUE) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(
    pvals = signif(pvals, digits = 2),
    pvals = ifelse(pvals < 0.001, "<0.001", pvals), 
    pvals = ifelse(test == "None", 
                   paste0("<sup>1</sup>", pvals, sig), 
                   paste0("<sup>2</sup>", pvals, sig))) %>% 
  dplyr::group_by(Taxa, Abundance, !!ref) %>% 
  dplyr::summarise(pvals = paste0(pvals, collapse = "<br>"))


## assign colours for plot
pal <- c("#946666", "#66948e")


## construct viobox
gt_viobox_JN <- ggplot2::ggplot(
  data = p2_JN %>%
    dplyr::filter(Abundance > 0), 
  aes(x = !!ref,
      y = Abundance,
      fill = !!ref)) +
  
  ## Add jitter plot
  ggplot2::geom_jitter(aes(colour = !!ref),
                       show.legend = F, 
                       size = 5,
                       alpha = 0.8) +
  
  ## half violin on the left
  gghalves::geom_half_violin(
    mapping = aes(colour = !!ref),
    trim = F, 
    side = "l",
    draw_quantiles = c(0.25, 0.5, 0.75), 
    scale = "area",
    #colour = "grey70",
    size = .5,
    alpha = 0.7) +
  
  ## half box on the right
  gghalves::geom_half_boxplot(
    mapping = aes(colour = !!ref),
    side = "r", 
    nudge = .05,
    size = .5, 
    #colour = "grey70",
    outlier.shape = NA, 
    alpha = 0.7,
    show.legend = F) +
  
  ## facet by Taxa
  ggplot2::facet_grid(. ~ Taxa, scales = "free_y") +
  ggtext::geom_richtext(
    data = p2_lab_JN, 
    aes(x = !!ref,
        y = 20,
        label = pvals),
    colour = "grey60",
    size = 6,
    alpha = 0.85,
    inherit.aes = F, 
    show.legend = FALSE
  ) +
  
  ## add colours
  ggplot2::scale_fill_manual(values = pal, name = "") +
  ggplot2::scale_colour_manual(values = pal, name = "")  +
  
  ## log10 transform y-axis
  ggplot2::scale_y_log10() +
  
  ## add point for median
  ggplot2::stat_summary(
    fun = "mean", 
    geom = "point", 
    colour = "#c1e0dd", 
    size = 5, 
    show.legend = FALSE) +
  
  ## adjust labels and theme
  ggplot2::xlab("") +
  ggplot2::ylab("Relative abundance (log10)") +
  ggplot2::theme(
    axis.text.x = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(1, "cm"), 
    plot.title = element_text(colour = "grey50")) +
  ggplot2::guides(fill = guide_legend(ncol = 2))


gt_viobox_JN



### richness and diversity metrics ----
rh_div <- plt_dt_JN_cln %>% 
  dplyr::group_by(sub_grp) %>% 
  dplyr::distinct(ID, .keep_all = TRUE) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(desc(sub_grp), ID) %>%
  dplyr::mutate(sub_grp = factor(
    sub_grp, levels = unique(sub_grp)
  )) %>% 
  dplyr::ungroup()


pal <- c("#946666", "#66948e")


plt_rh <- ggpubr::ggpaired(
  rh_div, 
  x = "sub_grp", 
  y = "Richness", 
  #fill = "sub_grp", 
  add = "mean", 
  color = "sub_grp",
  line.color = "grey80",
  palette = pal, 
  ggtheme = plot_theme(angle = 0, 
                       axis.line.x = element_line(
                         colour = "grey50"), 
                       axis.line.y = element_line(
                         colour = "grey50")),
  add.params = list(color = "grey95"), 
  xlab = "Time of sampling", 
  ylab = "Richness", 
  id = "ID") +
  ggpubr::rremove("legend") +
  ggpubr::stat_compare_means(paired = TRUE)


plt_rh

plt_div <- ggpubr::ggpaired(
  rh_div, 
  x = "sub_grp", 
  y = "Diversity", 
  add = "mean", 
  color = "sub_grp",
  line.color = "grey80",
  palette = pal, 
  ggtheme = plot_theme(angle = 0, 
                       axis.line.x = element_line(
                         colour = "grey50"), 
                       axis.line.y = element_line(
                         colour = "grey50")),
  add.params = list(color = "grey95"), 
  xlab = "Time of sampling", 
  ylab = "Diversity", 
  id = "ID") +
  ggpubr::rremove("legend") +
  ggpubr::stat_compare_means(paired = TRUE)


plt_div


plot_theme()

plt_A <- cowplot::plot_grid(
  stack_comp_JN, 
  plt_rh,
  plt_div,
  ncol = 3, 
  rel_widths = c(1, 0.4, 0.4), 
  labels = "auto",
  label_colour = "grey30", 
  label_size = 24
)

plt_A

plt_B <- cowplot::plot_grid(
  stack_bar_summ_JN, 
  gt_viobox_JN, 
  nrow = 2, 
  rel_heights = c(0.4, 0.8), 
  labels = c("d", "e"),
  label_colour = "grey30", 
  label_size = 24
  )

plt_B



### compile and save infographic ----
comb_plot_JN <- cowplot::plot_grid(
  plt_A,
  plt_B, 
  nrow = 2
)

comb_plot_JN

## compile and save infographic 
ggplot2::ggsave(
  path = here::here(pth_figrs),
  plot = comb_plot_JN,
  filename = paste0(
    "Fig3.pdf"), 
  device = "pdf", 
  dpi = 1200,
  width = 45,
  height = 45,
  units = "cm", 
  limitsize = F)


ggplot2::ggsave(
  path = here::here(pth_figrs),
  plot = comb_plot_JN,
  filename = paste0(
    "Fig3.eps"), 
  device = cairo_ps, 
  dpi = 1200,
  width = 45,
  height = 45,
  units = "cm",
  limitsize = F)


## Supplement ----
### Supplement Figure 1 ----
alluvdt <- phy_HEMI_rel %>% 
  phyloseq::psmelt() %>%
  dplyr::arrange(Sample) %>%
  dplyr::rename("taxa" = "OTU") %>%
  dplyr::filter(!is.na(prgcy)) %>%
  dplyr::rowwise() %>% 
  #dplyr::filter(
  #  any(sub_grp == "N" & prgcy == 1) | sub_grp == "J"
  #) %>% 
  dplyr::filter(sub_grp == "J") %>% 
  dplyr::mutate(
    prgcy = ifelse(
      prgcy == 1,
      "pregnancy",
      "non-pregnancy"
    ),
    birth = ifelse(
      birth == 1, 
      "live birth", 
      "no live birth")) %>% 
  dplyr::select(ID, taxa, Abundance, Sample, birth, prgcy, sub_grp) %>% 
  dplyr::group_by(Sample) %>% 
  dplyr::filter(Abundance == max(Abundance)) %>% 
  #dplyr::filter(Abundance >= 0.50) %>% 
  dplyr::ungroup() %>% 
  #dplyr::group_by(taxa) %>% 
  dplyr::mutate(dom = ifelse(Abundance >= 0.50, "a", "b"), 
                tax_ord = ifelse(stringr::str_detect(taxa, 
                                                     "Lactobacillus"), 1, 0
                )) %>% 
  dplyr::group_by(Sample) %>% 
  dplyr::add_count() %>% 
  dplyr::group_by(prgcy) %>% 
  dplyr::add_tally(name = "prgcy_c") %>%
  dplyr::group_by(birth) %>% 
  dplyr::add_tally(name = "birth_c") %>% 
  dplyr::mutate(prgcy = paste0(prgcy, " (", prgcy_c, ")"), 
                birth = paste0(birth, " (", birth_c, ")")) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(desc(tax_ord), 
                 desc(Abundance)) %>% 
  dplyr::mutate(prgcy = factor(prgcy, levels = c("pregnancy (29)", "non-pregnancy (46)")), 
                birth = factor(birth, levels = c("live birth (25)", "no live birth (50)")), 
                taxa = factor(taxa, levels = unique(taxa)), 
                taxa = forcats::fct_relevel(taxa, "Lactobacillus crispatus", after = 0L)) %>% 
  #tidyr::pivot_longer(cols = c("prgcy", "birth"), 
  #                    names_to = "var", values_to = "vals") %>%
  dplyr::arrange(taxa) %>% 
  dplyr::mutate(ID = factor(ID, levels = unique(ID)))


tally <- alluvdt %>% 
  dplyr::group_by(prgcy, birth) %>% 
  dplyr::tally()


clrs <- pick_colrs(alluvdt$taxa)
plot_theme(text_size = 20)

alluv <- ggplot2::ggplot(
  alluvdt, 
  aes(y = n, 
      axis1 = taxa, 
      axis2 = birth, 
      axis3 = prgcy, 
      axis4 = ID)) +
  ggalluvial::stat_alluvium(aes(fill = taxa), alpha = 0.8) +
  ggalluvial::geom_lode(aes(fill = "taxa"), 
                        stat = "stratum", 
                        colour = "grey95", 
                        fill = "grey20", 
                        lwd = 1) +
  #ggalluvial::stat_alluvium(geom = "text", aes(label = sub_grp)) +
  ggfittext::geom_fit_text(stat = "stratum", 
                           aes(label = after_stat(stratum)), 
                           colour = "grey90", 
                           size = 25) +
  ggplot2::scale_fill_manual(values = clrs, name = "Species") +
  ggplot2::scale_x_continuous(
    breaks = 1:4, 
    labels = c("Most abundant species", 
               "Birth status", 
               "Pregnancy status", 
               "IDs")
  ) +
  ggplot2::theme(
    legend.position = "bottom"
  ) +
  ggplot2::guides(fill = guide_legend(nrow = 5)) +
  ggplot2::coord_flip()


ggplot2::ggsave(
  plot = alluv, 
  filename = "Supplement_Fig_1.pdf", 
  path = pth_figrs, 
  device = "pdf", 
  dpi = 1200,
  height = 60, 
  width = 70, 
  units = "cm"
)

ggplot2::ggsave(
  plot = alluv, 
  filename = "Supplement_Fig_1.eps", 
  path = pth_figrs, 
  device = cairo_ps, 
  dpi = 1200,
  height = 60, 
  width = 70, 
  units = "cm"
)


### Supplement Figure 2 ----
prgcy_rh_div <- plt_dt_plr_cln %>%
  dplyr::distinct(ID, .keep_all = TRUE) %>% 
  dplyr::select(ID, Richness, Diversity, prgcy) %>% 
  dplyr::arrange(prgcy)

pr1 <- ggpubr::ggboxplot(
  prgcy_rh_div, 
  x = "prgcy", 
  y = "Diversity", 
  fill = "prgcy", 
  color = "grey20",
  add = "mean", 
  palette = pal,
  ggtheme = plot_theme(angle = 0, 
                       axis.line.x = element_line(
                         colour = "grey50"), 
                       axis.line.y = element_line(
                         colour = "grey50")),
  add.params = list(color = "grey95")) +
  ggpubr::rremove("legend") +
  ggpubr::stat_compare_means()

pr1


pr2 <- ggpubr::ggboxplot(
  prgcy_rh_div, 
  x = "prgcy", 
  y = "Richness", 
  fill = "prgcy", 
  color = "grey20",
  add = "mean", 
  palette = pal,
  ggtheme = plot_theme(angle = 0, 
                       axis.line.x = element_line(
                         colour = "grey50"), 
                       axis.line.y = element_line(
                         colour = "grey50")),
  add.params = list(color = "grey95")) +
  ggpubr::rremove("legend") +
  ggpubr::stat_compare_means()

pr2


### compile and save infographic ----
comb_plot_pr <- cowplot::plot_grid(
  pr1,
  pr2, 
  nrow = 1,
  labels = "auto",
  label_colour = "grey30", 
  label_size = 24
)

comb_plot_pr

## compile and save infographic 
ggplot2::ggsave(
  path = here::here(pth_figrs),
  plot = comb_plot_pr,
  filename = paste0(
    "Supplement_Fig_2.pdf"), 
  device = "pdf", 
  dpi = 1200,
  width = 25,
  height = 15,
  units = "cm", 
  limitsize = F)


ggplot2::ggsave(
  path = here::here(pth_figrs),
  plot = comb_plot_pr,
  filename = paste0(
    "Supplement_Fig_2.eps"), 
  device = cairo_ps, 
  dpi = 1200,
  width = 25,
  height = 15,
  units = "cm",
  limitsize = F)
