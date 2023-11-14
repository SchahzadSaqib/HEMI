# Session parameters and packages ----
old_path <- Sys.getenv("PATH")
old_path
Sys.setenv(PATH = paste(old_path, "~/opt/anaconda3/bin", sep = ":"))

## load packages ----
pkg_to_load <- list.files(
  here::here(
    "Scripts",
    "to_load"
  ),
  pattern = "^1_|^3_",
  full.names = T
)

pkg_to_load <- purrr::map(
  pkg_to_load,
  source
)

## Set NCBI key ----
rentrez::set_entrez_key("insertNCBIKEY")
Sys.getenv("ENTREZ_KEY")

## define major directory paths ----
path_to_data <- here::here(
  "Data"
)

path_to_anlys <- here::here(
  "Data",
  "Analysis"
)

path_to_NCBI_db <- here::here(
  "~",
  "Documents",
  "NCBI_databases"
)

path_to_alt_db <- here::here(
  "~",
  "Documents",
  "Annot_sets"
)

path_to_ecolst <- here::here(
  "~",
  "Documents",
  "Eco_lists"
)

path_to_rslt <- here::here(
  "Results",
  "16S",
  "MiSeq47"
)


# 16S rRNA pipeline ----
## Remove Ns from the sequences ----
fnFs_MiSeq47 <- sort(
  list.files(
    here::here(
      path_to_data,
      "Raw",
      "MiSeq47"
    ),
    pattern = "_R1",
    full.names = T
  )
)

fnRs_MiSeq47 <- sort(
  list.files(
    here::here(
      path_to_data,
      "Raw",
      "MiSeq47"
    ),
    pattern = "_R2",
    full.names = T
  )
)

sample.names <- c(sapply(
  strsplit(basename(fnFs_MiSeq47), "_"),
  function(x) paste(x[1:3], collapse = "_")
))
sample.names

# Filter Ns and place in separate directory
Nrmv_Fs_MiSeq47 <- file.path(
  here::here(
    path_to_data,
    "rmv_Ns",
    "MiSeq47"
  ),
  paste0(sample.names, "_R1_Nrmv.fastq.gz")
)

Nrmv_Rs_MiSeq47 <- file.path(
  here::here(
    path_to_data,
    "rmv_Ns",
    "MiSeq47"
  ),
  paste0(sample.names, "_R2_Nrmv.fastq.gz")
)

names(Nrmv_Fs_MiSeq47) <- sample.names
names(Nrmv_Rs_MiSeq47) <- sample.names

outN_MiSeq47 <- dada2::filterAndTrim(
  fwd = fnFs_MiSeq47,
  filt = Nrmv_Fs_MiSeq47,
  rev = fnRs_MiSeq47,
  filt.rev = Nrmv_Rs_MiSeq47,
  #trimLeft = c(0, 2),
  maxN = 0,
  compress = T,
  multithread = T,
  verbose = T
)

head(outN_MiSeq47)
rownames(outN_MiSeq47) <- sample.names

Nrmv_Fs_MiSeq47 <- sort(
  list.files(
    here::here(
      path_to_data,
      "rmv_Ns",
      "MiSeq47"
    ),
    pattern = "_R1_",
    full.names = T
  )
)

Nrmv_Rs_MiSeq47 <- sort(
  list.files(
    here::here(
      path_to_data,
      "rmv_Ns",
      "MiSeq47"
    ),
    pattern = "_R2_",
    full.names = T
  )
)

sample.names <- c(sapply(
  strsplit(basename(Nrmv_Fs_MiSeq47), "_"),
  function(x) paste(x[1:3], collapse = "_")
))
sample.names

# Quality profiles
Qualp_Nrmv_R1_MiSeq47 <- dada2::plotQualityProfile(
  Nrmv_Fs_MiSeq47,
  aggregate = T
) # Quality of forward reads
Qualp_Nrmv_R1_MiSeq47

Qualp_Nrmv_R2_MiSeq47 <- dada2::plotQualityProfile(
  Nrmv_Rs_MiSeq47,
  aggregate = T
) # Quality of reverse reads
Qualp_Nrmv_R2_MiSeq47


## Cutadapt to split libraries ----
system.time(
  split.res <- splitLibrariesByPrimers(
    fnFs = list.files(
      path = here::here(
        path_to_data,
        "rmv_Ns",
        "MiSeq47"
      ),
      pattern = "R1_Nrmv.fastq.gz",
      full.names = T
    ), # Path for raw data
    fnRs = list.files(
      path = here::here(
        path_to_data,
        "rmv_Ns",
        "MiSeq47"
      ),
      pattern = "R2_Nrmv.fastq.gz",
      full.names = T
    ),
    path.16S = here::here(
      path_to_data,
      "trim",
      "16S",
      "MiSeq47"
    ),
    path.ITS = here::here(
      path_to_data,
      "trim",
      "ITS",
      "MiSeq47"
    )
  )
)

Cutadpt_Fs_MiSeq47 <- sort(
  list.files(
    here::here(
      path_to_data,
      "trim",
      "16S",
      "MiSeq47"
    ),
    pattern = "_R1_",
    full.names = T
  )
)

Cutadpt_Rs_MiSeq47 <- sort(
  list.files(
    here::here(
      path_to_data,
      "trim",
      "16S",
      "MiSeq47"
    ),
    pattern = "_R2_", full.names = T
  )
)

sample.names <- c(sapply(
  strsplit(basename(Cutadpt_Fs_MiSeq47), "_"),
  function(x) paste(x[1:3], collapse = "_")
))
sample.names


## Filter, trim, and place in separate directory ----
QC_trim_Fs_MiSeq47 <- file.path(
  here::here(
    path_to_data,
    "clean",
    "16S",
    "MiSeq47"
  ),
  paste0(sample.names, "_R1_16S_QC_trim.fastq.gz")
)

QC_trim_Rs_MiSeq47 <- file.path(
  here::here(
    path_to_data,
    "clean",
    "16S",
    "MiSeq47"
  ),
  paste0(sample.names, "_R2_16S_QC_trim.fastq.gz")
)

names(QC_trim_Fs_MiSeq47) <- sample.names
names(QC_trim_Rs_MiSeq47) <- sample.names

out_QC_trim_MiSeq47 <- dada2::filterAndTrim(
  fwd = Cutadpt_Fs_MiSeq47,
  filt = QC_trim_Fs_MiSeq47,
  rev = Cutadpt_Rs_MiSeq47,
  filt.rev = QC_trim_Rs_MiSeq47,
  maxN = 0,
  truncQ = c(2, 4),
  maxEE = c(2, 4),
  truncLen = c(240, 220),
  minLen = 150,
  rm.phix = T,
  compress = T,
  multithread = T,
  matchIDs = T,
  verbose = T
)

head(out_QC_trim_MiSeq47)

QC_trim_Fs_MiSeq47 <- sort(
  list.files(
    here::here(
      path_to_data,
      "clean",
      "16S",
      "MiSeq47"
    ),
    pattern = "_R1_",
    full.names = T
  )
)

QC_trim_Rs_MiSeq47 <- sort(
  list.files(
    here::here(
      path_to_data,
      "clean",
      "16S",
      "MiSeq47"
    ),
    pattern = "_R2_",
    full.names = T
  )
)

sample.names <- c(sapply(
  strsplit(basename(QC_trim_Fs_MiSeq47), "_"),
  function(x) paste(x[1:3], collapse = "_")
))
sample.names

rownames(out_QC_trim_MiSeq47) <- sample.names

# Check Primer removal
Qualp_QC_trim_R1_MiSeq47 <- dada2::plotQualityProfile(
  QC_trim_Fs_MiSeq47,
  aggregate = T
)

Qualp_QC_trim_R1_MiSeq47

Qualp_QC_trim_R2_MiSeq47 <- dada2::plotQualityProfile(
  QC_trim_Rs_MiSeq47,
  aggregate = T
)

Qualp_QC_trim_R2_MiSeq47

gridExtra::grid.arrange(
  nrow = 2, Qualp_Nrmv_R1_MiSeq47,
  Qualp_Nrmv_R2_MiSeq47,
  Qualp_QC_trim_R1_MiSeq47,
  Qualp_QC_trim_R2_MiSeq47
)

g <- gridExtra::arrangeGrob(
  Qualp_Nrmv_R1_MiSeq47,
  Qualp_Nrmv_R2_MiSeq47,
  Qualp_QC_trim_R1_MiSeq47,
  Qualp_QC_trim_R2_MiSeq47,
  nrow = 2,
  ncol = 2
)

ggplot2::ggsave(g,
  path = here::here("Results"),
  filename = paste("MiSeq47_QC_aggregate",
    ".pdf",
    sep = ""
  ),
  device = "pdf",
  width = 40,
  height = 30,
  units = "cm",
  dpi = 320
)


## Learn the Error rates ----
errors_F <- dada2::learnErrors(
  QC_trim_Fs_MiSeq47,
  randomize = T,
  multithread = T,
  verbose = T
)

errors_R <- dada2::learnErrors(
  QC_trim_Rs_MiSeq47,
  randomize = T,
  multithread = T,
  verbose = T
)

plot_err_F <- dada2::plotErrors(
  errors_F,
  nominalQ = T
)
plot_err_F

ggplot2::ggsave(
  path = here::here("Results"),
  filename = "MiSeq47_error_rates_F.pdf",
  plot_err_F,
  device = "pdf",
  dpi = 320
)

plot_err_R <- dada2::plotErrors(
  errors_R,
  nominalQ = T
)
plot_err_R

ggplot2::ggsave(
  path = here::here("Results"),
  filename = "MiSeq47_error_rates_R.pdf",
  plot_err_R,
  device = "pdf",
  dpi = 320
)
# Black: Observed errors
# Red: Expected error rate


## Sample inference dada2 ----
SInf_Fs <- dada2::dada(
  QC_trim_Fs_MiSeq47,
  err = errors_F,
  multithread = T,
  verbose = T
)

SInf_Rs <- dada2::dada(
  QC_trim_Rs_MiSeq47,
  err = errors_R,
  multithread = T,
  verbose = T
)

## Merge paired reads ----
mergers <- dada2::mergePairs(
  SInf_Fs,
  QC_trim_Fs_MiSeq47,
  SInf_Rs,
  QC_trim_Rs_MiSeq47,
  verbose = T
)
head(mergers[[1]])

## Construct sequence table ----
seqtab_MiSeq47 <- dada2::makeSequenceTable(mergers)
dim(seqtab_MiSeq47)
table(nchar(getSequences((seqtab_MiSeq47))))


## Remove Chimeras ----
seqtab.nochim_MiSeq47 <- dada2::removeBimeraDenovo(
  seqtab_MiSeq47,
  minFoldParentOverAbundance = 8,
  method = "consensus",
  multithread = T,
  verbose = T
)

dim(seqtab.nochim_MiSeq47)
sum(seqtab.nochim_MiSeq47) / sum(seqtab_MiSeq47)
rownames(seqtab.nochim_MiSeq47) <- sample.names

## Tracking reads through the pipeline ----
getN <- function(x) sum(getUniques(x))
track <- outN_MiSeq47 %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "samples") %>%
  dplyr::inner_join(
    out_QC_trim_MiSeq47 <- data.frame(out_QC_trim_MiSeq47) %>%
      tibble::rownames_to_column(var = "samples"),
    by = "samples"
  ) %>%
  dplyr::inner_join(
    SInf_Fs <- data.frame(sapply(SInf_Fs, getN)) %>%
      dplyr::mutate(samples = sample.names),
    by = "samples"
  ) %>%
  dplyr::inner_join(
    SInf_Rs <- data.frame(sapply(SInf_Rs, getN)) %>%
      dplyr::mutate(samples = sample.names),
    by = "samples"
  ) %>%
  dplyr::inner_join(
    mergers <- data.frame(sapply(mergers, getN)) %>%
      dplyr::mutate(samples = sample.names),
    by = "samples"
  ) %>%
  dplyr::inner_join(
    nochim <- data.frame(rowSums(seqtab.nochim_MiSeq47)) %>%
      dplyr::mutate(samples = sample.names),
    by = "samples"
  ) %>%
  purrr::set_names(c(
    "samples",
    "inputN",
    "filteredN",
    "inputF",
    "filteredF",
    "denoisedF",
    "denoisedR",
    "merged",
    "nonchim"
  )) %>%
  dplyr::arrange(desc(nonchim)) %>%
  dplyr::mutate(Pct_filt = nonchim / filteredF * 100)


# Taxonomic annotations - 16S ----
## taxminer ----
### Alignments 
Blast_hits_C <- taxminer::txm_align(
  seq_in = seqtab.nochim_MiSeq47,
  db_path = here::here(
    path_to_NCBI_db,
    "bac",
    "c"
  ),
  db_name = "bacC",
  tab_path = path_to_rslt,
  tab_out = "MiSeq47",
  alt_annot = TRUE,
  alt_path = path_to_alt_db,
  asbd_tbl_lge = here::here(
    path_to_ecolst,
    "bac_lineage.fst"
  ),
  run_blst = TRUE,
  threads = 3,
  batches = 3,
  chunks = 3,
  qcvg = 98,
  pctidt = 98,
  max_out = 5000
)

Blast_hits_UC <- taxminer::txm_align( 
  seq_in = seqtab.nochim_MiSeq47,
  db_path = here::here(
    path_to_NCBI_db,
    "bac",
    "uc"
  ),
  db_name = "bacUC",
  tab_path = path_to_rslt,
  tab_out = "MiSeq47_UC",
  asbd_tbl_lge = here::here(
    path_to_ecolst,
    "bac_lineage.fst"
  ),
  run_blst = TRUE,
  alt_annot = FALSE,
  threads = 3,
  batches = 3,
  chunks = 3,
  qcvg = 98,
  pctidt = 98,
  max_out = 5000
)

## Text-mining based filtration of hits
hits_w_ecosrcC <- Blast_hits_C %>%
  taxminer::txm_ecosrc(
    filt_host = "human",
    filt_site = c("vagina+FRS", "gut+oral+clinical"),
    filt_negt = "non_human",
    asbd_tbl = here::here(
      path_to_ecolst,
      "bac_ecosys.fst"
    ),
    alt_tbl_path = path_to_rslt,
    alt_tbl_name = "MiSeq47",
    do_filt = TRUE,
    add_scrs = TRUE
  )
  

Annot_MiSeq47_clnC <- hits_w_ecosrcC %>% 
  dplyr::group_by(ID) %>%
  dplyr::mutate(
    TaxID = as.numeric(TaxID), 
    ord = ifelse(
      str_starts(
      species, "unclassified"), 
      0, 1)) %>% 
  dplyr::filter(ord == max(ord)) %>% 
  dplyr::filter(Evalue == min(Evalue)) %>%
  dplyr::filter(bitscore == max(bitscore)) %>%
  dplyr::filter(score == max(score)) %>%
  dplyr::filter(Pct == max(Pct)) %>%
  dplyr::filter(qcovs == max(qcovs)) %>%
  dplyr::distinct(species, .keep_all = T) %>%
  dplyr::mutate(across(where(is.list), as.character)) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  base::as.data.frame() 

hits_w_ecosrcUC <- Blast_hits_UC %>%
  dplyr::filter(!ID %in% Annot_MiSeq47_clnC$ID) %>%
  taxminer::txm_ecosrc(
    filt_host = "human",
    filt_site = c("vagina+FRS", "gut+oral+clinical"),
    filt_negt = "non_human",
    asbd_tbl = here::here(
      path_to_ecolst,
      "bac_ecosys.fst"
    ),
    alt_tbl_path = path_to_rslt,
    alt_tbl_name = "MiSeq47",
    do_filt = TRUE, 
    add_scrs = TRUE
  )


Annot_MiSeq47_clnUC <- hits_w_ecosrcUC %>% 
  dplyr::group_by(ID) %>%
  dplyr::mutate(TaxID = as.numeric(TaxID)) %>%
  dplyr::filter(bitscore == max(bitscore)) %>%
  dplyr::filter(score == max(score)) %>%
  dplyr::filter(Pct == max(Pct)) %>%
  dplyr::filter(qcovs == max(qcovs)) %>%
  dplyr::filter(Evalue == min(Evalue)) %>%
  dplyr::distinct(species, .keep_all = T) %>%
  dplyr::mutate(across(where(is.list), as.character)) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  base::as.data.frame()

## Combining annotations
Annot_MiSeq47 <- Annot_MiSeq47_clnC %>%
  dplyr::bind_rows(Annot_MiSeq47_clnUC) %>%
  dplyr::arrange(nchar(ID))

saveRDS(
  Annot_MiSeq47, 
  here::here(
    path_to_anlys, 
    "1_MiSeq47_HEMI_wEco.rds"
  ))

## score summaries
score_summ <- Annot_MiSeq47 %>%
  dplyr::select(species:score) %>%
  dplyr::group_by(species) %>%
  dplyr::filter(
    silva_sp == max(silva_sp) |
      RDP_sp == max(RDP_sp)
  ) %>%
  dplyr::filter(
    silva == max(silva) |
      RDP == max(RDP)
  ) %>%
  dplyr::mutate(
    qual = case_when(
      silva_sp == 2 & RDP_sp == 2 ~ "a",
      silva == 6 | RDP == 6 ~ "b",
      silva < 6 & RDP < 6 ~ "c",
      TRUE ~ "d"
    )
  ) %>%
  dplyr::distinct(species, .keep_all = T)

saveRDS(
  score_summ,
  here::here(
    path_to_anlys,
    "1_MiSeq47_HEMI_scr_summ.rds"
  )
)

## Create taxonomy table
taxonomy <- Annot_MiSeq47 %>%
  tibble::column_to_rownames(var = "Seq") %>%
  dplyr::select(
    superkingdom,
    phylum,
    class,
    order,
    family,
    genus,
    species
  ) %>%
  dplyr::rename_with(tolower) %>%
  as.matrix() 

tax_chck <- taxonomy %>%
  data.frame() %>% 
  dplyr::group_by(across(superkingdom:genus)) %>% 
  dplyr::distinct(species, .keep_all = TRUE) %>% 
  dplyr::group_by(species) %>% 
  dplyr::add_tally(name = "sp")

if (any(tax_chck$sp > 1)) {
  stop("Potential duplicate taxa")
}

# Filter out sequences that were removed
MiSeq47_seqtab_filt <- seqtab.nochim_MiSeq47 %>%
  base::data.frame() %>%
  dplyr::select(Annot_MiSeq47$Seq) %>%
  base::as.matrix()

base::all.equal(
  rownames(taxonomy),
  colnames(MiSeq47_seqtab_filt)
)

saveRDS(taxonomy,
  file = here::here(
    path_to_anlys,
    "1_MiSeq47_HEMI_taxonomy.rds"
  )
)

saveRDS(MiSeq47_seqtab_filt,
  file = here::here(
    path_to_anlys,
    "1_MiSeq47_HEMI_seqtab_filt.rds"
  )
)

Align_test <- Blast_hits_C %>%
  dplyr::bind_rows(
    dplyr::filter(
      Blast_hits_UC,
      !ID %in% Blast_hits_C$ID
    )
  ) %>%
  dplyr::arrange(nchar(ID))

seq_prefilt <- seqtab.nochim_MiSeq47 %>%
  base::data.frame() %>%
  dplyr::select(Align_test$Seq) %>%
  dplyr::mutate(reads_aligned = rowSums(.)) %>%
  tibble::rownames_to_column("samples") %>%
  dplyr::select(samples, reads_aligned)

## Add annotated read percentage to tracking
track_cmpl <- data.frame(MiSeq47_seqtab_filt) %>%
  dplyr::mutate(reads_annotated = rowSums(.)) %>%
  tibble::rownames_to_column("samples") %>%
  dplyr::select(samples, reads_annotated) %>%
  dplyr::inner_join(track) %>%
  dplyr::inner_join(seq_prefilt) %>%
  dplyr::relocate(reads_annotated, .after = last_col()) %>%
  dplyr::mutate(
    Pct_aligned = reads_aligned / nonchim * 100,
    Pct_annot = reads_annotated / reads_aligned * 100
  ) %>%
  dplyr::arrange(
    str_detect(
      samples,
      "\\+"
    ),
    desc(Pct_annot)
  )

writexl::write_xlsx(track_cmpl,
  path = here::here(
    path_to_rslt,
    "MiSeq47_readtrack_annot.xlsx"
  )
)


# save image ----
save.image(here::here(
  path_to_anlys,
  "1_MiSeq47_HEMI_preprs_taxannot.RData"
))
