# Create phylogenetic trees ----
## https://f1000research.com/articles/5-1492/v2
mk_trees <- function(x) {
print(paste0("creating tree"))
  
## Constructing a phylogenetic tree
seqs <- phyloseq::taxa_names(x)
names(seqs) <- phyloseq::taxa_names(x)
fitGTR <- DECIPHER::AlignSeqs(DNAStringSet(seqs),
                              anchor = NA
) %>%
  base::as.matrix() %>%
  phangorn::phyDat(type = "DNA") %>%
  base::assign("align", ., envir = .GlobalEnv) %>%
  phangorn::dist.ml() %>%
  phangorn::NJ() %>%
  phangorn::pml(data = align) %>%
  stats::update(k = 4, inv = 0.2) %>%
  phangorn::optim.pml(
    model = "GTR",
    optInv = T,
    optGamma = T,
    rearrangement = "stochastic",
    control = phangorn::pml.control(trace = 0)
  )
}


# Verify read counts ----
check_reads <- function(x) {
  if (identical(
    as.vector(rowSums(data.frame(x@otu_table))),
    as.vector(x@sam_data$ReadCounts)
  )) {
    print(paste0("reads match"))
  } else {
    stop(paste0("reads are not identical"))
  }
}

