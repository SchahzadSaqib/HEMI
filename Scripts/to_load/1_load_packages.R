if (!require("pacman")) install.packages("pacman")
if (!require("remotes")) install.packages("remotes")


pacman::p_load(phyloseq, 
               here, 
               writexl, 
               tidyverse,
               readr,
               vegan,
               dada2,
               metacoder,
               rentrez,
               summarytools,
               ggpubr,
               ggExtra,
               ggnewscale,
               ggExtra,
               gghalves,
               shades,
               geomtextpath,
               DECIPHER,
               phangorn,
               scales,
               MASS,
               cowplot,
               scico,
               PNWColors,
               shades,
               viridis, 
               KEGGREST)

pacman::p_install_gh("SchahzadSaqib/taxminer")
pacman::p_install_gh("AllanCameron/geomtextpath")
