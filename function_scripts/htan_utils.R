

###################################################################
#
#   FUNCTIONS FOR LOADING SOURCE DATA FOR HTAN HR+ MANUSCRIPT
#
###################################################################

library(data.table)
library(dplyr)



# Function to load sample meta table 
load_meta <- function(fn = 'meta.csv', fn.dir = '') {
  
  # Factor order for sample IDs
  lvls <- c("HTA9-1_Bx1", "HTA9-1_Bx2", "HTA9-1_Bx4", "HTA9-1_Bx5", 
            "HTA9-2_Bx1", "HTA9-2_Bx2", 
            "HTA9-3_Bx1", "HTA9-3_Bx2", 
            "HTA9-14_Bx1", "HTA9-14_Bx2", "HTA9-14_Bx3", 
            "HTA9-15_Bx1", "HTA9-15_Bx2")
  
  # Read in meta table
  df <- fread(paste(fn.dir, fn, sep = '/')) %>%
    mutate(Sample = factor(Sample, levels = lvls)) %>%
    arrange(Sample) %>%
    data.frame()
  rownames(df) <- df$Sample
  
  return(df)
  
}

# Function to load CNV calls 
load_cnvs <- function(meta, fn = 'cnvs.csv', fn.dir = '') {
  
  # Load CNV calls
  df <- fread(paste(fn.dir, fn, sep = '/')) %>%
    data.frame(row.names = 1, check.names = FALSE) 
  df <- df[,meta$Sample]
  
  return(df)
  
}

# Function to load SNV calls 
load_snvs <- function(meta, fn = 'snvs.csv', fn.dir = '') {
  
  # Load SNV calls
  df <- fread(paste(fn.dir, fn, sep = '/')) %>%
    data.frame(row.names = 1, check.names = FALSE) 
  df <- df[,meta$Sample]
  
  return(df)
  
}

# Function to load GSVA pathway activity (RNA-seq) enrichment scores
load_gsva <- function(meta, fn = 'gsva_scores.csv', fn.dir = '') {
  
  # Load pathway scores
  df <- fread(paste(fn.dir, fn, sep = '/')) %>%
    data.frame(row.names = 1, check.names = FALSE) 
  df <- df[,meta$Sample]
  
  return(df)
  
}

# Function to load proteomic pathway activity (RPPA) scores
load_ppws <- function(meta, fn = 'protein_pathway_scores.csv', fn.dir = '') {
  
  # Load pathway scores
  df <- fread(paste(fn.dir, fn, sep = '/')) %>%
    data.frame(row.names = 1, check.names = FALSE) 
  
  select_samples <- as.character(meta$Sample[meta$Sample %in% colnames(df)])
  
  df <- df[,select_samples]
  
  return(df)
  
}


# Function to load gene expression (select categories, pre-scaled to background MMTERT cohort)
load_gene_expression <- function(meta, fn = 'gene_expression_scaled_categories.csv', fn.dir = '') {
  
  # Load expression
  df <- fread(paste(fn.dir, fn, sep = '/')) %>%
    select(-Category) %>%
    data.frame(row.names = 1, check.names = FALSE)
  df <- df[,meta$Sample]
  
  return(df)
  
}

# Function to load viper activity (select categories, pre-scaled to background MMTERT cohort)
load_viper <- function(meta, fn = 'viper_scaled_categories.csv', fn.dir = '') {
  
  # Load viper
  df <- fread(paste(fn.dir, fn, sep = '/')) %>%
    select(-Category) %>%
    data.frame(row.names = 1, check.names = FALSE)
  df <- df[,meta$Sample]
  
  return(df)
  
}

# Function to load RPPA protein abundace (pre-scaled to background MMTERT cohort)
load_htan_rppa <- function(meta, fn = 'rppa_scaled.csv', fn.dir = '') {
  
  # Load rppa
  df <- fread(paste(fn.dir, fn, sep = '/')) %>%
    data.frame(row.names = 1, check.names = FALSE)
  
  select_samples <- as.character(meta$Sample[meta$Sample %in% colnames(df)])
  
  df <- df[,select_samples]
  
  return(df)
  
}


load_protein_to_rna <- function(fn = 'protein_to_rna.csv', fn.dir = '') {
  
  fread(file.path(fn.dir, fn)) %>%
    return()
  
}


load_merged_rna_to_protein_names <- function(fn = 'merged_rna_protein_names.csv', fn.dir = '') {
  
  fread(file.path(fn.dir, fn)) %>%
    return()
  
}
