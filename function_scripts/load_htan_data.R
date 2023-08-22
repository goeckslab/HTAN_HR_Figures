


###################################################################
#
#   FUNCTIONS FOR LOADING SOURCE DATA FOR HTAN HR+ MANUSCRIPT
#
###################################################################

library(data.table)



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






