


###################################################################
#
#   FUNCTIONS FOR LOADING SOURCE DATA FOR HTAN HR+ MANUSCRIPT
#
###################################################################

library(data.table)



# Function to load sample meta table 
load_meta <- function(fn = 'meta.csv', fn.dir = '') {
  
  df <- fread(paste(fn.dir, fn, sep = '/')) %>%
    mutate(Sample = factor(Sample, levels = c("HTA9-1_Bx1", "HTA9-1_Bx2", "HTA9-1_Bx4", "HTA9-1_Bx5", "HTA9-2_Bx1", 
                                              "HTA9-2_Bx2", "HTA9-3_Bx1", "HTA9-3_Bx2", "HTA9-14_Bx1", "HTA9-14_Bx2", 
                                              "HTA9-14_Bx3", "HTA9-15_Bx1", "HTA9-15_Bx2"))) %>%
    arrange(Sample) %>%
    data.frame()
  rownames(df) <- df$Sample
  
  return(df)
  
}


# Function to load CNV calls 
load_cnvs <- function(meta, fn.cnvs = 'cnvs.csv', fn.dir = '') {
  
  # Load CNV calls
  df <- fread(paste(fn.dir, fn.cnvs, sep = '/')) %>%
    data.frame(row.names = 1, check.names = FALSE) 
  df <- df[,meta$Sample]
  
  return(df)
  
}

# Function to load SNV calls 
load_snvs <- function(meta, fn.snvs = 'snvs.csv', fn.dir = '') {
  
  # Load snv calls
  df <- fread(paste(fn.dir, fn.snvs, sep = '/')) %>%
    data.frame(row.names = 1, check.names = FALSE) 
  df <- df[,meta$Sample]
  
  return(df)
  
}



