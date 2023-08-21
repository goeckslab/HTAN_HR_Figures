


###################################################################
#
#   FUNCTIONS FOR LOADING SOURCE DATA FOR HTAN HR+ MANUSCRIPT
#
###################################################################

library(data.table)



# Function to load sample meta table 
load_meta <- function(fn = 'meta.csv', fn.dir = '') {
  
  df <- fread(paste(fn.dir, fn, sep = '/')) %>%
    data.frame()
  rownames(df) <- df$Sample
  
  return(df)
  
}


# Function to load CNV and SNV calls 
load_htan_variants <- function(meta, fn.cnvs = 'cnvs.csv', fn.svs = 'snvs.csv', fn.dir = '') {
  
  # Load CNV calls
  df.cnvs <- fread(paste(fn.dir, fn.cnvs, sep = '/')) 
    
  # Load SNV calls
  df.snvs <- fread(paste(fn.dir, fn.svs, sep = '/'))
  
  return(list('cnvs' = df.cnvs,
              'snvs' = df.snvs))
    
}





