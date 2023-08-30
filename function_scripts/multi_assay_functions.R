








# Function for integrating change across multiple RNA and protein assays
merge_assays <- function(select_samples, 
                         meta, 
                         df.rna = NULL, 
                         df.viper = NULL, 
                         df.rppa = NULL,
                         protein_rna_table = NULL, 
                         select_gene_cats = NULL, 
                         fill_all_assays = FALSE,
                         Zchange = TRUE, 
                         select_assays = NULL) {
  
  
  # Subset meta table
  meta <- meta %>% 
    select(Sample, Patient, HTAN, HtanID, Date, Biopsy) %>% 
    distinct() %>% 
    arrange(Date)
  
  
  
  # Store assay tables here
  df.list <- list()
  
  # Convert to long format and compute change; return difference using second sample
  if (!is.null(df.rna)) {
    
    if (Zchange) {
      
      df.rna <- df.rna %>% as.matrix() %>% melt() %>% 
        setNames(c('Gene', 'Sample', 'Zscore')) %>% 
        filter(Sample %in% select_samples,
               Gene %in% select_gene_cats$Gene) %>% 
        left_join(meta) %>% 
        group_by(Patient, Gene) %>%
        arrange(Date, .by_group = TRUE) %>% 
        mutate(Count = n()) %>% 
        filter(Count >= 2) %>%
        mutate(Zchange = Zscore - lag(Zscore)) %>%
        ungroup() %>%
        filter(!is.na(Zchange)) %>%
        select(Sample, Gene, Zchange) %>% 
        mutate(Assay = 'RNA') 
      
    } else {
      
      df.rna <- df.rna %>% as.matrix() %>% melt() %>% 
        setNames(c('Gene', 'Sample', 'Zscore')) %>% 
        filter(Sample %in% select_samples) %>% 
        select(Sample, Gene, Zscore) %>% 
        mutate(Assay = 'RNA') 
      
    }
    
    df.list[['RNA']] <- df.rna
    
  }
  
  # Format Viper
  if (!is.null(df.viper)) {
    
    # Scale viper activity
    df.viper <- t(scale(t(df.viper),scale = TRUE, center = TRUE))
    
    if (Zchange) {
      
      df.viper <- df.viper %>% as.matrix() %>% melt() %>% 
        setNames(c('Gene', 'Sample', 'Zscore')) %>% 
        filter(Sample %in% select_samples,
               Gene %in% select_gene_cats$Gene) %>% 
        left_join(meta) %>% 
        group_by(Patient, Gene) %>%
        arrange(Date, .by_group = TRUE) %>% 
        mutate(Count = n()) %>% 
        filter(Count >= 2) %>%
        mutate(Zchange = Zscore - lag(Zscore)) %>%
        ungroup() %>%
        filter(!is.na(Zchange)) %>%
        select(Sample, Gene, Zchange) %>% 
        mutate(Assay = 'Viper') 
      
    } else {
      
      df.viper <- df.viper %>% as.matrix() %>% melt() %>% 
        setNames(c('Gene', 'Sample', 'Zscore')) %>% 
        filter(Gene %in% select_gene_cats$Gene,
               Sample %in% select_samples) %>% 
        select(Sample, Gene, Zscore) %>% 
        mutate(Assay = 'Viper') 
      
    }
    
    df.list[['Viper']] <- df.viper
    
  }
  
  
  # Scale RPPA and split protein and phosphoprotein
  if (!is.null(df.rppa)) {
    
    # Select proteins available
    select_proteins <- protein_rna_table %>% 
      filter(RNA %in% select_gene_cats$Gene) %>% 
      pull(Protein) %>% as.character()  
    
    # Subset down to select proteins
    df.rppa <- df.rppa[rownames(df.rppa) %in% select_proteins,]
    
    # Scale RPPA data
    df.rppa <- t(scale(t(df.rppa),scale = TRUE, center = TRUE)) 
    
    
    # Split proteins and phosphoproteins
    df.protein <- df.rppa[grep('_p', rownames(df.rppa), invert = TRUE),,drop = FALSE]
    df.phospho <- df.rppa[grep('_p', rownames(df.rppa), invert = FALSE),,drop = FALSE]
    
    if (nrow(df.protein) > 0) {
      
      # Melt to long, then merge with protein/RNA table to get gene names
      df.protein <- df.protein %>% as.matrix() %>% melt() %>% 
        setNames(c('Protein', 'Sample', 'Zscore')) %>%
        left_join(protein_rna_table) %>% 
        filter(RNA != "") %>% distinct() %>%
        select(Sample, RNA, Protein, Zscore)
      
      # Compute change between pairs  
      if (Zchange) {
        
        df.protein <- df.protein %>% 
          filter(Sample %in% select_samples) %>% 
          left_join(meta) %>% 
          group_by(Patient, Protein, RNA) %>%
          arrange(Date, .by_group = TRUE) %>% 
          mutate(Count = n()) %>% 
          filter(Count >= 2) %>%
          mutate(Zchange = Zscore - lag(Zscore)) %>%
          ungroup() %>%
          filter(!is.na(Zchange)) %>%
          select(Sample, Protein, RNA, Zchange) %>% 
          setNames(c('Sample', 'Protein', 'Gene', 'Zchange')) %>%
          mutate(Assay = 'Protein') 
        
        # Otherwise just show value
      } else {
        
        df.protein <- df.protein %>% 
          filter(Sample %in% select_samples) %>% 
          select(Sample, Protein, RNA, Zscore) %>% 
          setNames(c('Sample', 'Protein', 'Gene', 'Zscore')) %>% 
          mutate(Assay = 'Protein') 
        
      }
      
      df.list[['Protein']] <- df.protein
      
    }
    
    # Now phosphoproteins
    if (nrow(df.phospho) > 0) {
      
      # Melt to long, then merge with protein/RNA table to get gene names
      df.phospho <- df.phospho %>% as.matrix() %>% melt() %>% 
        setNames(c('Protein', 'Sample', 'Zscore')) %>% 
        left_join(protein_rna_table) %>% 
        filter(RNA != "") %>% distinct() %>%
        select(Sample, RNA, Protein, Zscore)
      
      # Compute change between pairs
      if (Zchange) {
        
        df.phospho <- df.phospho %>% 
          filter(Sample %in% select_samples) %>% 
          left_join(meta) %>% 
          group_by(Patient, Protein, RNA) %>%
          arrange(Date, .by_group = TRUE) %>% 
          mutate(Count = n()) %>% 
          filter(Count >= 2) %>%
          mutate(Zchange = Zscore - lag(Zscore)) %>%
          ungroup() %>%
          filter(!is.na(Zchange)) %>%
          select(Sample, Protein, RNA, Zchange) %>% 
          setNames(c('Sample', 'Phospho', 'Gene', 'Zchange')) %>%
          mutate(Assay = 'Phospho') 
        
        
        # Otherwise use value  
      } else {
        
        df.phospho <- df.phospho %>% 
          filter(Sample %in% select_samples) %>% 
          select(Sample, Protein, RNA, Zscore) %>% 
          setNames(c('Sample', 'Phospho', 'Gene', 'Zscore')) %>% 
          mutate(Assay = 'Phospho') 
        
      }
      
      df.list[['Phospho']] <- df.phospho
      
    }
    
  }
  
  
  
  
  
 
  
  
  
  
  
  
  
  ##########################################
  
  # NEW: CONVERT GENE NAMES TO MERGED GENE/PROTEIN/PHOSPHOPROTEIN NAMES
  
  if ( 'Protein' %in% names(df.list) ) {
    
    # Compute merged gene/protein/phosphoprotein names from protein_rna_table
    #  First, filter protein/rna table for only proteins in sample data
    if (!is.null(df.list[['Protein']])) {
      
      select_proteins <- df.list[['Protein']] %>%
        pull(Protein) %>% as.character() %>%
        unique() %>% sort()
      
    } else {
      
      select_proteins <- c()
      
    }
    
    if (!is.null(df.list[['Phospho']])) {
      
      select_phosphos <- df.list[['Phospho']] %>%
        pull(Phospho) %>% as.character() %>%
        unique() %>% sort()
      
    } else {
      
      select_phosphos <- c()
      
    }
    
    protein_rna_tbl <- protein_rna_tbl %>% 
      filter(RNA %in% select_gene_cats$Gene) %>%
      filter(Protein %in% c(select_proteins, select_phosphos))
    
    merged_names <- merge_rna_protein_table(protein_rna_tbl) %>% arrange(Protein, Gene)
    
    # Update all gene names from each assay dataframe to new merged gene/protein names
    if (Zchange) {value.var <- 'Zchange'} else {value.var <- 'Zscore'}
    select_cols <- c('Sample', 'FullName', value.var, 'Assay')
    
    for (n in names(df.list)) {
      
      if (nrow(df.list[[n]])) {
        
        df.list[[n]] <- df.list[[n]] %>% 
          left_join(merged_names) %>%
          mutate(FullName = case_when(is.na(FullName) ~ Gene,
                                      TRUE ~ FullName)) %>%
          select(all_of(select_cols)) %>%
          setNames(c('Sample', 'Gene', value.var, 'Assay')) 
        
      }
      
    }
    
    # Get new names for select gene categories
    select_gene_cats <- select_gene_cats %>% 
      left_join(merged_names) %>%
      mutate(FullName = case_when(is.na(FullName) ~ Gene,
                                  TRUE ~ FullName)) %>%
      select(FullName, Category) %>%
      setNames(c('Gene', 'Category')) %>%
      data.frame()
    
  }
  
  
  
  ##########################################
  
  # Merge dataframe
  df.merged <- do.call(rbind, df.list)
  
  
  # Now create gene by assay matrices for each sample and store in list
  assay_mats <- list()
  
  # Get samples
  select_samples <- as.character(unique(df.merged$Sample))
  
  # Get genes (only genes that are found in at least one sample)
  all_genes <- as.character(unique(select_gene_cats$Gene))
  
  df.merged <- df.merged %>% filter(Gene %in% all_genes)
  
  all_genes <- all_genes[all_genes %in% as.character(df.merged$Gene)]
  
  # Order for acast?
  if (!is.null(select_assays)) {
    all_assays <- select_assays
  } else {
    all_assays <- c('RNA', 'Viper', 'Protein', 'Phospho')
  }
  
  df.merged$Assay <- factor(df.merged$Assay, levels = all_assays)
  
  for (s in select_samples) {
    
    sub.df <- df.merged %>% filter(Sample == s) 
    
    # Add to list of matrices
    if (Zchange) {value.var <- 'Zchange'} else {value.var <- 'Zscore'}
    assay_mats[[s]] <- acast(sub.df, Gene ~ Assay, value.var = value.var)
    
    if (fill_all_assays) {
      
      present <- colnames(assay_mats[[s]])
      no_dat <- all_assays[!all_assays %in% present]
      na_mat <- matrix(NA, nrow = nrow(assay_mats[[s]]), ncol = length(no_dat), dimnames = list(rownames(assay_mats[[s]]), no_dat))
      assay_mats[[s]] <- cbind(assay_mats[[s]], na_mat)[,all_assays]
      
    }
    
  }
  
  return(assay_mats)
  
}  
