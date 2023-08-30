








# Function for integrating change across multiple RNA and protein assays
merge_assays <- function(select_samples, 
                         meta_table, 
                         df.rna = NULL, 
                         df.viper = NULL, 
                         df.rppa = NULL,
                         protein_rna_table = NULL, 
                         select_gene_cats = NULL, 
                         fill_all_assays = FALSE,
                         Zchange = TRUE, 
                         select_assays = NULL) {
  
  
  
  
  
  
  # Subset meta table
  meta_table <- meta_table %>% 
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
        left_join(meta_table) %>% 
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
        left_join(meta_table) %>% 
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
          left_join(meta_table) %>% 
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
          left_join(meta_table) %>% 
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
  
  
  
  # Scale clinicalIHC
  # NO LONGER FUNTIONAL!!!
  if (!is.null(df.clinicalIHC)) {
    
    df.clinicalIHC <- df.clinicalIHC[rownames(df.clinicalIHC) %in% select_gene_cats$Gene,,drop=FALSE]
    
    df.clinicalIHC <- t(scale(t(df.clinicalIHC),scale = TRUE, center = TRUE))
    
    if (Zchange) {
      
      df.clinicalIHC <- df.clinicalIHC %>% as.matrix() %>% melt() %>% 
        mutate_at(vars(value), ~replace(., is.nan(.), 0)) %>% 
        setNames(c('Gene', 'Sample', 'Zscore')) %>% 
        filter(Gene %in% select_gene_cats$Gene) %>%
        filter(Sample %in% select_samples) %>% left_join(meta_table) %>% 
        group_by(Patient) %>% filter(n() == 2*nrow(df.clinicalIHC)) %>% ungroup() %>%
        group_by(Patient, Gene) %>% arrange(Date) %>% mutate(Zchange = diff(Zscore)) %>% 
        dplyr::slice(n()) %>% ungroup() %>% select(Sample, Gene, Zchange) %>% mutate(Assay = 'cIHC') 
      
    } else {
      
      df.clinicalIHC <- df.clinicalIHC %>% as.matrix() %>% melt() %>% 
        mutate_at(vars(value), ~replace(., is.nan(.), 0)) %>% 
        setNames(c('Gene', 'Sample', 'Zscore')) %>% 
        filter(Gene %in% select_gene_cats$Gene) %>%
        filter(Sample %in% select_samples) %>% 
        select(Sample, Gene, Zscore) %>% mutate(Assay = 'cIHC') 
      
    }
    
    df.list[['cIHC']] <- df.clinicalIHC
    
  }
  
  # Convert to long format and compute change; return difference using second sample
  # NO LONGER FUNTIONAL!!!
  if (!is.null(df.cosmx)) {
    
    df.cosmx <- df.cosmx[rownames(df.cosmx) %in% select_gene_cats$Gene,]
    
    # NOT TESTED
    if (Zchange) {
      
      df.cosmx <-  df.cosmx %>% as.matrix() %>% melt() %>% setNames(c('Gene', 'Sample', 'Zscore')) %>% 
        filter(Sample %in% select_samples) %>% left_join(meta_table) %>% 
        group_by(Patient) %>% filter(n() == 2*nrow(df.cosmx)) %>% ungroup() %>%
        group_by(Patient, Gene) %>% arrange(Date) %>% mutate(Zchange = diff(Zscore)) %>% 
        dplyr::slice(n()) %>% ungroup() %>% select(Sample, Gene, Zchange) %>% mutate(Assay = 'CosMX_RNA') 
      
    } else {
      
      df.cosmx <- df.cosmx %>% as.matrix() %>% melt() %>% setNames(c('Gene', 'Sample', 'Zscore')) %>% 
        filter(Sample %in% select_samples) %>% select(Sample, Gene, Zscore) %>% mutate(Assay = 'CosMX_RNA') 
      
    }
    
    df.list[['CosMX_RNA']] <- df.cosmx
    
  }
  
  # cycIF
  # NO LONGER FUNTIONAL!!!
  if (!is.null(df.cycIF)) {
    
    select_cycIF_markers <- cycIF_rna_tbl %>% filter(RNA %in% select_gene_cats$Gene) %>% 
      pull(Protein) %>% as.character()
    
    df.cycIF <- df.cycIF[rownames(df.cycIF) %in% select_cycIF_markers,,drop=FALSE]
    
    df.cycIF <- t(scale(t(df.cycIF),scale = TRUE, center = TRUE))
    
    if (Zchange) {
      
      df.cycIF <- df.cycIF %>% as.matrix() %>% melt() %>% 
        mutate_at(vars(value), ~replace(., is.nan(.), 0)) %>%
        setNames(c('Protein', 'Sample', 'Zscore')) %>% 
        left_join(cycIF_rna_tbl) %>%
        select(Sample, RNA, Zscore) %>% distinct() %>%
        setNames(c('Sample', 'Gene', 'Zscore')) %>%
        filter(Sample %in% select_samples) %>% left_join(meta_table) %>% 
        group_by(Patient) %>% filter(n() == 2*nrow(df.cycIF)) %>% ungroup() %>%
        group_by(Patient, Gene) %>% arrange(Date) %>% mutate(Zchange = diff(Zscore)) %>% 
        dplyr::slice(n()) %>% ungroup() %>% select(Sample, Gene, Zchange) %>% mutate(Assay = 't-CycIF') 
      
    } else {
      
      df.cycIF <- df.cycIF %>% as.matrix() %>% melt() %>% 
        mutate_at(vars(value), ~replace(., is.nan(.), 0)) %>%
        setNames(c('Protein', 'Sample', 'Zscore')) %>% 
        left_join(cycIF_rna_tbl) %>%
        select(Sample, RNA, Zscore) %>% distinct() %>%
        setNames(c('Sample', 'Gene', 'Zscore')) %>%
        filter(Sample %in% select_samples) %>% 
        select(Sample, Gene, Zscore) %>% mutate(Assay = 't-CycIF') 
      
    }
    
    df.list[['t-CycIF']] <- df.cycIF
    
  }
  
  # mIHC
  # NO LONGER FUNTIONAL!!!
  if (!is.null(df.mIHC)) {
    
    select_mIHC_markers <- mIHC_rna_tbl %>% filter(RNA %in% select_gene_cats$Gene) %>% 
      pull(Protein) %>% as.character()
    
    df.mIHC <- df.mIHC[rownames(df.mIHC) %in% select_mIHC_markers,,drop=FALSE]
    
    df.mIHC <- t(scale(t(df.mIHC),scale = TRUE, center = TRUE))
    
    if (Zchange) {
      
      df.mIHC <- df.mIHC %>% as.matrix() %>% melt() %>% 
        mutate_at(vars(value), ~replace(., is.nan(.), 0)) %>%
        setNames(c('Protein', 'Sample', 'Zscore')) %>% 
        left_join(mIHC_rna_tbl, by = 'Protein') %>%
        select(Sample, RNA, Zscore) %>% distinct() %>%
        setNames(c('Sample', 'Gene', 'Zscore')) %>%
        filter(Sample %in% select_samples) %>% left_join(meta_table) %>% 
        group_by(Patient) %>% filter(n() == 2*nrow(df.mIHC)) %>% ungroup() %>%
        group_by(Patient, Gene) %>% arrange(Date) %>% mutate(Zchange = diff(Zscore)) %>% 
        dplyr::slice(n()) %>% ungroup() %>% select(Sample, Gene, Zchange) %>% mutate(Assay = 'mIHC') 
      
    } else {
      
      df.mIHC <- df.mIHC %>% as.matrix() %>% melt() %>% 
        mutate_at(vars(value), ~replace(., is.nan(.), 0)) %>%
        setNames(c('Protein', 'Sample', 'Zscore')) %>% 
        left_join(mIHC_rna_tbl, by = 'Protein') %>%
        select(Sample, RNA, Zscore) %>% distinct() %>%
        setNames(c('Sample', 'Gene', 'Zscore')) %>%
        filter(Sample %in% select_samples) %>% 
        select(Sample, Gene, Zscore) %>% mutate(Assay = 'mIHC') 
      
    }
    
    df.list[['mIHC']] <- df.mIHC
    
  }
  
  
  # Mutation scores
  # NO LONGER FUNTIONAL!!!
  if (!is.null(df.mutation_scores)) {
    
    df.mutation_scores <- df.mutation_scores[rownames(df.mutation_scores) %in% select_gene_cats$Gene,,drop = FALSE]
    
    if (Zchange) {
      
      df.mutation_scores <- df.mutation_scores %>% as.matrix() %>% melt() %>% 
        setNames(c('Gene', 'Sample', 'Zscore')) %>% 
        filter(Sample %in% select_samples) %>% left_join(meta_table) %>% 
        group_by(Patient) %>% filter(n() == 2*nrow(df.mutation_scores)) %>% ungroup() %>%
        group_by(Patient, Gene) %>% arrange(Date) %>% mutate(Zchange = diff(Zscore)) %>% 
        dplyr::slice(n()) %>% ungroup() %>% select(Sample, Gene, Zchange) %>% mutate(Assay = 'Mutation Signature') 
      
    } else {
      
      df.mutation_scores <- df.mutation_scores %>% as.matrix() %>% melt() %>% 
        setNames(c('Gene', 'Sample', 'Zscore')) %>% 
        filter(Sample %in% select_samples) %>% 
        select(Sample, Gene, Zscore) %>% mutate(Assay = 'Mutation Signature') 
      
    }
    
    df.list[['Mutation Signature']] <- df.mutation_scores
    
  }
  
  
  # Pancancer mutation scores
  # NO LONGER FUNTIONAL!!!
  if (!is.null(df.pancancer)) {
    
    df.pancancer <- df.pancancer[rownames(df.pancancer) %in% select_gene_cats$Gene,,drop = FALSE]
    
    if (Zchange) {
      
      df.pancancer <- df.pancancer %>% as.matrix() %>% melt() %>% 
        setNames(c('Gene', 'Sample', 'Zscore')) %>% 
        filter(Sample %in% select_samples) %>% left_join(meta_table) %>% 
        group_by(Patient) %>% filter(n() == 2*nrow(df.pancancer)) %>% ungroup() %>%
        group_by(Patient, Gene) %>% arrange(Date) %>% mutate(Zchange = diff(Zscore)) %>% 
        dplyr::slice(n()) %>% ungroup() %>% select(Sample, Gene, Zchange) %>% mutate(Assay = 'pancancer') 
      
    } else {
      
      df.pancancer <- df.pancancer %>% as.matrix() %>% melt() %>% 
        setNames(c('Gene', 'Sample', 'Zscore')) %>% 
        filter(Sample %in% select_samples) %>% 
        select(Sample, Gene, Zscore) %>% mutate(Assay = 'pancancer') 
      
    }
    
    df.list[['pancancer']] <- df.pancancer
    
  }
  
  # BRCA mutation scores
  # NO LONGER FUNTIONAL!!!
  if (!is.null(df.brca)) {
    
    df.brca <- df.brca[rownames(df.brca) %in% select_gene_cats$Gene,,drop = FALSE]
    
    if (Zchange) {
      
      df.brca <- df.brca %>% as.matrix() %>% melt() %>% 
        setNames(c('Gene', 'Sample', 'Zscore')) %>% 
        filter(Sample %in% select_samples) %>% left_join(meta_table) %>% 
        group_by(Patient) %>% filter(n() == 2*nrow(df.brca)) %>% ungroup() %>%
        group_by(Patient, Gene) %>% arrange(Date) %>% mutate(Zchange = diff(Zscore)) %>% 
        dplyr::slice(n()) %>% ungroup() %>% select(Sample, Gene, Zchange) %>% mutate(Assay = 'BRCA') 
      
    } else {
      
      df.brca <- df.brca %>% as.matrix() %>% melt() %>% 
        setNames(c('Gene', 'Sample', 'Zscore')) %>% 
        filter(Sample %in% select_samples) %>% 
        select(Sample, Gene, Zscore) %>% mutate(Assay = 'BRCA') 
      
    }
    
    df.list[['BRCA']] <- df.brca
    
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
    all_assays <- c('RNA', 'Viper', 'pancancer', 'BRCA', 'Protein', 'Phospho', 'cIHC', 'mIHC', 't-CycIF', 'CosMX_RNA')
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
