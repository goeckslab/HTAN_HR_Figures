

library(dplyr)


# Function that merges names regardless if all features are available
#  Similar to function below, but was needed for some reason for plotting
#  Might be able to merge the two?
make_merged_rna_protein_categories <- function(assay_mats, gene_cats, protein_rna_table) {
  
  # Get all genes/proteins across assay mattrices
  all_features <- Reduce(union, lapply(assay_mats, FUN = row.names)) %>% unique() %>% sort()
  
  # Find names when both protein and phosphoprotein available
  rna_protein_phospho <- gene_cats %>% 
    left_join(merge_rna_protein_table(protein_rna_tbl))
  
  # Find names when only protein available
  rna_protein <- gene_cats %>% 
    left_join(merge_rna_protein_table(protein_rna_tbl %>% 
                                        filter(!grepl('_p', Protein))))
  
  # Find names when only phosphoprotein available
  rna_phospho <- gene_cats %>% 
    left_join(merge_rna_protein_table(protein_rna_tbl %>% 
                                        filter(grepl('_p', Protein))))
  
  # Merge all together
  gene_cats <- rbind(rbind(rna_protein_phospho, rna_protein), rna_phospho) %>% 
    mutate(FullName = case_when(is.na(FullName) ~ Gene,
                                TRUE ~ FullName)) %>%
    filter(FullName %in% all_features) %>%
    select(FullName, Category) %>%
    setNames(c('Gene', 'Category')) %>%
    distinct()
  
  rownames(gene_cats) <- gene_cats$Gene
  
  return(gene_cats)
  
}

# Function to create new merged gene/protein names when integrating 
#  genes, proteins, and phosphoproteins across assays
merge_rna_protein_table <- function(protein_rna_tbl) {
  
  # Remove proteins/phosphoproteins without a matching gene 
  protein_rna_tbl <- protein_rna_tbl %>% 
    setNames(c('Protein', 'RNA', 'Category')) %>%
    filter(RNA != '') %>% 
    distinct()
  
  # Separate proteins and phosphoproteins
  protein.tbl <- protein_rna_tbl[grep('_p', protein_rna_tbl$Protein, invert = TRUE),,drop = FALSE]
  phospho.tbl <- protein_rna_tbl[grep('_p', protein_rna_tbl$Protein, invert = FALSE),,drop = FALSE] %>% 
    setNames(c('Phospho', 'RNA', 'Category'))
  
  # Merge back together with proteins and phosphoproteins now as separate columns
  merge_gene_names <- full_join(protein.tbl, phospho.tbl) %>%
    distinct() %>% 
    select(Category, RNA, Protein, Phospho) %>%
    arrange(RNA, Phospho)
  
  # Separate phosphoprotein protein names and actual phospho groups
  merge_gene_names <- merge_gene_names %>%
    separate(col = 'Phospho', into = c('Phospho', 'Phosphogroup'), sep = '_', extra = 'merge')
  
  # Loop through and create new merged names
  merge_gene_names$FullName <- ''
  
  for (i in 1:nrow(merge_gene_names)) {
    
    gene <- merge_gene_names[i,2]
    protein <- merge_gene_names[i,3]
    phospho <- merge_gene_names[i,4]
    pGroup <- merge_gene_names[i,5]
    
    # Merge genes and proteins without phospho data
    if (is.na(phospho)) {
      
      if (gsub('-', '', str_to_lower(gene)) == gsub('-', '', str_to_lower(protein))) {
        
        merge_gene_names[i,'FullName'] <- gene
        
      } else {
        
        merge_gene_names[i,'FullName'] <- paste(gene, protein, sep = '/')
        
      }
      
      # Merge genes/phosphoproteins withoutprotein data    
    } else if (is.na(protein)) {
      
      if (gsub('-', '', str_to_lower(gene)) == gsub('-', '', str_to_lower(phospho))) {
        
        merge_gene_names[i,'FullName'] <- paste(gene, pGroup, sep = '/')
        
      } else {
        
        merge_gene_names[i,'FullName'] <- paste0(gene, '/', phospho, '_', pGroup)
        
      }
      
      # Merge genes/proteins with phospho data    
    } else {
      
      if (gsub('-', '', str_to_lower(gene)) == gsub('-', '', str_to_lower(protein)) &
          gsub('-', '', str_to_lower(gene)) == gsub('-', '', str_to_lower(phospho))) {
        
        merge_gene_names[i,'FullName'] <- paste(gene, pGroup, sep = '/')
        
      } else if (gsub('-', '', str_to_lower(protein)) == gsub('-', '', str_to_lower(phospho))) {
        
        merge_gene_names[i,'FullName'] <- paste0(gene, '/', protein, '_', pGroup)
        
      } else {
        
        merge_gene_names[i,'FullName'] <- paste0(gene, '/', protein, '/', phospho, '_', pGroup)
        
      }
      
    }
    
  }
  
  # Manually fill some others
  merge_gene_names[merge_gene_names$FullName == 'AURKA/Aurora-A/Aurora-ABC_pT288_pT232_pT198', 'FullName'] <- 'AURKA_pT288_pT232_pT198'
  merge_gene_names[merge_gene_names$FullName == 'BAK1/BAK','FullName'] <- 'BAK1'
  merge_gene_names[merge_gene_names$FullName == 'PRKCA/PKCALPHA/PKC-a-b-II_pT638_T641','FullName'] <- 'PRKCA/PKC-a-b-II_pT638_T641'
  merge_gene_names[merge_gene_names$FullName == 'STAT5A/STAT5ALPHA','FullName'] <- 'STAT5A'
  merge_gene_names[merge_gene_names$FullName == 'GZMB/Granzyme-B','FullName'] <- 'GZMB'
  merge_gene_names[merge_gene_names$FullName == 'RPS6KB1/P70S6K1/P70S6K_pT389','FullName'] <- 'RPS6KB1/P70S6K1_pT389'
  merge_gene_names[merge_gene_names$FullName == 'CCND1/CYCLIND1','FullName'] <- 'CCND1'
  merge_gene_names[merge_gene_names$FullName == 'CCND3/Cyclin-D3','FullName'] <- 'CCND3'
  merge_gene_names[merge_gene_names$FullName == 'CCNE1/CYCLINE1','FullName'] <- 'CCNE1'
  merge_gene_names[merge_gene_names$FullName == 'CCNB1/CYCLINB1','FullName'] <- 'CCNB1'
  merge_gene_names[merge_gene_names$FullName == 'CHEK1/CHK1_pS296','FullName'] <- 'CHK1_pS296'
  merge_gene_names[merge_gene_names$FullName == 'CHEK1/CHK1_pS345','FullName'] <- 'CHK1_pS345'
  merge_gene_names[merge_gene_names$FullName == 'CHEK2/CHK2_pT68','FullName'] <- 'CHK2_pT68'
  merge_gene_names[merge_gene_names$FullName == 'EIF4EBP1/X4EBP1_pS65','FullName'] <- 'EIF4EBP1/4EBP1_pS65'
  merge_gene_names[merge_gene_names$FullName == 'EIF4EBP1/X4EBP1_pT37T46','FullName'] <- 'EIF4EBP1/4EBP1_pT37T46'
  
  # Merge phosphoproteins back to original names
  merge_gene_names <- merge_gene_names %>% 
    mutate(Phospho = case_when(!is.na(Phospho) ~ paste(Phospho, Phosphogroup, sep = '_'),
                               TRUE ~ Phospho)) %>%
    select(RNA, Protein, Phospho, FullName) %>%
    setNames(c('Gene', "Protein", 'Phospho', 'FullName'))
  
  return(merge_gene_names)
  
}


# Function to subset assay matrix by select samples and features
# TODO: May move to separate functions script
subset_assay_data <- function(df, select_samples = NULL, select_features = NULL) {
  
  if (!is.null(select_samples)) {
    
    df <- df[,colnames(df) %in% select_samples,drop=FALSE]
    
  }
  
  if (!is.null(select_features)) {
    
    df <- df[rownames(df) %in% select_features,,drop=FALSE]
    
  }
  
  return(df)
  
}


# Function to convert wide matrix into long format for multi-assay analysis
add_assay <- function(df, assay, meta = NULL, 
                      feature_type = 'Gene', Zchange = TRUE, 
                      select_samples = NULL, select_features = NULL) {
  
  # Get samples
  if (is.null(select_samples)) {select_samples <- colnames(df)}
  
  # Get genes/proteins
  if (is.null(select_features)) {select_features <- rownames(df)}
  
  # Compuate change in value between paired samples
  if (Zchange) {
    
    df <- df %>% 
      compute_paired_change(meta, 
                            select_features = select_features,
                            select_samples = select_samples,
                            return_matrix = FALSE) %>%
      setNames(c('Sample', feature_type, 'Zchange')) %>%
      mutate(Assay = assay) 
    
    # Else, format to long using relative values  
  } else {
    
    df <- df %>% 
      as.matrix() %>% 
      melt() %>% 
      setNames(c(feature_type, 'Sample', 'Zscore')) %>% 
      select(all_of(c('Sample', feature_type, 'Zscore'))) %>% 
      filter(.data[[feature_type]] %in% select_features,
             Sample %in% select_samples) %>%
      mutate(Assay = assay) 
    
  }
  
  return(df)
  
}

# Function for integrating change across multiple RNA and protein assays
merge_assays <- function(meta, 
                         df.rna = NULL, 
                         df.viper = NULL, 
                         df.rppa = NULL,
                         scale_rna = FALSE,
                         scale_viper = FALSE,
                         scale_rppa = FALSE,
                         protein_rna_table = NULL, 
                         select_gene_cats = NULL,   # IS THERE ACTUALLY A POINT TO select_gene_cats ANYMORE? (maybe for gene/protein name merging?)
                         select_genes = NULL,
                         fill_all_assays = FALSE,
                         Zchange = TRUE, 
                         select_samples = NULL,
                         patient_column = 'Patient',
                         sample_column = 'Sample',
                         date_column = 'Date',
                         select_assays = c("RNA", 
                                           "Viper", 
                                           "Protein", 
                                           "Phospho")) {
  
  # Select columns from meta table and set names
  meta <- meta %>% 
    select(.data[[patient_column]], .data[[sample_column]], .data[[date_column]]) %>%
    setNames(c('Patient', 'Sample', 'Date'))
  
  # Select samples (might not be present/paired for every assay)
  if (is.null(select_samples)) {
    
    select_samples <- meta %>%
      pull(Sample) %>%
      as.character()
    
  } else {
    
    select_samples <- meta %>%
      filter(Sample %in% select_samples) %>%
      pull(Sample) %>%
      as.character()
    
  }
  
  # Subset meta to select samples
  meta <- meta %>% 
    filter(Sample %in% select_samples)
  
  # Subset assay data
  # TODO: Handle NULL assays
  if (is.null(select_genes)) {
    
    if (is.null(select_gene_cats)) {
      
      select_genes <- rownames(df.rna)
      
    } else {
      
      select_genes <- select_gene_cats$Gene
      
    }
    
  }
  
  # Format Protein to RNA mapping table
  if (!is.null(protein_rna_table)) {
    
    protein_rna_table <- protein_rna_table %>%
      setNames(c('Protein', 'Gene', 'Category')) %>%
      filter(Gene %in% select_genes,
             Protein %in% rownames(df.rppa))
    
  }
  
  # Select proteins and phosphoproteins
  # TODO: Handle NULL assays
  if (!is.null(select_genes)) {
    
    if (!is.null(protein_rna_table)) {
      
      select_proteins <- protein_rna_table %>%
        filter(Gene %in% select_genes) %>%
        pull(Protein)
      
    } else {
      
      select_proteins <- rownames(df.rppa)
      
    } 
    
  } else {
    
    select_proteins <- rownames(df.rppa)
    
  }
  
  # Use Zscore or compute delta between paired sample
  if (Zchange) {
    
    value.var <- 'Zchange'
    
  } else {
    
    value.var <- 'Zscore'
    
  }
  
  # Store assay tables here
  df.list <- list()
  
  # Format RNA expression
  if (!is.null(df.rna)) {
    
    df.list[['RNA']] <- add_assay(df = df.rna, 
                                  assay = 'RNA', 
                                  meta = meta, 
                                  feature_type = 'Gene', 
                                  Zchange = Zchange, 
                                  select_samples = select_samples, 
                                  select_features = select_genes)
    
  }
  
  # Format Viper
  if (!is.null(df.viper)) {
    
    # Scale viper activity
    if (scale_viper) {
      
      df.viper <- t(scale(t(df.viper),scale = TRUE, center = TRUE))
      
    }
    
    df.list[['Viper']] <- add_assay(df = df.viper, 
                                    assay = 'Viper', 
                                    meta = meta, 
                                    feature_type = 'Gene', 
                                    Zchange = Zchange, 
                                    select_samples = select_samples, 
                                    select_features = select_genes)
    
  }
  
  # Scale RPPA and split protein and phosphoprotein
  if (!is.null(df.rppa)) {
    
    # Scale RPPA data
    if (scale_rppa) {
      
      df.rppa <- t(scale(t(df.rppa),scale = TRUE, center = TRUE)) 
      
    }
    
    
    # Split proteins and phosphoproteins into separate matrices
    df.protein <- df.rppa[grep('_p', rownames(df.rppa), invert = TRUE),,drop = FALSE]
    df.phospho <- df.rppa[grep('_p', rownames(df.rppa), invert = FALSE),,drop = FALSE]
    
    if (nrow(df.protein) > 0) {
      
      df.list[['Protein']] <- add_assay(df = df.protein, 
                                        assay = 'Protein', 
                                        meta = meta, 
                                        feature_type = 'Protein', 
                                        Zchange = Zchange, 
                                        select_samples = select_samples, 
                                        select_features = select_proteins)
      
      # Store remaining proteins (incase of missing antibodies in paired samples)
      valid_proteins <- df.list[['Protein']] %>%
        pull(Protein) %>%
        as.character() %>%
        unique() %>%
        sort()
      
    }
    
    # Now phosphoproteins
    if (nrow(df.phospho) > 0) {
      
      df.list[['Phospho']] <- add_assay(df = df.phospho, 
                                        assay = 'Phospho', 
                                        meta = meta, 
                                        feature_type = 'Phospho', 
                                        Zchange = Zchange, 
                                        select_samples = select_samples, 
                                        select_features = select_proteins) 
      
      # Store remaining proteins (incase of missing antibodies in paired samples)
      valid_phosphos <- df.list[['Phospho']] %>%
        pull(Phospho) %>%
        as.character() %>%
        unique() %>%
        sort()
      
    }
    
  }
  

 
  
  
  ##########################################
  
  # Refilter protein_rna_table
  protein_rna_table <- protein_rna_table %>% 
    filter(Gene %in% select_gene_cats$Gene) %>%
    filter(Protein %in% c(valid_proteins, valid_phosphos))
  
  # Create new protein to rna table that includes merged names for final plotting
  merged_names <- merge_rna_protein_table(protein_rna_table) %>% 
    arrange(Protein, Gene)
  
  # Parse each assay and swap gene/protein name with merged full name
  for (n in names(df.list)) {
    
    if (nrow(df.list[[n]])) {
      
      df.list[[n]] <- df.list[[n]] %>% 
        left_join(merged_names) %>%
        mutate(FullName = case_when(is.na(FullName) ~ Gene,
                                    TRUE ~ FullName)) %>%
        mutate(Gene = FullName) %>% 
        select(all_of(c('Sample', 'Gene', value.var, 'Assay'))) 
      
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
  
  
  ##########################################
  
  # Merge dataframe
  df.merged <- do.call(rbind, df.list)
  
  # Set assay levels
  df.merged$Assay <- factor(df.merged$Assay, levels = select_assays)
  
  # Get genes (only genes that are found in at least one sample)
  all_genes <- as.character(intersect(select_gene_cats$Gene, 
                                      df.merged$Gene))
  
  # Update select samples to only those remaining if computing delta
  select_samples <- select_samples[select_samples %in% df.merged$Sample]
    
  # Now create gene by assay matrices for each sample and store in list
  assay_mats <- list()
  for (s in select_samples) {
    
    sub.df <- df.merged %>% 
      filter(Sample == s) 
    
    assay_mats[[s]] <- acast(sub.df, Gene ~ Assay, value.var = value.var)
    
    if (fill_all_assays) {
      
      present <- colnames(assay_mats[[s]])
      no_dat <- select_assays[!select_assays %in% present]
      na_mat <- matrix(NA, nrow = nrow(assay_mats[[s]]), 
                       ncol = length(no_dat), 
                       dimnames = list(rownames(assay_mats[[s]]), no_dat))
      assay_mats[[s]] <- cbind(assay_mats[[s]], na_mat)[,select_assays]
      
    }
    
  }
  
  return(assay_mats)
  
}  




sub_assay_anno <- function(meta, column_groups, anno_spec = NULL,
                           assay_mat_sample_column = 'Sample',
                           border = TRUE, gap = unit(3, "points"),
                           show_annotation_name = FALSE,
                           show_annotation_legend = TRUE,
                           spacer_below = TRUE) {
  
  # If not annotation specs, no annotation
  if (is.null(anno_spec)) {
    
    return(NULL)
    
  }
  
  # Fill with meta values for each instance of each sample in group
  idx <- match(column_groups, meta[,assay_mat_sample_column])
  
  print(idx)
  
  # Build annotation vectors as a named list of columns
  anno_cols <- lapply(names(anno_spec), function(meta_col) {
    as.character(meta[idx, meta_col])
  })
  
  print(anno_cols)
  
  # Rename columns to the display names (anno_name)
  display_names <- vapply(anno_spec, `[[`, character(1), "anno_name")
  names(anno_cols) <- display_names
  
  # Data frame passed to HeatmapAnnotation
  df <- as.data.frame(anno_cols, check.names = FALSE, stringsAsFactors = FALSE)
  
  # Build `col=` mapping: must be named by the *display names* (df colnames)
  col_list <- setNames(
    lapply(anno_spec, `[[`, "anno_colors"),
    display_names
  )
  
  
  # Annotation legend size paramters (can make argument later)
  anno_lgd_params <- list(nrow = 2, 
                          labels_gp = gpar(fontsize = 12), 
                          title_gp = gpar(fontsize = 12, 
                                          fontface = 'bold'))
  
  df <- cbind(df, df)
  #colnames(df) <- c('ER1', )
  print(df)
  
  # Build HeatmapAnnotation object
  if (spacer_below) {
    
    anno <- HeatmapAnnotation(
      df = df,
      #foo = anno_empty(border = FALSE, 
      #                 height = unit(0.00001, "points")),
      col = col_list,
      border = border,
      gap = gap,
      show_annotation_name = show_annotation_name,
      show_legend = show_annotation_legend,
      annotation_legend_param = anno_lgd_params
    )
    
  } else {
    
    anno <- HeatmapAnnotation(
      #foo = anno_empty(border = FALSE, 
      #                 height = unit(0.00001, "points")),
      df = df,
      col = col_list,
      border = border,
      gap = gap,
      show_annotation_name = show_annotation_name,
      show_legend = show_annotation_legend,
      annotation_legend_param = anno_lgd_params
    )
    
  }
    
  
  return(anno)
  
}




sub_assay_heatmap <- function(assay_mats, 
                              meta, 
                              ht_group,
                              
                              group_heatmaps_by = 'Sample', 
                              assay_mat_sample_column = 'Sample',
                              
                              show_annotation_name = FALSE,
                              show_annotation_legend = TRUE,
                              patient_plots = FALSE,
                              pre = NULL,
                              fn = NULL,
                              
                              
                              # Update/remove these after fixing top/bottom annotation code blocks
                              make_anno = 'Treatment',
                              include_ERi = TRUE, 
                              include_pam50 = TRUE,
                              annotate_assay_types = FALSE,
                              
                              # Legend parameters
                              lgdFntSize = 12, 
                              ht_lgd_length = 2,
                              value.var = 'Zscore',
                              
                              # Heatmap parameters 
                              cluster_rows = FALSE, 
                              col_fun = NULL,
                              row_title_rot = 0,
                              
                              # Heatmap dimensions
                              ht_height = unit(17, 'in'), 
                              ht_width = unit(3, 'in'),
                              
                              
                              # Not sure if all of these will be needed
                              group_order = TRUE,
                              order_rows = NULL,
                              row_split = NULL,
                              row_gap = NULL
                          
                              
                              ) {
  
  
  # Get available samples/assays for patient
  ht_samples <- meta %>%
    filter(.data[[group_heatmaps_by]] == ht_group) %>%
    filter(.data[[assay_mat_sample_column]] %in% names(assay_mats)) %>%
    pull(.data[[assay_mat_sample_column]]) %>%
    as.character()
  
  
  # Select all assay matrices for current group and merge
  group_mats <- do.call(cbind, assay_mats[ht_samples])
  
  # TODO: BELOW COULD BE MUCH MROE GENERALIZED
  #  HOW DOES column_groups WORK IF MULTIPLE SAMPLES EXIST FOR SAME GROUP
  
  # Label groupings for heatmaps
  column_groups <- c()
  for (s in ht_samples) { column_groups <- c(column_groups, rep(s, ncol(assay_mats[[s]]))) }
  
  
  print(column_groups)
  
  anno_spec <- list(
    #CDKi = list(anno_name = "CDK4/6i", anno_colors = colors.treatment),
    ERi  = list(anno_name = "ERi",     anno_colors = colors.treatment)
    # pam50 = list(anno_name = "PAM50", anno_colors = colors.pam)
  )
  
  anno_spec <- list(
    #CDKi = list(anno_name = "CDK4/6i", anno_colors = colors.treatment),
    ERi  = list(anno_name = "ERi",     anno_colors = all_extra_drug_cols)
    # pam50 = list(anno_name = "PAM50", anno_colors = colors.pam)
  )
  
  top_anno <- sub_assay_anno(meta = meta,
                             column_groups = column_groups,
                             anno_spec = anno_spec)
  
  btm_anno <- sub_assay_anno(meta = meta,
                             column_groups = column_groups,
                             anno_spec = anno_spec)
  
  
  # Biopsies for each column group
  # WHY IS THIS NEEDED AND BEING USED IN HEATMAP CALL?
  cols.biop <- as.character(meta[match(column_groups, meta$Sample),'Biopsy'])
  
  print(cols.biop)
  
  use_old_annos <- FALSE
  if (use_old_annos) {
  
  
  # Biopsies for each column group
  cols.biop <- as.character(meta[match(column_groups, meta$Sample),'Biopsy'])
  
  # HtanIDs for each column group
  cols.hta <- as.character(meta[match(column_groups, meta$Sample),'HtanID'])
  
  # Treatment for each column group
  cols.treat <- as.character(meta[match(column_groups, meta$Sample),'Treatment'])
  
  # ER therapy
  cols.er <- as.character(meta[match(column_groups, meta$Sample),'ERi'])
  
  # PAM50 (change)
  cols.pam <- as.character(meta[match(column_groups, meta$Sample),'pamChange'])
  
  # Patient response for each column group
  cols.intrinsic <- as.character(meta[match(column_groups, meta$Sample),'PatientResponse'])
  
  # onProgression for each column group
  cols.prog <- as.character(meta[match(column_groups, meta$Sample),'ProgressionStage'])
  
  # Patient
  cols.pat <- as.character(meta[match(column_groups, meta$Sample),'HTAN'])
  
  # Treatment legend that also indicates on-progression biopsies
  onProg.idx <- rep(NA, length(cols.prog))
  onProg.idx[which(cols.prog == 'OnProgression')] <- 8
  
  
  print(cols.treat)
  print(cols.er)
  cat('\n')
  cat('\n')
  
  
  
  
  
  
  
  
  
  # Annotation legend size paramters
  anno_lgd_params <- list(nrow = 2, 
                          labels_gp = gpar(fontsize = 12), 
                          title_gp = gpar(fontsize = 12, 
                                          fontface = 'bold'))
  
  if (make_anno == 'Treatment') {
    
    all_drug_cols <- all_drug_cols
    names(all_drug_cols) <- str_to_title(names(all_drug_cols))
    all_extra_drug_cols <- all_extra_drug_cols
    names(all_extra_drug_cols) <- str_to_title(names(all_extra_drug_cols))
    
    cols.treat <- str_to_title(cols.treat)
    cols.er <- str_to_title(cols.er)
    
    if (include_ERi) {
      
      top_anno <- HeatmapAnnotation('CDK4/6i' = cols.treat,
                                    'ERi' = cols.er,
                                    foo = anno_empty(border = FALSE, 
                                                     height = unit(0.00001, 'points')),
                                    col = list('CDK4/6i' = all_drug_cols,
                                               'ERi' = all_extra_drug_cols), 
                                    border = TRUE, 
                                    show_annotation_name = show_annotation_name,
                                    gap = unit(3, 'points'),
                                    show_legend = show_annotation_legend,
                                    annotation_legend_param = anno_lgd_params)
      
    } else {
      
      top_anno <- HeatmapAnnotation('CDK4/6i' = cols.treat,
                                    foo = anno_empty(border = FALSE, 
                                                     height = unit(0.00001, 'points')),
                                    col = list('CDK4/6i' = all_drug_cols), 
                                    border = TRUE, 
                                    show_annotation_name = show_annotation_name,
                                    gap = unit(3, 'points'),
                                    show_legend = show_annotation_legend,
                                    annotation_legend_param = anno_lgd_params)
      
    }
    
  } else if (make_anno == 'PatientResponse') {
    
    top_anno <- HeatmapAnnotation(PatientResponse = cols.intrinsic,
                                  col = list(PatientResponse = intrinsic_cols), 
                                  show_legend = show_annotation_legend,
                                  border = TRUE, 
                                  show_annotation_name = show_annotation_name,
                                  annotation_legend_param = anno_lgd_params)
    
  } else {
    
    top_anno <- FALSE
    
  }
  
  # Color annotate assay types
  if (annotate_assay_types) {
    
    btm_anno <- HeatmapAnnotation(' ' = anno_empty(border = FALSE, 
                                                   height = unit(0.00001, 'points')),
                                  'PAM50' = cols.pam,
                                  Assay = colnames(group_mats), 
                                  col = list(Assay = assay_cols,
                                             'PAM50' = pam_change_colors),
                                  show_annotation_name = show_annotation_name, 
                                  border = TRUE,
                                  show_legend = show_annotation_legend,
                                  gap = unit(3, 'points'),
                                  annotation_legend_param = anno_lgd_params)
    
  } else {
    
    if (include_pam50) {
      
      btm_anno <- HeatmapAnnotation(' ' = anno_empty(border = FALSE, 
                                                     height = unit(0.00001, 'points')),
                                    'PAM50' = cols.pam,
                                    col = list(Assay = assay_cols,
                                               'PAM50' = pam_change_colors),
                                    show_annotation_name = show_annotation_name, 
                                    border = TRUE,
                                    show_legend = show_annotation_legend,
                                    gap = unit(3, 'points'),
                                    annotation_legend_param = anno_lgd_params)
      
    } else {
      
      btm_anno <- NULL
      
    }
    
  }
  
  }
  
  
  print(head(group_mats))
  
  # Order genes by mean across assays
  if (group_order) {
    
    group_mats <- group_mats[order_rows,]
    
  } else {
    
    group_mats <- group_mats[order(rowMeans(group_mats,na.rm = TRUE)),]
    if (!is.null(row_split)) {
      
      row_split <- gene_cats[rownames(group_mats),'Category']
      
    }
    
  }
  
  
  
  # Set heatmap color function
  if (is.null(col_fun)) {
    
    # Set color intensities
    maxHT <- max(c(2, abs(quantile(group_mats, 0.01, na.rm = TRUE)), 
                   abs(quantile(group_mats, 0.99, na.rm = TRUE))))
    col_fun <- colorRamp2(c(-2,0,2), c("#015D5F",'white','magenta4'))
    ht_cols = c('cyan', 'black', 'magenta')
    col_fun <- colorRamp2(c(-2,0,2), ht_cols)
    
  }
  
  # Build current heatmap object
  ht <- Heatmap(group_mats, 
                name = value.var, 
                cluster_columns = FALSE, 
                column_split = cols.biop, 
                column_gap = unit(0, "mm"),
                column_title = ht_group, 
                cluster_rows = cluster_rows, 
                row_title_rot = row_title_rot, 
                row_split = row_split, 
                row_gap = row_gap,
                width = ht_width, 
                height = ht_height, 
                border = TRUE, 
                top_annotation = top_anno, 
                bottom_annotation = btm_anno,
                heatmap_legend_param = list(direction = 'horizontal', 
                                            title_gp = gpar(fontsize = lgdFntSize, 
                                                            fontface = "bold"), 
                                            labels_gp = gpar(fontsize = lgdFntSize), 
                                            legend_width = unit(ht_lgd_length, 'in'), 
                                            title_position = 'topcenter'),
                col = col_fun, 
                column_names_rot = 45
  )
  
  
  try(dev.off(), silent = TRUE)
  
  # Save individual patient plots to file
  if (patient_plots) {
    
    # First create faux plot to capture dimensions of all heatmap objects
    pdf(NULL)
    dht <- draw(ht) # Length needs to match number of components (total annotations plus heatmap) unit(c(0.2,0,0.2, 0.0001), "mm")
    try(dev.off(), silent = TRUE)
    
    # Measure object dimensions and determine proper figure size
    wh <- calc_ht_size(dht) # wh[1] = width, wh[2] = height
    
    # Now save to file
    fn.ht <- paste0(pre, '/assay_heatmaps/', ht_group, '.', fn, '.heatmap.png')
    png(filename =  fn.ht, width = wh[1], height = (wh[2]+wh[2]*.03), units = 'in', res = wh[1]*wh[2]*2)
    draw(ht, heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom', merge_legend = TRUE)
    try(dev.off(), silent = TRUE)
    
  }
  
  return(ht)
  
  
}


# Function for creating heatmaps of genome assays by patient
multi_assay_heatmap <- function(assay_mats,
                                meta, 
                                
                                # Sample meta parameters
                                assay_mat_sample_column = 'Sample',
                                select_samples = NULL,
                                group_heatmaps_by = 'Sample',
                                select_patients = NULL,  # Replace with select_groups
                                select_groups = NULL, 
                                
                                patientID = 'HTAN', 
                                patient_plots = FALSE, 
                                
                                # Features and feature tables
                                order_rows = NULL,
                                group_order = TRUE,
                                gene_cats = NULL, 
                                protein_rna_tbl = NULL, 
                                category_table = NULL,
                                sub_sep = c(' ', '/', '_'),
                                
                                # Heatmap annotations
                                show_annotation_legend = TRUE,
                                annotate_assay_types = TRUE, 
                                
                                # Remove/update these after fixing annotation code blocks
                                make_anno = 'Treatment',
                                include_ERi = TRUE, 
                                include_pam50 = TRUE,
                                
                                # Legend parameters
                                lgdFntSize = 12, 
                                ht_lgd_length = 2,
                                value.var = 'Zscore',
                                
                                # Heatmap parameters 
                                pre = "/Users/eggerj/Documents/CompBio/HRplus_Project/smmart_hrplus_figures/htan",
                                fn = NULL,
                                cluster_rows = FALSE, 
                                col_fun = NULL,
                                row_title_rot = 0,
                                
                                # Heatmap dimensions
                                add_height = 0, 
                                ht_height = NULL, 
                                add_width = 0,
                                ht_width = unit(3, 'in'),
                                res = NULL) {
  
  # Select samples (intersect of assay_mat names and Sample meta column)
  if (is.null(select_samples)) {
    
    select_samples <- meta %>%
      filter(.data[[assay_mat_sample_column]] %in% names(assay_mats)) %>%
      pull(.data[[assay_mat_sample_column]]) %>%
      as.character()
    
  } else {
    
    select_samples <- meta %>%
      filter(.data[[assay_mat_sample_column]] %in% names(assay_mats)) %>%
      filter(.data[[assay_mat_sample_column]] %in% select_samples) %>%
      pull(.data[[assay_mat_sample_column]]) %>%
      as.character()
    
  }
  
  # Select groups to plot (default: "Sample")
  if (is.null(select_groups)) {
    
    select_groups <- meta %>%
      filter(.data[[assay_mat_sample_column]] %in% select_samples) %>%
      pull(.data[[group_heatmaps_by]]) %>%
      unique() %>%
      sort() %>% # sorting by levels (change?)
      as.character()
    
  } else {
    
    select_groups <- meta %>%
      filter(.data[[assay_mat_sample_column]] %in% select_samples) %>%
      filter(.data[[group_heatmaps_by]] %in% select_groups) %>%
      pull(.data[[group_heatmaps_by]]) %>%
      unique() %>%
      sort() %>% # sorting by levels (change?)
      as.character()
    
  }
  
  # Subset meta to available samples/groups
  meta <- meta %>% 
    filter(.data[[assay_mat_sample_column]] %in% select_samples) %>%
    filter(.data[[group_heatmaps_by]] %in% select_groups)
  
  
  
  # TODO: Update ht_height and ht_width to use cell sizes for different grouping options
  
  # Set heatmap height if not specified
  if (is.null(ht_height)) {
    
    if (nrow(assay_mats[[1]]) >= 100) {
      ht_height <- unit(20, 'in')
    } else if (nrow(assay_mats[[1]]) >= 50) {
      ht_height <- unit(17, 'in')
    } else {
      ht_height <- unit(8, 'in')
    }
    
  }
  
  # Cluster rows
  if (is.logical(cluster_rows)) {
    if (cluster_rows == TRUE) {
      cluster_rows <- create_dendrogram(do.call(cbind, assay_mats))
    }
  }
  
  # Get all genes/proteins across assay matrices
  all_features <- Reduce(union, lapply(assay_mats, FUN = row.names)) %>% 
    unique() %>% 
    sort()
  
  # Sort alphabetically
  if (is.null(order_rows)) { 
    
    order_rows <- sort(all_features) 
    
  } else {
    
    order_rows <- order_rows[order_rows %in% all_features]
    
  }
  
  # Split rows by gene set categories
  if (!is.null(category_table)) {
    
    # Format gene/protein categories
    for (ss in sub_sep) {
      
      category_table <- category_table %>%
        mutate(Category = gsub(ss, '\n', Category))
      
    }
    
    rownames(category_table) <- category_table$Gene
    row_split <- category_table[order_rows,'Category']
    row_gap <- unit(2, 'mm')
    
  } else if (!is.null(gene_cats)) {
    
    # Format gene/protein categories
    for (ss in sub_sep) {
      
      gene_cats <- gene_cats %>%
        mutate(Category = gsub(ss, '\n', Category))
      
    }
    
    # Replace gene names with merged gene/protein/phospho names in gene cats table
    if (!is.null(protein_rna_tbl)) {
      
      gene_cats <- make_merged_rna_protein_categories(assay_mats, 
                                                      gene_cats, 
                                                      protein_rna_tbl)
      
    } 
    
    rownames(gene_cats) <- gene_cats$Gene
    row_split <- gene_cats[order_rows,'Category']
    row_gap <- unit(2, 'mm')
    
  } else {
    
    row_split = NULL
    row_gap = NULL
    
  }
  
  
  
  # Build heatmaps from left to right
  merge_ht <- NULL
  
  for (i in 1:length(select_groups)) {
    
    #p <- select_groups[i]
    ht_group <- select_groups[i]
    
    # Place top/bottom annotation names on right most heatmap
    if (i == length(select_groups)) {
      
      show_annotation_name = TRUE
      
    } else {
      
      show_annotation_name = FALSE
      
    }
    
    
    
    use_old <- FALSE
    if (use_old) {
    
    # Get available samples/assays for patient
    select_samples <- intersect(meta[meta[,group_heatmaps_by] == p,'Sample'],names(assay_mats))
    
    
    
    # Merge matrices
    assay_mat.patient <- do.call(cbind, assay_mats[select_samples])
    
    
    
    
    # Label groupings for heatmaps
    column_groups <- c()
    for (s in select_samples) { column_groups <- c(column_groups, rep(s, ncol(assay_mats[[s]]))) }
    
    
    
    # Biopsies for each column group
    cols.biop <- as.character(meta[match(column_groups, meta$Sample),'Biopsy'])
    
    # HtanIDs for each column group
    cols.hta <- as.character(meta[match(column_groups, meta$Sample),'HtanID'])
    
    # Treatment for each column group
    cols.treat <- as.character(meta[match(column_groups, meta$Sample),'Treatment'])
    
    # ER therapy
    cols.er <- as.character(meta[match(column_groups, meta$Sample),'ERi'])
    
    # PAM50 (change)
    cols.pam <- as.character(meta[match(column_groups, meta$Sample),'pamChange'])
    
    # Patient response for each column group
    cols.intrinsic <- as.character(meta[match(column_groups, meta$Sample),'PatientResponse'])
    
    # onProgression for each column group
    cols.prog <- as.character(meta[match(column_groups, meta$Sample),'ProgressionStage'])
    
    # Patient
    cols.pat <- as.character(meta[match(column_groups, meta$Sample),'HTAN'])
    
    # Treatment legend that also indicates on-progression biopsies
    onProg.idx <- rep(NA, length(cols.prog))
    onProg.idx[which(cols.prog == 'OnProgression')] <- 8
    
    
    # Annotation legend size paramters
    anno_lgd_params <- list(nrow = 2, 
                            labels_gp = gpar(fontsize = 12), 
                            title_gp = gpar(fontsize = 12, 
                                            fontface = 'bold'))
    
    if (make_anno == 'Treatment') {
      
      all_drug_cols <- all_drug_cols
      names(all_drug_cols) <- str_to_title(names(all_drug_cols))
      all_extra_drug_cols <- all_extra_drug_cols
      names(all_extra_drug_cols) <- str_to_title(names(all_extra_drug_cols))
      
      cols.treat <- str_to_title(cols.treat)
      cols.er <- str_to_title(cols.er)
      
      if (include_ERi) {
        
        top_anno <- HeatmapAnnotation('CDK4/6i' = cols.treat,
                                      'ERi' = cols.er,
                                      foo = anno_empty(border = FALSE, height = unit(0.00001, 'points')),
                                      col = list('CDK4/6i' = all_drug_cols,
                                                 'ERi' = all_extra_drug_cols), 
                                      border = TRUE, show_annotation_name = show_annotation_name,
                                      gap = unit(3, 'points'),
                                      show_legend = show_annotation_legend,
                                      annotation_legend_param = anno_lgd_params)
        
      } else {
        
        top_anno <- HeatmapAnnotation('CDK4/6i' = cols.treat,
                                      foo = anno_empty(border = FALSE, height = unit(0.00001, 'points')),
                                      col = list('CDK4/6i' = all_drug_cols), 
                                      border = TRUE, show_annotation_name = show_annotation_name,
                                      gap = unit(3, 'points'),
                                      show_legend = show_annotation_legend,
                                      annotation_legend_param = anno_lgd_params)
        
      }
      
    } else if (make_anno == 'PatientResponse') {
      
      top_anno <- HeatmapAnnotation(PatientResponse = cols.intrinsic,
                                    col = list(PatientResponse = intrinsic_cols), 
                                    show_legend = show_annotation_legend,
                                    border = TRUE, show_annotation_name = show_annotation_name,
                                    annotation_legend_param = anno_lgd_params)
      
    } else {
      
      top_anno <- FALSE
      
    }
    
    # Color annotate assay types
    if (annotate_assay_types) {
      
      btm_anno <- HeatmapAnnotation(' ' = anno_empty(border = FALSE, height = unit(0.00001, 'points')),
                                    'PAM50' = cols.pam,
                                    Assay = colnames(assay_mat.patient), 
                                    col = list(Assay = assay_cols,
                                               'PAM50' = pam_change_colors),
                                    show_annotation_name = show_annotation_name, border = TRUE,
                                    show_legend = show_annotation_legend,
                                    gap = unit(3, 'points'),
                                    annotation_legend_param = anno_lgd_params)
      
    } else {
      
      if (include_pam50) {
        
        btm_anno <- HeatmapAnnotation(' ' = anno_empty(border = FALSE, height = unit(0.00001, 'points')),
                                      'PAM50' = cols.pam,
                                      col = list(Assay = assay_cols,
                                                 'PAM50' = pam_change_colors),
                                      show_annotation_name = show_annotation_name, border = TRUE,
                                      show_legend = show_annotation_legend,
                                      gap = unit(3, 'points'),
                                      annotation_legend_param = anno_lgd_params)
        
      } else {
        
        btm_anno <- NULL
        
      }
      
    }
    
    
    # Order genes by mean across assays
    if (group_order) {
      
      assay_mat.patient <- assay_mat.patient[order_rows,]
      
    } else {
      
      assay_mat.patient <- assay_mat.patient[order(rowMeans(assay_mat.patient,na.rm = TRUE)),]
      if (!is.null(row_split)) {
        
        row_split <- gene_cats[rownames(assay_mat.patient),'Category']
        
      }
      
    }
    
    
    
    # Set heatmap color function
    if (is.null(col_fun)) {
      
      # Set color intensities
      maxHT <- max(c(2, abs(quantile(assay_mat.patient, 0.01, na.rm = TRUE)), 
                     abs(quantile(assay_mat.patient, 0.99, na.rm = TRUE))))
      col_fun <- colorRamp2(c(-2,0,2), c("#015D5F",'white','magenta4'))
      ht_cols = c('cyan', 'black', 'magenta')
      col_fun <- colorRamp2(c(-2,0,2), ht_cols)
      
    }
    
    # Build current heatmap object
    ht <- Heatmap(assay_mat.patient, 
                  name = value.var, 
                  cluster_columns = FALSE, 
                  column_split = cols.biop, 
                  column_gap = unit(0, "mm"),
                  column_title = p, 
                  cluster_rows = cluster_rows, 
                  row_title_rot = row_title_rot, 
                  row_split = row_split, 
                  row_gap = row_gap,
                  width = ht_width, 
                  height = ht_height, 
                  border = TRUE, 
                  top_annotation = top_anno, 
                  bottom_annotation = btm_anno,
                  heatmap_legend_param = list(direction = 'horizontal', 
                                              title_gp = gpar(fontsize = lgdFntSize, 
                                                              fontface = "bold"), 
                                              labels_gp = gpar(fontsize = lgdFntSize), 
                                              legend_width = unit(ht_lgd_length, 'in'), 
                                              title_position = 'topcenter'),
                  col = col_fun, 
                  column_names_rot = 45
    )
    
    
    try(dev.off(), silent = TRUE)
    
    # Save individual patient plots to file
    if (patient_plots) {
      
      # First create faux plot to capture dimensions of all heatmap objects
      pdf(NULL)
      dht <- draw(ht) # Length needs to match number of components (total annotations plus heatmap) unit(c(0.2,0,0.2, 0.0001), "mm")
      try(dev.off(), silent = TRUE)
      
      # Measure object dimensions and determine proper figure size
      wh <- calc_ht_size(dht) # wh[1] = width, wh[2] = height
      
      # Now save to file
      fn.ht <- paste0(pre, '/assay_heatmaps/', p, '.', fn, '.heatmap.png')
      png(filename =  fn.ht, width = wh[1], height = (wh[2]+wh[2]*.03), units = 'in', res = wh[1]*wh[2]*2)
      draw(ht, heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom', merge_legend = TRUE)
      try(dev.off(), silent = TRUE)
      
    }
    
    }
    
    
    # Merge heatmaps  
    #merge_ht <- merge_ht + ht
    
    merge_ht <- merge_ht + sub_assay_heatmap(assay_mats = assay_mats, 
                                             meta = meta, 
                                             ht_group = ht_group,
                                             #p = p, 
                                             group_heatmaps_by = group_heatmaps_by, 
                                             assay_mat_sample_column = assay_mat_sample_column,
                                             
                                             show_annotation_name = show_annotation_name,
                                             show_annotation_legend = show_annotation_legend,
                                             patient_plots = patient_plots,
                                             pre = pre,
                                             fn = fn,
                                             
                                             
                                             # Update/remove these after fixing top/bottom annotation code blocks
                                             make_anno = make_anno,
                                             include_ERi = include_ERi, 
                                             include_pam50 = include_pam50,
                                             annotate_assay_types = annotate_assay_types,
                                             
                                             # Legend parameters
                                             lgdFntSize = lgdFntSize, 
                                             ht_lgd_length = ht_lgd_length,
                                             value.var = value.var,
                                             
                                             # Heatmap parameters 
                                             cluster_rows = cluster_rows, 
                                             col_fun = col_fun,
                                             row_title_rot = row_title_rot,
                                             
                                             # Heatmap dimensions
                                             ht_height = ht_height, 
                                             ht_width = ht_width,
                                             
                                             
                                             # Not sure if all of these will be needed
                                             group_order = group_order,
                                             order_rows = order_rows,
                                             row_split = row_split,
                                             row_gap = row_gap
                                             
                                             
                                             
                                             )
    
  }
  
  try(dev.off(), silent = TRUE)
  
  # First create faux plot to capture dimensions of all heatmap objects
  pdf(NULL)
  dht <- draw(merge_ht, heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom', merge_legend = TRUE) 
  try(dev.off(), silent = TRUE)
  
  # Measure object dimensions and determine proper figure size
  wh <- calc_ht_size(dht) # wh[1] = width, wh[2] = height
  w <- wh[1]
  h <- wh[2]
  if (is.null(res)) {res <- wh[1]*wh[2]*2}
  
  # Add extra width/height if heatmap is being cropped
  w <- w + add_width
  h <- h + add_height
  
  # Now save to file
  #png(filename =  paste0(pre, '/assay_heatmaps/', fn, '.heatmap.png'), width = w, height = h, units = 'in', res = res)
  
  ht_filename <- paste0(pre, '/',  fn, '.heatmap.png')
  print(ht_filename)
  
  png(filename = ht_filename, width = w, height = h, units = 'in', res = res)
  
  
  draw(merge_ht, heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom', merge_legend = TRUE)
  try(dev.off(), silent = TRUE)
  
  return(merge_ht)
  
}

