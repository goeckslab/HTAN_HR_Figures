

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



# Function to create new merged gene/protein names when integrating 
#  genes, proteins, and phosphoproteins across assays
merge_rna_protein_names <- function(protein_rna_tbl, additional_RNA = NULL) {
  
  # Add additional genes to table to include in multi assay heatmap
  if (!is.null(additional_RNA)) {
    
    protein_rna_tbl <- tibble(RNA = additional_RNA) %>%
      full_join(protein_rna_tbl)
    
  }
  
  # Update name if "Gene" is used instead of "RNA"
  col.names <- colnames(protein_rna_tbl)
  if ("Gene" %in% col.names) {
    
    protein_rna_tbl <- protein_rna_tbl %>%
      mutate(RNA = Gene) 
      
  }
  
  # Set missing RNA to NA
  protein_rna_tbl <- protein_rna_tbl %>% 
    mutate(RNA = case_when(RNA == '' ~ NA_character_,
                           TRUE ~ RNA)) %>%
    filter(!is.na(RNA)) %>%
    select(RNA, Protein) %>%
    distinct()
  
  # Split Protein and Phosphoproteins into separate columns
  protein_only <- protein_rna_tbl %>%
    filter(!grepl('_p', Protein))
  
  phospho_only <- protein_rna_tbl %>%
    filter(grepl('_p', Protein)) %>%
    setNames(c('RNA', 'Phospho'))
  
  # Merge back together with proteins and phosphoproteins now as separate columns
  protein_rna_tbl <- protein_only %>%
    full_join(phospho_only) %>%
    distinct() %>% 
    select(RNA, Protein, Phospho) 
  
  # Separate phosphoprotein protein names and actual phospho groups
  protein_rna_tbl <- protein_rna_tbl %>%
    separate(col = 'Phospho', 
             into = c('Phospho', 'Phosphogroup'), 
             sep = '_', 
             extra = 'merge')
  
  # Merge RNA, protein, and phosphoproteins into new name
  protein_rna_tbl <- protein_rna_tbl %>%
    
    mutate(
      
      # Clean up names to find matches
      gene_clean    = gsub("-", "", str_to_lower(RNA)),
      protein_clean = gsub("-", "", str_to_lower(Protein)),
      phospho_clean = gsub("-", "", str_to_lower(Phospho)),
      
      # Create name
      FullName = case_when(
        
        # Gene + protein, no phospho
        is.na(Phospho) & !is.na(Protein) &
          gene_clean == protein_clean ~
          RNA,
        
        is.na(Phospho) & !is.na(Protein) ~
          paste(RNA, Protein, sep = "/"),
        
        # Gene + phospho, no protein
        is.na(Protein) & !is.na(Phospho) &
          gene_clean == phospho_clean ~
          paste(RNA, Phosphogroup, sep = "/"),
        
        is.na(Protein) & !is.na(Phospho) ~
          paste0(RNA, "/", Phospho, "_", Phosphogroup),
        
        # Gene only
        is.na(Protein) & is.na(Phospho) ~
          RNA,
        
        # Gene + protein + phospho
        gene_clean == protein_clean &
          gene_clean == phospho_clean ~
          paste(RNA, Phosphogroup, sep = "/"),
        
        protein_clean == phospho_clean ~
          paste0(RNA, "/", Protein, "_", Phosphogroup),
        
        TRUE ~
          paste0(RNA, "/", Protein, "/", Phospho, "_", Phosphogroup)
      )
    ) %>%
    select(-gene_clean, -protein_clean, -phospho_clean)
  

  # Manually swap certain names
  full_name_map <- c(
    "AURKA/Aurora-A/Aurora-ABC_pT288_pT232_pT198" = "AURKA_pT288_pT232_pT198",
    "BAK1/BAK" = "BAK1",
    "CCNB1/CYCLINB1" = "CCNB1",
    "CCND1/CYCLIND1" = "CCND1",
    "CCND3/Cyclin-D3" = "CCND3",
    "CCNE1/CYCLINE1" = "CCNE1",
    "CHEK1/CHK1_pS296" = "CHK1_pS296",
    "CHEK1/CHK1_pS345" = "CHK1_pS345",
    "CHEK2/CHK2_pT68" = "CHK2_pT68",
    "EIF4EBP1/X4EBP1_pS65" = "EIF4EBP1/4EBP1_pS65",
    "EIF4EBP1/X4EBP1_pT37T46" = "EIF4EBP1/4EBP1_pT37T46",
    "GZMB/Granzyme-B" = "GZMB",
    "PRKCA/PKCALPHA/PKC-a-b-II_pT638_T641" = "PRKCA/PKC-a-b-II_pT638_T641",
    "RPS6KB1/P70S6K1/P70S6K_pT389" = "RPS6KB1/P70S6K1_pT389",
    "STAT5A/STAT5ALPHA" = "STAT5A"
  )
  
  protein_rna_tbl <- protein_rna_tbl %>%
    mutate(
      FullName = recode(FullName, !!!full_name_map)
    )
  
  # Merge phosphoproteins back to original names
  protein_rna_tbl <- protein_rna_tbl %>% 
    mutate(Phospho = case_when(!is.na(Phospho) ~ paste(Phospho, Phosphogroup, sep = '_'),
                               TRUE ~ Phospho)) %>%
    select(RNA, Protein, Phospho, FullName) %>%
    setNames(c('Gene', "Protein", 'Phospho', 'FullName')) %>%
    arrange(Gene, Protein, Phospho, FullName)
  
  return(protein_rna_tbl)
  
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
                         
                         # Merge RNA/Protein
                         merged_rna_protein_table = NULL,
                         select_merged_names = NULL,
                         
                         # If no merged table available, use to make it
                         protein_rna_table = NULL, 
                         
                        
                         # Not functional
                         select_genes = NULL,
                         select_proteins = NULL,
                         
                         
                         full_matches_only = TRUE,
                         
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
  
  
  # ~~~ SAMPLE DATA ~~~ #
  
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
  
  
  # ~~~ FEATURE (RNA/PROTEIN) DATA ~~~ #
  
  # Check merged RNA/Protein table exists
  if (is.null(merged_rna_protein_table)) {
    
    if (is.null(protein_rna_table)) {
      
      if (full_matches_only) {
        
        stop('Unable to merge RNA and Protein names without mapping table')
        
      }
      
    } else {
      
      merged_rna_protein_table <- merge_rna_protein_names(protein_rna_table)
      
    }
    
  }
  
  # Reduce to select genes/RNAs
  if (is.null(select_genes)) {
    
    all_rna <- union(rownames(df.rna), rownames(df.viper)) 
    
    # Gene must be present in merged naming table
    if (full_matches_only) {
      
      select_genes <- all_rna[all_rna %in% merged_rna_protein_table$Gene]
      
    } else {
      
      select_genes <- all_rna
      
    }
    
  }
  
  # Reduce to select proteins/phosphoproteins
  if (is.null(select_proteins)) {
    
    all_protein <- rownames(df.rppa)
    
    # Protein/phosphoprotein must be present in merged naming table
    if (full_matches_only) {
      
      all_matched_protein <- c(merged_rna_protein_table$Protein,
                              merged_rna_protein_table$Phospho)
      
      select_proteins <- all_protein[all_protein %in% all_matched_protein]
    
    } else {
      
      select_proteins <- all_protein
      
    }
    
  }
  
  
  # Reduce merged naming table to only genes/proteins with matches
  if (full_matches_only) {
    
    # Filter to select features based on merged RNA/protein names
    if (!is.null(select_merged_names)) {
      
      merged_rna_protein_table <- merged_rna_protein_table %>%
        filter(FullName %in% select_merged_names)
      
    }
    
    # Make sure features exist in RNA/protein data
    merged_rna_protein_table <- merged_rna_protein_table %>%
      filter(Gene %in% select_genes | 
               Protein %in% select_proteins | 
               Phospho %in% select_proteins)
    
    # Update select RNA/protein
    #  THIS ISN"T REALLY FUNCTIONAL ANYMORE 
    #    NEED TO FIGURE OUT HOW WE WANT TO FILTER GENES TO PROTEINS
    select_genes <- merged_rna_protein_table$Gene
    select_proteins <- c(merged_rna_protein_table$Protein, 
                         merged_rna_protein_table$Phospho)
    
  }
  
  
  
  # TODO: Handle NULL assays
  
  
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
    
    # Format protein data
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
  
  
  # TODO: WHAT TO DO ABOUT DROP OUT IN RPPA DATA?
  
  # Refilter merged naming table incase of RPPA feature dropout
  #merged_rna_protein_table <- merged_rna_protein_table %>%
  #  filter(Gene %in% select_genes) %>%
  #  filter(Protein %in% valid_proteins | Phospho %in% valid_phosphos)
  
    
  # Parse each assay and swap gene/protein name with merged full name
  for (n in names(df.list)) {
    
    if (nrow(df.list[[n]])) {
      
      df.list[[n]] <- df.list[[n]] %>% 
        left_join(merged_rna_protein_table) %>%
        
        # Still needed?
        mutate(FullName = case_when(is.na(FullName) ~ Gene,
                                    TRUE ~ FullName)) %>%
        #mutate(Gene = FullName) %>% 
        #select(all_of(c('Sample', 'Gene', value.var, 'Assay'))) 
        select(all_of(c('Sample', 'FullName', value.var, 'Assay'))) 
      
    }
    
  }
  

  
  ##########################################
  
  # Merge dataframe
  df.merged <- do.call(rbind, df.list)
  
  # Set assay levels
  df.merged$Assay <- factor(df.merged$Assay, levels = select_assays)
  
  # Get genes (only genes that are found in at least one sample)
  all_features <- as.character(intersect(merged_rna_protein_table$FullName, df.merged$FullName))
  
  # Update select samples to only those remaining if computing delta
  select_samples <- select_samples[select_samples %in% df.merged$Sample]
    
  # Now create gene by assay matrices for each sample and store in list
  assay_mats <- list()
  for (s in select_samples) {
    
    sub.df <- df.merged %>% 
      filter(Sample == s) 
    
    assay_mats[[s]] <- acast(sub.df, FullName ~ Assay, value.var = value.var)
    
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

# Function to rename heatmap annotation names if present in anno_specs list
get_display_names <- function(anno_spec) {
  
  # Rename columns to the display names (anno_name); if missing, use the outer name
  display_names <- vapply(names(anno_spec), function(meta_col) {
    
    
    nm <- anno_spec[[meta_col]]$anno_name
    
    # Use if listed
    if (is.null(nm) || nm == "") {
      
      meta_col 
      
    # Else, use original column name  
    } else { 
      
      nm
      
    }  
    
  }, character(1))
  
  return(display_names)
  
}


sub_assay_anno <- function(assay_mats,
                           meta, 
                           ht_samples, 
                           anno_spec = NULL,
                           assay_mat_sample_column = 'Sample',
                           border = TRUE, gap = unit(3, "points"),
                           show_annotation_name = FALSE,
                           show_annotation_legend = TRUE,
                           spacer_below = TRUE,
                           spacer_gap = unit(0.00001, "points")) {
  
  # If no annotation specs, no annotation
  if (is.null(anno_spec) || length(anno_spec) == 0) {
    
    return(NULL)
    
  }
  
  # For each sample in group, mulitply by number of coloumns in assay mat for that sample
  column_samples <- rep(ht_samples, vapply(assay_mats[ht_samples], ncol, integer(1)))
  
  # Fill with meta values for each instance of each sample in group
  idx <- match(column_samples, meta[,assay_mat_sample_column])
  
  # Rename columns to the display names (anno_name); if missing, use the outer name
  display_names <- get_display_names(anno_spec)
 
  
  # Sub meta for sub heatmap annotations
  df <- meta %>%
    slice(idx) %>%
    select(all_of(names(anno_spec))) %>%
    setNames(display_names)
  
  # Make color list for annotations (make sure to match new display names)
  color_list <- anno_spec %>%
    lapply(`[[`, "anno_colors") %>%
    rlang::set_names(display_names)
  
  
  # Annotation legend size paramters (can make argument later)
  anno_lgd_params <- list(
    nrow = 2, 
    labels_gp = gpar(fontsize = 12), 
    title_gp = gpar(fontsize = 12, 
                    fontface = 'bold')
  )
  
  # Build spacer to reduce text cluttering
  spacer <- anno_empty(border = FALSE, 
                       height = spacer_gap)
  
  # Build HeatmapAnnotation object with spacer above or below
  if (spacer_below) {
    
    anno <- HeatmapAnnotation(
      df = df,
      spacer = spacer,
      col = color_list,
      border = border,
      gap = gap,
      show_annotation_name = show_annotation_name,
      show_legend = show_annotation_legend,
      annotation_legend_param = anno_lgd_params
    )
    
  } else {
    
    anno <- HeatmapAnnotation(
      spacer = spacer,
      df = df,
      col = color_list,
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
                              plot_each_group = FALSE,
                              pre = NULL,
                              fn = NULL,
                              
                              # New for annotations (list of specs)
                              top_anno = NULL,
                              btm_anno = NULL,
                              
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
                              add_height = 0,
                              add_width = 0,
                              
                              
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
  
  
  # Select all assay matrices for current group and merge for heatmap
  group_mats <- do.call(cbind, assay_mats[ht_samples])
  
  # Build top annotations using top_anno_specs
  top_anno <- sub_assay_anno(assay_mats = assay_mats, 
                             meta = meta,
                             ht_samples = ht_samples,
                             anno_spec = top_anno,
                             spacer_below = TRUE,
                             show_annotation_name = show_annotation_name,
                             show_annotation_legend = show_annotation_legend)
  
  # Build bottom annotations using btm_anno_specs
  btm_anno <- sub_assay_anno(assay_mats = assay_mats, 
                             meta = meta,
                             ht_samples = ht_samples,
                             anno_spec = btm_anno,
                             spacer_below = FALSE,
                             show_annotation_name = show_annotation_name,
                             show_annotation_legend = show_annotation_legend)
  
  
 
  
  # Order genes by mean across assays
  if (group_order) {
    
    group_mats <- group_mats[order_rows,]
    
  } else {
    
    group_mats <- group_mats[order(rowMeans(group_mats,na.rm = TRUE)),]
    if (!is.null(row_split)) {
      
      row_split <- gene_cats[rownames(group_mats),'Category']
      
    }
    
  }
  
  

  
  # Build current heatmap object
  ht <- Heatmap(group_mats, 
                
                # Column features
                top_annotation = top_anno, 
                bottom_annotation = btm_anno,
                cluster_columns = FALSE, 
                column_title = ht_group, 
                column_names_rot = 45,
                
                # Row features
                cluster_rows = cluster_rows, 
                row_title_rot = row_title_rot, 
                row_split = row_split, 
                row_gap = row_gap,
                
                # Heatmap size
                width = ht_width, 
                height = ht_height, 
                border = TRUE, 
                
                # Legend
                name = value.var,
                col = col_fun, 
                heatmap_legend_param = list(direction = 'horizontal', 
                                            title_gp = gpar(fontsize = lgdFntSize, 
                                                            fontface = "bold"), 
                                            labels_gp = gpar(fontsize = lgdFntSize), 
                                            legend_width = unit(ht_lgd_length, 'in'), 
                                            title_position = 'topcenter')
  )
  
  
  # Save individual patient plots to file
  if (plot_each_group) {
    
    # Make group filename (make sure group string is valid for filenaming)
    fn.ht <- paste0(pre, '/',  make.names(ht_group), '.', fn, '.heatmap.png')
    
    # Save as png
    save_htan_heatmap(ht_objects = list(ht, 
                                        NULL), 
                      fn = fn.ht, 
                      ht_gap = unit(2, "mm"), # default is 2mm?
                      
                      add_height = add_height, 
                      add_width = add_width,
                      
                      extend_w = 0 # 1.5 used for make_heatmap()
                      )
    
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
                                plot_each_group = FALSE, 
                                
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
                                
                                # New for annotations (list of specs)
                                top_anno = NULL,
                                btm_anno = NULL,
                                
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
      as.character()
    
  } else {
    
    select_groups <- meta %>%
      filter(.data[[assay_mat_sample_column]] %in% select_samples) %>%
      filter(.data[[group_heatmaps_by]] %in% select_groups) %>%
      mutate(!!group_heatmaps_by := factor(.data[[group_heatmaps_by]], 
                                           levels = select_groups)) %>%
      pull(.data[[group_heatmaps_by]]) %>%
      unique() %>%
      sort() %>%
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
  
  # Set heatmap color function
  if (is.null(col_fun)) {
    
    col_fun <- make_heatmap_colors(NULL, minHt = -2, maxHt = 2)
   
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
    
    # Get current group
    ht_group <- select_groups[i]
    
    # Place top/bottom annotation names on right most heatmap
    if (i == length(select_groups)) {
      
      show_annotation_name = TRUE
      
    } else {
      
      show_annotation_name = FALSE
      
    }
    
    
    # Add heatmap
    merge_ht <- merge_ht + 
      sub_assay_heatmap(assay_mats = assay_mats, 
                        meta = meta, 
                        
                        ht_group = ht_group,
                        group_heatmaps_by = group_heatmaps_by, 
                        assay_mat_sample_column = assay_mat_sample_column,
                        
                        
                        top_anno = top_anno,
                        btm_anno = btm_anno,
                        show_annotation_name = show_annotation_name,
                        show_annotation_legend = show_annotation_legend,
                        
                        
                        plot_each_group = plot_each_group,
                        pre = pre,
                        fn = fn,
                        
                        
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
                        add_height = add_height,
                        add_width = add_width,
                        
                        # Not sure if all of these will be needed
                        group_order = group_order,
                        order_rows = order_rows,
                        row_split = row_split,
                        row_gap = row_gap
      )
    
  }
  
  # Filename
  fn.ht <- paste0(pre, '/',  fn, '.heatmap.png')
  
  # Save merged heatmaps as png
  save_htan_heatmap(ht_objects = list(merge_ht, NULL), 
                    fn = fn.ht, 
                    
                    ht_gap = unit(2, "mm"), # default is 2mm?
                    
                    add_height = add_height, 
                    add_width = add_width,
                    extend_w = 0 # 1.5 used for make_heatmap()
                    
  )
  
  
  
  return(merge_ht)
  
}

