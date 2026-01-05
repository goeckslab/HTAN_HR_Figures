



# Function to perform a permutation test for enrichment scores
perm_test_module <- function(s1, s2, n_perm = 10000, alternative = 'greater') {
  
  #observed_stat <- mean(s1)
  observed_stat <- median(s1)
  
  perm_stats <- replicate(n_perm, {
    
    combined <- sample(c(s1, s2))  # Permute combined vector
    group1 <- combined[1:length(s1)]
    #mean(group1)
    median(group1)
    
  })
  
  
  if (alternative == 'greater') {
    
    #pval <- mean(perm_stats >= observed_stat)
    pval <- (1 + sum(perm_stats >= observed_stat)) / (n_perm + 1)
    
  } else if (alternative == 'less') {
    
    #pval <- mean(perm_stats <= observed_stat)
    pval <- (1 + sum(perm_stats <= observed_stat)) / (n_perm + 1)
  } 
  else {  # two-sided
    
    #pval <- mean(abs(perm_stats - mean(c(s1, s2))) >= abs(observed_stat - mean(c(s1, s2))))
    center <- mean(c(s1, s2))
    pval <- (1 + sum(abs(perm_stats - center) >= abs(observed_stat - center))) / (n_perm + 1)
    
  }
  
  return(pval)
}


get_sciATAC_module_permpvals <- function(df, 
                                         select_pathways, 
                                         meta = NULL,
                                         
                                         #sample_group = 'HtanID', 
                                         sample_group = 'Sample',
                                         
                                         alt_group = 'rest', 
                                         n_perm = 10000, 
                                         alternative = 'greater') {
  
  # Merge with meta table (needed if testing group column not sciATAC enrichment score table)
  if (!is.null(meta)) {
    
    df <- df %>%
      left_join(meta)
    
  }
  
  # Get all sample IDs for sciATAC
  samples <- df %>% 
    pull(.data[[sample_group]]) %>% 
    as.character() %>% 
    unique()
  
  # Make sure pathways are available
  select_pathways <- select_pathways[select_pathways %in% colnames(df)]
  
  # Dataframe to store p-values
  pvals <- data.frame(Module = character(), 
                      Sample = character(), 
                      Pvalue = numeric())
  
  # Parse each module
  for (s in samples) {
    
    # Parse each pathway (module)
    for (m in select_pathways) {
      
      s1 <- df %>%
        filter(.data[[sample_group]] == s) %>%
        pull(.data[[m]])
      
      # Get all other cells enrichment scores for current pathway 
      s2 <- if (alt_group == 'rest') {
        
        df %>% 
          filter(.data[[sample_group]] != s) %>%
          pull(.data[[m]])
        
        # Or all cells plus current sample (more conservative)    
      } else {
        
        df %>% 
          pull(.data[[m]])
        
      }
      
      # Run test and add pvalues to dataframe
      rez <- perm_test_module(s1, s2, n_perm = n_perm, alternative = alternative)
      #pvals <- rbind(pvals, data.frame(Module = m, Sample = s, Pvalue = rez))
      pvals <- rbind(pvals, data.frame(m, s, rez))
      
    }
  }
  
  # Perform multiple correction
  pvals <- pvals %>% 
    setNames(c('Pathway', sample_group, 'pvalue')) %>% 
    mutate(FDR = p.adjust(pvalue, method = 'fdr')) %>% 
    mutate(log10 = -log10(FDR))
  
  return(pvals)
  
}



# Function to compute statistical testing on module scores from pathway genesets
#  Computes test on sciATAC cells from each sample by testing all cells in that sample to
#   all cells from all other samples (default) OR all single cell samples in cohort
get_sciATAC_pathway_pvals <- function(df, 
                                      select_pathways, 
                                      meta = NULL,
                                      
                                      #sample_group = 'HtanID', 
                                      sample_group = 'Sample',
                                      
                                      alt_group = 'all', 
                                      alternative = 'greater') {
  
  # Merge with meta table (needed if testing group column not sciATAC enrichment score table)
  if (!is.null(meta)) {
    
    df <- df %>%
      left_join(meta)
    
  }
  
  # Get all sample IDs for single cell
  samples <- df %>% 
    pull(.data[[sample_group]]) %>% 
    as.character() %>% 
    unique()
  
  # Make sure select_pathways present in dataframe of module scores
  #   Each row is cell and columns contain meta elements or pathways
  select_pathways <- intersect(select_pathways, colnames(df))
  
  # Start data frame to store pvalues for each Pathway and sample
  pvals <- tibble(Pathway = character(), 
                  Sample = character(), 
                  Pvalue = numeric())
  
  # Parse each sample
  for (s in samples) {
    
    # Parse each Pathway
    for (pw in select_pathways) {
      
      # Get all cells (rows) for current sample
      s1 <- df %>% 
        filter(.data[[sample_group]] == s) %>%
        pull(.data[[pw]])
      
      # Get all cells for all other samples
      if (alt_group == 'rest') {
        
        s2 <- df %>% 
          filter(.data[[sample_group]] != s) %>%
          pull(.data[[pw]])
        
        # Or, get all cells from dataset    
      } else {
        
        s2 <- df %>% 
          pull(.data[[pw]])
        
      }
      
      # Perform test and add to dataframe
      rez <- wilcox.test(s1, s2, alternative = alternative)$p.value
      pvals <- bind_rows(
        pvals, 
        tibble(Pathway = pw, 
               Sample = s, 
               Pvalue = rez)
      )
      
    }
    
  }
  
  # Format p-values dataframe and correct for multiple testing
  pvals <- pvals %>% 
    setNames(c('Pathway', sample_group, 'pvalue')) %>% 
    mutate(FDR = p.adjust(pvalue, method = 'fdr')) %>% 
    mutate(log10 = -log10(FDR))
  
  return(pvals)
  
}

# Legacy calls
get_sciATAC_module_pvals <- get_sciATAC_pathway_pvals


# Function for plotting sciATAC pathway enrichment distributions and statistical testing
plot_sciATAC_pathway_activity <- function(df, 
                                          meta,
                                          select_pathways, 
                                          
                                          x = 'Sample', 
                                          col_by = 'Sample', 
                                          
                                          celltypes = NULL, 
                                          select_samples = NULL, 
                                          sample_order = NULL,
                                          test_all_samples = TRUE,
                                          ncol = NULL, 
                                          category_table = NULL, 
                                          hide_legend = FALSE,
                                          hide_x_axis = FALSE, 
                                          facet_scale = 'fixed',
                                          order_by_sample = FALSE,
                                          show_cell_count = FALSE,
                                          fixed_y_scale = TRUE, 
                                          facet_formula = c('Category ~ Pathway'),
                                          return_plot = TRUE,
                                          fill_collins = NULL,
                                          alt_group = 'all', 
                                          fn = NULL, 
                                          width = NULL, 
                                          height = NULL) {
  
  # Merge with meta table
  df <- df %>%
    left_join(meta)
  
  # Factor level?
  if (is.null(sample_order)) {
    
    sample_order <- df %>%
      pull(.data[[x]]) %>%
      unique()
    
  }
    
  # Set factor levels for ordering
  df <- df %>%
    mutate(!!x := factor(.data[[x]], levels = sample_order))
    
  
  # Get available pathways
  select_pathways <- select_pathways[select_pathways %in% colnames(df)]
  
  # Filter to select cell types (e.g. tumor only)
  if (!is.null(celltypes)) { 
    
    df <- df %>% 
      filter(Celltype %in% celltypes) 
    
  }
  
  # Subset down to select samples before statistical testing
  if (!is.null(select_samples) & !test_all_samples) {
    
    df <- df %>% 
      filter(Sample %in% select_samples)
    
  }
  
  # Merge sample and cluster
  df <- df %>% 
    mutate(SampleCluster = paste(Sample, Cluster, sep = '_'))
  
  # Compute pvalues for each module
  cat('Computing hypothesis testing and gathering p-values...\n')
  module_pvals <- get_sciATAC_pathway_pvals(df, select_pathways = select_pathways, sample_group = x, alt_group = 'rest')
  cat('\nDone\n')
  
  
  # Subset down to select samples after statistical testing
  if (!is.null(select_samples) & test_all_samples) {
    
    df <- df %>% 
      filter(Sample %in% select_samples)
    
  }
  
  
  # Subset down to needed columns and merge with module pvals
  select_cols <- unique(c('Barcode', x, col_by, 'SampleCluster', 'Cluster', 'Celltype'))
  select_cols <- select_cols[select_cols %in% colnames(df)]
  select_cols <- c(select_cols, c('Pathway', 'Score'))
  df <- df %>%
    pivot_longer(cols = select_pathways,
                 names_to = 'Pathway',
                 values_to = 'Score') %>%
    select(all_of(select_cols)) %>%
    left_join(module_pvals) %>%
    group_by(.data[[x]], Pathway) %>%
    mutate(median_score = median(Score),
           Count = n()) 
  
  
  # Add cell counts to legend
  if (show_cell_count) {
    
    df <- df %>%
      mutate(lgd := paste0(.data[[x]], ' (', Count, ')')) %>%
      ungroup()
    
  } else {
    
    df <- df %>%
      mutate(lgd := .data[[x]]) %>%
      ungroup()
    
  }
  
  # Set order of x axis
  if (order_by_sample) {
    
    lvls <- df %>% 
      arrange(.data[[x]]) %>% 
      pull(.data[[x]]) %>%
      unique()
    
  } else {
    
    # Set factor levels based on group median (highest to lowest)
    lvls <- df %>%
      arrange(-median_score) %>%
      pull(.data[[x]]) %>%
      as.character() %>%
      unique()
    
  }
  
  # Update levels in table
  df <- df %>%
    mutate(!!x := factor(.data[[x]], levels = lvls))
  
  
  
  # Format to scientific notation
  expSup <- function(w, digits=0) {
    tryCatch(sprintf(paste0("%.", digits, "f*x*10^%d"), w/10^floor(log10(abs(w))), floor(log10(abs(w)))), error = function(e) {0})
  }
  
  
  # Get dataframe for plotting FDR above each boxplot
  pvals <- df %>%
    select(c(x, 'FDR', Pathway)) %>%
    mutate(Sig = case_when(FDR < 0.01 ~ '*',
                           TRUE ~ '')) %>%
    distinct()
  
  
  # So counts are only in legend
  if (col_by == x) {
    
    col_by <- 'lgd'
    lgd_name <- x
    
  } else {
    
    x <- 'lgd'
    lgd_name <- col_by
    pvals <- pvals %>% 
      left_join(df)
    
  }
  
  
  # Drop outliers before plotting
  df <- df %>%
    group_by(Pathway) %>%
    mutate(q1 = quantile(Score, 0.01)) %>%
    mutate(q99 = quantile(Score, 0.99)) %>%
    filter(Score > q1) %>%
    filter(Score < q99) %>%
    mutate(pwMax = round_any(max(Score), 0.25, f = ceiling)) %>%
    ungroup()
  
  yMax <- round_any(max(df %>% pull(Score)), 0.25, f = ceiling)
  yMin <- round_any(min(df %>% pull(Score)), 0.25, f = floor)
  
  
  if (fixed_y_scale) {
    
    ast_pos <- 'yMax'
    
  } else {
    
    pvals <- df %>%
      select(Pathway, pwMax) %>%
      right_join(pvals)
    
    ast_pos <- 'pwMax'
    
  }
  
  # If legend shows Sample.Drug, change to "Sample" for title
  if (lgd_name == 'Sample.Drug') {lgd_name <- 'Sample'}
  
  # Fill colors
  if (is.null(fill_collins)) {
    
    fill_collins <- htan_sample_cols
    
  }
  
  # Merge categories if provided
  if (!is.null(category_table)) {
    
    df <- df %>% 
      left_join(category_table)
    
    g <- ggplot(df) +
      geom_violin(aes_string(x = x, y = 'Score', fill = col_by), 
                  draw_quantiles = c(0.25, 0.75)) + 
      geom_text(data = pvals, 
                aes_string(x = x, y = ast_pos, label = 'Sig'), 
                position = position_dodge(width = 1), 
                col = 'red',
                size = 10.5,
                vjust = 1) + 
      facet_grid(as.formula(facet_formula), scales = facet_scale) + #, ncol = 4) +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            text = element_text(size = 15)) +
      stat_summary(data = df, 
                   aes_string(x = x, y = 'Score'),
                   fun.y="mean", 
                   geom="point", 
                   size = 2, 
                   alpha = 0.6) +
      scale_fill_manual(name = lgd_name, 
                        values = fill_collins) +
      xlab(NULL) + 
      ylab('Chromatin Accessibility\nEnrichment')
    
  }
  
  else {
    
    g <- ggplot(df) +
      geom_violin(aes_string(x = x, y = 'Score', fill = col_by), 
                  draw_quantiles = c(0.25, 0.75)) + 
      geom_text(data = pvals, 
                aes_string(x = x, y = ast_pos, label = 'Sig'), 
                parse = FALSE,
                position = position_dodge(width = 1), 
                col = 'red',
                size = 10.5,
                vjust = 1) + 
      facet_wrap(~ Pathway, ncol = ncol, scales = facet_scale) +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            text = element_text(size = 15)) +
      stat_summary(data = df, 
                   aes_string(x = x, y = 'Score'),
                   fun.y="mean", 
                   geom="point", 
                   size = 2, 
                   alpha = 0.6) +
      scale_fill_manual(name = lgd_name, 
                        values = fill_collins) +
      xlab(NULL) + 
      ylab('Chromatin Accessibility\nEnrichment')
    
  }
  
  # Fix y scale
  if (fixed_y_scale) {
    
    g <- g + scale_y_continuous(breaks = seq(yMin, yMax, 0.25), 
                                limits = c(yMin, yMax))
    
  }
  
  # Hide legend\
  if (hide_legend) { 
    
    g <- g + theme(legend.position = "none") 
    
  }
  
  # Hide x axis labels
  if (hide_x_axis) {
    
    g <- g + theme(axis.text.x=element_blank(),
                   axis.ticks.x=element_blank())
    
  }
  
  # make 1 row shorter, 3 rows wider
  
  if (!is.null(fn)) {
    
    numPWs <- df %>% 
      pull(Pathway) %>% 
      unique() %>% 
      length()
    
    if (numPWs < 4) {
      
      w <- 14
      h <- 5
      
    } else if (numPWs > 9) {
      
      w <- 18
      h <- 8
      
    } else {
      
      w <- 14
      h <- 8
      
    }
    
    if (is.null(width)) {width <- w}
    if (is.null(height)) {height <- h}
    
    ggsave(g, filename = fn, width = width, height = height)
    
  }
  
  if (return_plot) {
    
    return(g)
    
  }
  
}

