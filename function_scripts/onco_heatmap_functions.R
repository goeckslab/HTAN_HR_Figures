

# --------------------------------------------------------------------------
#
#     FUNCTIONS FOR CREATING ONCOPLOTS USING COMPLEXHEATMAP
# 
# --------------------------------------------------------------------------

library(ComplexHeatmap)
library(rlang)
library(dplyr)

# Set color scheme for variant types
variant_cols <- c('Gain' = "#E31A1C", 
                  'Loss' = "#1F78B4", 
                  "Frameshift Deletion" = "#A6CEE3", 
                  "Frameshift Insertion" = "#FFFF99", 
                  "In Frame Deletion" = "#B2DF8A", 
                  "In Frame Insertion" = "#6A3D9A", 
                  "In Frame Deletion/insertion (indels)" = 'grey15', #change (and probably not keep?)
                  "In Frame Indel" = "#9F4D23FF", #change (and probably not keep?)
                  "Intron Substitution" = "#FB9A99", 
                  "Missense Substitution" = "#33A02C", 
                  "Nonsense Substitution" = "#FF7F00", 
                  "Promoter Substitution" = "#FDBF6F", 
                  #"Splice Site Substitution" = "#CAB2D6", 
                  "Splice Site Substitution" = "#B1A0BA", 
                  "Splice Site Deletion" = "blue", # change?
                  "Nonsense Insertion" = 'darkred', # change?
                  'Multi-hit' = 'black',
                  # Simpified SNVs
                  "Substitution" = "#33A02C", 
                  "Insertion/Deletion" = "#FF7F00", 
                  'SNV' = "black", 
                  #'background' = "grey80"
                  'background' = "grey82"
                  )

# Oncoplot cell and variant sizes
size.cell <- 1
size.cnv <- unit(1, 'pt')
size.snv <- .33

# Oncoplot alteration function
alter_fun = list(
  "Gain" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(size.cell, "pt"), h - size.cnv, 
              gp = gpar(fill = variant_cols["Gain"], col = NA))
  }, 
  "Loss" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(size.cell, "pt"), h - size.cnv, 
              gp = gpar(fill = variant_cols["Loss"], col = NA))
  }, 
  "Frameshift Deletion" = function (x, y, w, h) { 
    grid.rect(x, y, w - unit(size.cell, "pt"), h*size.snv, 
              gp = gpar(fill = variant_cols["Frameshift Deletion"], col = NA))
  }, 
  "Frameshift Insertion" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(size.cell, "pt"), h*size.snv, 
              gp = gpar(fill = variant_cols["Frameshift Insertion"], col = NA))
  }, 
  "In Frame Deletion" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(size.cell, "pt"), h*size.snv, 
              gp = gpar(fill = variant_cols["In Frame Deletion"], col = NA))
  }, 
  "In Frame Insertion" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(size.cell, "pt"), h*size.snv, 
              gp = gpar(fill = variant_cols["In Frame Insertion"], col = NA))
  }, 
  "In Frame Indel" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(size.cell, "pt"), h*size.snv, 
              gp = gpar(fill = variant_cols["In Frame Indel"], col = NA))
  }, 
  "Nonsense Insertion" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(size.cell, "pt"), h*size.snv, 
              gp = gpar(fill = variant_cols["Nonsense Insertion"], col = NA))
  }, 
  "Intron Substitution" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(size.cell, "pt"), h*size.snv, 
              gp = gpar(fill = variant_cols["Intron Substitution"], col = NA))
  }, 
  "Missense Substitution" = function (x, y, w, h)  {
    grid.rect(x, y, w - unit(size.cell, "pt"), h*size.snv, 
              gp = gpar(fill = variant_cols["Missense Substitution"], col = NA))
  }, 
  "Nonsense Substitution" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(size.cell, "pt"), h*size.snv, 
              gp = gpar(fill = variant_cols["Nonsense Substitution"], col = NA))
  }, 
  "Promoter Substitution" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(size.cell, "pt"), h*size.snv, 
              gp = gpar(fill = variant_cols["Promoter Substitution"], col = NA))
  }, 
  "Splice Site Substitution" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(size.cell, "pt"), h*size.snv, 
              gp = gpar(fill = variant_cols["Splice Site Substitution"], col = NA))
  },
  "Splice Site Deletion" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(size.cell, "pt"), h*size.snv, 
              gp = gpar(fill = variant_cols["Splice Site Deletion"], col = NA))
  },
  "Substitution" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(size.cell, "pt"), h*size.snv, 
              gp = gpar(fill = variant_cols["Substitution"], col = NA))
  }, 
  "Insertion/Deletion" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(size.cell, "pt"), h*size.snv, 
              gp = gpar(fill = variant_cols["Insertion/Deletion"], col = NA))
  }, 
  "SNV" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(size.cell, "pt"), h*size.snv, 
              gp = gpar(fill = variant_cols["SNV"], col = NA))
  }, 
  "Multi-hit" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(size.cell, "pt"), h*size.snv, 
              gp = gpar(fill = variant_cols["Multi-hit"], col = NA))
  }, 
  "ND_CNV" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(size.cell, "pt"), h*size.snv, 
              gp = gpar(fill = variant_cols["ND_CNV"], col = NA))
  }, 
  "ND_SNV" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(size.cell, "pt"), h*size.snv, 
              gp = gpar(fill = variant_cols["ND_SNV"], col = NA))
  }, 
  "background" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(size.cell, "pt"), h - unit(size.cell, 'pt'), 
              gp = gpar(fill = variant_cols["background"], col = NA))
  }
)


# Function for saving oncoplot as png image
# TODO: Could probably just use save heatmap function if formmated correctly
draw_oncoprint <- function(ht, fn, pd, 
                           column_title = NULL) {
  
  # Draw dummy plot to determine figure size
  pdf(NULL)
  dht <- draw(ht,heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom',
              annotation_legend_list = pd, column_title = column_title)
  wh <- calc_ht_size(dht) # wh[1] = width, wh[2] = height
  try(dev.off(), silent = TRUE)
  
  # Set width, height, and resolution
  w <- wh[1] 
  h <- wh[2] 
  res <- (w*h)
  
  # Save as png
  png(filename = fn, width = w, height = h, units = 'in', res = res)
  draw(ht,heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom',
       annotation_legend_list = pd, column_title = column_title)
  dev.off()
  
}

# Function to sort columns using mutual exclusivity
meSort <- function(m, row_order = NULL) {
  
  # Order rows decreasing by frequency
  if (is.null(row_order)) {
    
    row_order <- sort(rowSums(m), decreasing=TRUE, index.return=TRUE)$ix
    
  }
  
  # Compute score for each olumn
  scoreCol <- function(x) {
    
    score <- 0
    
    for (i in 1:length(x)) {
      
      if (x[i]) {
        
        score <- score + 2^(length(x)-i);
        
      }
    }
    
    return(score)
    
  }
  
  # Get column scores
  scores <- apply(m[row_order, ], 2, scoreCol)
  
  # Set new order and return
  col_order <- sort(scores, decreasing=TRUE, index.return=TRUE)$ix
  
  return(col_order)
  
}




# Function for creating bar annotation showing column counts of variant types
onco_column_anno <- function(var_list) {
  
  # Compute column counts
  var_counts <- sapply(var_list, function(x) apply(x, 2, sum))
  
  # Make annotation object and return
  HeatmapAnnotation(cbar = anno_barplot(var_counts, height = unit(20, 'mm'), border = FALSE,
                                        gp = gpar(fill = variant_cols[names(var_list)], 
                                                  col = variant_cols[names(var_list)])), 
                    show_annotation_name = FALSE) %>%
    return()
  
}

# Function for creating bar annotation showing row counts of variant types
onco_row_anno <- function(var_list) {
  
  # Compute row counts
  var_counts <- sapply(var_list, function(x) apply(x, 1, sum))
  
  # Make annotation object and return
  rowAnnotation(rbar = anno_barplot(var_counts, border = TRUE, width = unit(3, "cm"),
                                    gp = gpar(fill = variant_cols[colnames(var_counts)]),
                                    axis_param = list(side = 'top', labels_rot = 0)),
                show_annotation_name = FALSE) %>%
    return()
  
}


# Function to make binary versions of each variant type
#   * Note samples with no CNV or SNV calls are not removed
# TODO: Make a more generalizable function that parses through alteration types from CNV and SNV matrices
binary_variant_matrices <- function(cnvs.dat, snvs.dat, merge_snvs = FALSE) {
  
  # Get all CNV types
  cnv_types <- cnvs.dat %>% as.matrix() %>% as.vector() %>% unique() 
  cnv_types <- cnv_types[!is.na(cnv_types)]
  cnv_types <- cnv_types[cnv_types != 'ND-CNV']
  
  # Get all SNV types
  snv_types <- snvs.dat %>% as.matrix() %>% as.vector() %>% unique() 
  snv_types <- snv_types[!is.na(snv_types)]
  snv_types <- snv_types[snv_types != 'ND-SNV']
  
  # Create matrices to store each variant type
  var_list <- list()
  
  # CNVs
  for (ct in cnv_types) {
    m <- matrix(0, nrow = dim(cnvs.dat)[1], ncol = dim(cnvs.dat)[2], dimnames = list(rownames(cnvs.dat), colnames(cnvs.dat)))
    m[as.matrix(cnvs.dat) %in% c(ct)] <- 1
    var_list[[ct]] <- t(m)
  }
  
  # SNVs
  if (merge_snvs) { # Merge into single matrix
    
    m <- matrix(0, nrow = dim(snvs.dat)[1], ncol = dim(snvs.dat)[2], dimnames = list(rownames(snvs.dat), colnames(snvs.dat)))
    m[as.matrix(snvs.dat) %in% snv_types] <- 1
    var_list[['SNV']] <- t(m)
    
  } else { # Else split by mutation type
    
    for (st in snv_types) { 
      m <- matrix(0, nrow = dim(snvs.dat)[1], ncol = dim(snvs.dat)[2], dimnames = list(rownames(snvs.dat), colnames(snvs.dat)))
      m[as.matrix(snvs.dat) %in% c(st)] <- 1
      var_list[[st]] <- t(m)
    }
    
  }
  
  
  # Samples missing CNV calls
  nd.cnvs.dat <- matrix(0, nrow = dim(cnvs.dat)[1], ncol = dim(cnvs.dat)[2], dimnames = list(rownames(cnvs.dat), colnames(cnvs.dat)))
  nd.cnvs.dat[as.matrix(cnvs.dat) %in% c('ND-CNV')] <- 1
  var_list[['ND-CNV']] <- t(nd.cnvs.dat)
  
  # Samples missing SNV calls
  nd.snvs.dat <- matrix(0, nrow = dim(snvs.dat)[1], ncol = dim(snvs.dat)[2], dimnames = list(rownames(snvs.dat), colnames(snvs.dat)))
  nd.snvs.dat[as.matrix(snvs.dat) %in% c('ND-SNV')] <- 1
  var_list[['ND-SNV']] <- t(nd.snvs.dat)
  
  # Return just SNV if merging
  if (merge_snvs) {snv_types <- c('SNV')}
  
  return(list(var_list, cnv_types, snv_types))
  
}

# Function for filtering variant matrices by min samples
#  Can be numeric to filter variants in min samples, but can also use named list to filter by
#   variants found in samples that have data for a specified column in a meta table, or
#    can be a named list with each item containing a named vector for specifiying a min number
#     of samples for each group
#  NOTE: NAMES LIST OPTION IS CURRENTLY SET UP AS AN "OR" OPTION --> VARIANTS MEETING ANY
#          OF SUCH THRESHOLDS ARE SELECTED
filter_variants <- function(var_list, min_vars = 2, meta.sub = NULL) {
  
  # Collapse list of variant matrices
  if (length(var_list) > 0) {
    
    var_mat <- do.call(pmax,var_list)
    
    # If min_vars is numeric, apply to all matrices in var_list
    if (is.numeric(min_vars)) { 
      
      select_vars <- rownames(var_mat[rowSums(var_mat) >= min_vars,,drop=FALSE])
      
    } else {
      
      # Start variant list
      select_vars <- c()
      
      # Parse named list
      for (c in names(min_vars)) {
        
        # Else if min_vars is a named list pointing towards a numeric value, select variants meeting threshold for meta column (exclude NA)
        if (is.null(names(min_vars[[c]]))) {
          
          thres <- min_vars[[c]]
          select_samples <- unique(intersect(rownames(meta.sub[!(is.na(meta.sub[,c])),]), colnames(var_mat)))
          var_sub <- var_mat[,select_samples,drop=FALSE]
          select_vars <- c(select_vars, rownames(var_sub[rowSums(var_sub) >= thres,]))
          
          # Else if named list values contain named vectors, use named vectors to filter based on group  
        } else {
          
          # Parse named vector
          for (g in names(min_vars[[c]])) {
            
            thres <- min_vars[[c]][g]
            select_samples <- unique(intersect(rownames(meta.sub[!(is.na(meta.sub[,c])) & meta.sub[,c] == g,]), colnames(var_mat)))
            var_sub <- var_mat[,select_samples,drop=FALSE]
            select_vars <- c(select_vars, rownames(var_sub[rowSums(var_sub) >= thres,,drop=FALSE]))
            
          }
        }
      }
    }
    
    # Remove duplicates
    select_vars <- sort(unique(select_vars))
    
    # Now update each matrix of list to only include select variants
    for (n in names(var_list)) {var_list[[n]] <- var_list[[n]][select_vars,,drop=FALSE] }
    
  }
  
  return(var_list)
  
}


# Function to create lists of matrices for each variant type including matrices indicating potentially missing data
variant_matrix_lists <- function(biMats, cnv_types, snv_types) { 
  
  # Get any samples with missing data
  no_cnv_samples <- colnames(biMats[['ND-CNV']])[colSums(biMats[['ND-CNV']]) > 0]
  no_snv_samples <- colnames(biMats[['ND-SNV']])[colSums(biMats[['ND-SNV']]) > 0]
  no_var_samples <- union(no_cnv_samples, no_snv_samples)
  
  # All alteration types
  var_list.all <- list()
  for (v in c(cnv_types, snv_types)) { var_list.all[[v]] <- biMats[[v]][,!colnames(biMats[[v]]) %in% no_var_samples,drop=FALSE] }
  var_list.all <- unify_mat_list(var_list.all)
  
  # CNVs only
  var_list.cnvs <- list()
  for (ct in cnv_types) { var_list.cnvs[[ct]] <- biMats[[ct]][,!colnames(biMats[[ct]]) %in% no_cnv_samples,drop=FALSE] }
  var_list.cnvs <- unify_mat_list(var_list.cnvs)
  
  # SNVs only
  var_list.snvs <- list()
  for (st in snv_types) { var_list.snvs[[st]] <- biMats[[st]][,!colnames(biMats[[st]]) %in% no_snv_samples,drop=FALSE] }
  var_list.snvs <- unify_mat_list(var_list.snvs)
  
  # Include all missing data
  var_list.missing <- unify_mat_list(biMats)
  
  return(list(var_list.all, var_list.cnvs, var_list.snvs, var_list.missing))
  
}


# A convienence function for processing variant data into lists of 
#   binary matrices suitable for oncoprint and other functions
var_list_pipeline <- function(cnvs.dat, snvs.dat, meta, min_vars = 1, 
                              select_samples = NULL, select_variants = NULL) {
  
  # Select only specific variants if specified
  if (!is.null(select_variants)) {
    
    cnvs.dat <- cnvs.dat[,colnames(cnvs.dat) %in% select_variants,drop = FALSE]
    snvs.dat <- snvs.dat[,colnames(snvs.dat) %in% select_variants,drop = FALSE]
    
  }   
  
  # Subset down to select samples
  if (!is.null(select_samples)) {
    
    cnvs.dat <- cnvs.dat[select_samples,,drop=FALSE]
    snvs.dat <- snvs.dat[select_samples,,drop=FALSE]
    meta.sub <- meta[select_samples,,drop=FALSE]
    
  } else {
    
    meta.sub <- meta
    
  }
  
  # Create binary matrices for each variant type
  bvm.return <- binary_variant_matrices(cnvs.dat, snvs.dat)
  biMats <- bvm.return[[1]]
  cnv_types <- bvm.return[[2]]
  snv_types <- bvm.return[[3]]
  
  # Create unified lists for all variant types (each plot type)
  var_lists <- variant_matrix_lists(biMats, cnv_types, snv_types)
  var_list.all <- var_lists[[1]]
  var_list.cnvs <- var_lists[[2]]
  var_list.snvs <- var_lists[[3]]
  var_list.missing <- var_lists[[4]]
  
  # Filter variant matrices using thresholds
  var_list.all <- filter_variants(var_list.all, min_vars = min_vars, meta.sub)
  var_list.cnvs <- filter_variants(var_list.cnvs, min_vars = min_vars, meta.sub)
  var_list.snvs <- filter_variants(var_list.snvs, min_vars = min_vars, meta.sub)
  
  # Only use variants from joint CNV/SNV for missing matrices
  if (length(var_list.all) >= 1) {
    
    keep_vars <- rownames(var_list.all[[1]])
    var_list.missing_sub <- list()
    
    for (vt in names(var_list.missing)) {  var_list.missing[[vt]] <- var_list.missing[[vt]][keep_vars,,drop=FALSE]  }
    
    # Change to underscore for missing data
    names(var_list.missing) <- gsub('-', '_', names(var_list.missing))
    
  } else {
    
    var_list.missing <- list()
    
  }
  
  return(list(var_list.all, var_list.cnvs, var_list.snvs, var_list.missing, meta.sub))
  
}






# Function for gathering all components of oncoprint and setting oncoprint parameters before drawing
build_oncoprint <- function(var_list, meta, select_samples = NULL, fn = NULL, heatmap_title = '', 
                            
                            top_anno = NULL, 
                            bottom_anno = NULL, 
                            
                            cluster_columns = FALSE, 
                            column_order = NULL, 
                            fix_order = FALSE, 
                            keep_original_column_order = TRUE, 
                            column_barplot = TRUE,
                            
                            column_split = NULL, 
                            column_split_order = NULL,
                            column_title = NULL,
                            show_column_split_titles = TRUE,
                            column_split_fill_cols = NULL,
                            column_split_title_border = NULL,
                            column_title_side = 'top',
                            column_title_gp = gpar(fontsize = 14),
                            column_gap = unit(2.5, "mm"), # library default is unit(1, "mm"),
                            ht_border = FALSE, 
                            gap_border = TRUE,
                            
                            category_table = NULL, cat_order = NULL,
                            
                            cluster_rows = FALSE, 
                            row_barplot = TRUE,
                            show_pct = TRUE, 
                            pct_side = 'right', 
                            row_names_side = 'left',
                            
                            lgd_fontsize = 12, 
                            lgd_gridsize = 4,
                            lgd_rows = 3,
                            lgd_gap = unit(2, 'mm'),
                            lgd_pack_gap = unit(2, 'mm'),
                          
                            ht_gap = unit(4, "mm"),
                            
                            ht_height = NULL, ht_width = NULL) {
  
  # Compute binary matrix for row and column ordering
  biMat <- do.call(pmax,var_list)
  
  # Compute summation matrix for row ordering tie breaker
  sumMat <- Reduce('+', var_list)
  
  # Make sure at least one variant remains after subsetting/filtering step
  if (nrow(biMat) > 0) {
    
    
    ########### LEGEND PARAMETERS ##################
    
    # TODO: Creating oncoLgd list first no longer needed?
    # Create oncoplot legend (passed to oncoPrint function)
    oncoLgd = list(title = "Alterations", at = names(variant_cols), 
                   
                   labels = names(variant_cols), nrow = lgd_rows,
                   labels_gp = gpar(fontsize = lgd_fontsize),
                   
                   grid_height = unit(lgd_gridsize, 'mm'), 
                   grid_width = unit(lgd_gridsize, 'mm'), 
                   gap = lgd_gap,
                   
                   title_gp = gpar(fontsize = lgd_fontsize, fontface = "bold"))
    
    # Same legend but for drawing png image
    varLgd <- Legend(title = oncoLgd[['title']], labels = names(var_list), nrow = oncoLgd[['nrow']],
                     legend_gp = gpar(fill = variant_cols[names(var_list)]), 
                     
                     grid_height = oncoLgd[['grid_height']], grid_width = oncoLgd[['grid_width']],
                     gap = lgd_gap,
                     
                     title_gp = oncoLgd[['title_gp']], labels_gp = oncoLgd[['labels_gp']])
    
    # Organize all legends into single list and pack
    lgd_list <- make_legend_list(row_lgd_list = NULL, top_lgd_list = top_anno[[2]], 
                                 ht_lgd = varLgd, btm_lgd_list = bottom_anno[[2]])
    pd <- packLegend(list = lgd_list, direction = 'horizontal', 
                     column_gap = lgd_pack_gap, 
                     
                     max_width = unit(13, 'in'))
    
    ######### END LEGEND PARAMETERS #################
    
    
    
    ######### ROW PARAMETERS ###########
    
    # Set row order first by (binary) percent altered, then by total alt types per gene, then by name
    order_rows <- data.frame(Gene = rownames(biMat),
                             Percent = rowSums(biMat), 
                             Total = rowSums(sumMat)) %>%
      mutate(Gene = as.character(Gene)) %>%
      arrange(-Percent, -Total, Gene) %>%
      pull(Gene)
    
    # Split by categories if table provided
    if (!is.null(category_table)) {
      
      var.df <- data.frame(Gene = rownames(biMat))
      category_table <- left_join(var.df, category_table)
      rownames(category_table) <- category_table$Gene
      category_table <- category_table[rownames(biMat),]
      split_rows <- category_table[,2]
      ht_border <- TRUE
      
      # Order categories by order parameter or by highest count
      if (is.null(cat_order)) {
        
        cat_order <- unique(category_table[order_rows,2])
        
      } 
      
      split_rows <- factor(split_rows, levels = cat_order) 
      
    } else {
      
      split_rows <- NULL
      
    }
  
    # Build row annotations for right side of oncoplot
    if (row_barplot) {
      
      right_anno <- onco_row_anno(var_list)  
      
    }  else {
      
      right_anno <- NULL
      
    }
  
    
    # Cluster rows if specified
    if (cluster_rows & nrow(biMat) >= 2) { 
      
      cluster_rows <- create_dendrogram(biMat)
   
    } else {
      
      cluster_rows <- NULL
      
    }
    
    
    ########## END ROW PARAMETERS ##################
    
    
    ############ COLUMN PARAMETERS ###################
    
    # Subset top_annotations if provided
    meta$AnnoIndex <- 1:nrow(meta)
    top_anno <- update_annotations(meta[colnames(biMat),], top_anno)
    bottom_anno <- update_annotations(meta[colnames(biMat),], bottom_anno)
    
    # Split oncoplot by phenotype
    if (!is.null(column_split)) {
      
      column_split <- meta[colnames(biMat),column_split]
      
      # Set order if provided, else sort alphabetically
      if (!is.null(column_split_order)) {
        
        column_split <- factor(column_split, levels = column_split_order)
        
      } else {
        
        column_split_order <- column_split %>% sort() %>% unique()
        column_split <- factor(column_split, levels = unique(column_split_order))
        
      }
      
      # Set split titles using column_title parameter
      if (show_column_split_titles) {
        
        column_title = unique(column_split_order)
        
        # Set colors for column splits
        if (is.null(column_split_fill_cols)) {
          
          column_split_fill_cols <- 'white'
          if (is.null(column_split_title_border)) {column_split_title_border <- FALSE}
          
          
        } else {
          
          # Alt black and white
          #column_split_fill_cols <- rep(c('black', 'white'), length.out = length(column_title))
          column_split_fill_cols <- column_split_fill_cols[column_title]
          if (is.null(column_split_title_border)) {column_split_title_border <- TRUE}
          
        }
        
        # Set column split header parameters
        column_title_gp = gpar(fill = column_split_fill_cols,
                               col = sapply(column_split_fill_cols, FUN = isDark),
                               border = column_split_title_border)
         
      } 

      # Add gap between splits
      if (gap_border) {ht_border = TRUE}
      
    }
    
    ######## COLUMN ORDERING, SPLITTING, AND CLUSTERING ##########
    
    # SET COLUMN ORDERING
    if (fix_order | cluster_columns) {
      
      if (is.null(select_samples)) {select_samples <- colnames(biMat)}
      col_order <- select_samples
      
    } else {
      
      # If splitting columns by group, ME sort each group 
      #  independently while keeping row order fixed
      if (!is.null(column_split)) {
        
        col_order <- c()
        column_groups <- unique(column_split)
        for (cg in sort(column_groups)) {
          
          biMat.cg <- biMat[,column_split == cg]
          cg.order <- colnames(biMat.cg)[meSort(biMat.cg, 
                                                #row_order = NULL
                                                row_order = order_rows
                                                )]
          col_order <- c(col_order, cg.order)
          
        }
        
      } else {
        
        # Order columns using mutually exclusive sort
        col_order <- colnames(biMat)[meSort(biMat)]
        
      }
      
    }
    
    # Build empty top heatmap object for ordering, clustering, and splitting
    top_anno.pre <- Heatmap(biMat, 
                            col = c('white', 'white'),
                            show_heatmap_legend = FALSE,
                            show_row_names = FALSE,
                            cluster_rows = FALSE,
                            height = unit(0, 'in'),
                            heatmap_width = ht_width, 
                            column_order = col_order,
                            column_split = column_split,
                            column_title = column_title,
                            column_title_gp = column_title_gp,
                            column_title_side = column_title_side,
                            column_gap = column_gap,
                            show_parent_dend_line = FALSE,
                            cluster_columns = cluster_columns) 
    
    # Make column barplot annotation
    if (column_barplot) {
      
      top_anno.bar <- onco_column_anno(var_list)  
      
    } else {
      
      top_anno.bar <- NULL
      
    }
    
    
    # Stack annotations
    # TODO: Make option for putting split titles below barplot
    if (column_title_side == 'top') {
      
      top_anno.final <- top_anno.pre %v% top_anno.bar %v% top_anno[[1]]
      
    } else {
      
      top_anno.final <- top_anno.bar %v% top_anno[[1]] %v% top_anno.pre
      
    }
    
 
    ############### END COLUMN PARAMETERS ####################
    
 
    # Set dimensions of main heatmap if not specified
    if (is.null(ht_height)) {
      
      if (nrow(biMat) > 15) {
        
        ht_height <- unit(11, 'in')
        
      } else {
        
        ht_height <- unit(4, 'in')
        
      }
      
    }
    
    if (is.null(ht_width)) {
      
      ht_width <- unit(18, 'in')
      
    }
    
    # So column barplot isn't created automatically by oncoPrint function
    faux_top_anno <- HeatmapAnnotation(foo = anno_empty(border = FALSE, height = unit(.000001, "mm")))
    
    # Create initial set of arguments for oncoPrint
    onco_args <- list(mat = var_list, 
                      alter_fun = alter_fun, 
                      show_heatmap_legend = FALSE,
                      remove_empty_columns = FALSE, 
                      remove_empty_rows = FALSE, 
                      show_column_names = FALSE,
                      column_names_side = 'bottom',
                      
                      border = ht_border,
                      
                      top_annotation = faux_top_anno, 
                      right_annotation = right_anno,
                      
                      row_order = order_rows,
                      cluster_rows = cluster_rows,
                      
                      row_split = split_rows, 
                      row_title_rot = 0, 
                      row_gap = unit(2, "mm"), 
                      show_pct = show_pct,
                      pct_side = pct_side,
                      row_names_side = row_names_side,
                      
                      height = ht_height,
                      width = ht_width,
                      col = variant_cols)
    
    # Set cluster_rows parameter (for some reason it breaks if not set after creating initial list?)
    onco_args[['cluster_rows']] <- cluster_rows
 
    
    # TEST PROVIDING REAL ANNOTATION INSTEAD OF FAUX ANNOTATION
    #print('TEST TOP ANNOTATION')
    #onco_args[['top_annotation']] <- onco_column_anno(var_list)
    
    
    # Create oncoPrint heatmap object
    #oncoHT <- purrr::invoke(oncoPrint, .args = onco_args)
    oncoHT <- rlang::invoke(oncoPrint, .args = onco_args)
    
    # Append top/bottom annotations 
    ht <- top_anno.final %v% oncoHT %v% bottom_anno[[1]]
    
    # TODO: Merge draw_oncoprint and save_htan_heatmap
    #        remove hardcode gaps and create parameter
    
    # Draw and save oncoplot to file
    #if (!is.null(fn)) {draw_oncoprint(ht, fn, pd, heatmap_title)}
    
    if (!is.null(fn)) {save_htan_heatmap(list(ht, lgd_list), fn, 
                                         extend_w = 0, res_factor = 1,
                                         ht_gap = ht_gap
                                         )}
    
    # Return a second oncoPrint object containing just the 
    #   heatmap (and legend object) without annotation objects
    onco_args[['show_pct']] <- FALSE
    onco_args[['right_annotation']] <- rowAnnotation(foo = anno_empty(height = unit(.000001, "mm"),
                                                                      border = FALSE))
    
    # Keep original column order when returning (so can be concatenated with other objects later)
    if (keep_original_column_order) {
      
      onco_args[['column_order']] <- 1:ncol(var_list[[1]])
    
    }
    
    # Hide row dendrogram
    onco_args[['show_row_dend']] <- FALSE
    
    # Run oncoPrint() with list of parameters
    #oncoHT <- purrr::invoke(oncoPrint, .args = onco_args)
    oncoHT <- rlang::invoke(oncoPrint, .args = onco_args)
    
    return(list(oncoHT,varLgd))
    
  } else {
    
    return(list(NULL,NULL))
    
  }
  
}


# Main function for formatting alteration data and generating 
#   oncoplots using oncoPrint() from ComplexHeatmap
make_oncoplots <- function(cnvs.dat, snvs.dat, meta, select_samples = NULL, 
                           select_variants = NULL, min_vars = 1, 
                           top_anno = NULL, bottom_anno = NULL, 
                           ht_height = NULL, ht_width = NULL, 
                           heatmap_title = NULL, 
                           
                           ht_gap = unit(2, "mm"),
                           lgd_fontsize = 12, 
                           lgd_gridsize = 4,
                           lgd_rows = 3,
                           lgd_gap = unit(2, 'mm'),
                           lgd_pack_gap = unit(2, 'mm'),
                           
                           column_barplot = TRUE,
                           column_split = NULL, 
                           column_split_order = NULL, 
                           show_column_split_titles = TRUE,
                           column_title_side = 'top',
                           column_split_fill_cols = NULL,
                           column_split_title_border = NULL,
                           gap_border = TRUE,
                           cluster_columns = FALSE, 
                           
                           fix_order = FALSE, keep_original_column_order = TRUE,
                           
                           all_variants = TRUE, cnvs_only = FALSE, snvs_only = FALSE,
                           
                           cluster_rows = FALSE, row_names_side = 'left', 
                           show_pct = TRUE, pct_side = 'right', row_barplot = TRUE,
                           category_table = NULL, cat_order = NULL,
                           pre = NULL, return_objects = FALSE) {
  
  # Format gene categories if table provided
  if (!is.null(category_table)) {
    
    category_table <- category_table %>%
      mutate(Category = gsub('_', '\n', Category)) %>%
      mutate(Category = gsub(' ', '\n', Category)) %>%
      mutate(Category = gsub('-', '\n', Category))
    
  }
  
  # Process all data into unified lists for oncoprint
  var_lists <- var_list_pipeline(t(cnvs.dat), t(snvs.dat), meta, min_vars = min_vars, 
                                 select_samples = select_samples, select_variants = select_variants)
  var_list.all <- var_lists[[1]]
  var_list.cnvs <- var_lists[[2]]
  var_list.snvs <- var_lists[[3]]
  meta.sub <- var_lists[[5]]
  
  # CNVS AND SNVS
  if (length(var_list.all) >= 1 & all_variants) {
    
    if (nrow(do.call(pmax,var_list.all)) >= 1) {
      
      if (!is.null(pre)) {fn <- paste(pre, ".png", sep = '')} else {fn <- NULL}
      
      oncoHT.all <- build_oncoprint(var_list.all, meta, meta.sub, select_samples = select_samples, 
                                    fn = fn, 
                                    heatmap_title = heatmap_title, 
                                    top_anno = top_anno, bottom_anno = bottom_anno, 
                                    
                                    lgd_fontsize = lgd_fontsize, 
                                    lgd_gridsize = lgd_gridsize,
                                    lgd_rows = lgd_rows,
                                    lgd_gap = lgd_gap,
                                    lgd_pack_gap = lgd_pack_gap, 
                                    
                                    fix_order = fix_order, cluster_columns = cluster_columns, cluster_rows = cluster_rows, 
                                    
                                    column_barplot = column_barplot,
                                    column_split = column_split, column_split_order = column_split_order, 
                                    show_column_split_titles = show_column_split_titles, gap_border = gap_border,
                                    column_split_fill_cols = column_split_fill_cols, column_title_side = column_title_side,
                                    column_split_title_border = column_split_title_border,
                                    
                                    keep_original_column_order = keep_original_column_order,
                                    category_table = category_table, cat_order = cat_order,
                                    show_pct = show_pct, pct_side = pct_side, row_names_side = row_names_side,
                                    row_barplot = row_barplot,
                                    
                                    ht_gap = ht_gap,
                                    ht_height = ht_height, ht_width = ht_width) 
      
    } 
    
  } else {
    
    oncoHT.all <- NULL
    
  }
  
  # CNVS ONLY
  if (length(var_list.cnvs) >= 1 & cnvs_only) {
    
    if (nrow(do.call(pmax,var_list.cnvs)) >= 1) {
      
      if (!is.null(pre)) {fn <- paste(pre, "_cnvs.png", sep = '')} else {fn <- NULL}
      
      
      oncoHT.cnvs <- build_oncoprint(var_list.cnvs, meta, meta.sub, select_samples = select_samples, 
                                     fn = fn, heatmap_title = heatmap_title, top_anno = top_anno, 
                                     bottom_anno = bottom_anno, 
                                     
                                     lgd_fontsize = lgd_fontsize, 
                                     lgd_gridsize = lgd_gridsize,
                                     lgd_rows = lgd_rows,
                                     lgd_gap = lgd_gap,
                                     lgd_pack_gap = lgd_pack_gap, 
                                     
                                     fix_order = fix_order, cluster_columns = cluster_columns, cluster_rows = cluster_rows, 
                                     
                                     column_barplot = column_barplot,
                                     column_split = column_split, column_split_order = column_split_order, 
                                     show_column_split_titles = show_column_split_titles, gap_border = gap_border,
                                     column_split_fill_cols = column_split_fill_cols, column_title_side = column_title_side,
                                     column_split_title_border = column_split_title_border,
                                     
                                     keep_original_column_order = keep_original_column_order,
                                     category_table = category_table, cat_order = cat_order,
                                     show_pct = show_pct, pct_side = pct_side, row_names_side = row_names_side,
                                     row_barplot = row_barplot,
                                     ht_gap = ht_gap,
                                     ht_height = ht_height, ht_width = ht_width) 
    } 
    
  } else {
    
    oncoHT.cnvs <- NULL
    
  }
  
  # SNVS ONLY
  if (length(var_list.snvs) >= 1 & snvs_only) {
    
    if (nrow(do.call(pmax,var_list.snvs)) >= 1) {
      
      if (!is.null(pre)) {fn <- paste(pre, "_snvs.png", sep = '')} else {fn <- NULL}
      oncoHT.snvs <- build_oncoprint(var_list.snvs, meta, meta.sub, select_samples = select_samples, 
                                     
                                     fn = fn, heatmap_title = heatmap_title, top_anno = top_anno, 
                                     bottom_anno = bottom_anno, 
                                     
                                     lgd_fontsize = lgd_fontsize, 
                                     lgd_gridsize = lgd_gridsize,
                                     lgd_rows = lgd_rows,
                                     lgd_gap = lgd_gap,
                                     lgd_pack_gap = lgd_pack_gap, 
                                     
                                     fix_order = fix_order, cluster_columns = cluster_columns, cluster_rows = cluster_rows, 
                                     
                                     column_barplot = column_barplot,
                                     column_split = column_split, column_split_order = column_split_order, 
                                     show_column_split_titles = show_column_split_titles, gap_border = gap_border,
                                     column_split_fill_cols = column_split_fill_cols, column_title_side = column_title_side,
                                     column_split_title_border = column_split_title_border,
                                     
                                     keep_original_column_order = keep_original_column_order,
                                     category_table = category_table, cat_order = cat_order,
                                     show_pct = show_pct, pct_side = pct_side, row_names_side = row_names_side,
                                     row_barplot = row_barplot,
                                     ht_gap = ht_gap,
                                     ht_height = ht_height, ht_width = ht_width) 
    } 
    
  } else {
    
    oncoHT.snvs <- NULL
    
  }
  
  # Return oncoplot objects
  if (return_objects) {
    
    return(list(oncoHT.all, 
                oncoHT.cnvs, 
                oncoHT.snvs))
    
  }
  
}

