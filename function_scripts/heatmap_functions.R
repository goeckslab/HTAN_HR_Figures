

######################################################
#
#   HEATMAP FUNCTIONS FOR HTAN HR+ MANUSCRIPT
#
######################################################

library(ComplexHeatmap)
library(circlize)
library(tidyr)
library(dplyr)
library(yarrr)


# Annotation colors
colors.treatment <- c("Abemaciclib" = "#EE0011FF", 
                      "Palbociclib" = "#0C5BB0FF", 
                      "Everolimus" = "#7F8624FF", 
                      "Fulvestrant" = "#FA6B09FF", 
                      "Letrozole" = "#15983DFF", 
                      "Tamoxifen" = "maroon", 
                      "Post-Abemaciclib" = "#FFC6C7FF", 
                      "Post-Palbociclib" = "#B7CBFFFF", 
                      "Post-Ribociclib" = "#87EF9AFF",
                      "Post-Everolimus" = "#EBF09C", 
                      "Post-Fulvestrant" = "#FFE7D4", 
                      "Post-Letrozole" = "#8EFFBA", 
                      "Post-Tamoxifen" = "#FFD3E5",
                      'Pre-Palbociclib' = 'white',
                      'Pre-Abemaciclib' = 'white',
                      'Pre-Ribociclib' = 'white',
                      'Pre-Fulvestrant' = 'white',
                      'Pre-Letrozole' = 'white',
                      'Pre-Tamoxifen' = 'white')
                      
colors.response <- c("Responder" = "#11776CFF", 
                     "Non-responder" = "#94220EFF")  

colors.pam <- c("Basal" = "#A3DA4BFF", 
                "Her2" = "#E9D738FF", 
                "LumA" = "#4D709CFF", 
                "LumB" = "#B91226FF",
                "LumA:LumB" = "#8AB8CFFF", 
                "LumB:LumA" = "tomato", 
                "LumA:Her2" = "#7F8624FF")

colors.site <- c("Bone" = "#C4BB90FF", 
                 "Liver" = "#5FB233FF", 
                 "Lymph Node" = "#F19C1FFF")

colors.patient <- c("9-1" = "#F57206FF", 
                    "9-2" = "#8F2F8BFF", 
                    "9-3" = "red", 
                    "9-14" = "darkcyan", 
                    "9-15" = "#434159FF")

colors.gsva <- c("Cellular Component" = "#EC579AFF", 
                 "Development" = "#7F8624FF", 
                 "DNA Damage" = "#0C5BB0FF", 
                 "Immune Cells" = "navy", 
                 "Immune Signaling" = "#EE0011FF", 
                 "JAK/STAT" = "turquoise3",
                 "Metabolic" = "#149BEDFF", 
                 "Pathway" = "#FA6B09FF", 
                 "Proliferation" = "#FEC10BFF", 
                 "Receptor Signaling" = "#15983DFF", 
                 "Signaling" = "purple")

# Function to determine darkness of background color for overlaying black and white text
isDark <- function(colr) { 
  
  if (sum(col2rgb(colr) * c(299, 587,114))/1000 < 123) {
    
    return("White")
    
  } else {
    
    return("Black")
  }
  
}

# Function to assign colors to vector of items using yarrr piratepal()
set_colors <- function(x, col_name){
  
  my_colors <- setNames(piratepal(palette = col_name, length.out = length(unique(x))),unique(x))
  return(my_colors[sort(names(my_colors))])
  
}

# Function for building heatmap annotation legend objects
make_heatmap_legends <- function(meta, select_samples = NULL, lgd_rows = 2,
                                 treatment_title = ' ', strip_names = TRUE, rm_duplicate_pair = TRUE,
                                 lgd_fontsize = 12, lgd_gridsize = 4) {
  
  if (is.null(select_samples)) {
    
    select_samples <- meta %>% pull(Sample)
    
  }
  
  # All color schemes subsetted down to select patients
  if (!is.null(select_samples)) {
    
    colors.treatment <- colors.treatment[names(colors.treatment) %in% c(meta[select_samples, 'CDKi'], 
                                                                        meta[select_samples, 'ERi'])]
    colors.response <- colors.response[names(colors.response) %in% meta[select_samples, 'PatientResponse']]
    colors.pam <- colors.pam[names(colors.pam) %in% meta[select_samples, 'PAM50']]
    colors.pam_change <- colors.pam[names(colors.pam) %in% meta[select_samples, 'pamChange']]
    colors.site <- colors.site[names(colors.site) %in% meta[select_samples, 'BiopsySite']]
    colors.patient <- colors.patient[names(colors.patient) %in% meta[select_samples, 'Patient']]
    
  }
  
  # Legend for pre-treatment and on-progression with on-progression labeled with asterisk
  opAstrLgd <- Legend(title = 'Treatment', type = 'points', pch = c(NA,8), 
                      at = c('PreTreatment', 'OnProgression'), background = c('white'),
                      labels_gp = gpar(fontsize = lgd_fontsize), 
                      grid_height = unit(lgd_gridsize, 'mm'), grid_width = unit(lgd_gridsize, 'mm'),
                      title_gp = gpar(fontsize = lgd_fontsize, fontface = "bold"))
  
  # Legend for CDK4/6i and ERi treatment (merged with above)
  onProgLgd <- Legend(title = treatment_title, labels = names(colors.treatment), nrow = lgd_rows,
                      legend_gp = gpar(fill = colors.treatment),  labels_gp = gpar(fontsize = lgd_fontsize), 
                      grid_height = unit(lgd_gridsize, 'mm'), grid_width = unit(lgd_gridsize, 'mm'), 
                      title_gp = gpar(fontsize = lgd_fontsize, fontface = "bold"))
  
  # Patient response legend
  responseLgd <- Legend(title = 'Patient Response', labels = names(colors.response), nrow = lgd_rows, 
                        legend_gp = gpar(fill = colors.response), labels_gp = gpar(fontsize = lgd_fontsize), 
                        grid_height = unit(lgd_gridsize, 'mm'), grid_width = unit(lgd_gridsize, 'mm'),
                         title_gp = gpar(fontsize = lgd_fontsize, fontface = "bold"))
  
  # PAM50 annotation legend
  pamLgd <- Legend(title = 'PAM50', labels = names(colors.pam), nrow = lgd_rows, 
                   legend_gp = gpar(fill = colors.pam), labels_gp = gpar(fontsize = lgd_fontsize), 
                   grid_height = unit(lgd_gridsize, 'mm'), grid_width = unit(lgd_gridsize, 'mm'),
                   title_gp = gpar(fontsize = lgd_fontsize, fontface = "bold"))
  
  # PAM50 change annotation legend
  pamChangeLgd <- Legend(title = 'PAM50', labels = names(colors.pam_change), nrow = lgd_rows, 
                         legend_gp = gpar(fill = colors.pam_change), labels_gp = gpar(fontsize = lgd_fontsize), 
                         grid_height = unit(lgd_gridsize, 'mm'), grid_width = unit(lgd_gridsize, 'mm'),
                         title_gp = gpar(fontsize = lgd_fontsize, fontface = "bold"))
  
  # Metastatic site legend
  siteLgd <- Legend(title = 'Biopsy Site', labels = names(colors.site), nrow = lgd_rows, 
                    legend_gp = gpar(fill = colors.site), labels_gp = gpar(fontsize = lgd_fontsize), 
                    grid_height = unit(lgd_gridsize, 'mm'), grid_width = unit(lgd_gridsize, 'mm'),
                    title_gp = gpar(fontsize = lgd_fontsize, fontface = "bold"))
  
  # Patient legend
  patientLgd <- Legend(title = 'Patient ID', labels = names(colors.patient), nrow = lgd_rows,
                       legend_gp = gpar(fill = colors.patient), labels_gp = gpar(fontsize = lgd_fontsize), 
                       grid_height = unit(lgd_gridsize, 'mm'), grid_width = unit(lgd_gridsize, 'mm'),
                       title_gp = gpar(fontsize = lgd_fontsize, fontface = "bold"))
  
  return(list('opAstrLgd' = opAstrLgd, 
              'onProgLgd' = onProgLgd,
              'responseLgd' = responseLgd, 
              'pamLgd' = pamLgd, 
              'pamChangeLgd' = pamChangeLgd,
              'patientLgd' = patientLgd, 
              'siteLgd' = siteLgd))
  
}

# Function to subset annotation objects to select samples
update_annotations <- function(meta, annotation_objects) {
  
  # Pull out annotation from list and leave legend
  anno <- annotation_objects[[1]]
  
  # Make sure annotation object exists
  if (!is.null(anno)) {
    
    # If only one annotation, just subset
    if (length(anno) == 1) {
      
      anno_update <- anno[meta$AnnoIndex,]
      
    }
    
    # If multiple annotations, start new object, iteratively subset and concatenate 
    else {
      
      anno_update <- NULL
      
      for (i in 1:length(anno)) {
        
        anno_update <- anno_update %v% anno[i][,meta$AnnoIndex]
        
      }
      
    }
    
    # Replace with subsetted annotations
    annotation_objects[[1]] <- anno_update
    
  }
  
  return(annotation_objects)
  
}

# Function for creating barplot column annotations for heatmaps
# TODO: Update axis parameters
make_barplot_annotation <- function(meta, anno, anno_height = unit(1, 'in'), 
                                    labs = seq(-0.5, 1, .5), ylim = c(-.8, 1.25),
                                    lSize = 0.5) {
  
  meta[,anno] <- as.numeric(meta[,anno])
  
  # Compute axis parameters if not provided
  if (is.null(labs) | is.null(ylim)) {
    
    m.ax <- ceiling(max(meta[,anno], na.rm = TRUE)*2) / 2
    m.in <- min(0, (floor(min(meta[,anno], na.rm = TRUE)*2) / 2))
    
    if (is.null(labs)) {
      
      labs <- seq(m.in,m.ax,lSize)
      
    }
    
    if (is.null(ylim)) {
      
      ylim <- c(m.in,m.ax)
      
    }
    
  }
  
  
  # Create barplot anotation
  # TODO: Recode this to make generalizable
  if (anno == 'g1Score') {
    
    bar_anno <- HeatmapAnnotation('G1 Arrest Score' = anno_barplot(meta[,anno],
                                                                   height = anno_height,
                                                                   axis_param = list(at = labs,
                                                                                     labels = labs),
                                                                   ylim = ylim))

  } else if (anno == 'CytoChange') {
    
    bar_anno <- HeatmapAnnotation('Cytolytic Change' = anno_barplot(meta[,anno],
                                                                   height = anno_height,
                                                                   axis_param = list(at = labs,
                                                                                     labels = labs),
                                                                   ylim = ylim))
    
  # TODO: Add cytolytic activity (no change)  
  } else if (anno == 'CytolyticActivity') {
    
    bar_anno <- HeatmapAnnotation('Cytolytic Activity' = anno_barplot(meta[,anno],
                                                                      height = anno_height,
                                                                      axis_param = list(at = labs,
                                                                                        labels = labs),
                                                                      ylim = ylim))
    
  } else {
    
    bar_anno <- HeatmapAnnotation(ha = anno_barplot(meta[,anno],
                                                    height = anno_height,
                                                    axis_param = list(at = labs,
                                                                      labels = labs),
                                                    ylim = ylim))
    
    names(bar_anno) <- anno
    
  }
  
  return(bar_anno)
  
}



# Function to make column ID annotations for heatmaps
make_heatmap_columnIDs <- function(meta) {
  
  # Create annotation for pointing to column names
  pointer <- HeatmapAnnotation(Sample = anno_simple(rep('0', nrow(meta)), pch = rep('|', nrow(meta)), 
                                                    pt_gp = gpar(col = "black", fill = 'white'), 
                                                    pt_size = unit(5, "mm"),col = c('0' = 'white')),
                               show_annotation_name = FALSE)
  
  # Sample IDs
  sampleIDs <- as.character(meta[,'Sample'])
  sampleIDs <- gsub('HTA', '', sampleIDs)
  sampleIDs <- gsub('_', ' ', sampleIDs)
  anno.sIDs <- HeatmapAnnotation(cn = anno_text(sampleIDs, rot = 50))
  
  # Patient IDs with biopsy pairs
  pIDpairs <- as.character(meta[,'BiopsyChange.Drug'])
  anno.pIDpairs <- HeatmapAnnotation(cn = anno_text(pIDpairs, rot = 50))
  
  # Return annotations
  return(list('pointer' = pointer, 
              'anno.sIDs' = anno.sIDs, 
              'anno.pIDpairs' = anno.pIDpairs))
  
}

# Convenience function to create all cohort heatmap annotations (except assay availability annotations)
make_heatmap_annotations <- function(meta) {
  
  # Treatment legend that also indicates on-progression biopsies
  # New approach to avoid splitting bug
  onProg.idx <- rep(8, nrow(meta))
  #onProg.idx[which(meta$ProgressionStage == 'OnProgression')] <- 8
  onProg.col <- colors.treatment[as.character(meta$Treatment)]
  onProg.col[which(meta$ProgressionStage == 'OnProgression')] <- 'black'
  onProgAnno <- HeatmapAnnotation('CDK4/6i' = anno_simple(meta$CDKi, pch = onProg.idx, 
                                                          pt_gp = gpar(col = onProg.col),
                                                          col = colors.treatment, border = TRUE), 
                                  show_legend = FALSE, border = TRUE)
  
  # Make legend for ERi therapy
  erAnno <- HeatmapAnnotation(ERi = meta$ERi, show_legend = FALSE, border = TRUE, 
                              col = list(ERi = colors.treatment))
  
  # Make response annotation object and legend list
  responseAnno <- HeatmapAnnotation(PatientResponse = meta$PatientResponse, show_legend = FALSE, 
                                    border = TRUE, col = list(PatientResponse = colors.response))
  
  # PAM50 annotations
  pamAnno <- HeatmapAnnotation(PAM50 = meta$PAM50, show_legend = FALSE, na_col = 'grey',
                               col = list(PAM50 = colors.pam), border = TRUE)
  
  # PAM50 change annotations
  pamChangeAnno <- HeatmapAnnotation(PAM50 = meta$pamChange, show_legend = FALSE, na_col = 'grey', 
                                     col = list(PAM50 = colors.pam), border = TRUE)
  
  # Color code biopsies by patient
  patientAnno <- HeatmapAnnotation(PatientID = meta$Patient, show_legend = FALSE, border = TRUE, 
                                   col = list(PatientID = c(colors.patient)))
  
  # Metastatic site annotations
  siteAnno <- HeatmapAnnotation(BiopsySite = meta$BiopsySite, show_legend = FALSE, na_col = 'grey', 
                                border = TRUE, col = list(BiopsySite = colors.site))
  
  # Make custom column name annotations for working with HTAN patients
  colAnno <- make_heatmap_columnIDs(meta)
  pointers <- colAnno[['pointer']]
  anno.sIDs <- colAnno[['anno.sIDs']]
  anno.pIDpairs <- colAnno[['anno.pIDpairs']]
  
  # Return list of annotation objects  
  return(list('onProgAnno' = onProgAnno, 
              'erAnno' = erAnno,
              'responseAnno' = responseAnno, # change to "responseAnno"
              'pamAnno' = pamAnno, 
              'pamChangeAnno' = pamChangeAnno,
              'siteAnno' <- siteAnno,
              'patientAnno' = patientAnno,
              'htanPointer' = pointers,
              'sampleAnno' = anno.sIDs,
              'biopPairAnno' = anno.pIDpairs))
        
}

# Function to make color row annotations for heatmaps
make_row_annotations <- function(mat, category_table, category_colors, 
                                 lgd_rows = 2, fntSize = 10, gridSize = 4) {
  
  # Only use features present in matrix
  category_table <- category_table[category_table[,c(1)] %in% rownames(mat),]
  
  # Set category for features in multiple categories to NA
  multiRows <-  names(table(category_table[,c(1)])[table(category_table[,c(1)])>1])
  category_table[category_table[,c(1)] %in% multiRows,c(2)] <- NA
  category_table <- distinct(category_table)
  
  # Match feature category table with matrix
  rownames(category_table) <- as.character(category_table[,c(1)])
  category_table <- category_table[rownames(mat),]
  
  # Set any categories without available color annotation assignments to NA
  category_table[!(category_table[,c(2)] %in% names(category_colors)),c(2)] <- NA
  
  # Create row annotation and legend if any remaining features can be color annotated
  if (sum(!(is.na(category_table[,2]))) > 0) {
    
    row_anno <- rowAnnotation(Category = category_table[,c(2)], col = list(Category = category_colors), 
                              show_legend = FALSE, na_col = 'white', show_annotation_name = FALSE)
    
    row_lgd <- list(Legend(labels = names(category_colors)[names(category_colors) %in% category_table[,c(2)]], 
                           legend_gp = gpar(fill = category_colors[names(category_colors) %in% category_table[,c(2)]]),
                           title_gp = gpar(fontsize = fntSize, fontface = "bold"), labels_gp = gpar(fontsize = fntSize), 
                           title = 'Category', nrow = lgd_rows, grid_height = unit(gridSize, 'mm'), grid_width = unit(gridSize, 'mm')))
  
  } else {
    
    row_anno = NULL
    row_lgd = NULL
    
  }
  
  return(list(row_anno, row_lgd))
  
}

# Convenience function to clean up names when splitting heatmaps by category
format_cats <- function(categories) {
  
  categories <- gsub('_', '\n', categories)
  categories <- gsub(' ', '\n', categories)
  categories <- gsub('-', '\n', categories)
  
  return(categories)
  
}


# Function to create heatmap annotation based on results of cutree 
#   Currently only set up for column clustering
make_cutree_anno <- function(meta, 
                             mat = NULL, 
                             ct = NULL, 
                             compute_change = FALSE,
                             assign_to_pair = TRUE,
                             select_features = NULL,
                             select_samples = NULL,
                             dist_func = 'euclidean', 
                             clust_type = 'complete',
                             cluster_columns = TRUE,
                             cluster_name = 'Cluster',
                             cluster_colors = NULL,
                             cluster_values = c('High', 'Low'),
                             pirate_pallete = 'southpark',
                             lgd_rows = 2,
                             lgd_fontsize = 12, 
                             lgd_gridsize = 4,
                             k = 2) {
  
  # Select samples
  if (is.null(select_samples)) {
    
    select_samples <- colnames(mat)
    
  }
  
  # Select features
  if (is.null(select_features)) {
    
    select_features <- rownames(mat)
    
  }
  
  # Compute delta between paired samples
  if (compute_change) {
    
    # Subset matrix to select samples
    mat <- mat[,select_samples,drop=FALSE]
    
    # Compute change
    mat.change <- compute_paired_change(mat, meta)
    
    # Update sample list to 2nd of all sample pairs
    select_samples <- select_samples[select_samples %in% colnames(mat.change)]
    
    # Update matrix for heatmap
    mat <- mat.change[rownames(mat),select_samples,drop=FALSE]
    
  } else {
    
    # Subset matrix to select samples
    mat <- mat[,select_samples,drop=FALSE]
    
  }
  
  # Run clustering and cut tree if no clustering assignment already provided
  if (is.null(ct)) {
    
    ct <- cut_dendro(mat, 
                     select_features = select_features, 
                     select_samples = select_samples,
                     dist_func = dist_func, 
                     clust_type = clust_type,
                     cluster_columns = cluster_columns,
                     k = k)
    
  }
  
  
  # Merge clustering assignments with meta table to build annotation
  meta.rownames <- rownames(meta)
  
  # Set name of cluster annotation for meta table
  ct <- ct %>%
    as.matrix() %>%
    melt() %>%
    select(Var1, value) %>%
    setNames(c('Sample', cluster_name)) 
  
  # Swap cluster assignment names
  if (!is.null(cluster_values) & length(cluster_values) == k) {
    
    # If using high/low scheme, assign names by cluster mean
    if (k == 2 & 'High' %in% cluster_values & 'Low' %in% cluster_values) {
      
      ct <- mat %>%
        as.matrix() %>%
        melt() %>%
        setNames(c('Feature', 'Sample', 'Value')) %>%
        filter(Feature %in% select_features) %>%
        left_join(ct) %>%
        group_by(.data[[cluster_name]]) %>%
        mutate(Mean = mean(Value)) %>%
        ungroup() %>%
        mutate(High = max(Mean)) %>%
        select(Sample, .data[[cluster_name]], Mean, High) %>%
        distinct() %>%
        mutate(!!cluster_name := case_when(Mean == High ~ 'High',
                                           TRUE ~ 'Low')) %>%
        select(Sample, .data[[cluster_name]]) %>%
        distinct()
      
    } else {
      
      names(cluster_values) <- 1:k
      
      ct <- ct %>%
        mutate(!!cluster_name := cluster_values[.data[[cluster_name]]])
      
    }
    
  }
  
  # Remove from meta if already present from previous run
  if (cluster_name %in% colnames(meta)) {
    
    meta <- meta %>% select(-c(cluster_name))
    
  }
  
  # Merge with meta table  
  meta <- ct %>%
    right_join(meta) %>%
    data.frame(check.names = FALSE)
  rownames(meta) <- meta$Sample
  meta <- meta[meta.rownames,]
  
  
  # Assign cluster assignments to NA values in meta table when NA comes 
  #   from first of two samples when computing delta between pairs
  if (compute_change & assign_to_pair) {
    
    meta <- meta %>%
      group_by(Patient) %>%
      mutate(!!cluster_name := case_when(is.na(.data[[cluster_name]]) & !is.na(lead(.data[[cluster_name]])) ~ lead(.data[[cluster_name]]),
                                         TRUE ~ .data[[cluster_name]])) %>%
      ungroup() %>%
      data.frame(check.names = FALSE)
    rownames(meta) <- meta$Sample
    meta <- meta[meta.rownames,]
    
  }
  
  
  print(cluster_colors)
  
  # Set annotation colors
  if (is.null(cluster_colors)) {
    
    cluster_colors <- set_colors(meta[,cluster_name], pirate_pallete)
    
  }
  

  
  # Build heatmap annotation
  ctAnno = HeatmapAnnotation(Cluster = meta[,cluster_name], 
                             show_legend = FALSE, 
                             na_col = 'grey',
                             col = list(Cluster = cluster_colors), 
                             border = TRUE)
  names(ctAnno) <- cluster_name
  
  # Build legend
  ctLgd <- Legend(title = cluster_name, labels = names(cluster_colors), nrow = lgd_rows,
                  legend_gp = gpar(fill = cluster_colors),  labels_gp = gpar(fontsize = lgd_fontsize), 
                  grid_height = unit(lgd_gridsize, 'mm'), grid_width = unit(lgd_gridsize, 'mm'), 
                  title_gp = gpar(fontsize = lgd_fontsize, fontface = "bold"))
  
  
  return(list(ctAnno = ctAnno,
              ctLgd = ctLgd,
              meta = meta))
  
}


# Convienence function to add new heatmap column annotation without recreating list
# TODO: Add option to add new annotation at specified position
add_ht_annotation <- function(anno, anno_list, pos = 'last', swap = TRUE) {
  
  # Get name of new annotation
  anno.name <- names(anno)
  
  # Separate annotations and legends
  annos <- anno_list[[1]]
  lgds <- anno_list[[2]]
  
  # Get names of current annotation lsit
  annos.names <- names(annos)
  
  # Get positions in annotation list
  idx <- c(1:length(annos))
  names(idx) <- annos.names
  
  # Swap if already present
  if (swap & (anno.name %in% annos.names)) {
    
    # Get current position
    to_swap <- which(names(idx) == anno.name)
    
    # Remove from current list
    annos <- annos[-c(to_swap)]
    
  }
  
  # Append new annotation to list
  annos <- annos %v% anno
  
  return(list(annos, lgds))
  
}



# Function that organizes list of legends in a set order for plotting with heatmaps
make_legend_list <- function(row_lgd_list, top_lgd_list, ht_lgd, btm_lgd_list) {
  
  # Row annotation legend goes on left followed by heatmap legend
  if (!(is.null(row_lgd_list))) {
    
    lgd_list <- row_lgd_list
    l <- length(lgd_list)
    lgd_list[[l + 1]] <- ht_lgd
    
    # Then add top legends to list
    l <- length(lgd_list)
    if (!(is.null(top_lgd_list))) {
      
      for (i in 1:length(top_lgd_list)) {
        
        lgd_list[[l + i]] <- top_lgd_list[[i]]
        
      }  
      
    }
    
    # Then add bottom legends to list
    l <- length(lgd_list)
    if (!(is.null(btm_lgd_list))) {
      
      for (i in 1:length(btm_lgd_list)) {
        
        lgd_list[[l + i]] <- btm_lgd_list[[i]]
        
      }  
      
    }
    
  } 
  
  # If no row annotation, top annotation goes on left followed by heatmap legend
  else { 
    
    if (!(is.null(top_lgd_list))) {
      
      lgd_list <- top_lgd_list
      l <- length(lgd_list)
      lgd_list[[l + 1]] <- ht_lgd
      
    } 
    
    # If not top annotation, start with heatmap legend
    else {
      
      lgd_list <- list(ht_lgd)
      
    } 
    
    # Then add bottom legends to list
    l <- length(lgd_list)
    if (!(is.null(btm_lgd_list))) {
      
      for (i in 1:length(btm_lgd_list)) {
        
        lgd_list[[l + i]] <- btm_lgd_list[[i]]
        
      }  
      
    }
    
  }
  
  return(lgd_list)
  
}




# Function to create box around column of samples to highlight in heatmap
box_samples <- function(mat, samples, htName, dend, box_col = 'white', 
                        order_columns = NULL,
                        box_width = 1, merge_overlaps = FALSE) {
  
  # Find columns on heatmap
  if (class(dend) == 'dendrogram') {
    idx <- which(order.dendrogram(dend) %in% which(colnames(mat) %in% samples))
  } else if (!is.null(order_columns)) {
    idx <- which(order_columns %in% samples)
  } else {
    idx <- which(colnames(mat) %in% samples)
  }
  print(idx)
  rBound <- idx / ncol(mat)
  lBound <- rBound - (1/ncol(mat))
  
  # Merge overlapping boundaries
  if (merge_overlaps) {
    
  df <- round(data.frame(lBound, rBound), digits = 6) %>% 
    mutate(indx = c(0, cumsum(as.numeric(lead(lBound)) > cummax(as.numeric(rBound)))[-n()])) %>%
    group_by(indx) %>% 
    summarise(lBound = min(lBound), rBound = max(rBound)) %>% 
    data.frame()
  lBound <- df$lBound
  rBound <- df$rBound
  
  }
  
  # Overlay on top of heatmap
  decorate_heatmap_body(htName, { 
    
    for (i in 1:length(lBound)) {  
      
      # Vertical lines (between cells)
      grid.lines(c(lBound[i], lBound[i]), c(0, 1), gp = gpar(lty = 1, lwd = box_width, col = box_col))
      grid.lines(c(rBound[i], rBound[i]), c(0, 1), gp = gpar(lty = 1, lwd = box_width, col = box_col))
      
      # Horizontal lines (top and bottom of heatmap)
      grid.lines(c(lBound[i], rBound[i]), c(0, 0), gp = gpar(lty = 1, lwd = 0.5, col = box_col))
      grid.lines(c(lBound[i], rBound[i]), c(1, 1), gp = gpar(lty = 1, lwd = 0.5, col = box_col))
      
    }
    
  })
  
}

# Function to capture size of final heatmap and all of its components
#   to ensure proper image size during saving
calc_ht_size = function(ht, unit = "inch") {
  
  w = ComplexHeatmap:::width(ht)
  w = convertX(w, unit, valueOnly = TRUE)
  h = ComplexHeatmap:::height(ht)
  h = convertY(h, unit, valueOnly = TRUE)
  
  return(c(w, h))

}

# TODO: Rename function; also create separate parameter for legend instead of passing list
save_htan_heatmap <- function(ht_objects, fn, ht_gap = unit(4, "mm"), res = NULL, 
                              pointsize = 12, lgd_direction = 'horizontal', 
                              max_width = NULL, add_anno_title = NULL, 
                              add_width = 0, extend_w = 1.5, res_factor = 2,
                              add_height = 0, lgd_gap = unit(2, 'mm'),
                              mark_samples = NULL, mat = NULL, order_columns = NULL,
                              col_dend = NULL, box_col = NULL, box_width = NULL) {
  
  # Pull out heatmap and legend objects
  ht <- ht_objects[[1]]
  lgd <- ht_objects[[2]]
  
  # If no heatmap provided, just make empty one (used if just saving legend to png)
  if (is.null(ht)) {ht <- Heatmap(matrix(0, nrow = 0, ncol = 1))}
  
  # Pack list of legends
  pd <- packLegend(list = lgd, direction = lgd_direction,
                   max_width = max_width, gap = lgd_gap)
  
  # Draw dummy plot to determine figure size
  pdf(NULL)
  dht <- draw(ht, heatmap_legend_side = 'bottom', 
              annotation_legend_side = 'bottom', 
              annotation_legend_list = pd, ht_gap = ht_gap)
  wh <- calc_ht_size(dht) # wh[1] = width, wh[2] = height
  try(dev.off(), silent = TRUE)
  
  # Set width, height, and resolution
  w <- wh[1] + extend_w
  h <- wh[2]
  if (is.null(res)) {res <- (w*h)*res_factor}
  
  # Extend width and height if needed
  w <- w + add_width
  h <- h + add_height
  
  # Save as png
  png(filename = fn, width = w, height = h, units = 'in', 
      res = res, pointsize = pointsize)
  draw(ht, annotation_legend_list = pd, ht_gap = ht_gap, 
       heatmap_legend_side = 'bottom', 
       annotation_legend_side = 'bottom')
  
  # Create box around select samples if specified
  if (!(is.null(mark_samples))) {
    
    box_samples(mat, mark_samples, 'heat', col_dend, 
                order_columns = order_columns,
                box_col = box_col, box_width = box_width)
    
  }
  try(dev.off(), silent = TRUE)
  
  
  #dev.off()
  
}

# Function for creating dendrogram object used for clustering heatmaps
create_dendrogram <- function(mat, dist_func = 'euclidean', 
                              clust_type = 'complete', 
                              reorder_dend = TRUE) {
  
  # Only cluster if at least 2 rows, else return FALSE 
  if (nrow(mat) >= 2) {
    
    # Cluster rows/columns
    dend <- as.dendrogram(hclust(dist(mat, method = dist_func), method = clust_type))
    
    # Reorder dendrogram by row/column means
    if (reorder_dend) {
      
      dOrder = -rowMeans(mat, na.rm = T)
      dend <- reorder(dend, dOrder, agglo.FUN = mean)
      
    }
    
  } else {
    
    dend <- FALSE
    
  }
  
  return(dend)
  
}

# Function to create color scheme object for heatmap
make_heatmap_colors <- function(mat, htColors = NULL, numColors = 3, minHt = NULL, maxHt = NULL) {
  
  # Get color scheme and number of intensity points
  if (is.null(htColors)) {
    
    htColors = list(c("navy", "yellow"), c("cyan", "black", "magenta"), 
                    c('navy', 'darkseagreen1', 'goldenrod1', 'red4'))[[numColors-1]]
    
  } else {
    
    numColors <- length(htColors)
    
  }
  
  # Get min/max intesnity limits
  if (is.null(minHt)) {
    
    minHt <- max(abs(quantile(mat, 0.01, na.rm = TRUE)), 
                 quantile(mat, 0.99, na.rm = TRUE)) * (-1)
    
  } 
  
  if (is.null(maxHt)) {
    
    maxHt <- max(abs(quantile(mat, 0.01, na.rm = TRUE)), 
                 quantile(mat, 0.99, na.rm = TRUE))
    
  }
  
  # Set heatmap legend limits (ensure 0 in middle)
  if (numColors == 2) {
    
    htLims <- c(minHt, maxHt)
    
  } else if (numColors == 3) {
    
    htLims <- c(minHt, 0, maxHt)
    
  } else {
    
    htLims <- seq(minHt, maxHt, length = numColors)
    
  }
  
  # Build heatmap color function and return
  return(colorRamp2(htLims, htColors))
  
}


# Function to compute change across pairs of patient biopsies 
#   Needs meta with date specifying biopsy timepoint
# TODO: May store in separate library script if more space needed
compute_paired_change <- function(df, meta, return_matrix = TRUE, 
                                  select_features = NULL) {
  
  # Melt to long format
  df <- df %>% 
    as.matrix() %>%
    melt() %>%
    setNames(c('Feature', 'Sample', 'Value')) 
  
  # Subset to select features
  if (!is.null(select_features)) {
    
    df <- df %>% filter(Feature %in% select_features)
    
  }
  
  # Merge with meta table to get pairs and compute change
  df <- df %>%
    left_join(meta) %>%
    select(Patient, Sample, Date, Feature, Value) %>% 
    group_by(Patient, Feature) %>% 
    arrange(Date, by_group = TRUE) %>%
    filter(n() >= 2) %>%
    mutate(Change = Value - lag(Value)) %>%
    ungroup() %>%
    filter(!is.na(Change)) %>%
    select(Sample, Feature, Change)
  
  # Cast to matrix
  if (return_matrix) {
    
    df <- acast(df, Feature ~ Sample, value.var = 'Change')
    
  }
  
  return(df)
  
}  


# Main function for building omics heatmaps using ComplexHeatmap 
make_heatmap <- function(mat, 
                         meta, 
                         select_samples = NULL, 
                         select_features = NULL, 
                         
                         compute_change = FALSE,
                         
                         lgd_name = 'Activity', 
                         ht_cols = c('cyan', 'black', 'magenta'),
                         col_func = NULL,
                         show_heatmap_legend = TRUE,
                         lgdFntSize = 12, 
                         ht_lgd_length = 2,
                         lgd_gridsize = 4,
                         lgd_rows = 2, 
                         max_lgd_width = NULL, 
                         lgd_gap = unit(2, 'mm'),
                         
                         heatmap_width = unit(6, 'in'), add_width = 0, 
                         heatmap_height = unit(8, 'in'), add_height = 0,
                         cell_width = NULL, cell_height = NULL,
                         res = NULL, rect_gp = gpar(col = NA), 
                         
                         top_anno = NULL, btm_anno = NULL,
                         
                         top_anno.ht = NULL,
                         
                         show_column_annotation_legend = TRUE,
                         cluster_columns = TRUE, 
                         split_column_by_pheno = NULL, 
                         column_split_order = NULL,
                         split_column_by_dendrogram = NULL,
                         order_columns = NULL,
                         order_columns_by = NULL,
                         
                         cluster_rows = FALSE,
                         keep_row_order = TRUE,
                         category_table = NULL, 
                         split_by_cat = TRUE, 
                         cat_order = NULL, 
                         category_colors = NULL, annotate_categories = TRUE,
                         show_row_annotation_legend = TRUE,
                         show_row_names = NULL,
                         row_title_rot = 0,
                         
                         bar_anno = NULL,
                         bar_labs = seq(-0.5, 1, .5), 
                         bar_ylim = c(-.8, 1.25),
                         
                         
                         mark_samples = NULL,
                         box_col = 'white', 
                         box_width = 1,
                         
                         fn = NULL, 
                         return_heat_objects = FALSE) {
  
  

  # ~~~~~~~~~~~~~~~~~~~ ROW PARAMETERS ~~~~~~~~~~~~~~~~~~~~~ #
  
  # Select features to plot from matrix
  if (is.null(select_features)) {
    
    select_features <- rownames(mat)
    
  } else {
    
    select_features <- select_features[select_features %in% intersect(select_features, rownames(mat))]
    
  }
  
  # Sort unless keeping original order
  if (!keep_row_order) {select_features <- sort(select_features)}
  
  # Subset matrix to select features
  mat <- mat[select_features,,drop = FALSE]
  
  # Hide row names if too many features
  if (is.null(show_row_names)) {
    
    if (nrow(mat) > 70) {
      
      show_row_names <- FALSE
      
    } else {
      
      show_row_names <- TRUE
      
    }
    
  }
  
  # Set parameters for row splitting and/or annotating
  if (!is.null(category_table)) {
    
    # Subset category table down to select features
    select_features <- select_features[select_features %in% category_table[[1]]]
    #category_table <- category_table[category_table[[1]] %in% select_features,]
    category_table <- category_table[select_features,]
    mat <- mat[category_table[[1]],,drop = FALSE]
    
    # Color annotate categories
    if (annotate_categories) {
      
      # Make row annotation object and legend
      rowAnno.items <- make_row_annotations(mat, category_table = category_table, 
                                            fntSize = lgdFntSize, gridSize = lgd_gridsize,
                                            category_colors = category_colors, lgd_rows = lgd_rows)
      
      rowAnno <- rowAnno.items[[1]]
      rowLgd <- rowAnno.items[[2]]
      
    } else {
      
      rowAnno <- NULL
      rowLgd <- NULL
      
    }
    
    # Split rows by category
    if (split_by_cat) {
      
      # Format gene/protein categories
      category_table <- category_table %>%
        mutate(Category = format_cats(Category)) 
      
      # Set specific order if provided (reformat first)
      if (!is.null(cat_order)) {
        
        category_table <- category_table %>%
          mutate(Category = factor(Category, levels = format_cats(cat_order)))
        
      } 
      
      # Pull out row split for heatmap
      row_split <- category_table %>% pull(Category)
      
    } else {
      
      row_split <- NULL
    
    }
    
  } else {
    
    row_split <- NULL
    rowAnno <- NULL
    rowLgd <- NULL
  
  }
  
 
  # Set heatmap height based on cell height if provided
  if (!is.null(cell_height)) {heatmap_height <- cell_height*nrow(mat)}
  
  # ~~~~~~~~~~~~~~~~~~~ COLUMN PARAMETERS ~~~~~~~~~~~~~~~~~~~~~ #
  
  # Select samples for heatmap
  if (is.null(select_samples)) {
    
    select_samples <- colnames(mat)
    
  } else {
    
    select_samples <- intersect(select_samples, colnames(mat))
    
  }
  
  # Compute delta between paired samples
  if (compute_change) {
    
    # Subset matrix to select samples
    mat <- mat[,select_samples,drop=FALSE]
    
    # Compute change
    mat.change <- compute_paired_change(mat, meta)
    
    # Update sample list to 2nd of all sample pairs
    select_samples <- select_samples[select_samples %in% colnames(mat.change)]
    
    # Update matrix for heatmap
    mat <- mat.change[rownames(mat),select_samples,drop=FALSE]
    
  } else {
    
    # Subset matrix to select samples
    mat <- mat[,select_samples,drop=FALSE]
    
  }
  
  # Add barplot annotation above column annotations
  if (!is.null(bar_anno)) {
    
    if (!is.null(top_anno)) {
      
      top_anno[[1]] <- make_barplot_annotation(meta, bar_anno,
                                               labs = bar_labs, ylim = bar_ylim) %v% top_anno[[1]]
      
    } else {
      
      top_anno <- list(make_barplot_annotation(meta, bar_anno,
                                               labs = bar_labs, ylim = bar_ylim), NULL)
      
    }
    
  }
  
  # Subset column annotations to select samples
  meta$AnnoIndex <- 1:nrow(meta)
  if (!is.null(top_anno)) {top_anno <- update_annotations(meta[select_samples,], top_anno)}
  if (!is.null(btm_anno)) {btm_anno <- update_annotations(meta[select_samples,], btm_anno)}
  
  
  # Subset annotation for direct to heatmp
  # TODO: may need to update to include legend as 2nd item in list
  #        (currently just a single anno object)
  if (!is.null(top_anno.ht)) {
    
    idx <- meta[select_samples,'AnnoIndex']
    top_anno.ht <- top_anno.ht[idx,]
    
    # Use heatmap function to perform clustering 
    # instead of creating dendrogram ahead of time
    cluster_columns.ht <- cluster_columns
    
  } else {
    
    cluster_columns.ht <- FALSE
    
  }
  
  
  # Split columns by specified phenotype
  if (!is.null(split_column_by_pheno)) {
    
    col_split <- meta[select_samples,split_column_by_pheno]
    #cluster_columns.ht <- cluster_columns
    cluster_columns <- FALSE
    
    # Get order of column groups
    if (is.null(column_split_order)) {
      
      column_split_order <- unique(col_split)
      
    }
    
    # Assign order as factors 
    col_split <- factor(col_split, column_split_order)
    
  } else {
    
    col_split <- NULL
    
  }
  
  # Split by dendrogram using number of clusters specified
  if (!is.null(split_column_by_dendrogram)) {
    
    col_split <- split_column_by_dendrogram
    column_title <- NULL
    cluster_columns <- TRUE
    if (!is.null(top_anno.ht)) {
      
      cluster_columns.ht <- TRUE
      
    }
    
 
  } else {
    
    col_split <- col_split
    column_title <- unique(sort(col_split))
  
  }
  
  # Cluster columns if specified
  if (cluster_columns) {cluster_columns <- create_dendrogram(t(mat))}
  
  # Set heatmap width based on cell width if provided
  if (!is.null(cell_width)) {heatmap_width <- cell_width*ncol(mat)}
  
  
  
  # Set column order
  if (!is.null(order_columns_by)) {
    
    order_columns <- names(sort(mat[order_columns_by,], decreasing = TRUE, na.last = TRUE))
    
  }
  
  
  # Update parameters when passing top annotation directly into heamtap function
  col_split.ht <- NULL
  if (!is.null(top_anno.ht)) {
    
    if (!is.null(col_split)) {
      
      col_split.ht <- col_split
      
    } 

  } 
  

  # Create top heatmap object to control column clustering and splitting
  mat.top <- matrix(0, nrow = 0, ncol = ncol(mat))
  colnames(mat.top) <- colnames(mat)
  ht.top <- Heatmap(mat.top,
                    cluster_columns = cluster_columns,
                    column_order = order_columns,
                    column_split = col_split,
                    column_title = column_title,
                    heatmap_width = heatmap_width)
  
  
  # ~~~~~~~~~~~~~~ BUILD HEATMAP OBJECTS ~~~~~~~~~~~~~~~~ #
  
  # Create heatmap color function
  if (is.null(col_func)) {
    
    col_func <- make_heatmap_colors(mat, htColors = ht_cols)
    
  }
  
  # Create heatmap legend
  ht_lgd = Legend(title = lgd_name, direction = "horizontal", 
                  legend_width = unit(ht_lgd_length, 'in'), 
                  legend_height = unit(ht_lgd_length, 'in'), 
                  title_gp = gpar(fontsize = lgdFntSize, fontface = "bold"), 
                  labels_gp = gpar(fontsize = lgdFntSize), 
                  title_position = "topcenter", col_fun = col_func)
  
  # Main heatmap body
  ht <- Heatmap(mat, col = col_func, show_heatmap_legend = FALSE, 
                name = 'heat',
                show_row_names = show_row_names,
                left_annotation = rowAnno, rect_gp = rect_gp,
                row_split = row_split, 
                row_title_rot = row_title_rot,
                cluster_rows = cluster_rows, 
                
                # TEST
                #cluster_columns = FALSE, 
                #cluster_columns = TRUE,
                #column_split = col_split,
                #column_title = NULL,
                #top_annotation = top_anno.ht,
                
                cluster_columns = cluster_columns.ht,
                column_split = col_split.ht,
                column_title = NULL,
                top_annotation = top_anno.ht,
                
                
                width = heatmap_width, height = heatmap_height)
  
  # Stack heatmap objects
  if (is.null(top_anno.ht)) {
  
  ht <- ht.top %v% 
    top_anno[[1]] %v% 
    ht %v% 
    btm_anno[[1]]
  
  } else {
    
    ht <- ht %v% 
      btm_anno[[1]]
    
  }
  
  
 
  
  # ~~~~~~~~~~~~~~ BUILD LEGEND OBJECTS ~~~~~~~~~~~~~~~~ #
  
  # Remove legends not to be included
  if (!show_heatmap_legend) {ht_lgd <- NULL}
  if (!show_row_annotation_legend) {rowLgd <- NULL}
  if (!show_column_annotation_legend) {top_anno <- NULL}
  if (!show_column_annotation_legend) {btm_anno <- NULL}
  
  # Pack legends into list
  lgd_list <- make_legend_list(row_lgd_list = rowLgd, 
                               top_lgd_list = top_anno[[2]], 
                               ht_lgd = ht_lgd, 
                               btm_lgd_list = btm_anno[[2]])
  
  # Save heatmap to file
  if (!is.null(fn)) {
    
    if (is.null(max_lgd_width)) {max_lgd_width <- heatmap_width}
    
    save_htan_heatmap(list(ht, lgd_list), fn, ht_gap = unit(1, "mm"), 
                      add_height = add_height, add_width = add_width,
                      res = res, max_width = max_lgd_width, lgd_gap = lgd_gap,
                      mark_samples = mark_samples, order_columns = order_columns,
                      mat = mat, col_dend = cluster_columns,
                      box_width = box_width, box_col = box_col)
    
  }
  
  # Return heatmap and legend objects
  if (return_heat_objects) { return(list(ht = ht, lgd = lgd_list))  }
  
}







