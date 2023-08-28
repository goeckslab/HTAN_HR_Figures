

######################################################
#
#   HEATMAP FUNCTIONS FOR HTAN HR+ MANUSCRIPT
#
######################################################

library(ComplexHeatmap)


# Annotation colors
colors.treatment <- c("Abemaciclib" = "#EE0011FF", 
                      "Palbociclib" = "#0C5BB0FF", 
                      "Everolimus" = "#7F8624FF", 
                      "Fulvestrant" = "#FA6B09FF", 
                      "Letrozole" = "#15983DFF", 
                      "Tamoxifen" = "maroon", 
                      "Post-Abemaciclib" = "#FFC6C7FF", 
                      "Post-Palbociclib" = "#B7CBFFFF", 
                      "Post-Everolimus" = "#EBF09C", 
                      "Post-Fulvestrant" = "#FFE7D4", 
                      "Post-Letrozole" = "#8EFFBA", 
                      "Post-Tamoxifen" = "#FFD3E5",
                      'Pre-Palbociclib' = 'white',
                      'Pre-Abemaciclib' = 'white',
                      'Pre-Fulvestrant' = 'white',
                      'Pre-Letrozole' = 'white')
                      
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


# Function for building heatmap annotation legend objects
make_heatmap_legends <- function(meta, select_samples = NULL, lgd_rows = 2,
                                 treatment_title = ' ', strip_names = TRUE, rm_duplicate_pair = TRUE,
                                 lgd_fontsize = 12, lgd_gridsize = 4) {
  
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

# Function for creating barplot column annotations for heatmaps
make_barplot_annotation <- function(meta, anno, anno_height = unit(1, 'in'), 
                                    labs = seq(-0.5, 1, .5), ylim = c(-.8, 1.25)) {
  
  # Compute axis parameters if not provided
  if (is.null(labs) | is.null(ylim)) {
    
    m.ax <- ceiling(max(meta[,anno], na.rm = TRUE)*2) / 2
    m.in <- min(0, (floor(min(meta[,anno], na.rm = TRUE)*2) / 2))
    
    if (is.null(labs)) {
      
      labs <- seq(m.in,m.ax,0.5)
      
    }
    
    if (is.null(ylim)) {
      
      ylim <- c(m.in,m.ax)
      
    }
    
  }
  
  
  # Create barplot anotation
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
    
    bar_anno <- NULL
    
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
  onProg.idx <- rep(NA, ncol(meta))
  onProg.idx[which(meta$ProgressionStage == 'OnProgression')] <- 8
  onProgAnno <- HeatmapAnnotation('CDK4/6i' = anno_simple(meta$CDKi, pch = onProg.idx, 
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
  pointers <- colAnno[['pointers']]
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


# TODO: Rename function; also create separate parameter for legend instead of passing list
save_htan_heatmap <- function(ht_objects, fn, ht_gap = unit(4, "mm"),  res = NULL, 
                              pointsize = 12, lgd_direction = 'horizontal', 
                              max_width = NULL, add_anno_title = NULL, add_width = 0, 
                              add_height = 0, lgd_gap = unit(2, 'mm')) {
  
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
  w <- wh[1] + 1.5
  h <- wh[2]
  if (is.null(res)) {res <- (w*h)*2}
  
  # Extend width and height if needed
  w <- w + add_width
  h <- h + add_height
  
  # Save as png
  png(filename = fn, width = w, height = h, units = 'in', 
      res = res, pointsize = pointsize)
  draw(ht, annotation_legend_list = pd, ht_gap = ht_gap, 
       heatmap_legend_side = 'bottom', 
       annotation_legend_side = 'bottom')
  dev.off()
  
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




# Main function for building omics heatmaps using ComplexHeatmap 
make_heatmap <- function(mat, 
                         meta, 
                         select_samples = NULL, 
                         select_features = NULL, 
                         
                         compute_change = FALSE,
                         
                         lgd_name = 'Activity', 
                         ht_cols = c('cyan', 'black', 'magenta'),
                         force_col_fun = NULL,
                         show_heatmap_legend = TRUE,
                         lgdFntSize = 12, 
                         ht_lgd_length = 2,
                         lgd_gridsize = 4,
                         max_lgd_width = NULL, 
                         
                         heatmap_width = unit(6, 'in'), add_width = 0, 
                         heatmap_height = unit(8, 'in'), add_height = 0,
                         res = NULL, rect_gp = gpar(col = NA), 
                         
                         top_anno = NULL, btm_anno = NULL,
                         show_column_annotation_legend = TRUE,
                         cluster_columns = FALSE, 
                         split_column_by_pheno = NULL, 
                         split_column_by_dendrogram = NULL,
                         
                         cluster_rows = FALSE,
                         category_table = NULL, 
                         split_by_cat = TRUE, 
                         cat_order = NULL, 
                         category_colors = NULL, annotate_categories = TRUE,
                         show_row_annotation_legend = TRUE,
                         show_row_names = NULL,
                         row_title_rot = 0,
                         
                         compute_cyto = NULL,
                         g1_score = NULL,
                         
                         fn = NULL, 
                         return_heat_objects = FALSE) {
  
  
  # Subset down to select features
  if (is.null(select_features)) {select_features <- rownames(mat)}
  mat <- mat[select_features,,drop = FALSE ]
  
  
  # Subset category table down to select features
  if (!is.null(category_table)) {
    
    category_table <- category_table[category_table[[1]] %in% c(select_features,rownames(mat)),]
    mat <- mat[category_table[[1]],,drop = FALSE]
    
    # Color annotate categories
    if (annotate_categories) {
      rowAnno.items <- make_row_annotations(mat, category_table = category_table, 
                                            fntSize = lgdFntSize, gridSize = lgd_gridsize,
                                            category_colors = category_colors, lgd_rows = 2)
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
        mutate(Category = gsub('_', '\n', Category)) %>%
        mutate(Category = gsub(' ', '\n', Category)) %>%
        mutate(Category = gsub('-', '\n', Category))
      
      row_split <- category_table[[2]]
      
      # Set specific order if provided
      if (!is.null(cat_order)) {
        
        cat_order <- gsub('_', '\n', cat_order)
        cat_order <- gsub(' ', '\n', cat_order)
        cat_order <- gsub('-', '\n', cat_order)
        
        row_split <- factor(row_split, levels = cat_order)
        
      }
      
    } else {
      row_split <- NULL
    }
    
  } else {
    row_split <- NULL
    rowAnno <- NULL
    rowLgd <- NULL
  }
  
  
  
  
  
  
  
  
  # Select samples for heatmap
  if (is.null(select_samples)) {select_samples <- colnames(mat)}
  
  # Subset down to select samples
  mat <- mat[,select_samples,drop = FALSE]
  
  # Compute cytolytic activity
  # TODO: Need to compute cyto, then compute change, then create annotation object
  if (!is.null(compute_cyto)) {
    cyto_act <- get_cyto_activity(compute_cyto[,select_samples])
    #cyto_act <- get_cyto_activity(compute_cyto)
    cyto_anno <- cyto_annotation(cyto_act)
    
    
    
  } else {
    cyto_anno <- NULL
  }
  
  
  
  # Compute change between samples
  if (compute_change) {
    
    mat.change <- mat %>% melt() %>% 
      setNames(c('Feature', 'Sample', 'Value')) %>%
      left_join(meta) %>%
      select(Patient, Sample, Date, Feature, Value) %>% 
      group_by(Patient, Feature) %>% 
      arrange(Date, by_group = TRUE) %>%
      filter(n() >= 2) %>%
      mutate(Change = Value - lag(Value)) %>%
      ungroup() %>%
      filter(!is.na(Change)) 
    
    mat.change <- acast(mat.change, Feature ~ Sample, value.var = 'Change')
    
    select_samples <- select_samples[select_samples %in% colnames(mat.change)]
    
    # Create heatmap colors
    df.change <- calc_change(mat, meta[colnames(mat),], c("BxE", "Bx1", "Bx2", "Bx3", "Bx4", "Bx5", "Bx", "BxA"))[[1]]
    col_func <- make_heatmap_colors(df.change, htColors = ht_cols)
    #col_func <- make_heatmap_colors(mat.change, htColors = ht_cols)
    
    mat <- mat.change[rownames(mat),select_samples]
    
    # Get G1 score if provided
    if (!is.null(g1_score)) {
      
      g1_score <- g1_score %>%
        filter(Sample %in% colnames(mat),
               Signature == 'G1-arrest') %>%
        select(Sample, Score) %>%
        distinct() %>%
        data.frame()
      rownames(g1_score) <- g1_score$Sample
      g1_score <- g1_score[colnames(mat),]
      
      # Create barplot anotation
      g1_anno <- HeatmapAnnotation('G1 Arrest Score' = anno_barplot(g1_score$Score, 
                                                                    height = unit(1, "in"), 
                                                                    axis_param = list(
                                                                      at = seq(-0.5, 1, .5),
                                                                      labels = seq(-0.5, 1, .5)),
                                                                    ylim = c(-.8,1.25))) 
      
    }
    
    if (!is.null(compute_cyto)) {
      
      # Get cytolytic activity
      cyto_act <- get_cyto_activity(compute_cyto)
      cyto_act$Sample <- rownames(cyto_act)
      
      # Compute change across pairs
      cyto_act <- cyto_act %>%  
        left_join(meta) %>%
        select(Patient, Sample, Date, Score) %>% 
        group_by(Patient) %>% 
        arrange(Date, by_group = TRUE) %>%
        filter(n() >= 2) %>%
        mutate(Change = Score - lag(Score)) %>%
        ungroup() %>%
        filter(!is.na(Change)) %>%
        data.frame()
      
      rownames(cyto_act) <- cyto_act$Sample
      cyto_act <- cyto_act[colnames(mat),]
      
      # Create barplot annotation
      cyto_anno <- HeatmapAnnotation('Cytolytic Change' = anno_barplot(cyto_act$Change, 
                                                                       height = unit(1, "in"), 
                                                                       axis_param = list(#at = c(-1,-0.5,0,0.5,1,1.5), 
                                                                         at = seq(-0.5, 1, .5),
                                                                         #labels = c(-1,-0.5,0,0.5,1,1.5)), 
                                                                         labels = seq(-0.5, 1, .5)),
                                                                       ylim = c(-.8,1.25))) 
      
    } else {
      cyto_anno <- NULL
    }
    
  }
  
  
  
  
  # Subset annotations
  meta$AnnoIndex <- 1:nrow(meta)
  if (!is.null(top_anno)) {top_anno <- update_annotations(meta[select_samples,], top_anno)}
  if (!is.null(btm_anno)) {btm_anno <- update_annotations(meta[select_samples,], btm_anno)}
  
  if (!is.null(cyto_anno)) {
    top_anno[[1]] <- cyto_anno %v% top_anno[[1]]
  }
  
  # Add G1 arrest score (only if change computed)
  if (!is.null(g1_score)) {
    top_anno[[1]] <- g1_anno %v% top_anno[[1]]
  }
  
  
  # Split columns by specified phenotype
  if (!is.null(split_column_by_pheno)) {
    
    if (split_column_by_pheno == 'HTAN') {
      
      meta[,'HTAN'] <- gsub('_p2', '', meta[,'HTAN'])
      meta[,'HTAN'] <- gsub('HTA', '', meta[,'HTAN'])
      meta[,'HTAN'] <- factor(meta[,'HTAN'], levels = c("9-1", "9-2", "9-3", "9-14", "9-15"))
    
    }
    
    col_split <- meta[select_samples,split_column_by_pheno]
    cluster_columns <- FALSE
  
  } else {
    
    col_split <- NULL
    
  }
  
  
  # Split by dendrogram
  if (!is.null(split_column_by_dendrogram)) {
    
    col_split <- split_column_by_dendrogram
    cluster_columns <- TRUE
    column_title <- NULL
  
  } else {
    
    col_split <- col_split
    column_title <- unique(col_split)
  
  }
  
  
  
  
  
  
  
  
  # Use supplied heatmap color function 
  if (!is.null(force_col_fun)) { col_func <- force_col_fun }
  
  
  # Create heatmap colors
  # TODO: shouldn't this be swapped with above????
  col_func <- make_heatmap_colors(mat, htColors = ht_cols)
  
  # Create heatmap legend
  ht_lgd = Legend(title = lgd_name, direction = "horizontal", 
                  legend_width = unit(ht_lgd_length, 'in'), legend_height = unit(ht_lgd_length, 'in'), 
                  title_gp = gpar(fontsize = lgdFntSize, fontface = "bold"), labels_gp = gpar(fontsize = lgdFntSize), 
                  title_position = "topcenter", col_fun = col_func)
  
  
  
  
  
  # Hide row names
  if (is.null(show_row_names)) {
    
    if (nrow(mat) > 70) {
      
      show_row_names <- FALSE
      
    } else {
      
      show_row_names <- TRUE
      
    }
    
  }
  
  
  # Main heatmap body
  ht <- Heatmap(mat, col = col_func, show_heatmap_legend = FALSE, 
                left_annotation = rowAnno, row_split = row_split,
                row_title_rot = row_title_rot,
                cluster_rows = cluster_rows, cluster_columns = FALSE, 
                show_row_names = show_row_names,
                rect_gp = rect_gp,
                width = heatmap_width, height = heatmap_height)
  
  
  # Cluster columns if specified
  if (cluster_columns) {cluster_columns <- create_dendrogram(t(mat))}
  
  
  
  # Create empty matrix to store plot title
  ht.title <- Heatmap(matrix(0, nrow = 0, ncol = length(select_samples)), 
                      #column_order = order_columns, 
                      cluster_columns = cluster_columns,
                      column_split = col_split,
                      column_title = column_title,
                      #column_title = ht_title, 
                      heatmap_width = heatmap_width)
  
  # Stack heatmap objects
  ht <- ht.title %v% top_anno[[1]] %v% ht %v% btm_anno[[1]]
  
  
  
  
  
  
  
  
  
  # Merge legends
  # TODO: CAN WE JUST USE make_legend_list() function?
  lgd_list <- list()
  
  # Add row annotation legend
  if (!is.null(rowLgd) & show_row_annotation_legend) {
    n <- 1
    lgd_list[[n]] <- rowLgd[[1]]
  } else {
    n <- 0
  }
  
  # Add top annotations
  if (!is.null(top_anno) & show_column_annotation_legend) {
    for (i in 1:length(top_anno[[2]])) { lgd_list[[i+n]] <- top_anno[[2]][[i]]  }
  } else {
    i <- 0
  }
  
  # Add heatmap legends
  if (show_heatmap_legend) {
    n <- n + i + 1
    lgd_list[[n]] <- ht_lgd
    names(lgd_list)[n] <- 'ht_lgd'
  }
  
  # Add bottom annotations
  if (!is.null(btm_anno) & show_column_annotation_legend) {
    for (i in 1:length(btm_anno[[2]])) {  lgd_list[[i+n]] <- btm_anno[[2]][[i]]  }
  }
  
  
  
  
  
  # Save heatmap to file
  if (!is.null(fn)) {
    
    if (length(lgd_list) == 0) {lgd_list <- list(NULL)}
    
    if (is.null(max_lgd_width)) {max_lgd_width <- heatmap_width}
    
    save_htan_heatmap(list(ht, lgd_list), fn, ht_gap = unit(1, "mm"), 
                      add_height = add_height, add_width = add_width,
                      res = res, max_width = max_lgd_width)
    
  }
  
  
  
  
  # Return heatmap and legend objects
  if (return_heat_objects) { return(list(ht = ht, lgd = lgd_list))  }
  
}







