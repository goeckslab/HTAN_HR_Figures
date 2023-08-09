

######################################################
#
#   HEATMAP FUNCTIONS FOR HTAN HR+ MANUSCRIPT
#
######################################################

library(ComplexHeatmap)


# Function for building heatmap annotation legend objects
make_heatmap_legends <- function(samples_table, meta, select_samples = NULL, lgd_rows = 2,
                                 treatment_title = ' ', strip_names = TRUE, rm_duplicate_pair = TRUE,
                                 lgd_fontsize = 12, lgd_gridsize = 4) {
  
  # All color schemes subsetted down to select patients
  if (!is.null(select_samples)) {
    
    drug_cols <- drug_cols[names(drug_cols) %in% meta[select_samples, 'Treatment']]
    
    # TODO: Change to response colors
    intrinsic_cols <- intrinsic_cols[names(intrinsic_cols) %in% meta[select_samples, 'PatientResponse']]
    
    pam_colors <- pam_colors[names(pam_colors) %in% meta[select_samples, 'PAM50']]
    pam_change_colors <- pam_change_colors[names(pam_change_colors) %in% meta[select_samples, 'pamChange']]
    site_colors <- site_colors[names(site_colors) %in% meta[select_samples, 'BiopsySite']]
    
    htan_cols <- htan_cols[names(htan_cols) %in% meta[select_samples, 'HTAN']]
    
  }
  
  # Legend for pre-treatment and on-progression with on-progression labeled with asterisk
  opAstrLgd <- Legend(title = 'Treatment', type = 'points', pch = c(NA,8), at = c('PreTreatment', 'OnProgression'), background = c('white'),
                      labels_gp = gpar(fontsize = lgd_fontsize), grid_height = unit(lgd_gridsize, 'mm'), grid_width = unit(lgd_gridsize, 'mm'),
                      title_gp = gpar(fontsize = lgd_fontsize, fontface = "bold"))
  
  # Legend for CDK4/6i and ERi treatment
  # TODO: Either bring in function or isoloate here
  onProgLgd <- make_treatment_heatmap_annotations(meta, select_samples = select_samples, lgd_rows = lgd_rows, 
                                                  treatment_title = treatment_title, onProgression = TRUE)[[2]]
  
  # Make response annotation object and legend list
  # TODO: Change to responseLgd
  intrinsicLgd <- Legend(labels = names(intrinsic_cols), legend_gp = gpar(fill = intrinsic_cols), title = 'Patient Response', nrow = lgd_rows, 
                         labels_gp = gpar(fontsize = lgd_fontsize), grid_height = unit(lgd_gridsize, 'mm'), grid_width = unit(lgd_gridsize, 'mm'),
                         title_gp = gpar(fontsize = lgd_fontsize, fontface = "bold"))
  
  # PAM50 annotation legend
  pamLgd <- Legend(labels = names(pam_colors), legend_gp = gpar(fill = pam_colors), title = 'PAM50', nrow = lgd_rows, 
                   labels_gp = gpar(fontsize = lgd_fontsize), grid_height = unit(lgd_gridsize, 'mm'), grid_width = unit(lgd_gridsize, 'mm'),
                   title_gp = gpar(fontsize = lgd_fontsize, fontface = "bold"))
  
  # PAM50 change annotation legend
  pamChangeLgd <- Legend(labels = names(pam_change_colors), legend_gp = gpar(fill = pam_change_colors), title = 'PAM50', nrow = lgd_rows, 
                         labels_gp = gpar(fontsize = lgd_fontsize), grid_height = unit(lgd_gridsize, 'mm'), grid_width = unit(lgd_gridsize, 'mm'),
                         title_gp = gpar(fontsize = lgd_fontsize, fontface = "bold"))
  
  
  # Metastatic site annotations
  siteLgd <- Legend(labels = names(site_colors), legend_gp = gpar(fill = site_colors), title = 'Biopsy Site', nrow = lgd_rows, 
                    labels_gp = gpar(fontsize = lgd_fontsize), grid_height = unit(lgd_gridsize, 'mm'), grid_width = unit(lgd_gridsize, 'mm'),
                    title_gp = gpar(fontsize = lgd_fontsize, fontface = "bold"))
  
  
  
  # Make HTAN patient color legend
  if (strip_names) { names(htan_cols) <- gsub('HTA', '', names(htan_cols))}
  if (rm_duplicate_pair) {htan_cols <- htan_cols[!names(htan_cols) %in% c('HTA9-1_p2', '9-1_p2')]}
  patientLgd <- Legend(labels = names(htan_cols), legend_gp = gpar(fill = htan_cols), title = 'Patient ID', nrow = lgd_rows, 
                       labels_gp = gpar(fontsize = lgd_fontsize), grid_height = unit(lgd_gridsize, 'mm'), grid_width = unit(lgd_gridsize, 'mm'),
                       title_gp = gpar(fontsize = lgd_fontsize, fontface = "bold"))
  
  
  return(list('opAstrLgd' = opAstrLgd, 
              'onProgLgd' = onProgLgd,
              'intrinsicLgd' = intrinsicLgd, # change to responseLgd
              'pamLgd' = pamLgd, 
              'pamChangeLgd' = pamChangeLgd,
              'patientLgd' = patientLgd, 
              'siteLgd' = siteLgd))
  
}

# Convenience function to create all cohort heatmap annotations (except assay availability annotations)
make_heatmap_annotations <- function(meta) {
  
  # Treatment legend that also indicates on-progression biopsies
  onProg.idx <- rep(NA, ncol(meta))
  onProg.idx[which(meta$ProgressionStage == 'OnProgression')] <- 8
  onProgAnno <- HeatmapAnnotation('CDK4/6i' = anno_simple(meta$Treatment, pch = onProg.idx, col = all_drug_cols, border = TRUE), show_legend = FALSE, border = TRUE)
  
  # Make legend for ERi therapy
  erAnno <- HeatmapAnnotation(ERi = meta$ERi, show_legend = FALSE, border = TRUE, col = list(ERi = all_extra_drug_cols))
  
  # Make response annotation object and legend list
  intrinsicAnno <- HeatmapAnnotation(PatientResponse = meta$PatientResponse, show_legend = FALSE, border = TRUE, col = list(PatientResponse = intrinsic_cols))
  
  # PAM50 annotations
  pamAnno <- HeatmapAnnotation(PAM50 = meta$PAM50, show_legend = FALSE, col = list(PAM50 = pam_colors), na_col = 'grey', border = TRUE)
  
  # PAM50 change annotations
  pamChangeAnno <- HeatmapAnnotation(PAM50 = meta$pamChange, show_legend = FALSE, col = list(PAM50 = pam_change_colors), na_col = 'grey', border = TRUE)
  
  # Color code biopsies by patient
  # TODO: UPDATE TO USE HTAN PATIENT IDs and HTAN COLORS
  patientAnno <- HeatmapAnnotation(PatientID = meta$Patient, show_legend = FALSE, col = list(PatientID = c(patient_cols)), border = TRUE)
  
  # Metastatic site annotations
  siteAnno <- HeatmapAnnotation(BiopsySite = meta$BiopsySite, show_legend = FALSE, col = list(BiopsySite = site_colors), na_col = 'grey', border = TRUE)
  
  # Make custom column name annotations for working with HTAN patients
  htanAnno <- make_htan_pointers(meta)
  htanPointer <- htanAnno[['pointer']]
  htanNames <- htanAnno[['colNames']]
  htanPIDs = htanAnno[['patientNames']] # Might not need anymore
  htanBiopChange = htanAnno[['biopChange']]
  
  # Return list of annotation objects  
  return(list('onProgAnno' = onProgAnno, 
              'erAnno' = erAnno,
              'intrinsicAnno' = intrinsicAnno, # change to "responseAnno"
              'pamAnno' = pamAnno, 
              'pamChangeAnno' = pamChangeAnno,
              'siteAnno' <- siteAnno,
              'patientAnno' = patientAnno,
              'htanPointer' = htanPointer,
              'htanNames' = htanNames,
              'biopChangeHTAN' = htanBiopChange))
        
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

