

######################################################
#
#   HEATMAP FUNCTIONS FOR HTAN HR+ MANUSCRIPT
#
######################################################

library(ComplexHeatmap)


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

