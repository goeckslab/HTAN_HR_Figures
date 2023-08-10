

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
                      "Post-Tamoxifen" = "#FFD3E5")
                      
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
make_heatmap_legends <- function(samples_table, meta, select_samples = NULL, lgd_rows = 2,
                                 treatment_title = ' ', strip_names = TRUE, rm_duplicate_pair = TRUE,
                                 lgd_fontsize = 12, lgd_gridsize = 4) {
  
  # All color schemes subsetted down to select patients
  if (!is.null(select_samples)) {
    
    colors.treatment <- colors.treatment[names(colors.treatment) %in% c(meta[select_samples, 'CDKi'], meta[select_samples, 'ERi'])]
    colors.response <- colors.response[names(colors.response) %in% meta[select_samples, 'PatientResponse']]
    colors.pam <- colors.pam[names(colors.pam) %in% meta[select_samples, 'PAM50']]
    colors.pam_change <- colors.pam[names(colors.pam) %in% meta[select_samples, 'pamChange']]
    colors.site <- colors.site[names(colors.site) %in% meta[select_samples, 'BiopsySite']]
    colors.patient <- colors.patient[names(colors.patient) %in% meta[select_samples, 'Patient']]
    
  }
  
  # Legend for pre-treatment and on-progression with on-progression labeled with asterisk
  opAstrLgd <- Legend(title = 'Treatment', type = 'points', pch = c(NA,8), at = c('PreTreatment', 'OnProgression'), background = c('white'),
                      labels_gp = gpar(fontsize = lgd_fontsize), grid_height = unit(lgd_gridsize, 'mm'), grid_width = unit(lgd_gridsize, 'mm'),
                      title_gp = gpar(fontsize = lgd_fontsize, fontface = "bold"))
  
  # Legend for CDK4/6i and ERi treatment (merged with above)
  onProgLgd <- Legend(labels = names(colors.treatment), legend_gp = gpar(fill = colors.treatment), title = treatment_title, 
                      nrow = lgd_rows, labels_gp = gpar(fontsize = lgd_fontsize), grid_height = unit(lgd_gridsize, 'mm'), 
                      grid_width = unit(lgd_gridsize, 'mm'), title_gp = gpar(fontsize = lgd_fontsize, fontface = "bold"))
  
  
  # Patient response legend
  responseLgd <- Legend(labels = names(colors.response), legend_gp = gpar(fill = colors.response), title = 'Patient Response', nrow = lgd_rows, 
                         labels_gp = gpar(fontsize = lgd_fontsize), grid_height = unit(lgd_gridsize, 'mm'), grid_width = unit(lgd_gridsize, 'mm'),
                         title_gp = gpar(fontsize = lgd_fontsize, fontface = "bold"))
  
  # PAM50 annotation legend
  pamLgd <- Legend(labels = names(colors.pam), legend_gp = gpar(fill = colors.pam), title = 'PAM50', nrow = lgd_rows, 
                   labels_gp = gpar(fontsize = lgd_fontsize), grid_height = unit(lgd_gridsize, 'mm'), grid_width = unit(lgd_gridsize, 'mm'),
                   title_gp = gpar(fontsize = lgd_fontsize, fontface = "bold"))
  
  # PAM50 change annotation legend
  pamChangeLgd <- Legend(labels = names(colors.pam_change), legend_gp = gpar(fill = colors.pam_change), title = 'PAM50', nrow = lgd_rows, 
                         labels_gp = gpar(fontsize = lgd_fontsize), grid_height = unit(lgd_gridsize, 'mm'), grid_width = unit(lgd_gridsize, 'mm'),
                         title_gp = gpar(fontsize = lgd_fontsize, fontface = "bold"))
  
  # Metastatic site legend
  siteLgd <- Legend(labels = names(colors.site), legend_gp = gpar(fill = colors.site), title = 'Biopsy Site', nrow = lgd_rows, 
                    labels_gp = gpar(fontsize = lgd_fontsize), grid_height = unit(lgd_gridsize, 'mm'), grid_width = unit(lgd_gridsize, 'mm'),
                    title_gp = gpar(fontsize = lgd_fontsize, fontface = "bold"))
  
  # Patient legend
  patientLgd <- Legend(labels = names(colors.patient), legend_gp = gpar(fill = colors.patient), title = 'Patient ID', nrow = lgd_rows, 
                       labels_gp = gpar(fontsize = lgd_fontsize), grid_height = unit(lgd_gridsize, 'mm'), grid_width = unit(lgd_gridsize, 'mm'),
                       title_gp = gpar(fontsize = lgd_fontsize, fontface = "bold"))
  
  
  return(list('opAstrLgd' = opAstrLgd, 
              'onProgLgd' = onProgLgd,
              'responseLgd' = responseLgd, 
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
  onProgAnno <- HeatmapAnnotation('CDK4/6i' = anno_simple(meta$CDKi, pch = onProg.idx, col = colors.treatment, border = TRUE), show_legend = FALSE, border = TRUE)
  
  # Make legend for ERi therapy
  erAnno <- HeatmapAnnotation(ERi = meta$ERi, show_legend = FALSE, border = TRUE, col = list(ERi = colors.treatment))
  
  # Make response annotation object and legend list
  responseAnno <- HeatmapAnnotation(PatientResponse = meta$PatientResponse, show_legend = FALSE, border = TRUE, col = list(PatientResponse = colors.response))
  
  # PAM50 annotations
  pamAnno <- HeatmapAnnotation(PAM50 = meta$PAM50, show_legend = FALSE, col = list(PAM50 = colors.pam), na_col = 'grey', border = TRUE)
  
  # PAM50 change annotations
  pamChangeAnno <- HeatmapAnnotation(PAM50 = meta$pamChange, show_legend = FALSE, col = list(PAM50 = colors.pam), na_col = 'grey', border = TRUE)
  
  # Color code biopsies by patient
  patientAnno <- HeatmapAnnotation(PatientID = meta$Patient, show_legend = FALSE, col = list(PatientID = c(colors.patient)), border = TRUE)
  
  # Metastatic site annotations
  siteAnno <- HeatmapAnnotation(BiopsySite = meta$BiopsySite, show_legend = FALSE, col = list(BiopsySite = colors.site), na_col = 'grey', border = TRUE)
  
  # Make custom column name annotations for working with HTAN patients
  htanAnno <- make_htan_pointers(meta)
  htanPointer <- htanAnno[['pointer']]
  htanNames <- htanAnno[['colNames']]
  htanPIDs = htanAnno[['patientNames']] # Might not need anymore
  htanBiopChange = htanAnno[['biopChange']]
  
  # Return list of annotation objects  
  return(list('onProgAnno' = onProgAnno, 
              'erAnno' = erAnno,
              'responseAnno' = responseAnno, # change to "responseAnno"
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

