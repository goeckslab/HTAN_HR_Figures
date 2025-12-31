

source('/Users/eggerj/Documents/CompBio/HRplus_Project/manuscript_repo/HTAN_HR_Figures/function_scripts/heatmap_functions.R')


#########################################################################
#
#   LOAD AND SET ANNOTATION OBJECTS FOR HEATMAPS AND ONCOPLOTS
#
#########################################################################

# Load list of all heatmap and oncoplot annotations
annotations.htan <- make_heatmap_annotations(meta.htan)

# List of legends for heatmaps showing all (paired) samples
lgds.htan <- make_heatmap_legends(meta.htan, select_samples = htan.paired)

# List of legends for heatmaps showing delta between pairs during therapy
lgds.change <- make_heatmap_legends(meta.htan, select_samples = htan.onProgression)



#  -- Oncoplot annotations (single time point) -- #
onco_annotations.htan <- list(
  
  # Annotations
  annotations.htan$onProgAnno %v% 
    annotations.htan$erAnno %v%
    annotations.htan$responseAnno %v%
    annotations.htan$pamAnno %v%
    annotations.htan$patientAnno %v% 
    annotations.htan$htanPointer %v% 
    annotations.htan$sampleAnno, 
  
  # Legends
  list(lgds.htan$opAstrLgd,
       lgds.htan$onProgLgd,
       lgds.htan$responseLgd,
       lgds.htan$pamLgd)
)


# -- Annotations for all samples with patient ID annotation -- #
top_annotations.htan <- list(
  
  # Annotations
  annotations.htan$onProgAnno %v% 
    annotations.htan$erAnno %v%
    annotations.htan$responseAnno, 
  
  # Legends
  list(lgds.htan$opAstrLgd,
       lgds.htan$onProgLgd, 
       lgds.htan$responseLgd)
)


btm_annotations.htan <- list(
  
  # Ananotations
  annotations.htan$pamAnno %v%
    annotations.htan$patientAnno %v% 
    annotations.htan$htanPointer %v% 
    annotations.htan$sampleAnno, 
  
  # Legends
  list(lgds.htan$pamLgd)
  
)



# -- Delta Heatmap Annotation -- #
top_annotations.change.htan <- list(
  
  # Annotations
  annotations.htan$onProgAnno %v% 
    annotations.htan$erAnno %v%
    annotations.htan$responseAnno, 
  
  # Legends
  list(lgds.change$opAstrLgd,
       lgds.change$onProgLgd, 
       lgds.change$responseLgd)
  
)


btm_annotations.change.htan <- list(
  
  # Heatmaps
  annotations.htan$pamChangeAnno %v%
    annotations.htan$htanPointer %v% 
    annotations.htan$biopPairAnno, 
  
  # Legends
  list(lgds.change$pamChangeLgd)
  
)



# Annotations for heatmaps split by patient (doesn't need patient ID annotation)
top_annotations.split.htan <- list(
  
  # Heatmaps
  annotations.htan$onProgAnno %v% 
    annotations.htan$erAnno %v%
    annotations.htan$responseAnno, 
  
  # Legends
  list(lgds.htan$opAstrLgd,
       lgds.htan$onProgLgd, 
       lgds.htan$responseLgd)
  
)

btm_annotations.split.htan <- list(
  
  # Heatmaps
  annotations.htan$pamAnno %v%
    annotations.htan$htanPointer %v% 
    annotations.htan$sampleAnno, 
  
  # Legends
  list(lgds.htan$pamLgd)
  
)


# Top annotations for multi-assay heatmaps
top_annotations.multiassay.change.htan <- list(
  CDKi = list(anno_name = "CDK4/6i", 
              anno_colors = colors.treatment
  ),
  ERi  = list(anno_colors = colors.treatment
  )
)

# Bottom annotations for multi-assay heatmaps
btm_annotations.multiassay.change.htan <- list(
  pamChange = list(anno_name = "PAM50", 
                   anno_colors = colors.pam
  )
)

