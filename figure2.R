

#######################################################################################
#
#   FIGURE 2: MULTIMODAL PROFILING REVEALS TUMOR-INTRINSIC ADAPTATIONS TO CDK4/6i
#
#######################################################################################

# Source function libraries






###################################
#
#   LOAD DATA
#
###################################

# Meta table


# Genomic alterations (DNA-seq) 


# GSVA enrichment scores (RNA-seq)


# Gene expression (RNA-seq)


# Viper regulator activity (RNA-seq)


# Protein/phosphoprotein abundance (RPPA)


# Proteomic pathway scores (RPPA)


# Chromatin accessibility enrichment (sciATAC-seq)


# Tile density (CycIF)


# Ripley's K (CycIF)



#########################################################################
#
#   LOAD AND CONSTRUCT ANNOTATION OBJECTS FOR HEATMAPS AND ONCOPLOTS
#
#########################################################################

# Load list of all annotations
annotations.hrplus <- make_cohort_annotations(samples.hrplus, meta.hrplus)

# List of heatmap legends specifically for htan samples
lgds.htan <- make_cohort_legends(samples.hrplus, meta.hrplus, select_samples = htan.paired)

# For heatmaps showing change in value (only need on-progression labels)
lgds.change <- make_cohort_legends(samples.hrplus, meta.hrplus, select_samples = htan.onProgression)



# Oncoplot annotations
onco_annotations.htan <- list(annotations.hrplus$onProgAnno %v% 
                                annotations.hrplus$erAnno %v%
                                annotations.hrplus$intrinsicAnno %v%
                                annotations.hrplus$pamAnno %v%
                                annotations.hrplus$patientAnno %v% 
                                annotations.hrplus$htanPointer %v% 
                                annotations.hrplus$htanNames, 
                              list(lgds.htan$opAstrLgd,
                                   lgds.htan$onProgLgd,
                                   lgds.htan$intrinsicLgd,
                                   lgds.htan$pamLgd))


# Annotations for all samples with patient ID annotation
top_annotations.htan <- list(annotations.hrplus$onProgAnno %v% 
                               annotations.hrplus$erAnno %v%
                               annotations.hrplus$intrinsicAnno, 
                             list(lgds.htan$opAstrLgd,
                                  lgds.htan$onProgLgd, 
                                  lgds.htan$intrinsicLgd))
btm_annotations.htan <- list(annotations.hrplus$pamAnno %v%
                               annotations.hrplus$patientAnno %v% 
                               annotations.hrplus$htanPointer %v% 
                               annotations.hrplus$htanNames, 
                             list(lgds.htan$pamLgd))



# Annotations for delta heatmaps
top_annotations.change.htan <- list(annotations.hrplus$onProgAnno %v% 
                                      annotations.hrplus$erAnno %v%
                                      annotations.hrplus$intrinsicAnno, 
                                    list(lgds.change$opAstrLgd,
                                         lgds.change$onProgLgd, 
                                         lgds.change$intrinsicLgd))
btm_annotations.change.htan <- list(annotations.hrplus$pamChangeAnno %v%
                                      annotations.hrplus$htanPointer %v% 
                                      annotations.hrplus$biopChangeHTAN, 
                                    list(lgds.change$pamChangeLgd))



# Select top and bottom annotations for HTAN samples when splitting by patient (doesn't need patient ID annotation)
top_annotations.split.htan <- list(annotations.hrplus$onProgAnno %v% 
                                     annotations.hrplus$erAnno %v%
                                     annotations.hrplus$intrinsicAnno, 
                                   list(lgds.htan$opAstrLgd,
                                        lgds.htan$onProgLgd, 
                                        lgds.htan$intrinsicLgd))
btm_annotations.split.htan <- list(annotations.hrplus$pamAnno %v%
                                     annotations.hrplus$htanPointer %v% 
                                     annotations.hrplus$htanNames, 
                                   list(lgds.htan$pamLgd))


############################################################
#
#   LOAD TABLES FOR GENE AND PATHWAY SETS
#
############################################################







##################################################################
#
#  FIGURE 2A: 
#
##################################################################



##################################################################
#
#  FIGURE 2B: 
#
##################################################################



##################################################################
#
#  FIGURE 2C: 
#
##################################################################



##################################################################
#
#  FIGURE 2D: 
#
##################################################################



##################################################################
#
#  FIGURE 2E: 
#
##################################################################



##################################################################
#
#  FIGURE 2F: 
#
##################################################################






