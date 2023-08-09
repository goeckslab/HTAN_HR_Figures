

#######################################################################################
#
#   FIGURE 2: MULTIMODAL PROFILING REVEALS TUMOR-INTRINSIC ADAPTATIONS TO CDK4/6i
#
#######################################################################################

# Source function libraries
source('function_scripts/heatmap_functions.R')





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
#   DEFINE SAMPLE SETS USING HTAN IDS
#
#########################################################################


htan.paired <- c("HTA9-1_Bx1", "HTA9-1_Bx2", "HTA9-1_Bx4", "HTA9-1_Bx5", 
                 "HTA9-2_Bx1", "HTA9-2_Bx2", 
                 "HTA9-3_Bx1", "HTA9-3_Bx2", 
                 "HTA9-14_Bx1", "HTA9-14_Bx2", "HTA9-14_Bx3", 
                 "HTA9-15_Bx1", "HTA9-15_Bx2")


htan.onProgression <- c("HTA9-1_Bx2", "HTA9-1_Bx5", 
                        "HTA9-2_Bx2", 
                        "HTA9-3_Bx2", 
                        "HTA9-14_Bx2", "HTA9-14_Bx3", 
                        "HTA9-15_Bx2")



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


source('~/Documents/CompBio/HRplus_Project/manuscript_repo/HTAN_HR_Figures/function_scripts/load_marker_sets.R')




##################################################################
#
#  FIGURE 2A: 
#
##################################################################


# Oncoplot with fixes column order
make_oncoprints(cnvs.hrplus, snvs.hrplus, 
                meta.hrplus %>% mutate(HTAN = gsub('_p2', '', HTAN)), 
                pre = paste0(results_dir.htan,'/dna_figures/oncoplots/select_cats_fixed_'), min_vars = 1,
                select_variants = select_dna_cats$Gene,
                category_table = select_dna_cats,
                #select_samples = htan.hrplus,
                select_samples = htan.paired,
                show_pct = FALSE,
                cluster_columns = FALSE, cluster_rows = FALSE,
                fix_order = TRUE,
                #column_split = 'HTAN',
                ht_width = unit(10.5, 'in'),
                ht_height = unit(7.7, 'in'),
                bottom_anno = list(onco_annotations.htan[[1]], NULL))



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






