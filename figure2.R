

#######################################################################################
#
#   FIGURE 2: MULTIMODAL PROFILING REVEALS TUMOR-INTRINSIC ADAPTATIONS TO CDK4/6i
#
#######################################################################################

# Source function libraries

source('function_scripts/heatmap_functions.R')
source('function_scripts/load_marker_sets.R')
source('function_scripts/onco_heatmap_functions.R')




###################################
#
#   LOAD DATA
#
###################################

source_dir <- '/Users/eggerj/OneDrive - Oregon Health & Science University/SMMART_HR+/Manuscript/Supp-Tables-Data/source_data'

# Meta table
meta.htan <- load_meta(fn.dir = source_dir)

# Copy number alterations (DNA-seq) 
cnvs.htan <- load_cnvs(meta.htan, fn.dir = source_dir)

# Single nucleotide variants (DNA-seq) 
snvs.htan <- load_snvs(meta.htan, fn.dir = source_dir)

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
annotations.htan <- make_heatmap_annotations(meta.htan)

# List of legends for heatmaps showing all (paired) samples
lgds.htan <- make_heatmap_legends(meta.htan, select_samples = htan.paired)

# List of legends for heatmaps showing delta between pairs during therapy
lgds.change <- make_heatmap_legends(meta.htan, select_samples = htan.onProgression)



# Oncoplot annotations
onco_annotations.htan <- list(annotations.htan$onProgAnno %v% 
                                annotations.htan$erAnno %v%
                                annotations.htan$responseAnno %v%
                                annotations.htan$pamAnno %v%
                                annotations.htan$patientAnno %v% 
                                annotations.htan$htanPointer %v% 
                                annotations.htan$sampleAnno, 
                              list(lgds.htan$opAstrLgd,
                                   lgds.htan$onProgLgd,
                                   lgds.htan$responseLgd,
                                   lgds.htan$pamLgd))


# Annotations for all samples with patient ID annotation
top_annotations.htan <- list(annotations.htan$onProgAnno %v% 
                               annotations.htan$erAnno %v%
                               annotations.htan$responseAnno, 
                             list(lgds.htan$opAstrLgd,
                                  lgds.htan$onProgLgd, 
                                  lgds.htan$responseLgd))
btm_annotations.htan <- list(annotations.htan$pamAnno %v%
                               annotations.htan$patientAnno %v% 
                               annotations.htan$htanPointer %v% 
                               annotations.htan$sampleAnno, 
                             list(lgds.htan$pamLgd))



# Annotations for delta heatmaps
top_annotations.change.htan <- list(annotations.htan$onProgAnno %v% 
                                      annotations.htan$erAnno %v%
                                      annotations.htan$responseAnno, 
                                    list(lgds.change$opAstrLgd,
                                         lgds.change$onProgLgd, 
                                         lgds.change$responseLgd))
btm_annotations.change.htan <- list(annotations.htan$pamChangeAnno %v%
                                      annotations.htan$htanPointer %v% 
                                      annotations.htan$biopPairAnno, 
                                    list(lgds.change$pamChangeLgd))



# Annotations for heatmaps split by patient (doesn't need patient ID annotation)
top_annotations.split.htan <- list(annotations.htan$onProgAnno %v% 
                                     annotations.htan$erAnno %v%
                                     annotations.htan$responseAnno, 
                                   list(lgds.htan$opAstrLgd,
                                        lgds.htan$onProgLgd, 
                                        lgds.htan$responseLgd))
btm_annotations.split.htan <- list(annotations.htan$pamAnno %v%
                                     annotations.htan$htanPointer %v% 
                                     annotations.htan$sampleAnno, 
                                   list(lgds.htan$pamLgd))


############################################################
#
#   LOAD TABLES FOR GENE AND PATHWAY SETS
#
############################################################


source('~/Documents/CompBio/HRplus_Project/manuscript_repo/HTAN_HR_Figures/function_scripts/load_marker_sets.R')




####################################################################################
#
#  FIGURE 2A: ONCOPLOT OF GENOMIC ALTERATIONS FROM HIGH DEPTH TARGETED DNA-SEQ
#
####################################################################################

# TODO: Test when select_samples = NULL

# Oncoplot of select variants from relevant pathways with ordered columns
make_oncoplots(cnvs.htan, 
               snvs.htan,
               meta.htan,
               select_samples = htan.paired,
               select_variants = dna_cats.htan$Gene,
               category_table = dna_cats.htan,
               
               #pre = paste0(results_dir.htan,'/dna_figures/oncoplots/select_cats_fixed_'), 
               pre = paste0(results_dir.htan,'/dna_figures/oncoplots/test_'), 
               
               min_vars = 1,
               
               show_pct = FALSE,
               cluster_columns = FALSE, 
               cluster_rows = FALSE,
               fix_order = TRUE,
               
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






