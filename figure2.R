#!/usr/bin/env Rscript

#######################################################################################
#
#   FIGURE 2: MULTIMODAL PROFILING REVEALS TUMOR-INTRINSIC ADAPTATIONS TO CDK4/6i
#
#######################################################################################

# Source function libraries
source('function_scripts/load_htan_data.R')
source('function_scripts/heatmap_functions.R')
source('function_scripts/load_marker_sets.R')
source('function_scripts/onco_heatmap_functions.R')

# Temp: Store test figures here
results_dir.test <- "/Users/eggerj/Documents/CompBio/HRplus_Project/manuscript_repo/test_figures"


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
gsva.htan <- load_gsva(meta.htan, fn.dir = source_dir)

# Gene expression (RNA-seq)


# Viper regulator activity (RNA-seq)


# Protein/phosphoprotein abundance (RPPA)


# Proteomic pathway scores (RPPA)
ppws.htan <- load_ppws(meta.htan, fn.dir = source_dir)

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

# Oncoplot of select variants from relevant pathways with ordered columns
# TODO: Change file naming
make_oncoplots(cnvs.htan, 
               snvs.htan,
               meta.htan,
               select_samples = htan.paired,
               select_variants = dna_cats.htan$Gene,
               category_table = dna_cats.htan,
               
               pre = paste0(results_dir.test,'/oncoplot_test'), 
               
               min_vars = 1,
               
               show_pct = FALSE,
               cluster_columns = FALSE, 
               cluster_rows = FALSE,
               fix_order = TRUE,
               ht_width = unit(10.5, 'in'),
               ht_height = unit(7.7, 'in'),
               bottom_anno = list(onco_annotations.htan[[1]], NULL))



#################################################################################
#
#  FIGURE 2B: TRANSCRIPTIONAL GENE SET VARIATION ANALYSIS (GSVA) OF INTRINSIC
#             PATHWAYS DURING CDK4/6I THERAPY AND PROGRESSION
#
#################################################################################

# Intrinsic pathways from Mann-Whitney test (p < 0.1)
gsva_pws.intrinsic <- c("E2F_TARGETS", "G2M_CHECKPOINT", 
                        "KEGG_DNA_REPLICATION", "MTORC1_SIGNALING", 
                        "MYC_TARGETS_V1", "OXIDATIVE_PHOSPHORYLATION", 
                        "PROTEIN_SECRETION", "REACTOME_CELL_CYCLE", 
                        "REACTOME_REPLICATION_STRESS", "REACTOME_S_PHASE", 
                        "UNFOLDED_PROTEIN_RESPONSE")


# CLUSTER HEATMAP USING CHANGE ACROSS CDK4/6i
ht.fn <- paste(results_dir.test, "gsva_test_heatmap.png", sep = '/')

# TODO: Change function name
make_heatmap(gsva.htan,
             meta.htan,
             select_samples = htan.onProgression,
             top_anno = top_annotations.change.htan,  
             btm_anno = btm_annotations.change.htan,
             category_table = gsva_cats.htan, 
             cat_order = c("Cell Cycle", "Replication Stress", "PI3K/AKT/mTOR", 
                           "Cellular Process", "Metabolic"),
             bar_anno = 'g1Score',
             cluster_columns = TRUE,
             select_features = gsva_pws.intrinsic,
             show_column_annotation_legend = FALSE,
             heatmap_width = unit(5.5, 'in'),
             heatmap_height = unit(6.25, 'in'),
             add_width = -0.35,
             split_column_by_dendrogram = 2,
             compute_change = TRUE,
             lgd_name = 'Activity Change',
             split_by_cat = TRUE, 
             fn = ht.fn
)



###########################################################################
#
#  FIGURE 2C: INTEGRATED HEATMAP OF INTRINSIC RNA AND PROTEIN 
#             CHANGES DURING CDK4/6I THERAPY AND PROGRESSION
#
###########################################################################



###########################################################################
#
#  FIGURE 2D: INTRINSIC PROTEOMIC PATHWAY SIGNALING 
#             DURING CDK4/6I THERAPY AND PROGRESSION
#
###########################################################################

# Select intrinsic protein pathways to match intrinsic GSVA gene sets
ppws.intrinsic <- c("Cell_cycle_progression", "G0_G1", "G1_S", "G2_M",
                        "G2M_Checkpoint", "TSC_mTOR")



# Heatmap splitting by patient
ht.fn <- paste(results_dir.test, "rppa_pathways_test_heatmap.png", sep = '/')

make_heatmap(ppws.htan, 
             meta.htan,
             select_samples = htan.paired,
             top_anno = top_annotations.split.htan,  
             btm_anno = btm_annotations.split.htan,
             category_table = ppw_cats.htan, 
             select_features = ppws.intrinsic,
             split_column_by_pheno = 'Patient',
             cluster_columns = FALSE,
             heatmap_width = unit(4.5, 'in'),
             heatmap_height = unit(5.5, 'in'),
             add_width = -1.5,
             compute_change = FALSE,
             show_column_annotation_legend = FALSE,
             lgd_name = 'Activity',
             fn = ht.fn)

#############################################################################################
#
#  FIGURE 2E: VIOLIN PLOTS SHOWING DISTRIBUTIONS OF PER-CELL CHROMATIC ACCESSIBILITY 
#             ENRICHMENT SCORES FOR INTRINSIC PATHWAYS FROM SCIATAC-SEQ
#
#############################################################################################



####################################################################################
#
#  FIGURE 2F: QUANTIFIACTION OF SPATIAL HETEROGENEITY OF PROLIFERATING TUMOR 
#             CELLS DURING CDK4/6I THERAPY AND PROGRESSION USING CYCIF
#
####################################################################################






