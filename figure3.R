
#######################################################################################
#
#   FIGURE 3: MULTIMODAL PROFILING REVEALS TUMOR-INTRINSIC ADAPTATIONS TO CDK4/6i
#
#######################################################################################

# Source function libraries
source('function_scripts/load_htan_data.R')
source('function_scripts/heatmap_functions.R')
source('function_scripts/load_marker_sets.R')




###################################
#
#   LOAD DATA
#
###################################

source_dir <- '/Users/eggerj/OneDrive - Oregon Health & Science University/SMMART_HR+/Manuscript/Supp-Tables-Data/source_data'

# Meta table
meta.htan <- load_meta(fn.dir = source_dir)

# GSVA enrichment scores (RNA-seq)
gsva.htan <- load_gsva(meta.htan, fn.dir = source_dir)

# Gene expression (RNA-seq)


# Viper regulator activity (RNA-seq)


# Protein/phosphoprotein abundance (RPPA)


# mIHC



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

# List of legends for heatmaps showing delta between pairs during therapy
lgds.change <- make_heatmap_legends(meta.htan, select_samples = htan.onProgression)



# Annotations for delta heatmaps
# TODO: May not have to define separate annotations if all legends will be merged 
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


############################################################
#
#   LOAD TABLES FOR GENE AND PATHWAY SETS
#
############################################################


source('~/Documents/CompBio/HRplus_Project/manuscript_repo/HTAN_HR_Figures/function_scripts/load_marker_sets.R')




####################################################################################
#
#  FIGURE 3A: TRANSCRIPTIONAL GENE SET VARIATION ANALYSIS (GSVA) OF EXTRINSIC
#             PATHWAYS DURING CDK4/6I THERAPY AND PROGRESSION
#
####################################################################################

# Extrinsic pathways from Mann-Whitney test (p < 0.1)
gsva_pws.extrinsic <- c("Activated CD8 T cell", "ALLOGRAFT_REJECTION", "Antigen Presentation", 
                        "COMPLEMENT", "Gamma delta T cell", "iDC", "IL2_STAT5_SIGNALING", 
                        "IL6_JAK_STAT3_SIGNALING", "INFLAMMATORY_RESPONSE", "INTERFERON_ALPHA_RESPONSE", 
                        "INTERFERON_GAMMA_RESPONSE", "KEGG_JAK_STAT_SIGNALING_PATHWAY", 
                        "NK CD56dim cells", "REACTOME_PD_1_SIGNALING", "Regulatory T cell", 
                        "T Cell Inflamed GEP", "Tem cells")


# CLUSTER HEATMAP USING CHANGE ACROSS CDK4/6i
ht.fn <- paste(results_dir.htan, "rna_figures/heatmaps/gsva_test_heatmap.png", sep = '/')

# TODO: Change function name
make_heatmap(gsva.htan,
             meta.htan,
             select_samples = htan.onProgression,
             top_anno = top_annotations.change.htan,  
             btm_anno = btm_annotations.change.htan,
             category_table = gsva_cats.htan, 
             bar_anno = 'CytoChange',
             cluster_columns = TRUE,
             select_features = gsva_pws.extrinsic,
             show_column_annotation_legend = FALSE,
             heatmap_width = unit(7.5, 'in'),
             heatmap_height = unit(6, 'in'),
             add_width = 0.5,
             res = 100,
             split_column_by_dendrogram = 2,
             compute_change = TRUE,
             lgd_name = 'Activity Change',
             split_by_cat = TRUE, 
             fn = ht.fn
)





##########################################################################
#
#  FIGURE 3B: INTEGRATED HEATMAP OF EXTRINSIC RNA AND PROTEIN 
#             CHANGES DURING CDK4/6I THERAPY AND PROGRESSION
#
##########################################################################




###############################################################################################
#
#  FIGURE 3C: Cell densities reported for all CD45+ immune cells and individual immune 
#             cell type populations across pre-treatment (Bx1, white) and on-progression 
#             (Bx2, dark grey) biopsies for all biopsy pairs using mIHC
#
###############################################################################################



#####################################################################################
#
#  FIGURE 3D: Psuedocolored mIHC whole slide image and select regions of interest 
#             of 9-15 Bx2 showing spatial heterogeneity of CD8 marker (magenta)
#
#####################################################################################



####################################################################################
#
#  FIGURE 3E: Spatial heterogeneity and cytotoxic function of 
#             CD8+ T cells during CDK4/6i therapy and progression
#
####################################################################################


