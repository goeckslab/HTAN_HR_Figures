#!/usr/bin/env Rscript

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

# NEW ANNOTATIONS AFTER PERFORMING PARK VALIDATION

# USE BITS OF CODE BELOW


# Make G1 arrest barplot annotation function
g1Anno.merged <- make_barplot_annotation(meta.merged, anno = 'g1Score', return_bar_anno = TRUE)

# TCI delta barplot annotation function
tcDeltaAnno.merged <- make_barplot_annotation(meta.merged, anno = 'TCIDelta', 
                                              ylim = c(-2,2),
                                              return_bar_anno = TRUE)

# MC cluster annotation function
mcAnno.merged <- make_anno_simple(meta = meta.merged, anno = 'MC Cluster', 
                                  anno_colors = colors.intrinsic, return_anno_simple = TRUE)

# TIME cluster annotation function
timeAnno.merged <- make_anno_simple(meta = meta.merged, anno = 'TIME Cluster', 
                                    anno_colors = colors.extrinsic, return_anno_simple = TRUE)




# Make separate annotation of extrinsic clusters just for MMTERT set
# TODO: change name
ct.extrinsic.hrplus <- make_anno_simple(meta.merged[rownames(meta.hrplus),'TIME Cluster'], 
                                        anno_name = 'TIME Cluster', anno_colors = colors.extrinsic)


# USE FOR FIGURE 3A


top.celltypes.htan <- list(annotations.hrplus$tciDeltaAnno %v%
                             ct.extrinsic.hrplus,
                           list(ct.extrinsic.merged$ctLgd))


# Move rest of annotations to bottom of heatmap
btm.celltypes.htan <- list(annotations.hrplus$onProgAnno %v%
                             annotations.hrplus$erAnno %v%
                             annotations.hrplus$intrinsicAnno %v%
                             annotations.hrplus$pamChangeAnno %v%
                             annotations.hrplus$htanPointer %v%
                             annotations.hrplus$biopChangeHTAN,
                           list(lgds.change$opAstrLgd,
                                lgds.change$onProgLgd,
                                lgds.change$intrinsicLgd, # still using? new name?
                                lgds.change$pamChangeLgd,
                                lgds.change$cohortLgd))




# USE IMMUNE CELL TYPES INSTEAD OF BELOW

# Extrinsic pathways from Mann-Whitney test (p < 0.1)
gsva_pws.extrinsic <- c("Activated CD8 T cell", "ALLOGRAFT_REJECTION", "Antigen Presentation", 
                        "COMPLEMENT", "Gamma delta T cell", "iDC", "IL2_STAT5_SIGNALING", 
                        "IL6_JAK_STAT3_SIGNALING", "INFLAMMATORY_RESPONSE", "INTERFERON_ALPHA_RESPONSE", 
                        "INTERFERON_GAMMA_RESPONSE", "KEGG_JAK_STAT_SIGNALING_PATHWAY", 
                        "NK CD56dim cells", "REACTOME_PD_1_SIGNALING", "Regulatory T cell", 
                        "T Cell Inflamed GEP", "Tem cells")


# CLUSTER HEATMAP USING CHANGE ACROSS CDK4/6i
ht.fn <- paste(results_dir.test, "gsva_test_heatmap.png", sep = '/')

# TODO: Change function name
make_heatmap(gsva.htan,
             meta.htan,
             select_samples = htan.onProgression,
             
             # TODO: SWAP WITH NEW ANNOTATIONS
             top_anno = top_annotations.change.htan,  
             btm_anno = btm_annotations.change.htan,
             
             
             #category_table = gsva_cats.htan, 
             #select_features = gsva_pws.extrinsic,
             
             select_features = pws.celltypes,
             
             bar_anno = 'CytoChange',
             cluster_columns = TRUE,
             
             
             
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



merged_cats.time <- c("Immune", "JAK/STAT")


# Malignant cell RNA and protein marker sets
merged_rna_protein.time <- merged_rna_protein_cats.main %>%
  filter(Category %in% merged_cats.time) %>%
  pull(MergedName)


# Multi-assay heatmap of delta values for malignant cell markers in RNA and protein modalities
multi_assay_heatmap(
  multi_modal.change.htan, 
  meta.htan, 
  
  pre = results_dir.test,
  fn = 'figure3B',
  
  
  # Heatmap groupings (pre-treatment to on-progression deltas)
  group_heatmaps_by = 'BiopsyChange.Drug', # Default: Sample
  
  top_anno = top_annotations.multiassay.change.htan,
  btm_anno = btm_annotations.multiassay.change.htan,
  
  # Features and annotations
  select_features = merged_rna_protein.time,
  category_table = merged_rna_protein_cats.main,
  sub_sep = c(' '),
  
  # Heatmap parameters
  # TODO: UPDATE SIZSE TO MATCH ORIGINAL
  add_width = .6,
  ht_width = unit(1.75, 'in'),
  ht_height = unit(11.2,'in'),
  annotate_assay_types = FALSE, # REMOVE
  show_annotation_legend = FALSE,
  value.var = 'Zchange'
)



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


