#!/usr/bin/env Rscript

#######################################################################################
#
#   FIGURE 2: MULTIMODAL PROFILING REVEALS TUMOR-INTRINSIC ADAPTATIONS TO CDK4/6i
#
#######################################################################################

# Source function libraries
#source('function_scripts/load_htan_data.R')
#source('function_scripts/heatmap_functions.R')
#source('function_scripts/load_marker_sets.R')
#source('function_scripts/onco_heatmap_functions.R')

# Load functions
# (temp paths)
source('~/Documents/CompBio/HRplus_Project/manuscript_repo/HTAN_HR_Figures/function_scripts/htan_utils.R')
source('/Users/eggerj/Documents/CompBio/HRplus_Project/manuscript_repo/HTAN_HR_Figures/function_scripts/heatmap_functions.R')
source('/Users/eggerj/Documents/CompBio/HRplus_Project/manuscript_repo/HTAN_HR_Figures/function_scripts/onco_heatmap_functions.R')
source('/Users/eggerj/Documents/CompBio/HRplus_Project/manuscript_repo/HTAN_HR_Figures/function_scripts/multi_assay_functions.R')



# Load sample data, gene/pathways sets, and annotations
source('/Users/eggerj/Documents/CompBio/HRplus_Project/manuscript_repo/HTAN_HR_Figures/function_scripts/load_htan_data.R')
source('/Users/eggerj/Documents/CompBio/HRplus_Project/manuscript_repo/HTAN_HR_Figures/function_scripts/load_marker_sets.R')
source('~/Documents/CompBio/HRplus_Project/manuscript_repo/HTAN_HR_Figures/function_scripts/htan_annotations.R')


# Temp: Store test figures here
results_dir.test <- "/Users/eggerj/Documents/CompBio/HRplus_Project/manuscript_repo/test_figures"

# To retest outputs
results_dir.test <- "/Users/eggerj/Documents/CompBio/HRplus_Project/manuscript_repo/test_figures2"




####################################################################################
#
#  FIGURE 2A: ONCOPLOT OF GENOMIC ALTERATIONS FROM HIGH DEPTH TARGETED DNA-SEQ
#
####################################################################################

# Oncoplot of select variants from relevant pathways with ordered columns
make_oncoplots(
  
  # CNV/SNV Calls
  cnvs.htan, 
  snvs.htan,
  
  # Sample meta data and annotations
  meta.htan,
  select_samples = htan.paired,
  bottom_anno = list(onco_annotations.htan[[1]], NULL),
  
  # Select variants and annotations
  select_variants = dna_cats.main$Gene,
  category_table = dna_cats.main,
  min_vars = 1,
  
  # Oncoplot arguments
  pre = paste0(results_dir.test,'/figure2A'), 
  show_pct = FALSE,
  cluster_columns = FALSE, 
  cluster_rows = FALSE,
  fix_order = TRUE,
  ht_width = unit(10.5, 'in'),
  ht_height = unit(7.7, 'in')
  
)



#################################################################################
#
#  FIGURE 2B: TRANSCRIPTIONAL GENE SET VARIATION ANALYSIS (GSVA) OF INTRINSIC
#             PATHWAYS DURING CDK4/6I THERAPY AND PROGRESSION
#
#################################################################################

# Note: make sure "REACTOME_REPLICATION_STRESS" is named appropriately



# Select malignant cell pathways with average delta > 0.5 between G1 arrest and G1 entry tumors
pws.mc.htan <- gsva.change.htan %>%
  
  # Convert to long format and merge with meta data and pathway categories
  melt() %>%
  setNames(c('Pathway', 'Sample', 'Delta')) %>%
  left_join(meta.htan) %>%
  left_join(gsva_cats.main) %>%
  filter(Sample %in% htan.onProgression,
         !grepl('Immune', Category)) %>%
  
  # Compute mean delta by group per pathway
  group_by(g1Pheno, Pathway) %>%
  mutate(MeanDelta = mean(Delta)) %>%
  ungroup() %>%
  select(g1Pheno, Pathway, MeanDelta) %>%
  distinct() %>%
  
  # Compute difference in mean delta between groups per pathway
  group_by(Pathway) %>%
  mutate(GroupDiff = diff(MeanDelta)) %>%
  ungroup() %>%
  
  # Select significant pathways for heatmap
  arrange(Pathway) %>%
  filter(abs(GroupDiff) > 0.5) %>%
  pull(Pathway) %>%
  as.character() %>%
  unique()
  
  

# Fixed order for groups
pw_order.mc <- c(
  "Cell Cycle", 
  "Replication Stress", 
  "PI3K/AKT/mTOR", 
  "Metabolic"
)


# Heatmap of delta activity for intrinsic expression pathways
make_heatmap(
  
  # GSVA scores
  gsva.change.htan,
  
  # Sample meta data and annotations
  meta.htan,
  select_samples = htan.onProgression,
  top_anno = top_annotations.change.htan,  
  btm_anno = btm_annotations.change.htan,
  bar_anno = 'g1Score',
  
  # Select pathways and annotations
  select_features = pws.mc.htan,
  category_table = gsva_cats.main, 
  cat_order = pw_order.mc,
  split_by_cat = TRUE, 
  
  # Heatmap arguments
  fn = file.path(results_dir.test, 'figure2B.png'),
  lgd_name = 'Activity Change',
  cluster_columns = TRUE,
  split_column_by_dendrogram = 2,
  show_column_annotation_legend = FALSE,
  heatmap_width = unit(5.5, 'in'),
  heatmap_height = unit(6.25, 'in'),
  add_width = -0.35
  
)



###########################################################################
#
#  FIGURE 2C: INTEGRATED HEATMAP OF INTRINSIC RNA AND PROTEIN 
#             CHANGES DURING CDK4/6I THERAPY AND PROGRESSION
#
###########################################################################

merged_cats.mc <- c("G0", 
                    "G1",
                    "G1/S",
                    "G2/M", 
                    "G2/M CHECKPOINT", 
                    "mTORC1", 
                    "mTORC2", 
                    "PI3K/AKT")


# Malignant cell RNA and protein marker sets
merged_rna_protein.mc <- merged_rna_protein_cats.main %>%
  filter(Category %in% merged_cats.mc) %>%
  pull(MergedName)


# Multi-assay heatmap of delta values for malignant cell markers in RNA and protein modalities
multi_assay_heatmap(
  multi_modal.change.htan, 
  meta.htan, 
  
  pre = results_dir.test,
  fn = 'figure2C',
  
  
  # Heatmap groupings (pre-treatment to on-progression deltas)
  group_heatmaps_by = 'BiopsyChange.Drug', # Default: Sample
  
  top_anno = top_annotations.multiassay.change.htan,
  btm_anno = btm_annotations.multiassay.change.htan,
  
  # Features and annotations
  select_features = merged_rna_protein.mc,
  category_table = merged_rna_protein_cats.main,
  sub_sep = c(' '),
  
  # Heatmap parameters
  add_width = .6,
  ht_width = unit(1.75, 'in'),
  ht_height = unit(11.2,'in'),
  annotate_assay_types = FALSE, # REMOVE
  show_annotation_legend = FALSE,
  value.var = 'Zchange'
)




###########################################################################
#
#  FIGURE 2D: INTRINSIC PROTEOMIC PATHWAY SIGNALING 
#             DURING CDK4/6I THERAPY AND PROGRESSION
#
###########################################################################

# Select intrinsic protein pathways to match intrinsic GSVA gene sets
ppws.intrinsic <- c("Cell_cycle_progression", 
                    "G0_G1", 
                    "G1_S", 
                    "G2_M", 
                    "G2M_Checkpoint", 
                    "TSC_mTOR")

ppws.mc <- c("Cell_cycle_progression", 
             "G0_G1", 
             "G1_S", 
             "G2_M",
             "G2M_Checkpoint", 
             "TSC_mTOR")



make_heatmap(
  ppws.htan, 
  meta.htan,
  select_samples = htan.paired,
  top_anno = top_annotations.split.htan,  
  btm_anno = btm_annotations.split.htan,
  category_table = ppw_cats.main, 
  select_features = ppws.intrinsic,
  split_column_by_pheno = 'Patient',
  cluster_columns = FALSE,
  heatmap_width = unit(4.5, 'in'),
  heatmap_height = unit(5.5, 'in'),
  add_width = -1.5,
  compute_change = FALSE,
  show_column_annotation_legend = FALSE,
  lgd_name = 'Activity',
  fn = file.path(results_dir.test, "rppa_pathways_test_heatmap.png")
)


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






