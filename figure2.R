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
source('/Users/eggerj/Documents/CompBio/HRplus_Project/manuscript_repo/HTAN_HR_Figures/function_scripts/sciATAC_functions.R')



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
  
  # Multi-assay delta z-scores
  multi_modal.change.htan, 
  
  # Sample data and annotations
  meta = meta.htan, 
  top_anno = top_annotations.multiassay.change.htan,
  btm_anno = btm_annotations.multiassay.change.htan,
  group_heatmaps_by = 'BiopsyChange.Drug',
  
  # Features and annotations
  select_features = merged_rna_protein.mc,
  category_table = merged_rna_protein_cats.main,
  sub_sep = c(' '),
  
  # Heatmap parameters
  pre = results_dir.test,
  fn = 'figure2C',
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


# Protein pathwey activity heatmap of intrinsic pathways
make_heatmap(
  
  # Protein pathway activity
  ppws.htan, 
  
  # Sample meta data and annotations
  meta.htan,
  select_samples = htan.paired,
  top_anno = top_annotations.split.htan,  
  btm_anno = btm_annotations.split.htan,
  split_column_by_pheno = 'Patient',
  
  # Features and annotations
  category_table = ppw_cats.main, 
  select_features = ppws.intrinsic,
  
  # Heatmap parameters
  cluster_columns = FALSE,
  show_column_annotation_legend = FALSE,
  heatmap_width = unit(4.5, 'in'),
  heatmap_height = unit(5.5, 'in'),
  add_width = -1.5,
  lgd_name = 'Activity',
  fn = file.path(results_dir.test, "figure2D.png")
  
)


#############################################################################################
#
#  FIGURE 2E: VIOLIN PLOTS SHOWING DISTRIBUTIONS OF PER-CELL CHROMATIC ACCESSIBILITY 
#             ENRICHMENT SCORES FOR INTRINSIC PATHWAYS FROM SCIATAC-SEQ
#
#############################################################################################





# Malignant cell sciATAC pathways
scipws.mc.htan <- c(
  "E2F_TARGETS",
  "KEGG_DNA_REPLICATION",
  "MYC_TARGETS_V1",
  "REACTOME_CELL_CYCLE",
  "REACTOME_S_PHASE",
  "PROTEIN_SECRETION",
  "UNFOLDED_PROTEIN_RESPONSE",
  "OXIDATIVE_PHOSPHORYLATION",
  "MTORC1_SIGNALING",
  "G2M_CHECKPOINT",
  "REACTOME_REPLICATION_STRESS"
)





# test getting p-values ahead of time
test_pvals.sciATAC <- get_sciATAC_pathway_pvals(sciATAC.htan, 
                                                select_pathways = scipws.mc.htan, 
                                                meta = meta.htan,
                                                
                                                sample_group = 'Sample.Drug', 
                                                
                                                alt_group = 'rest')

# TODO: Make computing p-values optional in plotting function so can be skipped if pvals table already provided
#   make sure all filtering steps are also added to p-values function


# TODO: Can we make this into generic function

plot_list <- list()

for (c in unique(sciATAC_cats.main$Category)) {
  
  sciATAC_subPws <- sciATAC_cats.main %>%
    filter(Category == c,
           Pathway %in% scipws.mc.htan) %>%
    pull(Pathway)
  
  if (length(sciATAC_subPws)) {
    
    plot_list[[c]] <- plot_sciATAC_pathway_activity(
      sciATAC.htan, 
      meta = meta.htan,
      
      select_pathways = sciATAC_subPws,
      category_table = sciATAC_cats.main,
      
      x = 'Sample.Drug',
      col_by = 'Sample.Drug',
      celltypes = c('Tumor'), 
      
      hide_legend = TRUE,
      #hide_legend = FALSE,
      hide_x_axis = TRUE,
      order_by_sample = TRUE,
      
      alt_group = 'rest'
    )
    
  }
  
}



layout_mat <- t(matrix(c(1,1,1,1,
                         2,2,NA, NA,
                         3,3,NA,NA,
                         4,NA,NA,NA,
                         5,NA,NA,NA), 
                       nrow = 4))

g <- grid.arrange(plot_list$`Cell Cycle`, 
                  plot_list$`Replication Stress`,
                  plot_list$`Cellular Process`, 
                  plot_list$`PI3K/AKT/mTOR`, 
                  plot_list$Metabolic,
                  nrow = 4,
                  widths = c(3,2.2,2.2,2.2),
                  layout_matrix = layout_mat)



# May 19, 2025 Figure 2?
ggsave(g, width = 14, height = 12,
       filename = file.path(results_dir.test, 'figure2E.png'))







# Save legend
plot_sciATAC_pathway_activity(sciATAC_scores.htan, 
                           c('G2M_CHECKPOINT'),
                           category_table = sciATAC_cats.main,
                           hide_legend = FALSE,
                           hide_x_axis = TRUE,
                           order_by_sample = TRUE,
                           facet_formula = c('Category ~ Pathway'),
                           #facet_formula = c('Pathway ~ Category'),
                           celltypes = c('Tumor'), 
                           
                           alt_group = 'rest') %>% 
  
  # where is this?
  get_legend() %>%
  
  ggsave(width = 8, height = 8,
         filename = file.path(results_dir.test, 'figure2E_legend.png'))


####################################################################################
#
#  FIGURE 2F: QUANTIFIACTION OF SPATIAL HETEROGENEITY OF PROLIFERATING TUMOR 
#             CELLS DURING CDK4/6I THERAPY AND PROGRESSION USING CYCIF
#
####################################################################################






