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




# Compare change in activity between G1 arrest and G1 entry patients
pvals.G1score.htan <- two_sample_test(gsva.change.htan, 
                                      meta.htan, 
                                      select_samples = htan.onProgression,
                                      features = gsva_cats.main %>% 
                                        filter(!grepl('Immune', Category)) %>%
                                        pull(Pathway),
                                      pheno = 'g1Pheno', testType = 'wilcox', 
                                      paired_test = FALSE)

two_sample_test(gsva.change.htan, 
                meta.htan, 
                select_samples = htan.onProgression,
                features = gsva_cats.main %>% 
                  filter(!grepl('Immune', Category)) %>%
                  pull(Pathway),
                pheno = 'g1Pheno', testType = 'wilcox', 
                paired_test = FALSE) %>%
  arrange(-abs(MeanChange)) %>%
  filter(abs(MeanChange) > 0.5) %>%
  pull(Feature) %>%
  as.character()


# Select top malignant cell pathways by average absolute delta
pws.mc.htan <- pvals.G1score.htan %>% 
  arrange(-abs(MeanChange)) %>%
  filter(abs(MeanChange) > 0.5) %>%
  pull(Feature) %>%
  as.character()


# Intrinsic pathways from Mann-Whitney test (p < 0.1)
gsva_pws.intrinsic <- c(
  "E2F_TARGETS",
  "G2M_CHECKPOINT",
  "KEGG_DNA_REPLICATION",
  "MTORC1_SIGNALING",
  "MYC_TARGETS_V1",
  "MYC_TARGETS_V2",
  "OXIDATIVE_PHOSPHORYLATION",
  "REACTOME_CELL_CYCLE",
  "REACTOME_REPLICATION_STRESS",
  "REACTOME_S_PHASE"
)

# Fixed order for groups
gsva_cat.order.intrinsic <- c(
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
  select_features = gsva_pws.intrinsic,
  category_table = gsva_cats.main, 
  cat_order = gsva_cat.order.intrinsic,
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


# Replicating figures created in htan_scripts/htan_rna_integration.R

# TODO: RENAME RETURNED MATRICES

assays.test.htan <- merge_assays(meta.htan, 
                                 df.rna = exp.scaled.htan, 
                                 df.viper = viper.scaled.htan, 
                                 df.rppa = rppa.scaled.htan,
                                 
                                 scale_rna = FALSE,
                                 scale_viper = FALSE,
                                 scale_rppa = FALSE,
                                 
                                 
                                 protein_rna_table = protein_to_rna.htan, 
                                 
                                 merged_rna_protein_table = merged_rna_protein_names.htan,
                                 select_merged_names = merged_rna_protein_cats.main$MergedName,
                                 
                                 #select_genes = gene_cats.main$Gene[1:40],
                                 
                                 select_samples = htan.paired,
                                 patient_column = 'Patient.Drug',
                                 
                                 fill_all_assays = TRUE,
                                 
                                 
                                 Zchange = TRUE)







merged_rna_protein.intrinsic <- merged_rna_protein_cats.main %>%
  filter(Category %in% c("G0", "G1", "G1/S", "G2/M", "G2/M CHECKPOINT", #"Immune", "JAK/STAT", 
                          "mTORC1", "mTORC2", "PI3K/AKT")) %>%
  pull(MergedName)


multi_assay_heatmap(assays.test.htan, 
                    meta.htan, 
                    
                    pre = results_dir.test,
                    fn = 'test_merged_assays',
                    
                    
                    # Use to select how heatmaps are grouped (default "Sample", but "BiopsyChange.Drug" for paired delta values)
                    group_heatmaps_by = 'BiopsyChange.Drug', # Default: Sample
                    #group_heatmaps_by = 'pamChange', # Default: Sample
                    
                    
                    
                    top_anno = top_annotations.multiassay.change.htan,
                    btm_anno = btm_annotations.multiassay.change.htan,
                    
                    
                    # Contains merged names?
                    category_table = merged_rna_protein_cats.main,
                    sub_sep = c(' '),
                    
                    # TODO: change parameter name
                    order_rows = merged_rna_protein.intrinsic,
                    
                    
                    add_width = .6,
                    ht_width = unit(1.75, 'in'),
                    ht_height = unit(11.2,'in'),
                    annotate_assay_types = FALSE,
                    show_annotation_legend = FALSE,
                    value.var = 'Zchange')




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




# Heatmap splitting by patient
ht.fn <- paste(results_dir.test, "rppa_pathways_test_heatmap.png", sep = '/')

make_heatmap(ppws.htan, 
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






