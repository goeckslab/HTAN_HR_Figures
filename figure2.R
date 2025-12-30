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

# TEMP
source('/Users/eggerj/Documents/CompBio/HRplus_Project/manuscript_repo/HTAN_HR_Figures/function_scripts/load_htan_data.R')
source('/Users/eggerj/Documents/CompBio/HRplus_Project/manuscript_repo/HTAN_HR_Figures/function_scripts/heatmap_functions.R')
source('/Users/eggerj/Documents/CompBio/HRplus_Project/manuscript_repo/HTAN_HR_Figures/function_scripts/load_marker_sets.R')
source('/Users/eggerj/Documents/CompBio/HRplus_Project/manuscript_repo/HTAN_HR_Figures/function_scripts/onco_heatmap_functions.R')
source('/Users/eggerj/Documents/CompBio/HRplus_Project/manuscript_repo/HTAN_HR_Figures/function_scripts/multi_assay_functions.R')

# Temp: Store test figures here
results_dir.test <- "/Users/eggerj/Documents/CompBio/HRplus_Project/manuscript_repo/test_figures"

# To retest outputs
results_dir.test <- "/Users/eggerj/Documents/CompBio/HRplus_Project/manuscript_repo/test_figures2"


###################################
#
#   LOAD MULTI-OMIC DATASETS
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

# GSVA paired delta values (RNA-seq)
gsva.change.htan <- compute_paired_change(gsva.htan, meta.htan)


# Gene expression (RNA-seq)
exp.scaled.htan <- load_gene_expression(meta.htan, fn.dir = source_dir)

# Viper regulator activity (RNA-seq)
viper.scaled.htan <- load_viper(meta.htan, fn.dir = source_dir)

# Protein/phosphoprotein abundance (RPPA)
rppa.scaled.htan <- load_htan_rppa(meta.htan, fn.dir = source_dir)

# Proteomic pathway scores (RPPA)
ppws.htan <- load_ppws(meta.htan, fn.dir = source_dir)

# Chromatin accessibility enrichment (sciATAC-seq)


# Tile density (CycIF)


# Ripley's K (CycIF)


# Accessory tables
protein_to_rna.htan <- load_protein_to_rna(fn.dir = source_dir)
merged_rna_protein_names.htan <- load_merged_rna_to_protein_names(fn.dir = source_dir)


#########################################################################
#
#   DEFINE SAMPLE SETS USING HTAN IDS
#
#########################################################################

# Paired tumor groups (for single time point)
htan.paired <- c(
  # HTA9-1
  "HTA9-1_Bx1",
  "HTA9-1_Bx2",
  "HTA9-1_Bx4",
  "HTA9-1_Bx5",
  
  # HTA9-2
  "HTA9-2_Bx1",
  "HTA9-2_Bx2",
  
  # HTA9-3
  "HTA9-3_Bx1",
  "HTA9-3_Bx2",
  
  # HTA9-14
  "HTA9-14_Bx1",
  "HTA9-14_Bx2",
  "HTA9-14_Bx3",
  
  # HTA9-15
  "HTA9-15_Bx1",
  "HTA9-15_Bx2"
)


# On-progression (for delta values)
htan.onProgression <- c(
  # HTA9-1
  "HTA9-1_Bx2",
  "HTA9-1_Bx5",
  
  # HTA9-2
  "HTA9-2_Bx2",
  
  # HTA9-3
  "HTA9-3_Bx2",
  
  # HTA9-14
  "HTA9-14_Bx2",
  "HTA9-14_Bx3",
  
  # HTA9-15
  "HTA9-15_Bx2"
)



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
  fn = file.path(results_dir.test, 'figure2B'),
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






