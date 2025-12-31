

source('~/Documents/CompBio/HRplus_Project/manuscript_repo/HTAN_HR_Figures/function_scripts/htan_utils.R')

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





##########################################################################
#
#   LOAD MULTI-OMIC DATASETS
#
##########################################################################

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



# Multi-modal matrices (Pre-treatment to On-Progression Deltas)

multi_modal.change.htan <- merge_assays(
  meta.htan, 
  df.rna = exp.scaled.htan, 
  df.viper = viper.scaled.htan, 
  df.rppa = rppa.scaled.htan,
  
  scale_rna = FALSE,
  scale_viper = FALSE,
  scale_rppa = FALSE,
  
  
  protein_rna_table = protein_to_rna.htan, 
  
  merged_rna_protein_table = merged_rna_protein_names.htan,
  select_merged_names = merged_rna_protein_cats.main$MergedName,
  
  select_samples = htan.paired,
  patient_column = 'Patient.Drug',
  
  fill_all_assays = TRUE,
  
  
  Zchange = TRUE
)

