
##########################################################
#
#    GENE AND PATHWAY SETS FOR HTAN HR+ MANUSCRIPT
#
##########################################################


# Helper function to make feature category table from named list
make_category_table <-  function(feature_sets, feature_name = 'Gene') {
  
  # Melt list of features into table 
  cat_tbl <- feature_sets %>%
    melt() %>%
    setNames(c(feature_name, 'Category')) %>%
    mutate(!!feature_name := as.character(.data[[feature_name]]),
           Category = as.character(Category)) %>%
    arrange(Category, .data[[feature_name]]) %>%
    data.frame()
  rownames(cat_tbl) <- cat_tbl[,feature_name]
  
  return(cat_tbl)

}

########################################
#
#   GENE ALTERATION CATAGORIES (DNA)
#
########################################

dna_cats.htan <- list("Cell Cycle" = c("AURKA", 
                                       "BRD4", 
                                       "CCND1", 
                                       "CDKN1B",  
                                       "CDKN2A", 
                                       "CDKN2B", 
                                       "RB1"), 
                      "DNA Damage" = c("ARID1A",
                                       "BRCA2",  
                                       "MSH6", 
                                       "PARP1", 
                                       "SLFN11"), 
                      "MAPK" = c("BRAF"), 
                      "P53" = c("MDM4", 
                                "TP53"), 
                      "PI3K/AKT/mTOR" = c("AKT3", 
                                          "PIK3CA",  
                                          "PRKCA",  
                                          "RICTOR", 
                                          "RPS6KB1"), 
                      "Receptor" = c("ERBB2",  
                                     "ESR1", 
                                     "FGFR1",
                                     "GATA3",  
                                     "KMT2C"), 
                      "Replication Stress" = c("ATM", 
                                               "ATR", 
                                               "CHEK1", 
                                               "CHEK2"))

# Make into table and return
dna_cats.htan <- make_category_table(dna_cats.htan, 'Gene')


########################################
#
#   GENES AND REGULATORS (RNA)
#
########################################

gene_cats.htan <- list("G0" = c("CDKN1A", 
                                "CDKN1B", 
                                "RB1"), 
                       "G1" = c("CCND1", 
                                "CCND3", 
                                "CDK4", 
                                "CDK6"), 
                       "G1/S" = c("BRD4", 
                                  "CCNE1", 
                                  "CDK2", 
                                  "E2F1", 
                                  "E2F2", 
                                  "PCNA"), 
                       "G2/M" = c("AURKA", 
                                  "CCNB1", 
                                  "CDC25C", 
                                  "CDK1", 
                                  "CDK7", 
                                  "PLK1"), 
                       "G2/M CHECKPOINT" = c("ATM", 
                                             "ATR", 
                                             "CHEK1", 
                                             "CHEK2", 
                                             "HIST1H3A", 
                                             "RPA2", 
                                             "WEE1"), 
                       "Immune" = c("CD274", 
                                    "PD-L1", # alt ID
                                    "CD276", 
                                    "CD38", 
                                    "CD4", 
                                    "CD68", 
                                    "GZMB", 
                                    "HLA-DRA", 
                                    "IDO1", 
                                    "IRF1", 
                                    "IRF3", 
                                    "LCK", 
                                    "PDCD1", #alt ID
                                    "PD-1", 
                                    "TNFRSF4", 
                                    "VTCN1", 
                                    "ZAP70"), 
                       "JAK/STAT" = c("IL6", 
                                      "JAK2", 
                                      "MUC1",   
                                      "PTPRC", 
                                      "STAT1", 
                                      "STAT3", 
                                      "STAT5A"),  
                       "mTORC1" = c("AKT1S1", 
                                    "EEF2K", 
                                    "EIF4E", 
                                    "EIF4EBP1", 
                                    "RPS6", 
                                    "RPS6KB1", 
                                    "RPTOR"), 
                       "mTORC2" = c("NDRG1", 
                                    "PRKCA", 
                                    "RICTOR", 
                                    "SGK1"), 
                       "PI3K/AKT" = c("AKT1", 
                                      "GSK3A", 
                                      "MTOR", 
                                      "PDPK1", 
                                      "PIK3CA", 
                                      "PTEN", 
                                      "TSC1", 
                                      "TSC2"), 
                       "Receptor" = c("AR", 
                                      "EGFR", 
                                      "ERBB2", 
                                      "ERBB3", 
                                      "ESR1", 
                                      "FOXM1", 
                                      "GATA3", 
                                      "MET", 
                                      "PGR"))


# Make into table and return
gene_cats.htan <- make_category_table(gene_cats.htan, 'Gene') %>%
  arrange(Category, Gene)



#######################################################
#
#   RNA PATHWAY CATEGORIES FOR MAIN FIGURES
#
#######################################################


gsva_cats.htan <- list('Cell Cycle' = c("E2F_TARGETS", 
                                      "MYC_TARGETS_V1", 
                                      "KEGG_DNA_REPLICATION",  
                                      "REACTOME_S_PHASE", 
                                      "REACTOME_CELL_CYCLE"),
                     'EMT' = c("EPITHELIAL_MESENCHYMAL_TRANSITION"),
                     'Replication Stress' = c("G2M_CHECKPOINT", 
                                              "DNA_REPAIR",  
                                              "REACTOME_REPLICATION_STRESS",
                                              "REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS",
                                              "REACTOME_G2_M_CHECKPOINTS", 
                                              'REACTOME_G2_M_DNA_DAMAGE_CHECKPOINT'),
                     'PI3K/AKT/mTOR' = c("MTORC1_SIGNALING", 
                                         "PI3K_AKT_MTOR_SIGNALING",
                                         "KEGG_MTOR_SIGNALING_PATHWAY"),
                     'Cellular Process' = c("PROTEIN_SECRETION", 
                                            "UNFOLDED_PROTEIN_RESPONSE"),
                     'Metabolic' = c("OXIDATIVE_PHOSPHORYLATION"),
                     'Estrogen Signaling' = c("ESTROGEN_RESPONSE_EARLY", 
                                              "ESTROGEN_RESPONSE_LATE", 
                                              "ANDROGEN_RESPONSE"),
                     'JAK/STAT' = c("BIOCARTA_STAT3_PATHWAY", 
                                    "REACTOME_STAT5_ACTIVATION",
                                    "IL6_JAK_STAT3_SIGNALING",
                                    "IL2_STAT5_SIGNALING",
                                    "KEGG_JAK_STAT_SIGNALING_PATHWAY"),
                     'Immune Signaling' = c('Antigen Presentation',
                                            'Chemokine',
                                            'T Cell Inflamed GEP',
                                            "ALLOGRAFT_REJECTION", 
                                            "COAGULATION", 
                                            "COMPLEMENT", 
                                            "INFLAMMATORY_RESPONSE", 
                                            "INTERFERON_ALPHA_RESPONSE", 
                                            "INTERFERON_GAMMA_RESPONSE",
                                            "REACTOME_PD_1_SIGNALING"),
                     'Immune Cell Activity' = c("Activated CD8 T cell", 
                                                "B cells", 
                                                "Cytotoxic cells", 
                                                "Eosinophils", 
                                                "Gamma delta T cell", 
                                                "iDC", 
                                                "Macrophages", 
                                                "Mast cells", 
                                                "Neutrophils", 
                                                "NK CD56bright cells", 
                                                "NK CD56dim cells", 
                                                "Regulatory T cell", 
                                                "T helper cells", 
                                                "Tcm cells", 
                                                "Tem cells", 
                                                "Tfh cells"))

gsva_cats.htan <- make_category_table(gsva_cats.htan, 'Pathway')


###############################################################
#
#   RNA PATHWAY CATEGORIES FOR SUPPLEMENTAL FIGURES
#
###############################################################


# Merged categories
gsva_cats.merged <- list("Cellular Component" = c("APICAL_JUNCTION", 
                                                "APICAL_SURFACE", 
                                                "PEROXISOME"), 
                       "Development" = c("ADIPOGENESIS", 
                                         "ANGIOGENESIS", 
                                         "EPITHELIAL_MESENCHYMAL_TRANSITION", 
                                         "MYOGENESIS", 
                                         "PANCREAS_BETA_CELLS", 
                                         "SPERMATOGENESIS"), 
                       "DNA Damage" = c("DNA_REPAIR", 
                                        "UV_RESPONSE_DOWN", 
                                        "UV_RESPONSE_UP"), 
                       "Immune Cell Activity" = c("B cells", 
                                                  "Eosinophils", 
                                                  "Macrophages", 
                                                  "Mast cells",
                                                  "NK CD56bright cells", 
                                                  "NK CD56dim cells", 
                                                  "Neutrophils", 
                                                  "T helper cells", 
                                                  "Tcm cells", 
                                                  "Tem cells", 
                                                  "Tfh cells", 
                                                  "iDC", "Activated CD8 T cell",
                                                  "Gamma delta T cell", 
                                                  "Regulatory T cell", 
                                                  "Cytotoxic cells"), 
                       "Immune Signaling" = c("ALLOGRAFT_REJECTION", 
                                              "COAGULATION", 
                                              "COMPLEMENT", 
                                              "INFLAMMATORY_RESPONSE", 
                                              "INTERFERON_ALPHA_RESPONSE", 
                                              "INTERFERON_GAMMA_RESPONSE", 
                                              "REACTOME_PD_1_SIGNALING",
                                              "Antigen Presentation", 
                                              "Chemokine",
                                              'T Cell Inflamed GEP'), 
                       "JAK/STAT" = c("IL2_STAT5_SIGNALING", 
                                      "IL6_JAK_STAT3_SIGNALING", 
                                      "BIOCARTA_STAT3_PATHWAY", 
                                      "REACTOME_STAT5_ACTIVATION", 
                                      "KEGG_JAK_STAT_SIGNALING_PATHWAY"), 
                       "Metabolic" = c("BILE_ACID_METABOLISM", 
                                       "CHOLESTEROL_HOMEOSTASIS", 
                                       "FATTY_ACID_METABOLISM", 
                                       "GLYCOLYSIS", 
                                       "HEME_METABOLISM", 
                                       "OXIDATIVE_PHOSPHORYLATION", 
                                       "XENOBIOTIC_METABOLISM"), 
                       "Pathway" = c("APOPTOSIS", 
                                     "HYPOXIA", 
                                     "PROTEIN_SECRETION", 
                                     "REACTIVE_OXYGEN_SPECIES_PATHWAY", 
                                     "UNFOLDED_PROTEIN_RESPONSE"), 
                       "Proliferation" = c("E2F_TARGETS", 
                                           "G2M_CHECKPOINT", 
                                           "MITOTIC_SPINDLE", 
                                           "MYC_TARGETS_V1", 
                                           "MYC_TARGETS_V2", 
                                           "P53_PATHWAY", 
                                           "KEGG_DNA_REPLICATION", 
                                           "REACTOME_S_PHASE", 
                                           "REACTOME_REPLICATION_STRESS", 
                                           "REACTOME_CELL_CYCLE"), 
                       "Receptor Signaling" = c("ANDROGEN_RESPONSE", 
                                                "ESTROGEN_RESPONSE_EARLY", 
                                                "ESTROGEN_RESPONSE_LATE"), 
                       "Signaling" = c("HEDGEHOG_SIGNALING", 
                                       "KRAS_SIGNALING_DN", 
                                       "KRAS_SIGNALING_UP", 
                                       "MTORC1_SIGNALING", 
                                       "NOTCH_SIGNALING", 
                                       "PI3K_AKT_MTOR_SIGNALING", 
                                       "TGF_BETA_SIGNALING", 
                                       "TNFA_SIGNALING_VIA_NFKB", 
                                       "WNT_BETA_CATENIN_SIGNALING"))

gsva_cats.merged <- make_category_table(gsva_cats.merged, 'Pathway')

# Color scheme for row annotations
gsva_colors <- c("Cellular Component" = "#EC579AFF", 
                 "Development" = "#7F8624FF", 
                 "DNA Damage" = "#0C5BB0FF", 
                 "Metabolic" = "#149BEDFF", 
                 "Pathway" = "#FA6B09FF", 
                 "Proliferation" = "#FEC10BFF", 
                 "Signaling" = "purple", 
                 "JAK/STAT" = "turquoise3", 
                 "Immune Signaling" = "#EE0011FF", 
                 "Receptor Signaling" = "#15983DFF"
)

###############################################################
#
#   PROTEIN PATHWAY CATEGORIES FOR SUPPLEMENTAL FIGURES
#
###############################################################

ppw_cats.htan <- list('Cell Cycle' = c("Cell_cycle_progression", 
                                      "G0_G1", 
                                      "G1_S", 
                                      "G2_M"),
                     'Replication Stress' = c("G2M_Checkpoint"),
                     'Receptor Signaling' = c("Hormone_receptor", 
                                              "Hormone_signaling_Breast", 
                                              "RTK"),
                     'Immune Signaling' = c("Immune", 
                                            "Immune_Checkpoint", 
                                            "JAK/STAT"),
                     'Apoptosis' = c('Apoptosis', 
                                     'BH3_Balance'),
                     'Pathway Signaling' = c("Notch", 
                                             "RAS_MAPK"),
                     'Tumor' = c('Tumor_Content'),
                     'PI3K/AKT/mTOR' = c("PI3K_Akt", 
                                         "TSC_mTOR"))

ppw_cats.htan <- make_category_table(pw_cats.rppa, 'Pathway')
