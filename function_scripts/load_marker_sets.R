
##########################################################
#
#    GENE AND PATHWAY SETS FOR HTAN HR+ MANUSCRIPT
#
##########################################################

library(dplyr)

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

dna_cats.main <- list("Cell Cycle" = c("AURKA", 
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
dna_cats.main <- make_category_table(dna_cats.main, 'Gene')


########################################
#
#   GENES AND REGULATORS (RNA)
#
########################################

gene_cats.main <- list("G0" = c("CDKN1A", 
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
gene_cats.main <- make_category_table(gene_cats.main, 'Gene') %>%
  arrange(Category, Gene)



#######################################################
#
#   RNA PATHWAY CATEGORIES FOR MAIN FIGURES
#
#######################################################


gsva_cats.main <- list('Cell Cycle' = c("E2F_TARGETS", 
                                      "MYC_TARGETS_V1", 
                                      "MYC_TARGETS_V2",
                                      "MITOTIC_SPINDLE",
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
                                            #"REACTOME_PD_1_SIGNALING",
                                            "INTERFERON_ALPHA_RESPONSE", 
                                            "INTERFERON_GAMMA_RESPONSE"
                                            ),
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

gsva_cats.main <- make_category_table(gsva_cats.main, 'Pathway')


###############################################################
#
#   RNA PATHWAY CATEGORIES FOR SUPPLEMENTAL FIGURES
#
###############################################################


# Merged categories
gsva_cats.extend <- list("Cellular Component" = c("APICAL_JUNCTION", 
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
                                              #"REACTOME_PD_1_SIGNALING",
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
                                           "REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS",
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

gsva_cats.extend <- make_category_table(gsva_cats.extend, 'Pathway')

# Mix of categories
gsva_cats.plus <- rbind(gsva_cats.main,
                        gsva_cats.extend %>%
                          filter(!Pathway %in% gsva_cats.main$Pathway))

# Add Reactome adapative immune sets
#reactome_cats <- list('Adaptive Immune' = c("Antigen Presentation: Folding, assembly and peptide loading of class I MHC", 
#                                            "Antigen activates B Cell Receptor (BCR) leading to generation of second messengers", 
#                                            "Antigen processing-Cross presentation", "Antigen processing: Ubiquitination & Proteasome degradation", 
#                                            "CD22 mediated BCR regulation", "CD28 co-stimulation", "CTLA4 inhibitory signaling", 
#                                            "Downstream TCR signaling", "Downstream signaling events of B Cell Receptor (BCR)", 
#                                            "Generation of second messenger molecules", "PD-1 signaling", 
#                                            "Phosphorylation of CD3 and TCR zeta chains", "Translocation of ZAP-70 to Immunological synapse"
#))




add_reactome <- FALSE
add_reactome <- TRUE
if (add_reactome) {
reactome_cats <- list('Adaptive Immune' = c("REACTOME_TCR_SIGNALING",
                                            "REACTOME_CD28_CO_STIMULATION", 
                                            "REACTOME_CTLA4_INHIBITORY_SIGNALING", 
                                            "REACTOME_PD_1_SIGNALING", 
                                            #"REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION", 
                                            #"REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC", 
                                            "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION"),
                      
                      'Cytokine Signaling' = c("REACTOME_INTERLEUKIN_7_SIGNALING", 
                                               "REACTOME_INTERLEUKIN_1_FAMILY_SIGNALING", 
                                               "REACTOME_INTERLEUKIN_12_FAMILY_SIGNALING", 
                                               "REACTOME_INTERLEUKIN_17_SIGNALING", 
                                               "REACTOME_OTHER_INTERLEUKIN_SIGNALING", 
                                               "REACTOME_INTERLEUKIN_2_FAMILY_SIGNALING", 
                                               "REACTOME_INTERLEUKIN_3_INTERLEUKIN_5_AND_GM_CSF_SIGNALING", 
                                               "REACTOME_INTERLEUKIN_6_FAMILY_SIGNALING", 
                                               "REACTOME_INTERLEUKIN_10_SIGNALING", 
                                               "REACTOME_INTERLEUKIN_4_AND_INTERLEUKIN_13_SIGNALING", 
                                               "REACTOME_INTERLEUKIN_20_FAMILY_SIGNALING", 
                                               "REACTOME_ANTIVIRAL_MECHANISM_BY_IFN_STIMULATED_GENES", 
                                               "REACTOME_INTERFERON_GAMMA_SIGNALING", 
                                               "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING", 
                                               #"REACTOME_STAT5_ACTIVATION", 
                                               "REACTOME_NEGATIVE_REGULATION_OF_FLT3", 
                                               "REACTOME_FLT3_SIGNALING_THROUGH_SRC_FAMILY_KINASES"
                      )
                      )


split_reactome_cats <- FALSE
split_reactome_cats <- TRUE

if (split_reactome_cats) {
reactome_cats <- list('Adaptive Immune' = reactome_cats[['Adaptive Immune']],
                      
                      'Interferon Signaling' = c("REACTOME_ANTIVIRAL_MECHANISM_BY_IFN_STIMULATED_GENES", 
                                                 "REACTOME_INTERFERON_GAMMA_SIGNALING", 
                                                 "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING"),
                      
                      'Interleukin Signaling' = c("REACTOME_INTERLEUKIN_7_SIGNALING", 
                                                  "REACTOME_INTERLEUKIN_1_FAMILY_SIGNALING", 
                                                  "REACTOME_INTERLEUKIN_12_FAMILY_SIGNALING", 
                                                  "REACTOME_INTERLEUKIN_17_SIGNALING", 
                                                  "REACTOME_OTHER_INTERLEUKIN_SIGNALING", 
                                                  "REACTOME_INTERLEUKIN_2_FAMILY_SIGNALING", 
                                                  "REACTOME_INTERLEUKIN_3_INTERLEUKIN_5_AND_GM_CSF_SIGNALING", 
                                                  "REACTOME_INTERLEUKIN_6_FAMILY_SIGNALING", 
                                                  "REACTOME_INTERLEUKIN_10_SIGNALING", 
                                                  "REACTOME_INTERLEUKIN_4_AND_INTERLEUKIN_13_SIGNALING", 
                                                  "REACTOME_INTERLEUKIN_20_FAMILY_SIGNALING"),
                      
                      'FLT3 Signaling' = c("REACTOME_NEGATIVE_REGULATION_OF_FLT3", 
                                           "REACTOME_FLT3_SIGNALING_THROUGH_SRC_FAMILY_KINASES")
                      
                      )
}


# Just adaptive immune
reactome_cats <- list('Adaptive Immune' = c("REACTOME_TCR_SIGNALING",
                                            "REACTOME_CD28_CO_STIMULATION", 
                                            "REACTOME_CTLA4_INHIBITORY_SIGNALING", 
                                            "REACTOME_PD_1_SIGNALING", 
                                            #"REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION", 
                                            #"REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC", 
                                            "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION"))

reactome_cats <- make_category_table(reactome_cats, 'Pathway')
gsva_cats.plus <- rbind(gsva_cats.plus, reactome_cats)

}

# Hallmarks 
hallmarks <- c("ADIPOGENESIS", "ALLOGRAFT_REJECTION", "ANDROGEN_RESPONSE", 
               "ANGIOGENESIS", "APICAL_JUNCTION", "APICAL_SURFACE", "APOPTOSIS", 
               "BILE_ACID_METABOLISM", "CHOLESTEROL_HOMEOSTASIS", "COAGULATION", 
               "COMPLEMENT", "DNA_REPAIR", "E2F_TARGETS", "EPITHELIAL_MESENCHYMAL_TRANSITION", 
               "ESTROGEN_RESPONSE_EARLY", "ESTROGEN_RESPONSE_LATE", "FATTY_ACID_METABOLISM", 
               "G2M_CHECKPOINT", "GLYCOLYSIS", "HEDGEHOG_SIGNALING", "HEME_METABOLISM", 
               "HYPOXIA", "IL2_STAT5_SIGNALING", "IL6_JAK_STAT3_SIGNALING", 
               "INFLAMMATORY_RESPONSE", "INTERFERON_ALPHA_RESPONSE", "INTERFERON_GAMMA_RESPONSE", 
               "KRAS_SIGNALING_DN", "KRAS_SIGNALING_UP", "MITOTIC_SPINDLE", 
               "MTORC1_SIGNALING", "MYC_TARGETS_V1", "MYC_TARGETS_V2", "MYOGENESIS", 
               "NOTCH_SIGNALING", "OXIDATIVE_PHOSPHORYLATION", "P53_PATHWAY", 
               "PANCREAS_BETA_CELLS", "PEROXISOME", "PI3K_AKT_MTOR_SIGNALING", 
               "PROTEIN_SECRETION", "REACTIVE_OXYGEN_SPECIES_PATHWAY", "SPERMATOGENESIS", 
               "TGF_BETA_SIGNALING", "TNFA_SIGNALING_VIA_NFKB", "UNFOLDED_PROTEIN_RESPONSE", 
               "UV_RESPONSE_DOWN", "UV_RESPONSE_UP", "WNT_BETA_CATENIN_SIGNALING", 
               "XENOBIOTIC_METABOLISM")

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
#   PATHWAY LISTS
#
###############################################################



# Immune cell types
pws.celltypes <- gsva_cats.main %>%
  filter(Category == 'Immune Cell Activity') %>%
  pull(Pathway)



###############################################################
#
#   PROTEIN PATHWAY CATEGORIES FOR SUPPLEMENTAL FIGURES
#
###############################################################

ppw_cats.main <- list('Cell Cycle' = c("Cell_cycle_progression", 
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

ppw_cats.main <- make_category_table(ppw_cats.main, 'Pathway')
