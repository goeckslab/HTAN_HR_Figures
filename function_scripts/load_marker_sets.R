
##########################################################
#
#    GENE AND PATHWAY SETS FOR HTAN HR+ MANUSCRIPT
#
##########################################################


###############################
#
#   GENES AND REGULATORS
#
###############################

select_gene_cats <- list("G0" = c("CDKN1A", "CDKN1B", "CDKN1C", "CDKN2A", "CDKN2B", 
                                  "CDKN2C", "CDKN2D", "RB1", "RBL1"), 
                         "G1" = c("CCND1", "CCND2", 
                                  "CCND3", "CDK4", "CDK6"), 
                         "G1/S" = c("BRD4", "CCNA1", "CCNE1", 
                                    "CDC25A", "CDC6", "CDK2", 
                                    "E2F1", "E2F2", "E2F3", "MKI67", "MYC", 
                                    "PCNA", "TK1"), 
                         "G2/M" = c("AURKA", "CCNB1", "CDC25C", "CDK1", 
                                    "CDK7", "PLK1"), 
                         "G2/M CHECKPOINT" = c("ATM", "ATR", "CHEK1", 
                                               "CHEK2", "H2AFX", "HIST1H3A", 
                                               "RPA1", "RPA2", "WEE1"), 
                         "Immune" = c("CD274", "CD276", "CD38", "CD3D", "CD4", 
                                      "CD68", "CD8A", "CD8B", "CTLA4", 
                                      "FOXP3", "GZMA", "GZMB", "HLA-DRA", 
                                      "IDO1", "IFGR1", "IFGR2", 
                                      "IFNG", "IL10", "IL2RA", "IRF1", 
                                      "IRF3", "ITGAM", "LCK", "NFATC3", 
                                      "NFATC4", "PD-1", "PD-L1", "PDCD1", 
                                      "PRF1", "TNFRSF4", "VTCN1", 
                                      "ZAP70"), 
                         "JAK/STAT" = c("IL6", "IL6R", "JAK1", "JAK2", 
                                        "MUC1",   "PIAS1", "PTPRC", "SOCS1", 
                                        "STAT1", "STAT3", "STAT4", "STAT5A", 
                                        "STAT5B"), 
                         "mTORC1" = c("AKT1S1", "EEF2K", "EIF4E", "EIF4EBP1", 
                                      "RPS6", "RPS6KB1", "RPTOR"), 
                         "mTORC2" = c("MAPKAP1", "NDRG1", "PRKCA", 
                                      "RICTOR", "SGK1"), 
                         "PI3K/AKT" = c("AKT1", "AKT3", "GSK3A", "MTOR", 
                                        "PDPK1", "PIK3CA", "PIP", "PTEN", "TSC1", "TSC2"), 
                         "Receptor" = c("AR", 
                                        "EGFR", "ERBB2", "ERBB3", "ESR1", "FOXA1", "FOXM1", "GATA3", 
                                        "IGF1R", "INSR", "KMT2C", "MET", "NR3C1", "PGR"))


# Make into table and return
select_gene_cats <- make_category_table(select_gene_cats, 'Gene') %>%
  arrange(Category, Gene)



##############################################
#
#   PATHWAY CATEGORIES FOR MAIN FIGURES
#
##############################################



pw_cats.htan <- list('Cell Cycle' = c("E2F_TARGETS", 
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

pw_cats.htan <- make_category_table(pw_cats.htan, 'Pathway')


######################################################
#
#   PATHWAY CATEGORIES FOR SUPPLEMENTAL FIGURES
#
######################################################


# Merged categories
pw_cats.merged <- list("Cellular Component" = c("APICAL_JUNCTION", 
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
                                              "Chemokine"), 
                       "JAK/STAT" = c("IL2_STAT5_SIGNALING", "IL6_JAK_STAT3_SIGNALING", 
                                      "BIOCARTA_STAT3_PATHWAY", "REACTOME_STAT5_ACTIVATION", "KEGG_JAK_STAT_SIGNALING_PATHWAY"
                       ), 
                       "Metabolic" = c("BILE_ACID_METABOLISM", "CHOLESTEROL_HOMEOSTASIS", 
                                       "FATTY_ACID_METABOLISM", "GLYCOLYSIS", "HEME_METABOLISM", "OXIDATIVE_PHOSPHORYLATION", 
                                       "XENOBIOTIC_METABOLISM"), 
                       "Pathway" = c("APOPTOSIS", "HYPOXIA", 
                                     "PROTEIN_SECRETION", "REACTIVE_OXYGEN_SPECIES_PATHWAY", "UNFOLDED_PROTEIN_RESPONSE"
                       ), 
                       "Proliferation" = c("E2F_TARGETS", "G2M_CHECKPOINT", "MITOTIC_SPINDLE", 
                                           "MYC_TARGETS_V1", "MYC_TARGETS_V2", "P53_PATHWAY", "KEGG_DNA_REPLICATION", 
                                           "REACTOME_S_PHASE", "REACTOME_REPLICATION_STRESS", "REACTOME_CELL_CYCLE"
                       ), 
                       "Receptor Signaling" = c("ANDROGEN_RESPONSE", "ESTROGEN_RESPONSE_EARLY", 
                                                "ESTROGEN_RESPONSE_LATE"), 
                       "Signaling" = c("HEDGEHOG_SIGNALING", 
                                       "KRAS_SIGNALING_DN", "KRAS_SIGNALING_UP", "MTORC1_SIGNALING", 
                                       "NOTCH_SIGNALING", "PI3K_AKT_MTOR_SIGNALING", "TGF_BETA_SIGNALING", 
                                       "TNFA_SIGNALING_VIA_NFKB", "WNT_BETA_CATENIN_SIGNALING"))

pw_cats.merged <- make_category_table(pw_cats.merged, 'Pathway')


