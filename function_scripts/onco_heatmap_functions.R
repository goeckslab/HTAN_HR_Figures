





variant_cols <- c('Gain' = "#E31A1C", 
                  'Loss' = "#1F78B4", 
                  "Frameshift Deletion" = "#A6CEE3", 
                  "Frameshift Insertion" = "#FFFF99", 
                  "In Frame Deletion" = "#B2DF8A", 
                  "In Frame Insertion" = "#6A3D9A", 
                  "Intron Substitution" = "#FB9A99", 
                  "Missense Substitution" = "#33A02C", 
                  "Nonsense Substitution" = "#FF7F00", 
                  "Promoter Substitution" = "#FDBF6F", 
                  "Splice Site Substitution" = "#CAB2D6", 
                  'SNV' = "black", 
                  'ND_CNV' = "grey65", 
                  'ND_SNV' = "grey45", 
                  'background' = "grey80")


onco_box_size <- 1
size.cnv <- unit(onco_box_size, 'pt')
size.snv <- .33


alter_fun = list(
  "Gain" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(onco_box_size, "pt"), h - size.cnv, 
              gp = gpar(fill = variant_cols["Gain"], col = NA))
  }, 
  "Loss" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(onco_box_size, "pt"), h - size.cnv, 
              gp = gpar(fill = variant_cols["Loss"], col = NA))
  }, 
  "Frameshift Deletion" = function (x, y, w, h) { 
    grid.rect(x, y, w - unit(onco_box_size, "pt"), h*size.snv, 
              gp = gpar(fill = variant_cols["Frameshift Deletion"], col = NA))
  }, 
  "Frameshift Insertion" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(onco_box_size, "pt"), h*size.snv, 
              gp = gpar(fill = variant_cols["Frameshift Insertion"], col = NA))
  }, 
  "In Frame Deletion" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(onco_box_size, "pt"), h - unit(onco_box_size,  "pt"), 
              gp = gpar(fill = variant_cols["In Frame Deletion"], col = NA))
  }, 
  "In Frame Insertion" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(onco_box_size, "pt"), h*size.snv, 
              gp = gpar(fill = variant_cols["In Frame Insertion"], col = NA))
  }, 
  "Intron Substitution" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(onco_box_size, "pt"), h*size.snv, 
              gp = gpar(fill = variant_cols["Intron Substitution"], col = NA))
  }, 
  "Missense Substitution" = function (x, y, w, h)  {
    grid.rect(x, y, w - unit(onco_box_size, "pt"), h*size.snv, 
              gp = gpar(fill = variant_cols["Missense Substitution"], col = NA))
  }, 
  "Nonsense Substitution" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(onco_box_size, "pt"), h*size.snv, 
              gp = gpar(fill = variant_cols["Nonsense Substitution"], col = NA))
  }, 
  "Promoter Substitution" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(onco_box_size, "pt"), h*size.snv, 
              gp = gpar(fill = variant_cols[v], col = NA))
  }, 
  "Splice Site Substitution" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(onco_box_size, "pt"), h*size.snv, 
              gp = gpar(fill = variant_cols["Splice Site Substitution"], col = NA))
  },
  "SNV" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(onco_box_size, "pt"), h*size.snv, 
              gp = gpar(fill = variant_cols["SNV"], col = NA))
  }, 
  "ND_CNV" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(onco_box_size, "pt"), h*size.snv, 
              gp = gpar(fill = variant_cols["ND_CNV"], col = NA))
  }, 
  "ND_SNV" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(onco_box_size, "pt"), h*size.snv, 
              gp = gpar(fill = variant_cols["ND_SNV"], col = NA))
  }, 
  "background" = function (x, y, w, h) {
    grid.rect(x, y, w - unit(onco_box_size, "pt"), h - unit(onco_box_size, 'pt'), 
              gp = gpar(fill = variant_cols["background"], col = NA))
  }
)



