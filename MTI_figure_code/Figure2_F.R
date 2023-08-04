library(ggplot2)
library(tidyverse)
library(gridExtra)
library(cowplot)

### TO-DO: INGEST DATA VIA SYNAPSE CLIENT 
### FOR NOW, UNCOMMENT BELOW & SET TO PATH WITH FIGURE DATA 
# setwd('')

### --------------------------------------------------------------------------
### Code for generating Figure 2 Panel F
### Figure 2F : Proliferating tumor Ripley's K and tile densities 
### --------------------------------------------------------------------------

# read proliferating tumor cell data: tile densities and WSI Ripley's K 
tile_data <- read.csv('cycIF_proliferating_tumor_tile_density.csv')
ripley_data <- read.csv('cycIF_proliferating_tumor_ripleys.csv')

# set factor levels 
tile_data$patient <- factor(tile_data$patient, levels = c('9-1P','9-2','9-3','9-14','9-15'))
ripley_data$patient <- factor(ripley_data$patient, levels = c('9-1P','9-2','9-3','9-14','9-15'))

tile_data$biopsy <- factor(tile_data$biopsy, levels = c('Bx1','Bx2'))
ripley_data$biopsy <- factor(ripley_data$biopsy, levels = c('Bx1','Bx2'))

# tile distribution violin plot
tileplot <- tile_data %>% 
  filter(prolif_tumor_tile_density > 0) %>%
  ggplot(aes(y = prolif_tumor_tile_density, x = biopsy, fill = biopsy)) +
  geom_violin(alpha = 0.4, draw_quantiles = c(0.25,0.75), scale = 'width') +
  scale_fill_discrete(limits = c('Bx1','Bx2')) +
  facet_wrap(.~patient, ncol= 5, scales = 'free_x', strip.position = 'bottom') +
  labs(color="Biopsy") +
  labs(fill="Biopsy") +
  xlab('') +
  ylab('Per-tile cell density') +
  theme_bw() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(plot.margin=unit(c(1,0.1,1,0.1), "pt")) +
  stat_summary(fun.y="mean", geom="point", size = 2, alpha = 0.6) +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  ylim(c(0,6500))

# plot of Ripley's K values for proliferating tumor cells at r = 50 microns 
ripleyplot <- ripley_data %>%
  ggplot(aes(x = biopsy, y = prolif_tumor_k_r50, fill = biopsy)) +
  geom_col(alpha = 0.4, color = 'black', position = 'dodge') +
  scale_fill_discrete(limits = c('Bx1','Bx2')) +
  facet_wrap(.~patient, ncol = 5, strip.position = 'bottom') +
  xlab('') +
  ylab("Ripley's K, R = 50") +
  theme_bw() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(plot.margin=unit(c(1,1,1,0.1), "pt")) +
  theme(legend.position="none") +
  labs(fill="Biopsy") +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

# grab the legend from the min. distance plot and use for all three subplots 
legend_b <- get_legend(tileplot + theme(legend.position="bottom"))

# arrange plots on a single row 
prow <- plot_grid( ripleyplot + theme(legend.position="none"),
                   tileplot + theme(legend.position="none"),
                   nrow = 2,
                   rel_heights = c(0.6,1)
)

# add a title 
title <- ggdraw() + 
  draw_label(
    "Spatial heterogeneity of proliferating tumor cells",
    x = 0.13,
    hjust = 0)

# combine plot row with legend and title 
p <- plot_grid(title, prow, legend_b, ncol = 1, rel_heights = c(.1, 1, .1)) +
  theme(plot.background = element_rect(fill = "white", colour = NA)) +
  theme(plot.margin=unit(c(5,5,5,10), "pt"))

# save
ggsave('Figure_2F.png', p, dpi = 500, width = 5, height = 6, units = "in", device='png')







