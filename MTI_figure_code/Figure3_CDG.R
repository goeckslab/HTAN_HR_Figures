library(ggplot2)
library(tidyverse)
library(gridExtra)
library(cowplot)

### TO-DO: INGEST DATA VIA SYNAPSE CLIENT 
### FOR NOW, UNCOMMENT BELOW & SET TO PATH WITH FIGURE DATA 
# setwd('')

### --------------------------------------------------------------------------
### Code for generating Figure 3 Panels C, D, G, and Supp. Figure 1
### Figure 3C        : Broad phenotype density line plots
### Figure 3D        : Immune phenotype %CD45 line plots
### Figure 3G        : CD8 T cell tile densities, min. distances, and %GRZB+
### Supp. Figure 1   : Immune phenotype density line plots
### --------------------------------------------------------------------------

### Figure 3: Panel C -----------------------
data <- read.csv('mIHC_broad_phenotype_densities.csv')

# set factor levels/axis order 
data$patient <- factor(data$patient, levels = rev(c('9-1P','9-1A','9-2','9-3','9-14','9-15')))
data$broad_phenotype <- factor(data$broad_phenotype, levels = c('Neoplastic Epithelial','Immune','Other Stroma'))

# these columns allow for biopsy timepoint-specific density values while keeping long format data necessary for plots 
data <- data %>%
  select('patient','biopsy','broad_phenotype','density') %>%
  mutate(Bx1 = case_when(biopsy == 'Bx1' ~ density)) %>%
  mutate(Bx2 = case_when(biopsy == 'Bx2' ~ density))

# need a column that states whether the metric is increasing or decreasing from Bx1 to Bx2 
# pivot wider so can compare Bx1 to Bx2
data_wide <- data %>% 
  select('patient','biopsy','broad_phenotype','density') %>%
  pivot_wider(
    id_cols = c('patient','broad_phenotype'),
    names_from = 'biopsy',
    values_from = c('density')) %>%
  mutate(density_direction = case_when(
    Bx1 < Bx2 ~ 'Increase', 
    Bx1 > Bx2 ~ 'Decrease')) %>%
  select('patient','broad_phenotype','density_direction')

# merge the direction columns back into original dataframe 
data <- merge(data, data_wide, by = c('patient','broad_phenotype'))

# plot broad phenotype cell density line plots
p <- ggplot(data, aes(x = density, y = patient)) + 
  geom_line(aes(color = density_direction), linewidth = 1.5) +
  scale_color_manual(values=c('cyan3', 'magenta3')) +
  labs(color='\u0394 Bx1 to Bx2') +
  geom_point(color='black', shape=21, size=3, aes(fill=biopsy)) +
  scale_fill_manual(values=c('white', 'darkgrey')) +
  facet_wrap(.~broad_phenotype, nrow = 3, scales = 'free_x') +
  ylab('') +
  xlab('Cell Density') +
  labs(fill="Biopsy") +
  theme_bw() +
  guides(fill = guide_legend(order = 1), color = guide_legend(order = 2)) +
  theme(legend.position = "bottom", 
        legend.box = 'horizontal',
        legend.justification = c(1,0))

# save 
ggsave('Figure_3C.png', p, width = 5, height = 6, units = "in", device='png', dpi = 500)

### Clear environment -----------------------
rm(list = ls())

### Figure 3: Panel D -----------------------
data <- read.csv('mIHC_phenotype_composition.csv')

# set factor levels/axis order
data$patient <- factor(data$patient, levels = rev(c('9-1P','9-1A','9-2','9-3','9-14','9-15')))

# need a Percent column for the figure axis, then Bx1 and Bx2 specific percent columns 
data <- data %>%
  select('patient','biopsy','phenotype','percent_CD45','density') %>%
  mutate(Bx1 = case_when(biopsy == 'Bx1' ~ percent_CD45)) %>%
  mutate(Bx2 = case_when(biopsy == 'Bx2' ~ percent_CD45))

# need a column that states whether the metric is increasing or decreasing from Bx1 to Bx2 
# pivot wider so can compare Bx1 to Bx2
data_wide <- data %>% 
  pivot_wider(
    id_cols = c('patient','phenotype'),
    names_from = 'biopsy',
    values_from = c('percent_CD45','density')) %>%
  mutate(percent_CD45_direction = case_when(
    percent_CD45_Bx1 < percent_CD45_Bx2 ~ 'Increase', 
    percent_CD45_Bx1 > percent_CD45_Bx2 ~ 'Decrease')) %>%
  mutate(density_direction = case_when(
    density_Bx1 < density_Bx2 ~ 'Increase', 
    density_Bx1 > density_Bx2 ~ 'Decrease')) %>%
  select('patient','phenotype','percent_CD45_direction','density_direction')

# merge the direction columns back into original dataframe 
data <- merge(data, data_wide, by = c('patient','phenotype'))

# percent of CD45+ population immune line plot
p <- ggplot(data, aes(x = percent_CD45, y = patient)) + 
  geom_line(aes(color = percent_CD45_direction), linewidth = 1.5) +
  scale_color_manual(values=c('cyan3', 'magenta3')) +
  labs(color='\u0394 Bx1 to Bx2') +
  geom_point(color='black', shape=21, size=3, aes(fill=biopsy)) + 
  scale_fill_manual(values=c('white', 'darkgrey')) +
  facet_wrap(.~phenotype, nrow = 6, scales = 'free') + 
  ylab('') +
  xlab('% of CD45+ Population') +
  labs(fill = 'Biopsy') +
  theme_bw() +
  guides(fill = guide_legend(order = 1), color = guide_legend(order = 2)) +
  theme(legend.position = "bottom", 
        legend.box = 'horizontal',
        legend.justification = c(1,0))

ggsave('Figure_3D.png', p, width = 5, height = 12, units = "in", device='png')

# immune cell density line plot (supplemental)
p2 <- ggplot(data, aes(x = density, y = patient)) + 
  geom_line(aes(color = density_direction), linewidth = 1.5) +
  scale_color_manual(values=c('cyan3', 'magenta3')) +
  labs(color='\u0394 Bx1 to Bx2') +
  geom_point(color='black', shape=21, size=3, aes(fill=biopsy)) + 
  scale_fill_manual(values=c('white', 'darkgrey')) +
  facet_wrap(.~phenotype, nrow = 6, scales = 'free') + 
  ylab('') +
  xlab('Cell Density') +
  labs(fill = 'Biopsy') +
  theme_bw() +
  guides(fill = guide_legend(order = 1), color = guide_legend(order = 2)) +
  theme(legend.position = "bottom", 
        legend.box = 'horizontal',
        legend.justification = c(1,0))

ggsave('Supplemental_Figure_1.png', p2, width = 5, height = 12, units = "in", device='png')

### Clear environment -----------------------
rm(list = ls())

### Figure 3: Panel G -----------------------
dist_df <- read.csv('mIHC_cell_phenotype_min_distances.csv') # CD8 T cell min. distances
tile_df <- read_csv('mIHC_tile_densities.csv')               # CD8 T cell tile densities
grzb_df <- read_csv('mIHC_CD8_GRZB_percent.csv')             # CD8 T cell %GRZB+ 

# set factor levels/axis order
dist_df$biopsy <- factor(dist_df$biopsy, levels = c('Bx2','Bx1'))
tile_df$biopsy <- factor(tile_df$biopsy, levels = c('Bx2','Bx1'))
grzb_df$biopsy <- factor(grzb_df$biopsy, levels = c('Bx2','Bx1'))

dist_df$patient <- factor(dist_df$patient, levels = c('9-1A','9-2','9-3','9-14','9-15'))
tile_df$patient <- factor(tile_df$patient, levels = c('9-1A','9-2','9-3','9-14','9-15'))
grzb_df$patient <- factor(grzb_df$patient, levels = c('9-1A','9-2','9-3','9-14','9-15'))

# get descriptive statistics table for minimum distances
distance_summary_stats <- dist_df %>%
  # subset to phenotype of interest for individual plotting
  filter(query_cell_phenotype %in% c("CD8+ T cells")) %>%
  dplyr::group_by(patient, biopsy) %>%
  dplyr::summarise(mean_dist = mean(PANCK), sd_dist = sd(PANCK))

# get descriptive statistics table for tiles 
tile_summary_stats <- tile_df %>%
  # density threshold to exclude tiles with very few CD8 T cells 
  filter(CD8_density > 100) %>%
  dplyr::group_by(patient, biopsy) %>%
  dplyr::summarise(mean_CD8_density = mean(CD8_density), sd_CD8_density = sd(CD8_density))

# violin plot of minimum distances from CD8s to tumor cells 
distanceplot <- dist_df %>%
  # subset to phenotype of interest for individual plotting
  filter(query_cell_phenotype %in% c('CD8+ T cells')) %>%
  ggplot(aes(x = PANCK, y = biopsy, fill = biopsy)) +
  geom_violin(alpha = 0.4, scale = "width", draw_quantiles = c(0.25,0.75)) +
  scale_x_continuous(trans='log10') +
  scale_fill_discrete(limits = c('Bx1','Bx2')) +
  facet_wrap(.~patient, nrow= 5, scales = 'free_y', strip.position = 'left') +
  labs(color="Biopsy") +
  labs(fill="Biopsy") +
  ylab('') +
  xlab('Min. Distance to Tumor Cells') +
  theme_bw() +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(plot.margin=unit(c(1,0.1,1,0.1), "pt")) +
  stat_summary(fun.y="mean", geom="point", size = 2, alpha = 0.6) +
  guides(fill = guide_legend(override.aes = list(shape = NA)))

# violin plot of CD8 tile densities 
tileplot <- tile_df %>%
  # remove tiles below a threshold of CD8 T cell presence
  filter(CD8_density > 100) %>%
  ggplot(aes(x = CD8_density, y = biopsy, fill = biopsy)) +
  geom_violin(alpha = 0.4, scale = "width", draw_quantiles = c(0.25,0.75)) +
  scale_x_continuous(trans='log10') +
  scale_fill_discrete(limits = c('Bx1','Bx2')) +
  facet_wrap(.~patient, nrow= 5, scales = 'free_y', strip.position = 'left') +
  ylab('') +
  xlab('Per-tile cell density') +
  theme_bw() +
  theme(legend.position="none") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(plot.margin=unit(c(1,0.1,1,0.3), "pt")) +
  stat_summary(fun.y="mean", geom="point", size = 2, alpha = 0.6)

# bar plot of %GRZB+ 
grzbplot <- grzb_df %>%
  ggplot(aes(y = biopsy, x = CD8_percent_GRZB_pos, fill = biopsy)) +
  geom_col(alpha = 0.4, color = 'black') +
  scale_fill_discrete(limits = c('Bx1','Bx2')) +
  facet_wrap(.~patient, nrow = 5, strip.position = 'left') +
  ylab('') +
  xlab('% GRZB+') +
  theme_bw() +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(plot.margin=unit(c(1,1,1,0.1), "pt")) +
  theme(legend.position="none") +
  labs(fill="Biopsy")

# grab the legend from the min. distance plot and use for all three subplots
legend_b <- get_legend(distanceplot + theme(legend.position="bottom"))

# arrange plots on a single row
prow <- plot_grid( tileplot + theme(legend.position="none"),
                   distanceplot + theme(legend.position="none"),
                   grzbplot + theme(legend.position="none"),
                   nrow = 1,
                   rel_widths = c(1,1,0.6)
)

# combine plot row with legend
p <- plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1)) +
  theme(plot.background = element_rect(fill = "white", colour = NA)) +
  theme(plot.margin=unit(c(5,5,5,5), "pt"))

# save 
ggsave('Figure_3G.png', p, dpi = 500, width = 8, height = 5, units = "in", device='png')

### Clear environment -----------------------
rm(list = ls())

