library(ggplot2)
library(tidyverse)
library(gridExtra)
library(cowplot)

### TO-DO: INGEST DATA VIA SYNAPSE CLIENT 
### FOR NOW, UNCOMMENT BELOW & SET TO PATH WITH FIGURE DATA 
# setwd('')

### --------------------------------------------------------------------------
### Code for generating Figure 4 Panels A, C, E, F, G
### Figure 4A : CD163+ Macrophage tile density and minimum distance plots
### Figure 4C : SpatialScore violin plot
### Figure 4E : %PD-L1+ by phenotype line plot
### Figure 4F : PD1 to PDL1 minimum distance violin plot
### Figure 4G : CD8 differentiation (PD1/EOMES) bubble plot 
### --------------------------------------------------------------------------


### Figure 4A ----------------------
dist_df <- read.csv('mIHC_cell_phenotype_min_distances.csv')
tile_df <- read_csv('mIHC_tile_densities.csv')

# set axis orders 
dist_df$biopsy <- factor(dist_df$biopsy, levels = c('Bx2','Bx1'))
tile_df$biopsy <- factor(tile_df$biopsy, levels = c('Bx2','Bx1'))

dist_df$patient <- factor(dist_df$patient, levels = c('9-1A','9-2','9-3','9-14','9-15'))
tile_df$patient <- factor(tile_df$patient, levels = c('9-1A','9-2','9-3','9-14','9-15'))

# get descriptive statistics table for minimum distances
distance_summary_stats <- dist_df %>% 
  # subset to phenotype of interest for individual plotting 
  filter(query_cell_phenotype %in% c("CD163+ Macrophages/Monocytes")) %>% 
  dplyr::group_by(patient, biopsy) %>%
  dplyr::summarise(mean_dist = mean(PANCK), sd_dist = sd(PANCK))

# get descriptive statistics table for tile densities 
tile_summary_stats <- tile_df %>% 
  filter(CD163_density > 100) %>%
  dplyr::group_by(patient, biopsy) %>%
  dplyr::summarise(mean_CD163_density = mean(CD163_density), sd_CD163_density = sd(CD163_density))

# violin plot of minimum distances between CD163+ cells and tumor cells 
distanceplot <- dist_df %>% 
  # subset to phenotype of interest for individual plotting 
  filter(query_cell_phenotype %in% c("CD163+ Macrophages/Monocytes")) %>% 
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

# violin plot of CD163+ tile densities 
tileplot <- tile_df %>% 
  # remove tiles below a threshold of CD163 T cell presence 
  filter(CD163_density > 100) %>%
  ggplot(aes(x = CD163_density, y = biopsy, fill = biopsy)) +
  geom_violin(alpha = 0.4, scale = "width", draw_quantiles = c(0.25,0.75)) +
  scale_x_continuous(trans='log10') +
  scale_fill_discrete(limits = c('Bx1','Bx2')) +
  facet_wrap(.~patient, nrow= 5, scales = 'free_y', strip.position = 'left') +
  ylab('') +
  xlab('Per-tile cell density') +
  theme_bw() +
  theme(legend.position="none") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(plot.margin=unit(c(1,5,1,0.3), "pt")) +
  stat_summary(fun.y="mean", geom="point", size = 2, alpha = 0.6) +
  guides(fill = guide_legend(override.aes = list(shape = NA)))

# grab the legend from the min. distance plot and use for all three subplots 
legend_b <- get_legend(distanceplot + theme(legend.position="bottom"))

# arrange plots on a single row 
prow <- plot_grid( tileplot + theme(legend.position="none"),
                   distanceplot + theme(legend.position="none"),
                   nrow = 1,
                   rel_widths = c(1,1)
)

# combine plot row with legend 
p <- plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1)) +
  theme(plot.background = element_rect(fill = "white", colour = NA)) +
  theme(plot.margin=unit(c(5,5,5,5), "pt"))

# save
ggsave('Figure_4A.png', p, dpi = 500, width = 8, height = 5, units = "in", device='png')


# FIGURE 4C USES SAME INPUT FILE AS 4A SO DON'T CLEAR ENVIRONMENT HERE 

### Figure 4C ----------------------
# calculate the SpatialScore as the ratio of min dist(Helper T cells - Tumor) to min dist (Helper T - T-reg)
dist_df <- dist_df %>%
  dplyr::filter(query_cell_phenotype == "Helper T cells") %>%
  mutate(SpatialScore = PANCK / Regulatory)

# violin plot of the SpatialScore
p2 <- dist_df %>%
  ggplot(aes(y = biopsy, x = SpatialScore, fill = biopsy)) +
  geom_violin(alpha = 0.4, scale = "width", draw_quantiles = c(0.25,0.75)) +
  scale_fill_discrete(limits = c('Bx1','Bx2')) +
  scale_x_continuous(trans='log10') +
  facet_wrap(.~patient, nrow= 5, scales = 'free_y', strip.position = 'left') +
  ylab('') +
  xlab(expression(italic('SpatialScore'))) +
  theme_bw() +
  theme(legend.position="bottom") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  stat_summary(fun.y="mean", geom="point", size = 2, alpha = 0.6) +
  guides(fill = guide_legend(override.aes = list(shape = NA)))

# save 
ggsave('Figure_4C.png', p2, dpi = 500, width = 5, height = 6, units = "in", device='png')

### Clear environment -----------------------
rm(list = ls())

### Figure 4E ----------------------
# read data
data <- read.csv('mIHC_PDL1_phenotype_percentages.csv')

# set factor levels/axis order 
data$biopsy <- factor(data$biopsy, levels = c('Bx1','Bx2'))
data$patient <- factor(data$patient, levels = rev(c("9-1P", "9-1A", "9-2",  "9-3",  "9-14", "9-15")))
data$phenotype <- factor(data$phenotype, levels = c('Tumor','DC','CD163 Mac/Mono'))

# need a percent column for axis and Bx1/Bx2 specific columns 
data <- data %>%
  select('patient','biopsy','phenotype','percent_pos_PDL1') %>%
  mutate(Bx1 = case_when(biopsy == 'Bx1' ~ percent_pos_PDL1)) %>%
  mutate(Bx2 = case_when(biopsy == 'Bx2' ~ percent_pos_PDL1))

# need a column that states whether the metric is increasing or decreasing from Bx1 to Bx2 
# pivot wider so can compare Bx1 to Bx2
data_wide <- data %>% 
  select('patient','biopsy','phenotype','percent_pos_PDL1') %>%
  pivot_wider(
    id_cols = c('patient','phenotype'),
    names_from = 'biopsy',
    values_from = c('percent_pos_PDL1')) %>%
  mutate(percent_direction = case_when(
    Bx1 < Bx2 ~ 'Increase',
    Bx1 > Bx2 ~ 'Decrease')) %>%
  select('patient','phenotype','percent_direction')

# merge the direction columns back into original dataframe
data <- merge(data, data_wide, by = c('patient','phenotype'))

# plot
p <- ggplot(data, aes(x = percent_pos_PDL1, y = patient)) + 
  geom_line(aes(color = percent_direction), linewidth = 1.5) +
  scale_color_manual(values=c('cyan3', 'magenta3')) +
  labs(color='\u0394 Bx1 to Bx2') +
  geom_point(color='black', shape=21, size=3, aes(fill=biopsy)) + 
  scale_fill_manual(values=c('white', 'darkgrey')) +
  facet_wrap(.~phenotype, nrow = 3) + 
  ylab('') +
  xlab('% PDL1+ Cells') +
  labs(fill = 'Biopsy') +
  theme_bw() +
  guides(fill = guide_legend(order = 1), color = guide_legend(order = 2)) +
  theme(legend.position = "bottom", 
        legend.box = 'horizontal',
        legend.justification = c(1,0))

ggsave('Figure_4E.png', p, width = 5, height = 7, units = "in", device='png')

### Clear environment -----------------------
rm(list = ls())

### Figure 4F ----------------------
# read data
dist_df <- read.csv('mIHC_PD1_PDL1_min_distance.csv')

# set axis orders
dist_df$patient <- factor(dist_df$patient, levels = c('9-1A','9-2','9-3','9-14','9-15'))
dist_df$biopsy <- factor(dist_df$biopsy, levels = c('Bx2','Bx1'))

# get descriptive stats 
distance_summary_stats <- dist_df %>%
  dplyr::group_by(patient, biopsy) %>%
  dplyr::summarise(mean_dist = mean(min_dist_to_PDL1), sd_dist = sd(min_dist_to_PDL1))

# violin plot of PD1-PDL1 minimum distances 
p <- dist_df %>%
  ggplot(aes(x = min_dist_to_PDL1, y = biopsy, fill = biopsy)) +
  geom_violin(alpha = 0.4, scale = "width", draw_quantiles = c(0.25,0.75)) +
  scale_fill_discrete(limits = c('Bx1','Bx2'), name = 'Biopsy') +
  scale_x_continuous(trans='log10') +
  ylab('') +
  xlab('PD1+ cells: Minimum distance to PD-L1+ cells') +
  facet_wrap(.~patient, nrow= 6, scales = 'free_y', strip.position = 'left') +
  theme_bw() +
  theme(legend.position="bottom") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  stat_summary(fun.y="mean", geom="point", size = 2, alpha = 0.6) +
  guides(fill = guide_legend(override.aes = list(shape = NA)))

ggsave('Figure_4F.png', p, width = 4, height = 5, units = 'in', dpi = 500)

### Clear environment -----------------------
rm(list = ls())

### Figure 4G ----------------------
# Read CD8 T cell differentiation (PD1/EOMES) data 
data <- read.csv('mIHC_CD8_PD1_EOMES.csv')

# Sample order for plotting 
data$sample_ID <- factor(data$sample_ID, 
                         levels = rev(c("9-1P Bx1", "9-1P Bx2",
                                        "9-1A Bx1", "9-1A Bx2",
                                        "9-2 Bx1",  "9-2 Bx2",
                                        "9-3 Bx1",  "9-3 Bx2",
                                        "9-14 Bx1", "9-14 Bx2",
                                        "9-15 Bx1", "9-15 Bx2")))

p <- data %>%
  # bubble plot where cell type density is color, bubble size is proportion 
  ggplot(aes(y=sample_ID, x=phenotype, size=proportion, color=density)) +
  geom_point(alpha=0.8) +
  # setting custom colorscale to match manuscript heatmaps 
  scale_colour_gradient2(
    low = 'cyan',
    mid = "black",
    high = 'magenta',
    midpoint = 1,
    trans = 'log10') + 
  # setting bubble size range 
  scale_size(range = c(.1, 10), name="% of total \nCD8+ T cells") +
  theme_bw() +
  ylab('') +
  xlab('') +
  theme(axis.text.x = element_blank()) +               
  theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  # All the rest of this is formatting the axis to be a +/- matrix for PD1 and EOMES
  coord_cartesian(clip = "off", ylim = c(1,12)) +
  annotate(geom = "text", x = 1:4, y = 0, label = c('-','-','+','+'), vjust = 1.3, size = 6) +
  annotate(geom = "text", x = 0, y = 0, label = 'PD1 -', vjust = 1.8, hjust = 1, size = 4) +
  annotate(geom = "text", x = 1:4, y = 0, label = c('-','+','-','+'), vjust = 2.9, size = 6) +
  annotate(geom = "text", x = 0, y = 0, label = 'EOMES -', vjust = 4, hjust = 1, size = 4) +
  annotate("rect", xmin = 0, xmax = 5, ymin = -1.8, ymax = 0, colour="black", alpha=0)

# save; above annotation coordinates work for these W x H dimensions specifically 
ggsave('Figure_4G.png', p, width = 5, height = 5, units = "in", device='png', dpi = 500)

### Clear environment -----------------------
rm(list = ls())
























