###plotDAbeeswarm

####plotDAbeeswarm
library(ggbeeswarm)
alpha = 0.05
group.by = "major_cluster"
da.res <- mutate(da_results, group_by = da_results[, group.by ])
da.res <- mutate(da.res, group_by = factor(group_by, levels = unique(group_by)))
head(da.res)
beeswarm_pos <- ggplot_build(da.res %>% mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) 
                             %>% arrange(group_by) %>% ggplot(aes(group_by, logFC)) + 
                               geom_quasirandom()
                             ) # build the beeswarm
pos_x <- beeswarm_pos$data[[1]]$x # get the pos from beeswarm_pos
pos_y <- beeswarm_pos$data[[1]]$y # get the pos from beeswarm_pos
n_groups <- unique(da.res$group_by) %>% length()

###########
#plot the figure
###########
da.res %>% mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>% 
  mutate(logFC_color = ifelse(is_signif == 1, logFC, NA)) %>% 
  arrange(group_by) %>% mutate(Nhood = factor(Nhood, levels = unique(Nhood))) %>% 
  mutate(pos_x = pos_x, pos_y = pos_y) %>% 
  ggplot(aes(pos_x, pos_y, color = logFC_color)) + scale_color_gradient2() + 
  guides(color = "none") + xlab(group.by) + ylab("Log Fold Change") + 
  scale_x_continuous(breaks = seq(1, n_groups), 
                     labels = setNames(levels(da.res$group_by), seq(1, n_groups))) + 
  geom_point() + coord_flip() + 
  theme_bw(base_size = 22) + theme(strip.text.y = element_text(angle = 0))
