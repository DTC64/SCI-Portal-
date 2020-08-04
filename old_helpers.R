
###### LR Functions for Shiny #######

cell <- list()
cell$Myeloid <- c('H-Microglia','DA-Microglia','Foamy-Microglia', 'Div-Microglia','Apop-Microglia',
                  'Monocyte','Trans-Monocyte','Foamy-Macrophage','Apoe-Macrophage', 'Div-Myeloid',
                  'U1-Myeloid','U2-Myeloid','U3-Myeloid','U4-Myeloid','U5-Myeloid',
                  'BA-Macrophage','Dendritic Cell',
                  'Neutrophil','Ph-Neutrophil','B Cell','T Cell')
cell$Vascular <- c('A-Endothelial','C-Endothelial','V-Endothelial',
                   'Pericyte','SMC',
                   'Fibroblast',
                   'IFN-Vascular', 'U1-Vascular','U2-Vascular')
cell$Macroglia <- c('H1-Ependymal','H2-Ependymal','Trans-Ependymal','Astroependymal','Astrocyte',
                    'Oligodendrocyte','OPC','Div-OPC','Pre-Oligo','Neuron')
times <- c('Uninjured','1dpi','3dpi','7dpi')

# Take subset of LR data.frame for plotting
selectLR <- function(
  score.results.df,
  timepoints = NULL,
  receptor.cell = NULL,
  receptor = NULL,
  ligand.cell = NULL,
  ligand = NULL,
  pair.name = NULL,
  significant.only = FALSE,
  organize = 'Ligand_Cell',
  sort.order = 'Ligand'
) {
  
  params <- list()
  params$data <- score.results.df
  # Functionality to add: DotPlot to add p-val info (low priority)
  params$style <- list('dims' = 1,
                       'p' = 'TilePlotLR')
  
  if(length(ligand.cell) == 0) {
    ligand.cell <- unlist(cell)
  } else if(any(ligand.cell %in% names(cell))) {
    ligand.cell <- union(ligand.cell[-which(ligand.cell %in% names(cell))],
                         unlist(cell[ligand.cell[ligand.cell %in% names(cell)]]))
  }
  if(length(receptor.cell) == 0) {
    receptor.cell <- unlist(cell)
  } else if(any(receptor.cell %in% names(cell))) {
    receptor.cell <- union(receptor.cell[-which(receptor.cell %in% names(cell))],
                         unlist(cell[receptor.cell[receptor.cell %in% names(cell)]]))
  }
  if(length(timepoints) == 0) {
    stop('No timepoint selected. Please choose a timepoint.')
  }
  
  if(all(c(length(pair.name), length(ligand), length(receptor)) == 0)) {
    pair.name.tmp <- levels(params$data$Pair_name)
  } else {
    if(length(ligand) == 0) {ligand <- levels(score.results.df$Ligand)}
    if(length(receptor) == 0) {receptor <- levels(score.results.df$Receptor)}
    pair.name.tmp <- union(pair.name, 
                           intersect(params$data$Pair_name[params$data$Ligand %in% ligand],
                                     params$data$Pair_name[params$data$Receptor %in% receptor]))
  }

  params$data <- params$data %>% 
    {if(significant.only) filter(., Pval < 0.05) else .} %>%
    filter(Ligand_Cell %in% ligand.cell) %>%
    filter(Receptor_Cell %in% receptor.cell) %>%
    filter(Pair_name %in% pair.name.tmp) %>%
    filter(Time %in% timepoints)
  
  if(!nrow(params$data) > 0) {stop('No interaction score available.')}

  params$data$Pair_name <- factor(params$data$Pair_name, levels = unique(params$data$Pair_name[order(params$data[[sort.order]])]))

  
  if(all(c(length(ligand.cell), length(receptor.cell)) > 1)) {
    params$style$dims <- 2
    params$x.axis <- 'Pair_name'
    params$y.axis <- c('Ligand_Cell','Receptor_Cell')[which(c('Ligand_Cell','Receptor_Cell') != organize)]
    params$y.facet <- organize
    params$switch <- list('y.position' = ifelse(organize == 'Ligand_Cell', 'right', 'left'),
                          'switch' = switch((organize == 'Ligand_Cell')+1, NULL, 'y'))
  } else if(length(ligand.cell) == 1) {
    params$x.axis <- 'Pair_name'
    params$y.axis <- 'Receptor_Cell'
    params$y.facet <- 'Ligand_Cell'
    params$switch <- list('y.position' = 'right',
                          'switch' = 'y')
  } else if(length(receptor.cell) == 1) {
    params$x.axis <- 'Pair_name'
    params$y.axis <- 'Ligand_Cell'
    params$y.facet <- 'Receptor_Cell'
    params$switch <- list('y.position' = 'left',
                          'switch' = NULL)
  }

  return(params)
}

# Dotplot of LR, where size of dot corresponds to log_Pval. 
# Only needed when {significant.only} == TRUE
# DotPlotLR <- function(
#   score.results.df
# ) {
#   m = median(score.results.df$Score)
#   lim_high = quantile(score.results.df$Score, 0.99)
#   lim_low = min(score.results.df$Score)
#   
#   score.results.df %>%
#     ggplot(mapping = aes(x = Pair_name, y = Ligand_Cell)) +
#     geom_point(aes_string(color = Score, size = log_Pval)) +
#     scale_color_gradient2(low = "blue3", high = "red3", midpoint = m, limits = c(lim_low, lim_high), na.value = 'red3') +
#     scale_size(range = c(0.1,4)) +
#     facet_grid(Receptor_Cell + Time ~ .) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 8),
#           axis.text.y = element_text(size = 8),
#           panel.background = element_rect(fill = "gray60"),
#           panel.grid.major = element_line(color = "black", size = 0.5)) -> tmp
#   return(tmp)
# }


TilePlotLR <- function(
  score.results.df,
  x.axis = 'Pair_name',
  y.axis = 'Receptor_Cell',
  y.facet = 'Ligand_Cell',
  switch = list('y.position' = 'right', 'switch' = 'y'),
  dims = 1
) {
  mid = floor(median(score.results.df$Score))
  hi = quantile(score.results.df$Score, 0.99)
  lo = min(score.results.df$Score)
  
  if(dims == 1) {
    y.facet.new <- paste(y.facet, 'Time', sep = '+')
  } else {
    plot.list <- list()
    for(t in unique(score.results.df$Time)) {
      plot.list[[t]] <- score.results.df %>%
        filter(Time %in% t) %>%
        ggplot(aes_string(y = y.axis, x = x.axis, fill = 'Score')) +
        geom_tile(color = 'black', size = 0.5) +
        labs(title = t) +
        scale_y_discrete(position = switch$y.position) +
        facet_grid(reformulate('.', y.facet), switch = switch$switch) +
        xlab(label = 'Ligand_Receptor Pairs') +
        ylab(label = y.axis) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.35),
              panel.grid.major.y = element_line(color = 'grey50'),
              panel.grid.major.x = element_line(color = 'grey85'),
              panel.background = element_rect(fill = 'grey85'),
              strip.text.y = element_text(size = 14),
              axis.text = element_text(size = 10),
              axis.title = element_text(size = 12),
              title = element_text(size = 20)) +
        scale_fill_gradient2(low = 'blue3', high = 'red3', na.value = 'red3',
                             midpoint = mid, limits = c(lo, hi))
    }
    return(plot_grid(plotlist = plot.list, ncol = 1))
  }

  score.results.df %>%
    ggplot(aes_string(y = y.axis, x = x.axis, fill = 'Score')) +
    geom_tile(color = 'black', size = 0.5) +
    facet_grid(reformulate('.', y.facet.new), switch = switch[[2]]) +
    scale_y_discrete(position = switch[[1]]) +
    xlab(label = 'Ligand_Receptor Pairs') +
    ylab(label = y.axis) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.35),
          panel.grid.major.y = element_line(color = 'grey50'),
          panel.grid.major.x = element_line(color = 'grey85'),
          panel.background = element_rect(fill = 'grey85'),
          strip.text.y = element_text(size = 14),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12)) +
    scale_fill_gradient2(low = 'blue3', high = 'red3', na.value = 'red3',
                         midpoint = mid, limits = c(lo, hi)) -> tmp
  
  return(tmp)
}


selectUMAP <- function(
  df,
  color.by = 'Cluster'
) {
  return(df %>% select(Time, Cluster, UMAP_1, UMAP_2, color.by))
}


UMAP_LR <- function(
  df,
  color.by = 'Cluster',
  split.time = FALSE
){
  p <- df %>%
    ggplot(mapping = aes_string(x = 'UMAP_1', y = 'UMAP_2', color = color.by)) +
    geom_point(size = 1.5) +
    theme(panel.background = element_rect(fill = 'grey95'),
          axis.line = element_line(size = 1),
          panel.grid = element_blank())
  if(!color.by %in% c('Cluster','Time')) {
    p <- p + scale_color_gradient(low = 'grey80', high = 'red1')
  } else if (color.by == 'Cluster') {
    p <- p + geom_text_repel(
      data = df %>% 
        group_by(Cluster) %>%
        mutate(median_umap_1 = median(UMAP_1)) %>%
        mutate(median_umap_2 = median(UMAP_2)) %>%
        select(median_umap_1, median_umap_2, Cluster) %>%
        distinct(),
      aes(x = median_umap_1, y = median_umap_2, label = Cluster),
      inherit.aes = FALSE,
      size = 4) +
      guides(fill = guide_legend(ncol = 2))
  }
  if(split.time) {
    p <- p + 
      facet_wrap(Time~., ncol = length(levels(df$Time))/2) +
      theme(strip.text = element_text(size = 16))
  }
  p.legend <- get_legend(p)
  p <- p + theme(legend.position = 'none', legend.text = element_text(size = 12), legend.title = element_text(size = 16))
  return(plot_grid(p, p.legend, ncol = 1, rel_heights = c(1, 0.3)))
}
  

# Precious scratch work...
# 
#   if(split.time) {
#     # Basic plot
#     p <- lapply(levels(df$Time), function(t) {
#       df %>%
#         filter(Time %in% t) %>%
#         ggplot(mapping = aes_string(x = 'UMAP_1', y = 'UMAP_2', color = color.by)) +
#         geom_point(size = 1) +
#         theme(panel.background = element_rect(fill = 'grey95'),
#               axis.line = element_line(size = 1),
#               panel.grid = element_blank(),
#               plot.title = element_text(size = 16)) +
#         labs(title = t)
#     })
#     # For gene expression levels
#     if(!color.by %in% c('Cluster', 'Time')) {
#       p <- lapply(1:length(p), function(t, p) {
#         p[[t]] + scale_color_viridis(
#           alpha = 0.8,
#           rescaler = function(x, to = c(0,1), from = NULL) {
#             scales::rescale(df[[color.by]], to = c(0,1))[df$Time %in% levels(df$Time)]
#           }
#         )
#       }, p)
#     } 
#     [df$Time %in% levels(df$Time)[t]]
#     # For cluster identity
#     if(color.by == 'Cluster') {
#       p <- lapply(1:length(p), function(t, p) {
#         p[[t]] + geom_text_repel(
#           data = df %>% 
#             filter(Time %in% levels(df$Time)[t]) %>%
#             mutate(median_umap_1 = median(UMAP_1)) %>%
#             mutate(median_umap_2 = median(UMAP_2)) %>%
#             select(median_umap_1, median_umap_2, Cluster) %>%
#             distinct(),
#           aes(x = median_umap_1, y = median_umap_2, label = Cluster),
#           inherit.aes = FALSE,
#           size = 4
#         )
#       }, p)
#     }
#     # Retrieve legend
#     # p.legend <- get_legend(plot = )
#   }
#   
#   return(plot_grid(plotlist = p))
# }
# 
# 
# 
# 
#   } else {
#     p[[1]] <- df %>%
#       ggplot(mapping = aes_string(x = 'UMAP_1', y = 'UMAP_2', color = color.by)) +
#       geom_point(size = 1) +
#       theme(panel.background = element_rect(fill = 'grey95'),
#             axis.line = element_line(size = 1),
#             panel.grid = element_blank()) +
#       labs('All Cells, All Times')
#   }
# 
# 
#   if(!color.by %in% c('Cluster','Time')) {
#     if(split.time) {
#       t <- names(p)
#       p <- lapply(seq_len(length(p)), function(x, p, t) {
#         p[[x]] <- p[[x]] + scale_color_viridis_c(
#           alpha = 0.8,
#           rescaler = function(x, to = c(0,1), from = NULL) {
#             scales::rescale(df[[color.by]], to = c(0,1)[df$Time %in% t])
#           }
#         )
#       }, p, t)
#     } else if(color.by == 'Cluster') {
#       l <- df %>%
#         group_by(Cluster) %>%
#         filter()
#         mutate(median_umap_1 = median(UMAP_1)) %>%
#         mutate(median_umap_2 = median(UMAP_2)) %>%
#         select(median_umap_1, median_umap_2, Cluster) %>%
#         distinct()
#       x + geom_text_repel(data = l,
#                           aes(x = median_umap_1, y = median_umap_2, label = Cluster),
#                           inherit.aes = FALSE,
#                           size = 4)
#     }
# 
#     p$legend <- p.legend
#   } else if (color.by == 'Cluster') {
#     p <- lapply(p, function(x) {
#       l <- df %>%
#         group_by(Cluster) %>%
#         mutate(median_umap_1 = median(UMAP_1)) %>%
#         mutate(median_umap_2 = median(UMAP_2)) %>%
#         select(median_umap_1, median_umap_2, Cluster) %>%
#         distinct()
#       x + geom_text_repel(data = l,
#                           aes(x = median_umap_1, y = median_umap_2, label = Cluster),
#                           inherit.aes = FALSE,
#                           size = 4)
#     })
#     if(split.time) {
#       p.legend <- get_legend(
#         plot = df %>%
#           filter(Time %in% t) %>%
#           ggplot(mapping = aes_string(x = 'UMAP_1', y = 'UMAP_2', color = color.by)) +
#           geom_point(size = 1) +
#           theme(panel.background = element_rect(fill = 'grey95'),
#                 axis.line = element_line(size = 1),
#                 panel.grid = element_blank())
#       )
#     }
#     p$legend <- p.legend
#   } else {
#     p.legend <- get_legend(
#       plot = df %>%
#         filter(Time %in% t) %>%
#         ggplot(mapping = aes_string(x = 'UMAP_1', y = 'UMAP_2', color = color.by)) +
#         geom_point(size = 1) +
#         theme(panel.background = element_rect(fill = 'grey95'),
#               axis.line = element_line(size = 1),
#               panel.grid = element_blank())
#     )
#     p$legend <- p.legend
#   }
# 
#   p.out <- plot_grid(plotlist = p[1:(length(p)-1)], ncol = ceiling((length(p)-1)/2))
#   p.out <- plot_grid(p.out, p$legend, ncol = 1, rel_heights = c(1,0.5))
# 
#   return(p.out)
# 
#   
#   
#   # 
#   # if(split.time) {
#   #   plot.list <- list()
#   #   for(t in levels(df$Time)) {
#   #     plot.list[[t]] <- df %>%
#   #       filter(Time %in% t) %>%
#   #       ggplot(mapping = aes_string(x = 'UMAP_1', y = 'UMAP_2', color = color.by)) +
#   #       geom_point(size = 1) +
#   #       theme(panel.background = element_rect(fill = 'grey95'),
#   #             axis.line = element_line(size = 1),
#   #             panel.grid = element_blank(),
#   #             legend.position = 'none')
#   #     if(!color.by %in% c('Cluster', 'Time')) {
#   #       tmp.legend <- get_legend(
#   #         plot = df %>%
#   #       )
#   #       plot.list[[t]] <- plot.list[[t]] +
#   #         scale_color_viridis_c(alpha = 0.75)
#   #     }
#   #   }
#   #   tmp.legend <- get_legend(
#   #     plot = ggplot(mapping = aes_string(x = 'UMAP_1', y = 'UMAP_2', color = color.by)) +
#   #       geom_point(size = 1) +
#   #       theme(panel.background = element_rect(fill = 'grey95'),
#   #             axis.line = element_line(size = 1),
#   #             panel.grid = element_blank()) + scale
#   #   )
#   #   
#   #   
#   #   return(plot_grid(plotlist = plot.list, ncol = 2))
#   # }
#   # 
#   # 
#   # scale_color_viridis_c(alpha = 0.75) +
#   #   
#   
# }
