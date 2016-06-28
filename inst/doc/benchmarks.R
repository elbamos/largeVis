## ----setupbenchmark,eval=T,echo=F,warning=F,error=F,message=F------------
# Note to reader:  Please don't steal the semi-distinctive visual style I spent several minutes creating for myself.
require(ggplot2, 
        quietly = TRUE)
require(RColorBrewer, 
        quietly = TRUE)
require(wesanderson, 
        quietly = TRUE)
require(dplyr, quietly = TRUE)
knitr::opts_chunk$set(collapse = TRUE, 
                      comment = "#>",
                      fig.width = 7, 
                      fig.height = 5)
colors_discrete <- function(x) rep(wes_palette("Darjeeling", 
                                               n = min(x, 5)), 
                                   2)[1:x]
colors_divergent_discrete <- function(x) 
  grDevices::colorRampPalette(RColorBrewer::brewer.pal(x, "Spectral"))
colors_continuous <-  function(x) wes_palette(name = "Zissou",
                                              n = x, 
                                              type = "continuous")

nacol <- colors_discrete(4)[4]
theme_set(
  theme_bw() %+replace%
  theme(
    legend.key.size = unit(4, "mm"), 
    legend.title = element_text(size = rel(0.8),
                              face = "bold"),
    legend.margin = unit(0, "cm"),
    legend.position = "bottom",
    legend.key.size = unit(0.5, "lines"),
    legend.text=element_text(size = unit(8, "points")), 
    axis.title.y = element_text(angle = 90),
    axis.text = element_text(size = rel(0.7)),
    plot.margin = unit(c(0, 0.5, 1, 0), "lines"), 
    axis.title = element_text(size = rel(0.8),
                              face = "bold"),
    title = element_text(size = rel(0.9))
  ) 
)

## ----plotpeformance,echo=F,fig.align='center',warning=FALSE,message=FALSE----
load(system.file("extdata", "benchmark.Rda", package = "largeVis"))
benchmark %>% 
  filter(machine != 'Large Server',
         machine == 'Workstation' | K == 50) %>%
  mutate(facet = precision, 
         facet = ifelse(facet < 0.95, '', 'Closeup'), 
         facet = factor(facet)) %>%
  ggplot(aes( y = time, 
              x = precision, 
              group = series, 
              fill = series, 
              shape = series)) +
  geom_point(size = 1.5, alpha = 0.7, color = "grey80") + 
  scale_y_log10(name = "Speed, log (nodes / seconds)") + 
  scale_x_continuous("Precision", 
                breaks = c(0, 0.2, 0.4, 0.6, 0.8, 0.925, 0.95, 0.975, 1.0)) +
  facet_grid(K + machine ~ facet, scales = "free") +
  scale_fill_manual(name = "Method & n. iter.", 
    values = colors_divergent_discrete(nlevels(benchmark$series))(nlevels(benchmark$series))) +
  scale_shape_manual(name = "Method & n. iter.", 
                     values = c(21, 21, 21, 21, 23)) +
 # guides(color = guide_legend(nrow=3)) +
  ggtitle(expression(
    atop("Precision-Performance tradeoff, RcppAnnoy and largeVis",
         atop(italic("(n = 10000; Upper Right is Better)"))
         )
    ))

## ----constn,echo=F,warning=F---------------------------------------------
bench <- benchmark %>%
  filter(method == 'largeVis', machine == 'Large Server') %>%
  mutate(nn = threshold * n_trees) %>% 
  group_by(max_iters, nn)  %>% 
  filter(n() > 2) %>% 
  mutate(series = paste(max_iters, ", ", nn, sep = " "))
bench$facet <- factor(ifelse(bench$n_trees >= 4, "", "n. trees < 10"))
bench %>%
  ggplot(aes(y = time, 
           x = precision, 
           fill = series, 
           group = series,
           color = factor(n_trees))) + 
  geom_point(size = 1.5, alpha = 0.8, shape = 21) +
  scale_fill_manual("n. iter, tth", values = colors_divergent_discrete(6)(6)) +
  scale_color_grey("n. trees", start = 0.8, end = 0 ) +
#  guides(color = FALSE) +
 # scale_shape(name = "Iterations", solid = FALSE) +
  facet_grid(machine ~ .) +
  scale_y_log10(name = "Speed, log (nodes / second)", limits = c(1e2,1e5)) + 
  scale_x_continuous("Precision", 
                breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +  
  ggtitle(expression(
    atop("Precision-Performance tradeoff, n_trees and tree_threshold",
         atop(italic("(100-NN precision, n = 10000; Upper Right is Better)")))))

## ----tree_threshold,echo=F-----------------------------------------------
bench <- benchmark %>%
  filter(method == 'largeVis',
         machine != 'Large Server') %>%
  mutate(label = ifelse(threshold == 128, "128", "Other"), 
         label = factor(label), 
         facet = precision, 
         facet = ifelse(facet < 0.85, '', 'Closeup'))
bench$facet <- factor(bench$facet)
bench %>% 
  arrange(nn) %>%
  mutate(max_iters = factor(max_iters)) %>%
  ggplot(aes(y = time, 
             x = precision , 
             color = max_iters, 
             group = max_iters)) + 
#  geom_path(size = 0.5, alpha =0.8, arrow = arrow(length = unit(0.05, "inches"))) +
  geom_point(size = 1, alpha = 0.8, shape = 16) +
  facet_grid(K + machine ~ facet, scales = 'free') +
  scale_y_log10(name = "Speed, log (nodes / second)") + 
  scale_x_continuous("Precision", 
                breaks = c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0)) +
 # scale_shape_discrete(name = "", solid = FALSE) + 
 #   guides(color = FALSE) +
   scale_color_manual("n. iter", values = colors_discrete(4)) +    
  ggtitle(expression(
    atop("Precision-Performance tradeoff, effect of increasing tth vs. max_iters",
         atop(italic("(n = 10000; Upper Right is Better)")))))

