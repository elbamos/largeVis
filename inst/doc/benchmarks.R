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
                      comment = "#>")
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
rebuild <- FALSE
if (!exists("buildManifolds")) buildManifolds <- rebuild

## ----performance,echo=F,eval=rebuild-------------------------------------
#  # benchmark <- readr::read_csv(system.file("extdata", "results.csv", package="largeVis"),
#  #                              col_names = FALSE)
#  benchmark1 <- readr::read_csv("../inst/results.csv", col_names = FALSE)
#  benchmark1$machine <- 1
#  benchmark2 <- readr::read_csv("../inst/nelsonresults.csv", col_names = FALSE)
#  benchmark2$machine <- 2
#  
#  benchmark <- rbind(benchmark1, benchmark2, benchmark3) %>%
#    set_colnames(c("time",
#      "precision",
#      "n_trees",
#      "max_iters",
#      "threshold",
#      "method",
#      "tree_type",
#      "search_type",
#      "eps",
#      "machine")
#      ) %>%
#    select(-tree_type, -search_type, - eps) %>%
#    mutate(method = ifelse(grepl("largeVis", method), "largeVis", method),
#           series = ifelse(grepl("largeVis", method),
#                           paste(method, max_iters, "iter."),
#                           method),
#           series = factor(series),
#           method = factor(method)) %>%
#    group_by(machine) %>%
#    mutate(time = time / min(time))
#  
#  #benchmark[(benchmark$n_trees == 10 & benchmark$method == 'RcppAnnoy'),]
#  
#  # benchmark$time <- benchmark$time / ifelse(benchmark$machine == 1, 1.444920, 1.535254)
#  

## ----plotpeformance,echo=F,fig.width=6,fig.height=6,fig.align='center',warning=FALSE,message=FALSE----
load(system.file("extdata", "benchmark.Rda", package = "largeVis"))
benchmark %>% 
  dplyr::filter(series != 'RcppAnnoy') %>%
  mutate(facet = precision / 100, 
         facet = ifelse(facet < 0.9, '', 'Closeup'), 
         facet = factor(facet)) %>%
  ggplot(aes(y = 1 / time, 
                      x = (precision / 100), 
                      group = series, 
                      color = series, 
                      shape = method)) +
  geom_point(size = 1.5, alpha = 0.8) + 
 # geom_line(size = 0.5) + 
  # geom_text(size = 3, 
  #           nudge_x = 10) +
  scale_y_log10(name = "Time, log (speed relative to RcppAnnoy with 10 trees)") + 
  scale_x_continuous("Precision", 
               # limits = c(0,1), 
                breaks = c(0, 0.2, 0.4, 0.6, 0.8, 0.925, 0.95, 0.975, 1.0)) +
  facet_grid(. ~ facet, scales = "free_x") +
  scale_color_manual(values =      colors_divergent_discrete(nlevels(benchmark$series))(nlevels(benchmark$series))) +
  guides(color = guide_legend(nrow=3)) +
  ggtitle(expression(
    atop("Time vs. Error Rate (K = 100, n = 10000)",
         atop(italic("(Upper Right is Better)"))
         )
    ))

## ----n_trees,echo=F,fig.width=5,fig.height=7-----------------------------
bench <- benchmark %>%
  filter(method == 'largeVis') %>%
  select(-series) %>%
  mutate(seriesa = paste(machine, method, threshold, max_iters))
bench$series <- factor(bench$seriesa)
bench %>% 
  group_by(series) %>% 
  filter(n() > 1, 
         max_iters < 2) %>%
  arrange(n_trees) %>%
  ungroup() %>%
  ggplot(aes(y = 1 / time, 
             x = (precision / 100), 
             color = series, 
             group = series)) + 
  geom_path(size = 0.5, 
            arrow = arrow(length = unit(0.05, "inches"))) +
  geom_point(size = 0.5) +
  facet_grid(threshold ~ max_iters) +
  scale_y_log10(name = "Time, log (speed relative to RcppAnnoy with 10 trees)") + 
  scale_x_continuous("Precision", 
                limits = c(0,1), 
                breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  guides(color = FALSE) +
  # scale_color_manual(values =      colors_divergent_discrete(nlevels(benchmark$series))(nlevels(benchmark$series))) +
  ggtitle(expression(
    atop("Time vs. Precision Rate (K = 100, n = 10000)",
         atop(italic("(Upper Right is Better)"))
         )
    ))

## ----tree_threshold,echo=F-----------------------------------------------
bench <- benchmark %>%
  filter(method == 'largeVis') %>%
  select(-series) %>%
  mutate(seriesa = paste(machine, n_trees, max_iters))
bench$series <- factor(bench$seriesa)
bench %<>%
  select(-n_trees, -max_iters, -method, -machine) %>%
  group_by(series) %>% 
  filter(n() > 1) %>%
  mutate(label = ifelse(threshold == 128, "128", "Other"), 
         label = factor(label), 
         facet = precision / 100, 
         facet = ifelse(facet < 0.9, '', 'Closeup'))
bench$facet <- factor(bench$facet)
bench %>% 
  arrange(threshold) %>%
  ggplot(aes(y = 1/time, 
             x = precision / 100, 
             color = series, 
             group = series, 
             shape = label)) + 
  geom_path(size = 0.5, alpha =0.8, arrow = arrow(length = unit(0.05, "inches"))) +
  geom_point(size = 1.5) +
  facet_grid(. ~ facet, scales = 'free_x') +
  scale_y_log10(name = "Time, log (speed relative to RcppAnnoy with 10 trees)") + 
  scale_x_continuous("Precision", 
                breaks = c(0, 0.2, 0.4, 0.6, 0.8, 0.925, 0.95, 0.975, 1.0)) +
  scale_shape_discrete(name = "", solid = FALSE) + 
    guides(color = FALSE) +
  # scale_color_manual(values =      colors_divergent_discrete(nlevels(benchmark$series))(nlevels(benchmark$series))) +
  ggtitle(expression(
    atop("Time vs. Precision (K = 100, n = 10000)",
         atop(italic("(Upper Right is Better)"))
         )
    ))

## ----max_iters,echo=F----------------------------------------------------
bench <- benchmark %>%
  filter(method == 'largeVis') %>%
  select(-series) %>%
  mutate(series = paste(machine, n_trees, threshold))
bench$series <- factor(bench$series)
bench %>%
  select(-n_trees, -threshold, -method, -machine) %>%
  group_by(series) %>% 
  filter(n() > 1) %>%
  mutate(facet = precision / 100, 
         facet = ifelse(facet < 0.9, '', 'Closeup'), 
         facet = factor(facet), 
         shape = ifelse(max_iters == 1, '1', 'Other'), 
         shape = ifelse(max_iters == 0, '0', shape), 
         shape = ifelse(max_iters > 1, '> 1', shape), 
         shape = factor(shape)) %>%
  arrange(max_iters) %>%
  ggplot(aes(y = 1 / time, 
             x = precision / 100, 
             color = series, 
             group = series,
             shape = shape)) + 
  geom_line(size = 0.2, arrow = arrow(length = unit(0.05, "inches"))) +
  geom_point(size = 1.5, alpha = 0.8) +
  facet_grid(. ~ facet, scales = 'free_x') +
  guides(color = FALSE) +
  scale_shape(name = "Iterations", solid = FALSE) +
  scale_y_log10(name = "Time, log (speed relative to RcppAnnoy with 10 trees)") + 
  scale_x_continuous("Precision", 
                breaks = c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.975, 1.0)) +  ggtitle(expression(
    atop("Time vs. Precision (K = 100, n = 10000)",
         atop(italic("(Upper Right is Better)"))
         )
    ))

