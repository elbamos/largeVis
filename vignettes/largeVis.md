---
title: "largeVis"
author: "Amos Elberg"
date: "2016-05-19"
output: 
  rmarkdown::html_vignette:
    fig_caption: yes
  rmarkdown::pdf_document:
    fig_caption: yes
vignette: >
  %\VignetteIndexEntry{LargeVis}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


This Vingette provides an overview of the largeVis package.  

### Manifold Embeddings
Manifold embeddings attempt to find a low-dimensional (i.e., 2-D) representation of high-dimensional data that preserves the relationships among examples in the dataset. Examples that are similar to each other in high-dimensional space should appear close to each other in the low-dimensional space. 

LargeVis is intended to overcome a significant limitation of prior methods, that they are impractical for large datasets because they run in $O(n^2)$ or, at best, $O(n \log n)$. LargeVis should run in $O(n)$.

The neighboring nodes are used to estimate a probility distribution of each node's neighbors. 

Then, a lower dimensional embedding is estimated, first by random initialization.  Each node is assigned a position chosen at radnom. The embeddings are then updated to move similar nodes close together (and distant nodes further apart), according to a probabilstic functon that realates their distance in the low-dimensional space to the probability that the nodes are neighbors in the high-dimensional space. 

### Observations

t-SNE embedding has proven effecitve in wide variety of data contexts. 

However, the t-SNE algorithm depends on the identifying each node's nearest neighbors, an $O(N^2)$ operation at best. 

## Vignette Info


Note the various macros within the `vignette` section of the metadata block above. These are required in order to instruct R how to build the vignette. Note that you should change the `title` field and the `\VignetteIndexEntry` to match the title of your vignette.

## Styles

The `html_vignette` template includes a basic CSS theme. To override this theme you can specify your own CSS in the document metadata as follows:

    output: 
      rmarkdown::html_vignette:
        css: mystyles.css

## Effect of Hyperparameters

The following grid illustrates the effect of the $\alpha$ and $\gamma$ hyperparameters on the `wiki` dataset:


```
## Loading required package: largeVis
```

```
## Warning: Removed 11372 rows containing missing values (geom_point).
```

![plot of chunk wikihyperparameters](figure/wikihyperparameters-1.png)

($\alpha = 0$ uses the alternate probabilistic distance function, $p(e_{ij} = 1) = f(1 / (1 + \exp(||\vec(y_1) - \vec(y_2)||^2)))$).

## Figures




The figure sizes have been customised so that you can easily put two images side-by-side. 


```r
plot(1:10)
plot(10:1)
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png)![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-2.png)

You can enable figure captions by `fig_caption: yes` in YAML:

    output:
      rmarkdown::html_vignette:
        fig_caption: yes

Then you can use the chunk option `fig.cap = "Your figure caption."` in **knitr**.

## More Examples

You can write math expressions, e.g. $Y = X\beta + \epsilon$, footnotes^[A footnote here.], and tables, e.g. using `knitr::kable()`.


|                  |  mpg| cyl|  disp|  hp| drat|    wt|  qsec| vs| am| gear| carb|
|:-----------------|----:|---:|-----:|---:|----:|-----:|-----:|--:|--:|----:|----:|
|Mazda RX4         | 21.0|   6| 160.0| 110| 3.90| 2.620| 16.46|  0|  1|    4|    4|
|Mazda RX4 Wag     | 21.0|   6| 160.0| 110| 3.90| 2.875| 17.02|  0|  1|    4|    4|
|Datsun 710        | 22.8|   4| 108.0|  93| 3.85| 2.320| 18.61|  1|  1|    4|    1|
|Hornet 4 Drive    | 21.4|   6| 258.0| 110| 3.08| 3.215| 19.44|  1|  0|    3|    1|
|Hornet Sportabout | 18.7|   8| 360.0| 175| 3.15| 3.440| 17.02|  0|  0|    3|    2|
|Valiant           | 18.1|   6| 225.0| 105| 2.76| 3.460| 20.22|  1|  0|    3|    1|
|Duster 360        | 14.3|   8| 360.0| 245| 3.21| 3.570| 15.84|  0|  0|    3|    4|
|Merc 240D         | 24.4|   4| 146.7|  62| 3.69| 3.190| 20.00|  1|  0|    4|    2|
|Merc 230          | 22.8|   4| 140.8|  95| 3.92| 3.150| 22.90|  1|  0|    4|    2|
|Merc 280          | 19.2|   6| 167.6| 123| 3.92| 3.440| 18.30|  1|  0|    4|    4|

Also a quote using `>`:

> "He who gives up [code] safety for [code] speed deserves neither."
([via](https://twitter.com/hadleywickham/status/504368538874703872))
