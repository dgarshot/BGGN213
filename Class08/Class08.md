Class 08\_Principal Component Analysis(PCA)
================
Danielle Garshott
2/6/2019

### Clustering

``` r
# Generate some example data for clustering
tmp <- c(rnorm(30,-3), rnorm(30,3))
x <- cbind(x=tmp, y=rev(tmp))
plot(x)
```

![](Class08_files/figure-markdown_github/unnamed-chunk-1-1.png)

``` r
# k means
# kmeans(x, center=3, nstart = 20)
# nstart (how many times dou you want to start this calculation over)
# Euclidean calculation d= sqr(x^2 + y^2 + z^2)
```

tmp &lt;- c(rnorm(30,-3), rnorm(30,3)) -&gt; 30 points centered around -3, and 30 around 3 (60 points total)
============================================================================================================

cbind -&gt; putting the columns together
========================================

``` r
km <- kmeans(x, 2, nstart = 20)
km
```

    ## K-means clustering with 2 clusters of sizes 30, 30
    ## 
    ## Cluster means:
    ##           x         y
    ## 1 -2.731096  2.858829
    ## 2  2.858829 -2.731096
    ## 
    ## Clustering vector:
    ##  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2
    ## [36] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
    ## 
    ## Within cluster sum of squares by cluster:
    ## [1] 71.99792 71.99792
    ##  (between_SS / total_SS =  86.7 %)
    ## 
    ## Available components:
    ## 
    ## [1] "cluster"      "centers"      "totss"        "withinss"    
    ## [5] "tot.withinss" "betweenss"    "size"         "iter"        
    ## [9] "ifault"

``` r
km$size
```

    ## [1] 30 30

``` r
km$cluster
```

    ##  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2
    ## [36] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

``` r
plot(x, col=km$cluster)
points(km$centers, col="blue", pch=15, cex=2)
```

![](Class08_files/figure-markdown_github/unnamed-chunk-5-1.png)

Hierarchical clustering
-----------------------

``` r
# First we need to calculate point (dis)similarity
# as the Euclidean distance between observations
dist_matrix <- dist(x)
# The hclust() function returns a hierarchical
# clustering model
hc <- hclust(d = dist_matrix)
# the print method is not so useful here
hc 
```

    ## 
    ## Call:
    ## hclust(d = dist_matrix)
    ## 
    ## Cluster method   : complete 
    ## Distance         : euclidean 
    ## Number of objects: 60

``` r
plot(hc)
abline(h=6, col="red")
```

![](Class08_files/figure-markdown_github/unnamed-chunk-7-1.png)

``` r
grp2 <- cutree(hc, h=6)
```

``` r
plot(x, col=grp2)
```

![](Class08_files/figure-markdown_github/unnamed-chunk-8-1.png)

``` r
plot(hc)
abline(h=2.5, col="blue")
```

![](Class08_files/figure-markdown_github/unnamed-chunk-9-1.png)

``` r
grp6 <- cutree(hc, h=2.5)
table(grp6)
```

    ## grp6
    ##  1  2  3  4  5  6  7  8 
    ##  5 15  7  3  7 15  5  3

We can also use k= groups as an argument to cutree! Instead of h= height
========================================================================

### Using different Linkage methods in R

``` r
# Using different hierarchical clustering methods
hc.complete <- hclust(dist_matrix, method="complete")
plot(hc.complete)
```

![](Class08_files/figure-markdown_github/unnamed-chunk-10-1.png)

``` r
hc.average <- hclust(dist_matrix, method="average")
plot(hc.average)
```

![](Class08_files/figure-markdown_github/unnamed-chunk-10-2.png)

``` r
hc.single <- hclust(dist_matrix, method="single")
plot(hc.single)
```

![](Class08_files/figure-markdown_github/unnamed-chunk-10-3.png)

Made up overlapping data - a bit more real life like
====================================================

``` r
# Step 1. Generate some example data for clustering
x <- rbind(
 matrix(rnorm(100, mean=0, sd = 0.3), ncol = 2), # c1
 matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2), # c2
 matrix(c(rnorm(50, mean = 1, sd = 0.3), # c3
 rnorm(50, mean = 0, sd = 0.3)), ncol = 2))
colnames(x) <- c("x", "y")
# Step 2. Plot the data without clustering
plot(x)
```

![](Class08_files/figure-markdown_github/unnamed-chunk-11-1.png)

``` r
# Step 3. Generate colors for known clusters
# (just so we can compare to hclust results)
col <- as.factor( rep(c("c1","c2","c3"), each=50) )
plot(x, col=col)
```

![](Class08_files/figure-markdown_github/unnamed-chunk-11-2.png)

``` r
dist_realLife <- dist(x)

hc_realLife <- hclust(dist_realLife)
plot(hc_realLife)
```

![](Class08_files/figure-markdown_github/unnamed-chunk-12-1.png)

``` r
 gp2 <- cutree(hc_realLife, k=2)
 plot(x, col=gp2)
```

![](Class08_files/figure-markdown_github/unnamed-chunk-12-2.png)

``` r
 gp3 <- cutree(hc_realLife, k=3)
 plot(x, col=gp3)
```

![](Class08_files/figure-markdown_github/unnamed-chunk-12-3.png)

``` r
 gp4 <- cutree(hc_realLife, k=9)
 plot(x, col=gp4)
```

![](Class08_files/figure-markdown_github/unnamed-chunk-12-4.png)

### Principal Component Analysis (PCA)

We have so many different points it compares all the data in a way that is more easily interpreted
==================================================================================================

``` r
## You can also download this file from the class website!
mydata <- read.csv("https://tinyurl.com/expression-CSV",
 row.names=1)
# head command (first 6 lines of the data table)
head(mydata) 
```

    ##        wt1 wt2  wt3  wt4 wt5 ko1 ko2 ko3 ko4 ko5
    ## gene1  439 458  408  429 420  90  88  86  90  93
    ## gene2  219 200  204  210 187 427 423 434 433 426
    ## gene3 1006 989 1030 1017 973 252 237 238 226 210
    ## gene4  783 792  829  856 760 849 856 835 885 894
    ## gene5  181 249  204  244 225 277 305 272 270 279
    ## gene6  460 502  491  491 493 612 594 577 618 638

**NOTE**: prcomp() expects the samples to be rows and genes to be columns so we need to first transpose the matrix with the t() function!

``` r
## lets do PCA
# pca <- prcomp(t(mydata), scale=TRUE) 
# Take the transpose as the input data 

head( t(mydata))
```

    ##     gene1 gene2 gene3 gene4 gene5 gene6 gene7 gene8 gene9 gene10 gene11
    ## wt1   439   219  1006   783   181   460    27   175   658    121    337
    ## wt2   458   200   989   792   249   502    30   182   669    116    337
    ## wt3   408   204  1030   829   204   491    37   184   653    134    330
    ## wt4   429   210  1017   856   244   491    29   166   633    117    322
    ## wt5   420   187   973   760   225   493    34   180   657    133    313
    ## ko1    90   427   252   849   277   612   304   255   628    931    100
    ##     gene12 gene13 gene14 gene15 gene16 gene17 gene18 gene19 gene20 gene21
    ## wt1    214    789    458    551    390    900    951    436    244    119
    ## wt2    194    738    490    555    400    970    991    414    266     87
    ## wt3    213    807    493    527    403    905    991    388    228     87
    ## wt4    192    768    446    552    402    850    983    418    223     88
    ## wt5    207    820    496    503    401    834    984    410    240     93
    ## ko1     97    293    694    712    755    353    217    162    540    914
    ##     gene22 gene23 gene24 gene25 gene26 gene27 gene28 gene29 gene30 gene31
    ## wt1    156     89    570    788   1007    937    224    809    624    218
    ## wt2    170     97    567    796    972    876    232    869    598    259
    ## wt3    150     96    563    766    977    901    231    815    587    213
    ## wt4    167     97    587    778   1003    958    238    788    552    204
    ## wt5    155     82    563    825   1027    957    226    781    592    213
    ## ko1    346    788    424    456    945    414    850    482    956     69
    ##     gene32 gene33 gene34 gene35 gene36 gene37 gene38 gene39 gene40 gene41
    ## wt1    906    262    155    100    117    286    321    388    606    379
    ## wt2    798    291    172    104    147    262    353    372    576    377
    ## wt3    828    258    173     94    120    260    334    345    558    362
    ## wt4    874    271    173    114    147    270    340    373    581    346
    ## wt5    890    279    192     90    145    293    316    359    574    354
    ## ko1    541    534    643    212    353    360    642     50    415    991
    ##     gene42 gene43 gene44 gene45 gene46 gene47 gene48 gene49 gene50 gene51
    ## wt1    471    592    755     35    758     24    100    809    955    453
    ## wt2    492    615    733     40    734     25    113    825    994    419
    ## wt3    473    602    775     28    704     12    136    833    994    443
    ## wt4    470    602    687     25    761     13    117    800    975    459
    ## wt5    471    655    776     32    672     22    103    776    973    469
    ## ko1    401    514    255    947    567    324    912    538    175    174
    ##     gene52 gene53 gene54 gene55 gene56 gene57 gene58 gene59 gene60 gene61
    ## wt1    327    657    678    304    659    673    785    501    232    928
    ## wt2    320    669    638    325    687    668    772    513    228    936
    ## wt3    324    631    676    312    659    694    817    462    193   1015
    ## wt4    321    701    683    327    667    699    766    484    247    971
    ## wt5    318    647    671    320    639    726    784    504    231    964
    ## ko1    489    246    259    819    109     18    467     37    997    428
    ##     gene62 gene63 gene64 gene65 gene66 gene67 gene68 gene69 gene70 gene71
    ## wt1    159    336    968    339     35     27     80    744    766    672
    ## wt2    169    344    888    335     32     28     69    685    739    736
    ## wt3    163    372    907    373     45     25     87    733    751    672
    ## wt4    151    389    914    338     37     35     87    693    720    715
    ## wt5    166    357    883    328     38     27     81    746    738    693
    ## ko1    869    664    886    275    765    200    693    745    645    839
    ##     gene72 gene73 gene74 gene75 gene76 gene77 gene78 gene79 gene80 gene81
    ## wt1    526    627    468    986    348    719    883    837    666    804
    ## wt2    553    650    466    945    333    714    899    883    657    735
    ## wt3    534    664    477   1006    344    734    868    864    719    771
    ## wt4    511    622    469   1020    321    693    873    807    656    763
    ## wt5    529    606    494   1024    296    682    882    854    638    813
    ## ko1    922    805    703    359    770    620    803    210    549    613
    ##     gene82 gene83 gene84 gene85 gene86 gene87 gene88 gene89 gene90 gene91
    ## wt1    476    438    938     29    810    575    451    174    158    371
    ## wt2    494    430    934     29    830    579    471    170    122    367
    ## wt3    521    477    976     30    760    567    494    205    138    369
    ## wt4    494    457    965     19    796    565    447    175    159    339
    ## wt5    482    481    960     21    807    576    470    179    128    360
    ## ko1    183    466    904    618    486    352    540    298    863    103
    ##     gene92 gene93 gene94 gene95 gene96 gene97 gene98 gene99 gene100
    ## wt1    853    208    555    527    589    396     33    321      25
    ## wt2    798    214    584    573    607    384     27    343      34
    ## wt3    866    200    574    548    579    382     39    349      34
    ## wt4    843    196    599    548    536    399     42    367      36
    ## wt5    823    206    581    552    583    401     33    343      32
    ## ko1    934    409    292    686    497    460    977    949     661

``` r
## lets do PCA
pca <- prcomp(t(mydata), scale=TRUE)
summary(pca)
```

    ## Importance of components:
    ##                           PC1    PC2     PC3     PC4     PC5     PC6
    ## Standard deviation     9.6237 1.5198 1.05787 1.05203 0.88062 0.82545
    ## Proportion of Variance 0.9262 0.0231 0.01119 0.01107 0.00775 0.00681
    ## Cumulative Proportion  0.9262 0.9493 0.96045 0.97152 0.97928 0.98609
    ##                            PC7     PC8     PC9      PC10
    ## Standard deviation     0.80111 0.62065 0.60342 3.348e-15
    ## Proportion of Variance 0.00642 0.00385 0.00364 0.000e+00
    ## Cumulative Proportion  0.99251 0.99636 1.00000 1.000e+00

Make our first PCA plot
=======================

``` r
## A basic PC1 vs PC2 2-D plot
dim(pca$x)
```

    ## [1] 10 10

``` r
plot(pca$x[,1], pca$x[,2], xlab="PC1", ylab="PC2")
```

![](Class08_files/figure-markdown_github/unnamed-chunk-16-1.png)

``` r
## Variance captured per PC
pca.var <- pca$sdev^2 
# round command rounds up to a specified # of significant figures (i.e. 1 (shown))
round(pca.var/sum(pca.var)*100, 1) 
```

    ##  [1] 92.6  2.3  1.1  1.1  0.8  0.7  0.6  0.4  0.4  0.0

``` r
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
```

``` r
barplot(pca.var.per, main="Scree Plot",
 xlab="Principal Component", ylab="Percent Variation")
```

![](Class08_files/figure-markdown_github/unnamed-chunk-18-1.png)

Make our PCA plot nice
======================

``` r
## A vector of colors for wt and ko samples
colvec <- as.factor( substr( colnames(mydata), 1, 2) )
plot(pca$x[,1], pca$x[,2], col=colvec, pch=16,
 xlab=paste0("PC1 (", pca.var.per[1], "%)"),
 ylab=paste0("PC2 (", pca.var.per[2], "%)")) 
```

![](Class08_files/figure-markdown_github/unnamed-chunk-19-1.png)

``` r
# Takes a string of characters or numerals and strings it out (i.e. 1,2 arguments call for position 1 and 2 of each colname)
substr( colnames(mydata), 1, 2)
```

    ##  [1] "wt" "wt" "wt" "wt" "wt" "ko" "ko" "ko" "ko" "ko"

### Principal Component Analysis (PCA)

UK Foods
========

``` r
x <- read.csv("UK_foods.csv")
head(x) 
```

    ##                X England Wales Scotland N.Ireland
    ## 1         Cheese     105   103      103        66
    ## 2  Carcass_meat      245   227      242       267
    ## 3    Other_meat      685   803      750       586
    ## 4           Fish     147   160      122        93
    ## 5 Fats_and_oils      193   235      184       209
    ## 6         Sugars     156   175      147       139

``` r
dim(x)
```

    ## [1] 17  5

``` r
#17 rows, and 5 columns
```

``` r
# Note how the minus indexing works
# Note how the minus indexing works: the -1 includes everything but the first column (this is what the neg sign does)
rownames(x) <- x[,1]
x <- x[,-1]
head(x)
```

    ##                England Wales Scotland N.Ireland
    ## Cheese             105   103      103        66
    ## Carcass_meat       245   227      242       267
    ## Other_meat         685   803      750       586
    ## Fish               147   160      122        93
    ## Fats_and_oils      193   235      184       209
    ## Sugars             156   175      147       139

``` r
dim(x)
```

    ## [1] 17  4

**Side-note**: An alternative approach to setting the correct row-names in this case would be to read the data filie again and this time set the row.names argument of read.csv() to be the first column (i.e. use argument setting row.names=1), see below:

x &lt;- read.csv("data/UK\_foods.csv", row.names=1)
===================================================

Assign row names from the first col of the data upon reading. This is a safer approach.
=======================================================================================

``` r
x <- read.csv("UK_foods.csv", row.names=1)
head(x)
```

    ##                England Wales Scotland N.Ireland
    ## Cheese             105   103      103        66
    ## Carcass_meat       245   227      242       267
    ## Other_meat         685   803      750       586
    ## Fish               147   160      122        93
    ## Fats_and_oils      193   235      184       209
    ## Sugars             156   175      147       139

``` r
tail(x)
```

    ##                   England Wales Scotland N.Ireland
    ## Fresh_fruit          1102  1137      957       674
    ## Cereals              1472  1582     1462      1494
    ## Beverages              57    73       53        47
    ## Soft_drinks          1374  1256     1572      1506
    ## Alcoholic_drinks      375   475      458       135
    ## Confectionery          54    64       62        41

``` r
# Generating barplots
# *Note* beside defaults to FALSE (which stacks bars)
barplot(as.matrix(x), beside=T, col=rainbow(nrow(x)))
```

![](Class08_files/figure-markdown_github/unnamed-chunk-26-1.png)

``` r
pairs(x, col=rainbow(10), pch=16)
```

![](Class08_files/figure-markdown_github/unnamed-chunk-27-1.png)

``` r
# i.e. England is the y-axis for the first row, then each subsequent graph has the other countries on the x-axis
```

PCA to the rescue!
==================

``` r
# Use the prcomp() PCA function (note: t is the transpose)
pca <- prcomp( t(x) )
summary(pca)
```

    ## Importance of components:
    ##                             PC1      PC2      PC3       PC4
    ## Standard deviation     324.1502 212.7478 73.87622 4.189e-14
    ## Proportion of Variance   0.6744   0.2905  0.03503 0.000e+00
    ## Cumulative Proportion    0.6744   0.9650  1.00000 1.000e+00

``` r
# Plot PC1 vs PC2
plot(pca$x[,1], pca$x[,2], xlab="PC1", ylab="PC2", xlim=c(-270,500))
text(pca$x[,1], pca$x[,2], colnames(x), col=c("orange", "red", "blue", "darkgreen"))
```

![](Class08_files/figure-markdown_github/unnamed-chunk-29-1.png)

Loading values
==============

``` r
## Lets focus on PC1 as it accounts for > 90% of variance 
par(mar=c(10, 3, 0.35, 0))
barplot( pca$rotation[,1], las=2 )
```

![](Class08_files/figure-markdown_github/unnamed-chunk-30-1.png)
