# Project 1

``` r
knitr::opts_chunk$set(message = FALSE, warning = FALSE)
```

# Part 1 - Ready!

Next-Generation Sequencing (NGS) entails several high-throughput
techniques that generate big amount of data. The main two approaches
involve long-read sequencing with nanopore technology and short-read
sequencing with an Illumina machine. In the latter case, Polymerase
Chain Reaction (PCR) is performed by the sequence to amplify the library
and quantify gene expression.

In the ideal situation, we would expect to have a doubling of each RNA
strand in a linear fashion, so that each molecule is duplicated and the
ratio between the expression levels of the different genes is
maintained. However, this process introduces different types of errors
due to factors such as the composition of the sequence (for example high
fraction of G or C) or the secondary structure, stochasticity and
mutation that make it impossible to deduce the absolute gene expression,
especially for low-copy genes.

In addition, the sequencing process involves breaking the mRNA sequence
into several pieces because the machine cannot sequence a full-length
genome. Thus, multiple sequences may originate from the same sequence,
so we may think there are several copies of the same gene, but in
reality it is one broken copy. Therefore, we may infer a wrong
expression level. Instead, sequencing data provides a compositional
measure of gene expression, where every observation (or gene) is
quantified in relation with the rest of the pool. Thus, it is essential
to normalise reads from different samples when drawing between-sample
comparisons.

# Part 2 - Ready!

``` r
# Import package
library(pasilla)
```

``` r
# Load data
pasCts <- system.file("extdata", "pasilla_gene_counts.tsv",
                      package = "pasilla", mustWork = TRUE)

pasAnno <- system.file("extdata", "pasilla_sample_annotation.csv",
                       package="pasilla", mustWork = TRUE)

cts <- as.matrix(read.csv(pasCts, sep = "\t", row.names = "gene_id"))

coldata <- read.csv(pasAnno, row.names = 1)
coldata <- coldata[ , c("condition", "type")]

rownames(coldata) <- sub("fb", "", rownames(coldata))
cts <- cts[, rownames(coldata)]
```

The matrix `cts` contains the gene expression data in terms of the
number of reads for 14599 genes from 7 samples.

``` r
# View head of counts assay
head(cts, 10)
```

                treated1 treated2 treated3 untreated1 untreated2 untreated3
    FBgn0000003        0        0        1          0          0          0
    FBgn0000008      140       88       70         92        161         76
    FBgn0000014        4        0        0          5          1          0
    FBgn0000015        1        0        0          0          2          1
    FBgn0000017     6205     3072     3334       4664       8714       3564
    FBgn0000018      722      299      308        583        761        245
    FBgn0000022        0        0        0          0          1          0
    FBgn0000024       10        7        5         10         11          3
    FBgn0000028        0        1        1          0          1          0
    FBgn0000032     1698      696      757       1446       1713        615
                untreated4
    FBgn0000003          0
    FBgn0000008         70
    FBgn0000014          0
    FBgn0000015          2
    FBgn0000017       3150
    FBgn0000018        310
    FBgn0000022          0
    FBgn0000024          3
    FBgn0000028          0
    FBgn0000032        672

``` r
# Store assay dimensions
mat_dim <- dim(cts)
# View assay dimensions
mat_dim
```

    [1] 14599     7

Samples differ from one another in the library size. Therefore, read
counts should be normalised.

``` r
# View raw library sizes
colSums(cts)
```

      treated1   treated2   treated3 untreated1 untreated2 untreated3 untreated4 
      18670279    9571826   10343856   13972512   21911438    8358426    9841335 

Normalisation is performed by dividing the reads from every sample by
the corresponding library size and multiplying by the library size of
the first sample.

``` r
# Normalise all samples to library size of first sample
norm_cts <- apply(cts, 2, function(col) col / sum(col) * sum(cts[ , 1]))
# Check normalised library sizes
colSums(norm_cts)
```

      treated1   treated2   treated3 untreated1 untreated2 untreated3 untreated4 
      18670279   18670279   18670279   18670279   18670279   18670279   18670279 

# Part 3

``` r
find_log_stats <- function(assay) {
  
  # Compute row means and vars
  control_means <- rowMeans(assay)
  control_vars <- rowVars(assay)

  # Take log of means and vars with pseudocount
  control_log_means <- log(control_means + 1)
  control_log_vars <- log(control_vars + 1)

  # Combine data into a data.frame
  df <- data.frame(Var = control_log_vars, Mean = control_log_means)
  return(df)
}
```

``` r
# Subset data to untreated samples
control_cts <- norm_cts[ , grep("^untreated", colnames(norm_cts))]
# Find log means and vars for genes in control samples
control_stats <- find_log_stats(control_cts)
```

The Poisson distribution fits the data relatively well for lower numbers
of reads, but it creates bias for higher numbers of reads where the
variance becomes larger than the mean. Such phenomenon is called
*overdispersion*.

``` r
# Visualise fit of Poisson distribution to data
plot(control_stats$Mean, control_stats$Var, xlab = "Log Means", ylab = "Log Vars")
abline(a = 0, b = 1, col = "red")
legend("topleft", "Poisson Fit", lty = 1, col = "red")
```

![](project1.markdown_strict_files/figure-markdown_strict/unnamed-chunk-10-1.png)

``` r
# Fit dispersion coefficient of negative binomial distribution to data
nls_fit <- nls(Var ~ Mean + a * Mean^2, data = control_stats, start = list(a = 0))
summary(nls_fit)
```


    Formula: Var ~ Mean + a * Mean^2

    Parameters:
      Estimate Std. Error t value Pr(>|t|)    
    a 0.052201   0.000206   253.3   <2e-16 ***
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    Residual standard error: 0.8978 on 14598 degrees of freedom

    Number of iterations to convergence: 1 
    Achieved convergence tolerance: 1.408e-11

Unlike the Poisson distribution, the negative binomial distribution fits
well the data both for lower and higher numbers of reads.

``` r
# Define var-mean relationship of negative binomial distribution
disper_coef <- summary(nls_fit)$coefficients["a", "Estimate"]
neg_binom <- function(x) x + disper_coef * x^2

# Visualise fit of negative binomial distribution to data
plot(control_stats$Mean, control_stats$Var, xlab = "Log Means", ylab = "Log Vars")
curve(neg_binom, col = "red", add = TRUE)
legend("topleft", "Negative Binomial Fit", lty = 1, col = "red")
```

![](project1.markdown_strict_files/figure-markdown_strict/unnamed-chunk-12-1.png)

``` r
# Subset data to treated samples
treated_cts <- norm_cts[ , grep("^treated", colnames(norm_cts))]
# Find log means and vars for genes in treated samples
treated_stats <- find_log_stats(treated_cts)
```

``` r
# Visualise fit of Poisson distribution to data
plot(control_stats$Mean, treated_stats$Var, xlab = "Log Means", ylab = "Log Vars")
abline(a = 0, b = 1, col = "red")
legend("topleft", "Poisson Fit", lty = 1, col = "red")
```

![](project1.markdown_strict_files/figure-markdown_strict/unnamed-chunk-14-1.png)

``` r
# Fit dispersion coefficient of negative binomial distribution to data
nls_fit <- nls(Var ~ Mean + a * Mean^2, data = treated_stats, start = list(a = 0))
summary(nls_fit)
```


    Formula: Var ~ Mean + a * Mean^2

    Parameters:
       Estimate Std. Error t value Pr(>|t|)    
    a 0.0439520  0.0002798   157.1   <2e-16 ***
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    Residual standard error: 1.227 on 14598 degrees of freedom

    Number of iterations to convergence: 1 
    Achieved convergence tolerance: 7.173e-11

``` r
# Define var-mean relationship of negative binomial distribution
disper_coef <- summary(nls_fit)$coefficients["a", "Estimate"]
neg_binom <- function(x) x + disper_coef * x^2

# Visualise fit of negative binomial distribution to data
plot(treated_stats$Mean, treated_stats$Var, xlab = "Log Means", ylab = "Log Vars")
curve(neg_binom, col = "red", add = TRUE)
legend("topleft", "Negative Binomial Fit", lty = 1, col = "red")
```

![](project1.markdown_strict_files/figure-markdown_strict/unnamed-chunk-16-1.png)
