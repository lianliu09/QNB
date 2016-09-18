QNB

Version: 1.0

Date: 2016-6-29

Author: Lian Liu liulian19860905@163.com, Jia Meng jia.meng@hotmail.com

Maintainer: Lian Liu liulian19860905@163.com

The package is developed for differential methylation analysis. Estimate variance and mean dependence in count data from MeRIP-seq and test for differential
methylation based on a model using quadratic-negative-binomial distribution.


Depends: DESeq

Installation

QNB depends on DESeq and please make sure install them before installing QNB The package can be installed as:


source("http://bioconductor.org/biocLite.R")
biocLite("DESeq")

In addtion, R package devtools is required for QNB to be installed from GitHub.

install.packages("devtools")
At last, QNB can be installed as:

library("devtools")
install_github("lianliu09/QNB")
Toy Example

# using the data included in the package with replicates
library(QNB)
f1 <- system.file("extdata", "meth1.txt", package="QNB")
f2 <- system.file("extdata", "meth2.txt", package="QNB")
f3 <- system.file("extdata", "unmeth1.txt", package="QNB")
f4 <- system.file("extdata", "unmeth2.txt", package="QNB")

meth1 <- read.table(f1,header=TRUE)
meth2 <- read.table(f2,header=TRUE)
unmeth1 <- read.table(f3,header=TRUE)
unmeth2 <- read.table(f4,header=TRUE)

#When there are replicates under two conditions, we could select "mode=per-conditon" or 
#"mode=pooled" to estimate the dispersion. The default is "per-condition". 
result <- qnbtest(meth1, meth2,unmeth1,unmeth2)

#without replicates
f1 <- system.file("extdata", "no_rep_meth1.txt", package="QNB")
f2 <- system.file("extdata", "no_rep_meth2.txt", package="QNB")
f3 <- system.file("extdata", "no_rep_unmeth1.txt", package="QNB")
f4 <- system.file("extdata", "no_rep_unmeth2.txt", package="QNB")

no_rep_meth1 <- read.table(f1,header=TRUE)
no_rep_meth2 <- read.table(f2,header=TRUE)
no_rep_unmeth1 <- read.table(f3,header=TRUE)
no_rep_unmeth2 <- read.table(f4,header=TRUE)
head(no_rep_meth1)
head(no_rep_unmeth1)
result = qnbtest(no_rep_meth1, 
                 no_rep_meth2,
                 no_rep_unmeth1,
                 no_rep_unmeth2,
                 mode="blind")

#mode is auto
result = qnbtest(meth1, meth2,unmeth1,unmeth2,mode="auto")
