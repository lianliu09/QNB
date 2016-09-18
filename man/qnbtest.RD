\name{qnbtest}
\alias{qnbtest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
This is the main function of differential methylation analysis for MeRIP-seq.
}
\description{
This function tests for differential methylation between two conditions.
}
\usage{
qnbtest(meth1, meth2, unmeth1, unmeth2, mode = "per-condition", plot.dispersion = TRUE, pvals.only = TRUE,output.dir=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{meth1}{
        a  reads count matrix for IP samples of one condition
}
  \item{meth2}{
        a  reads count matrix for IP samples of two condition
}
  \item{unmeth1}{
        a  reads count matrix for input samples of one condition
}
  \item{unmeth2}{
        a  reads count matrix for input samples of two condition
}
  \item{mode}{
There are four ways how the empirical dispersion can be computed:

\itemize{
  \item{pooled}{
Use the samples from all conditions with replicates to estimate a single pooled empirical dispersion value, called "pooled", and assign it to all samples.
}
  \item{per-condition}{
For each condition with replicates, compute an empirical dispersion value by considering the data from samples for this condition. The default is per-conditon.
}
  \item{blind}{
Ignore the sample labels and compute an empirical dispersion value as if all samples were replicates of a single condition. This can be done even if there are no biological 
replicates.}
  \item{auto}{
select mode according to the size of samples automaticly.
}
}

}
  \item{plot.dispersion}{
the default is TRUE. If plot.dispersion = FALSE, it will not save the dispersion figure. 
}
  \item{pvals.only}{
get pvalue only. If pvals.only = TRUE, get pvalue only in result matrix.If 
pvals.only = FALSE,get pvalue and FDR in result matrix.
}
  \item{output.dir}{
The saved file path. The default is NA. If output.dir = NA, the path is the current path.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Results will get a matrix including 8 columns(pvalue,mean,variance,fold-change,IP 
fold-change,input fold-change,expression,FDR(if pvals.only=FALSE))


  \item{pvalue}{
Indicate the significance of the methylation site as an RNA differential methylation site
}
  \item{log2.fc}{
  Fold-change. log 2 (fold enrichment)  within the peak in the IP sample compared with the input sample.
  }
  \item{q}{
  The expression of each methylation site.
  }
  \item{FDR}{
  fdr of the methylation site, indicating the significance of the peak as an RNA 
  differential methylation site after multiple hypothesis correction
  }
%%  ~Describe the value returned
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Lian Liu <liulian19860905@163.com>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## The function is currently defined as
# load library and specify the parameters
library("QNB")

path=system.file("extdata", "MeRIP_Data.RData", package="QNB")
load(path)
#with replicates, and the number of replicate for meth and unmeth must be the same, but meth1 and #meth2 may be different
res=QNB(meth1,meth2,unmeth1,unmeth2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ differential methylation }
\keyword{ quadratic-negative-binomial model }% __ONLY ONE__ keyword per line
