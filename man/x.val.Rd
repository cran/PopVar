% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cross_validation.R
\name{x.val}
\alias{x.val}
\title{Estimate genome-wide prediction accuracy using cross-validation}
\usage{
x.val(
  G.in = NULL,
  y.in = NULL,
  min.maf = 0.01,
  mkr.cutoff = 0.5,
  entry.cutoff = 0.5,
  remove.dups = TRUE,
  impute = "EM",
  frac.train = 0.6,
  nCV.iter = 100,
  nFold = NULL,
  nFold.reps = 1,
  return.estimates = FALSE,
  CV.burnIn = 750,
  CV.nIter = 1500,
  models = c("rrBLUP", "BayesA", "BayesB", "BayesC", "BL", "BRR"),
  saveAt = tempdir()
)
}
\arguments{
\item{G.in}{\code{Matrix} of genotypic data. First row contains marker names and the first column contains entry (taxa) names. Genotypes should be coded as follows: \itemize{
\item \code{1}: homozygous for minor allele
\item \code{0}: heterozygous
\item \code{-1}: homozygous for major allele
\item \code{NA}: missing data
\item Imputed genotypes can be passed, see \code{impute} below for details
}
TIP - Set header=\code{FALSE} within \code{\link{read.table}} or \code{\link{read.csv}} when importing a tab-delimited file containing data for \code{G.in}.}

\item{y.in}{\code{Matrix} of phenotypic data. First column contains entry (taxa) names found in \code{G.in}, regardless of whether the entry has a phenotype for any or all traits. Additional columns contain phenotypic data; column names should reflect the trait name(s). TIP - Set header=\code{TRUE} within \code{read.table} or \code{read.csv} when importing a tab-delimited file containing dat}

\item{min.maf}{Optional \code{numeric} indicating a minimum minor allele frequency (MAF) when filtering \code{G.in}. Markers with an MAF < \code{min.maf} will be removed. Default is \code{0.01} to remove monomorphic markers. Set to \code{0} for no filtering.}

\item{mkr.cutoff}{Optional \code{numeric} indicating the maximum missing data per marker when filtering \code{G.in}. Markers missing > \code{mkr.cutoff} data will be removed. Default is \code{0.50}. Set to \code{1} for no filtering.}

\item{entry.cutoff}{Optional \code{numeric} indicating the maximum missing genotypic data per entry allowed when filtering \code{G.in}. Entries missing > \code{entry.cutoff} marker data will be removed. Default is \code{0.50}. Set to \code{1} for no filtering.}

\item{remove.dups}{Optional \code{logical}. If \code{TRUE} duplicate entries in the genotype matrix, if present, will be removed. This step may be necessary for missing marker imputation (see \code{impute}). Default is \code{TRUE}.}

\item{impute}{Options include \code{c("EM", "mean", "pass")}. By default (i.e. \code{"EM"}), after filtering missing genotypic data will be imputed via the EM algorithm implemented in \code{\link[rrBLUP]{A.mat}}  (\cite{Endelman, 2011}; \cite{Poland et al., 2012}). If \code{"mean"} missing genotypic data will be imputed via the 'marker mean' method, also implemented in \code{\link[rrBLUP]{A.mat}}. Enter \code{"pass"} if a pre-filtered and imputed genotype matrix is provided to \code{G.in}.}

\item{frac.train}{Optional \code{numeric} indicating the fraction of the TP that is used to estimate marker effects (i.e. the prediction set) under cross-validation (CV) method 1 (see \code{Details}). The remaining \eqn{(1-frac.trait)} of the TP will then comprise the prediction set.}

\item{nCV.iter}{Optional \code{integer} indicating the number of times to iterate \emph{CV method 1} described in \code{Details}. Default is \code{100}.}

\item{nFold}{Optional \code{integer}. If a number is provided, denoting the number of "folds", then CV will be conducted using \emph{CV method 2} (see \code{Details}). Default is \code{NULL}, resulting in the default use of the \emph{CV method 1}.}

\item{nFold.reps}{Optional \code{integer} indicating the number of times \emph{CV method 2} is repeated. The CV accuracy returned is the average \emph{r} of each rep. Default is \code{1}.}

\item{return.estimates}{Optional \code{logical}. If \code{TRUE} additional items including the marker effect and beta estimates from the selected prediction model (i.e. highest CV accuracy) will be returned.}

\item{CV.burnIn}{Optional \code{integer} argument used by \code{\link[BGLR]{BGLR}} when fitting Bayesian models. Default is \code{750}.}

\item{CV.nIter}{Optional \code{integer} argument used by \code{\link[BGLR]{BGLR}} (\cite{de los Compos and Rodriguez, 2014}) when fitting Bayesian models. Default is \code{1500}.}

\item{models}{Optional \code{character vector} of the regression models to be used in CV and to estimate marker effects. Options include \code{rrBLUP, BayesA, BayesB, BayesC, BL, BRR}, one or more may be included at a time. By default all models are tested.}

\item{saveAt}{When using models other than "rrBLUP" (i.e. Bayesian models), this is a path and prefix for saving temporary files
the are produced by the \code{\link[BGLR]{BGLR}} function.}
}
\value{
A list containing: \itemize{
\item \code{CVs} A \code{dataframe} of CV results for each trait/model combination specified
\item If \code{return.estimates} is \code{TRUE} the additional items will be returned: \itemize{
\item \code{models.used} A \code{list} of the models chosen to estimate marker effects for each trait
\item \code{mkr.effects} A \code{vector} of marker effect estimates for each trait generated by the respective prediction model used
\item \code{betas} A \code{list} of beta values for each trait generated by the respective prediction model used
}
}
}
\description{
\code{x.val} performs cross-validation (CV) to estimate the accuracy of genome-wide prediction (otherwise known as genomic selection) for a specific training population (TP), i.e. a set of individuals for which phenotypic and genotypic data is available. Cross-validation can be conducted via one of two methods within \code{x.val}, see \code{Details} for more information.

\if{html}{\out{<div class="sourceCode">}}\preformatted{         NOTE - \code{x.val}, specifically \code{\link[BGLR]{BGLR}} writes and reads files to disk so it is highly recommended to set your working directory
}\if{html}{\out{</div>}}
}
\details{
Two CV methods are available within \code{PopVar}: \itemize{
\item \code{CV method 1}: During each iteration a training (i.e. model training) set will be \strong{randomly sampled} from the TP of size \eqn{N*(frac.train)}, where \emph{N} is the size of the TP, and the remainder of the TP is assigned to the validation set. The accuracies of individual models are expressed as average Pearson's correlation coefficient (\emph{r}) between the genome estimated breeding value (GEBV) and observed phenotypic values in the validation set across all \code{nCV.iter} iterations. Due to its amendibility to various TP sizes, \emph{CV method 1} is the default CV method in \code{\link{pop.predict}}.
\item \code{CV method 2}: \code{nFold} \strong{independent} validation sets are sampled from the TP and predicted by the remainder. For example, if \eqn{nFold = 10} the TP will be split into 10 equal sets, each containing \eqn{1/10}-th of the TP, which will be predicted by the remaining \eqn{9/10}-ths of the TP. The accuracies of individual models are expressed as the average (\emph{r}) between the GEBV and observed phenotypic values in the validation set across all \code{nFold} folds. The process can be repeated \code{nFold.reps} times with \code{nFold} new independent sets being sampled each replication, in which case the reported prediction accuracies are averages across all folds and replications.
}
}
\examples{
\donttest{

## The following examples only use the model 'rrBLUP' for the sake of testing. Functions
## BGLR package write temporary files to the disk.

## CV using method 1 with 25 iterations
CV.mthd1 <- x.val(G.in = G.in_ex, y.in = y.in_ex, nCV.iter = 25, models = "rrBLUP")
CV.mthd1$CVs

## CV using method 2 with 5 folds and 3 replications
x.val(G.in = G.in_ex, y.in = y.in_ex, nFold = 5, nFold.reps = 3, models = "rrBLUP")
}
}
