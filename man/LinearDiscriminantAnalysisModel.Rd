% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LinearDiscriminantAnalysisModel.r
\name{LinearDiscriminantAnalysisModel}
\alias{LinearDiscriminantAnalysisModel}
\title{Linear discriminant analysis model}
\arguments{
\item{SplitRatio}{Train and test split ratio. A value greater or equal than zero and less than one.}

\item{CV}{The number of folds for cross validation. It should be greater than one.}

\item{NCores}{The number of cores for parallel computing.}

\item{NewData}{A methylation beta values input from the ReadMethylFile function.}
}
\value{
A list
}
\description{
A function to train a linear discriminant analysis model to classify medulloblastoma subgroups using the DNA methylation beta values (Illumina Infinium HumanMethylation450). Prediction is followed by training if new data is provided.
}
\examples{
set.seed(123)
lda <- LinearDiscriminantAnalysisModel(SplitRatio = 0.8,
                                       CV = 2,
                                       NCores = 1,
                                       NewData = NULL)
}
