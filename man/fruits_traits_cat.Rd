% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fruits_traits_cat.R
\docType{data}
\name{fruits_traits_cat}
\alias{fruits_traits_cat}
\title{Dataset: Fruits Traits Informations}
\format{
A data frame with 8 rows (traits) and the following three columns:
\describe{
\item{trait_name}{a character giving the trait name}
\item{trait_type}{a character giving the trait type (\strong{O} for Ordinal
trait, \strong{N} for Nominal trait, \strong{Q} for Quantitative trait, and \strong{F}
for Fuzzy-coded trait)}
\item{fuzzy_name}{a character giving the name of fuzzy-coded traits
(i.e. \code{Use}) to  which 'sub-traits' (i.e. \code{raw}, \code{pastry}, and \code{jam})
belongs}
}
}
\usage{
fruits_traits_cat
}
\description{
This dataset summarizes information about the 6 traits used in the
\code{fruits_traits} dataset.
}
\note{
If your dataset does not contain fuzzy trait, the column \code{fuzzy_name}
can be ignored but the first two columns are mandatory.

Traits in this dataset correspond to columns (traits) of the \code{fruits_traits}
dataset.
}
\examples{
\dontrun{
# Load Traits Information
data("fruits_traits_cat", package = "mFD")
fruits_traits_cat

# Load Species x Traits Data Frame
data("fruits_traits", package = "mFD")

# Summarize Species x Traits Data
mFD::sp.tr.summary(tr_cat = fruits_traits_cat, sp_tr = fruits_traits)
}
}
\seealso{
\code{fruits_traits}, \code{baskets_fruits_weights}
}
\keyword{datasets}
