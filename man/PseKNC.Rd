% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PseKNC.R
\name{PseKNC}
\alias{PseKNC}
\title{PseKNC}
\usage{
PseKNC(Seq)
}
\arguments{
\item{Seq}{the Sequence}
}
\value{
A vector contains the PseKNC, its formula is like \{0,0,1,1,0,1,0,0.5 \}.
   each 4 number denotes a nucleutide in the Seq.
}
\description{
The function returns a vector which return the PseKNC's nucleatiode part
PseKNC_each is the daughter-function for the funtion which is designed to lapply
}
