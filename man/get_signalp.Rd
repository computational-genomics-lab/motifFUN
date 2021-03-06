% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signalp.R
\name{get_signalp}
\alias{get_signalp}
\title{Predicts the secretory proteins}
\usage{
get_signalp(signalp.path = NULL, signalp.version = NULL,
  input_file = NULL, org.type = NULL)
}
\arguments{
\item{signalp.path}{Localpath of folder containing signalp binary excutable file or the the the excutable file itself, if path not specified, then signalp must be in program search path}

\item{signalp.version}{please provide signalp version to use}

\item{input_file}{Please provide fasta input file}

\item{org.type}{Provide organism type [euk,gram+,gram-]}
}
\value{
Provide  dataframe its contain summary of the signalp results file

 @examples

 signalp <- get_signalp(signalp.path = NULL, signalp.version = NULL, input_file = NULL, org.type = NULL)
}
\description{
This function predicts the presence and location of signal peptide cleavage sites in amino acid sequences from different organisms(eukaryotes,gram+,gram-).
}
