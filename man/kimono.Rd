% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kimono.R
\name{kimono}
\alias{kimono}
\title{Run Kimono - Knowledge-guIded Multi-Omic Netowrk inference}
\usage{
kimono(
  input_list,
  mapping_list,
  metainfo,
  main_layer = 1,
  min_features = 2,
  sel_iterations = 0,
  core = 1,
  ...
)
}
\arguments{
\item{input_list}{- list of omics data. First list element will be used as predictor}

\item{mapping_list}{- list of mappings between each data type one}

\item{main_layer}{- which input data set represents the main layer (default = 1)}

\item{min_features}{- autoexclude models with less than 2 features (default = 2)}

\item{core}{- if core != 1 kimono will perform an parallell computation}
}
\value{
a network in form of an edge table
}
\description{
Run Kimono - Knowledge-guIded Multi-Omic Netowrk inference
}
