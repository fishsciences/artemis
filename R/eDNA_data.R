#' eDNA_data
#' 
#' A subset of the data collected in Delta Smelt eDNA live car experiments in the Sacramento River Delta.
#' 
#' This sample dataset is designed to be representative of the type of data collected in a semi-controlled eDNA survey study.  The data consists of 732 technical replicates extracted from samples filtered at different volumes and taken at different distances from the 'source' of eDNA, which was a small underwater enclosure ("live car") containing 22.5 deceased individual (hatchery-origin) Delta Smelt.  All samples were filtered from water at the surface (depth < 1.0m).
#' 
#' @format A data frame of 732 observations and 7 variables associated with individaul eDNA samples.
#'  \describe{
#'    \item{Date}{Date of sample, yyyy-mm-dd}
#'    \item{SampleID}{Represents a single Sterivex filter, associated with 8-12 extractions (technical replicates)}
#'    \item{TechnicalRep}{Index of technical replicate from a given SampleID.}
#'    \item{FilterNumber}{Index (from 1-15) of filters associated with a particular volume/distance combination.  For example, on 2019-01-24, fifteen replicate filters were taken at 50mL, 0m from source, while only 5 replicate filters are associated with this volume/distance combo on 2019-03-21.}
#'    \item{Distance}{Distance (in meters) from the 'source' of eDNA, which was a small underwater enclosure containing 22.5 deceased individual (hatchery-origin) Delta Smelt.}
#'    \item{Volume}{Volume (in milliliters) of water pulled through the Sterivex filter in the sample.}
#'    \item{Cq}{Quantification cycle of qPCR - lower values represent higher concentration of eDNA.  Cq == 40.0 indicates no eDNA present (the cutoff value for amplification cycles during qPCR).}
#'    }
"eDNA_data"