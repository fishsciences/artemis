#' eDNA_data
#' 
#' Data from Delta Smelt (\emph{Hypomesus transpacificus}) eDNA experiments, collected at the Central Valley Project facility in California.
#' 
#' This sample dataset is designed to be representative of the type of data collected in a semi-controlled eDNA survey study.  The data consists of 180 technical replicates processed from samples filtered at different distances from the 'source' of eDNA, which was a small underwater enclosure containing 100 live individual (cultured) Delta Smelt.  All samples were filtered from water flowing unidirectionally at the surface (depth < 1.0m). Filtered volume was fixed at either 50mL or 200mL. From distances of 10-50m, three  filters were taken in series every 10m at 50mL and 200mL, sampled from near to far relative to live car. Note that the live car itself (Distance_m = ~0) was not actually sampled. Full details on this data are presented in Espe et al. 2021 (in prep).
#' 
#' The standard curve equation associated with this data is: Cq_star = -1.529*ln[eDNA concentration] + 21.168, where Cq_star is the Cq value estimated by the standard curve for the assay.
#' 
#' @format A data frame of 180 observations and 9 variables associated with individual qPCR replicates within filters.
#'  \describe{
#'    \item{Date}{Date of sample, yyyy-mm-dd}
#'    \item{FilterID}{Represents a single Sterivex filter.  Each unique FilterID is associated with 3 qPCR reactions (technical replicates)}
#'    \item{TechRep}{Index of technical replicate from a given SampleID.}
#'    \item{Cq}{Quantification cycle of qPCR.  For this assay, the limit of detection was set at 40.0}
#'    \item{Distance_m}{Distance (in meters) from the 'source' of eDNA, which was a small underwater enclosure containing 100 live individual (cultured) Delta Smelt.}
#'    \item{Volume_mL}{Volume (in milliliters) of water pulled through the Sterivex filter in the sample.}
#'    \item{Biomass_N}{The number of individual Delta Smelt carcasses present in the submerged cage.}
#'    \item{StdCrvAlpha_lnForm}{The intercept of the standard curve equation associated with these filters, in ln form}
#'    \item{StdCrvBeta_lnForm}{The slope of the standard curve equation associated with these filters, in ln form}
#'    }
"eDNA_data"
