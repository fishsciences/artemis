#' eDNA_data
#' 
#' Data from Delta Smelt eDNA experiments, collected at the Central Valley Project facility in California.
#' 
#' This sample dataset is designed to be representative of the type of data collected in a semi-controlled eDNA survey study.  The data consists of 624 technical replicates extracted from samples filtered at different volumes and taken at different distances from the 'source' of eDNA, which was a small underwater enclosure containing varying numbers of dead individual (hatchery-origin) Delta Smelt.  All samples were filtered from water flowing unidirectionally at the surface (depth < 1.0m).
#' 
#' The standard curve equation associated with this data is: Cq_star = -1.529*ln[eDNA concentration] + 21.168, where Cq_star is the Cq value estimated by the standard curve for the assay.
#' 
#' @format A data frame of 624 observations and 7 variables associated with individual qPCR replicates within filters.
#'  \describe{
#'    \item{Date}{Date of sample, yyyy-mm-dd}
#'    \item{FilterID}{Represents a single Sterivex filter.  Each unique FilterID is associated with 3 qPCR reactions (technical replicates)}
#'    \item{TechnicalRep}{Index of technical replicate from a given SampleID.}
#'    \item{Cq}{Quantification cycle of qPCR.  For this assay, the limit of detection was set at 40.0}
#'    \item{Distance_m}{Distance (in meters) from the 'source' of eDNA, which was a small underwater enclosure containing 22.5 deceased individual (hatchery-origin) Delta Smelt.}
#'    \item{Volume_mL}{Volume (in milliliters) of water pulled through the Sterivex filter in the sample.}
#'    \item{N_fish}{The number of individual Delta Smelt carcasses present in the submerged cage.}
#'    }
"eDNA_data"
