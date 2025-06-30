#' Microexon-tag clusters information in plants
#' 
#' A dataset containing comprehensive information about microexon-tag clusters 
#' in plant genomes, including cluster features, exon blocks, and DNA consensus matrices.
#' 
#' @usage data(MEPdata)
#' 
#' @format A list with three components:
#' \describe{
#'   \item{cluster}{A data frame with 11 variables describing microexon clusters:
#'     \describe{
#'       \item{cluster}{Microexon-tag cluster identifier}
#'       \item{size}{Microexon size in nucleotides}
#'       \item{phase}{Microexon phase (0, 1, or 2)}
#'       \item{motif}{Coding motif pattern}
#'       \item{exons}{Number of exons in cluster}
#'       \item{me_order}{Microexon order in the tag}
#'       \item{rep_complex}{Representative sequence complexity score}
#'       \item{string}{Conserved sequence string}
#'       \item{L_complex}{Left border sequence complexity score}
#'       \item{R_complex}{Right border sequence complexity score}
#'       \item{low_copy}{Logical indicating low-copy cluster (TRUE/FALSE)}
#'     }
#'   }
#'   \item{blocks}{An IntegerList of exon blocks for each cluster}
#'   \item{matrix}{A list of DNA consensus matrices for each cluster}
#' }
#' 
#' @examples
#' # Load data
#' data(MEPdata)
#' 
#' # Explore data structure
#' names(MEPdata)
#' sapply(MEPdata, class)
#' 
#' # Examine cluster information
#' head(MEPdata$cluster)
#' 
#' # View exon blocks for first 5 clusters
#' head(MEPdata$blocks, 5)
#' 
#' # View matrix names
#' head(names(MEPdata$matrix))
#' 
#' # Analyze specific cluster (e.g., cluster 8)
#' MEPdata$cluster[8, ]
#' MEPdata$blocks[[8]]
#' MEPdata$matrix[[8]]
#' 
#' # Visualize sequence logo (requires ggseqlogo)
#' if(requireNamespace("ggseqlogo", quietly = TRUE)) {
#'   ggseqlogo::ggseqlogo(MEPdata$matrix[[8]][, 53:57])
#' }
#' 
"MEPdata"
