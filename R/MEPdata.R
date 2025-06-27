#' Information of microexon-tag clusters in plant
#' 
#' MEPdata contains information of microexon-tag clusters in plants.
#' 
#' @usage data(MEPdata)
#' 
#' @format A list containing three elements:
#' \itemize{
#' \item MEPdata$cluster: a data.frame containing microexon-tag cluster 
#' (cluster), microexon size (zise) and phase (phase), the coding motif (motif),
#'  number of exons (exons) and the order of microexon (me_order).
#' \item MEPdata$blocks: an IntegerList containing exon blocks in each cluster.
#' \item MEPdata$matrix: a list of containng the DNA consensus 
#' matrix for ecach microexon cluster.
#' }
#' 
#' @examples
#' data(MEPdata)
#' names(MEPdata)
#' sapply(MEPdata, class)
#' head(MEPdata$cluster)
#' head(MEPdata$blocks)
#' names(MEPdata$matrix)
#' 
#' # Microexon cluster 8
#' MEPdata$cluster[8,]
#' MEPdata$blocks[[8]]
#' MEPdata$matrix[[8]]
#' ggseqlogo::ggseqlogo(MEPdata$matrix[[8]][,53:57])
#' 
"MEPdata"
