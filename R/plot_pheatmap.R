#' Plot Method for pheatmap Objects
#'
#' @description
#' This S3 method allows pheatmap objects (returned by \code{Celltype_Calculate()})
#' to be plotted using the generic \code{plot()} function. Without this method,
#' attempting to use \code{plot()} on a pheatmap object results in an error.
#'
#' @param x A pheatmap object, typically from \code{cluster_results$Heatmap_plot}
#' @param ... Additional arguments (currently ignored)
#'
#' @return Invisibly returns the input pheatmap object after displaying it
#'
#' @details
#' Pheatmap objects contain a gtable component that needs to be drawn using
#' grid graphics. This method handles that automatically when \code{plot()} is called.
#'
#' Alternative ways to display pheatmaps:
#' \itemize{
#'   \item \code{print(pheatmap_object)} - Works natively
#'   \item \code{plot(pheatmap_object)} - Works after loading SlimR
#'   \item \code{grid::grid.draw(pheatmap_object$gtable)} - Direct access
#' }
#'
#' @examples
#' \dontrun{
#' # After running Celltype_Calculate()
#' cluster_results <- Celltype_Calculate(
#'     seurat_obj = sce,
#'     gene_list = Markers_list,
#'     species = "Human"
#' )
#'
#' # Now both of these work:
#' print(cluster_results$Heatmap_plot)
#' plot(cluster_results$Heatmap_plot)
#' }
#'
#' @export
#' @method plot pheatmap
#' @importFrom grid grid.newpage grid.draw
plot.pheatmap <- function(x, ...) {
  # Validate input
  if (!inherits(x, "pheatmap")) {
    stop("Input must be a pheatmap object")
  }
  
  # Check if gtable component exists
  if (is.null(x$gtable)) {
    warning("pheatmap object does not contain a gtable component")
    return(invisible(x))
  }
  
  # Draw the heatmap using grid graphics
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  
  # Return invisibly for chaining
  invisible(x)
}
