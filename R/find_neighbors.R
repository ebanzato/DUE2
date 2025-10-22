#' Finding neighbors
#'
#' Finding neighbors of nodes in a graph.
#'
#' @param g A symmetric adjacency matrix representing the undirected graph.
#'
#' @return a list. Each list element represents a graph node and contains its neighboring nodes.
#'
#' @examples
#'
#' set.seed(1)
#' g <- as.matrix(igraph::as_adjacency_matrix(igraph::sample_gnp(6, 0.7)))
#' colnames(g) <- rownames(g) <- paste0('X',1:6)
#'
#' ng <- find_neighbors(g)
#'

find_neighbors = function(graph){

  if(!isSymmetric(graph)){
    stop('The adjacency matrix is not symmetric.')
  }
  g_nodes = colnames(graph)

  out = lapply(1:length(g_nodes), function(g){

    gcol = graph[, which(colnames(graph) == g_nodes[g])]
    gnei = g_nodes[which(gcol == 1)]
    gnei

  })

  names(out) = g_nodes

  return(out)
}
