#' Finding neighbors
#'
#' Finding neighbors of nodes in a graph.
#'
#' @param graph An igraph object.
#'
#' @return a list. Each element of the list corresponds to a node in the graph, with
#' the corresponding neighboring nodes.
#'
#' @examples
#'
#' g <- igraph::make_graph(c(1,2,2,3,3,1,3,4,4,5,5,6,6,3,6,4),directed=F)
#'
#' ng <- find_neighbors(g)
#'

find_neighbors = function(graph){

  if(!igraph::is.igraph(graph)){
    stop('Not an igraph object.')
  }

  g_nodes = igraph::V(graph)
  out = lapply(g_nodes, function(node) igraph::neighbors(graph, node))
  return(out)
}
