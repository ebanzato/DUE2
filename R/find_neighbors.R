#' Finding neighbors
#'
#' Finding neighbors of nodes in a graph.
#'
#' @param g An undirected graph represented as an igraph object.
#'
#' @return a list. Each list element represents a graph node and contains its neighboring nodes.
#'
#' @examples
#'
#' g <- igraph::make_graph(c(1,2,2,3,3,1,3,4,4,5,5,6,6,3,6,4), directed=F)
#'
#' ng <- find_neighbors(g)
#'

find_neighbors = function(g){

  if(!igraph::is_igraph(g)){
    stop('Not an igraph object.')
  }

  if(igraph::is_directed(g)){
    stop('Not an undirected graph.')
  }

  g_nodes = igraph::V(g)
  out = lapply(g_nodes, function(node) igraph::neighbors(g, node))
  return(out)
}
