
#' DUE2 function
#'
#' Function to identify differences in networks between two conditions.
#'
#' 1. Screening stage
#'
#' 2. Confirmatory stage
#'
#' @param graph A symmetric adjacency matrix representing the undirected graph. Node names must match the variable names in the data.
#'
#' @param data A data frame with the observed variables. Variable names in the data must match the node names.
#'
#' @param group A vector of group labels for comparison, one entry per observation.
#'
#' @param glm.family Description of the distribution and link function to be used in the model.
#'
#' @param alpha The alpha level for conducting the test. Default is 0.05.
#'
#' @param method.FDR Method for adjusting p-values in the first stage to control the FDR. Default is 'BH'.
#'
#' @param method.FWER Method for adjusting p-values in the first stage to control the FWER. Default is 'holm'.
#'
#' @param cell.group A vector of cell group labels, used as an extra covariate in the model and not for comparison, with length equal to the sample size.
#'
#' @param size.factor A vector of size factors with length equal to the sample size. It is used as an offset and to scale the covariates in the model.
#'
#' @param progressbar Logical, whether to display the progress bar. Default is TRUE.
#'
#' @return a list with four elements.
#' 1. diff: A matrix of dimension p x (p + 2), where p is the number of nodes, indicating differences in the network.
#' The columns are as follows:
#' 1st column: First screening stage: 0 for equality, 1 for nodes with a difference.
#' 2nd column: Difference in the mean parameter of the node in the row: 0 for no difference, 1 for a difference, NA if the first column is 0.
#' 3-(p+2)th columns: Differences in the association parameters between the node in the row and nodes in the columns: 0 for no difference, 1 for a difference, NA if not tested.
#'
#' 2. delta_pvalues: corrected p-values of the two stages.
#'
#' 3. delta_coef: coefficients of the difference parameters.
#'
#' 4. base_coef: coefficients of the reference condition.
#'
#' @examples
#' #'
#' set.seed(1)
#' graph <- as.matrix(igraph::as_adjacency_matrix(igraph::sample_gnp(6, 0.7)))
#' colnames(graph) <- rownames(graph) <- paste0('X',1:6)
#'
#' data1 <- data.frame(matrix(rpois(600, 3), ncol=6)); colnames(data1) = paste0('X',1:6)
#' data2 <- data.frame(matrix(rpois(600, 5), ncol=6)); colnames(data2) = paste0('X',1:6)
#' data  <- rbind(data1, data2)
#'
#' group <- rep(c(0,1), each=100)
#'
#' DUE2(graph, data, group, glm.family='poisson', alpha=0.05, method.FDR='BH', method.FWER='holm')
#'

DUE2 = function(graph, data, group, glm.family, alpha=0.05, method.FDR='BH', method.FWER='holm', cell.group=NULL, size.factor=NULL, progressbar=TRUE){

  # Check if adjm names match
  if(!identical(colnames(graph), rownames(graph))){
    stop('No match between columns names and rows names of \'graph\'.')
  }

  # Check if 'graph' is symmetric
  if(!isSymmetric(graph)){
    warning('The adjacency matrix is not symmetric.')
  }

  # Check if data contains all the variables in graph
  if(!sum(colnames(graph) %in% colnames(data)) == length(colnames(graph))){
    stop('Not all variables in \'graph\' are in \'data\' OR no match between nodes and variables names')
  }

  # Group and cell group as factors
  if(!is.null(group)){
    group = as.factor(group)
  }
  if(!is.null(cell.group)){
    cell.group = as.factor(cell.group)
  }

  # Tests
  test.res = test_known(graph, data, group, glm.family, cell.group, size.factor, progressbar)

  # Two stages alg
  res = two_stages_pval(graph, test.pval = test.res$p.mat, method.FDR, method.FWER, alpha)

  # OUT

  out = list('g.diff' = res, 'delta_pvalues' = test.res$p.mat,
             'delta_coef' = test.res$c.mat, 'base_coef'= test.res$b.mat)
  return(out)
}

