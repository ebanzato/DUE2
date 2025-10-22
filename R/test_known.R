#' Regressions
#'
#' Function to estimate the models to compare.
#'
#' @param graph An igraph object.
#'
#' @param data A dataset with the columns name corresponding to the names of the nodes.
#'
#' @param group a vector indicating the groups to compare
#'
#' @param glm.family a description of the error distribution and link function to be used in the model.
#'
#' @param cell.group vector specifying the groups of observations.
#'
#'
#' @return a list with two objects: a matrix of pvalues and a matrix of coefficients.
#'
#'
#' @examples
#'
#' set.seed(1)
#' g <- igraph::sample_gnp(6, 0.7); igraph::V(g) = paste0('X',1:6)
#'
#' g.df <- data.frame(matrix(rpois(600, 3), ncol=6)); colnames(g.df) = 1:6
#'
#' g.gr <- rep(c(0,1), each=50)
#'
#' test.g <- test_known(g, g.df, g.gr, glm.family='Poisson', cell.group=NULL)
#'

###### SISTEMARE ############
#
# 1. esempio da sistemare
#
# 2. aggiungere la possibilità di ottenere la matrice dei coefficienti

test_known = function(graph, data, group, glm.family, cell.group){

  if(glm.family=='quasipoisson'){
    test.aov = 'F'
  }else{
    test.aov = 'Chisq'
  }

  if(sum(names(igraph::V(graph)) %in% colnames(data)) >= length(igraph::V(graph))){

    neigh = find_neighbors(graph)
    n_nodes = names(igraph::V(graph))

    results = lapply(1:ncol(data), function(j) {

      n_neigh = neigh[[j]]

      # out
      res = rep(NA, length(n_nodes)+2)
      names(res) = c('LRT','group',n_nodes)

      if(is.null(cell.group)){
        # define formula
        formula_h0 = paste0(names(V(graph)[j]), ' ~ ', paste(names(neigh[[j]]), collapse='+', sep=''))
        formula_h1 = paste0(names(V(graph)[j]), ' ~ (', paste(names(neigh[[j]]), collapse='+', sep=''),')*group')
        # models
        mod_h0 = glm(formula_h0, data.frame(data,'group'=group), family=glm.family)  # model H0
        mod_h1 = glm(formula_h1, data.frame(data,'group'=group), family=glm.family)  # model H1
        cg = 0
      }else{
        # define formula
        formula_h0 = paste0(names(V(graph)[j]), ' ~ cgroup +', paste(names(neigh[[j]]), collapse='+', sep=''))
        formula_h1 = paste0(names(V(graph)[j]), ' ~ cgroup + (', paste(names(neigh[[j]]), collapse='+', sep=''),')*group')
        # models
        mod_h0 = glm(formula_h0, data.frame(data,'group'=group,'cgroup'=cell.group), family=glm.family)  # model H0
        mod_h1 = glm(formula_h1, data.frame(data,'group'=group,'cgroup'=cell.group), family=glm.family)  # model H1
        cg = length(levels(as.factor(cell.group)))-1
      }

      # Test
      test = anova(mod_h0, mod_h1, test=test.aov)
      if(test.aov=='Chisq'){
        res[1] = test$`Pr(>Chi)`[2]
      }
      if(test.aov=='F'){
        res[1] = test$`Pr(>F)`[2]
      }
      pvalbeta = summary(mod_h1)$coefficients[-c(1:(1+cg+length(n_neigh))),4]
      # out
      res[2] = pvalbeta[1]
      res[-c(1:2)][which(n_nodes %in% names(neigh[[j]]))] = pvalbeta[-1]
      return(res)
    })

    out = t(sapply(results, c))
    rownames(out) = n_nodes

    # out
    return(out)
  }else{
    print('No match between nodes and variables names')
  }
}
