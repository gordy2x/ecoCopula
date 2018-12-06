

#' Plot graph of species intearactions.
#'
#' @param obj is a mvabund object, e.g. from output of saint
#' @return a plot of species interactions, positive/negative interactions are blue/pink,
#' @export
#' @examples
#' data(spider)
#' abund <- mvabund(spider$abund)
#' X <- spider$x
#' spider.mod=manyglm(abund~X)
#' spid.graph=saint(spider.mod)
#' plot(spid.graph)
#' 
#' 

plot.saint=function(obj,P=NULL) 
{
  labs=colnames(obj$best.graph$Y)
  
  Theta = -cov2cor(obj$best.graph$prec)
  Graph = (Theta != 0) * 1
  
  posneg = Theta
  diag(posneg) = 0
  posneg[posneg > 0] = "light blue"
  posneg[posneg < 0] = "pink"
  
  
  if(is.null(P)){
    P = sna::gplot.layout.fruchtermanreingold(Graph, list())
  }
  
  sna::gplot(Graph, gmode = "graph", label = labs, coord = P, 
             vertex.col = "blue", edge.col = posneg,
             label.cex = 0.8, edge.lwd = 4)

  invisible(P)
}

