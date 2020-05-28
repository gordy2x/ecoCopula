

#' Plot graph of direct species associations.
#'
#' @param obj is a cgr object, e.g. from output of cgr.
#' @param P locations of graph nodes, if NULL (default) these are generated with a Fruchterman Reingold algorithm.
#' @param ...	other parameters to be passed through to plotting gplot, in particular \code{pad}, the amount to pad the plotting range is useful if labels are being clipped.  
#' @return a plot of species associations after accounting for the effect of all other species, positive/negative are blue/pink.
#' The matrix of node positions (\code{P}) is returned silently.
#' @export
#' @examples
#' data(spider)
#' abund <- mvabund(spider$abund)
#' X <- spider$x
#' spider_mod=manyglm(abund~1)
#' spid_graph=cgr(spider_mod)
#' plot(spid_graph)

plot.cgr = function(obj, P = NULL,plot.raw=FALSE, ...) {
    

    #Partial correlations
    Theta = -cov2cor(obj$best_graph$prec)
    
    #graph of partial correlations
    Graph = (Theta != 0) * 1
    
    #color of lines from sign of partial correlations
    posneg = Theta
    diag(posneg) = 0
    posneg[posneg > 0] = "light blue"
    posneg[posneg < 0] = "pink"
    
    #use algorithm to arrange graph nodes
    if (is.null(P)) {
        P = sna::gplot.layout.fruchtermanreingold(Graph, list())
    }
    
    #plot graph
    sna::gplot(Graph, gmode = "graph", label = colnames(obj$best_graph$Y), coord = P, 
               vertex.col = "blue", edge.col = posneg, 
               label.cex = 0.8, edge.lwd = 4, ...)
    
    #raturn location of nodes
    invisible(P)
}

