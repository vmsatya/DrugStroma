nRows <- 10
adMatrix <- matrix (round (runif (nRows * nRows)), ncol=nRows)
adMatrix <- adMatrix * upper.tri (adMatrix, diag=T)
colnames (adMatrix) <- 1:nRows
adMatrix <- tmp2
g <- new ("graphAM", adjMat = adMatrix, edgemode="directed")

nodeList <- nodes (g)
nodeListLength <- length (nodeList)
from <- c ()
to <- c ()
nodeInteger <- c ()
nodeChar <- c ()
nodeFloat <- c()
edgeInteger <- c ()
edgeChar <- c ()
edgeFloat <- c ()

for (i in seq (1, nodeListLength)) {
  node = nodeList [i]
  nodeInteger <- c (nodeInteger, round (runif (1, max = 500)))
  chars = paste (sample (letters), collapse="")
  nodeChar <- c (nodeChar, chars)
  nodeFloat <- c (nodeFloat, runif (1, max = 500))
  for (j in seq (i, nodeListLength)) {
    node2 = nodeList [j]
    if (adMatrix [i, j] == 1) {
      from <- c (from, node)
      to <- c (to, node2)
      edgeInteger <- c (edgeInteger, round (runif (1, max = 500)))
      chars <- paste (sample (letters), collapse="")
      edgeChar <- c (edgeChar, chars)
      edgeFloat <- c (edgeFloat, runif (1, max = 500))
    } # if
  } # for j
} # for i

g <- initNodeAttribute (g, "nodeInteger", "integer", 0)
nodeData (g, nodes(g), "nodeInteger") <- nodeInteger
g <- initNodeAttribute (g, "nodeChar", "char", '')
nodeData (g, nodes(g), "nodeChar") <- nodeChar
g <- initNodeAttribute (g, "nodeFloat", "numeric", 0.0)
nodeData (g, nodes(g), "nodeFloat") <- nodeFloat

g <- initEdgeAttribute (g, "edgeInteger", "integer", 0)
edgeData (g, from, to, "edgeInteger") <- edgeInteger
g <- initEdgeAttribute (g, "edgeChar", "char", '')
edgeData (g, from, to, "edgeChar") <- edgeChar
g <- initEdgeAttribute (g, "edgeFloat", "numeric", 0.0)
edgeData (g, from, to, "edgeFloat") <- edgeFloat

window.name <- 'Random ADMatrix Graph'
if (window.name %in% as.character (getWindowList (cy)))
  destroyWindow (cy, window.name)
cw <- CytoscapeWindow (window.name, g)
displayGraph (cw)
layoutNetwork (cw, "jgraph-circle")
redraw (cw)

cw2 <- existing.CytoscapeWindow (window.name, copy.graph.from.cytoscape.to.R=TRUE)

tmp <- which(adjmat[which(rownames(mydata) == "1019-S") ,]!=0)

# select all the row names first and then column names
tmp1 <- adjmat[tmp,]
tmp2 <- tmp1[,tmp]


