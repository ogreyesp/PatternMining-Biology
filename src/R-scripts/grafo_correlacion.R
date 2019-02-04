rm(list = ls())
pacman::p_load("edgeR", "SummarizedExperiment", "igraph")

setwd("/home/antonio/Escritorio/Riñoncete")
condicion <- factor(c(rep("kich", 66), rep("kirc", 529), rep("kirp", 289)))
load("KICH.rda")
load("KIRC.rda")
load("KIRP.rda")
dfinal <- data.frame(kichfinal, kircfinal, kirpfinal)
rm(kichfinal)
rm(kirpfinal)
rm(kircfinal)
dfinal <- dfinal[-1,]

indx <- sapply(dfinal, is.factor)
dfinal[indx] <- lapply(dfinal[indx], function(x) as.numeric(as.character(x)))

dge <- DGEList(dfinal, group = condicion)
keep <- rowSums(cpm(dge)>1) >= 15 # Filtro para aquellos genes que tengan al menos cpm >= 1 en la mitad de las muestras
dge <- dge[keep,]
dge$samples$lib.size <- colSums(dge$counts) # Recalculo el tamaño de libreria

dge.n <- calcNormFactors(dge)

design <- model.matrix( ~ condicion)
colnames(design) <- c("kich", "kirc", "kirp")


design
dge.c <- estimateCommonDisp(dge.n)
dge.t <- estimateTagwiseDisp(dge.c)

fit <- glmQLFit(dge.t, design)
qlf <- glmQLFTest(fit, coef = 2:3)

topTags(qlf, n = 100)

genes <- topTags(qlf, n = 100)
genes <- rownames(genes$table)

# Voy a crear un grafo de coexpresión génica para ver si sale parecido a 
# minería de datos
significativos <- dfinal[genes,]

g <- graph.adjacency(as.matrix(as.dist(cor(t(significativos), method = "pearson"))), mode = "undirected", weighted = TRUE, diag = FALSE)

#Simplfy the adjacency object
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE)

#Colour negative correlation edges as blue
E(g)[which(E(g)$weight<0)]$color <- "darkblue"

#Colour positive correlation edges as red
E(g)[which(E(g)$weight>0)]$color <- "darkred"

#Convert edge weights to absolute values
E(g)$weight <- abs(E(g)$weight)

#Change arrow size
#For directed graphs only
#E(g)$arrow.size <- 1.0

#Remove edges below absolute Pearson correlation 0.8
g <- delete_edges(g, E(g)[which(E(g)$weight<0.8)])

#Assign names to the graph vertices (optional)
V(g)$name <- V(g)$name

#Change shape of graph vertices
V(g)$shape <- "sphere"

#Change colour of graph vertices
V(g)$color <- "skyblue"

#Change colour of vertex frames
V(g)$vertex.frame.color <- "white"

#Scale the size of the vertices to be proportional to the level of expression of each gene represented by each vertex
#Multiply scaled vales by a factor of 10
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
vSizes <- (scale01(apply(significativos, 1, mean)) + 1.0) * 10

#Amplify or decrease the width of the edges
edgeweights <- E(g)$weight * 2.0

#Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
mst <- mst(g, algorithm="prim")

#Plot the tree object
plot(
  mst,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.size=vSizes,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.6,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main="My first graph"
)

# PASO 3 IDENTIFICAR GRUPOS EN EL ARBOL BASANDOSE EN LA DISTANCIA ENTRE EJES
mst.communities <- edge.betweenness.community(mst, weights=NULL, directed=FALSE)
mst.clustering <- make_clusters(mst, membership=mst.communities$membership)
V(mst)$color <- mst.communities$membership + 1

par(mfrow=c(1,2))
plot(
  mst.clustering, mst,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.size=vSizes,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.6,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main="Cluster"
)

plot(
  mst,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.size=vSizes,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.6,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main="Genes coloreados"
)

