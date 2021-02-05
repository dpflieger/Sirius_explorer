library(tidyverse)
library(data.table)

json_test <- "~/Projects/Gaquerel/Elser/Sirius/4_Tissue+Exudates-sirius301120_38/191_Tissue+Exudates-sirius301120_1328/trees/C23H38N2O2_[M+H]+.json"
res <- jsonlite::fromJSON(json_test, flatten = TRUE)

fragments_mz <- res$fragments$mz
names(fragments_mz) <- res$fragments$id

fragments_mz["0"]

setDT(res$losses)

names(res$losses)
names(res$fragments)

res$losses[, source_mz := fragments_mz[as.character(source)]]
res$losses[, targets_mz := fragments_mz[as.character(target)]]

res$losses[, delta_mz := source_mz - targets_mz]


test <- fread("/Volumes/4TB/Users/dpflieger/Projects/Gaquerel/Elser/Sirius/4_Tissue+Exudates-sirius301120_38/formula_identifications.tsv")
colnames(test)

# 

library(treemap)
library(data.tree)

test <- res$losses

test$pathString <- paste("losses",
                         test$source,
                         test$target,
                         sep = "/"
                         )

print(as.Node(test), "molecularFormula")
plot(as.Node(test))

library(igraph)
as.igraph(as.Node(test))


library(collapsibleTree)

collapsibleTree(as.Node(test), nodeSize = )


collapsibleTree(df = res$losses, hierarchy = c("source", "target"))

species <- read.csv("https://apps.fs.usda.gov/fia/datamart/CSV/REF_SPECIES_GROUP.csv")
print(res$losses)

FromDataFrameTable()


library(ggtree)


# Igraph network ----------------------------------------------------------
library(igraph)
library(d3r)

nodes = res$fragments
edges = res$losses[, c("source", "target", "molecularFormula")]
network <- graph_from_data_frame(d = edges, directed = T, vertices = nodes)

E(network)$label.color <- rainbow(nrow(nodes) + 1)

set.seed(1)

l <- layout_as_tree(network)
#ll <- layout_(network, with_dh(weight.edge.lengths = edge_density(network)/1000)))

par(mar=c(0,0,0,0))
plot(network, 
     vertex.shape = "rectangle", vertex.size = 25, vertex.size2 = 10,
     vertex.color = "lightblue", 
     vertex.label=V(network)$molecularFormula, 
     vertex.label.cex = 0.8,
     vertex.label.font = 2,
     arrow.mode = 0,
     #edge.color = rainbow(n = 31), 
     edge.label = res$losses$molecularFormula,
     edge.width = 1,
     edge.arrow.size = 0.2,
     #edge.arrow.width = 0.5,
     edge.label.cex = 0.8,
     edge.lty = 2,
     layout = l, rescale=T)

# Test all layout igraph --------------------------------------------------

# See all layout
layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1] 
# Remove layouts that do not apply to our graph.
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama", layouts)]

par(mfrow=c(3,3), mar=c(1,1,1,1))
for (layout in layouts) {
    print(layout)
    l <- do.call(layout, list(network)) 
    plot(network, edge.arrow.mode=0, layout=l, main=layout,  vertex.shape = "vrectangle", 
         vertex.color = "lightblue", 
         edge.color = "black", 
         vertex.label=V(network)$molecularFormula) 
}


# visNetwork --------------------------------------------------------------
library("visNetwork")

nodes <- setDT(res$fragments[, c("id", "label")])

setnames(nodes, c("molecularFormula"), c("label"))
edges <- setDT(res$losses[, c("source", "target", "molecularFormula")])
setnames(edges, c("source", "target", "molecularFormula"), c("from", "to", "label"))

t <- visNetwork(nodes, edges) %>% 
    visIgraphLayout(layout = "layout_as_tree")
    #visHierarchicalLayout(direction = "UD") # to have always the same network     
t$x$nodes$y <- -t$x$nodes$y
t$x$nodes$shape <- "box"
t
