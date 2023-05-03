# Function to plot adjacency graph from ligand receptor interaction tables, obtained from LIANA
own_interaction_graph <- function(liana_trunc){
  lrp <- liana_trunc[,c(3,4)]
  ligands <- unique(lrp$ligand.complex)
  receptors <- unique(lrp$receptor.complex)
  adj_matrix <- matrix(0
                       ,nrow = length(ligands) + length(receptors)
                       ,ncol = length(ligands) + length(receptors)
  )
  dimnames(adj_matrix) <- list(c(ligands,receptors)
                               ,c(ligands,receptors)
  )
  for(i in 1:nrow(lrp)){
    idx_row <- which(lrp$ligand.complex[i] == rownames(adj_matrix))
    idx_col <- which(lrp$receptor.complex[i] == colnames(adj_matrix))
    adj_matrix[idx_row, idx_col] <- 1
  }
  
  graph <- igraph::graph.adjacency(adj_matrix,mode = "undirected")
  
  V(graph)$color <- {
    color <- c(rep(c("#C3C6C2", "#231F20")[1]
                   ,length(ligands))
               ,rep(c("#C3C6C2", "#231F20")[2]
                    ,length(receptors)))
    names(color) <- c(ligands
                      ,receptors)
    color <- color[colnames(adj_matrix)]
    color
  }
  
  V(graph)$label.color <- {
    color <- c(rep(c("black", "white")[1]
                   ,length(ligands))
               ,rep(c("black", "white")[2]
                    ,length(receptors)))
    names(color) <- c(ligands
                      ,receptors)
    color <- color[colnames(adj_matrix)]
    color
  }
  
  set.seed(42)
  # pdf(save_path, 9, 9)
  igraph::plot.igraph(graph, vertex.size = 13, edge.width = 3, 
                      vertex.frame.color = NA, vertex.label.family = "Helvetica", 
                      layout=layout_with_fr)
  # dev.off()
}
