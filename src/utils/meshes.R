## utils

import_mesh_data <- function(path) {
  mesh_data <- list(
    nodes = as.matrix(read.csv(paste(path, "points.csv", sep = ""))[, -1]),
    edges = as.matrix(read.csv(paste(path, "edges.csv", sep = ""))[, -1]),
    elements = as.matrix(read.csv(paste(path, "elements.csv", sep = ""))[, -1]),
    neigh = as.matrix(read.csv(paste(path, "neigh.csv", sep = ""))[, -1]),
    boundary = as.matrix(read.csv(paste(path, "boundary.csv", sep = ""))[, -1])
  )
  return(mesh_data)
}


## mesh is an fdaPDE1 mesh
evaluate_field <- function(grid, f_at_nodes, mesh) {
  
  FEMbasis = create.FEM.basis(mesh)
  
  FEMfunction = FEM(f_at_nodes, FEMbasis)
  f_at_grid <- eval.FEM(FEMfunction, grid)
  
  return(f_at_grid)
}


## mesh generators

unit_segment <- function(n_nodes) {
  x <- seq(0, n_nodes-1)/(n_nodes-1)
  mesh_data <- list(
    nodes = as.matrix(x),
    elements = matrix(rep(0, 2*n_nodes), ncol = 2),
    boundary = as.matrix(rep(0, n_nodes))
  )
  return(mesh_data)
}

unit_square <- function(n_nodes) {
  x <- y <- seq(0, 1, length = sqrt(n_nodes))
  grid <- meshgrid(x, y)
  nodes <- data.frame(as.vector(grid$X), as.vector(grid$Y))
  mesh <- fdaPDE::create.mesh.2D(nodes)
  mesh_data <- list(
    nodes = as.matrix(mesh$nodes),
    edges = as.matrix(mesh$edges),
    elements = as.matrix(mesh$triangles),
    neigh = as.matrix(mesh$neighbors),
    boundary = as.matrix(as.numeric(mesh$nodesmarkers))
  )
  return(mesh_data)
}
