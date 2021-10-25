#' loading mapping files
#'
#' @param df data frame or mapping file
#' @param layers - default NA, if no colnames please specify your layer names here
#' @return data frame containing all mappings
load_mapping <- function(df, layers=NA){
  if(length(layers) != 2 ){
    layers <- colnames(df)
  }
  map <- as.data.table(df)
  map <- cbind(map, layers[1], layers[2])
  colnames(map) <- c('A','B','layer_A','layer_B')
  map
}

#' creating a prior network on basis of the prior_map
#'
#' @param prior_map data frame - edge list -wir column names  'A','B','layer_A','layer_B'
#' @return igraph network
create_prior_network <- function(prior_map) {

  #cleaning
  idx_rm <- as.character(prior_map$A) %in% c('',NA) |
    as.character(prior_map$B) %in% c('',NA)
  prior_map <- prior_map[!idx_rm,]

  #rm duplicated mapings
  prior_map <- distinct(prior_map)

  #nodes
  A <- distinct(prior_map,A,layer_A)
  colnames(A) <- c('id_name','layer')
  B <- distinct(prior_map,B,layer_B)
  colnames(B) <- c('id_name','layer')
  nodes <- distinct(rbind(A,B))
  nodes <- cbind('id' = paste(nodes$layer,nodes$id_name,sep = '___'),
                 'name' = paste(nodes$layer,nodes$id_name,sep = '___'),
                 nodes)

  #links
  links <- data.table( from = paste(prior_map$layer_A,prior_map$A,sep='___'),
                       to = paste(prior_map$layer_B,prior_map$B,sep='___'),
                       relation = paste(prior_map$layer_A,prior_map$layer_B,sep='___')

  )

  graph_from_data_frame(links, directed = FALSE, vertices = nodes)
}
