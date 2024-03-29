#' generate igraph
#'
#' @param network edge list generated by kimono
#' @param directed if false the method will collaps all edges and remove free nodes as well as single directional edges.
#' @return igraph
to_igraph <- function(network, directed = TRUE){

  if(sum(c('target','predictor','value','r_squared','mse','predictor_layer','target_layer') %in% colnames(network)) != 7){
    if(sum(c("target","predictor","mean_value","sd_value","mean_rsq","sd_rsq","mean_mse","sd_mse","sel_freq","predictor_layer","target_layer") %in% colnames(network)) != 11){
      warning('data format unknown')
      return()
    }
  }

  #we only consider r_squared
  colnames(network)[colnames(network) %in% c('mean_value','value')] <- 'value'
  colnames(network)[colnames(network) %in% c('mean_rsq','r_squared')] <- 'r_squared'

  #filter intercept and 0 edges
  tmp_network <- network  %>%
    filter(value !=0 ) %>%  #filter edge effect size = 0
    filter(predictor != '(Intercept)') #filter intercepts and intercept only models

  #create a unique vertex list
  tmp_network$predictor <- paste0(tmp_network$predictor_layer,"__",tmp_network$predictor)
  tmp_network$target <- paste0(tmp_network$target_layer,"__",tmp_network$target)

  #preprocess lable
  id <- unique(c(tmp_network$predictor, tmp_network$target ) )
  id_tmp <- do.call(rbind, strsplit(id,"__")) %>% data.table
  id <- cbind(id_tmp,id)
  tmp <- unique(tmp_network[,c('target','r_squared')])
  colnames(tmp)[1] <- 'id'
  id <- id[tmp, on = 'id', r_squared := i.r_squared ]

  #node information
  nodes <- data.frame(id = id$id,
                      label = id$V2,
                      data_layer = id$V1,
                      r_squared = id$r_squared )


  #generate igraph
  ig_kimono <- graph_from_data_frame(tmp_network[,c("predictor","target","value")], vertices=nodes, directed = TRUE)

  if(!directed){
    ig_kimono <- as.undirected(ig_kimono,mode = c("mutual") )
    isolated = which(degree(ig_kimono)==0)
    ig_kimono <- delete.vertices(ig_kimono, isolated)
  }

  ig_kimono
}

#' generate igraph
#'
#' @param ig_kimono igraph of kimono network
#' @param title plot title
#' @return igraph
plot_kimono <- function(ig_kimono,title=''){

  layers <- V(ig_kimono)$data_layer
  color <- rainbow(length(layers), s = 0.4)

  plot(ig_kimono,
       main = title,
       edge.curved=0,
       vertex.color = color[layers %>% as.factor %>% as.numeric ],
       vertex.frame.color="white",
       vertex.label = V(ig_kimono)$label,
       vertex.label.color='black',
       vertex.label.cex=.7,
       rescale=T
  )

  legend(x=-1.5, y=-1.1, layers %>% unique, pch=21,
         col="#777777", pt.bg= color[layers %>% unique%>% as.factor %>% as.numeric ] , pt.cex=2, cex=.8, bty="n", ncol=1)
}



