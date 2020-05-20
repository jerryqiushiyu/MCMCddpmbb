#' show the posterior probability for a unique voter pair to be in the same cluster.
#' @param voteCube_t the vote matrix in a period t
#' @import ggplot2
#' @importFrom reshape2 melt 
#' @importFrom stats hclust
#' @export
pair_probability_same_cluster_heatmap=function(voteCube_t,clusters_mat,statesName_t_iso, n_voter, n_item, font_size=5){
  # n_item=sum(rowSums(is.na(voteCube_t))!=ncol(voteCube_t))
  # n_voter=sum(colSums(is.na(voteCube_t))!=nrow(voteCube_t))
  
  
  clusters_mat=clusters_mat[,1:n_voter]
  clusters=vector(mode = 'list', length = nrow(clusters_mat))
  for(i in 1:length(clusters)){
    clusters[[i]]=clusters_mat[i,]
  }
  
  corr_store=matrix(NA, nrow=length(clusters[[1]]),ncol=length(clusters[[1]]) )
  
  
  rownames(corr_store)=statesName_t_iso
  colnames(corr_store)=statesName_t_iso
  for(i in 1:nrow(corr_store)){
    for(j in 1:ncol(corr_store)){
      same_cluster_counter=0
      for(t in 1:length(clusters)){
        if(clusters[[t]][i]==clusters[[t]][j]) same_cluster_counter=same_cluster_counter+1
      }
      corr_store[i,j]=same_cluster_counter/length(clusters)
    }
  }
  
  #corr=corr_store*2-1
  corr=corr_store
  #library(reshape2)
  cor_mat=corr
  melted_cormat <- melt(cor_mat)
  
  #corr
  #reorder
  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd,method = "ward.D2")
    return(cormat[rev(hc$order), rev(hc$order)])
    # cormat <-cormat[rev(hc$order), rev(hc$order)]
  }
  # # Get lower triangle of the correlation matrix
  # get_lower_tri<-function(cormat){
  #   cormat[upper.tri(cormat)] <- NA
  #   return(cormat)
  # }
  # # Get upper triangle of the correlation matrix
  # get_upper_tri <- function(cormat){
  #   cormat[lower.tri(cormat)]<- NA
  #   return(cormat)
  # }
  cormat=corr
  # Reorder the correlation matrix
  cormat <- reorder_cormat(cormat)
  #upper_tri <- get_upper_tri(cormat)
  # Melt the correlation matrix
  #melted_cormat <- melt(upper_tri, na.rm = TRUE)
  melted_cormat <- melt(cormat, na.rm = TRUE)
  # library(ggplot2)
  # Create a ggheatmap
  ggheatmap <- ggplot(melted_cormat, aes(Var1, Var2, fill = value))+
    geom_tile(color = "white")+
    # scale_fill_gradient2(low = "white", high = "red", mid = "white", 
    #                      midpoint = 0, limit = c(-1,1), space = "Lab", 
    #                      name="Same Cluster") +
    scale_fill_gradient2(low = "white", high = "red", mid="orange",
                         midpoint = 0.5, limit = c(0,1), space = "Lab",
                         name="Probability") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     size = font_size, hjust = 0.5),
          axis.text.y = element_text(size = font_size),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())+
    coord_fixed()
  # Print the heatmap
  print(ggheatmap)
}
