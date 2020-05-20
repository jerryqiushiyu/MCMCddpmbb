#' show the mode of the numbers of clusters across iterations.
#' @export
Modecluster<-function(nclust_store=dynResult$nclust_store_data){
  nclust_mode=rep(NA, dim(nclust_store)[2])
  for(t in 1:dim(nclust_store)[2]){
    nclust_mode[t]=as.numeric(names(table(nclust_store[,t]))[table(nclust_store[,t])==max(table(nclust_store[,t]))])
  }
  return(nclust_mode)
}
