#' Extract the Maximum A Posteriori cluster estimator from the posterior sample.
#'
#' @export
MAPcluster<-function(voteCube=voteCube, nclust_store=dynResult$nclust_store_data, clusterSample=dynResult$c_store, dataStructure=dataStructure){
  T<-nrow(dataStructure)
  MAPcluster<-vector(mode = "list",length = T)
  nclust_mode<-Modecluster(nclust_store=nclust_store)
  for(t in 1:T){
    n_state<-dataStructure[t,1]
    n_item<-dataStructure[t,2]
    voteCube_t<-voteCube[1:n_item,1:n_state,t]
    cluster_t<-clusterSample[,1:n_state,t]
    cluster_t<-cluster_t[nclust_store[,t]==nclust_mode[t],]
    iterations<-nrow(cluster_t)
    loglike_t<-rep(NA, iterations )
    for(iter in 1:iterations){
      cluster_unique<-unique(cluster_t[iter,])
      cluster_marginal_log<-rep(NA,length(cluster_unique))
      cluster_unique_prob<-matrix(NA, nrow=n_item, ncol=length(cluster_unique))
      
      for(clus in 1:length(cluster_unique)){
        cluster_marginal_log[clus]<-log(sum(cluster_t[iter,]==cluster_unique[clus])/n_state)
        #The line below computes the mean of the posterior beta distributions with beta(1,1) prior and clustered vote records.
        cluster_unique_prob[,clus]<-rowMeans(cbind(voteCube_t[,cluster_t[iter,]==cluster_unique[clus]],1,0), na.rm = T)
        
      }
      #compute the posterior likelihood as the cluster marginal * bernoulli likelihood with cluster wise probability vectors computed above
      loglike_t_iter<-0
      for(ii in 1:n_state){
        vote_ii<-voteCube_t[,ii]
        cluster_prob_ii<-cluster_unique_prob[,cluster_t[iter, ii]==cluster_unique]
        loglike_t_iter<-loglike_t_iter+cluster_marginal_log[cluster_t[iter, ii]==cluster_unique]+sum(vote_ii*log(cluster_prob_ii),na.rm = T)+sum((1-vote_ii)*log((1-cluster_prob_ii)),na.rm = T)
      }
      loglike_t[iter]<-loglike_t_iter
      
    }
    if(length(which(loglike_t==max(loglike_t)))==1){
      MAPcluster[[t]]<-cluster_t[which(loglike_t==max(loglike_t)),]
      
    }else{
      MAPcluster[[t]]<-cluster_t[which(loglike_t==max(loglike_t))[1],]
      
    }
    
  }
  return(MAPcluster)
  
}


