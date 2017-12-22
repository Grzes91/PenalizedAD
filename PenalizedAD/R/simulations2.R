#' @title MAES algorithm
#' @description Method to obtain parameter estimates for the MAES algorithm. Slelection of regularization parameters is performed
#' automatically based on the modified BIC criteria
#' @importFrom mclust adjustedRandIndex
#' @import mvtnorm
#' @param data numeric matrix of size n and dimension P. The data on which the model is constructed.
#' @param J numeric(1). Number of Gaussian components used in the mixture to fit the model.
#' @param labs numeric vector of length equal to size of the data, default all labels are equal to 0. Used for computation of the classification results.
#' @return resulting list with respective entries
#' \itemize{
#' \item means - matrix (JxP) of component means by row
#' \item eigenvectors - list of J component matrices containing eigenvectors by column
#' \item eigenvalues - list of J component diagonal matrices containing eigenvalues
#' \item proportions - component proportions \eqn{\pi}
#' \item gamma1 - \eqn{\gamma_1} mean regularization parameter
#' \item gamma2 - \eqn{\gamma_2} eigenvalue regularization parameter
#' \item ARI - adjusted Rand Index for classification results
#' \item pred_classes - model predicted classes
#' }
#' @references Working paper on a Statistical Learning method for Model-Independent searches for New Physics - AMVA4NewPhysics ITN
#' @examples
#' data(bcg_data)
#'
#' #Obtain model parameter estimates (might take about a minute)
#' f <- MAES(bcg_data, 2, attr(bcg_data, "label"))
#' @export
MAES=function(data,J,labs = 0){
  n=dim(data)[1]
  rand=0
  iter=1
  n_iter=50
  p=dim(data)[2]
  #fit models with different penalties for eigenvalues shrinkage
  possible_lam=c(4,2,1,0.2,0.02,0.002)
  s=sapply(possible_lam, function (lam) bic_mine(EM_bis(data,J,49,lam,lambda_mean=0),n))
  E=EM_bis(data,J,n_iter,lambda = possible_lam[which.min(s)],lambda_mean=0)

  #Run multiply times for the selected parameters
  for(i in 1:10){G=EM_bis(data,J,n_iter,possible_lam[which.min(s)],lambda_mean=0)
  if(G[[4]][n_iter]>E[[4]][n_iter])E=G}

  #Again try to find best eigenvalue regularization param, but using adaptive penalty
  lamev=c(4,2,1,0.2,0.02,0.002)
  pv=matrix(0,nrow=2,ncol=length(lamev))
  Warm_start=E
  for(lam in 1:length(lamev)){
    Fd=EM_adaptive(data,J,n_iter,lamev[lam],lambda_mean = 0,adapt = E[[3]][[3]],pi=Warm_start[[3]][[4]],mean=Warm_start[[3]][[1]],Qm=Warm_start[[3]][[2]],Dm=Warm_start[[3]][[3]])
    pv[1,lam]=bic_mine(Fd,n)
    pv[2,lam]=sum(diag(Fd[[3]][[3]][[2]])==Fd[[3]][[3]][[2]][p,p])
    Warm_start=Fd
  }

  #Based on the weights E the adaptive models are fitted
  #a grid search for best lambda values is performed
  #where a grid values are selected based on the above pre-fit
  pv0=pv
  lambda_eigen=lamev[which.min(pv0[1,])]*c(0.5,0.8,1,1.2,1.4,1.8,2.2,3.2)
  lambda_mean=c(0,10,20,25,30,35,40,45)
  pv=matrix(0,nrow=7,ncol=length(lambda_eigen)*length(lambda_mean))
  Warm_start=E
  E_best=E
  bic_min=Inf
  ni=n_iter
  for(l1 in 1:length(lambda_eigen)){
    for(l2 in 1:length(lambda_mean)){
      if(l2!=lambda_mean[1])ni=15 else ni=n_iter
      Fd=EM_adaptive(data,J,ni,lambda = lambda_eigen[l1],lambda_mean = lambda_mean[l2],pi=E[[3]][[4]],adapt = E[[3]][[3]],mean=Warm_start[[3]][[1]],Qm=Warm_start[[3]][[2]],Dm=Warm_start[[3]][[3]])
      pv[1,(l1-1)*length(lambda_mean)+l2]=bic_mine(Fd,n)
      pv[2,(l1-1)*length(lambda_mean)+l2]=sum(diag(Fd[[3]][[3]][[1]])==Fd[[3]][[3]][[1]][p,p])
      pv[3,(l1-1)*length(lambda_mean)+l2]=sum(diag(Fd[[3]][[3]][[2]])==Fd[[3]][[3]][[2]][p,p])
      pv[4,(l1-1)*length(lambda_mean)+l2]=sum(Fd[[3]][[1]]==0)
      pv[5,(l1-1)*length(lambda_mean)+l2]=lambda_eigen[l1]
      pv[6,(l1-1)*length(lambda_mean)+l2]=lambda_mean[l2]
      pv[7,(l1-1)*length(lambda_mean)+l2]=Fd[[5]][1]
      Warm_start=Fd
      if(pv[1,(l1-1)*length(lambda_mean)+l2]<bic_min & Fd[[5]][1]>0){
        E_best=Fd
        bic_min=pv[1,(l1-1)*length(lambda_mean)+l2]
      }
    }
    Warm_start=E
  }

  rand=miss_error(data,J,E=E_best,labs)
  E_best[[7]]=rand
  E_best[[8]]=classes(data,J,E=E_best,labs)
  E_best[[1]] <- E_best[[3]][[1]]
  E_best[[2]] <- E_best[[3]][[2]]
  E_best[[4]] <- E_best[[3]][[4]]
  E_best[[3]] <- E_best[[3]][[3]]
  names(E_best) <- c("means", "eigenvectors","eigenvalues","proportions","gamma2","gamma1", "ARI", "pred_classes")
  return(E_best)
}

#' @title Penalized anomaly detection
#' @description Implementation of the PAD algorithm that enables to find the anomalous component in the experimental data.
#' @importFrom mclust adjustedRandIndex
#' @import mvtnorm
#' @param data_exp numeric matrix, experimental data which is tested if contain anomaly
#' @param data data.frame. background data against which anomalies are searched
#' @param E_best MAES result. The background model - starting point for the algorithm, possible output model if signal is not found
#' @param minres numeric(1) I don't know
#' @param gamma_mean numeric(1) regularization parameter \eqn{\gamma_1} for the mean parameters penalty
#' @param gamma_eigen numeric(1) regularization parameter \eqn{\gamma_2} for the component covariance eigenvalue penalty
#' @param n_iter numeric(1) default equal to 20, maximum number of iteration for the algorithm outer loop
#' @param minbic numeric(1) the minimal BIC criteria reached so far.
#' @references Working paper on a Statistical Learning method for Model-Independent searches for New Physics - AMVA4NewPhysics ITN
#' @return resulting list with respective entries
#' \itemize{
#' \item means - matrix (JxP) of component means by row
#' \item eigenvectors - list of J component matrices containing eigenvectors by column
#' \item eigenvalues - list of J component diagonal matrices containing eigenvalues
#' \item proportions - component proportions \eqn{\pi}
#' \item gamma1 - \eqn{\gamma_1} mean regularization parameter
#' \item gamma2 - \eqn{\gamma_2} eigenvalue regularization parameter
#' }
#'
#' @examples
#' #Load the datasets
#' data(bcg_data)
#' data(experiment_data)
#'
#' #1 fit the background model as a possible result if no anomaly is found
#' E_bcg=MAES(bcg_data,2,attr(bcg_data, "label"))
#'
#' #2 fit the signal
#'
#' #Iterate between fiting the bacground on data and fitting the signal on data2
#' minbic <- Inf
#' sig <- PAD(E_bcg, experiment_data,gamma_eigen = 1,
#'             gamma_mean = 4,minbic,bcg_data)
#'
#' #Resulting component proportions
#' sig[[4]]
#' @export
PAD=function(E_best,data_exp,gamma_eigen,gamma_mean,n_iter=20,minbic,data)#allows to find a signal in data_exp set
{
  minres=E_best
  E_best[[3]] <- list(E_best[[1]],E_best[[2]],E_best[[3]],E_best[[4]])
  n=nrow(data_exp)
  Fd=E_best
  G=E_best
  J=length(Fd[[3]][[3]])
  p=dim(data_exp)[2]
  Fd0=list()
  add=data_exp[n-1,]
  probs=prob_m1(data_exp,E_best[[3]][[1]],E_best[[3]][[2]],E_best[[3]][[3]],E_best[[3]][[4]])#2xn matrix
  Fd[[3]][[4]]=c(0.9*Fd[[3]][[4]],0.1)
  Fd[[3]][[3]][[J+1]]=diag(0.6,p)
  Fd[[3]][[2]][[J+1]]=diag(1,p)
  minbic=minbic
  minres=minres
  for(i in 1:10){
    ad=apply(cbind(diag(Fd[[3]][[3]][[1]]),diag(Fd[[3]][[3]][[2]])),1,mean)
    res=EM2_bis(data=data_exp,J=3,n_iter = n_iter,lambda = gamma_eigen,lambda_mean=gamma_mean,means =  Fd[[3]][[1]], Qm=Fd[[3]][[2]],Dm=Fd[[3]][[3]],pi =  Fd[[3]][[4]],adapt=ad,add=add,probs = probs)
    #model selection
    Fd0[[1]]=res[[1]][1:2,]
    Fd0[[2]]=list()
    Fd0[[3]]=list()
    Fd0[[2]][[1]]=res[[2]][[1]]
    Fd0[[2]][[2]]=res[[2]][[2]]
    Fd0[[3]][[1]]=res[[3]][[1]]
    Fd0[[3]][[2]]=res[[3]][[2]]
    Fd0[[4]]=res[[4]][1:2]/sum(res[[4]][1:2])
    Fd=EM_adaptive(data,J,n_iter,lambda = E_best[[5]],lambda_mean = E_best[[6]],pi=Fd0[[4]],adapt = E_best[[3]][[3]],mean=Fd0[[1]],
                   Qm=Fd0[[2]],Dm=Fd0[[3]],add=res[[1]][J+1,])
    probs=prob_m1(data_exp,Fd[[3]][[1]],Fd[[3]][[2]],Fd[[3]][[3]],Fd[[3]][[4]])#2xn matrix
    Fd[[3]][[2]][[3]]=res[[2]][[3]]
    Fd[[3]][[3]][[3]]=res[[3]][[3]]
    Fd[[3]][[4]]=c((res[[4]][[1]]+res[[4]][[2]])*Fd[[3]][[4]],res[[4]][[3]])
    add=res[[1]][3,]
  }
  G[[3]]=res
  G[[4]]=res[[5]]
  if(bic_mine(G,n)<minbic){
    minbic=bic_mine(G,n)
    minres=G
    #print(minbic)
  }
  minbic<<-minbic
  minres[[5]] <- minres[[3]]
  minres[[6]] <- minres[[4]]
  minres[[1]] <- minres[[3]][[1]]
  minres[[2]] <- minres[[3]][[2]]
  minres[[4]] <- minres[[3]][[4]]
  minres[[3]] <- minres[[3]][[3]]
  minres[[7]] <- NULL
  minres[[8]] <- NULL
  names(minres) <- c("means", "eigenvectors","eigenvalues","proportions","gamma2","gamma1")
  return(minres)
}

##########
#Other important functions to fit a model


#Compute the adjusted rand index
miss_error=function(data,J,E,labs){
  class = classes(data,J,E,labs)
  adjustedRandIndex(labs,class)
}

classes=function(data,J,E,labs){
  class=sapply(1:J,function(j) E[[3]][[4]][j]*dmvnorm(data,E[[3]][[1]][j,],E[[3]][[2]][[j]]%*%E[[3]][[3]][[j]]%*%t(E[[3]][[2]][[j]])))
  apply(class,1,which.max)
}
