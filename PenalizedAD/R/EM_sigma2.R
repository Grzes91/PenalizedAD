########
#The functions required by the EM algorithm


#Compute the likelihood given the data and parameters
log_likelihood=function(data,tau,Q,D,pi){
  d1=sapply(1:dim(tau)[1],function(k) dmvnorm(data,tau[k,],Q[[k]]%*%D[[k]]%*%t(Q[[k]])))
  sum(log(d1%*%pi))
}

#Used to calculate probability of given the fixed backgroun parameters
pb=function(data,mu,Q,D,pi){
  d1=sapply(1:length(pi),function(k) dmvnorm(data,mu[k,],Q[[k]]%*%D[[k]]%*%t(Q[[k]]))*pi[k])
  apply(d1,1,sum)
}
####################list(best_res,worst_res,last_res,l,l2,l3)
##All the functions EM_* run the EM algortihm and return a listt
#[[1]] list of parameters giving highest likelihood
#[[2]] list of parameters giving lowest likelihood
#[[3]] list of parameters for the converged algorithm
#[[4]] vector of likelihood values
#[[5]] lambda-eigen
#[[6]] lambda-mean
#for list 1-3
#[[.]][[1]] mean parameters
#[[.]][[2]] eigenvecotrs parameters
#[[.]][[3]] eigenvalues parameters
#[[.]][[4]] pi


#Computes the E step
#Returns the matrix of posterior probabilities
#matrix nxJ pf probabilities
E_bis=function(data,mu,Q,D,pi){
  d1=sapply(1:length(pi),function(k) dmvnorm(data,mu[k,],Q[[k]]%*%D[[k]]%*%t(Q[[k]]))*pi[k])
  apply(d1,1,function(x) x/sum(x))
}

#The bagkround fit using EM alg
EM_bis=function(data,J,n_iter,lambda,lambda_mean=0,start=F,means=0,Q=0,D=0){
  if(lambda==0)return(EM_only_mean(data,J,n_iter,lambda_mean))
  p=dim(data)[2]
  n=dim(data)[1]
#  km=kmeans(data,J)
  if(start==F){
    means=matrix(rnorm(J*p,sd=1),nrow=J,ncol=p)
    means=apply(means,2,function(k) k+seq(-1,1,length.out = J))
    Q=list()
    D=list()
    for(i in 1:J){
      D[[i]]=diag(1,p)
      Q[[i]]=diag(1,p)
    }
  }
#initialize parameters
  l=1:n_iter
  l2=1:n_iter
  l3=1:n_iter
  pi=rep(1/J,J)

  best_res=list(means,Q,D,pi)
  worst_res=list(means,Q,D,pi)
  last_res=list(means,Q,D,pi)
  like_min=+Inf
  like_max=-Inf
  for(iter in 1:n_iter){#A main loop
    if(sum(pi<0.02)>0){reset=(which(pi<0.02))[1]#resetting small components
    pi=rep(1/J,J)
    means=data[sample(1:n,J),]
    D[[reset]]=diag(1.5,p)
    }
    posterior=E_bis(data,means,Q,D,pi)#matrix Jxn
    pi=apply(posterior,1,mean)
    denom=apply(posterior,1,sum)
    cov=list()
    for(j in 1:J)cov[[j]]=Q[[j]]%*%D[[j]]%*%t(Q[[j]])
    m_cov=cov[[1]]
    for(j in 2:J)m_cov=diag(pmax(diag(m_cov),diag(cov[[j]])))
    means_wave=posterior%*%data
    means_wave=t(sapply(1:J,function(k) means_wave[k,]/denom[k]))
    norm_mw=sqrt(colSums(apply(means,2,function(c) c^2)))#vector of norms of hat mu
    means_temp=matrix(0,nrow=J,ncol=dim(data)[2])#means to be
    limit=sqrt(apply((posterior%*%data)^2,2,sum))
    for(i in 1:p){
         for(j in 1:J)
            if(limit[i]>lambda_mean*m_cov[i,i] & norm_mw[i]>0){
              means_temp[j,i]=means_wave[j,i]-lambda_mean*m_cov[i,i]*means[j,i]/denom[j]/norm_mw[i]
            }
    }

    Sigma=lapply(1:J, function(j) apply(sapply(1:dim(data)[1],function(i) posterior[j,i]*(data[i,]-means[j,])%*%t(data[i,]-means[j,])),1,sum)/(denom[j]))
      for(i in 1:J){

      E=tryCatch({eigen(matrix(Sigma[[i]],nrow = p))},error=function(e){return(list(best_res,worst_res,last_res,l,-1,l3))})
      Q[[i]]=E$vectors
      #Normal shrinkage
      E$values[E$values<0.005]=0.005
      diagonala=(-denom[i]+sqrt(denom[i]^2+8*lambda*denom[i]*E$values))/4/lambda
      diagonala[diagonala<0.005]=0.005

    #Shrinkage to epsilon based on hypothesis test
      for(kol in 2:p){
            eps_1=mean(diagonala[kol:p])
            eps_2=diagonala[p]
            eps_3=diagonala[kol]
            if(eps_3/eps_1<1-sqrt(2/n)*qnorm(0.025/(kol)) | eps_2/eps_1>1+sqrt(2/n)*qnorm(0.025/(kol))){
                  diagonala[kol:p]=eps_1
                  D[[i]]=diag(diagonala)
         break
       }# else print (eps_3/eps_1)
     }
    }
    means=means_temp
    like=log_likelihood(data,means,Q,D,pi)#-1*sum(sapply(D, function(k) sum(k)))
    if(iter>10 & like_min>like){like_min=like
      worst_res=list(means,Q,D,pi)}
    if(iter>10 & like_max<like){like_max=like
    best_res=list(means,Q,D,pi)}
    l[iter]=like
    l2[iter]=D[[1]][p-1,p-1]
    l3[iter]=D[[2]][p-1,p-1]
  }
#  print(pi)
  last_res=list(means,Q,D,pi)
  list(best_res,worst_res,last_res,l,lambda,lambda_mean)
}

EM_only_mean=function(data,J,n_iter,lambda_mean=0){
  p=dim(data)[2]
  n=dim(data)[1]
  means=matrix(rnorm(J*p,sd=1),nrow=J,ncol=p)
  means=apply(means,2,function(k) k+seq(-1,1,length.out = J))
  Q=list()
  D=list()
  l=1:n_iter

  for(i in 1:J){
    D[[i]]=diag(1,p)
    Q[[i]]=diag(1,p)
  }

  pi=rep(1/J,J)

  best_res=list(means,Q,D,pi)
  worst_res=list(means,Q,D,pi)
  last_res=list(means,Q,D,pi)
  like_min=+Inf
  like_max=-Inf
  for(iter in 1:n_iter){
    if(sum(pi<0.02)>0){reset=which(pi<0.02)
    pi=rep(1/J,J)
    for(r in reset){
        means[r,]=data[sample(1:n,1),]
    D[[r]]=diag(1.5,p)}
    }
    posterior=E_bis(data,means,Q,D,pi)#matrix Jxn
    pi=apply(posterior,1,mean)
    denom=apply(posterior,1,sum)
    means_wave=posterior%*%data
    #   means=t(sapply(1:J,function(k) means_wave[k,]/denom[k]))
    cov=list()
    for(j in 1:J)cov[[j]]=Q[[j]]%*%D[[j]]%*%t(Q[[j]])
    m_cov=diag(pmax(diag(cov[[1]]),diag(cov[[2]])))
    means_wave=posterior%*%data
    means_wave=t(sapply(1:J,function(k) means_wave[k,]/denom[k]))
    norm_mw=sqrt(colSums(apply(means,2,function(c) c^2)))#vector of norms of hat mu
    means_temp=matrix(0,nrow=J,ncol=dim(data)[2])#means to be
    limit=sqrt(apply((posterior%*%data)^2,2,sum))
    for(i in 1:p){
      for(j in 1:J)
        if(limit[i]>lambda_mean*m_cov[i,i] & norm_mw[i]>0){
          means_temp[j,i]=means_wave[j,i]-lambda_mean*m_cov[i,i]*means[j,i]/denom[j]/norm_mw[i]
        }
    }

    Sigma=lapply(1:J, function(j) apply(sapply(1:dim(data)[1],function(i) posterior[j,i]*(data[i,]-means[j,])%*%t(data[i,]-means[j,])),1,sum)/(denom[j]))
    for(i in 1:J){
      D[[i]]=matrix(Sigma[[i]],nrow = p)
    }
    means=means_temp
    like=log_likelihood(data,means,Q,D,pi)#-1*sum(sapply(D, function(k) sum(k)))
    if(iter>10 & like_min>like){like_min=like
    worst_res=list(means,Q,D,pi)}
    if(iter>10 & like_max<like){like_max=like
    best_res=list(means,Q,D,pi)}
    l[iter]=like
  }
 # print(pi)
  last_res=list(means,Q,D,pi)
  list(best_res,worst_res,last_res,l,pi)
}

EM_adaptive=function(data,J,n_iter,lambda,lambda_mean=0,adapt,mean,Qm,Dm,pi,add=NA){
  if(lambda==0)return(EM_only_mean(data,J,n_iter,lambda_mean))
  p=dim(data)[2]
  n=dim(data)[1]
  means=mean
  Q=Qm
  D=Dm
  l=1:n_iter
  l2=1:n_iter
  l3=1:n_iter
  pi=pi
  pi_ad=pi#the weight for components mean weigthing

  best_res=list(means,Q,D,pi)
  worst_res=list(means,Q,D,pi)
  last_res=list(means,Q,D,pi)
  like_min=+Inf
  like_max=-Inf
  count=0
  for(iter in 1:n_iter){
    if(sum(pi<0.02)>0){reset=which(pi<0.02)
    count=count+1
    if(count==10){print("Too many resets")
      return(list(best_res,worst_res,last_res,l,-1,l3))}#-1 on the 5th variable is a sign of not convergence
    for(r in reset){
    pi=rep(1/J,J)
    means[r,]=data[sample(1:n,1),]
    D[[r]]=diag(1.2,p)}
  }

    posterior=E_bis(data,means,Q,D,pi)#matrix Jxn
    denom=rowSums(posterior)
    if(sum(denom==0)){r=(which(denom==0))[1]
      denom[r]=sum(denom)/100
      denom[-r]=denom[-r]*0.99}
    pi=denom/n
    #print(pi)
    #print(pi_ad)
    #means_wave=posterior%*%data
    #   means=t(sapply(1:J,function(k) means_wave[k,]/denom[k]))
    cov=list()
    for(j in 1:J)cov[[j]]=Q[[j]]%*%D[[j]]%*%t(Q[[j]])
    m_cov=cov[[1]]
    for(j in 2:J)m_cov=diag(pmax(diag(m_cov),diag(cov[[j]])))
    means_wave=posterior%*%data
    means_wave=t(sapply(1:J,function(k) means_wave[k,]/denom[k]))
    if(!is.na(add[1]))m=rbind(means,add) else m=means
    norm_mw=sqrt(colSums(apply(m,2,function(c) c^2)))#vector of norms of hat mu
    means_temp=matrix(0,nrow=J,ncol=dim(data)[2])#means to be
    limit=sqrt(apply((posterior%*%data)^2,2,sum))
    for(i in 1:p){
      for(j in 1:J)
        if(limit[i]>lambda_mean*m_cov[i,i]*pi_ad[j] & norm_mw[i]>0){
          means_temp[j,i]=means_wave[j,i]-lambda_mean*pi_ad[j]*m_cov[i,i]*means[j,i]/denom[j]/norm_mw[i]
        }#else if (iter==1)print(paste(i,j,limit[i],lambda_mean*m_cov[i,i],pi_ad[j],norm_mw[i],iter))
    }
    Sigma=lapply(1:J, function(j) apply(sapply(1:dim(data)[1],function(i) posterior[j,i]*(data[i,]-means[j,])%*%t(data[i,]-means[j,])),1,sum)/(denom[j]))
    for(i in 1:J){
      E=tryCatch({eigen(matrix(Sigma[[i]],nrow = p))},error=function(e){return(list(best_res,worst_res,last_res,l,-2,l3))})
      Q[[i]]=E$vectors
      E$values[E$values<0.005]=0.005
      # D[[i]]=diag((-denom[i]+sqrt(denom[i]^2+8*lambda*denom[i]*E$values/E$values[1]))/4/lambda*E$values[1])
      diagonala=(-denom[i]+sqrt(denom[i]^2+8*lambda/diag(adapt[[i]])*denom[i]*E$values))/4/lambda*diag(adapt[[i]])
      diagonala[diagonala<0.005]=0.005
        #Normal shrinkage
      #Shrinkage to epsilon based on hypothesis test
      for(kol in 2:p){
        eps_1=mean(diagonala[kol:p])
        eps_2=diagonala[p]
        eps_3=diagonala[kol]
        if(eps_3/eps_1<1-sqrt(2/n)*qnorm(0.025/(kol)) | eps_2/eps_1>1+sqrt(2/n)*qnorm(0.025/(kol))){
          diagonala[kol:p]=eps_1
          D[[i]]=diag(diagonala)
          break
        }
      }
    }
    means=means_temp
    like=log_likelihood(data,means,Q,D,pi)
    if(iter>10 & like_min>like){like_min=like
    worst_res=list(means,Q,D,pi)}
    if(iter>10 & like_max<like){like_max=like
    best_res=list(means,Q,D,pi)}
    l[iter]=like
    l2[iter]=D[[1]][p-1,p-1]
    l3[iter]=D[[2]][p-1,p-1]
  }
#  print(pi)
  last_res=list(means,Q,D,pi)
  list(best_res,worst_res,last_res,l,lambda,lambda_mean)
}


like=function(data,mean,Sigma,cols){
  p=dim(data)[2]
  #cols are colums for uuninformative variables
  if(sum(length(cols)==p)){
    return(sapply(1:dim(data)[1],function(k)
      dmvnorm(data[k,cols],mean,Sigma)))}
  sigma1=Sigma[-cols,-cols]
  sigma2=Sigma[cols,cols]
  sigma12=Sigma[-cols,cols]
  if(length(cols)==(p-1))sigma12=t(as.matrix(sigma12))
  if(length(cols)==1)sigma12=as.matrix(sigma12)
  m1=mean[-cols]
  m2=mean[cols]
  m_correct=sigma12%*%solve(sigma2)
  s_correct=sigma1-sigma12%*%solve(sigma2)%*%t(sigma12)
  if(length(cols)>1){
    sapply(1:dim(data)[1],function(k)
      dmvnorm(data[k,cols],m2,sigma2)*dmvnorm(data[k,-cols],m1+m_correct%*%(data[k,cols]-m2),s_correct))
  } else {
    sapply(1:dim(data)[1],function(k)
      dnorm(data[k,cols],m2,sqrt(sigma2))*dmvnorm(data[k,-cols],m1+m_correct%*%(data[k,cols]-m2),s_correct))
  }
}

Common_cov=function(data,E){
  J=length(E[[3]][[3]])
  p=dim(data)[2]
  pi=E[[3]][[4]]
  S=list()
  for(j in 1:J)
    S[[j]]=E[[3]][[2]][[j]]%*%E[[3]][[3]][[j]]%*%t(E[[3]][[2]][[j]])

  Sigma_bar=pi[1]*S[[1]]
  for(j in 2:J)Sigma_bar=Sigma_bar+pi[j]*S[[j]]
  cols=(1:p)[apply(E[[3]][[1]],2,sum)==0]
  k=length(cols)
  posterior=sapply(1:length(pi),function(i) dmvnorm(data,E[[3]][[1]][i,],S[[i]]))
  res0=-2*sum(log(posterior%*%pi))
  clist=list()
  clist[[1]]=cols
  if(k>1){
    for (ki in round(length(cols)*3/4-0.1):(k-1))clist=c(clist,combn(cols,ki,simplify = F))
  }
  #  if(sum(clist[[1]]==1:p)==p)clist[[1]]=1
  res=vector("integer",length = length(clist))
  if(length(cols)>0) for(ki in 1:length(clist)){
    R=list()
    for(j in 1:J){
      R[[j]]=Sigma_bar
      R[[j]][-clist[[ki]],-clist[[ki]]]=S[[j]][-clist[[ki]],-clist[[ki]]]
    }
    #   print(clist[[ki]])
    posterior2=sapply(1:length(pi),function(i) like(data,E[[3]][[1]][i,],R[[i]],clist[[ki]]))
    res[ki]=-2*sum(log(posterior2%*%pi))-log(n)*(length(clist[[ki]])*(2*p-length(clist[[ki]])))
  }
  ki=which.min(res)
  if (length(cols)<1){clist[[ki]]=2*p+1
  res[1]=Inf}
  #Old fit after dim red.
  clist[[ki]]
}

#Computes the E step for signal detection
#matrix nxQ pf probabilities
E2_bis=function(data,mu,Q,D,pi,probs){
  d1=sapply(dim(mu)[1],function(k) dmvnorm(data,mu[k,],Q[[k]]%*%D[[k]]%*%t(Q[[k]]))*pi[k])
  d1=cbind(probs,d1)#nx(q+1)matrix
  for(k in 1:length(Q))d1[,k]=d1[,k]*pi[k]
  d=apply(d1,1,function(x) x/sum(x))
  as.matrix(d)
}

EM2_bis=function(data,J,n_iter,lambda,lambda_mean,means,Qm,Dm,pi,adapt,add=NA,probs){
  #Fit a single signal component to J-1 background components
  p=dim(data)[2]
  n=dim(data)[1]
  l=1:n_iter
  #  km=kmeans(data,J)
 # means=rbind(means,mean+attributes(data)$`scaled:center`)
  means=rbind(means,add)
    Q=Qm
    D=Dm
    #D[[J]]=diag(1,p)
    #Q[[J]]=diag(1,p)
  #  means=matrix(rep(0,dim(data)[2]),nrow=1)
  #pi=c(0.95*pi,0.05)
  pi=pi
  #  pi_ad=pi
  like=-Inf
    best_res=list(means,Q,D,pi,like)
  for(iter in 1:n_iter){
    if(sum(pi<0.01)>0){reset=(which(pi<0.01))
    print("reset")
    for(r in reset){
      pi[3]=0.1
      pi[1:2]=pi[1:2]*(1-pi[3])/sum(pi[1:2])
      sam=sample(1:n,1)
    means[r,]=data[sam,]
   # print(sam)
    D[[r]]=diag(1.2,p)
    }
    }
    posterior=E2_bis(data,means,Q,D,pi,probs = probs)#matrix Jxn
        denom=rowSums(posterior)#Q+1 vector
      if(sum(is.na(denom)==TRUE)){print("error")
      print(denom)}
    if(sum(denom==0))denom=c(denom[-J]*0.99,sum(denom[-J])/100)
    pi[3]=denom[3]/dim(data)[1]
    pi[1:2]=pi[1:2]*(1-pi[3])/sum(pi[1:2])
    #    means_wave=posterior[J,]%*%data/denom[J]
    means_wave=posterior%*%data
    means_wave=t(sapply(1:J,function(k) means_wave[k,]/denom[k]))
    means[J,]=means_wave[J,]
    cov=Q[[J]]%*%D[[J]]%*%t(Q[[J]])
    norm_mw=sqrt(colSums(apply(means,2,function(c) c^2)))#vector of norms of hat mu
    if(is.na(norm_mw)[1]){print(means)
      print(pi)
      print(denom)
      print(iter)}
    means_temp=matrix(0,nrow=J,ncol=dim(data)[2])#means to be
    limit=sqrt(apply((posterior%*%data)^2,2,sum))
    for(i in 1:p){
      j=J
        if(limit[i]>lambda_mean*cov[i,i] & norm_mw[i]>0){
          means_temp[j,i]=means_wave[j,i]-lambda_mean*cov[i,i]*means[j,i]/denom[j]/norm_mw[i]
          if(sign(means_temp[j,i])!=sign(means_wave[j,i]))means_temp[j,i]=0
        }
    }
    Sigma=apply(sapply(1:dim(data)[1],function(i) posterior[J,i]*(data[i,]-means[J,])%*%t(data[i,]-means[J,])),1,sum)/(denom[J])
    i=J
      E=tryCatch({eigen(matrix(Sigma,nrow = p))},error=function(e){print("error")
        return(list(list(means,Q,D,pi,like,l),-1))})
      Q[[i]]=E$vectors
      E$values[E$values<0]=0.0005
     # shrunk=p+1-sum(diag(D[[i]])==D[[i]][p,p])
#      if(sum(is.na(sqrt(denom[i]^2+8*lambda/adapt*denom[i]*E$values)))){print(denom[i])
#        print(E$values)
#        print(denom[i]^2+8*lambda/adapt*denom[i]*E$values)}
      diagonala=(-denom[i]+sqrt(denom[i]^2+8*lambda/adapt*denom[i]*E$values))/4/lambda*adapt
      diagonala[diagonala<0.05]=0.05
      for(kol in 2:p){
        eps_1=mean(diagonala[kol:p])
        eps_2=diagonala[p]
        eps_3=diagonala[kol]
        if(eps_3/eps_1<1-sqrt(2/n)*qnorm(0.025/(kol)) | eps_2/eps_1>1+sqrt(2/n)*qnorm(0.025/(kol))){
          #kol=min(kol,shrunk)
          diagonala[kol:p]=eps_1
          #diagonala[kol:p]=mean(diagonala[kol:p])
          D[[i]]=diag(diagonala)
          break
        }
      }

    means[J,]=means_temp[J,]
    like2=log_likelihood(data,means,Q,D,pi)
    l[iter]=like2

#print(l[iter])
#       print(means[,8:14])
  }
    #print(pi)
    list(means,Q,D,pi,l)
}

#Compute the BIC given the reuslt of EM_ function
bic_mine=function(H,n){
  p=dim(H[[3]][[3]][[1]])[1]
  J=length(H[[3]][[4]])
  npar=p*J+J-1+J*(p+1)*p/2-sum(H[[3]][[1]]==0)
  for(j in 1:J){
    s=sum(diag(H[[3]][[3]][[j]])==H[[3]][[3]][[j]][p,p])
    npar=npar-(s-1)*s/2-s+1
  }
  -2*H[[4]][length(H[[4]])]+log(2*n)*npar # the minimum should be found
}

EM_BIC_mean=function(d,J,nit,lamev,lambda_mean,EM_best){
  bic_min=bic_mine(EM_best,n)
  pv=matrix(0,nrow=2,ncol=length(lambda_mean))
  for(lam in 1:length(lambda_mean)){
    Fd=EM_bis(d,J,nit,lamev,lambda_mean = lambda_mean[lam])
    pv[1,lam]=bic_mine(Fd,n)
    if(pv[1,lam]<bic_min){
      EM_best=Fd
      bic_min=bic_mine(EM_best,n)
    }
    pv[2,lam]=sum(Fd[[3]][[1]]==0)
  }
  list(pv,EM_best)
}

#Computes the probabilities of bck model on new data
prob_m1=function(data,mu,Q,D,pi){
  d1=sapply(1:dim(mu)[1],function(k) dmvnorm(data,mu[k,],Q[[k]]%*%D[[k]]%*%t(Q[[k]]))*pi[k])
  for(k in 1:length(Q))d1[,k]=d1[,k]*pi[k]
  as.matrix(d1)
}

EM_BIC_eigen=function(d,J,nit,lamev,lambda_mean){
  pv=matrix(0,nrow=2,ncol=length(lamev))
  for(lam in 1:length(lamev)){
    Fd=EM_adaptive(d,J,nit,lamev[lam],lambda_mean = lambda_mean,adapt = E[[3]][[3]],mean=E[[3]][[1]],Qm=E[[3]][[2]],Dm=E[[3]][[3]])

    pv[1,lam]=bic_mine(Fd,n)
    pv[2,lam]=sum(diag(Fd[[3]][[3]][[2]])==Fd[[3]][[3]][[2]][11,11])
  }
  pv
}
