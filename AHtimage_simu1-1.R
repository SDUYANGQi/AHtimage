library("corpcor")

K     = 3
s     = 2+K
re    = 1
n     = 500
nseed = 1000
dimage= 300
Arec  = array(NA,dim=c(nseed,n,s))
Srec  = array(NA,dim=c(nseed,n,s))
v     = function(u){View(cbind(u))}
na    = function(u){which(is.na(u))}
su    = function(u){summary(u)}
de    = function(u){dev.off()}
fm2   = function(x){x = array(x,dim = c(length(x),1)); return(x%*%t(x))}
fmat  = function(x){m = sqrt(length(x)); return(array(x,dim = c(m,m)))}

currdir = getwd()
newdir  = paste(currdir,"/result record" , sep = "")
dir.create(newdir)
cat(paste("The result is recorded in the new folder: \n", newdir, sep = ""))

for(re in 1:nseed){
  
  # data generation
  if(1){
    # true eigenvectors
    eivec = array(0,c(K,dimage*dimage))
    for(j in 1:K){
      atemp = array(0,dim = c(dimage,dimage))
      a1 = dimage/K; c1 = c(1:a1)
      if(j == 1){ atemp[c1,c1+2*a1]      = 1; eivec[j,] = c(atemp) }
      if(j == 2){ atemp[c1+1*a1,c1+1*a1] = 1; eivec[j,] = c(atemp) }
      if(j == 3){ atemp[c1+2*a1,c1]      = 1; eivec[j,] = c(atemp) }
    }
    if(1){
      for(j in 1:K){
        atemp = array(0,dim = c(dimage,dimage))
        a1 = dimage/K; c1 = c(1:a1)
        if(j == 1){ atemp[c1,c1] = 1; eivec[j,] = c(atemp) }
        if(j == 2){ atemp[c1+a1,c1+a1]     = 1; eivec[j,] = c(atemp) }
        if(j == 3){ atemp[c1+2*a1,c1+2*a1] = 1; eivec[j,] = c(atemp) }
      }
    }
    
    # define the coordinate framework for the image
    x_co = c(0,seq(dimage))
    y_co = c(0,seq(dimage))
    par(mfrow = c(1,K))
    
    # Plot the true eigenimages
    for(jj in 1:K){
      # This is to standarize the eigenimage
      eivec[jj,] = eivec[jj,]/sqrt(sum(eivec[jj,]^2))
    }
    
    # we generate the eigen_score for each subject based on the normal distribution, just follows the mixed effect representation
    eigen_score = array(0,c(n,K))
    
    # This is the variance for the noraml distributions generating the eigen_scores, mixed effects
    eigen_var = rep(0,K)
    eigen_var[1:K] = 0.5^((1:K) - 1)
    eigen_sd = sqrt(eigen_var)
    for(jj in 1:K){  eigen_score[,jj] = rnorm(n, mean = 0, sd = eigen_sd[jj]) }
    
    # generate the X
    true.funcs = eigen_score%*%eivec
    
    # demean the process, get the miu(t) out of the X
    mean.funcs = apply(true.funcs,2,mean)
    
    for(ii in 1:n){ true.funcs[ii,] = true.funcs[ii,] - mean.funcs }
    
    # fast svd to get eigenvalues and eigenvectors
    #library("corpcor")
    eigen_svd = fast.svd(t(true.funcs),tol = 0.0001)
    
    # estimated eigenimages
    eigenimage_est = t(eigen_svd$u)
    
    # whether to keep the direction of the eigenimages
    ind_rev = rep(0,K)
    for(jj in 1:K){
      # keep the direction be positive
      if(sum(eigenimage_est[jj,]*eivec[jj,])<0){
        eigenimage_est[jj,] = - eigenimage_est[jj,]
        ind_rev[jj] = 1
      }
    }
    xi = t(t(eigen_svd$v)*eigen_svd$d)
    #write.table(cbind(xi,eigen_score), file = "C:/Users/Administrator/Desktop/AHtim/xiei.txt")
    for(j in 1:K){xi[,j] = xi[,j]*(-2*ind_rev[j]+1)}
  }
  
  Z   = cbind(1, rbinom(n,3,.5), eigen_score)
  C   = runif(n,0,3)
  X   = array()
  Ti  = array(0); ran = array(0)
  A1 = function(t){t^2+1*t}; A2 = function(t){t^(1.5)/2};
  A3 = function(t){.5*t};    A4 = function(t){-.5*t^2};  A5 = function(t){.5*log(1+t)};
  
  for(i in 1:n){
    ran[i] = runif(1)
    ft     = function(t){ (A1(t) + A2(t)*Z[i,2] + A3(t)*Z[i,3] + A4(t)*Z[i,4] + A5(t)*Z[i,5])*1 + log(ran[i])}
    if(ft(5) > 0 ){ Ti[i] = uniroot(ft,c(0,5))$root}else{Ti[i] = 10}
  }
  su(Ti); sum(Ti == 10) # Ti[which(Ti <= 0)] = 10; # Ti[i]  = nleqslv(1,ft)$x
  
  for(i in 1:n){ X[i]  = min(Ti[i],C[i]) }
  delta   = as.numeric(Ti < C)
  cenrate = 1 - sum(delta)/n
  
  data = data.frame(cbind(Z[,c(1,2)],xi),X,delta)
  Z1   = data[,c(1:s)]; X = data[,s+1]; delta = data[,s+2]
  
  f = function(Z1,X,delta){
    
    n  = nrow(Z1); s = ncol(Z1)
    
    # sort order for X
    XZ  = as.matrix(cbind(X,Z1,delta))
    sor = XZ[order(XZ[,1]),]
    Xs  = sor[,1]
    Zs  = sor[,2:(s+1)]
    ds  = sor[,(s+2)]
    
    # Estimation of A(t) and its SEE 
    if(1){
      Hn    = array(0,dim = c(n,s,s))
      Hnpie = array(0,dim = c(n,s,s))
      for(k in 1:n){
        Hnpie[k,,] = Zs[k,]%*%t(Zs[k,])
      }
      for(i in n:1){
        if(i == n){
          Hn[i,,] = Hnpie[i,,]
        }else{
          Hn[i,,] = Hn[i+1,,] + Hnpie[i,,]
        }
      }
      
      Ahat = array(0,dim = c(n,s))
      Apie = array(0,dim = c(n,s))
      for(i in 1:(n-1)){
        temp = Hn[i,,]
        if(ds[i] == 1 & det(temp)>1e-8){
          Apie[i,] = solve(temp)%*%Zs[i,]
        }
      }
      for(k in 1:(n-1)){
        if(k == 1){
          Ahat[k,] = Apie[k,]
        } else{
          Ahat[k,] = Apie[k,] + Ahat[k-1,]
        }
      }
      # end of para esti
      
      # SEE estimation
      Gtpie = array(0,dim = c(n,s,s))
      Gt    = array(0,dim = c(n,s,s))
      for(i in 1:(n-1)){
        Gtpie[i,,] = Apie[i,]%*%t(Apie[i,])
        if(i == 1){
          Gt[i,,]  = Gtpie[i,,]
        }else{
          Gt[i,,]  = Gtpie[i,,] + Gt[i-1,,]
        }
      }
      SEEa = array(0,dim = c(n,s))
      for(i in 1:n){  SEEa[i,] = sqrt(diag(Gt[i,,]))  }
    }
    
    # 2. Smooth
    if(1){
      Xmax = 1
      int  = seq(0, Xmax,by = Xmax/(n-1))
      catx = array(1,n); lambda = array(1, n)
      for(i in 2:n){
        catx[i]   = min(which(sort(c(Xs,int[i]))==int[i]),n)
        if(catx[i] == 1){
          lambda[i] = 1
        }else{
          lambda[i] = (int[i] - Xs[catx[i]-1])/(Xs[catx[i]]-Xs[catx[i]-1])
        }
      }
      Ahats  = array(0,dim = c(n,s))
      SEEas  = array(0,dim = c(n,s))
      for(i in 2:(n-1)){
        catx[i] = max(2,catx[i])
        for(j in 1:s){
          Ahats[i,j] = Ahat[catx[i]-1,j] + lambda[i]*(Ahat[catx[i],j]-Ahat[catx[i]-1,j])
          SEEas[i,j] = SEEa[catx[i]-1,j] + lambda[i]*(SEEa[catx[i],j]-SEEa[catx[i]-1,j])
        }
      }
    }
    
    return(list(Ahats = Ahats, SEEas = SEEas, int = int))
  }
  
  TEMP = f(Z1,X,delta)
  Arec[re,,] = floor(TEMP$Ahats*1000)/1000
  Srec[re,,] = floor(TEMP$SEEas*1000)/1000
  int = TEMP$int
  #file = paste('/lustre/users/s1155099185/trans/mix-', num, ".txt", sep = "")
  if(re %% 10 == 0){print(paste("re =",re," & ", Sys.time()))}
  write.table(cbind(Arec[re,,],Srec[re,,],int), 
              file = paste(newdir,"/n1c1-", re, ".txt" ,sep = "") )
  
}

cat(paste("The result is recorded in the new folder: \n", newdir, sep = ""))
cat(paste("Please Use read_simu1-1.R for further result analysis"))

