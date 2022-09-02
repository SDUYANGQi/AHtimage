v     = function(u){View(cbind(u))}
na    = function(u){which(is.na(u))}
su    = function(u){summary(u)}
de    = function(u){dev.off()}
fm2   = function(x){x = array(x,dim = c(length(x),1)); return(x%*%t(x))}
fmat  = function(x){m = sqrt(length(x)); return(array(x,dim = c(m,m)))}

# Package loading
loadPkg  = function(pkg) {
  new.pkg = pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  suppressMessages(invisible(sapply(pkg, require, character.only = TRUE)))
}
#packages = c("SemiPar", "corpcor", "rje", "MASS","mvtnorm","nleqslv")
#loadPkg(packages)

color.bar = function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  dev.new(width=1.75, height=5)
  par(bg="white")
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title,col.axis="black")
  axis(2, ticks, las=1,col="black",col.axis = "black")
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}


fmm  = function(x){return(c(min(x),max(x)))}

fper = function(x){
  len = length(x)
  new = x
  sux = sum(x)
  stepx = array(0,len)
  for(i in 1:len){
    stepx[i] = ifelse(i == 1, x[i], x[i]+stepx[i-1])
    new[i]   = stepx[i]/sux
  }
  return(new)
}

fpca = function(Xim, K){
  
  # Centralization
  meanx = apply(Xim, 1, mean)
  Xim1  = Xim
  for(j in 1:n){Xim1[,j] = (Xim[,j] - meanx)}
  
  xx = t(Xim1)%*%Xim1
  
  # Produce Xim = UDV'
  library("corpcor")
  eigen_svd = fast.svd(Xim1,tol = 0.0001)
  
  BU = eigen_svd$u
  BD = eigen_svd$d
  BV = eigen_svd$v
  
  # eigenimages psi_k(v) are the kth column of BU
  psi = BU[,c(1:K)]
  # eigenscores \xi_ik are the kth column of DV'
  xi  = t(diag(BD)%*%t(BV))[,1:K]
  return(list(xi = xi, psik = psi))
}

fplot  = function(aa){
  ####define the coordinate framework for the image
  x_co = c(0,seq(nrow(aa)))
  y_co = c(0,seq(ncol(aa)))
  image(x_co,y_co,aa,axes = F, col = rainbow(150), bg= NULL,xlab = "",ylab = "")
  # axis(1,at = seq(0,nrow(aa), by = 10)) 
  # axis(2,at = seq(0,ncol(aa), by = 10))
  title(main = "axial")
}

fplot1 = function(aa,nam){
  ####define the coordinate framework for the image
  x_co = c(0,seq(nrow(aa)))
  y_co = c(0,seq(ncol(aa)))
  image(x_co,y_co,aa,axes = F, col = hcl.colors(100, "YlOrRd", rev = TRUE), bg= NULL,xlab = "",ylab = "")
  # axis(1,at = seq(0,nrow(aa), by = 10)) 
  # axis(2,at = seq(0,ncol(aa), by = 10))
  title(main = nam)
}

fmtx  = function(psi){
  newmat = array(0, dim = c(nrow(psi), ncol(psi)))
  for(j in 1:M){
    newmat[,j] = as.numeric(psi[,j])
  }
  return(newmat)
}

faht = function(Z,X,delta){
  
  n  = nrow(Z); s = ncol(Z)
  
  # sort order for X
  XZ  = as.matrix(cbind(X,Z,delta))
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
  if(0){
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
  
  return(list(Ahats = Ahat, SEEas = SEEa, int = Xs, ds = ds))
}

