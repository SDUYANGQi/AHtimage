# Please read results from the correct path !!!!!!!!!
rm(list = ls())

path = "C:/Users/YANG QI/Desktop/code revision" #getwd()
setwd(path)
source("myfunc.R")
options(scipen = 100)
n    = 585
repe = 1000
cut  = c(0.05,seq(0.1,0.7,by = 0.1))
len  = length(cut)
pat1 = "D:/AHtsimu2/compare"
imagepath = paste(path,"/slice/", sep = "")

fread  = function(M,cut){
  len  = length(cut)
  res  = array(0,dim = c(repe, len, M))
  res1 = array(0,dim = c(repe, n, 6))
  eps  = 0.05
  for(re in 1:repe){
    temp = read.table(paste(pat1, "/M", M, "re", re, ".txt", sep = ""))
    atem = as.matrix(temp)
    Xs   = temp$Xs
    ds   = temp$ds
    Ahat = atem[,1:(M+2)]
    
    for(j in 1:len){
      mn   = min(which(Xs > (cut[j] - eps)))
      mx   = max(which(Xs < (cut[j] + eps)))
      res[re,j,]  = colMeans(Ahat[mn:mx,1:M+2])
    }
    res1[re,,] = atem[,c(1:2,1:2+(M+2),1:2+2*(M+2))]
  }
  return(list(Gamt = res, Betat = res1))
}

fmean  = function(gamt){
  temp = gamt[1,,]
  for(i in 2:repe){
    temp = temp + gamt[i,,]
  }
  return(temp/repe)
}


# Slice reading and vectorization
if(1){
  # Slice reading
  filename = list.files(imagepath)
  aa       = as.matrix(read.table(paste(imagepath,filename[1], sep = "")))
  sizerow  = nrow(aa)
  sizecol  = ncol(aa)
  size     = sizerow*sizecol
  n        = length(filename)
  slice    = array(0, dim = c(n,nrow(aa),ncol(aa)))
  for(i in 1:n){
    slice[i,,] = as.matrix(read.table(paste(imagepath, filename[i], sep = "")))
  }
  
  # Find the common region (delete the contour)
  fshape = function(mat){   return(abs(sign(mat)))  }
  comm   = array(0, dim = c(nrow(aa),ncol(aa)))
  for(i in 1:n){  comm = comm + fshape(slice[i,,])  }
  shape = sign(comm)
  
  # Vectorize a slice
  shapevec = c(shape)
  ind0 = which(shapevec == 0)
  ind1 = which(shapevec == 1)
  fvec = function(a){  return(c(a)[ind1])  }
  finv = function(a){
    temp = array(NA,size)
    temp[ind1] = a
    return(array(temp, dim = c(sizerow,sizecol)))
  }
  
  # Design matrix Xim
  V   = length(ind1)
  Xim = array(NA, dim = c(V,n))
  for(j in 1:n){  Xim[,j] = fvec(slice[j,,])/1000  }
}

# True gamma(v,t) and reconstructed gamma(v,t)
library("oro.nifti")
ima1  = readNIfTI(paste(path, "/gamma1.nii.gz", sep = ""), reorient = FALSE)
ima2  = readNIfTI(paste(path, "/gamma2.nii.gz", sep = ""), reorient = FALSE)
gamma1 = ima1@.Data[,65,]
gamma2 = ima2@.Data[,65,]
fgamma = function(cut){
  A3  = function(t){t^5}; A4 = function(t){t^(0.2)}
  gt  = array(0, dim = c(length(cut),nrow(gamma1),ncol(gamma1)) )
  for(j in 1:length(cut)){  gt[j,,] = A3(cut[j])*gamma1 + A4(cut[j])*gamma2 }
  return(gt)
}
gt = fgamma(cut)/10


# result read
if(1){
  # Read the result from M = 50 and output M50.txt. 
  eimage = read.table(file = paste(path, "/psi.txt", sep = ""))
  reslis = list(); l = 1
  for(M in c(50)){
    temp   = fread(M,cut = cut)
    gamhat = fmean(temp$Gamt)
    psi    = fmtx(eimage[,1:M])
    
    V     = 17337
    gamvt = array(0, dim = c(len,V))
    for(j in 1:len){ gamvt[j,] = rowSums(psi%*%diag(gamhat[j,])) }
    
    reslis[[l]] = gamvt
    l = l+1
    write.table(t(gamvt), file = paste(path, "/M",M,".txt", sep = ""))
    print(paste("M = ", M, ", ", Sys.time(), sep = ""))
  }
  
  
  # Final plotting the images under M = 50
  if(1){
    par(mfrow = c(5,4), mar = rep(2.5,4))
    gamread = read.table(file = paste(path, "/M",50,".txt", sep = ""))
    nam     = paste("t = ", cut, sep = "")
    rms     = array(NA, 5)
    a1 = finv(fvec(gamma1)/10)
    a2 = finv(fvec(gamma2)/10)
    a1[1:10,1:10]  = a2[1:10,1:10]  =   0.15
    a1[11:20,1:10] = a2[11:20,1:10] =  -0.15
    fplot1(a1, nam = expression(tilde(X)[1](v)))
    fplot1(a2, nam = expression(tilde(X)[2](v)))
    fplot1(a1, nam = expression(tilde(X)[1](v)))
    fplot1(a2, nam = expression(tilde(X)[2](v)))
    
    #par(mfrow = c(4,4))
    for(j in c(1:8)){
      esti   = gamread[,j]
      estim  = finv(esti)
      ti     = fvec(gt[j,,])
      tim    = finv(ti)
      rms[j] = sqrt(mean((esti - ti)^2))
      
      estim[1:10,1:10]  = tim[1:10,1:10]  =   0.15
      estim[11:20,1:10] = tim[11:20,1:10] =  -0.15
      fplot1(tim,  nam = nam[j])
      #fplot1(estim,nam = nam[j])
    }
    for(j in c(1:8)){
      esti   = gamread[,j]
      estim  = finv(esti)
      ti     = fvec(gt[j,,])
      tim    = finv(ti)
      rms[j] = sqrt(mean((esti - ti)^2))
      
      estim[1:10,1:10]  = tim[1:10,1:10]  =   0.15
      estim[11:20,1:10] = tim[11:20,1:10] =  -0.15
      #fplot1(tim,  nam = nam[j])
      fplot1(estim,nam = nam[j])
    }
    color.bar(hcl.colors(100, "YlOrRd", rev = TRUE),-0.15,0.15, nticks = 5)
  }
  
  # plot of the first eight eigenimages
  # par(mfrow = c(2,4))
  # for(j in 1:8){ 
  #   temp = finv(eimage[,j])
  #   temp[1:10,1:10]   =   0.05
  #   temp[11:20,1:10]  =  -0.05
  #   fplot1(temp, nam = paste(j,"th eigenimage", sep = ""))
  # }
  # color.bar(hcl.colors(100, "YlOrRd", rev = TRUE),-0.05,0.05, nticks = 5)
}






 