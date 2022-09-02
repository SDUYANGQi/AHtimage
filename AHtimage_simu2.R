rm(list = ls())

path = getwd()
library("oro.nifti")
source("myfunc.R")
imagepath = paste(path,"/slice/", sep = "")

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

# FPCA implementation (K is the retained numver of eigenimages/scores)
result = fpca(Xim, K = 100)

eigenscore = result$xi
write.table(result$xi,   file = paste(path, "/xi.txt", sep = ""))
write.table(result$psik, file = paste(path, "/psi.txt", sep = ""))

# True gamma(v) and data generation
ima1  = readNIfTI(paste(path, "/gamma1.nii.gz", sep = ""), reorient = FALSE)
ima2  = readNIfTI(paste(path, "/gamma2.nii.gz", sep = ""), reorient = FALSE)
gamma1 = ima1@.Data[,65,]
gamma2 = ima2@.Data[,65,]
g1    = fvec(gamma1)/10
g2    = fvec(gamma2)/10

xi1 = array(NA,n)
xi2 = array(NA,n)
for(i in 1:n){
  xi1[i] = sum(g1*Xim[,i])
  xi2[i] = sum(g2*Xim[,i])
}


repe = 1000; M = 30; re = 1
for(M in c(50)){
  for(re in 1:repe){
    # data generation
    Z   = cbind(1, rbinom(n,3,.5), xi1,xi2)
    C   = runif(n,0,2)
    X   = array()
    Ti  = array(0); ran = array(0)
    A1 = function(t){t^2}; A2 = function(t){0.5*t^(1.5)}
    A3 = function(t){t^2}; A4 = function(t){t^(0.5)}
    
    for(i in 1:n){
      ran[i] = runif(1)
      ft     = function(t){ A1(t) + A2(t)*Z[i,2] + A3(t)*Z[i,3] + A4(t)*Z[i,4] + log(ran[i]) } # A5(t)*Z[i,3]
      if(ft(2) > 0 ){ Ti[i] = uniroot(ft,c(0,2))$root}else{Ti[i] = 10}
    }
    # su(Ti); sum(Ti == 10)
    
    for(i in 1:n){ X[i]  = min(Ti[i],C[i]) }
    delta   = as.numeric(Ti < C)
    cenrate = 1 - sum(delta)/n
    
    Zr   = cbind(Z[,c(1,2)], as.matrix(eigenscore[,c(1:M)]))#; su(Zr)
    TEMP = faht(Zr,X,delta)
    Ahat = as.matrix(TEMP$Ahats)
    SEEa = as.matrix(TEMP$SEEas)
    Xs   = int = as.numeric(TEMP$int)
    ds   = as.numeric(TEMP$ds)
    
    if(re%%50 == 0) print(paste("M = ", M, ", re = ", re, " & ",Sys.time(), sep = ""))
    write.table(cbind(Ahat,SEEa,Xs,ds), file = paste(path, "/M", M, "re", re, ".txt", sep = ""))
  }
}

