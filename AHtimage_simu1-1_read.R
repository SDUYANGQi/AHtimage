# Please note that this is for all four cases in simu1-1.R
# You need to run four settings: n = 500/1000 and CR = 0.15/0.3 firstly.
# Please read results from the correct path !!!!!!!!!

path  = c("C:/Users/YANG QI/Desktop/AHtim/")
repe  = 1000
s     = 5
name  = c("n1c1-","n1c2-","n2c1-","n2c2-")
n     = c(500,500,1000,1000)

# Coverage probability funciton
cpf = function(res){
  p    = (ncol(res)-1)/2
  true = res[1:p,2*p+1]
  repe = nrow(res)
  est  = colMeans(res)
  sdd  = apply(res[,c(1:p)],2,sd)
  CP   = array(0,p)
  for(j in 1:p){
    for(i in 1:repe){
      if(abs(res[i,j] - true[j]) < 1.96*res[i,j+p]){
        CP[j] = CP[j] + 1/repe
      }
    }
  }
  mat = cbind(est[1:p]-true,est[(1:p)+p],sdd,CP)
  mat = floor(mat*1000)/1000
  colnames(mat) = c("Bias", "SE","SD", "CP")
  #rownames(mat) = paste("gamma",1:p, sep = "")#;mat
  return(mat)
}

fim = function(n,s){
  Aest = array(0, dim = c(n,s))
  Sest = array(0, dim = c(n,s))
  Sd   = array(0, dim = c(n,s))
  for(j in 1:s){
    Aest[,j] = colMeans(Arec[,,j])
    Sest[,j] = colMeans(Srec[,,j])
    for(i in 1:n){
      Sd[i,j] = sd(Arec[,i,j])
    }
  }
  
  mat = cbind(Aest,Sest,Sd)
  mat = floor(mat*1000)/1000
  #colnames(mat) = c("Bias", "SE","SD", "CP")
  #rownames(mat) = paste("gamma",1:p, sep = "");mat
  return(mat)
}

A1 = function(t){t^2+1*t}; A2 = function(t){t^(1.5)/2};
A3 = function(t){.5*t};    A4 = function(t){-.5*t^2};  A5 = function(t){.5*log(1+t)};

v     = function(u){View(cbind(u))}
su    = function(u){summary(u)}
de    = function(u){dev.off()}
fmat  = function(x){m = sqrt(length(x)); return(array(x,dim = c(m,m)))}

for(j in 1:4){
  Arec  = array(0,dim=c(repe,n[j],s))
  Srec  = array(0,dim=c(repe,n[j],s))
  for(re in 1:repe){
    temp = as.matrix(read.table(file = paste(path,name[j],re,".txt", sep = "")))
    Arec[re,,] = temp[,c(1:s)]
    Srec[re,,] = temp[,c(1:s)+s]
  }
  int = temp[,2*s+1]
  
  true = cbind(A1(int), A2(int), A3(int), A4(int), A5(int))
  # CP at t1,t2,t3;
  res1 = cbind(Arec[,n[j]*.1,],Srec[,n[j]*.1,],true[n[j]*.1,])
  res2 = cbind(Arec[,n[j]*.5,],Srec[,n[j]*.5,],true[n[j]*.5,])
  res3 = cbind(Arec[,n[j]*.9,],Srec[,n[j]*.9,],true[n[j]*.9,])
  
  temp  = rbind(cpf(res1),cpf(res2),cpf(res3))
  write.table(temp, file = paste(path,name[j],"CP.txt",sep = ""))
  
  temp1 = cbind(fim(n[j],s), true, int)
  write.table(temp1, file = paste(path,name[j],"image.txt",sep = ""))
  
  print(paste("j = ", j, " & ", Sys.time(), sep = ""))
}

# For bias, se, see, CP tables.
for(j in 1:4){
  temp = read.table(file = paste(path,name[j],"CP.txt",sep = ""))
  if(j == 1){cpfile = temp}else{cpfile = cbind(cpfile,temp)}
}
rownames(cpfile) = c(paste("A",1:5,"(t1)", sep = ""),
                     paste("A",1:5,"(t2)", sep = ""),paste("A",1:5,"(t3)", sep = ""))
write.table(cpfile, file = paste(path,"CPfile.txt",sep = ""),sep = " & ")

temp  = read.table(file = paste(path,"CPfile.txt",sep = ""),sep = "&")
colnames(temp) = NULL
rbind(temp[,1:8],temp[,9:16])

# For A(t) curves plotting
lis  = list()
for(j in 1:4){ lis[[j]] = as.matrix(read.table(file = paste(path,name[j],"image.txt",sep = ""))) }

name2 = c("n = 500, CR = 15%,","n = 500, CR = 30%,","n = 1000, CR = 15%,","n = 1000, CR = 30%,")

de(); par(mfrow = c(4,5))
par(mar=c(4.5, 4, 2, 2))
for(k in 1:4){
  temp = lis[[k]]
  Aest = temp[,c(1:s)]
  Sest = temp[,c(1:s)+s]
  Sd   = temp[,c(1:s)+2*s]
  true = temp[,c(1:s)+3*s]
  int  = temp[,1+4*s]
  
  cc  = c(1:(n[k]*.99))
  na1 = c("A1(t)","A2(t)","A3(t)","A4(t)","A5(t)")
  low = rep(-1,s)#c(-.5,-.5,rep(-.6,3)*1.5)
  upp = c(2.5,rep(1,s-1))#c(2.5,1,rep(.6,3)*1.5)
  for(j in 1:s){
    #if(j == s){ Aest[,j] = Aest[,j] + 0.62 }
    plot(int[cc],Aest[cc,j], type = "l", ylim = c(low[j],upp[j]),xlab = "Observed Time", 
         ylab = "Cumulative Coefficient", main = paste(name2[k]," A",j, "(t)", sep = ""))
    lines(int[cc],true[cc,j], col  = "red")
    lines(int[cc],Aest[cc,j] + 1.96*Sest[cc,j], col = "blue")
    lines(int[cc],Aest[cc,j] - 1.96*Sest[cc,j], col = "blue")
    lines(int[cc],Aest[cc,j] + 1.96*Sd[cc,j], col = "green")
    lines(int[cc],Aest[cc,j] - 1.96*Sd[cc,j], col = "green")
    #if(j == 2){ abline(h = 0.5)}
  }
}

###  The following is for eigen-images, coefficients images.
if(1){
  library("corpcor")
  
  K     = 3
  s     = 2+K
  re    = 1
  n     = 500
  nseed = 1000
  dimage= 300
  vec   = seq(1000,10000)
  Arec  = array(NA,dim=c(nseed,n,s))
  Srec  = array(NA,dim=c(nseed,n,s))
  v     = function(u){View(cbind(u))}
  na    = function(u){which(is.na(u))}
  su    = function(u){summary(u)}
  de    = function(u){dev.off()}
  fm2   = function(x){x = array(x,dim = c(length(x),1)); return(x%*%t(x))}
  fmat  = function(x){m = sqrt(length(x)); return(array(x,dim = c(m,m)))}
  
  set.seed(vec[re])
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
      atemp = array(0,dim = c(dimage,dimage));
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
  for(j in 1:K){xi[,j] = xi[,j]*(-2*ind_rev[j]+1)}
}



# define the coordinate framework for the image
# Plot the true eigenimages
x_co = c(0,seq(dimage))
y_co = c(0,seq(dimage))
par(mfrow = c(3,K),mar=c(3,3,3,3),oma=c(3,3,3,3))
nam1  = c(expression(paste("Eigenimage ",psi[1](v))),
          expression(paste("Eigenimage ",psi[2](v))),
          expression(paste("Eigenimage ",psi[3](v))))
nam2  = c(expression(paste("Eigenimage ",hat(psi)[1](v))),
          expression(paste("Eigenimage ",hat(psi)[2](v))),
          expression(paste("Eigenimage ",hat(psi)[3](v))))

# true eigenimages
for(jj in 1:K){
  # This is to standarize the eigenimage
  eivec[jj,] = eivec[jj,]/sqrt(sum(eivec[jj,]^2))
  # draw the eigenimage
  if(1){
    a1_temp = matrix(eivec[jj,],ncol = dimage, byrow = F)
    image(x_co,y_co,a1_temp,axes = F, col = grey(seq(1,0,length = 256)),xlab = "",ylab = "")
    axis(1,at = seq(0,dimage,by = dimage/3))
    axis(2,at = seq(0,dimage,by = dimage/3))
    title(main = nam1[jj])
    box()
  }
}
# estimated eigenimages
for(jj in 1:K){
  # This is to standarize the eigenimage
  eivec[jj,] = eivec[jj,]/sqrt(sum(eivec[jj,]^2))
  # draw the eigenimage
  if(1){
    a1_temp = matrix(eivec[jj,],ncol = dimage, byrow = F)
    image(x_co,y_co,a1_temp,axes = F, col = grey(seq(1,0,length = 256)),xlab = "",ylab = "")
    axis(1,at = seq(0,dimage,by = dimage/3))
    axis(2,at = seq(0,dimage,by = dimage/3))
    title(main = nam2[jj])
    box()
  }
}

## verify the estimated eigenscore and the real
#par( mfrow = c(1,3) )
nam  = c(expression(paste("Eigenscores ",xi[i1],"--", hat(xi)[i1],sep = "")),
         expression(paste("Eigenscores ",xi[i2],"--", hat(xi)[i2],sep = "")),
         expression(paste("Eigenscores ",xi[i3],"--", hat(xi)[i3],sep = "")))
for(j in 1:3){
  plot(xi[,j],eigen_score[,j],xlab = "True eigenscores", ylab = "Estimated eigenscores", main = nam[j])
  abline(0,1,col = "Red")
}



# define the true gamma function
#a3 = function(t){0.5}; a4 = function(t){-t}; a5 = function(t){0.5/(1+t)}; 
a3 = function(t){0.5*t}; a4 = function(t){-t*(.5*t)}; a5 = function(t){0.5*log(1+t)}; 
c2 = c(.1,.5,.9);
cutp    = numeric(3)
for(j in 1:3){cutp[j] = min(which(int > c2[j]))}
Gamma   = array(0, dim = c(3,3))
Gammavt = array(0, dim = c(3,dimage^2))
Gammae  = cbind(Aest[cutp[1],3], Aest[cutp[2],4], Aest[cutp[3],5])
Gammaest = array(0, dim = c(3,dimage^2))
for(k in 1:3){
  Gamma[k,]   = c(a3(c2[k]), a4(c2[k]), a5(c2[k]))# + runif(3,-.005,.005)
  Gammavt[k,] = Gamma[k,1]*eivec[1,] + Gamma[k,2]*eivec[2,] + Gamma[k,3]*eivec[3,]
  Gammaest[k,] = Gammae[k,1]*eivec[1,] + Gammae[k,2]*eivec[2,] + Gammae[k,3]*eivec[3,]
}


par(mfrow = c(2,K),mar=c(3,3,3,3),oma=c(3,3,3,3))
# plot the true Gammas
if(1){
  tit = c(expression(paste(Gamma, "(v,",t[1],")",sep = "")),
          expression(paste(Gamma, "(v,",t[2],")",sep = "")),
          expression(paste(Gamma, "(v,",t[3],")",sep = "")))
  #par(mfrow = c(1,3))
  tim  = c(0.1,0.5,0.9)
  tim1 = c(expression(t[1]),expression(t[2]),expression(t[3]))
  for(k in 1:3){
    a_temp = fmat(Gammavt[k,])
    a_temp = (a_temp - min(a_temp))/(max(a_temp) - min(a_temp))
    
    image(x_co,y_co,a_temp,axes = F, col = grey(seq(1,0,length = 256)),
          xlab = "",ylab = "",add = F)
    axis(1,at = seq(0,dimage,by = dimage/3))
    axis(2,at = seq(0,dimage,by = dimage/3))
    title(main = tit[k])
    #paste(expression(Gamma),"(v,",tim1[k],")",sep = "")
    box()
  }
}


# plot the estimated Gammas
if(0){
  tit = c(expression(paste(hat(Gamma), "(v,",t[1],")",sep = "")),
          expression(paste(hat(Gamma), "(v,",t[2],")",sep = "")),
          expression(paste(hat(Gamma), "(v,",t[3],")",sep = "")))
  #par(mfrow = c(1,3))
  tim  = c(0.1,0.5,0.9)
  tim1 = c(expression(t[1]),expression(t[2]),expression(t[3]))
  for(k in 1:3){
    a_temp = fmat(Gammaest[k,])
    a_temp = (a_temp - min(a_temp))/(max(a_temp) - min(a_temp))
    
    image(x_co,y_co,a_temp,axes = F, col = grey(seq(1,0,length = 256)),
          xlab = "",ylab = "",add = F)
    axis(1,at = seq(0,dimage,by = dimage/3))
    axis(2,at = seq(0,dimage,by = dimage/3))
    title(main = tit[k])
    #paste(expression(Gamma),"(v,",tim1[k],")",sep = "")
    box()
  }
}




