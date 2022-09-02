# Please read results from the correct path !!!!!!!!!

path  = c("C:/Users/YANG QI/Desktop/AHtim-extra/result record/")
repe  = 1000
s     = 6
name  = c("extra-")
n     = c(600)

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
  #rownames(mat) = paste("beta",1:p, sep = "")#;mat
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
  #rownames(mat) = paste("beta",1:p, sep = "");mat
  return(mat)
}

A1 = function(t){t^2+1*t}; A2 = function(t){t^(1.5)/2}; A3 = function(t){t/2};
A4 = function(t){.5*t};    A5 = function(t){-.5*t^2};   A6 = function(t){.5*log(1+t)};

v     = function(u){View(cbind(u))}
su    = function(u){summary(u)}
de    = function(u){dev.off()}
fmat  = function(x){m = sqrt(length(x)); return(array(x,dim = c(m,m)))}


for(j in 1:1){
  Arec  = array(0,dim=c(repe,n[j],s))
  Srec  = array(0,dim=c(repe,n[j],s))
  for(re in 1:repe){
    temp = as.matrix(read.table(file = paste(path,name[j],re,".txt", sep = "")))
    Arec[re,,] = temp[,c(1:s)]
    Srec[re,,] = temp[,c(1:s)+s]
  }
  int = temp[,2*s+1]
  
  true = cbind(A1(int), A2(int), A3(int), A4(int), A5(int), A6(int))
  
  # CP at t1,t2,t3;
  res1 = cbind(Arec[,n[j]*.1,],Srec[,n[j]*.1,],true[n[j]*.1,])
  res2 = cbind(Arec[,n[j]*.3,],Srec[,n[j]*.3,],true[n[j]*.3,])
  res3 = cbind(Arec[,n[j]*.6,],Srec[,n[j]*.6,],true[n[j]*.6,])
  
  temp  = rbind(cpf(res1),cpf(res2),cpf(res3))
  write.table(temp, file = paste(path,name[j],"CP.txt",sep = ""))
  
  temp1 = cbind(fim(n[j],s), true, int)
  write.table(temp1, file = paste(path,name[j],"image.txt",sep = ""))
  
  print(paste("j = ", j, " & ", Sys.time(), sep = ""))
}
# temp;temp1

# For CP
# for(j in 1:1){
#   temp = read.table(file = paste(path,name[j],"CP.txt",sep = ""))
#   if(j == 1){cpfile = temp}else{cpfile = cbind(cpfile,temp)}
# }
# rownames(cpfile) = c(paste("A",1:5,"(t1)", sep = ""),
#                      paste("A",1:5,"(t2)", sep = ""),paste("A",1:5,"(t3)", sep = ""))
# write.table(cpfile, file = paste(path,"CPfile.txt",sep = ""),sep = " & ")
# 
# temp  = read.table(file = paste(path,"CPfile.txt",sep = ""),sep = "&")
# colnames(temp) = NULL
# rbind(temp[,1:8],temp[,9:16])


# For A(t) curves
lis  = list()
for(j in 1:1){ lis[[j]] = as.matrix(read.table(file = paste(path,name[j],"image.txt",sep = ""))) }

name2 = c("n = 600, CR = 60%,","n = 500, CR = 30%,","n = 1000, CR = 15%,","n = 1000, CR = 30%,")

de(); par(mfrow = c(2,3))
par(mar=c(4.5, 4, 2, 2)); k = 1
for(k in 1:1){
  temp = lis[[k]]
  Aest = temp[,c(1:s)]
  Sest = temp[,c(1:s)+s]
  Sd   = temp[,c(1:s)+2*s]
  true = temp[,c(1:s)+3*s]
  int  = temp[,1+4*s]
  
  cc  = c(1:(n[k]*.6))
  #na1 = c("B1(t)","B2(t)","B3(t)","B4(t)","B5(t)","B6(t)")
  low = rep(-1,s)#c(-.5,-.5,rep(-.6,3)*1.5)
  upp = c(2.5,rep(1,s))#c(2.5,1,rep(.6,3)*1.5)
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

