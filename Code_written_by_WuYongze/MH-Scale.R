#基于MH-MCMC的ALSCAL降维算法
library(MASS)
set.seed(43221)

#计算样本距离函数
deta = function(X,p){
  if (p==1){
    rows = length(X)
  }
  else{
    rows = dim(X)[1]
  }
  deta = matrix(0,rows,rows)
  for (i in seq(1:rows)){
    for (j in seq(1:rows)){
      if (p == 1){
        deta[i,j] = sqrt(sum((X[i]-X[j])^2))
      }
      else{
        deta[i,j] = sqrt(sum((X[i,]-X[j,])^2))
      }
    }
  }
#  print(dim(deta))
  return(deta)
}

#计算分布中的SSR统计量
SSR_old = function(D,deta){
  SSR = 0
  rows = dim(D)[1]
  for (i in seq(1:rows)){
    for (j in seq(1:rows)){
      if (i<=j){
        D[i,j] = 0
        deta[i,j] = 0
      }
    }
  }
  SSR = sum((D-deta)^2)
  return(SSR)
}
#计算分布中的SSR统计量
SSR_new = function(n,sigma2){
  SSR = 0
  D = array(0,dim = c(n,n))
  for (i in seq(1:n)){
    for (j in seq(1:n)){
      D[i,j] = rnorm(1,0,sqrt(sigma2))
    }
  }
  for (i in seq(1:n)){
    for (j in seq(1:n)){
      if (i<=j){
        D[i,j] = 0
      }
    }
  }
  SSR = sum((D)^2)
  return(SSR)
}

#计算后验联合概率密度函数
log_P = function(Y,X,sigma2,lamda,a0,b0,arph0,beta0,p){
  #sigma2的IG先验的参数
  a = a0 
  b = b0
  #lamda的IG先验的参数
  arph = arph0
  beta = beta0 #长度为p的向量
  
  n = dim(X)[1]
  m = n*(n-1)/2
  p_11 = (-(m/2)-a-1)*log(sigma2)
  p_12 = sum((-(n/2))*log(lamda))#p的第一部分
  p_13 = sum((-(arph+1))*log(lamda))
  p_1 = p_11+p_12+p_13
  
  D_ij = deta(Y,dim(Y)[2])
  deta_ij = deta(X,p)
#  SSR = SSR_old(D_ij,deta_ij)
  SSR = SSR_new(n,sigma2)
  p_21 = -1/(2*sigma2)*SSR
  
  log_norm = log(pnorm(deta_ij/sqrt(sigma2),0,1))
  rows = dim(log_norm)[1]
  for (i in seq(1:rows)){
    for (j in seq(1:rows)){
      if (i<=j){
        log_norm[i,j] = 0
      }
    }
  }
  p_22 = -sum(log_norm)
  
  p_23 = 0
  X_new = X^2
  rows_Xnew = dim(X_new)[1]
  cols_Xnew = dim(X_new)[2]
  for (i in seq(1:rows_Xnew)){
    for (j in seq(1:cols_Xnew)){
      X_new[i,j] = X_new[i,j]/lamda[j]
      p_23 = p_23 + X_new[i,j]
    }
  }
  p_23 = -p_23/2
  
  p_24 = -b/sigma2
  
  p_25 = -sum(beta/lamda)
  
  p_2 = p_21+p_22+p_23+p_24+p_25
  p_value = p_1+p_2
  #print(p_value)
  return(p_value)
}

#Y = t(matrix(c(1,2,3,2,3,4,4,5,6,5,5,5,6,6,6,9,8,7),c(6,3)))#原始样本
BMDS = function(Y,p,niter){
  Y = as.matrix(Y)
#  p = 2 #降低到2维
  p = p
  n = dim(Y)[1] #样本数n
  m = n*(n-1)/2
#  niter = 100 #马氏链长度
  
  #构建参数马氏链
  X_chain = array(0,dim = c(n,p,niter)) #样本数n=3,降到维度p=2,迭代1000次
  sigma2_chain = array(0,niter) #sigma平方的链
  lamda_chain = array(0,dim = c(p,niter))#lamda的链 p行niter列
  
  #sigma2 = 5
  #设置先验参数
  a0 = 3
  #SSR = SSR_new(n,sigma2)
  #b0 = SSR/(n*(n-1)/2)*(a0-1)
  b0 = 4
  arph0 = 0.15
  beta0 = array(1,dim = c(p))
  
  #sigma2 = 1/rgamma(1,a0,b0)
  #定义提议分布的方差
  #sigma2_sigma2 = var(1/rgamma(100,m/2+a0,SSR/2+b0)) #sigma平方可以直接通过逆Gamma分布的公式得到
  sigma2_sigma2 = 0.1
  #Xi_sigma2 = diag(p)*(sigma2/(n-1)) #Xi的提议分布N的协方差矩阵(p*p 元素为1)
  Xi_sigma2 = diag(p)
  
  #设置后验初始值
  lamda = array(1,dim = c(p))
  #sigma2 = 1/rgamma(1,a0,b0)
  sigma2 = 0.1
  X = array(0,dim = c(n,p))
  for (i in 1:n){
    X[i,] = mvrnorm(n=1, rep(X[i,]), Xi_sigma2)
  }
  SSR = 0
  #for (j in 1:p){
  #  beta0[j] = var(X[,j])
  #}
  
  rate_X = 0
  rate_sigma2 = 0
  
  #MH方法迭代MCMC抽样
  for (iter in 1:niter){
    lamda_a = arph0+n/2
    lamda_b = array(0,dim = c(p))
    for (j in 1:p){
      beta0[j] = var(X[,j])
      lamda_b[j] = beta0[j]+var(X[,j])*n/2
      lamda[j] = 1/rgamma(1,lamda_a,lamda_b[j])
    }
    lamda_chain[,iter] = lamda
    
    X_new =  array(1,dim = c(n,p))
    for (i in 1:n){
      X_new[i,] = mvrnorm(n=1, rep(X[i,]), Xi_sigma2)
    }
    log_X = log_P(Y,X_new,sigma2,lamda,a0,b0,arph0,beta0,p)-log_P(Y,X,sigma2,lamda,a0,b0,arph0,beta0,p)
    if(runif(1) < exp(log_X) & !is.na(log_X)){
      X = X_new
      rate_X = rate_X+1
    }
    X_chain[,,iter] = X
    
    
    sigma2_new = rnorm(n=1, sigma2, sqrt(sigma2_sigma2))
    while(sigma2_new<0){
      sigma2_new = rnorm(n=1, sigma2, sqrt(sigma2_sigma2))
    }
    log_sigma2 = log_P(Y,X,sigma2_new,lamda,a0,b0,arph0,beta0,p)-log_P(Y,X,sigma2,lamda,a0,b0,arph0,beta0,p)
    if(runif(1) < exp(log_sigma2) & !is.na(log_sigma2)){
      sigma2 = sigma2_new
      rate_sigma2 = rate_sigma2+1
    }
    sigma2_chain[iter] = sigma2
    
    Xi_sigma2 = diag(p)*(sigma2/(n-1)) #更新X采样的sigma2
    deta_ij = deta(X,p)
    D_ij = deta(Y,dim(Y)[2])
    SSR = SSR_new(n,sigma2)
  #  SSR = SSR_old(D_ij,deta_ij)
    b0 = SSR/(n*(n-1)/2)*(a0-1) #更新b0参数
    sigma2_sigma2 = var(1/rgamma(100,m/2+a0,SSR/2+b0))#更新sigma2的方差
  }
  return(list(X_chain,sigma2_chain,lamda_chain,a0,b0,arph0,beta0,SSR,rate_X,rate_sigma2))
}

#MDSIC代码
log_LRp = function(SSR1,SSR2,X1,X2){
  n = dim(X2)[1]
  m = n*(n-1)/2
  p = dim(X2)[2]-1
  r_1 = (m-2)*(log(SSR2)-log(SSR1))
  
  s_1 = array(0,p)
  s_2 = array(0,p+1)
  for (j in 1:p){
    if (p == 1){
      s_1[j] = sum(X1^2)
    }
    else{
      s_1[j] = sum(X1[,j]^2)
    }
  }
  for (j in 1:(p+1)){
    s_2[j] = sum(X2[,j]^2)
  }
  r = array(0,p)
  for (j in 1:p){
    r[j] = s_2[j]/s_1[j]
  }
  r_21 = log(r*(n+1))-log(r+n)
  r_2 = (n+1)*sum(r_21)
  
  r_3 = (n+1)*log(n+1)
  
  LRP = r_1+r_2+r_3
  return(LRP)
}

MDSIC = function(Y,niter,p){
  X = BMDS(Y,1,niter)[[1]]
  deta0 = deta(X[,,niter],1)
  D = deta(Y,dim(Y)[2])
  SSR1 = SSR_old(D,deta0)
  n = dim(Y)[1]
  m = n*(n-1)/2
  MDSIC1 = (m-2)*log(SSR1)
  MDSICp = MDSIC1
  if (p == 1){
    MDSICp = MDSIC1
  }else{
    for (i in 1:(p-1)){
   #   SSR1 = BMDS(Y,i,niter)[[8]]
   #   print(SSR1)
   #   SSR2 = BMDS(Y,i+1,niter)[[8]]
   #   print(SSR2)
   #   print(i)
      print(i)
      X1_chain = BMDS(Y,i,niter)[[1]]
      X1 = X1_chain[,,niter]
      deta1 = deta(X1,i)
      SSR1 = SSR_old(D,deta1)
   #   print(X1)
      X2_chain = BMDS(Y,i+1,niter)[[1]]
      X2 = X2_chain[,,niter]
      deta2 = deta(X2,i+1)
      SSR2 = SSR_old(D,deta2)
   #   print(X2)
      MDSICp = MDSICp+log_LRp(SSR1,SSR2,X1,X2)
    }
  }
  return(MDSICp)
}

#算法使用部分
Y = read.csv("/Users/wuyongze/Desktop/统计计算实验数据.csv",header=T)
niter = 100
p = 1
out = BMDS(Y,p,niter)
X_chain = out[[1]]
#print(deta(X_chain[,,niter]))
MDSICp =MDSIC(Y,niter,p)
print(MDSICp)
#绘制图
plot(X_chain[2,1,],type = 'l')
#plot(sigma2_chain,type = 'l')
#plot(lamda_chain[1,],type = 'l')
#rate_X = rate_X/niter
#rate_sigma2 = rate_sigma2/niter
#print(rate_X)
#print(rate_sigma2)
print(X_chain[2,1,niter])
