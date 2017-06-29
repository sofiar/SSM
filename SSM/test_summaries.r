### busco la mejor combinacion de summaries para el parametro scale

Q=matrix(NA,12,14)
for(t in 2:12)
{
  ms=6
  while(ms>2)
  {
  A=matrix(0,2,2)
  while (abs(det(A))<0.1 )
  {
  q=sample(1:13,t)
  Ssim=s[,q[1]]
  Stobs=ss[q[1]] 
  
  for (j in 2:t)
  {
    Ssim=cbind(Ssim,s[,q[j]])
    Stobs=c(Stobs,ss[q[j]])
    
  }
  A=cov(na.omit(Ssim))
  
  }
  #XX=(as.matrix(Ssim[2,]-Stobs)%*%solve(A)%*%as.matrix(t(Ssim[2,]-Stobs)))
  Mahal=apply(Ssim,1,mahalanobis,center=Stobs,cov=A)
  
  which=numeric(nsims)
  which[order(Mahal)[1:nbest]]=1
  ms=mean(s_scale[which==1])
  }
  dd=numeric(13)
  dd[q]=1
  Q[t,]=c(ms,dd)
  print(t)
  print(ms)
}

### Ahora supongamos que tengo fijos dos summaries y quiero ver cual puedo agregar
# fijos 6,7
Ssim=cbind(s[,6],s[,7],s[,9])
Stobs=c(ss[6],ss[7],ss[9])
qq=c(1:5,8,10:13)
QQ=matrix(NA,10,2)


  for(j in 1:10)
  {
      q=qq[j]
      
      Ssim=cbind(s[,6],s[,7],s[,9])
      Stobs=c(ss[6],ss[7],ss[9])
      
      Ssim=cbind(Ssim,s[,q])
      Stobs=cbind(Stobs,ss[q]) 
      A=cov(na.omit(Ssim))
      
      if (abs(det(A))>0.000000001 )
      {

    #XX=(as.matrix(Ssim[2,]-Stobs)%*%solve(A)%*%as.matrix(t(Ssim[2,]-Stobs)))
    Mahal=apply(Ssim,1,mahalanobis,center=Stobs,cov=A)
    
    which=numeric(nsims)
    which[order(Mahal)[1:nbest]]=1
    ms=mean(s_scale[which==1])
      }
      else {
        ms=NA
           }
     QQ[j,]=cbind(q,ms)
    
    
    }
  


