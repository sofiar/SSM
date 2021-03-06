#### Funciones auxiliares
library(splancs)
library(RcppArmadillo)

## de una trayectoria de moviemiento (x,y,t) da la obsevacion correspondiente a tiempo dt
observe <- function(x,y,t,dt){
  nts <- ceiling(max(t)/dt)
  st <- (0:(nts-1)) * dt
  sx <- numeric(nts)
  sy <- numeric(nts)
  for(i in 2:nts){
    tmp <- max(which(t <= st[i]))
    dx <- x[tmp+1] - x[tmp]
    dy <- y[tmp+1] - y[tmp]
    dto <- t[tmp+1] - t[tmp]
    sx[i] <- x[tmp] + ((st[i]-t[tmp]) * dx)/dto
    sy[i] <- y[tmp] + ((st[i]-t[tmp]) * dy)/dto
  }
  return(data.frame(sx,sy,st))
}


# funcion para estimar el angulo de giro y el largo del paso.
# Devuelve: largo de paso, giro, direccion, seno del angulo ,coseno del angulo
pathelements <- function(X,Y){
  n = length(X)
  adj = X[2:n]-X[1:n-1] 
  op = Y[2:n]-Y[1:n-1]
  step = (adj^2 + op^2)^0.5
  aca = which(step==0)
  step[aca]=0.0000001
  
  si<-sign(adj)
  si[si==0]<-1  #corrects for sign(0) == 0
  ang = si*(op<0)*pi+atan(adj/op)
  adif <- ang[2:length(ang)]-ang[1:(length(ang)-1)]
  
  ## corregimos para que quede entre -pi y pi
  adif[which(adif < -pi)] = adif[which(adif < -pi)] + 2*pi
  adif[which(adif > pi)] = adif[which(adif > pi)] - 2*pi  
  
  turns <- adif
  direction <- ang  
  co = adj/step # cosine
  si = op/step  # sine
  res = list(step,turns,direction,co,si)
  names(res) = c("steps","turns","direction","cosine","sine")
  return(res)
}


# feature: how squared displacement increases with time DCT

cdt<-function(x,y,t)
{R1=x^2+y^2
 tmp = fastLm(formula = R1 ~ t + 0)
 return(as.numeric(tmp$coefficient))
}


cdt2<-function(sl,t)
{
tmp = fastLm(formula =sl ~ t + 0)
return(as.numeric(tmp$coefficient))
}


si<-function(c,s,p,b)
{
  tmp = 2*(p*((1-c^2-s^2)/((1-c)^2+s^2)+b^2))^(-0.5)
  return(tmp)
}


nsd<-function(elemento)
{
  ee=na.omit(elemento)
  val=numeric(length(ee))
  val[1]=ee[1]
  for(i in 2:length(ee))
  {
    #val[i]=sum(abs(ee[1:i]))  
    val[i]=val[i-1]+abs(ee[i])
  }
  return(val)
}

distancia2<-function(a,b)
{
  salida<-(a[1]-b[1])^2+(a[2]-b[2])^2
return(salida)
}
it<-function(pasos,xx,yy)
{
  N=length(xx)
  x0=xx[1]
  y0=yy[1]
  xF=xx[N]
  yF=yy[N]
  
  de=distancia2(c(x0,y0),c(xF,yF))  

  salida<-sqrt(de)/sum(pasos)
  #salida<-sum(pasos)/cha(xx,yy)
  
  return(salida)
  
   }

rtn<-function(xx,yy,tt)
{
  li=c(xx[1],yy[1])
  cords=cbind(xx,yy)
  dist.or=apply(cords,1,distancia2,b=li)
  sd(dist.or)
  }

cha<-function(x,y){
  tmp <- which(is.finite(x))
  i <- chull(x[tmp],y[tmp])
  return(areapl(cbind(x[i],y[i])))
}


rr<-function(ov,sv)
{
  qq<-qqplot(sv,ov,plot.it = FALSE)
  r<-1-sum((qq$x-qq$y)^2)/sum((qq$y-mean(qq$y))^2)
  return(r)
}



ppp<-function(x,y,t)
{
  if(length(na.omit(x))>5)
  {
    #val=numeric(ceiling((length(na.omit(x))/5))-1)
    val2=numeric(ceiling((length(na.omit(x))/5))-1)
    #tt=numeric(ceiling((length(x)/5))-1)
    for (i in 1:ceiling(length(na.omit(x))/5)-1)
    {
      xx=x[((i-1)*5+1):(i*5)]
      yy=y[((i-1)*5+1):(i*5)]
     # val[i]=(min(xx)-max(xx))*(min(yy)-max(yy))
      val2[i]=cha(xx,yy)
      
      #  tt[i]=sum(t((i-1)*5+1):(i*5))
    }
    
    
    #val2=lm(formula=val~tt+0)$coefficient
  }
  else
  {
    #val=NA
    val2=NA
    
  }
  return(val2)
}

  
Summaries2<-function(dir,turn,steps,cos,sin,ox,oy,ot,dt)
{
  a=cdt(ox,oy,ot)
  vv=length(ot)
  b=fastLm(cos[2:vv]~cos[1:(vv-1)]+0)$coefficients
  ct=mean(cos)
  st=mean(sin)
  
  #c=acf(ps$direction,plot=FALSE)$acf[3]  
  bo=sd(steps)
  
  salida<-c(mean(steps),
              sd(turn),
              bo,
              sicpp(ct,st,mean(steps),bo),
              b[1],
              sum(steps)/ot[vv],
              a)

  return(salida)

}
  
  



