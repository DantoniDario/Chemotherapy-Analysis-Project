Omega.matrix = function(model) {
  if(class(model)[1]=="lmerModLmerTest"|class(model)[1]=="lmerMod"){
    sigma2.epsilon = sigma(model)^2
    Psi.star       = crossprod(getME(model,"Lambdat"))*sigma2.epsilon
    Z              = getME(model,"Z")
    Omega          = Z %*% Psi.star %*% t(Z) + sigma2.epsilon* Diagonal(nrow(Z))
    Omega
  }else{
    warning("Function only works on outcput of function lmer()\n")
  }
}

Omega.plot = function(Omega,legend=TRUE,axes=TRUE){
  corw = cov2cor(Omega)
  if(any(corw<0)){
    colw=c(gray(1),rainbow(197)[197:100],gray(.9),rainbow(197)[99:1],gray(0))
    zlim=c(-1,1)
  }else{
    colw=c(gray(.9),rainbow(98)[98:1],gray(0))
    zlim=c(0,1)
  }
  image(z=as.matrix(corw[nrow(corw):1,]),zlim=zlim,axes=FALSE,col=colw)
  if(legend){    
    valw = as.numeric(names(table(as.matrix(corw))))  
    posw = round(valw*length(colw))
    posw[posw==0] = 1
    posw[posw>length(colw)] = length(colw)    
    legend("topright",ncol=1,legend=format(round(valw,4)),
           col=colw[posw],pch=15,bg="light gray",
           title="Values",box.lwd=NA)
  }
  if(axes){
    axis(3,at=seq(0,1,length=nrow(Omega)),labels=FALSE)
    axis(2,at=seq(0,1,length=nrow(Omega)),labels=FALSE)
    axis(2,at=c(1,0),c(1,nrow(Omega)),las=2)
    axis(3,at=c(0,1),c(1,nrow(Omega)),las=1)
  }
}
mod9.ml$
mod5.ml.hetero = mod9.ml

Omega = Omega.matrix(mod5.ml)
Omega.plot(Omega,legend=TRUE,axes=TRUE)

mod5.ml.sum = summary(mod5.ml)
mod5.ml.sum$varcor
mod5.ml.sum$varcor$patient[(,"stddev")]

sigma_0 = attr(mod5.ml.sum$varcor$patient, "stddev")
attr(mod5.ml.sum$varcor$patient, "stddev")
sigma_res = mod5.ml.sum$sigma
sigma_0^2/(sigma_0^2 + sigma_res^2)

mod5.ml<-lmer(tumour~sensitivity+line + month + I(month^2) + 
                      (1|patient), data = df, REML = F)
mod5.ml.hetero = lme(fixed = tumour~sensitivity*line + month + 
                       I(month^2), random=~1|patient, data = df, weights=varExp(form=~sensitivity), method = "ML")
summary(mod5.ml.hetero)


