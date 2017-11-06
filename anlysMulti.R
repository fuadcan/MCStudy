# n<-20;Tm=50;frho=0.6;k=4; noCons<-T

nlyzA <- function(Tm,n,k,frho,noCons){
  
  nocStr   <- if(noCons){"noCons"} else {"withCons"}
  dir <- paste0("Output/",nocStr,"/multiClub/")
  filename <- paste0("Results_",n,"-",k,"-A-",nocStr,"Mlt/reportA-",Tm,"-",frho,".rda")
  
  rep <- get(load(paste0(dir,filename)))
  
  nlyz <- sum(rep)
  
  return(nlyz)
  
}


nlyzHF05 <- function(Tm,n,k,frho,noCons){
  
  nocStr   <- if(noCons){"noCons"} else {"withCons"}
  dir <- paste0("Output/",nocStr,"/multiClub/")
  filename <- paste0(dir,"Results_",n,"-",k,"-",nocStr,"Mlt_HF/reportHF05-",Tm,"-",frho,".rda")
  repHF<-get(load(filename))[,(!noCons)+1]
  
  nlyz <- sum(repHF)
  
  return(nlyz)
}

overallnlyz<- function(TmVec,n,k,noCons){
  
  nlyzA5<- lapply(c(.2,.6), function(frho) sapply(TmVec, function(Tm)  nlyzA(Tm,n,k,frho,1)))
  nlyzHF<- lapply(c(.2,.6), function(frho) sapply(TmVec, function(Tm)  nlyzHF05(Tm,n,k,frho,1)))
  
  nlyz <- cbind(unlist(nlyzA5),unlist(nlyzHF))
  
  colnames(nlyz)<- c("adf","HF")
  
  return(nlyz)  
  
  
}

overallRepMulti <- function(){
a1<- lapply(c(T,F), function(nc) lapply(c(2,3), function(k) overallnlyz(c(50,75,100),10,k,nc)))
a2<- lapply(c(T,F), function(nc) lapply(c(4,5), function(k) overallnlyz(c(50,75,100),20,k,nc)))
a3<- lapply(c(T,F), function(nc) lapply(c(5,6), function(k) overallnlyz(c(50,75,100),30,k,nc)))
a4<- lapply(c(T,F), function(nc) lapply(c(7,8), function(k) overallnlyz(c(50,75,100),40,k,nc)))

noConst   <- do.call(rbind,c(a1[[1]],a2[[1]],a3[[1]],a4[[1]]))
withConst <- do.call(rbind,c(a1[[1]],a2[[1]],a3[[1]],a4[[1]]))

rep <- cbind(noConst,withConst)
colnames(rep) <- c("adf_noInt","HF_noInt","adf_wInt","HF_wInt")

dtType  <- sapply(c("2-10","3-10","4-20","5-20","5-30","6-30","7-40","8-40"), function(kn) 
  sapply(c(.2,.6), function(frho) sapply(c(50,75,100), function(Tm) paste(kn,frho,Tm,sep = "-"))))

rownames(rep) <- c(dtType)

dir.create("Results",recursive = TRUE)
write.csv(rep, file="Results/multi_perfectScores.csv", sep = ";",row.names = T, col.names = T)


}