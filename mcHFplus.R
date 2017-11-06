
# Tm=100;n=10;k=3; frho=.6; ind=1
library("stringr")

mcHFplus <- function(Tm,n,k,frho,noCons){
  
  inDirHFm  <- if(noCons){"Data/noCons/multiClub/"} else {"Data/withCons/multiClub/"}
  outDirHFm <- if(noCons){"Output/noCons/multiClub/"} else {"Output/withCons/multiClub/"}
  
  cat("initializing parameters\n")
  
  dirname  <- inDirHFm
  
  nocStr   <- if(noCons){"-noConsMlt"} else {"-withConsMlt"}
  
  filename <- paste0("zZ_",Tm,"-",n,"-",k,"-",frho,nocStr,".rda")
  zz       <- get(load(paste0(dirname,filename)))
  lenz<-length(zz[[1]])
  list   <- zz[[2]]
  gammas <- zz[[3]]
  
  fsts<- sapply(list, function(x) x[1])  
  list<- list[order(fsts)]
  
  
  write.table(c(Tm,n),"hfcodes/dims.csv",row.names = FALSE,col.names = FALSE)
  
  
  
  repForDat <- function(ind){
    
  z      <- zz[[1]][[ind]]
    
    
    
    write.table(z,file = "hfcodes/datt.csv",row.names = FALSE,col.names = FALSE)
    
    cat("Analyzing\n")
    
    ########################################################
    tempHF <- shell(paste("C:/gauss10/tgauss -b ",'RUN hfcodes\\main05.gss'), intern=TRUE,wait=TRUE)
    tempHF <- tempHF[1:(which(tempHF=="GAUSS 10.0.3 (Dec 22 2009, 1345) 32-bit")-2)]
    
    aCrude<- strsplit(tempHF[1:((which(tempHF=="brkpnt"))-1)]," ")
    aCrude<-lapply(1:length(aCrude), function(x) aCrude[[x]] <- aCrude[[x]][aCrude[[x]]!=""])
    aCrude<-lapply(1:length(aCrude), function(x) as.numeric(str_replace_all(aCrude[[x]],"c","")))
    absList<- aCrude[sapply(1:length(aCrude), function(x) length(aCrude[[x]])!=1)]
    
    rCrude<- strsplit(tempHF[((which(tempHF=="brkpnt"))+1):length(tempHF)]," ")
    rCrude<-lapply(1:length(rCrude), function(x) rCrude[[x]] <- rCrude[[x]][rCrude[[x]]!=""])
    rCrude<-lapply(1:length(rCrude), function(x) as.numeric(str_replace_all(rCrude[[x]],"c","")))
    relList<- rCrude[sapply(1:length(rCrude), function(x) length(rCrude[[x]])!=1)]
    
    fstsABS <- sapply(absList, function(x) x[1]); fstsREL <- sapply(relList, function(x) x[1])  
    absList  <- absList[order(fstsABS)]; relList  <- relList[order(fstsREL)]
    
    if(length(absList)==0){absList=list(c(0,0,0),c(0,0,0))}
    if(length(relList)==0){relList=list(c(0,0,0),c(0,0,0))}
    
  
    ############################## Evaluation ##############################
    cat("Pre-Evaluation\n")
    listCode <- unlist(lapply(list, function(x) c(x,0)))
    
    sccHFabs05<- suppressWarnings(mean(unlist(lapply(absList, function(x) c(sort(x),0)))==listCode)==1)
    sccHFrel05<- suppressWarnings(mean(unlist(lapply(relList, function(x) c(sort(x),0)))==listCode)==1)
    
    repHF05<- matrix(c(sccHFabs05,sccHFrel05),1)
    gmmlHF05 <- t(matrix(rep("",2*n),n))
    for(i in 1:length(absList)){gmmlHF05[1,][absList[[i]]] <- paste0("c",i)} 
    for(i in 1:length(relList)){gmmlHF05[2,][relList[[i]]] <- paste0("c",i)}
    
    
    return(list(repHF05,gmmlHF05))
    ############################## END REPORT ##############################
  }
  
  
  
  consConcs <- lapply(1:lenz,function(x){cat(paste0(frho,"-",x,"\n")); repForDat(x)})
  cat("Consolidating Results\n")
  
  reportHF05 <- lapply(1:lenz, function(x) t(consConcs[[x]][[1]]))
  reportHF05 <- t(matrix(unlist(reportHF05),,lenz))
  
  
  gmmlsHF05  <- lapply(1:lenz, function(x) t(consConcs[[x]][[2]]))
  gmmlsHF05  <- rbind(gammas,t(matrix(unlist(gmmlsHF05),n)))
  
  
  savedir <- outDirHFm
  outdir  <- paste0(savedir,"Results_",n,"-",k,nocStr,"_HF/")
  filedir <- paste0(Tm,"-",frho,".rda")
  
  
  dir.create(outdir,recursive = TRUE)
  cat("Saving\n")
  save(reportHF05,file = paste0(outdir,"reportHF05-",filedir))
  save(gmmlsHF05,file = paste0(outdir,"gmmlHF05-",filedir))
  
  cat("Finished\n")
  
}