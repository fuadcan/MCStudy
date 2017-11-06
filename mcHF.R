# Tm=75;n=10;clsize=5; frho=.6; ind=1
#########################


library("stringr")
mcHF <- function(Tm,n,clsize,frho,noCons){

  inDirHF  <- if(noCons){"Data/noCons/singleClub/"} else {"Data/withCons/singleClub/"}
  outDirHF <- if(noCons){"Output/noCons/singleClub/"} else {"Output/withCons/singleClub/"}
  
  dirname  <- inDirHF
  nocStr   <- if(noCons){"-noCons"} else {"-withCons"}
  
  cat("Initializing variables\n")
  
  filename <- paste0("zZ_",Tm,"-",n,"-",clsize,"-",frho,nocStr,".rda")
  zz       <- get(load(paste0(dirname,filename)))
  lenz<-length(zz[[1]])
  list   <- zz[[2]]
  gammas <- zz[[3]][1:n]
  
  #########################
  write.table(c(Tm,n),"hfcodes/dims.csv",row.names = FALSE,col.names = FALSE)
  
  repForDat <- function(ind){
    
    z      <- zz[[1]][[ind]]
    
    
  write.table(z,file = "hfcodes/datt.csv",row.names = FALSE,col.names = FALSE)
  ########################################################

  cat("Analyzing")
    
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
  
    if(length(absList)==0){absList=list(c(0,0,0),c(0,0,0))}
    if(length(relList)==0){relList=list(c(0,0,0),c(0,0,0))}
    
    ############################## Pre-evaluation ##############################
    
    indAbs    <- sapply(1:length(absList), function(x) length(intersect(absList[[x]],list)))
    maxIndAbs <- sample(which(max(indAbs)==indAbs),1)
    maxlAbs   <- absList[[maxIndAbs]]
    indRel    <- sapply(1:length(relList), function(x) length(intersect(relList[[x]],list)))
    maxIndRel <- sample(which(max(indRel)==indRel),1)
    maxlRel   <- relList[[maxIndRel]]
    
    maxHFabs <- length(intersect(maxlAbs,list)); excHFabs<- length(setdiff(maxlAbs,list))
    maxHFrel <- length(intersect(maxlRel,list)); excHFrel<- length(setdiff(maxlRel,list))
    repHF <- matrix(c(maxHFabs,excHFabs,maxHFrel,excHFrel),1,4)
    
    gmmlHF <- t(matrix(rep("",2*n),n))
    gmmlHF[1,][maxlAbs] <- "c";  gmmlHF[2,][maxlRel] <- "c"
    
    return(list(repHF,gmmlHF))
    ############################## END REPORT ##############################
  }
  
  
  
  ############################# CONSOLIDATION ######################
  
  consConcs <- lapply(1:lenz,function(x){cat(paste0(frho,"-",x,"\n")); repForDat(x)})
  cat("Consolidating\n")
  
  reportHF05 <- lapply(1:lenz, function(x) t(consConcs[[x]][[1]]))
  reportHF05 <- t(matrix(unlist(reportHF05),4))
  
  gmmlsHF05  <- lapply(1:lenz, function(x) t(consConcs[[x]][[2]]))
  gmmlsHF05  <- rbind(gammas,t(matrix(unlist(gmmlsHF05),n)))
  
  
  savedir <- outDirHF
  outdir  <- paste0(savedir,"Results_",n,"-",clsize,nocStr,"_HF/")
  filedir <- paste0(Tm,"-",frho,".rda")
  
  
  dir.create(outdir,recursive = TRUE)
  cat("Saving\n")
  
  save(reportHF05,file = paste0(outdir,"reportHF05-",filedir))
  save(gmmlsHF05,file = paste0(outdir,"gmmlHF05-",filedir))
  
}