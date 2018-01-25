#######
## multi.select: code to select combinations of predictors based on
##               AUCs for multiple levels of the outcome
#######

## function to estimate bootstrap-corrected AUCs
bootAUC <- function(dataset, Dname, Mnames){
  bootsamp <- dataset[sample(1:nrow(dataset), nrow(dataset), replace=TRUE), ]
  K <- length(unique(dataset[,1]))
  Kboot <- length(unique(bootsamp[,1]))
  if(K != Kboot){
    warning("bootstrap sample does not have data from all K levels of the outcome")
    return(c(NA,NA))
  }else{
    ### fit model in bootsamp
    modelB <- glm(as.formula(paste("as.numeric(",Dname,"==",K,") ~ ",
                                   paste(Mnames,collapse=" + "),sep="")),
                  data=bootsamp, family=binomial)
    ### get AUC in bootsamp
    appAUCv1 <- Hmisc::somers2(predict.glm(modelB), as.numeric(bootsamp[,1]==K))[1]
    appAUCv2 <- Hmisc::somers2(predict.glm(modelB)[which(bootsamp[,1] < K)],
                               as.numeric(bootsamp[which(bootsamp[,1] < K), 1]==K-1))[1]
    ### apply AUC to original data
    origAUCv1 <- Hmisc::somers2(predict.glm(modelB,newdata=dataset), as.numeric(dataset[,1]==K))[1]
    origAUCv2 <- Hmisc::somers2(predict.glm(modelB,newdata=dataset)[which(dataset[,1] < K)],
                                as.numeric(dataset[which(dataset[,1] < K), 1]==K-1))[1]
    ### estimate optimism
    return(c(appAUCv1-origAUCv1, appAUCv2-origAUCv2))
  }
}

## function to do selection
multiselect <- function(data, size=2, Breps=40, nummod=10){
  if(!is.data.frame(data)){
    stop("data must be a data.frame")
  }
  if(min(sapply(data, is.numeric)) != 1){
    stop("columns of data must be numeric")
  }
  if(max(sort(unique(data[,1])) != c(1:length(unique(data[,1])))) > 0){
    stop("outcome (first column of data) must include all integers between 1 and K")
  }
  if(max(sort(unique(data[,1]))) < 3){
    stop("outcome must have at least 3 levels (K >= 3)")
  }
  if(max(is.na(data))==1){
    stop("missing values in data not allowed")
  }

  p <- ncol(data)-1 ## 1st column is the outcome, other columns = markers
  if(is.null(names(data))){
    names(data) <- c("D",paste("V",c(1:p),sep=""))
  }
  rsltsmat <- t(combn(1:p, size))
  allrslts <- cbind(rsltsmat, matrix(NA, nrow=nrow(rsltsmat), ncol=size+3))
  for(i in 1:nrow(allrslts)){
    moddata <- data.frame(data[,1], data[,(allrslts[i,c(1:size)]+1)])
    names(moddata) <- c("D",paste("M",c(1:size),sep=""))
    Dname <- names(moddata)[1]
    Mnames <- names(moddata)[-1]
    K <- length(unique(moddata[,1]))
    model <- glm(as.formula(paste("as.numeric(",Dname,"==",K,") ~ ",
                                  paste(Mnames,collapse=" + "),sep="")),
                 data=moddata, family=binomial)
    trainAUCv1 <- Hmisc::somers2(predict.glm(model), as.numeric(moddata[,1]==K))[1]
    trainAUCv2 <- Hmisc::somers2(predict.glm(model)[which(moddata[,1] < K)],
                                 as.numeric(moddata[which(moddata[,1] < K), 1]==K-1))[1]
    bootrslt <- replicate(Breps, bootAUC(dataset=moddata, Dname=Dname, Mnames=Mnames))
    estopt <- rowMeans(bootrslt,na.rm=T)
    allrslts[i,(size+1):(size+2)] <- c(trainAUCv1,trainAUCv2)-estopt
    allrslts[i,(size+3)] <- sum(is.na(bootrslt[1,]) | is.na(bootrslt[2,]))
    allrslts[i,(size+4):(2*size+3)] <- model$coef[-1]
  }

  rankv1 <- rank(allrslts[,size+1])
  rankv2 <- rank(allrslts[,size+2])

  bestv1 <- allrslts[which.max(rankv1),]
  bestv1mod <- c(names(data)[bestv1[1:size]+1], bestv1[(size+1):(2*size+3)])
  bestv1DF <- data.frame(matrix(bestv1mod[c(1:size)],ncol=size),
                         matrix(as.numeric(bestv1mod[(size+1):(2*size+3)]),nrow=1))
  names(bestv1DF) <- c(paste("Var",c(1:size),sep=""), "AUC1", "AUC2", "numNA", paste("Coef",c(1:size),sep=""))

  bestv2 <- allrslts[which.max(rankv1+rankv2),]
  bestv2mod <- c(names(data)[bestv2[1:size]+1], bestv2[(size+1):(2*size+3)])
  bestv2DF <- data.frame(matrix(bestv2mod[c(1:size)],ncol=size),
                         matrix(as.numeric(bestv2mod[(size+1):(2*size+3)]),nrow=1))
  names(bestv2DF) <- c(paste("Var",c(1:size),sep=""), "AUC1", "AUC2", "numNA", paste("Coef",c(1:size),sep=""))

  modrslt <- t(apply(allrslts, 1, function(x) c(names(data)[x[1:size]+1], x[(size+1):(2*size+3)])))
  modrslt.ord <- data.frame(modrslt[order(modrslt[,size+1], decreasing=TRUE),])
  modrslt.ordDF <- data.frame(modrslt.ord[,c(1:size)],
                              apply(modrslt.ord[,(size+1):(2*size+3)], 2, as.numeric))
  names(modrslt.ordDF) <- c(paste("Var",c(1:size),sep=""), "AUC1", "AUC2","numNA", paste("Coef",c(1:size),sep=""))

  if(length(which(rankv1==max(rankv1))) > 1){
    warning("tie(s) in AUC for D=K vs. D<K")
  }
  r1r2 <- rankv1+rankv2
  if(length(which(r1r2==max(r1r2))) > 1){
    warning("tie(s) in sum of ranks for AUCs")
  }

  return(list("Best.Single"=bestv1DF, "Best.Multi"=bestv2DF,
              "Ranked.Rslts"=modrslt.ordDF[1:nummod,]))
}
