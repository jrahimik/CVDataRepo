## Cross validation
# https://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1000/matrix/
library(RUVcorr)
library(glmnet)

Yind <- simulateGEdata(n=3000, m=1000, k=10, size.alpha=2,
                       corr.strength=5, g=NULL, Sigma.eps=0.1, 
                       nc=2000, ne=1000, intercept=TRUE, check=TRUE)

x <- Yind$Truth   
zeros<- which(apply(x,2,sd)==0)
x <- x[,-zeros]

x <- x[,1:100]
x <- scale(x,T,T)

write.table(x,
            file="geneExpressionData.csv"
            ,col.names=F
            ,row.names=F
            ,sep=" ")

y      <- 1*x[,60] + 0.9*x[,87]+rnorm(sd=1,n=dim(x)[1])

write.table(y,
            file="RA_Phen.csv"
            ,col.names=F
            ,row.names=F
            ,sep="")



#Create 10 equally size folds




ReportRep <- NULL
## Shuffling process


  
  
  



#Perform 10 fold cross validation




## Get the values for folds
GetFoldValue <- function(x,y,nfolds=10,lambdaVal=0.1){

suffle_id <- sample(nrow(x))
Xs        <- x[suffle_id,]
Ys        <- y[suffle_id,]  
reportFold <- NULL  
folds <- cut(seq(1,nrow(x)),breaks=10,labels=FALSE)  
  
for(i in 1:nfolds){
  
  testIndexes <- which(folds==i,arr.ind=TRUE)
  XsTest <- Xs[testIndexes,]
  YsTest <- Ys[testIndexes]
  
  XsTrain <- Xs[-testIndexes,]
  YsTrain <- Ys[-testIndexes]
  
  
  model<- glmnet(
    XsTrain,
    YsTrain,
    family = c("gaussian"),
    alpha = 1,
    lambda = lambdaVal)
  
  YsTesthat  <- predict(model,
                        newx = as.matrix(XsTest),
                        type="response")
  
  res        <- cbind(YsTesthat,YsTest)
  reportFold <- rbind(reportFold,res)
  
  
  
}

return(reportFold)

}



PlotData <- NULL
lambda_seq <- seq(0,2,0.1)

for(lambdaVal in lambda_seq){

res <- sapply(1:100,function(nothing){res <- GetFoldValue(x,y,nfolds=10,lambdaVal=lambdaVal);
return(cor(res[,1],res[,2]))})


PlotData <- cbind(PlotData,res)

}


colnames(PlotData) <- lambda_seq



boxplot(PlotData,ylab="spearman correlation",xlab=expression(lambda))


























