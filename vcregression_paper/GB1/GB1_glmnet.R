# install.packages("glmnet", repos = "http://cran.us.r-project.org")
# install.packages("Matrix", repos = "http://cran.us.r-project.org")


library(tidyverse)
library(broom)
library(glmnet)
library(Matrix)
library(glmnetUtils)
library(data.table)


bestpara <- function(output){
	number.of.alphas.tested <- length(output$alpha)
	cv.glmnet.dt <- data.table()
	for (i in 1:number.of.alphas.tested){
	  glmnet.model <- output$modlist[[i]]
	  min.mse <-  min(glmnet.model$cvm)
	  min.lambda <- glmnet.model$lambda.min
	  alpha.value <- output$alpha[i]
	  new.cv.glmnet.dt <- data.table(alpha=alpha.value,min_mse=min.mse,min_lambda=min.lambda)
	  cv.glmnet.dt <- rbind(cv.glmnet.dt,new.cv.glmnet.dt)
  return(cv.glmnet.dt[which.min(cv.glmnet.dt$min_mse)])
}}


wd="/Users/keitokiddo/Dropbox/VC/"
wd
setwd(paste(wd,"GB1/GB1_glmnet",sep=""))
#########

response=read.table(file = '../data/response.txt', sep = '\t', header = FALSE)[[1]]
seqlist=as.integer(unlist(read.table(file='../data/seqlist.txt',sep = '\t', header = FALSE)))


rules=read.table(file = 'rules4.tsv', sep = '\t', header = FALSE)
design <- sparseMatrix(rules[[1]],rules[[2]])
dim(design)

mks=c(1, 80, 2400, 32000, 160000)
nrep=3



k=2
modelname="elasticnet_pairwise"
start_time <- Sys.time()

for (i in 1:18){
	for (j in 1:3){

		tr=as.integer(unlist(read.table(file=paste('../data/tr_',i,"_",j,".tsv",sep=''), header = FALSE)))
				
		cva_out <- cva.glmnet(design[tr,1: sum(mks[1:(k+1)])], response[tr])
		
		pd <- predict(cva_out, design[,1: sum(mks[1:(k+1)])], alpha=bestpara(cva_out)$alpha)
			
		write.table(pd, file=paste("pd_",modelname,"_", i, "_", j, ".tsv", sep=""), quote=FALSE, sep='\t', col.names = NA)
		
	}
}
end_time <- Sys.time()

end_time - start_time


k=3
modelname="elasticnet_3way"
candidatelambdas = lambda=10^ seq(-6,-0,.5)

start_time <- Sys.time()

for (i in 1:18){
	for (j in 1:3){

		tr=as.integer(unlist(read.table(file=paste('../data/tr_',i,"_",j,".tsv",sep=''), header = FALSE)))
				
		cva_out <- cva.glmnet(design[tr,1: sum(mks[1:(k+1)])], response[tr],lambda=candidatelambdas)
		pd <- predict(cva_out, design[,1: sum(mks[1:(k+1)])], alpha=bestpara(cva_out)$alpha)
			
		write.table(pd, file=paste("pd_","modelname","_", i, "_", j, ".tsv", sep=""), quote=FALSE, sep='\t', col.names = NA)
		
	}
}
end_time <- Sys.time()

end_time - start_time
