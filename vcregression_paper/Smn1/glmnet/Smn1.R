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
	  new.cv.glmnet.dt <- data.table(alpha = alpha.value,min_mse = min.mse,min_lambda = min.lambda)
	  cv.glmnet.dt <- rbind(cv.glmnet.dt,new.cv.glmnet.dt)
  return(cv.glmnet.dt[which.min(cv.glmnet.dt$min_mse)])
}}


wd = "/Users/keitokiddo/Dropbox/VC/Smn1"

setwd(paste(wd,"/glmnet",sep = ""))
#########
response = read.table(file = 'response.tsv', sep = '\t', header = FALSE)[[1]]


mks = c(1, 32, 448, 3584, 17920, 57344, 114688, 131072, 65536)
nrep = 5


rules = read.table(file = 'rules.tsv', sep = '\t', header = FALSE)
design <- sparseMatrix(rules[[1]],rules[[2]])
dim(design)


## Three-way
modelname = "3way_EN"
k = 3

start_time <- Sys.time()

for (i in 1:11){
	for (j in 1:5){

		tr = as.integer(unlist(read.table(file = paste('../data/tr_',i,"_",j,".tsv",sep = ''), header = FALSE)))
				
		cva_out <- cva.glmnet(design[tr,1: sum(mks[1:(k+1)])], response[tr])
		best.params = bestpara(cva_out)
		
		pd <- predict(cva_out, design[,1: sum(mks[1:(k+1)])], alpha = best.params$alpha)
			
		write.table(as.matrix(pd), file = paste("pd_",modelname,"_", i, "_", j, ".tsv", sep = ""), quote = FALSE, sep = '\t', col.names = NA)
		
	}
}

end_time <- Sys.time()

end_time - start_time
#Imputation of full landscape using all data
tr = as.integer(unlist(read.table(file = paste('../data/seqlistSmn1',".tsv",sep = ''), header = FALSE)))
print(tr)			
cva_out <- cva.glmnet(design[tr,1: sum(mks[1:(k+1)])], response[tr])
best.params = bestpara(cva_out)
		
pd <- predict(cva_out, design[,1: sum(mks[1:(k+1)])], alpha = best.params$alpha)
			
write.table(as.matrix(pd), file = paste("pd_",modelname,"_", "ALL", ".tsv", sep = ""), quote = FALSE, sep = '\t', col.names = NA)


## Pairwise
modelname = "2way_EN"
k = 2

start_time <- Sys.time()

for (i in 1:11){
	for (j in 1:5){

		tr = as.integer(unlist(read.table(file = paste('../data/tr_',i,"_",j,".tsv",sep = ''), header = FALSE)))
				
		cva_out <- cva.glmnet(design[tr,1: sum(mks[1:(k+1)])], response[tr])
		best.params = bestpara(cva_out)
		
		pd <- predict(cva_out, design[,1: sum(mks[1:(k+1)])], alpha = best.params$alpha)
			
		write.table(as.matrix(pd), file = paste("pd_",modelname,"_", i, "_", j, ".tsv", sep = ""), quote = FALSE, sep = '\t', col.names = NA)
		
	}
}

end_time <- Sys.time()

end_time - start_time

