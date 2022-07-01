library(tidyverse)
library(broom)
library(glmnet)
library(Matrix)
library(glmnetUtils)
library(data.table)
library(hash)

bestpara <- function(output){
  number.of.alphas.tested <- length(output$alpha)
  cv.glmnet.dt <- data.table()
  for (i in 1:number.of.alphas.tested){
    glmnet.model <- output$modlist[[i]]
    min.mse <-  min(glmnet.model$cvm)
    min.lambda <- glmnet.model$lambda.min
    alpha.value <- output$alpha[i]
    new.cv.glmnet.dt <- data.table(alpha = alpha.value, min_mse = min.mse, min_lambda = min.lambda)
    cv.glmnet.dt <- rbind(cv.glmnet.dt,new.cv.glmnet.dt)
    return(cv.glmnet.dt[which.min(cv.glmnet.dt$min_mse)])
  }}

fl = '2_16_DEN_Exp_1'
wd = paste("/Users/keitokiddo/VC_revision/Simulations/", fl, sep='')
setwd(wd)


#########



mks = c(1, 16, 120, 560, 1820, 4368, 8008, 11440, 12870, 11440, 8008, 4368, 1820, 560, 120, 16, 1)
nrep = 3


rules = read.table(file = 'data/rules.tsv', sep = '\t', header = FALSE)
design <- sparseMatrix(rules[[1]],rules[[2]])
dim(design)


schemes = c('random', 'noised', 'mut')

trs = hash()
trs[['random']] <- 'tr_'
trs[['noised']] <- 'tr_'
trs[['mut']] <- 'tr_mut_'

fs = hash()
fs[['random']] <- 'f_'
fs[['noised']] <- 'f_noised_'
fs[['mut']] <- 'f_'


## Three-way
modelname = "3way"
k = 3

for (m in 3:3) {
  scheme = schemes[m]
  name = paste(fl, scheme, sep='_')
  print(name)
  for (j in 1:3){
    for (i in 1:3){
      outfilename = paste(name, modelname, j, i, sep="_")
      print(i, j)
      tr = as.integer(unlist(read.table(file = paste('data/', trs[[scheme]], i,
                                                     ".tsv",sep = ''), header = FALSE)))
      if (m==2) {
        response = read.table(file = paste('data/', fl, "_", fs[[scheme]],  j, "_", "0.2",  
                                           ".tsv", sep = ''), sep = '\t', header = FALSE)[[1]]
      } else{
        response = read.table(file = paste('data/', fl, "_", fs[[scheme]],  j, 
                                           ".tsv", sep = ''), sep = '\t', header = FALSE)[[1]]
      }
      
      cvfit <- cv.glmnet(design[tr,1: sum(mks[1:(k+1)])], response[tr], alpha = 0, parallel = TRUE)
      pd <- predict(cvfit, newx = design[,1: sum(mks[1:(k+1)])], s = "lambda.min")
      
      write.table(as.matrix(pd), file = paste('out/',outfilename, ".tsv", sep = ""),
                  row.names=FALSE, quote = FALSE, sep = '\t', col.names = FALSE)
    }
  }
}
end_time <- Sys.time()
end_time - start_time




## Pairwise
modelname = "2way"
k = 2
start_time <- Sys.time()

schemes


for (m in 3:3) {
  scheme = schemes[m]
  name = paste(fl, scheme, sep='_')
  print(name)
  for (j in 1:3){
    for (i in 1:3){
      outfilename = paste(name, modelname, j, i, sep="_")
      print(i, j)
      tr = as.integer(unlist(read.table(file = paste('data/', trs[[scheme]], i,
                                                     ".tsv",sep = ''), header = FALSE)))
      if (m==2) {
        response = read.table(file = paste('data/', fl, "_", fs[[scheme]],  j, "_", "0.2",  
                                           ".tsv", sep = ''), sep = '\t', header = FALSE)[[1]]
      } else{
        response = read.table(file = paste('data/', fl, "_", fs[[scheme]],  j, 
                                           ".tsv", sep = ''), sep = '\t', header = FALSE)[[1]]
      }
      
      cvfit <- cv.glmnet(design[tr,1: sum(mks[1:(k+1)])], response[tr], alpha = 0, parallel = TRUE)
      pd <- predict(cvfit, newx = design[,1: sum(mks[1:(k+1)])], s = "lambda.min")
      
      write.table(as.matrix(pd), file = paste('out/',outfilename, ".tsv", sep = ""),
                  row.names=FALSE, quote = FALSE, sep = '\t', col.names = FALSE)
    }
  }
}
end_time <- Sys.time()
end_time - start_time
