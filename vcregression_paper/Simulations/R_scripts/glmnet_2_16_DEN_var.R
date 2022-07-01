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

fl = '2_16_DEN'
wd = paste("/Users/keitokiddo/VC_revision/Simulations/", fl,"_var_ldas", sep='')
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
blist = c(10, 3.29584, 2.19722, 1.09861)

scheme = schemes[3]

subreps = c(1,2,8)

for (n in 1:length(blist)) {
  name = paste(fl, scheme, sep='_')
  for (j in subreps){
    for (i in 2:2){
      outfilename = paste(name, modelname, blist[n], j, i, sep="_")
      print(outfilename)
      tr = as.integer(unlist(read.table(file = paste('data/', trs[[scheme]], i,
                                                     ".tsv",sep = ''), header = FALSE)))
      
      response = read.table(file = paste('data/', fl, "_", "f_", as.character(blist[n]), "_",  j,".tsv", sep = '')
                            , sep = '\t', header = FALSE)[[1]]
      
      cvfit <- cv.glmnet(design[tr,1: sum(mks[1:(k+1)])], response[tr], alpha = 0, parallel = TRUE)
      pd <- predict(cvfit, newx = design[,1: sum(mks[1:(k+1)])], s = "lambda.min")

      write.table(as.matrix(pd), file = paste('out/',outfilename, ".tsv", sep = ""),
                  row.names=FALSE, quote = FALSE, sep = '\t', col.names = FALSE)
    }
  }
}




## Pairwise
modelname = "2way"
k = 2
start_time <- Sys.time()

for (n in 1:length(blist)) {
  name = paste(fl, scheme, sep='_')
  for (j in subreps){
    for (i in 2:2){
      outfilename = paste(name, modelname, blist[n], j, i, sep="_")
      print(outfilename)
      tr = as.integer(unlist(read.table(file = paste('data/', trs[[scheme]], i,
                                                     ".tsv",sep = ''), header = FALSE)))
      
      response = read.table(file = paste('data/', fl, "_", "f_", as.character(blist[n]), "_",  j,".tsv", sep = '')
                            , sep = '\t', header = FALSE)[[1]]
      
      cvfit <- cv.glmnet(design[tr,1: sum(mks[1:(k+1)])], response[tr], alpha = 0, parallel = TRUE)
      pd <- predict(cvfit, newx = design[,1: sum(mks[1:(k+1)])], s = "lambda.min")
      
      write.table(as.matrix(pd), file = paste('out/',outfilename, ".tsv", sep = ""),
                  row.names=FALSE, quote = FALSE, sep = '\t', col.names = FALSE)
    }
  }
}

