# Variableselection in metabolomics study #

#####################################################
## LASSO MODEL + bootstrap loop for ALL VARIABLES ####
######################################################

############
## Open data
#############
library(boot);library(caret); library(penalized); library(mice); library("missForest", lib.loc="~/R/win-library/3.2")

rm(list=ls())

df <- read.table("C:/dane/renata bujak/PAH project/PAH_matryca.txt", header=TRUE,sep="\t"); dim(df)
df[, c(3:840)] <- sapply(df[, c(3:840)], as.numeric)
table(df$pah)
df$pah <- as.factor(df$pah);
xx2 <- df[sample(nrow(df)),]# 
mydata = xx2
table(mydata[,2])

##############################
# Missing value imputation ###
##############################

# random forrest imputation ##
set.seed(81)
snp.imp <- missForest(mydata, verbose = TRUE)
snp.imp$OOBerror
done <- snp.imp$ximp

mydata = done

# SNP selected by LASSO
Results_d <- rep(NA, nrow = 838)

## creating the loop to calculate reproducibility
# SNP selected by LASSO #

mz1 <- 0
mz2 <- 0
mz3 <- 0
mz4 <- 0
mz5 <-0

## creating the loop to calculate reproducibility

for(i in 1:1000){
  print(i)
  set.seed(i+1)
  
  mydata <- na.omit(mydata)
  ## random permutation
  msk <- sample(c(rep(trunc(length(mydata[,2])))),replace=TRUE) # bootstrap resampling with replacement
  snp_ind_var.boot <- mydata[msk,]
  head(snp_ind_var.boot) #randomly shuffled ID
  
  prof.snp <- try(profL1(mydata[,2] ~ as.matrix(mydata[,3:840]), 
                         model="logistic",data = mydata,
                         standardize = T, fold = 5)); 
  
  if(class(prof.snp) == "try-error") next
  #names(prof.snp)
  #plotpath(prof.snp$fullfit)
  
  
  opt.snp <- try(optL1(mydata[,2]~ as.matrix(mydata[,3:840]), 
                       model="logistic",data = mydata,
                       standardize = F, fold = prof.snp$fold));
  
  if(class(opt.snp) == "try-error") next#if error then stop
  
  snp.lasso <- tryCatch(penalized(mydata[,2]~ as.matrix(mydata[,3:840]), 
                                  model="logistic",data = mydata, standardize = F,
                                  lambda1 = opt.snp$lambda)); snp.lasso
  res <- rep(NA,840)
  res <- as.vector(snp.lasso@penalized)
  
  
  # Results - SNPs LASSO
  
  if(res[30]!=0){mz1 = mz1 + 1}# if B of this snp is not equal to 0 count reproducibility
  if(res[186]!=0){mz2 = mz2 + 1}#
  if(res[262]!=0){mz3 = mz3 + 1}#
  if(res[836]!=0){mz4 = mz4 + 1}#
  if(res[837]!=0){mz5 = mz5 + 1}#
}

print(mz1);
print(mz2)
print(mz3) 
print(mz4) 
print(mz5) #control