library(boot)
library(lme4)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
usage_cmd = "USAGE: Rscprit H2_script.R input_file output_file"
if (length(args) < 2){
	print(usage_cmd)
}
## cmd = args[1]
inputF = args[1]
outputF = args[2]

#data <- read.table("D:\\CONFIDENTIAL\\LINUX_DATA\\Amitava\\Northwest\\H2_exampleData.csv", sep = ',', header = TRUE)
data <- read.table(inputF, sep = ',', header = TRUE)

########################
### define functions ###
########################
# Heritability
# data is data frame that contains strain and Value column
# indicies are used by the boot function to sample from the 'data' data.frame
H2.test.boot <- function(data, indicies){
  
  d <- data[indicies,]
  
  pheno <- as.data.frame(dplyr::select(d, Value))[,1]

  Strain <- as.factor(d$Strain)
  
  
  reffMod <- lme4::lmer(pheno ~ 1 + (1|Strain))
  
  Variances <- as.data.frame(lme4::VarCorr(reffMod, comp = "Variance"))
  
  Vg <- Variances$vcov[1]
  Ve <- Variances$vcov[2]
  H2 <- Vg/(Vg+Ve)
  
  # errors <- sqrt(diag(lme4::VarCorr(reffMod, comp = "Variance")$strain))
  
  return(H2)
}

# data is data frame that contains strain and Value column
H2.test <- function(data){
  
  pheno <- as.data.frame(dplyr::select(data, Value))[,1]
  Strain <- as.factor(data$Strain)
  
  reffMod <- lme4::lmer(pheno ~ 1 + (1|Strain))
  
  Variances <- as.data.frame(lme4::VarCorr(reffMod, comp = "Variance"))
  
  Vg <- Variances$vcov[1]
  Ve <- Variances$vcov[2]
  H2 <- Vg/(Vg+Ve)
  
  # errors <- sqrt(diag(lme4::VarCorr(reffMod, comp = "Variance")$strain))
  
  return(H2)
}

# df is data frame that contains strain and Value column
#data$Strain = na.omit(data$Strain)
#data$Value = na.omit(data$Value)
#df <- as.data.frame(data)
#df <- dplyr::select(df,Strain,Value)
H2.calc <- function(data, boot = T){
  df <- as.data.frame(data)
  df <- dplyr::select(df,Strain,Value)
  
  if(boot == T){
    # bootstrapping with 1000 replications
    results <- boot(data=df, statistic=H2.test.boot, R=10000) # can reduce value to save time (500 is resonable most of the time). if you Error in bca.ci(boot.out, conf, index[1L], L = L, t = t.o, t0 = t0.o,  : estimated adjustment 'a' is NA, then you need to increase R value.
    
    # get 95% confidence interval
    ci <- boot.ci(results, type="bca")
    
    H2_errors <- data.frame(H2 = ci$t0, ci_l = ci$bca[4], ci_r = ci$bca[5])
    
    return(H2_errors)
    
  } else {
    
    H2 <- data.frame(H2 = H2.test(data = df), ci_l = NA, ci_r = NA)
    return(H2)
  }
  
}

####Example run
#H2.test.boot <- function(data, indicies)
#H2.test <- function(data)
#H2.calc(data, boot = T)

x<-H2.calc(data, boot = T)
write.table(x, outputF, sep = '\t')

#exportResults(H2.calc(data, boot = T), file.path(outputF))
#exportResults(H2.test(data), file.path(outputF))

