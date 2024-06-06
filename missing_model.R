#Calculating missing model
library(dplyr)
library(randomForest)
#Input: Covariates, Missing indicator
#Output: mis.model = Missing model, weights = inverse probability of being observed for each individual, rowid = indicator of the complete rows.

mmis <- function(data,name.outcome,form=NULL){#Only outcome has missing. Input data has "r" indicating observed data.
  #data$mis <- 1-is.na(data[,name.outcome])
  num.del <- which(names(data) %in% c(name.outcome,"r"))
  #build model(RF CLASSIFICATION? OR PARAMETRIC MODEL?)
  if(is.null(form)){
  form <- as.formula(paste("factor(r)","~",
                           paste(names(data)[-num.del],collapse = "+"),collapse = ""))
  }
  mis.model <- glm(form,data = data,family=binomial(link = "logit"))
  #mis.model <- randomForest(form,data = data,)
  pis <- predict(mis.model,newdata = data,type = "response")
  pis <- ifelse(pis==0,pis+0.001,pis)#some of the weights are 0
  return(pis)
}




