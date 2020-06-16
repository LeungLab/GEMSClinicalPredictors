library(tidyverse)
library(lubridate)
library(pROC)
library(fields)
library(zoo)
library(ks)
library(KernSmooth)
library(ranger)
library(broom)
library(furrr)
library(cvAUC)
library(slider);

dat_joined <- readRDS("dat_afe.RDS") # loads data file with clean data and weather information included and etiology defined 

names=colnames(dat_joined)[c(5:8,10:92,94:134,136,138,141:165,189)] # specify column numbers of variables included in study (see supplement table for those included)

resp_var="viral_only" #define response variable 

# This calculates overall importance using random forest 
out=ranger(as.formula(paste(resp_var,'~',paste(names,collapse="+"),sep="")),data=dat_joined,num.trees=1000,importance="impurity")
imps=importance(out)
df_imps=data.frame(names=names(imps),var_red=as.numeric(imps)) %>% arrange(desc(var_red))
head(df_imps,10)



# partial dependency plots 
out=ranger(as.formula(paste(resp_var,'~',paste(df_imps$names[1:10],collapse="+"),sep="")),data=dat_joined,num.trees=1000,importance="impurity")
pd_plots=purrr::map(as.character(df_imps$names[1:10]),~out %>% partial(pred.var=c(.x),plot.engine = "ggplot2",plot=T,rug=T))
pd_plots=map2(pd_plots,
              c("Age (mo.)","Season","Blood in Stool","HAZ","Breastfed","Vomiting","MUAC","Resp. Rate (per min.)","Wealth Index","Temperature (°C)"), 
              function(x,y) x + theme_bw() + ylab("Viral Prediction"))
grid.arrange(grobs=pd_plots,top="Partial Dependecy Plots of Top Ten Predictors")

##Get importance for each site

outs=dat_joined %>% split(.$site) %>% purrr::map(~ranger(as.formula(paste(resp_var,'~',paste(names,collapse="+"),sep="")),data=.,num.trees=1000,importance="impurity"))
dfs_imps=purrr::map(outs,~importance(.)) %>% purrr::map(~data.frame(names=names(.),var_red=as.numeric(.)) %>% arrange(desc(var_red)))
bind_cols(dfs_imps %>% purrr::map(~.[1:11,1]))
###

#Define number of variables sets in cross validation 
nvars_opts=c(1:10,15,20,30,40,50)

result=data.frame(site=NA,iter=NA,nvar=NA,true=NA,pred_glm=NA,pred_RF=NA)


CV_func=function(x,y){
  train=dat_joined %>% sample_frac(.80,replace=F)
  test=dat_joined[-which(dat_joined$index %in% train$index),]
  
  out=ranger(as.formula(paste(resp_var,'~',paste(names,collapse="+"),sep="")),data=train,num.trees=1000,importance="impurity")
  df_imps=data.frame(names=names(ranger::importance(out)),imps=ranger::importance(out)) %>% arrange(desc(imps))
  nvars=x
  each=y
  out1=glm(as.formula(paste(resp_var,'~',paste(df_imps$names[1:nvars],collapse="+"),sep="")),data=train,family="binomial")
  out2=ranger(as.formula(paste(resp_var,'~',paste(df_imps$names[1:nvars],collapse="+"),sep="")),data=train,num.trees=1000)
  df=data.frame(site=test$site,iter=each,nvar=nvars,true=test[[resp_var]],pred_glm=as.numeric(predict(out1,newdata=test,type="response")),pred_RF=as.numeric(predict(out2,data=test,type="response")$predictions))
  df
}


plan(multiprocess)
result=cross2(1:100,nvars_opts) %>% future_map(~CV_func(.[[2]],.[[1]]),.progress=T)


AUCs=bind_rows(result) %>% split(.$nvar) %>% purrr::map(~ci.cvAUC(.$pred_glm,.$true,folds=.$iter))
AUCs2=bind_rows(result) %>% split(.$nvar) %>% purrr::map(~ci.cvAUC(.$pred_RF,.$true,folds=.$iter))


#Calibration
list_res=bind_rows(result) %>% split(.$nvar) %>% purrr::map(. %>% arrange(pred_RF))
#this takes a long time
plan(multiprocess)
resVobs=map2(list_res,list_res %>% purrr::map(~.$pred_RF), function(x,y) future_map(y, function(z)
  x %>% filter(pred_RF>z-.05,pred_RF<z+.05) %>% summarize(mean(true)),.progress=T))

list_res=map2(list_res,resVobs,function(x,y) x %>% mutate(observed=unlist(y)))

ggplot(bind_rows(list_res) %>% filter(nvar %in% 2:5),aes(x=pred_glm,y=observed,group=as.factor(nvar),color=as.factor(nvar))) + geom_line() + 
  geom_segment(aes(x=0,y=0,xend=1,yend=1),color="black",size=1.1)

calib=bind_rows(list_res %>% purrr::map(~summary((lm(observed~pred_glm,data=.)))$coefficients[1:2,1]))
calib$Parameter=c("Intercept","Slope")
calib_tbl=bind_rows(calib) %>% pivot_longer(cols=1:15) %>% pivot_wider(id_cols=2,names_from="Parameter")

# get logistic regression magnitudes of effect 
tbl=summary(glm(as.formula(paste(resp_var,'~',paste(df_imps$names[c(1:5)],collapse="+"),sep="")),data=dat_joined ,family="binomial"))

adj_odds=data.frame(Adjusted_Odds=exp(tbl$coefficients[,1]), CI=
                      t(bind_rows(map2(tbl$coefficients[,1],tbl$coefficients[,2], function(x,y) exp(x + c(-1,1)*y*qnorm(.975))))),
                    P_value=tbl$coefficients[,4])
adj_odds$CI=paste(round(adj_odds$CI.1,4),"-",round(adj_odds$CI.2,4))

coef=data.frame(Coef=(exp(tbl$coefficients[,1])), CI=
                  t(bind_rows(map2((tbl$coefficients[,1]),(tbl$coefficients[,2]), function(x,y) exp(x + c(-1,1)*y*qnorm(.975))))),
                P_value=tbl$coefficients[,4])
coef$CI=paste("(",round(coef$CI.1,4)," - ",round(coef$CI.2,4),")",sep="")
xtable(coef %>% select(Coef,CI,P_value),digits=4)
coef %>% select(Coef,CI, P_value)
adj_odds$CI=paste(round(adj_odds$CI.1,4),"-",round(adj_odds$CI.2,4))
adj_odds %>% select(Adjusted_Odds,CI,P_value)




