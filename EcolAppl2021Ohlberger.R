##========================================================================##
##                                                                        ##
##   Code for analyses in Ohlberger et al. 2021 Ecological Applications   ##
##                                                                        ##
##========================================================================##
pkgs<-c("readxl","tidyverse","visreg","MuMIn","nlme","vegan","performance", "akima","Hmisc") 
install.packages(setdiff(pkgs,rownames(installed.packages()))) 
invisible(lapply(pkgs,library,character.only=T))
'%!in%'<-function(x,y)!('%in%'(x,y))

##========================================================================##
##===================================================================## data
##========================================================================##
#cod data from ICES 2020 report of the Arctic Fisheries Working Group (AFWG)

##======================================================## age-specific data
ages<-seq(3,14,1)
years<-seq(1946,2020,1)
data<-data.frame(Year=years)
##---------------------------------------------------------## numbers-at-age
## AFWG 2020 Tab3-16 numbers at age (model estimate for age-3 in 2020)
NaA<-data.frame(read_excel("data_cod_number_at_age_2020.xlsx"))
NaA<-NaA %>% mutate(across(where(is.character),str_trim)) 
NaA<-mutate_all(NaA,function(x) as.numeric(as.character(x)))
NaA<-NaA %>% filter(Year %in% data$Year) %>% select(-Year,-plusgroup,-TOTAL)
data$Recruitment<-NaA$age3 ## recruits are age-3
##----------------------------------------------------------## mature-at-age
## AFWG 2020 Tab3-11 proportion mature at age
MaA<-data.frame(read_excel("data_cod_prop_mat_at_age_2020.xlsx"))
MaA<-MaA %>% mutate(across(where(is.character),str_trim)) 
MaA<-mutate_all(MaA,function(x) as.numeric(as.character(x)))
MaA<-MaA %>% filter(Year %in% data$Year) %>% select(-Year,-plusgroup)
##----------------------------------------------------------## weight-at-age
## AFWG 2020 Tab3-9 stock weights at age (kg)
WaA<-data.frame(read_excel("data_cod_weight_at_age_2020.xlsx"))
WaA<-WaA %>% mutate(across(where(is.character),str_trim)) 
WaA<-mutate_all(WaA,function(x) as.numeric(as.character(x)))
WaA<-WaA %>% filter(Year %in% data$Year) %>% select(-Year,-plusgroup)

##==============================================================## cannibals
## cannibals ages 3-6 (cod mature into spawning stock after age-6)
Cage<-paste0("age",c(3,4,5,6))
nY<-dim(data)[1]
data$CannibalsY<-NA
for(i in 1:nY) { data$CannibalsY[i]<-sum(NaA[i,names(NaA) %in% Cage]* WaA[i,names(WaA) %in% Cage]) }
##-------------------------------------------------------------## 1-year lag
## highest mortality on age-1 victims (Yaragina et al. 2009)
data$Cannibals<-NA
for(i in 1:(nY-1)) { data$Cannibals[i]<-sum(NaA[i+1,names(NaA) %in% Cage]* WaA[i+1,names(WaA) %in% Cage]) }

##==================================================## total spawner biomass
SpNaA<-round(NaA*MaA) ## spawner numbers-at-age
SpBaA<-round(NaA*MaA*WaA) ## spawner biomass-at-age
SpBio<-rowSums(SpBaA) ## total spawner biomass
data$TotSpBio<-SpBio

##==================================================## spawner age structure
## biomass-weighted mean age, weight and age diversity
##------------------------------------------------------------## mean weight
data$MeanWeight<-rowSums(SpBaA*WaA)/rowSums(SpBaA)
##---------------------------------------------------------------## mean age
data$MeanAge<-apply(SpBaA,1,function(x) sum(x*ages)/sum(x))
##----------------------------------------------------------## age diversity
data$AgeDiv<-apply(SpBaA,1,function(x) diversity(x,index="shannon"))
##---------------------------------------------## proportion repeat spawners
## without first year >> may need to re-run model selection if not included
data$pRepeat<-NA
for(y in 2:dim(MaA)[1]) {
  repeat_a<-NA;repeat_a[1]<-0
  for(a in 2:dim(MaA)[2]) {
    repeat_a[a]<-(MaA[y-1,a-1]/MaA[y,a])*SpBaA[y,a]
  }
  data$pRepeat[y]<-sum(repeat_a,na.rm=T)/SpBio[y]
} 

##===========================================================## recruits/SSB
## 'Recruitment' by sampling year; 'Recruits' by spawn year
data$Recruits<-c(data$Recruitment[-c(1:3)],rep(NA,3))
data$RpSSB<-c(data$Recruits/data$TotSpBio)
data$lnRpSSB<-log(data$RpSSB)
alldata<-data ## includes all years

##============================================================## temperature
## Kola section (0-200m average at stations 3-7, source: www.pinro.ru)
Kola<-data.frame(read_excel("data_Kola_temp_months.xlsx")) 
Kola<-mutate_all(Kola,function(x) as.numeric(as.character(x)))
Kola<-Kola %>% filter(Year %in% data$Year)
##-----------------------------------------------## interpolate missing data
## Aug-Dec data missing in 2016
rowind<-which(!complete.cases(Kola))
colind<-which(is.na(Kola[rowind,]))
df<-data.frame(expand_grid(Year=Kola$Year,Month=seq(12)))
df$Temp<-as.numeric(unlist(t(Kola[,-1])))
df<-df[complete.cases(df),]  
yrs<-Kola$Year[rowind]
interpolated<-interp(df$Year,df$Month,df$Temp,xo=yrs,yo=c(colind-1))
Kola[rowind,colind]<-interpolated$z
##----------------------------------------------## select temperature metric
Kola$MeanTemp<-rowMeans(Kola[,names(Kola) %in% c("Jul","Aug","Sep","Oct","Nov","Dec")]) ## larval and 0-group temps
Kola<-Kola %>% select(Year,MeanTemp)
data<-merge(data,Kola,by="Year",all=T)

##========================================================================##
##=============================================## multiple linear regression
##========================================================================##
data<-data %>% filter(!is.na(RpSSB)) ## drop years without recruits/SSB
data<-data %>% select(-pRepeat) ## without proportion repeat spawners?
data<-data[complete.cases(data),] ## only complete cases/years
nY<-dim(data)[1];nY

##--------------------------------------## center/standardize to mean=0/sd=1
data$lnRecruits<-log(data$Recruits) ## for test model
data_unscaled<-dplyr::select(data,-Year,-RpSSB,-lnRpSSB,-Recruits,-Recruitment,-lnRecruits) 
data_other<-dplyr::select(data,Year,RpSSB,lnRpSSB,Recruits,lnRecruits)
data_means<-sapply(data_unscaled,mean)
data_sds<-sapply(data_unscaled,sd)
data<-data.frame(cbind(data_other,scale(data_unscaled))) 

##-----------------------------------------------------## threshold function
## function to select threshold in mean age dummy used in interaction models
threshold.gls<-function(form=form,newdata,a=0.1,b=0.9,step=0.05){
  thd.vals<-seq(from=a,to=b,by=step)
  aic.vals<-rep(NA,length(thd.vals))
  mod.list<-list()
  length(mod.list)<-length(thd.vals)
  for(i in 1:length(thd.vals)){
    cut<-thd.vals[i]
    newdata$MeanAgeDummy[newdata$MeanAge<quantile(newdata$MeanAge,cut)]<-"low"
    newdata$MeanAgeDummy[newdata$MeanAge>=quantile(newdata$MeanAge,cut)]<-"high"
    m<-gls(form,data=newdata,correlation=corAR1(),method="ML")
    aic.vals[i]<-AICc(m)
    mod.list[[i]]<-m
  }
  thd.select<-thd.vals[aic.vals==min(aic.vals)]
  mod.select<-mod.list[[i]]
  return(list(thd.select=thd.select,mod.select=mod.select))
}

##=============================================## multiple linear regression
options(na.action="na.fail") 
##-------------------------------------------------------------## full model
# mod_form<-form<-formula(lnRpSSB~TotSpBio+MeanTemp+I(MeanTemp^2)+MeanAge+pRepeat+Cannibals+MeanTemp:MeanAgeDummy) ## proportion repeat spawners
mod_form<-form<-formula(lnRpSSB~TotSpBio+MeanTemp+I(MeanTemp^2)+MeanAge+ Cannibals+MeanTemp:MeanAgeDummy) ## all years
##-------------------------------------------## threshold for mean age dummy
test_mod<-threshold.gls(form=form,newdata=data)
cut<-test_mod$thd.select
data$MeanAgeDummy[data$MeanAge<quantile(data$MeanAge,cut)]<-"low" 
data$MeanAgeDummy[data$MeanAge>=quantile(data$MeanAge,cut)]<-"high"
##--------------------------------------------------------## model selection
mod_full<-lm(mod_form,data=data) 
mod_select<-dredge(mod_full,fixed=c("TotSpBio"),trace=F,rank="AICc")
mod_sel<-get.models(mod_select,subset=1)[[1]] ## dredge selection
car::vif(mod_sel) ## does not suggest any problems with multicollinearity 
# summary(mod_sel) 
##----------------------------------------------## check for autocorrelation
as.numeric(pacf(residuals(mod_sel),lag=15,plot=F)$acf) ## significant lag-1

##============================================## account for autocorrelation
mod_ar1_full<-gls(mod_form,data=data,correlation=corAR1(),method="ML")
mod_select_ar1<-dredge(mod_ar1_full,fixed=c("TotSpBio"),trace=F,rank="AICc")
mod_sel_ar1<-get.models(mod_select_ar1,subset=1)[[1]]
AICc(mod_sel,mod_sel_ar1)
compare_performance(mod_ar1_full,mod_sel_ar1) ## other performance metrics
# summary(mod_sel_ar1)

##=========================================================## model variance
mod_1<-update(mod_sel_ar1,weights=varFixed(~TotSpBio))
mod_2<-update(mod_sel_ar1,weights=varFixed(~MeanAge))
mod_3<-update(mod_sel_ar1,weights=varFixed(~MeanTemp))
mod_4<-update(mod_sel_ar1,weights=varFixed(~Cannibals))
mod_5<-update(mod_sel_ar1,weights=varExp(form=~TotSpBio))
mod_6<-update(mod_sel_ar1,weights=varExp(form=~MeanAge))
mod_7<-update(mod_sel_ar1,weights=varExp(form=~MeanTemp))
mod_8<-update(mod_sel_ar1,weights=varExp(form=~Cannibals))
##--------------------------------------------------------------## AIC table
aics<-AICc(mod_sel_ar1,mod_1,mod_2,mod_3,mod_4,mod_5,mod_6,mod_7,mod_8)
aics$delta<-round(aics$AICc-min(aics$AICc),4)
rownames(aics)[aics$delta==0] ## choose simpler model without weights 

##===================## AIC with penalty for selecting interaction threshold
## (Ottersen et al 2013)
nmod<-dim(data.frame(mod_select_ar1))[1]
myAICc<-NA
for(i in 1:nmod) { 
  my_mod<-get.models(mod_select_ar1,subset=i) 
  myterms<-strsplit(as.character(my_mod[[1]]$call)[2]," ")[[1]]
  myterms<-myterms[myterms %!in% c("+","~","1","lnRpSSB")]
  n<-length(myterms)
  if("MeanAgeDummy:MeanTemp" %in% myterms) { 
    my_AICc<-AICc(my_mod[[1]],k=2+1/n) 
  } else { 
    my_AICc<-AICc(my_mod[[1]],k=2) } 
  myAICc[i]<-my_AICc
}

##==============================================================## AIC table
aic_table<-data.frame(mod_select_ar1)
##---------------------------------------------## model forms
num_of_mods<-dim(aic_table)[1] 
mod_forms<-get.models(mod_select_ar1,subset=delta<100)
mod_forms<-mod_forms[1:num_of_mods]
mod_forms_topX<-NA
for(i in 1:num_of_mods) {
  use_mod<-mod_forms[[i]]	
  form<-as.character(use_mod$call)[2]
  mod_forms_topX[i]<-form
}
##---------------------------------------------## edit vars in model forms
mod_forms_topX<-unlist(mod_forms_topX)
mod_forms_topX<-gsub("lnRSSB","ln(R/SSB)",mod_forms_topX)
mod_forms_topX<-gsub("TotSpBio","S",mod_forms_topX)
mod_forms_topX<-gsub("Cannibals","C",mod_forms_topX)
mod_forms_topX<-gsub("MeanAge","A",mod_forms_topX)
mod_forms_topX<-gsub("AgeDiv","D",mod_forms_topX)
mod_forms_topX<-gsub("MeanTemp","T",mod_forms_topX)
mod_forms_topX<-gsub("ADummy:T","A:T",mod_forms_topX)
mod_forms_topX<-gsub("1","",mod_forms_topX)
mod_forms_topX<-gsub(" ","",mod_forms_topX)
mod_forms_topX<-gsub("\\++","+",mod_forms_topX)
mod_covars<-mod_forms_topX
mod_covars<-lapply(mod_covars,function(x) substr(x,9,nchar(x)))
mod_covars<-unlist(mod_covars)
##-----------------------------------------------------## add back to table
aic_table$mod_covars<-mod_covars
aic_table$mod_forms<-mod_forms_topX
##--------------------------------------------------------## new AICc values
aic_table<-aic_table[c(1:7,13)]
aic_table$AICc_penalized<-round(myAICc,2)
aic_table$delta_AICc<-aic_table$AICc_penalized-min(aic_table$AICc_penalized)
##-----------------------------------------------------------## order table
aic_table<-aic_table[order(aic_table$AICc_penalized),]
aic_table$delta<-round(aic_table$delta,1)
#------------------------------------------------------## cumulative weights
aic_table$weights<-round(Weights(aic_table$AICc_penalized),3)
aic_table$cum.weights<-cumsum(aic_table$weights)
# write.csv(aic_table,"table_AICs.csv")

##========================================================================##
##=================## model cross-validation using out-of-sample predictions
##========================================================================##

##---------------------------------------------------------## models to test
mod<-update(mod_sel_ar1,method="ML")
mod_select<-dredge(mod,trace=F,rank="AICc")
mod_list<-get.models(mod_select,subset=delta<1e3) 
nM<-length(mod_list)
test_forms<-list()
cnt<-1
for(i in 1:nM) {
  mymod<-mod_list[[i]]
  terms<-attributes(mymod)$namBetaFull[-1]
  if("TotSpBio" %in% terms){ ## only test models with TotSpBio as fixed term
    nterms<-length(terms)
    terms<-gsub("MeanAgeDummylow","MeanAgeDummy",terms)
    if(nterms!=0) my_mod<-formula("lnRpSSB~1") ## intercept only model
    if(nterms!=0) my_mod<-formula(paste("lnRpSSB~",paste0(terms,collapse="+")))
    test_forms[[cnt]]<-my_mod 
    cnt<-cnt+1
  }
}
nMod<-length(test_forms)
##------------------------------------## randomly sample train and test data
nS<-1e3 ## (~1 minute per 100 samples)
seeds<-sample(seq(1e6),1e4)
RMSE<-array(dim=c(nS,nMod))
formulas<-list()
start<-Sys.time()
cnt<-1
for(i in 1:nS) {
  set.seed(seeds[cnt])
  train<-sort(sample(seq(nY),round(0.75*nY),replace=F)) ## use 75% to train
  traindata<-data[train,] 
  test<-seq(nY)[seq(nY) %!in% train]
  testdata<-data[test,]
  ##------------------------------------------------## loop model structures
  for(j in 1:nMod) { 
    test_mod<-form<-test_forms[[j]]
    term_list<-strsplit(as.character(test_mod)," ")[[3]]
    term_list<-term_list[term_list!="+"]
    ##------------------------------## select optimal cut for dummy if in model
    if("MeanTemp:MeanAgeDummy" %in% term_list) {
      trainmod<-threshold.gls(form=form,newdata=traindata)
      cut<-trainmod$thd.select
      testdata$MeanAgeDummy[testdata$MeanAge<quantile(traindata$MeanAge,cut)]<-"low" 
      testdata$MeanAgeDummy[testdata$MeanAge>=quantile(traindata$MeanAge,cut)]<-"high"
      nlevels<-length(unique(testdata$MeanAgeDummy))
      if(nlevels==1) {
        cnt<-cnt+1
        next
      } ## skip test datasets that contain only one level of the categorical
      predicted<-predict(trainmod$mod.select,newdata=testdata,se.fit=T,type="response")
    } else {
      trainmod<-gls(test_mod,data=traindata,correlation=corAR1(),method="ML")
      predicted<-predict(trainmod,newdata=testdata,se.fit=T,type="response")
    }
    ##--------------------------------------------## root mean squared error
    pred<-as.numeric(predicted$fit)
    true<-as.numeric(testdata$lnRpSSB)
    RMSE[i,j]<-sqrt(sum((pred-true)^2)/nY)
    ##------------------------------------------------------## save formulas
    if("MeanTemp:MeanAgeDummy" %in% term_list) { 
      prev_mod<-trainmod$mod.select
    } else {
      prev_mod<-trainmod 
    }
    save_terms<-attributes(prev_mod)$namBetaFull[-1]
    allterms<-paste(save_terms,collapse="+") 
    new_mod<-formula(paste("~",allterms,sep=""))
    formulas[[j]]<-new_mod
  } ## end loop over models (j)
  cnt<-cnt+1
} ## end stochastic loop (i)
end<-Sys.time()
print(end-start)

##=============================## choose model with best predictive ability
med_RMSE<-apply(RMSE,2,function(x) median(x,na.rm=T))
##------------------------------------------## reorder based on median RMSE
indices<-order(med_RMSE)
RMSE<-RMSE[,indices]
models<-as.character(formulas)[indices]
##------------------------------------------------------------## best model
best_mod_terms<-models[1]
best_mod<-formula(paste("lnRpSSB",best_mod_terms,sep=""))
##--------------------------------------------------------## model names
full_model_names<-models
models<-gsub("Cannibals","C",models)
models<-gsub("MeanAgeDummylow","A",models)
models<-gsub("MeanAge","A",models)
models<-gsub("MeanTemp","T",models)
models<-gsub("TotSpBio","S",models)
models<-gsub(" ","",models)

##=============================================================## plot RMSE
pdf("plot_model_cross_validation.pdf",height=3+nMod*0.1,width=4) 
par(mar=c(4,6,1,1),mgp=c(2,0.5,0),cex.lab=1,tcl=-0.3,las=1)#mar=c(4,25,1,1)
cols<-c("gray60",rep("gray90",nMod-1))
boxplot(RMSE,axes=F,horizontal=T,range=1,outline=F,col=cols,lty=1,lwd=0.5, xlim=c(0.5,nMod+0.5),ylim=c(0.18,0.46),pars=list(boxwex=0.5),pch=16,cex=0.25, ylab="",xlab="RMSE");box()
abline(v=median(RMSE[,1]),lty=1,lwd=0.1)
axis(1,at=seq(0,1,0.1),labels=T,line=0)
models<-gsub("MeanAgeDummylow","MeanAgeD",models)
axis(2,at=seq(1,nMod,1),labels=models,line=0,cex.axis=0.8)
dev.off() ## none of the simpler models has a lower median RMSE 

##========================================================================##
##===========================## alternative metrics of spawner age structure
##========================================================================##
best_mod<-formula(lnRpSSB~Cannibals+MeanAge+MeanTemp+TotSpBio)
mod<-gls(best_mod,data=data,correlation=corAR1(),method="ML")
mod_alt_1<-update(mod, .~. -MeanAge+AgeDiv)
mod_alt_2<-update(mod, .~. -MeanAge+MeanWeight)
aics_alt<-AICc(mod,mod_alt_1,mod_alt_2)
alt_names<-c("MeanAge","AgeDiv","MeanWeight")
rownames(aics_alt)<-alt_names

##==============## cross-validation for alternative metrics of age structure
test_forms<-list(formula(lnRpSSB~MeanTemp+Cannibals+TotSpBio+MeanAge), formula(lnRpSSB~MeanTemp+Cannibals+TotSpBio+MeanWeight), formula(lnRpSSB~MeanTemp+Cannibals+TotSpBio+AgeDiv))
nMod<-length(test_forms)
nS<-1e3 ## number of runs 
seeds<-seq(nS) 
RMSE<-array(dim=c(nS,nMod))
formulas<-list()
cnt<-1
for(i in 1:nS) {
  set.seed(seeds[cnt])
  train<-sort(sample(seq(nY),round(0.75*nY),replace=F)) ## use 75% to train
  traindata<-data[train,] 
  test<-seq(nY)[seq(nY) %!in% train]
  testdata<-data[test,]
  ##------------------------------------------------## loop model structures
  for(j in 1:nMod) { 
    test_mod<-form<-test_forms[[j]]
    term_list<-strsplit(as.character(test_mod)," ")[[3]]
    term_list<-term_list[term_list!="+"]
    trainmod<-gls(test_mod,data=traindata,correlation=corAR1(),method="ML")
    predicted<-predict(trainmod,newdata=testdata,se.fit=T,type="response")
    ##--------------------------------------------## root mean squared error
    pred<-as.numeric(predicted$fit)
    true<-as.numeric(testdata$lnRpSSB)
    RMSE[i,j]<-sqrt(sum((pred-true)^2)/nY)
    ##------------------------------------------------------## save formulas
    prev_mod<-trainmod
    save_terms<-attributes(prev_mod)$namBetaFull[-1]
    allterms<-paste(save_terms,collapse="+") 
    new_mod<-formula(paste("~",allterms,sep=""))
    formulas[[j]]<-new_mod
  } ## end loop over models (j)
  cnt<-cnt+1
} ## end stochastic loop (i)
med_RMSE_alt<-apply(RMSE,2,function(x) median(x,na.rm=T))
names(med_RMSE_alt)<-alt_names

##========================================================================##
##============================================================## final model
##========================================================================##
mod<-gls(best_mod,data=data,correlation=corAR1(),method="REML")
phi_est<-coef(mod$modelStruct$corStruct,unconstrained=F)
sigma_est<-summary(mod)$sigma 
residuals<-residuals(mod,type="normalized")
fitted<-fitted(mod)
car::vif(mod) ## test for multicollinearity
out_lm<-summary(mod) 
coeffs_lm<-out_lm$coefficients
CIs<-intervals(mod)

##============================================================## model plots
pdf("plot_model_partial_effects.pdf",width=6,height=6)
layout(matrix(c(1:4),nrow=2,byrow=T))
par(mar=c(4,4,0.5,0.5),oma=c(0,0,1,0.5),mgp=c(2.2,0.5,0),tcl=-0.3,cex.axis=1,cex.lab=1.1)
ylab<-"Partial effect on ln(recruits/SSB)"
cexll<-1.2 ## letter size
ll<-c(-0.36,-0.16) ## inset for letters
##-------------------------------------------------------## spawner biomass
xtrans<-function(x) { (x*data_sds[names(data_sds)=="TotSpBio"]+data_means[names(data_means)=="TotSpBio"])*1e-6 }
visreg(mod,xvar="TotSpBio",xlab="SSB (million tons)",partial=T, ylab=ylab,scale="response",points.par=list(cex=0.8,pch=21,bg=alpha(1,0.25)), fill.par=list(col=alpha(1,0.1)),line.par=list(lwd=2,col=1),xtrans=xtrans)
legend("topleft","a",text.font=2,cex=cexll,bty="n",inset=ll,xpd=NA)
##--------------------------------------------------------------## mean age
xtrans<-function(x) { x*data_sds[names(data_sds)=="MeanAge"]+data_means[names(data_means)=="MeanAge"] }
visreg(mod,xvar="MeanAge",xlab="Mean age of spawners (yrs)",partial=T, ylab=ylab,scale="response",points.par=list(cex=0.8,pch=21,bg=alpha(1,0.25)), fill.par=list(col=alpha(1,0.1)),line.par=list(lwd=2,col=1),xtrans=xtrans)
legend("topleft","b",text.font=2,cex=cexll,bty="n",inset=ll,xpd=NA)
##-------------------------------------------------------------## cannibals
xtrans<-function(x) { (x*data_sds[names(data_sds)=="Cannibals"]+data_means[names(data_means)=="Cannibals"])*1e-6 }
visreg(mod,xvar="Cannibals",xlab="Biomass ages 3-6 (million kg)",partial=T, ylab=ylab,scale="response",points.par=list(cex=0.8,pch=21,bg=alpha(1,0.25)), fill.par=list(col=alpha(1,0.1)),line.par=list(lwd=2,col=1),xtrans=xtrans)
legend("topleft","c",text.font=2,cex=cexll,bty="n",inset=ll,xpd=NA)
##-----------------------------------------------------------## temperature
xtrans<-function(x) { x*data_sds[names(data_sds)=="MeanTemp"]+data_means[names(data_means)=="MeanTemp"] }
visreg(mod,xvar="MeanTemp",xlab="Mean temperature (°C)",partial=T, ylab=ylab,scale="response",points.par=list(cex=0.8,pch=21,bg=alpha(1,0.25)), fill.par=list(col=alpha(1,0.1)),line.par=list(lwd=2,col=1),xtrans=xtrans)
legend("topleft","d",text.font=2,cex=cexll,bty="n",inset=ll,xpd=NA)
##------------------------------------------------------------------## save
dev.off()

##========================================================## residual checks
pdf("plot_model_residuals.pdf",width=8,height=4)
par(mar=c(4,4,1,1),mgp=c(2.5,0.75,0),cex.axis=1,cex.lab=1.2,pch=16)
layout(matrix(c(1:8),nrow=2,byrow=T))
ylab="Residuals"
plot(fitted,residuals,xlab="Fitted values",ylab=ylab)
abline(h=0,lty=3) 
qqnorm(residuals,main="")
qqline(residuals)
pacf(residuals,lag=10,plot=T,lwd=2,ylim=c(-0.4,0.4))$acf
plot(data$Year,residuals,xlab="Year",ylab=ylab)
abline(h=0,lty=3) 
plot(data$TotSpBio,residuals,xlab="Spawning stock biomass",ylab=ylab)
abline(h=0,lty=3)
plot(data$MeanAge,residuals,xlab="Mean age of spawners",ylab=ylab)
abline(h=0,lty=3)
plot(data$MeanTemp,residuals,xlab="Mean temperature",ylab=ylab)
abline(h=0,lty=3)
plot(data$Cannibals,residuals,xlab="Biomass ages 3-6",ylab=ylab)
abline(h=0,lty=3)
dev.off()

##========================================================================##
##===========================================## model predicted and observed
##========================================================================##
predicted<-predict(mod,newdata=list(TotSpBio=data$TotSpBio,MeanTemp=data$MeanTemp,MeanAge=data$MeanAge,Cannibals=data$Cannibals), se.fit=T,interval="confidence",level=0.95,type="response")
pred<-as.vector(predicted$fit)
upper<-as.vector(predicted$fit+2*predicted$se.fit)
lower<-as.vector(predicted$fit-2*predicted$se.fit)

##===================================================================## plot
pdf("plot_lnRpSSB_prediction_and_observed.pdf",width=5,height=3.5)
par(mar=c(3.5,3.5,0.5,0.5),mgp=c(2,0.5,0),cex.axis=0.9,cex.lab=1.2,tcl=-0.3)
xlim<-c(min(data$Year)+1,max(data$Year)+3)
plot(NA,NA,xlab="Year",ylab="ln(recruits/SSB)",xlim=xlim,ylim=c(-3.1,2.4))
##-----------------------------------------------------------## predictions
X.Vec<-c(data$Year,tail(data$Year,1),rev(data$Year),data$Year[1])
Y.Vec<-c(lower,tail(upper,1),rev(upper),lower[1])
polygon(X.Vec,Y.Vec,col="grey90",border=NA) 
lines(data$Year,pred,lwd=1.5,lty=1,col="grey50")
lines(data$Year,upper,lwd=0.5,lty=1,col="grey50")
lines(data$Year,lower,lwd=0.5,lty=1,col="grey50")
##-------------------------------------------------## observed time-series
lines(data$Year,data_other$lnRpSSB,type="o",pch=16,cex=0.8,lwd=1,col=1)
dev.off()

##========================================================================##
##===================================## predicted stock-recruit relationship
##========================================================================##

##---------------------------------------------## spawner biomass scenarios
set_min<-0;set_max<-2.75e6;nS<-1e3
min_scaled<-(set_min-data_means[names(data_means)=="TotSpBio"])/data_sds[names(data_sds)=="TotSpBio"]
max_scaled<-(set_max-data_means[names(data_means)=="TotSpBio"])/data_sds[names(data_sds)=="TotSpBio"]
spawner_scen<-seq(min_scaled,max_scaled,length=nS)
spawner_scen_true<-spawner_scen*data_sds[names(data_sds)=="TotSpBio"]+data_means[names(data_means)=="TotSpBio"]
##-----------------------------------------------------## mean age scenarios
set_min<-round(min(alldata$MeanAge,na.rm=T)) ## 7 years
set_med<-round(median(alldata$MeanAge,na.rm=T)) ## 9 years
set_max<-round(max(alldata$MeanAge,na.rm=T)) ## 11 years
meanage_scen_true<-c(set_min,set_med,set_max)
nA<-length(meanage_scen_true)
##------------------------------------------## scaled values for predictions
min_scen<-(set_min-data_means[names(data_means)=="MeanAge"])/data_sds[names(data_sds)=="MeanAge"]
med_scen<-(set_med-data_means[names(data_means)=="MeanAge"])/data_sds[names(data_sds)=="MeanAge"]
max_scen<-(set_max-data_means[names(data_means)=="MeanAge"])/data_sds[names(data_sds)=="MeanAge"]
meanage_scen<-c(min_scen,med_scen,max_scen)
##---------------------------------------------------------## new data frame
## use range of SSB and median values for other covariates
newdat<-expand.grid(TotSpBio=spawner_scen,MeanTemp=median(data$MeanTemp), MeanAge=meanage_scen,Cannibals=median(data$Cannibals))
##-----------------------------------------------------------## unscaled SSB
newdat$unscaled<-TotSpBio_unscaled<-newdat$TotSpBio*data_sds[names(data_sds)=="TotSpBio"]+data_means[names(data_means)=="TotSpBio"]
##------------------------------------------------------------## predictions
predicted<-predict(mod,newdata=newdat,se.fit=T,se=T,interval="confidence", level=0.95,type="response")
newdat$pred<-data.frame(predicted$fit)[,1]
newdat$se<-data.frame(predicted$se.fit)[,1]
newdat$lower<-newdat$pred-2*newdat$se
newdat$upper<-newdat$pred+2*newdat$se
newdat$RpSSB_fit<-exp(newdat$pred) 	
newdat$RpSSB_low<-exp(newdat$lower) 	
newdat$RpSSB_upp<-exp(newdat$upper) 
newdat$recruits_fit<-newdat$RpSSB_fit*newdat$unscaled 
newdat$recruits_low<-newdat$RpSSB_low*newdat$unscaled 
newdat$recruits_upp<-newdat$RpSSB_upp*newdat$unscaled 

##===================================================================## plot
pdf("plot_stock_recruit_relationship_by_age_CIs.pdf",width=5,height=4)
par(mar=c(3.5,3.5,0.5,0.5),mgp=c(2,0.5,0),cex.axis=0.9,cex.lab=1,tcl=-0.3,xaxs="i",yaxs="i")
fac<-1e6
##------------------------------------------------------------## data points
ylim<-c(0,1.4*max(data_other$Recruits/fac));ylim<-c(0,3.6)
plot(NA,NA,xlab="Spawning stock biomass (million tons)",ylab="Recruits (billions)", xlim=c(0,1.05*max(data_unscaled$TotSpBio/fac)),ylim=ylim)
pcols<-rep(1,dim(data_other)[1])
points(data_unscaled$TotSpBio/fac,data_other$Recruits/fac,bg=alpha(pcols,0.5),pch=21,lwd=0.5,cex=0.9)
##------------------------------------------------------## model predictions
cols<-c("goldenrod1","forestgreen","darkorchid3")
ys<-array(NA,dim=c(nS,nA))
for(i in 1:nA){
  x<-TotSpBio_unscaled[((i-1)*nS+1):(i*nS)]/fac
  y1<-newdat$recruits_low[((i-1)*nS+1):(i*nS)]/fac
  y2<-ys[,i]<-newdat$recruits_fit[((i-1)*nS+1):(i*nS)]/fac
  y3<-newdat$recruits_upp[((i-1)*nS+1):(i*nS)]/fac
  lines(x,y1,col=cols[i],lwd=0.1)
  lines(x,y2,col=cols[i],lwd=1.5)
  lines(x,y3,col=cols[i],lwd=0.1)
  X.Vec<-c(x,tail(x,1),rev(x),x[1])
  Y.Vec<-c(y1,tail(y3,1),rev(y3),y1[1])
  polygon(X.Vec,Y.Vec,col=alpha(cols[i],0.1),border=NA)
}
legend("topright",c("    ","    ","    "),title="",lwd=6,cex=0.8,col=alpha(rev(cols),0.1),bty="n",inset=c(0,0.06),seg.len=2)
legend("topright",rev(paste0(meanage_scen_true)),title="spawner\nmean age", lwd=1.5,cex=0.8,col=rev(cols),bty="n",inset=c(0,0.06),seg.len=1.8)
dev.off()

##========================================================================##
##========================================================================##
##========================================================================##