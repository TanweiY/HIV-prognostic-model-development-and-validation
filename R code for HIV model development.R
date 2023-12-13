# load packages
library(mice)
library(riskRegression)
library(glmnet)
library(rlist)
library(ggplot2)
library(survival)
library(rms)
library(survminer)
library(concreg)
################################################# 1. Model development #################################################
################################################# 1.1 Multiple imputation #################################################
# Derivation cohort
# Prepare data for imputation

# Prepare imputation for baseline data
fixvari<-subset(gzbase, select = c("id", "gender", "age", "maritalstatuscat", "HIVdiagdate",
                                   "infectionroute", "Bweight", "Bheight", "Btuberculosis", "HCV",
                                   "BWHO",  "ARTinitidate", "death", "time"))

# Calculate cumulative hazard
HT1 <- summary(survival::survfit(Surv(time,death)~1,data=fixvari))
fixvari$haz_os <-  approx(c(0,HT1$time),-log(c(1,HT1$surv)),xout=fixvari$time,method="constant",f=0,rule=2)$y

fixvari$time<-NULL

# Prepare imputation for time updated data

# Extract follow-up variables
followtimeup<-subset(trainmerge, select = c("id", "FCD4", "FCD8", "Fviralload", "FWBC", "FPlatelet",
                                            "Fhaemoglobin", "FSCR", "FTG", "FTch", "FGLU", "FAST", "FALT",
                                            "FT.BIL", "time"))

# Extract baseline variables
basetimeup<-subset(gzbase, select = c("id",  "BCD4", "BCD8", "Bviralload", "WBC", "Platelet", "haemoglobin",
                                      "SCR", "TG", "Tch", "GLU", "AST", "ALT", "T.BIL"))


# Set baseline time = 0, all alive
basetimeup$time<-0
basetimeup$status<-0

# Rename baseline variable names
colnames(basetimeup)[2]<-'FCD4'
colnames(basetimeup)[3]<-'FCD8'
colnames(basetimeup)[4]<-'Fviralload'
colnames(basetimeup)[5]<-'FWBC'
colnames(basetimeup)[6]<-'FPlatelet'
colnames(basetimeup)[7]<-'Fhaemoglobin'
colnames(basetimeup)[8]<-'FSCR'
colnames(basetimeup)[9]<-'FTG'
colnames(basetimeup)[10]<-'FTch'
colnames(basetimeup)[11]<-'FGLU'
colnames(basetimeup)[12]<-'FAST'
colnames(basetimeup)[13]<-'FALT'
colnames(basetimeup)[14]<-'FT.BIL'

# Merge baseline and follow up variable together
timeup<-bind_rows(basetimeup, followtimeup)
timeup<-arrange(timeup, id, time)

# Merge all data for imputation together
gz.long<-merge(fixvari, timeup, by='id')

# Log transform laboratory data
gz.long<-within.data.frame(gz.long,{
  logCD4<-log(FCD4)
  logCD8<-log(FCD8)
  logviralload<-log(Fviralload)
  logWBC<-log(FWBC)
  logPlatelet<-log(FPlatelet)
  loghaemoglobin<-log(Fhaemoglobin)
  logSCR<-log(FSCR)
  logTG<-log(FTG)
  logTch<-log(FTch)
  logGLU<-log(FGLU)
  logAST<-log(FAST)
  logALT<-log(FALT)
  logT.BIL<-log(FT.BIL)
})

gz.long<-subset(gz.long, select=c("id", "gender", "age", "maritalstatuscat", "HIVdiagdate", "infectionroute",
                                  "Bweight", "Bheight", "Btuberculosis", "HCV", "BWHO", "ARTinitidate",
                                  "death", "haz_os",  "time",  "logT.BIL", "logALT", "logAST", "logGLU",
                                  "logTch", "logTG", "logSCR", "loghaemoglobin", "logPlatelet",
                                  "logWBC", "logviralload", "logCD8", "logCD4"))

# Add outcome to imputation
gz.long$death<- factor(gz.long$death)

# Imputation starts from here

data.impu <- gz.long
data.impu$int <- rep(1,nrow(data.impu))

# See all the default settings for imputation
ini2 <- mice(data.impu, maxit=0)
summary(ini2)

# Check predictor structure
pred <- quickpred(data.impu, exclude= c("id","time"))
pred

# Check imputation method
meth<-ini2$meth
meth

# Death and time will not be imputed
pred["time", ] <- 0
pred["death", ] <- 0
meth[c("time","death")] <- ""


# Change longitudinal variables to 2 level format, as this is a longitudinal variable
meth[c("logT.BIL", "logALT", "logAST",
       "logGLU", "logTch", "logTG", "logSCR", "loghaemoglobin", "logPlatelet",
       "logWBC", "logviralload", "logCD8", "logCD4")] <- "2l.norm"

# Tell R for which variables are longitudinal variables, which is the intercept, id variable and time variable

pred["logT.BIL", c("id","time","int")] <- c(-2,2,2)
pred["logALT", c("id","time","int")] <- c(-2,2,2)
pred["logAST", c("id","time","int")] <- c(-2,2,2)
pred["logGLU", c("id","time","int")] <- c(-2,2,2)
pred["logTch", c("id","time","int")] <- c(-2,2,2)
pred["logTG", c("id","time","int")] <- c(-2,2,2)
pred["logSCR", c("id","time","int")] <- c(-2,2,2)
pred["loghaemoglobin", c("id","time","int")] <- c(-2,2,2)
pred["logPlatelet", c("id","time","int")] <- c(-2,2,2)
pred["logWBC", c("id","time","int")] <- c(-2,2,2)
pred["logviralload", c("id","time","int")] <- c(-2,2,2)
pred["logCD8", c("id","time","int")] <- c(-2,2,2)
pred["logCD4", c("id","time","int")] <- c(-2,2,2)

# Time updated variables are seen as longitudinal
meth

# Remove duplicates
data.impute.noneduplicate<-subset(data.impu,duplicated(data.impu)=="FALSE")

# Set number of imputation
K <- 10

# Impute the missing values with mice
imputation_10 <- mice(data.impute.noneduplicate, maxit=5, m=K, seed = 1234, pred=pred, meth=meth, print=TRUE)

################### 1.2 Inspect the linear relationship between continuous variables and hazard##################################
# Prepare the imputed data
n_impu <- 10
data_imputated <- vector(K,mode="list")
time<-subset(gzbase, select = c('id', 'time'))
for (i in 1:n_impu) {
  data_impu[[i]] <- mice::complete(gzimputation_10, i)

  ###death
  data_impu[[i]]$death<-as.numeric(data_impu[[i]]$death)-1

  ##BMI
  data_impu[[i]]$BMI<-(data_impu[[i]]$Bweight/((data_impu[[i]]$Bheight)*(data_impu[[i]]$Bheight)))

  ##t=0

  base_impu[[i]]<-data_impu[[i]][data_impu[[i]]$time==0, ]

  base_impu[[i]]$time<-NULL

  base_impu[[i]]<-merge(base_impu[[i]], time, by = 'id')

}

# Stack the ten imputed datasets together
gzbasestack<-base_impu[[1]]
for (i in 1:9) {
  gzbasestack<-bind_rows(gzbasestack, base_impu[[i+1]])
}

# Inspect the linear relationship
gzb<-gzbasestack
dd <- datadist(gzb)
options(datadist="dd")

fit.BMI <- cph(Surv(time, death)~ rcs(BMI,3),data=gzb) ###10-30
surv2<-Predict(fit.BMI, BMI, ref.zero=TRUE)
summary(surv2$yhat)
summary(surv$age)
his2<-hist(gzb$BMI, freq = FALSE)
BMIx1 <- surv2$BMI
BMIy1 <- surv2$yhat
BMIx2 <- his2$mids
BMIy2 <- his2$density
twoord.plot(lx = BMIx1, ly = BMIy1, rx = BMIx2, ry = BMIy2, type=c("l","bar"), lcol
            = "#d00000", rcol = "NA", xlab = "Body mass index(kg/m2)",
            ylab = "log Relative Hazard", rylab = "Density",  halfwidth=2.5, lylim = c(-1,2), rylim
            = c(0, 0.2),lwd=2)

########### same procedure applied for other continous variables including age, CD4 cell counts,
#CD8 cell counts, white blood cell count, platelet, haemoglobin, creatinine,
#triglyceride, total cholesterol, plasma glucose, aspartate transaminase, alanine aminotransferase,
# and total bilirubin

############################################### 1.3 Further data transformation  #############################
# for relationships with U-shaped pattern, the turning points were found for further transformations.
# We retain the transformation that generate linear relationship (absolute distance from the turning point)

# logT.BIL
# Find the turning point
fit.logT.BIL <- cph(Surv(time, death)~ rcs(logT.BIL,3),data=gzb)
surv3<-Predict(fit.logT.BIL, logT.BIL,ref.zero=TRUE)

pT.BILori<-ggplot(Predict(fit.logT.BIL, logT.BIL,ref.zero=TRUE),
                  conf="fill", colfill='white')+
  scale_x_continuous(name='logT.BIL')

hist3<-ggplot(gzb, aes(x=logT.BIL, y=..density..), break.x.by=1,) +
  geom_histogram(binwidth=0.1, color="skyblue", fill="white")+xlab("Total bilirubin(μmol/L)")+
  theme(axis.text = element_text(size = 15), axis.title=element_text(size=15),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color='black'))

summary(surv3)
b<-surv3[surv3$yhat == min(surv3$yhat), ] ##### turning point: logT.BIL=2.4,

# Transformation 1:  absolute distance from the turning point:
gzb$absT.BIL<-abs(gzb$logT.BIL-2.4)
dd <- datadist(gzb)
options(datadist="dd")
fit.absT.BIL <- cph(Surv(time, death)~ rcs(absT.BIL,3),data=gzb)
pabsT.BIL<-ggplot(Predict(fit.absT.BIL, absT.BIL,ref.zero=TRUE),
                  conf="fill", colfill='white')+
  scale_x_continuous(name='|logT.BIL-2.4|')

# Transformation 2: square distance from the turning point:
gzb$T.BIL2<-(gzb$logT.BIL-2.4)^2
summary(gzb$T.BIL2)
dd <- datadist(gzb)
options(datadist="dd")

fit.T.BIL2 <- cph(Surv(time, death)~ rcs(T.BIL2,3),data=gzb)

pT.BIL2<-ggplot(Predict(fit.T.BIL2, T.BIL2,ref.zero=TRUE),
                conf="fill", colfill='white')+
  scale_x_continuous(name='(logT.BIL-2.4)^2')

# same procedure applied for other variables including white blood cell count, platelet, creatinine,
# triglyceride, plasma glucose, alanine aminotransferase, and total cholesterol

###################################################### 1.4 Predictor selection ####################################################################
# Data preparation
n_impu <- 10
data_impu <- vector(n_impu,mode="list")
time<-subset(gzbase, select = c('id', 'time'))

for (i in 1:n_impu) {
  data_impu[[i]] <- mice::complete(gzimputation_10, i)

  ###death
  data_impu[[i]]$death<-as.numeric(data_impu[[i]]$death)-1

  ##BMI
  data_impu[[i]]$BMI<-(data_impu[[i]]$Bweight/((data_impu[[i]]$Bheight)*(data_impu[[i]]$Bheight)))

  ##variable transformation

  data_impu[[i]]$absT.BIL<-abs( data_impu[[i]]$logT.BIL-2.4)

  data_impu[[i]]$absGLU<-abs(data_impu[[i]]$logGLU-1.66)

  data_impu[[i]]$absTch<-abs(data_impu[[i]]$logTch-1.48)

  data_impu[[i]]$absSCR<-abs(data_impu[[i]]$logSCR-4.29)

  data_impu[[i]]$absPlatelet<-abs(data_impu[[i]]$logPlatelet-5.27)

  data_impu[[i]]$absWBC<-abs(data_impu[[i]]$logWBC-1.67)

  data_impu[[i]]$absALT<-abs(data_impu[[i]]$logALT-3.03)

  data_impu[[i]]$absTG<-abs(data_impu[[i]]$logTG-0.62)

  ### in order to avoid group lasso in predictor selection, we convert multi-categorical variables to binary variable
  data_impu[[i]]$WHO2<-as.factor(ifelse(data_impu[[i]]$BWHO=='I'|data_impu[[i]]$BWHO=='II', 'I/II', 'III/IV'))
  data_impu[[i]]$infectionroute2<-as.factor(ifelse(data_impu[[i]]$infectionroute =='IVDU', 'IVDU', 'Other'))

  ##select t=time, merge end time

  data_impu[[i]]<-data_impu[[i]][data_impu[[i]]$time==0, ]

  data_impu[[i]]$time<-NULL

  data_impu[[i]]<-merge(data_impu[[i]], time, by = 'id')

}

# Selection starts from here
set.seed (1234)
CandidateVariables <- c("gender", "age", "maritalstatuscat",  "infectionroute2",
                        "Btuberculosis", "HCV", "WHO2",  "BMI",
                        "logAST","absT.BIL",
                        "loghaemoglobin",
                        "logCD8", "logCD4",
                         "absGLU", "absTch", "absSCR", "absPlatelet", "absWBC",
                        "absALT", "absTG")

coef_1s<-vector(n_impu,mode="list")

for (i in 1:n_impu) {
 cv.model <- cv.glmnet(model.matrix(~., data_impu[[i]][CandidateVariables]),
                        Surv(data_impu[[i]]$time, data_impu[[i]]$death),
                        family="cox", nlambda=50, alpha=1, standardize=TRUE)
 coef_1s[[i]]<-as.data.frame(as.matrix(coef(cv.model, s=cv.model$lambda.1se) ))

}

# Merge predictor selection results from 10 datasets
coef_1sm<-coef_1s[[1]]

for (i in 1:9) {

  coef_1sm<-bind_cols(coef_1sm, coef_1s[[i+1]])

}

coef1<-coef_1sm[[1]]

colnames(coef_1sm)<-paste0("lambda.1se", 1:10)

coefall<- cbind(names = rownames(coef_1sm), coef_1sm)
write_xlsx(coefall, path ='variselection2.xlsx')

# We selected predictors that were consistently retained in the ten models calculated from the ten imputed datasets

###################################### 1.5 Ajustment for optimism ######################################
SelectedVariables <- c( "age",  "infectionroute",
                        "Btuberculosis", "HCV",
                        "BMI",
                        "logAST",
                        "loghaemoglobin",
                        "logCD4",
                         "absGLU", "absPlatelet")

fml<-as.formula(paste0('Surv(time, death)~', paste0(SelectedVariables, collapse = '+')))

# Extract optimism adjustment and C statistics
name<-c("age", "infectionrouteheterosex", "infectionrouteIVDU",
        "infectionrouteOther", "Btuberculosistuberculosis", "HCV",
        "BMI",
        "logAST",
        "loghaemoglobin",
        "logCD4",
         "absGLU",  "absPlatelet")

optim_ad<-vector(n_impu,mode="list")

for (i in 1:n_impu) {
  final_model <- cph(fml, data=data_impu[[1]],x=TRUE,y=TRUE)
  optim_ad[i]<-validate(final_model, B=100)[3,3] # optimism adjustment
}

op<-data.frame(list.cbind(optim_ad))
range(op[1,])
optim_ad[i]

# Extract coefficient from each model
final_slope<-vector(n_impu,mode="list")

for (i in 1:n_impu) {
  final_model <- coxph(fml, data=data_impu[[i]])
  final_slope[i] <- as.data.frame(final_model$coefficients)
}

final_slope<-list.cbind(final_slope)

finalcoef<-data.frame(final_slope, row.names = name)

# Coefficients multiple optimism

adjustcof<-vector(n_impu,mode="list")
for (i in 1:n_impu) {

  adjustcof[[i]]<-as.numeric(optim_ad[[i]])*finalcoef[,i]
}

adjustcof<-list.cbind(adjustcof)

adjustcof<-data.frame(adjustcof, row.names = name)

adjustcof$mean<-row_means(adjustcof)

adjustcof$variable<-rownames(adjustcof)

c(adjustcof$variable, adjustcof$mean)

# Calculate Final Linear predictor
for (i in 1:n_impu) {
  data_impu[[i]]$score_infecroute <- ifelse(data_impu[[i]]$infectionroute=="IVDU",1.07363195299915,
                                            ifelse(data_impu[[i]]$infectionroute=="heterosex", 0.537387987475878,
                                                   ifelse(data_impu[[i]]$infectionroute=="Other",	0.928995957706017,0)))
  data_impu[[i]]$score_tuberculosis<-ifelse(data_impu[[i]]$Btuberculosis=="tuberculosis",0.454994655378687,0)
  data_impu[[i]]$score_HCV<-ifelse(data_impu[[i]]$HCV=="HCV",0.61112813797825,0)

  data_impu[[i]]$lp_nomogram<-data_impu[[i]]$score_infecroute + data_impu[[i]]$score_tuberculosis +
    data_impu[[i]]$score_HCV+0.0534869544053514*data_impu[[i]]$age-0.0697910578563831*data_impu[[i]]$BMI+
    0.237997613173825*data_impu[[i]]$logAST-0.662203221121992*data_impu[[i]]$loghaemoglobin-0.199524528915675*data_impu[[i]]$logCD4+
   0.777760863142926*data_impu[[i]]$absGLU+0.339954731518507*data_impu[[i]]$absPlatelet
}

###################################### 1.6 Internal validation ######################################

# Calculate originial and adjusted C statistics and their 95%CI
C_statistics<-0
C_statistics_var<-0
C_statistics_opm<-0
C_statistics_opm_var<-0
c_harrell_apparent <- 0
for (i in 1:n_impu) {
  final_model <- cph(fml, data=data_impu[[i]],x=TRUE,y=TRUE)
  C_statistics[i]<-(validate(final_model, B=100)[1,1]+1)/2 # original C statistics
  #C_statistics_var[i] <- vardiffC(data_impu[[i]]$time, data_impu[[i]]$death, data_impu[[i]]$lp_nomogram, data_impu[[i]]$lp_nomogram)$est.varCxy
  C_statistics_opm[i] <-(validate(final_model, B=100)[1,5]+1)/2 # adjusted C statistics
}

# Original C statistics and its 95% CI
C_statistics_final<-mean(C_statistics)
C_statistics_var_final<-mean(C_statistics_var) + (1+1/n_impu)*var(C_statistics_var)

C_statistics_final-qnorm(0.975)*C_statistics_var_final^0.5 # lower limit of 95%
C_statistics_final+qnorm(0.975)*C_statistics_var_final^0.5 # upper limit of 95%

# Adjusted C statistics
C_statistics_opm_final<-mean(C_statistics_opm)

# Internal calibration accuracy
# compare predicted versus observed survival probability
for (i in 1:n_impu) {
  # create risk groups
  data_impu[[i]]$group3<-cut(data_impu[[i]]$lp_nomogram, quantile(data_impu[[i]]$lp_nomogram, probs=c(0,(1/3),(2/3),1)), right=FALSE, labels=c(1:3))
  data_impu[[i]]$group3[data_impu[[i]]$lp_nomogram==max(data_impu[[i]]$lp_nomogram)] <- 3

}

# Calculate baseline risk
baseriskav<-vector(n_impu,mode="list")
basehazardtime<-vector(n_impu,mode="list")

for (i in 1:n_impu) {
  final_model <- cph(Surv(time, death)~lp_nomogram, data=data_impu[[i]],init=1, iter.max=0,x=TRUE,y=TRUE)
  basehazard<-basehaz(final_model, centered=FALSE)
  baseriskav[i]<-as.data.frame(basehazard$hazard)
  basehazardtime[i]<-as.data.frame(basehazard$time)
}

#Baseline hazard average

baserisk<-row_means(list.cbind(baseriskav))
basehazardtime<-row_means(list.cbind(basehazardtime))
baserisk<-data.frame(baserisk, basehazardtime)
baserisk$'Baseline Survival'<-exp(-baserisk$baserisk)
colnames(baserisk)[2]<-'Year'

gzbasestack<-data_impu[[1]]

for (i in 1:9) {

  gzbasestack<-bind_rows(gzbasestack, data_impu[[i+1]])

}
summary(gzbasestack)
calipre3<-subset(gzbasestack, select = c("lp_nomogram", "group3"))

for (i in 1:5555) {
  calipre3[[paste0(baserisk[i,2])]] <-exp(-baserisk[i,1])^(exp(calipre3$lp_nomogram))
}


predictc3<-calipre3 %>% group_by(group3) %>% summarise_all(funs(mean))
predictc3$lp_nomogram<-NULL
predictcl3 <- melt(setDT(predictc3), id.vars = c("group3"), variable.name = "time")
colnames(predictcl3)[3]<-'predictedsurv'
colnames(predictcl3)[2]<-'timepredict'
write_xlsx(predictcl3, path = "predictcl3.xlsx", col_names = T)

## Read the file 'predictcl3.xlsx' into R, with the column 'timepredict' changed from character to numeric

#Calculate observed survival

# group1
data_temp1 <- subset(gzbasestack, gzbasestack$group3==1)

fit_observe1 <- survfit(Surv(time,death) ~ 1, data=data_temp1, id=id, cluster=id)
fit_observe1$lower
survival_observe1<- data.frame(as.data.frame(fit_observe1$surv),
                               as.data.frame(fit_observe1$lower),
                               as.data.frame(fit_observe1$upper),
                               as.data.frame(fit_observe1$time))
colnames(survival_observe1)[1]<-'observesur'
colnames(survival_observe1)[2]<-'oblower'
colnames(survival_observe1)[3]<-'obupper'
colnames(survival_observe1)[4]<-'observetime'
survival_observe1$group3<-1

# group2
data_temp2 <- subset(gzbasestack, gzbasestack$group3==2)

fit_observe2 <- survfit(Surv(time,death) ~ 1, data=data_temp2, id=id, cluster=id)

survival_observe2<- data.frame(as.data.frame(fit_observe2$surv),
                               as.data.frame(fit_observe2$lower),
                               as.data.frame(fit_observe2$upper),
                               as.data.frame(fit_observe2$time))

colnames(survival_observe2)[1]<-'observesur'
colnames(survival_observe2)[2]<-'oblower'
colnames(survival_observe2)[3]<-'obupper'
colnames(survival_observe2)[4]<-'observetime'
survival_observe2$group3<-2

# group3
data_temp3 <- subset(gzbasestack, gzbasestack$group3==3)

fit_observe3 <- survfit(Surv(time,death) ~ 1, data=data_temp3, id=id, cluster=id)

survival_observe3<- data.frame(as.data.frame(fit_observe3$surv),
                               as.data.frame(fit_observe3$lower),
                               as.data.frame(fit_observe3$upper),
                               as.data.frame(fit_observe3$time))

colnames(survival_observe3)[1]<-'observesur'
colnames(survival_observe3)[2]<-'oblower'
colnames(survival_observe3)[3]<-'obupper'
colnames(survival_observe3)[4]<-'observetime'
survival_observe3$group3<-3

survival_observe3<-bind_rows(survival_observe1, survival_observe2, survival_observe3)

survival_observe3$group3<-as.factor(survival_observe3$group3)
predictcl3$group3<-as.factor(predictcl3$group3)

colnames(survival_observe3)[5]<-'Group'
colnames(predictcl3)[1]<-'Group'

# Plot calibration graph
four3<- ggplot() +
  geom_line(data=predictcl3, aes(x=timepredict, y=predictedsurv, group=Group, color=Group), linetype = "dashed", size=1)+
  geom_line(data=survival_observe3, aes(x=observetime, y=observesur, group=Group, color=Group))+
  geom_ribbon(data=survival_observe3, aes(ymin=oblower, ymax=obupper, x=observetime, y=observesur, group=Group,
                                          fill=Group), alpha = 0.15, show.legend=F)+
  scale_x_continuous(limits = c(0,10),name='Follow-up Year', )+
  scale_y_continuous(limits = c(0.5,1), name='Survival Probability')+
  theme(axis.text = element_text(size = 15), axis.title=element_text(size=15),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color='black'),
        axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
        axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm") ))+
  scale_fill_manual(values = c("#023e7d", "#d90429", "#80b918"))+
  scale_color_manual(values = c("#023e7d", "#d90429", "#80b918"))

# Calibration plot
sample_size=16481
n_impu=10
S2_nomogram <- seq(1:sample_size)
S5_nomogram <- seq(1:sample_size)
S10_nomogram <- seq(1:sample_size)

survival_predicted_2_combine <- seq(1:10)
survival_predicted_5_combine <- seq(1:10)
survival_predicted_10_combine <- seq(1:10)

survival_observed_2_combine <- seq(1:10)
survival_observed_5_combine <- seq(1:10)
survival_observed_10_combine <- seq(1:10)

survival_observed_2_var_combine <- seq(1:10)
survival_observed_5_var_combine <- seq(1:10)
survival_observed_10_var_combine <- seq(1:10)


for (i in 1:n_impu) {
  # calculate predicted survival probability
  data_impu[[i]]$S2_nomogram <- 0.868134811^exp(data_impu[[i]]$lp_nomogram)
  data_impu[[i]]$S5_nomogram <- 0.800661362^exp(data_impu[[i]]$lp_nomogram)
  data_impu[[i]]$S10_nomogram <- 0.708730609^exp(data_impu[[i]]$lp_nomogram)

  # evenly divide patients into 10 groups

  data_impu[[i]]$group10<-cut(data_impu[[i]]$lp_nomogram, quantile(data_impu[[i]]$lp_nomogram, seq(0,1,0.1)), right=FALSE, labels=c(1:10))
  data_impu[[i]]$group10[data_impu[[i]]$lp_nomogram==max(data_impu[[i]]$lp_nomogram)] <- 10

  survival_predicted_2 <- 0
  survival_predicted_5 <- 0
  survival_predicted_10 <- 0

  # 2-year predicted survival
  survival_predicted_2 <- aggregate(data_impu[[i]]$S2_nomogram, list(data_impu[[i]]$group10), mean)
  survival_predicted_2_combine <- data.frame(survival_predicted_2_combine,survival_predicted_2$x)

  # 5-year predicted survival
  survival_predicted_5 <- aggregate(data_impu[[i]]$S5_nomogram, list(data_impu[[i]]$group10), mean)
  survival_predicted_5_combine <- data.frame(survival_predicted_5_combine,survival_predicted_5$x)

  # 10-year predicted survival
  survival_predicted_10 <- aggregate(data_impu[[i]]$S10_nomogram, list(data_impu[[i]]$group10), mean)
  survival_predicted_10_combine <- data.frame(survival_predicted_10_combine,survival_predicted_10$x)

  # observed survival
  survival_observed_2 <- 0
  survival_observed_2_var <- 0

  survival_observed_5 <- 0
  survival_observed_5_var <- 0

  survival_observed_10 <- 0
  survival_observed_10_var <- 0

  for (j in 1:10) {

    data_temp <- subset(data_impu[[i]],data_impu[[i]]$group10==j)

    fit_calibration <- survfit(Surv(time,death) ~ 1, data=data_temp)

    survival_observed_2[j] <- min(fit_calibration$surv[fit_calibration$time <= 2])
    survival_observed_2_var[j] <- (tail(fit_calibration$std.err[fit_calibration$time <= 2],1))^2

    survival_observed_5[j] <- min(fit_calibration$surv[fit_calibration$time <= 5])
    survival_observed_5_var[j] <- (tail(fit_calibration$std.err[fit_calibration$time <= 5],1))^2

    survival_observed_10[j] <- min(fit_calibration$surv[fit_calibration$time <= 10])
    survival_observed_10_var[j] <- (tail(fit_calibration$std.err[fit_calibration$time <= 10],1))^2
  }


  survival_observed_2_combine <- data.frame(survival_observed_2_combine,survival_observed_2)
  survival_observed_5_combine <- data.frame(survival_observed_5_combine,survival_observed_5)
  survival_observed_10_combine <- data.frame(survival_observed_10_combine,survival_observed_10)

  survival_observed_2_var_combine <- data.frame(survival_observed_2_var_combine,survival_observed_2_var)
  survival_observed_5_var_combine <- data.frame(survival_observed_5_var_combine,survival_observed_5_var)
  survival_observed_10_var_combine <- data.frame(survival_observed_10_var_combine,survival_observed_10_var)

}

# plot calibration plot

# 2-year

survival_predicted_2_final <- exp(rowMeans(log(survival_predicted_2_combine[,-1])))
survival_observed_2_final <- exp(rowMeans(log(survival_observed_2_combine[,-1])))
survival_observed_2_var_final <- rowMeans(survival_observed_2_var_combine[,-1]) + (1+1/n_impu)*apply(survival_observed_2_var_combine[,-1], MARGIN=1, FUN=var, na.rm=TRUE)

survival_lower2_final<-exp(log(survival_observed_2_final) - qnorm(0.975)*survival_observed_2_var_final^0.5)
survival_upper2_final<-exp(log(survival_observed_2_final) + qnorm(0.975)*survival_observed_2_var_final^0.5)
survival_comparison2 <- data.frame(survival_predicted_2_final, survival_observed_2_final,
                                   survival_lower2_final,survival_upper2_final)

survival_comparison2$survival_upper2_final<-ifelse(survival_comparison2$survival_upper2_final>1, 1,survival_comparison2$survival_upper1_final)
survival_comparison2$underestimate<-(survival_comparison2$survival_observed_2_final-survival_comparison2$survival_predicted_2_final)/survival_comparison2$survival_observed_2_final

c1<-ggplot(data=survival_comparison2, aes(x=survival_predicted_2_final, y=survival_observed_2_final)) +
  geom_line(size=1, colour="dodgerblue3")+
  geom_errorbar(data=survival_comparison2, mapping=aes(x=survival_predicted_2_final, ymin=survival_lower2_final,
                                                       ymax=survival_upper2_final),
                colour="dodgerblue3", size=1,alpha=0.5, linetype=1)+
  geom_point(size=3, colour="dodgerblue3")+
  xlim(0.25,1)+
  ylim(0.25,1)+
  geom_abline(intercept = 0, slope = 1,lty=2)+
  labs(title="A. Calibration curve at 2 year",x="Predicted Survival Probability", y = "Observed Survival Probability")+
  theme( axis.text = element_text(size = 15), axis.title=element_text(size=15),
         panel.background = element_rect(fill = "white"),
         plot.title = element_text(size=15),
         axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
         axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm")),
         plot.margin = margin(1, 1, 1, 1, "cm"),
         axis.line = element_line(color='black'),
  )

# 5-year

survival_predicted_5_final <- exp(rowMeans(log(survival_predicted_5_combine[,-1])))
survival_observed_5_final <- exp(rowMeans(log(survival_observed_5_combine[,-1])))
survival_observed_5_var_final <- rowMeans(survival_observed_5_var_combine[,-1]) + (1+1/n_impu)*apply(survival_observed_5_var_combine[,-1], MARGIN=1, FUN=var, na.rm=TRUE)

survival_lower5_final<-exp(log(survival_observed_5_final) - qnorm(0.975)*survival_observed_5_var_final^0.5)
survival_upper5_final<-exp(log(survival_observed_5_final) + qnorm(0.975)*survival_observed_5_var_final^0.5)
survival_comparison5 <- data.frame(survival_predicted_5_final, survival_observed_5_final,
                                   survival_lower5_final,survival_upper5_final)

survival_comparison5$survival_upper5_final<-ifelse(survival_comparison5$survival_upper5_final>1, 1,survival_comparison5$survival_upper5_final)
survival_comparison5$underestimate<-(survival_comparison5$survival_observed_5_final-survival_comparison5$survival_predicted_5_final)/survival_comparison5$survival_observed_5_final

c2<-ggplot(data=survival_comparison5, aes(x=survival_predicted_5_final, y=survival_observed_5_final)) +
  geom_line(size=1, colour="dodgerblue3")+
  geom_errorbar(data=survival_comparison5, mapping=aes(x=survival_predicted_5_final, ymin=survival_lower5_final,
                                                       ymax=survival_upper5_final),
                colour="dodgerblue3", size=1,alpha=0.5, linetype=1)+
  geom_point(size=3, colour="dodgerblue3")+
  xlim(0.25,1)+
  ylim(0.25,1)+
  geom_abline(intercept = 0, slope = 1,lty=2)+
  labs(title="B. Calibration curve at 5 year",x="Predicted Survival Probability", y = "Observed Survival Probability")+
  theme( axis.text = element_text(size = 15), axis.title=element_text(size=15),
         panel.background = element_rect(fill = "white"),
         plot.title = element_text(size=15),
         axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
         axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm")),
         plot.margin = margin(1, 1, 1, 1, "cm"),
         axis.line = element_line(color='black'),
  )

# 10-year
survival_predicted_10_final <- exp(rowMeans(log(survival_predicted_10_combine[,-1])))
survival_observed_10_final <- exp(rowMeans(log(survival_observed_10_combine[,-1])))
survival_observed_10_var_final <- rowMeans(survival_observed_10_var_combine[,-1]) + (1+1/n_impu)*apply(survival_observed_10_var_combine[,-1], MARGIN=1, FUN=var, na.rm=TRUE)

survival_lower10_final<-exp(log(survival_observed_10_final) - qnorm(0.975)*survival_observed_10_var_final^0.5)
survival_upper10_final<-exp(log(survival_observed_10_final) + qnorm(0.975)*survival_observed_10_var_final^0.5)
survival_comparison10 <- data.frame(survival_predicted_10_final, survival_observed_10_final,
                                    survival_lower10_final,survival_upper10_final)

survival_comparison10$survival_upper10_final<-ifelse(survival_comparison10$survival_upper10_final>1, 1,survival_comparison10$survival_upper10_final)
survival_comparison10$underestimate<-(survival_comparison10$survival_observed_10_final-survival_comparison10$survival_predicted_10_final)/survival_comparison10$survival_observed_10_final

c3<-ggplot(data=survival_comparison10, aes(x=survival_predicted_10_final, y=survival_observed_10_final)) +
  geom_line(size=1, colour="dodgerblue3")+
  geom_errorbar(data=survival_comparison10, mapping=aes(x=survival_predicted_10_final, ymin=survival_lower10_final,
                                                        ymax=survival_upper10_final),
                colour="dodgerblue3", size=1,alpha=0.5, linetype=1)+
  geom_point(size=3, colour="dodgerblue3")+
  xlim(0.25,1)+
  ylim(0.25,1)+
  geom_abline(intercept = 0, slope = 1,lty=2)+
  labs(title="C. Calibration curve at 10 year",x="Predicted Survival Probability", y = "Observed Survival Probability")+
  theme(axis.text = element_text(size = 15), axis.title=element_text(size=15),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(size=15),
        axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
        axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm")),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        axis.line = element_line(color='black'),
  )

gzcal10<-ggarrange(c1, c2, c3, ncol = 3, nrow = 1)
save(gzcal10, file='gzcal10.RData')
################################################# 2. External validation ###############################################
# Same multiple imputation procedure in 1.1 was performed for validation cohort

# Data preparation
data_impu <- vector(n_impu,mode="list")
time<-subset(sybase, select = c('id', 'time'))

for (i in 1:n_impu) {
  data_impu[[i]] <- mice::complete(imputation_10, i)

  # death
  data_impu[[i]]$death<-as.numeric(data_impu[[i]]$death)-1

  # BMI
  data_impu[[i]]$BMI<-(data_impu[[i]]$Bweight/((data_impu[[i]]$Bheight)*(data_impu[[i]]$Bheight)))

  # variable transformation

  data_impu[[i]]$absT.BIL<-abs( data_impu[[i]]$logT.BIL-2.4)

  data_impu[[i]]$absGLU<-abs(data_impu[[i]]$logGLU-1.66)

  data_impu[[i]]$absTch<-abs(data_impu[[i]]$logTch-1.48)

  data_impu[[i]]$absSCR<-abs(data_impu[[i]]$logSCR-4.29)

  data_impu[[i]]$absPlatelet<-abs(data_impu[[i]]$logPlatelet-5.27)

  data_impu[[i]]$absWBC<-abs(data_impu[[i]]$logWBC-1.67)

  data_impu[[i]]$absALT<-abs(data_impu[[i]]$logALT-3.03)

  data_impu[[i]]$absTG<-abs(data_impu[[i]]$logTG-0.62)

  # select t=time, merge end time

  data_impu[[i]]<-data_impu[[i]][data_impu[[i]]$time==0, ]

  data_impu[[i]]$time<-NULL

  data_impu[[i]]<-merge(data_impu[[i]], time, by = 'id')

  # calculate LP
  data_impu[[i]]$score_infecroute <- ifelse(data_impu[[i]]$infectionroute=="IVDU",1.07363195299915,
                                            ifelse(data_impu[[i]]$infectionroute=="heterosex", 0.537387987475878,
                                                   ifelse(data_impu[[i]]$infectionroute=="Other",	0.928995957706017,0)))
  data_impu[[i]]$score_tuberculosis<-ifelse(data_impu[[i]]$Btuberculosis=="tuberculosis",0.454994655378687,0)
  data_impu[[i]]$score_HCV<-ifelse(data_impu[[i]]$HCV=="HCV",0.61112813797825,0)

  data_impu[[i]]$lp_nomogram<-data_impu[[i]]$score_infecroute + data_impu[[i]]$score_tuberculosis +
    data_impu[[i]]$score_HCV+0.0534869544053514*data_impu[[i]]$age-0.0697910578563831*data_impu[[i]]$BMI+
    0.237997613173825*data_impu[[i]]$logAST-0.662203221121992*data_impu[[i]]$loghaemoglobin-0.199524528915675*data_impu[[i]]$logCD4+
    0.777760863142926*data_impu[[i]]$absGLU+0.339954731518507*data_impu[[i]]$absPlatelet

  # create three risk groups according to the threshold in derivation cohort
  data_impu[[i]]$group3<-as.factor(ifelse(data_impu[[i]]$lp_nomogram< -3.056722, "1",
                                          ifelse(data_impu[[i]]$lp_nomogram>= -1.761427, "3", "2")))

}

################################################## 2.1 Discrimination  #################################################
# C statistics
c_harrell_apparent <- 0
c_harrell_apparent_var<-0

for (i in 1:n_impu) {
  final_model <- cph(Surv(time, death)~lp_nomogram, data=data_impu[[i]],x=TRUE,y=TRUE) ##
  c_harrell_apparent[i] <- (final_model$stats["Dxy"]+1)/2
  c_harrell_apparent_var[i] <- vardiffC(data_impu[[i]]$time, data_impu[[i]]$death, data_impu[[i]]$lp_nomogram, data_impu[[i]]$lp_nomogram)$est.varCxy
}
# C statistics and its 95%CI
C_statistics_final<-mean(c_harrell_apparent)
C_statistics_var_final<-mean(c_harrell_apparent_var) + (1+1/n_impu)*var(c_harrell_apparent_var)
C_statistics_final-qnorm(0.975)*C_statistics_var_final^0.5 # lower limit of 95%
C_statistics_final+qnorm(0.975)*C_statistics_var_final^0.5 # upper limit of 95%

################################################## 2.2 Calibration curve #################################################
# compare average predicted and observed curves
# calculate predict survival
sybasestack<-data_impu[[1]]

for (i in 1:9) {

  sybasestack<-bind_rows(sybasestack, data_impu[[i+1]])

}

calipre3<-subset(sybasestack, select = c("lp_nomogram", "group3"))

#r load base risk of derivation cohort

for (i in 1:5555) {
  calipre3[[paste0(baserisk[i,2])]] <-exp(-baserisk[i,1])^(exp(calipre3$lp_nomogram))
}
# average
predictc3<-calipre3 %>% group_by(group3) %>% summarise_all(funs(mean))
predictc3$lp_nomogram<-NULL
# wide to long
predictcl3 <- melt(setDT(predictc3), id.vars = c("group3"), variable.name = "time")
colnames(predictcl3)[3]<-'predictedsurv'
colnames(predictcl3)[2]<-'timepredict'
write_xlsx(predictcl3, path = "valipredictcl3.xlsx", col_names = T)
# read the file 'valipredictcl3.xlsx' into R, with the column 'timepredict' changed from character to numeric
predictcl3<-valipredictcl3
#calculate observed survival

# group1
data_temp1 <- subset(sybasestack, sybasestack$group3==1)

fit_observe1 <- survfit(Surv(time,death) ~ 1, data=data_temp1, id=id, cluster=id)
fit_observe1$lower
survival_observe1<- data.frame(as.data.frame(fit_observe1$surv),
                               as.data.frame(fit_observe1$lower),
                               as.data.frame(fit_observe1$upper),
                               as.data.frame(fit_observe1$time))
colnames(survival_observe1)[1]<-'observesur'
colnames(survival_observe1)[2]<-'oblower'
colnames(survival_observe1)[3]<-'obupper'
colnames(survival_observe1)[4]<-'observetime'
survival_observe1$group3<-1

# group2
data_temp2 <- subset(sybasestack, sybasestack$group3==2)

fit_observe2 <- survfit(Surv(time,death) ~ 1, data=data_temp2, id=id, cluster=id)

survival_observe2<- data.frame(as.data.frame(fit_observe2$surv),
                               as.data.frame(fit_observe2$lower),
                               as.data.frame(fit_observe2$upper),
                               as.data.frame(fit_observe2$time))
colnames(survival_observe2)[1]<-'observesur'
colnames(survival_observe2)[2]<-'oblower'
colnames(survival_observe2)[3]<-'obupper'
colnames(survival_observe2)[4]<-'observetime'
survival_observe2$group3<-2

# group3
data_temp3 <- subset(sybasestack, sybasestack$group3==3)
data_temp3$weights=0.1

fit_observe3 <- survfit(Surv(time,death) ~ 1, data=data_temp3, id=id, cluster=id)

survival_observe3<- data.frame(as.data.frame(fit_observe3$surv),
                               as.data.frame(fit_observe3$lower),
                               as.data.frame(fit_observe3$upper),
                               as.data.frame(fit_observe3$time))
colnames(survival_observe3)[1]<-'observesur'
colnames(survival_observe3)[2]<-'oblower'
colnames(survival_observe3)[3]<-'obupper'
colnames(survival_observe3)[4]<-'observetime'
survival_observe3$group3<-3

survival_observe3<-bind_rows(survival_observe1, survival_observe2, survival_observe3)

survival_observe3$group3<-as.factor(survival_observe3$group3)
predictcl3$group3<-as.factor(predictcl3$group3)

colnames(survival_observe3)[5]<-'Group'
colnames(predictcl3)[1]<-'Group'

# Plot calibration graph
valicurv<- ggplot() +
  geom_line(data=predictcl3, aes(x=timepredict, y=predictedsurv, group=Group, color=Group), linetype = "dashed", size=1)+
  geom_line(data=survival_observe3, aes(x=observetime, y=observesur, group=Group, color=Group))+
  geom_ribbon(data=survival_observe3, aes(ymin=oblower, ymax=obupper, x=observetime, y=observesur, group=Group,
                                          fill=Group),  alpha = 0.15, show.legend=F)+
  scale_x_continuous(limits = c(0,10),name='Follow-up Year', )+
  scale_y_continuous(limits = c(0.5,1), name='Survival Probability')+
  theme(axis.text = element_text(size = 15), axis.title=element_text(size=15),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color='black'),
        axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
        axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm") ))+
  scale_fill_manual(values = c("#023e7d", "#d90429", "#80b918"))+
  scale_color_manual(values = c("#023e7d", "#d90429", "#80b918"))

ggarrange(four3, sycali, nrow=1, ncol = 2, common.legend = T)


# Calibration slope and its 95%CI
n_impu=10

calibration_slope <- 0
calibration_slope_var <- 0

for (i in 1:n_impu) {
  final_model <- cph(Surv(time, death)~lp_nomogram, data=data_impu[[i]],x=TRUE,y=TRUE)
  calibration_slope[i] <- final_model$coefficients
  calibration_slope_var[i] <- final_model$var

}

calibration_slope_final <- mean(calibration_slope) # 1.012721
calibration_slope_var_final <- mean(calibration_slope_var) + (1+1/n_impu)*var(calibration_slope_var)
c(calibration_slope_final-qnorm(0.975)*calibration_slope_var_final^0.5,calibration_slope_final+qnorm(0.975)*calibration_slope_var_final^0.5)  # 95% CI

# Calibration plot
n_impu <- 10
sample_size=5751
S2_nomogram <- seq(1:sample_size)
S5_nomogram <- seq(1:sample_size)
S10_nomogram <- seq(1:sample_size)

survival_predicted_2_combine <- seq(1:10)
survival_predicted_5_combine <- seq(1:10)
survival_predicted_10_combine <- seq(1:10)

survival_observed_2_combine <- seq(1:10)
survival_observed_5_combine <- seq(1:10)
survival_observed_10_combine <- seq(1:10)

survival_observed_2_var_combine <- seq(1:10)
survival_observed_5_var_combine <- seq(1:10)
survival_observed_10_var_combine <- seq(1:10)


for (i in 1:n_impu) {
  # calculate predicted survival probability
  data_impu[[i]]$S2_nomogram <- 0.868134811^exp(data_impu[[i]]$lp_nomogram)
  data_impu[[i]]$S5_nomogram <- 0.800661362^exp(data_impu[[i]]$lp_nomogram)
  data_impu[[i]]$S10_nomogram <- 0.708730609^exp(data_impu[[i]]$lp_nomogram)

  # evenly divide patients into 10 groups

  data_impu[[i]]$group10<-cut(data_impu[[i]]$lp_nomogram, quantile(data_impu[[i]]$lp_nomogram, seq(0,1,0.1)), right=FALSE, labels=c(1:10))
  data_impu[[i]]$group10[data_impu[[i]]$lp_nomogram==max(data_impu[[i]]$lp_nomogram)] <- 10

  survival_predicted_2 <- 0
  survival_predicted_5 <- 0
  survival_predicted_10 <- 0

  # 2-year predicted survival
  survival_predicted_2 <- aggregate(data_impu[[i]]$S2_nomogram, list(data_impu[[i]]$group10), mean)
  survival_predicted_2_combine <- data.frame(survival_predicted_2_combine,survival_predicted_2$x)

  # 5-year predicted survival
  survival_predicted_5 <- aggregate(data_impu[[i]]$S5_nomogram, list(data_impu[[i]]$group10), mean)
  survival_predicted_5_combine <- data.frame(survival_predicted_5_combine,survival_predicted_5$x)

  # 10-year predicted survival
  survival_predicted_10 <- aggregate(data_impu[[i]]$S10_nomogram, list(data_impu[[i]]$group10), mean)
  survival_predicted_10_combine <- data.frame(survival_predicted_10_combine,survival_predicted_10$x)

  # observed survival
  survival_observed_2 <- 0
  survival_observed_2_var <- 0

  survival_observed_5 <- 0
  survival_observed_5_var <- 0

  survival_observed_10 <- 0
  survival_observed_10_var <- 0

  for (j in 1:10) {

    data_temp <- subset(data_impu[[i]],data_impu[[i]]$group10==j)

    fit_calibration <- survfit(Surv(time,death) ~ 1, data=data_temp)

    survival_observed_2[j] <- min(fit_calibration$surv[fit_calibration$time <= 2])
    survival_observed_2_var[j] <- (tail(fit_calibration$std.err[fit_calibration$time <= 2],1))^2

    survival_observed_5[j] <- min(fit_calibration$surv[fit_calibration$time <= 5])
    survival_observed_5_var[j] <- (tail(fit_calibration$std.err[fit_calibration$time <= 5],1))^2

    survival_observed_10[j] <- min(fit_calibration$surv[fit_calibration$time <= 10])
    survival_observed_10_var[j] <- (tail(fit_calibration$std.err[fit_calibration$time <= 10],1))^2
  }


  survival_observed_2_combine <- data.frame(survival_observed_2_combine,survival_observed_2)
  survival_observed_5_combine <- data.frame(survival_observed_5_combine,survival_observed_5)
  survival_observed_10_combine <- data.frame(survival_observed_10_combine,survival_observed_10)

  survival_observed_2_var_combine <- data.frame(survival_observed_2_var_combine,survival_observed_2_var)
  survival_observed_5_var_combine <- data.frame(survival_observed_5_var_combine,survival_observed_5_var)
  survival_observed_10_var_combine <- data.frame(survival_observed_10_var_combine,survival_observed_10_var)

}

# plot calibration plot

# 2-year

survival_predicted_2_final <- exp(rowMeans(log(survival_predicted_2_combine[,-1])))
survival_observed_2_final <- exp(rowMeans(log(survival_observed_2_combine[,-1])))
survival_observed_2_var_final <- rowMeans(survival_observed_2_var_combine[,-1]) + (1+1/n_impu)*apply(survival_observed_2_var_combine[,-1], MARGIN=1, FUN=var, na.rm=TRUE)

survival_lower2_final<-exp(log(survival_observed_2_final) - qnorm(0.975)*survival_observed_2_var_final^0.5)
survival_upper2_final<-exp(log(survival_observed_2_final) + qnorm(0.975)*survival_observed_2_var_final^0.5)
survival_comparison2 <- data.frame(survival_predicted_2_final, survival_observed_2_final,
                                   survival_lower2_final,survival_upper2_final)

survival_comparison2$survival_upper2_final<-ifelse(survival_comparison2$survival_upper2_final>1, 1,survival_comparison2$survival_upper1_final)
survival_comparison2$underestimate<-(survival_comparison2$survival_observed_2_final-survival_comparison2$survival_predicted_2_final)/survival_comparison2$survival_observed_2_final

c1<-ggplot(data=survival_comparison2, aes(x=survival_predicted_2_final, y=survival_observed_2_final)) +
  geom_line(size=1, colour="dodgerblue3")+
  geom_errorbar(data=survival_comparison2, mapping=aes(x=survival_predicted_2_final, ymin=survival_lower2_final,
                                                       ymax=survival_upper2_final),
                colour="dodgerblue3", size=1,alpha=0.5, linetype=1)+
  geom_point(size=3, colour="dodgerblue3")+
  xlim(0.25,1)+
  ylim(0.25,1)+
  geom_abline(intercept = 0, slope = 1,lty=2)+
  labs(title="A. Calibration curve at 2 year",x="Predicted Survival Probability", y = "Observed Survival Probability")+
  theme( axis.text = element_text(size = 15), axis.title=element_text(size=15),
         panel.background = element_rect(fill = "white"),
         plot.title = element_text(size=15),
         axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
         axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm")),
         plot.margin = margin(1, 1, 1, 1, "cm"),
         axis.line = element_line(color='black'),
  )

# 5-year

survival_predicted_5_final <- exp(rowMeans(log(survival_predicted_5_combine[,-1])))
survival_observed_5_final <- exp(rowMeans(log(survival_observed_5_combine[,-1])))
survival_observed_5_var_final <- rowMeans(survival_observed_5_var_combine[,-1]) + (1+1/n_impu)*apply(survival_observed_5_var_combine[,-1], MARGIN=1, FUN=var, na.rm=TRUE)

survival_lower5_final<-exp(log(survival_observed_5_final) - qnorm(0.975)*survival_observed_5_var_final^0.5)
survival_upper5_final<-exp(log(survival_observed_5_final) + qnorm(0.975)*survival_observed_5_var_final^0.5)
survival_comparison5 <- data.frame(survival_predicted_5_final, survival_observed_5_final,
                                   survival_lower5_final,survival_upper5_final)

survival_comparison5$survival_upper5_final<-ifelse(survival_comparison5$survival_upper5_final>1, 1,survival_comparison5$survival_upper5_final)
survival_comparison5$underestimate<-(survival_comparison5$survival_observed_5_final-survival_comparison5$survival_predicted_5_final)/survival_comparison5$survival_observed_5_final

c2<-ggplot(data=survival_comparison5, aes(x=survival_predicted_5_final, y=survival_observed_5_final)) +
  geom_line(size=1, colour="dodgerblue3")+
  geom_errorbar(data=survival_comparison5, mapping=aes(x=survival_predicted_5_final, ymin=survival_lower5_final,
                                                       ymax=survival_upper5_final),
                colour="dodgerblue3", size=1,alpha=0.5, linetype=1)+
  geom_point(size=3, colour="dodgerblue3")+
  xlim(0.25,1)+
  ylim(0.25,1)+
  geom_abline(intercept = 0, slope = 1,lty=2)+
  labs(title="B. Calibration curve at 5 year",x="Predicted Survival Probability", y = "Observed Survival Probability")+
  theme( axis.text = element_text(size = 15), axis.title=element_text(size=15),
         panel.background = element_rect(fill = "white"),
         plot.title = element_text(size=15),
         axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
         axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm")),
         plot.margin = margin(1, 1, 1, 1, "cm"),
         axis.line = element_line(color='black'),
  )

# 10-year
survival_predicted_10_final <- exp(rowMeans(log(survival_predicted_10_combine[,-1])))
survival_observed_10_final <- exp(rowMeans(log(survival_observed_10_combine[,-1])))
survival_observed_10_var_final <- rowMeans(survival_observed_10_var_combine[,-1]) + (1+1/n_impu)*apply(survival_observed_10_var_combine[,-1], MARGIN=1, FUN=var, na.rm=TRUE)

survival_lower10_final<-exp(log(survival_observed_10_final) - qnorm(0.975)*survival_observed_10_var_final^0.5)
survival_upper10_final<-exp(log(survival_observed_10_final) + qnorm(0.975)*survival_observed_10_var_final^0.5)
survival_comparison10 <- data.frame(survival_predicted_10_final, survival_observed_10_final,
                                    survival_lower10_final,survival_upper10_final)

survival_comparison10$survival_upper10_final<-ifelse(survival_comparison10$survival_upper10_final>1, 1,survival_comparison10$survival_upper10_final)
survival_comparison10$underestimate<-(survival_comparison10$survival_observed_10_final-survival_comparison10$survival_predicted_10_final)/survival_comparison10$survival_observed_10_final

c3<-ggplot(data=survival_comparison10, aes(x=survival_predicted_10_final, y=survival_observed_10_final)) +
  geom_line(size=1, colour="dodgerblue3")+
  geom_errorbar(data=survival_comparison10, mapping=aes(x=survival_predicted_10_final, ymin=survival_lower10_final,
                                                        ymax=survival_upper10_final),
                colour="dodgerblue3", size=1,alpha=0.5, linetype=1)+
  geom_point(size=3, colour="dodgerblue3")+
  xlim(0.25,1)+
  ylim(0.25,1)+
  geom_abline(intercept = 0, slope = 1,lty=2)+
  labs(title="C. Calibration curve at 10 year",x="Predicted Survival Probability", y = "Observed Survival Probability")+
  theme(axis.text = element_text(size = 15), axis.title=element_text(size=15),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(size=15),
        axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
        axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm")),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        axis.line = element_line(color='black'),
  )

sycal10<-ggarrange(c1, c2, c3, ncol = 3, nrow = 1)

save(sycal10, file='sycal10.RData')
ggarrange(gzcal10, sycal10, ncol = 1, nrow = 2)


################################################# 3. Landmarking analysis ###############################################
## Derivation cohort
# Data preparation
# Extract baseline variables, full follow-up time, and endpoint
data_impu_base <- data_impu
for (i in 1:n_impu){
  data_impu_base[[i]]$BMI<-(data_impu_base[[i]]$Bweight/((data_impu_base[[i]]$Bheight)*(data_impu_base[[i]]$Bheight)))
  data_impu_base[[i]]<-subset(data_impu_base[[i]], select=c('id','age', 'infectionroute','Btuberculosis','HCV','BMI','time','death' ))
}

# For original imputation data,  retain follow-up variables, variables transformation

data_impu <- vector(n_impu,mode="list")
for (i in 1:n_impu) {
  data_impu[[i]] <- mice::complete(gzimputation_10, i)

  # death
  data_impu[[i]]$death<-as.numeric(data_impu[[i]]$death)-1

  # variable tranformation

  data_impu[[i]]$absT.BIL<-abs( data_impu[[i]]$logT.BIL-2.4)

  data_impu[[i]]$absGLU<-abs(data_impu[[i]]$logGLU-1.66)

  data_impu[[i]]$absTch<-abs(data_impu[[i]]$logTch-1.48)

  data_impu[[i]]$absSCR<-abs(data_impu[[i]]$logSCR-4.29)

  data_impu[[i]]$absPlatelet<-abs(data_impu[[i]]$logPlatelet-5.27)

  data_impu[[i]]$absWBC<-abs(data_impu[[i]]$logWBC-1.67)

  data_impu[[i]]$absALT<-abs(data_impu[[i]]$logALT-3.03)

  data_impu[[i]]$absTG<-abs(data_impu[[i]]$logTG-0.62)

  data_impu[[i]]$lgCD4CD8<-data_impu[[i]]$logCD4 - data_impu[[i]]$logCD8

}

# 0.5 year after ART initiation
data_impu_l6 <- vector(n_impu,mode="list")
data_impu_6m <- vector(n_impu,mode="list")

for (i in 1:n_impu) {

  # determine population, original population with follow-up >0.5 year
  data_impu_base[[i]]<-subset(data_impu_base[[i]], time>0.5)
  # update age
  data_impu_base[[i]]$age<-data_impu_base[[i]]$age+0.5
  # update start time
  data_impu_base[[i]]$time<-data_impu_base[[i]]$time-0.5

  # find the latest biochemical parameters among those patients
  data_impu_l6[[i]]<-subset(data_impu[[i]], time<=0.5)
  data_impu_l6[[i]]<-data_impu_l6[[i]][(data_impu_l6[[i]]$id %in% data_impu_base[[i]]$id),]

  data_impu_l6[[i]]$timedf<-0.5-data_impu_l6[[i]]$time
  data_impu_l6[[i]]<-data_impu_l6[[i]] %>%
    group_by(id) %>%
    slice(which.min(timedf))

  data_impu_l6[[i]]<-subset(data_impu_l6[[i]], select=c('id', 'logAST', 'loghaemoglobin','lgCD4CD8',
                                                        'logCD4', 'absT.BIL',
                                                        'absGLU', 'absPlatelet'))

  data_impu_6m[[i]]<-merge(data_impu_l6[[i]], data_impu_base[[i]], by='id')

  # calculate LP
  data_impu_6m[[i]]$score_infecroute <- ifelse(data_impu_6m[[i]]$infectionroute=="IVDU",1.07363195299915,
                                               ifelse(data_impu_6m[[i]]$infectionroute=="heterosex", 0.537387987475878,
                                                      ifelse(data_impu_6m[[i]]$infectionroute=="Other",	0.928995957706017,0)))
  data_impu_6m[[i]]$score_tuberculosis<-ifelse(data_impu_6m[[i]]$Btuberculosis=="tuberculosis", 0.454994655378687,0)
  data_impu_6m[[i]]$score_HCV<-ifelse(data_impu_6m[[i]]$HCV=="HCV",0.61112813797825,0)
  data_impu_6m[[i]]$lp_nomogram<-data_impu_6m[[i]]$score_infecroute + data_impu_6m[[i]]$score_tuberculosis +
    data_impu_6m[[i]]$score_HCV+ 0.0534869544053514*data_impu_6m[[i]]$age-0.0697910578563831*data_impu_6m[[i]]$BMI+
    0.237997613173825*data_impu_6m[[i]]$logAST-0.662203221121992*data_impu_6m[[i]]$loghaemoglobin-0.199524528915675*data_impu_6m[[i]]$logCD4
    +0.777760863142926*data_impu_6m[[i]]$absGLU+0.339954731518507*data_impu_6m[[i]]$absPlatelet
}

# calculate  C statistics and 95% CI
# CD4
c_harrell_lgCD4 <- 0
Cse_logCD4 <- 0

for (i in 1:n_impu) {
  C_lgCD4<-rcorr.cens(data_impu_6m[[i]]$logCD4, Surv(data_impu_6m[[i]]$time, data_impu_6m[[i]]$death))
  c_harrell_lgCD4[i]<- C_lgCD4['C Index']
  Cse_logCD4[i]<-C_lgCD4['S.D.']/2
}

# CD4/CD8
c_harrell_lgCD4CD8 <- 0
Cse_logCD4CD8 <- 0
for (i in 1:n_impu) {
  C_lgCD4CD8<-rcorr.cens(data_impu_6m[[i]]$lgCD4CD8, Surv(data_impu_6m[[i]]$time, data_impu_6m[[i]]$death))
  c_harrell_lgCD4CD8[i]<- C_lgCD4CD8['C Index']
  Cse_logCD4CD8[i]<-C_lgCD4CD8['S.D.']/2

}

# Prediction model
c_harrell_lp <- 0
Cse_lp <- 0
for (i in 1:n_impu) {
  C_lp<-rcorr.cens(data_impu_6m[[i]]$lp_nomogram, Surv(data_impu_6m[[i]]$time, data_impu_6m[[i]]$death))
  c_harrell_lp[i]<- C_lp['C Index']
  Cse_lp[i]<-C_lp['S.D.']/2
}

landmark<-rep('6 months')
Type<-c('CD4', 'CD4/CD8', 'Prediction model')
Cstatistics<-c(mean(c_harrell_lgCD4), mean(c_harrell_lgCD4CD8), mean(c_harrell_lp))
Cse<-c(mean(Cse_logCD4), mean(Cse_logCD4CD8), mean(Cse_lp))

Comp3C<-data.frame(landmark, Type, Cstatistics, Cse)

Comp3C<-within.data.frame(Comp3C,{
  cll<-Cstatistics-qnorm(0.975)*Cse
  cul<-Cstatistics+qnorm(0.975)*Cse

})

Comp3C[3,3]<-1- Comp3C[3,3]
Comp3C[3,5]<-1-  Comp3C[3,5]
Comp3C[3,6]<-1-  Comp3C[3,6]

# Same procedure  performed for t = 1, 2, 3, 4, and 5 years after ART initiation
# Merge all time points
Comp3C_final<-bind_rows(Comp3C, Comp3C1y, Comp3C2y, Comp3C3y, Comp3C4y, Comp3C5y)
Comp3C_final$landmark<-as.factor(Comp3C_final$landmark)
# Plot figure
compare3Cgzplot<-ggplot(Comp3C_final, aes(landmark, Cstatistics)) +
  geom_errorbar(
    aes(ymin = cll, ymax = cul, color = Type),
    position = position_dodge(0.3), width = 0.2, size=0.6
  )+
  ylim(0,1)+
  geom_point(aes(color = Type), position = position_dodge(0.3), size=2.5) +
  scale_color_manual(values = c("#023e8a", "#43aa8b", "#dc2f02"))+
  labs(x="Landmark time (years)", y = "C statistics")+
  geom_hline(yintercept=0.5, linetype="dashed", color = "gray")+
  theme(axis.text = element_text(size = 15), axis.title=element_text(size=15),
        panel.background = element_rect(fill = "white"),
        legend.position="top",
        legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill = "white"),
        legend.key.size = unit(1.5, "cm"),
        legend.key.width = unit(0.5,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size=13),
        legend.spacing.x = unit(0.5, 'cm'),
        axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
        axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm") ))+
  theme(axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1))+
  scale_x_discrete(limits=c("6 months", "1 year", "2 years", "3 years", "4 years", "5 years"))

# add a table under figure
# make table
data_table <- ggplot(Comp3C_final, aes(x = landmark, y = Type,
                                       label = format(result, nsmall = 1), colour = Type)) +
  geom_text(size = 6) + theme_bw() +
  scale_y_discrete(limits = c("CD4/CD8", "CD4", "Prediction model"))+
  scale_color_manual(values = c("#023e8a", "#43aa8b", "#dc2f02"))+
  scale_x_discrete(limits=c("6 months", "1 year", "2 years", "3 years", "4 years", "5 years"))+
  theme(axis.text = element_text(size = 15), axis.title=element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position="none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(colour=c("#43aa8b", "#023e8a","#dc2f02")))

# combine table and figure
Layout <- grid.layout(nrow = 2, ncol = 1, heights = unit(c(2,0.35), c("null", "null")))
grid.show.layout(Layout)
vplayout <- function(...) {
  grid.newpage()
  pushViewport(viewport(layout = Layout))
}
subplot <- function(x, y) viewport(layout.pos.row = x,
                                   layout.pos.col = y)

mmplot <- function(a, b) {
  vplayout()
  print(a, vp = subplot(1, 1))
  print(b, vp = subplot(2, 1))
}

mmplot(compare3Cgzplot, data_table)

# Same procedure performed for Validation cohort


