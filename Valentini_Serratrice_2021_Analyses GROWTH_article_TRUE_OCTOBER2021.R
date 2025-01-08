library(tidyverse)

###DATA HANDLING ##################################################################################################

###data handling###########################################################################################
DATA_growth_NOV_2019_articleNEW <- read_csv("Valentini_Serratrice_2021_DATA_growth_NOV_2019_articleNEW.csv")
DATA_growth_NOV_2019_articleNEW->DataverySHORT_T3 # short data for growth models


###data very short###

DataverySHORT_T3$Subject <- factor(DataverySHORT_T3$Subject)
DataverySHORT_T3$School <- factor(DataverySHORT_T3$School)
DataverySHORT_T3$Gender <- factor(DataverySHORT_T3$Gender)
DataverySHORT_T3$READ_passages_T3 <- factor(DataverySHORT_T3$READ_passages_T3)
DataverySHORT_T3$SEN <- factor(DataverySHORT_T3$SEN)

summary(DataverySHORT_T3)


###no SEN###

DataverySHORT_T3_NEW<-subset(DataverySHORT_T3,DataverySHORT_T3$SEN=="NO")
summary(DataverySHORT_T3_NEW)

sd(DataverySHORT_T3_NEW$ENGexposureTOT, na.rm=TRUE)
sd(DataverySHORT_T3_NEW$ENGexposurePERC, na.rm=TRUE)


###scaling covariates###

DataverySHORT_T3_NEW$Age_T1<-scale(DataverySHORT_T3_NEW$Age_months_1,scale=TRUE,center=TRUE)
DataverySHORT_T3_NEW$Mother_EDU<-scale(DataverySHORT_T3_NEW$Mum_EDU,scale=TRUE,center=TRUE)
DataverySHORT_T3_NEW$SES<-scale(DataverySHORT_T3_NEW$higherSES,scale=TRUE,center=TRUE)
DataverySHORT_T3_NEW$INPUT<-scale(DataverySHORT_T3_NEW$EngCUMInput,scale=TRUE,center=TRUE)

DataverySHORT_T3_NEW$lenghEXPOSURE_ENG_months<-scale(DataverySHORT_T3_NEW$ENGexposureTOT,scale=TRUE,center=TRUE)
DataverySHORT_T3_NEW$CURRENT_ENGinputPERC<-scale(DataverySHORT_T3_NEW$ENGexposurePERC,scale=TRUE,center=TRUE)

DataverySHORT_T3_NEW$SESrev<-(DataverySHORT_T3_NEW$SES*(-1))



summary(DataverySHORT_T3_NEW)

########correlations##################################################################################################

corstarsl <- function(x){ 
  require(Hmisc) 
  x <- as.matrix(x) 
  R <- rcorr(x, type = "spearman")$r 
  p <- rcorr(x, type = "spearman")$P 
  
  ## define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ", " ")))
  
  ## trunctuate the matrix that holds the correlations to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1] 
  
  ## build a new matrix that includes the correlations with their apropriate stars 
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x)) 
  diag(Rnew) <- paste(diag(R), " ", sep="") 
  rownames(Rnew) <- colnames(x) 
  colnames(Rnew) <- paste(colnames(x), "", sep="") 
  
  ## remove lower triangle of correlation matrix
  Rnew <- as.matrix(Rnew)
  Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
  Rnew <- as.data.frame(Rnew)
  
  
  ## remove last column and return the matrix (which is now a data frame)
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  return(Rnew) 
}



corstarslPEARSON <- function(x){ 
  require(Hmisc) 
  x <- as.matrix(x) 
  R <- rcorr(x, type = "pearson")$r 
  p <- rcorr(x, type = "pearson")$P 
  
  ## define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ", " ")))
  
  ## trunctuate the matrix that holds the correlations to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1] 
  
  ## build a new matrix that includes the correlations with their apropriate stars 
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x)) 
  diag(Rnew) <- paste(diag(R), " ", sep="") 
  rownames(Rnew) <- colnames(x) 
  colnames(Rnew) <- paste(colnames(x), "", sep="") 
  
  ## remove lower triangle of correlation matrix
  Rnew <- as.matrix(Rnew)
  Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
  Rnew <- as.data.frame(Rnew)
  
  
  ## remove last column and return the matrix (which is now a data frame)
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  return(Rnew) 
}


###CORRELATION MATRIX

names(DataverySHORT_T3_NEW)

library(ggplot2)
library(moments)

my_dataSHORT <- DataverySHORT_T3_NEW[, c(30,35,29,33,34,15,17,16,18,20,22,21,23,25,27,26,28,13)]
# print the first 6 rows
head(my_dataSHORT, 6)
summary(my_dataSHORT)

###normality tests
shapiro.test(my_dataSHORT$Mother_EDU)
shapiro.test(my_dataSHORT$SESrev)
ks.test(my_dataSHORT$SESrev, pnorm)
shapiro.test(my_dataSHORT$Age_T1)
shapiro.test(my_dataSHORT$lenghEXPOSURE_ENG_months)
shapiro.test(my_dataSHORT$CURRENT_ENGinputPERC)
shapiro.test(my_dataSHORT$VOC_BPVS1)
shapiro.test(my_dataSHORT$VOC_mean_SYNOPP1)
shapiro.test(my_dataSHORT$GRAMM_WordStruct1)
shapiro.test(my_dataSHORT$GRAMM_TROGshort1)
shapiro.test(my_dataSHORT$VOC_BPVS_2)
shapiro.test(my_dataSHORT$VOC_meanSYNOPP_2)
shapiro.test(my_dataSHORT$GRAMM_WordStruct_2)
shapiro.test(my_dataSHORT$GRAMM_TROGshort_2)
shapiro.test(my_dataSHORT$VOC_BPVS_3)
shapiro.test(my_dataSHORT$VOC_meanSYNOPP_3)
shapiro.test(my_dataSHORT$GRAMM_WordStruct_3)
shapiro.test(my_dataSHORT$GRAMM_TROGshort_3)



skewness(my_dataSHORT, na.rm = TRUE)
kurtosis(my_dataSHORT, na.rm = TRUE)

###(the p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution. In other words, we can assume the normality.)


corstarsl(my_dataSHORT) 
corstarslPEARSON(my_dataSHORT)


########BPVS GROWTH###################################################################################################

# For each variable, a linear growth model comprising of both an intercept and a linear slope was estimated using the package lavaan (Rosseel, 2012) in R (R core Team, 2019).
# A full information maximum likelihood approach was used to deal with missing values for both endogenous and exogenous variables, by specifying exogenous variables as random;
# the software estimated a likelihood function for each participant based on the data collected.
# This model included all possible significant covariates: age at Time 1, mother's education, highest household occupation (henceforth occupation), length of exposure to English (English length exposure), and current amount of English input (English current input) were entered in the model if they significantly correlated with a given outcome measure.
# The significance of the linear growth was assessed by the comparison of this model with a model containing only an intercept but not a slope (no growth model) as suggested by Grimm, Ram, and Estabrook (2017).


summary(DataverySHORT_T3_NEW)

###no growth model###

ng.BPVS.lavaan_model <- '
# latent variable definitions
#intercept
eta_1 =~ 1*VOC_BPVS1
eta_1 =~ 1*VOC_BPVS_2
eta_1 =~ 1*VOC_BPVS_3

# factor variances
eta_1 ~~ eta_1

# covariances among factors 
#none (only 1 factor)

# factor means 
eta_1 ~ start(30)*1

# manifest variances (made equivalent by naming theta)
VOC_BPVS1 ~~ theta*VOC_BPVS1
VOC_BPVS_2 ~~ theta*VOC_BPVS_2
VOC_BPVS_3 ~~ theta*VOC_BPVS_3

# manifest means (fixed at zero)
VOC_BPVS1 ~ 0*1
VOC_BPVS_2 ~ 0*1
VOC_BPVS_3 ~ 0*1

# modification_1#
VOC_BPVS_2 ~ VOC_BPVS1
VOC_BPVS_3 ~ VOC_BPVS_2

#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1~Mother_EDU+SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Mother_EDU ~~ Mother_EDU
SES ~~ SES

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ Mother_EDU
lenghEXPOSURE_ENG_months ~~ SES
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
CURRENT_ENGinputPERC ~~ Mother_EDU
CURRENT_ENGinputPERC ~~ SES
Mother_EDU ~~ SES

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
Mother_EDU ~ 1
CURRENT_ENGinputPERC ~ 1
SES ~ 1
'
ng.BPVS.lavaan_fit <- sem(ng.BPVS.lavaan_model, 
                          data = DataverySHORT_T3_NEW,
                          meanstructure = TRUE,
                          estimator = "ML",
                          missing = "fiml")


summary(ng.BPVS.lavaan_fit, fit.measures=TRUE, standardized = T)


summary(ng.BPVS.lavaan_fit)
fitMeasures(ng.BPVS.lavaan_fit)

semPlot::semPaths(ng.BPVS.lavaan_fit, title = FALSE,layout = "tree", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(ng.BPVS.lavaan_fit, title = FALSE,layout = "circle", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(ng.BPVS.lavaan_fit, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)


parameterEstimates(ng.BPVS.lavaan_fit)
inspect(ng.BPVS.lavaan_fit, what="est")

###linear growth model###

lg.BPVS.lavaan_model <- '
# latent variable definitions
#intercept (note intercept is a reserved term)
eta_1 =~ 1*VOC_BPVS1
eta_1 =~ 1*VOC_BPVS_2
eta_1 =~ 1*VOC_BPVS_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*VOC_BPVS1
eta_2 =~ 1*VOC_BPVS_2
eta_2 =~ 2*VOC_BPVS_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2

# covariances among factors 
eta_1 ~~ eta_2

# factor means 
eta_1 ~ start(35)*1
eta_2 ~ start(4)*1

# manifest variances (made equivalent by naming theta)
VOC_BPVS1 ~~ theta*VOC_BPVS1
VOC_BPVS_2 ~~ theta*VOC_BPVS_2
VOC_BPVS_3 ~~ theta*VOC_BPVS_3

# manifest means (fixed at zero)
VOC_BPVS1 ~ 0*1
VOC_BPVS_2 ~ 0*1
VOC_BPVS_3 ~ 0*1


# modification_1#
VOC_BPVS_2 ~ VOC_BPVS1
VOC_BPVS_3 ~ VOC_BPVS_2

#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1~Mother_EDU+SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~Mother_EDU+SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Mother_EDU ~~ Mother_EDU
SES ~~ SES

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ Mother_EDU
lenghEXPOSURE_ENG_months ~~ SES
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
CURRENT_ENGinputPERC ~~ Mother_EDU
CURRENT_ENGinputPERC ~~ SES
Mother_EDU ~~ SES

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
Mother_EDU ~ 1
CURRENT_ENGinputPERC ~ 1
SES ~ 1
'

lg.BPVS.lavaan_model_fit <- sem(lg.BPVS.lavaan_model, 
                                data = DataverySHORT_T3_NEW,
                                meanstructure = TRUE,
                                estimator = "ML",
                                missing = "fiml")

summary(lg.BPVS.lavaan_model_fit)
summary(lg.BPVS.lavaan_model_fit, fit.measures=TRUE, standardized = T)
fitMeasures(lg.BPVS.lavaan_model_fit)

semPlot::semPaths(lg.BPVS.lavaan_model_fit, title = FALSE,layout = "tree", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(lg.BPVS.lavaan_model_fit, title = FALSE,layout = "circle", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(lg.BPVS.lavaan_model_fit, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)


anova(lg.BPVS.lavaan_model_fit,ng.BPVS.lavaan_fit)

##comparison with no growth
anova(lg.BPVS.lavaan_model_fit,ng.BPVS.lavaan_fit)


###quadratic growth checks ####

###quadratic growth model### use Quad_model3!

Quad_model <- '
# latent variable definitions
#intercept (note intercept is a reserved term)
eta_1 =~ 1*VOC_BPVS1
eta_1 =~ 1*VOC_BPVS_2
eta_1 =~ 1*VOC_BPVS_3

#linear slope (note intercept is a reserved term)
eta_2 =~ -0.5*VOC_BPVS1
eta_2 =~ 0*VOC_BPVS_2
eta_2 =~ 0.5*VOC_BPVS_3

#quadratic change factor loadings
eta_3 =~ 0.25*VOC_BPVS1
eta_3 =~0*VOC_BPVS_2
eta_3 =~0.25*VOC_BPVS_3

#manifest variances to be same over time
VOC_BPVS1 ~~ start(1)*theta*VOC_BPVS1
VOC_BPVS_2 ~~ theta*VOC_BPVS_2
VOC_BPVS_3 ~~ theta*VOC_BPVS_3

#latent variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2
eta_3 ~~ eta_3

#latent covariances
eta_1 ~~ eta_2
eta_1 ~~ eta_3
eta_2 ~~ eta_3

# factor means 
eta_1 ~ start(35)*1
eta_2 ~ start(4)*1

#manifest means fixed to 0
VOC_BPVS1 ~ 0*1
VOC_BPVS_2 ~ 0*1
VOC_BPVS_3 ~ 0*1

# modification_1#
VOC_BPVS_2 ~ VOC_BPVS1
VOC_BPVS_3 ~ VOC_BPVS_2

#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1~Mother_EDU+SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~Mother_EDU+SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Mother_EDU ~~ Mother_EDU
SES ~~ SES

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ Mother_EDU
lenghEXPOSURE_ENG_months ~~ SES
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
CURRENT_ENGinputPERC ~~ Mother_EDU
CURRENT_ENGinputPERC ~~ SES
Mother_EDU ~~ SES

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
Mother_EDU ~ 1
CURRENT_ENGinputPERC ~ 1
SES ~ 1
'
quad_fit <- sem(Quad_model, data = DataverySHORT_T3_NEW,
                meanstructure = TRUE,
                estimator = "ML",
                missing = "fiml")
summary(quad_fit, fit.measures = T, standardized = T)

##estimated variances are negative - changing model##

Quad_model2 <- '
# latent variable definitions
#intercept (note intercept is a reserved term)
eta_1 =~ 1*VOC_BPVS1
eta_1 =~ 1*VOC_BPVS_2
eta_1 =~ 1*VOC_BPVS_3

#linear slope (note intercept is a reserved term)
eta_2 =~ -0.5*VOC_BPVS1
eta_2 =~ 0*VOC_BPVS_2
eta_2 =~ 0.5*VOC_BPVS_3

#quadratic change factor loadings
eta_3 =~ 0.25*VOC_BPVS1
eta_3 =~0*VOC_BPVS_2
eta_3 =~0.25*VOC_BPVS_3

#manifest variances to be same over time
VOC_BPVS1 ~~ start(1)*theta*VOC_BPVS1
VOC_BPVS_2 ~~ theta*VOC_BPVS_2
VOC_BPVS_3 ~~ theta*VOC_BPVS_3

#latent variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2
eta_3 ~~ eta_3

#latent covariances
eta_1 ~~ eta_2
eta_1 ~~ 0*eta_3 #fixing the covariance to 0
eta_2 ~~ 0*eta_3 #fixing the covariance to 0

# factor means 
eta_1 ~ start(35)*1
eta_2 ~ start(4)*1

#manifest means fixed to 0
VOC_BPVS1 ~ 0*1
VOC_BPVS_2 ~ 0*1
VOC_BPVS_3 ~ 0*1

# modification_1#
VOC_BPVS_2 ~ VOC_BPVS1
VOC_BPVS_3 ~ VOC_BPVS_2

#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1~Mother_EDU+SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~Mother_EDU+SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Mother_EDU ~~ Mother_EDU
SES ~~ SES

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ Mother_EDU
lenghEXPOSURE_ENG_months ~~ SES
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
CURRENT_ENGinputPERC ~~ Mother_EDU
CURRENT_ENGinputPERC ~~ SES
Mother_EDU ~~ SES

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
Mother_EDU ~ 1
CURRENT_ENGinputPERC ~ 1
SES ~ 1
'
quad_fit2 <- sem(Quad_model2, data = DataverySHORT_T3_NEW,
                meanstructure = TRUE,
                estimator = "ML",
                missing = "fiml")
summary(quad_fit2, fit.measures = T, standardized = T)

##comparison with linear growth
anova(lg.BPVS.lavaan_model_fit,quad_fit2) ### no significant quadratic model compared to linear growth model


QuadModel3 <- '
# latent variable definitions
#intercept (note intercept is a reserved term)
eta_1 =~ 1*VOC_BPVS1
eta_1 =~ 1*VOC_BPVS_2
eta_1 =~ 1*VOC_BPVS_3

#linear slope (note intercept is a reserved term)
eta_2 =~ -0.5*VOC_BPVS1
eta_2 =~ 0*VOC_BPVS_2
eta_2 =~ 0.5*VOC_BPVS_3

#quadratic change factor loadings
eta_3 =~ 0.25*VOC_BPVS1
eta_3 =~0*VOC_BPVS_2
eta_3 =~0.25*VOC_BPVS_3

#manifest variances to be same over time
VOC_BPVS1 ~~ start(1)*theta*VOC_BPVS1
VOC_BPVS_2 ~~ theta*VOC_BPVS_2
VOC_BPVS_3 ~~ theta*VOC_BPVS_3

#latent variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2
eta_3 ~~  1*eta_3

#latent covariances
eta_1 ~~ eta_2
eta_1 ~~ 0*eta_3 #fixing the covariance to 0
eta_2 ~~ 0*eta_3 #fixing the covariance to 0

# factor means 
eta_1 ~ start(35)*1
eta_2 ~ start(4)*1

#manifest means fixed to 0
VOC_BPVS1 ~ 0*1
VOC_BPVS_2 ~ 0*1
VOC_BPVS_3 ~ 0*1

# modification_1#
VOC_BPVS_2 ~ VOC_BPVS1
VOC_BPVS_3 ~ VOC_BPVS_2

#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1~lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Mother_EDU ~~ Mother_EDU
SES ~~ SES

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ Mother_EDU
lenghEXPOSURE_ENG_months ~~ SES
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
CURRENT_ENGinputPERC ~~ Mother_EDU
CURRENT_ENGinputPERC ~~ SES
Mother_EDU ~~ SES

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
Mother_EDU ~ 1
CURRENT_ENGinputPERC ~ 1
SES ~ 1
'
quad_fit3 <- sem(QuadModel3, data = DataverySHORT_T3_NEW,
                 meanstructure = TRUE,
                 estimator = "ML",
                 missing = "fiml")
summary(quad_fit3, fit.measures = T, standardized = T)


anova(lg.BPVS.lavaan_model_fit,quad_fit3) ### no significant quadratic model compared to linear growth model



##covariates##  use lgCOV.BPVS.FINAL_model

lgCOV.BPVS.lavaan_model <- '
# latent variable definitions
#intercept (note intercept is a reserved term)
eta_1 =~ 1*VOC_BPVS1
eta_1 =~ 1*VOC_BPVS_2
eta_1 =~ 1*VOC_BPVS_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*VOC_BPVS1
eta_2 =~ 1*VOC_BPVS_2
eta_2 =~ 2*VOC_BPVS_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2

# covariances among factors 
eta_1 ~~ eta_2

# factor means 
eta_1 ~ start(35)*1
eta_2 ~ start(4)*1

# manifest variances (made equivalent by naming theta)
VOC_BPVS1 ~~ theta*VOC_BPVS1
VOC_BPVS_2 ~~ theta*VOC_BPVS_2
VOC_BPVS_3 ~~ theta*VOC_BPVS_3

# manifest means (fixed at zero)
VOC_BPVS1 ~ 0*1
VOC_BPVS_2 ~ 0*1
VOC_BPVS_3 ~ 0*1


# modification_1#
VOC_BPVS_2 ~ VOC_BPVS1
VOC_BPVS_3 ~ VOC_BPVS_2

#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1~Mother_EDU+SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~Mother_EDU+SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Mother_EDU ~~ Mother_EDU
SES ~~ SES

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ Mother_EDU
lenghEXPOSURE_ENG_months ~~ SES
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
CURRENT_ENGinputPERC ~~ Mother_EDU
CURRENT_ENGinputPERC ~~ SES
Mother_EDU ~~ SES

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
Mother_EDU ~ 1
CURRENT_ENGinputPERC ~ 1
SES ~ 1
'

lg.BPVS.COV_fit <- sem(lgCOV.BPVS.lavaan_model, 
                       data = DataverySHORT_T3_NEW,
                       meanstructure = TRUE,
                       estimator = "ML",
                       missing = "fiml")

summary(lg.BPVS.COV_fit, fit.measures=TRUE, standardized = T)
fitMeasures(lg.BPVS.COV_fit)

semPlot::semPaths(lg.BPVS.COV_fit, title = FALSE,layout = "tree", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(lg.BPVS.COV_fit, title = FALSE,layout = "circle", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(lg.BPVS.COV_fit, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)


anova(lg.BPVS.lavaan_model_fit,lg.BPVS.COV_fit)

lgCOV2.BPVS.lavaan_model <- '
# latent variable definitions
#intercept (note intercept is a reserved term)
eta_1 =~ 1*VOC_BPVS1
eta_1 =~ 1*VOC_BPVS_2
eta_1 =~ 1*VOC_BPVS_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*VOC_BPVS1
eta_2 =~ 1*VOC_BPVS_2
eta_2 =~ 2*VOC_BPVS_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2

# covariances among factors 
eta_1 ~~ eta_2

# factor means 
eta_1 ~ start(35)*1
eta_2 ~ start(4)*1

# manifest variances (made equivalent by naming theta)
VOC_BPVS1 ~~ theta*VOC_BPVS1
VOC_BPVS_2 ~~ theta*VOC_BPVS_2
VOC_BPVS_3 ~~ theta*VOC_BPVS_3

# manifest means (fixed at zero)
VOC_BPVS1 ~ 0*1
VOC_BPVS_2 ~ 0*1
VOC_BPVS_3 ~ 0*1


# modification_1#
VOC_BPVS_2 ~ VOC_BPVS1
VOC_BPVS_3 ~ VOC_BPVS_2

#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1~0*Mother_EDU+0*SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~0*Mother_EDU+0*SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Mother_EDU ~~ Mother_EDU
SES ~~ SES

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ Mother_EDU
lenghEXPOSURE_ENG_months ~~ SES
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
CURRENT_ENGinputPERC ~~ Mother_EDU
CURRENT_ENGinputPERC ~~ SES
Mother_EDU ~~ SES

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
Mother_EDU ~ 1
CURRENT_ENGinputPERC ~ 1
SES ~ 1
'

lg.BPVS.COV_fit2 <- sem(lgCOV2.BPVS.lavaan_model, 
                        data = DataverySHORT_T3_NEW,
                        meanstructure = TRUE,
                        estimator = "ML",
                        missing = "fiml")
summary(lg.BPVS.COV_fit2, fit.measures=TRUE, standardized = T)
semPlot::semPaths(lg.BPVS.COV_fit2, title = FALSE,layout = "tree", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)


anova(lg.BPVS.COV_fit2,lg.BPVS.COV_fit) ##model with influence of MUM edu and SES restricted to 0 not worse than full model


### final COV model BPVS####

lgCOV.BPVS.FINAL_model <- '
# latent variable definitions
#intercept (note intercept is a reserved term)
eta_1 =~ 1*VOC_BPVS1
eta_1 =~ 1*VOC_BPVS_2
eta_1 =~ 1*VOC_BPVS_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*VOC_BPVS1
eta_2 =~ 1*VOC_BPVS_2
eta_2 =~ 2*VOC_BPVS_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2

# covariances among factors 
eta_1 ~~ eta_2

# factor means 
eta_1 ~ start(35)*1
eta_2 ~ start(4)*1

# manifest variances (made equivalent by naming theta)
VOC_BPVS1 ~~ theta*VOC_BPVS1
VOC_BPVS_2 ~~ theta*VOC_BPVS_2
VOC_BPVS_3 ~~ theta*VOC_BPVS_3

# manifest means (fixed at zero)
VOC_BPVS1 ~ 0*1
VOC_BPVS_2 ~ 0*1
VOC_BPVS_3 ~ 0*1


# modification_1#
VOC_BPVS_2 ~ VOC_BPVS1
VOC_BPVS_3 ~ VOC_BPVS_2

#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1~lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC

#covariance of TIV covariates
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
'

lg.BPVS.COVFINAL_fit <- sem(lgCOV.BPVS.FINAL_model, 
                            data = DataverySHORT_T3_NEW,
                            meanstructure = TRUE,
                            estimator = "ML",
                            missing = "fiml")

summary(lg.BPVS.COVFINAL_fit, fit.measures=TRUE, standardized = T)
fitMeasures(lg.BPVS.COVFINAL_fit)

##use fixed.x =FALSE and fixed.x =TRUE to check

semPlot::semPaths(lg.BPVS.COVFINAL_fit, title = FALSE,layout = "tree", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(lg.BPVS.COVFINAL_fit, title = FALSE,layout = "circle", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(lg.BPVS.COVFINAL_fit, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)

lgCOV.BPVS.FINALcontr_model <- '
# latent variable definitions
#intercept (note intercept is a reserved term)
eta_1 =~ 1*VOC_BPVS1
eta_1 =~ 1*VOC_BPVS_2
eta_1 =~ 1*VOC_BPVS_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*VOC_BPVS1
eta_2 =~ 1*VOC_BPVS_2
eta_2 =~ 2*VOC_BPVS_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2

# covariances among factors 
eta_1 ~~ eta_2

# factor means 
eta_1 ~ start(35)*1
eta_2 ~ start(4)*1

# manifest variances (made equivalent by naming theta)
VOC_BPVS1 ~~ theta*VOC_BPVS1
VOC_BPVS_2 ~~ theta*VOC_BPVS_2
VOC_BPVS_3 ~~ theta*VOC_BPVS_3

# manifest means (fixed at zero)
VOC_BPVS1 ~ 0*1
VOC_BPVS_2 ~ 0*1
VOC_BPVS_3 ~ 0*1


# modification_1#
VOC_BPVS_2 ~ VOC_BPVS1
VOC_BPVS_3 ~ VOC_BPVS_2

#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1~0*lenghEXPOSURE_ENG_months+0*CURRENT_ENGinputPERC
eta_2~0*lenghEXPOSURE_ENG_months+0*CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
'

lg.BPVS.COVFINALcontr_fit <- sem(lgCOV.BPVS.FINALcontr_model, 
                                 data = DataverySHORT_T3_NEW,
                                 meanstructure = TRUE,
                                 estimator = "ML",
                                 missing = "fiml")

semPlot::semPaths(lg.BPVS.COVFINALcontr_fit, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)

anova(lg.BPVS.COVFINAL_fit,lg.BPVS.COVFINALcontr_fit) ##model with INPUT1 restricted to 0 worse than model with INPUT non restricted

lgCOV.BPVS.FINALcontr_model1 <- '
# latent variable definitions
#intercept (note intercept is a reserved term)
eta_1 =~ 1*VOC_BPVS1
eta_1 =~ 1*VOC_BPVS_2
eta_1 =~ 1*VOC_BPVS_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*VOC_BPVS1
eta_2 =~ 1*VOC_BPVS_2
eta_2 =~ 2*VOC_BPVS_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2

# covariances among factors 
eta_1 ~~ eta_2

# factor means 
eta_1 ~ start(35)*1
eta_2 ~ start(4)*1

# manifest variances (made equivalent by naming theta)
VOC_BPVS1 ~~ theta*VOC_BPVS1
VOC_BPVS_2 ~~ theta*VOC_BPVS_2
VOC_BPVS_3 ~~ theta*VOC_BPVS_3

# manifest means (fixed at zero)
VOC_BPVS1 ~ 0*1
VOC_BPVS_2 ~ 0*1
VOC_BPVS_3 ~ 0*1


# modification_1#
VOC_BPVS_2 ~ VOC_BPVS1
VOC_BPVS_3 ~ VOC_BPVS_2

#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1~lenghEXPOSURE_ENG_months+0*CURRENT_ENGinputPERC
eta_2~lenghEXPOSURE_ENG_months+0*CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
'

lg.BPVS.COVFINALcontr_fit1 <- sem(lgCOV.BPVS.FINALcontr_model1, 
                                 data = DataverySHORT_T3_NEW,
                                 meanstructure = TRUE,
                                 estimator = "ML",
                                 missing = "fiml")

semPlot::semPaths(lg.BPVS.COVFINALcontr_fit1, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)

anova(lg.BPVS.COVFINAL_fit,lg.BPVS.COVFINALcontr_fit1) ##model with INPUT2 restricted to 0 worse than model with INPUT non restricted

lgCOV.BPVS.FINALcontr_model2 <- '
# latent variable definitions
#intercept (note intercept is a reserved term)
eta_1 =~ 1*VOC_BPVS1
eta_1 =~ 1*VOC_BPVS_2
eta_1 =~ 1*VOC_BPVS_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*VOC_BPVS1
eta_2 =~ 1*VOC_BPVS_2
eta_2 =~ 2*VOC_BPVS_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2

# covariances among factors 
eta_1 ~~ eta_2

# factor means 
eta_1 ~ start(35)*1
eta_2 ~ start(4)*1

# manifest variances (made equivalent by naming theta)
VOC_BPVS1 ~~ theta*VOC_BPVS1
VOC_BPVS_2 ~~ theta*VOC_BPVS_2
VOC_BPVS_3 ~~ theta*VOC_BPVS_3

# manifest means (fixed at zero)
VOC_BPVS1 ~ 0*1
VOC_BPVS_2 ~ 0*1
VOC_BPVS_3 ~ 0*1


# modification_1#
VOC_BPVS_2 ~ VOC_BPVS1
VOC_BPVS_3 ~ VOC_BPVS_2

#Time invariant covariates
regression of time-invariant covariate on intercept and slope factors
eta_1~0*lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~0*lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC

#covariance of TIV covaraites
#lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
'

lg.BPVS.COVFINALcontr_fit2 <- sem(lgCOV.BPVS.FINALcontr_model2, 
                                  data = DataverySHORT_T3_NEW,
                                  meanstructure = TRUE,
                                  estimator = "ML",
                                  missing = "fiml")

semPlot::semPaths(lg.BPVS.COVFINALcontr_fit2, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)

anova(lg.BPVS.COVFINAL_fit,lg.BPVS.COVFINALcontr_fit2) ##model with INPUT restricted to 0 worse than model with INPUT non restricted


###final graph BPVS #####

summary(lg.BPVS.COVFINAL_fit, fit.measures=TRUE, standardized = T)

semPlot::semPaths(lg.BPVS.COVFINAL_fit, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram",label.cex=1.0, edge.label.cex=1, node.width = 2, curvePivot = TRUE)
semPlot::semPaths(lg.BPVS.COVFINAL_fit, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram",nodeLabels = c("VOC BREADTH T1", "VOC BREADTH T2", "VOC BREADTH T3", "Lenght Exp","Current Input","Intercept","Slope"), label.cex=1.0, edge.label.cex=1, node.width = 2, curvePivot = TRUE)



########VOC DEPTH GROWTH###########################################

summary(DataverySHORT_T3_NEW)

###no growth model###

ng.DEPTH.lavaan_model <- '
# latent variable definitions
#intercept
eta_1 =~ 1*VOC_mean_SYNOPP1
eta_1 =~ 1*VOC_meanSYNOPP_2
eta_1 =~ 1*VOC_meanSYNOPP_3

# factor variances
eta_1 ~~ eta_1

# covariances among factors 
#none (only 1 factor)

# factor means 
eta_1 ~ start(30)*1

# manifest variances (made equivalent by naming theta)
VOC_mean_SYNOPP1 ~~ theta*VOC_mean_SYNOPP1
VOC_meanSYNOPP_2 ~~ theta*VOC_meanSYNOPP_2
VOC_meanSYNOPP_3 ~~ theta*VOC_meanSYNOPP_3

# manifest means (fixed at zero)
VOC_mean_SYNOPP1 ~ 0*1
VOC_meanSYNOPP_2 ~ 0*1
VOC_meanSYNOPP_3 ~ 0*1

# modification_1#
VOC_meanSYNOPP_2 ~ VOC_mean_SYNOPP1
VOC_meanSYNOPP_3 ~ VOC_meanSYNOPP_2

#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1~Mother_EDU+SES+Age_T1+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Mother_EDU ~~ Mother_EDU
SES ~~ SES
Age_T1 ~~ Age_T1

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ Mother_EDU
lenghEXPOSURE_ENG_months ~~ SES
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
CURRENT_ENGinputPERC ~~ Mother_EDU
CURRENT_ENGinputPERC ~~ SES
Mother_EDU ~~ SES
Age_T1 ~~ lenghEXPOSURE_ENG_months

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
Mother_EDU ~ 1
CURRENT_ENGinputPERC ~ 1
SES ~ 1
Age_T1 ~ 1

'

ng.DEPTH.lavaan_fit <- sem(ng.DEPTH.lavaan_model, 
                           data = DataverySHORT_T3_NEW,
                           meanstructure = TRUE,
                           estimator = "ML",
                           missing = "fiml")


summary(ng.DEPTH.lavaan_fit, fit.measures=TRUE, standardized = T)


summary(ng.DEPTH.lavaan_fit)
fitMeasures(ng.DEPTH.lavaan_fit)

semPlot::semPaths(ng.DEPTH.lavaan_fit, title = FALSE,layout = "tree", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(ng.DEPTH.lavaan_fit, title = FALSE,layout = "circle", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(ng.DEPTH.lavaan_fit, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)


###linear growth model###

lg.DEPTH.lavaan_model <- '
# latent variable definitions
#intercept
eta_1 =~ 1*VOC_mean_SYNOPP1
eta_1 =~ 1*VOC_meanSYNOPP_2
eta_1 =~ 1*VOC_meanSYNOPP_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*VOC_mean_SYNOPP1
eta_2 =~ 1*VOC_meanSYNOPP_2
eta_2 =~ 2*VOC_meanSYNOPP_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2

# covariances among factors 
eta_1 ~~ eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
VOC_mean_SYNOPP1 ~~ theta*VOC_mean_SYNOPP1
VOC_meanSYNOPP_2 ~~ theta*VOC_meanSYNOPP_2
VOC_meanSYNOPP_3 ~~ theta*VOC_meanSYNOPP_3

# manifest means (fixed at zero)
VOC_mean_SYNOPP1 ~ 0*1
VOC_meanSYNOPP_2 ~ 0*1
VOC_meanSYNOPP_3 ~ 0*1


# modification_1#
VOC_meanSYNOPP_2 ~ VOC_mean_SYNOPP1
VOC_meanSYNOPP_3 ~ VOC_meanSYNOPP_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~Mother_EDU+SES+Age_T1+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~Mother_EDU+SES+Age_T1+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Mother_EDU ~~ Mother_EDU
SES ~~ SES
Age_T1 ~~ Age_T1

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ Mother_EDU
lenghEXPOSURE_ENG_months ~~ SES
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
CURRENT_ENGinputPERC ~~ Mother_EDU
CURRENT_ENGinputPERC ~~ SES
Mother_EDU ~~ SES
Age_T1 ~~ lenghEXPOSURE_ENG_months

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
Mother_EDU ~ 1
CURRENT_ENGinputPERC ~ 1
SES ~ 1
Age_T1 ~ 1

'


lg.DEPTH_model_fit <- sem(lg.DEPTH.lavaan_model, 
                          data = DataverySHORT_T3_NEW,
                          meanstructure = TRUE,
                          estimator = "ML",
                          missing = "fiml")

summary(lg.DEPTH_model_fit)
summary(lg.DEPTH_model_fit, fit.measures=TRUE, standardized = T)
fitMeasures(lg.DEPTH_model_fit)

semPlot::semPaths(lg.DEPTH_model_fit, title = FALSE,layout = "tree", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(lg.DEPTH_model_fit, title = FALSE,layout = "circle", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(lg.DEPTH_model_fit, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)


anova(lg.DEPTH_model_fit,ng.DEPTH.lavaan_fit)

###quadratic growth checks ####

###quadratic growth model###

Quad_modelDEPTH <- '
#intercept
eta_1 =~ 1*VOC_mean_SYNOPP1
eta_1 =~ 1*VOC_meanSYNOPP_2
eta_1 =~ 1*VOC_meanSYNOPP_3

#linear change factor loadings
eta_2 =~ -0.5*VOC_mean_SYNOPP1
eta_2 =~ 0*VOC_meanSYNOPP_2
eta_2 =~ 0.5*VOC_meanSYNOPP_3

#quadratic change factor loadings
eta_3 =~ 0.25*VOC_mean_SYNOPP1
eta_3 =~ 0*VOC_meanSYNOPP_2
eta_3 =~ 0.25*VOC_meanSYNOPP_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2
eta_3 ~~ eta_3

# covariances among factors 
eta_1 ~~ eta_2
eta_3 ~~ eta_2
eta_1 ~~ eta_3

# manifest variances (made equivalent by naming theta)
VOC_mean_SYNOPP1 ~~ theta*VOC_mean_SYNOPP1
VOC_meanSYNOPP_2 ~~ theta*VOC_meanSYNOPP_2
VOC_meanSYNOPP_3 ~~ theta*VOC_meanSYNOPP_3

# manifest means (fixed at zero)
VOC_mean_SYNOPP1 ~ 0*1
VOC_meanSYNOPP_2 ~ 0*1
VOC_meanSYNOPP_3 ~ 0*1


# modification_1#
VOC_meanSYNOPP_2 ~ VOC_mean_SYNOPP1
VOC_meanSYNOPP_3 ~ VOC_meanSYNOPP_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~Mother_EDU+SES+Age_T1+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~Mother_EDU+SES+Age_T1+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Mother_EDU ~~ Mother_EDU
SES ~~ SES
Age_T1 ~~ Age_T1

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ Mother_EDU
lenghEXPOSURE_ENG_months ~~ SES
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
CURRENT_ENGinputPERC ~~ Mother_EDU
CURRENT_ENGinputPERC ~~ SES
Mother_EDU ~~ SES
Age_T1 ~~ lenghEXPOSURE_ENG_months

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
Mother_EDU ~ 1
CURRENT_ENGinputPERC ~ 1
SES ~ 1
Age_T1 ~ 1

'
quad_fitDEPTH <- sem(Quad_modelDEPTH, data = DataverySHORT_T3_NEW,
                meanstructure = TRUE,
                estimator = "ML",
                missing = "fiml")
summary(quad_fitDEPTH, fit.measures = T, standardized = T)

#change for model with var-cov not positive defined

Quad_modelDEPTH2 <- '
#intercept
eta_1 =~ 1*VOC_mean_SYNOPP1
eta_1 =~ 1*VOC_meanSYNOPP_2
eta_1 =~ 1*VOC_meanSYNOPP_3

#linear change factor loadings
eta_2 =~ -0.5*VOC_mean_SYNOPP1
eta_2 =~ 0*VOC_meanSYNOPP_2
eta_2 =~ 0.5*VOC_meanSYNOPP_3

#quadratic change factor loadings
eta_3 =~ 0.25*VOC_mean_SYNOPP1
eta_3 =~ 0*VOC_meanSYNOPP_2
eta_3 =~ 0.25*VOC_meanSYNOPP_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2
eta_3 ~~ 1*eta_3

# covariances among factors 
eta_1 ~~ eta_2
eta_3 ~~ 0*eta_2
eta_1 ~~ 0*eta_3

# manifest variances (made equivalent by naming theta)
VOC_mean_SYNOPP1 ~~ theta*VOC_mean_SYNOPP1
VOC_meanSYNOPP_2 ~~ theta*VOC_meanSYNOPP_2
VOC_meanSYNOPP_3 ~~ theta*VOC_meanSYNOPP_3

# manifest means (fixed at zero)
VOC_mean_SYNOPP1 ~ 0*1
VOC_meanSYNOPP_2 ~ 0*1
VOC_meanSYNOPP_3 ~ 0*1


# modification_1#
VOC_meanSYNOPP_2 ~ VOC_mean_SYNOPP1
VOC_meanSYNOPP_3 ~ VOC_meanSYNOPP_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~0*Mother_EDU+0*SES+0*Age_T1+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~0*Mother_EDU+0*SES+0*Age_T1+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Mother_EDU ~~ Mother_EDU
SES ~~ SES
Age_T1 ~~ Age_T1

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ Mother_EDU
lenghEXPOSURE_ENG_months ~~ SES
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
CURRENT_ENGinputPERC ~~ Mother_EDU
CURRENT_ENGinputPERC ~~ SES
Mother_EDU ~~ SES
Age_T1 ~~ lenghEXPOSURE_ENG_months

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
Mother_EDU ~ 1
CURRENT_ENGinputPERC ~ 1
SES ~ 1
Age_T1 ~ 1

'
quad_fitDEPTH2 <- sem(Quad_modelDEPTH2, data = DataverySHORT_T3_NEW,
                     meanstructure = TRUE,
                     estimator = "ML",
                     missing = "fiml")
summary(quad_fitDEPTH2, fit.measures = T, standardized = T) ###some estimated lv variances are negative
##check for quadratic growth with final model because of model incorrect estimation


Quad_modelDEPTH4 <- '
#intercept
eta_1 =~ 1*VOC_mean_SYNOPP1
eta_1 =~ 1*VOC_meanSYNOPP_2
eta_1 =~ 1*VOC_meanSYNOPP_3

#quadratic change factor loadings
eta_3 =~ 0.25*VOC_mean_SYNOPP1
eta_3 =~ 0*VOC_meanSYNOPP_2
eta_3 =~ 0.25*VOC_meanSYNOPP_3

# factor variances
eta_1 ~~ eta_1
eta_3 ~~ eta_3

# covariances among factors 
eta_1 ~~ eta_3

# manifest variances (made equivalent by naming theta)
VOC_mean_SYNOPP1 ~~ theta*VOC_mean_SYNOPP1
VOC_meanSYNOPP_2 ~~ theta*VOC_meanSYNOPP_2
VOC_meanSYNOPP_3 ~~ theta*VOC_meanSYNOPP_3

# manifest means (fixed at zero)
VOC_mean_SYNOPP1 ~ 0*1
VOC_meanSYNOPP_2 ~ 0*1
VOC_meanSYNOPP_3 ~ 0*1


# modification_1#
VOC_meanSYNOPP_2 ~ VOC_mean_SYNOPP1
VOC_meanSYNOPP_3 ~ VOC_meanSYNOPP_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~0*Mother_EDU+0*SES+0*Age_T1+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_3~0*Mother_EDU+0*SES+0*Age_T1+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Mother_EDU ~~ Mother_EDU
SES ~~ SES
Age_T1 ~~ Age_T1

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ Mother_EDU
lenghEXPOSURE_ENG_months ~~ SES
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
CURRENT_ENGinputPERC ~~ Mother_EDU
CURRENT_ENGinputPERC ~~ SES
Mother_EDU ~~ SES
Age_T1 ~~ lenghEXPOSURE_ENG_months

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
Mother_EDU ~ 1
CURRENT_ENGinputPERC ~ 1
SES ~ 1
Age_T1 ~ 1

'
quad_fitDEPTH4 <- sem(Quad_modelDEPTH4, data = DataverySHORT_T3_NEW,
                      meanstructure = TRUE,
                      estimator = "ML",
                      missing = "fiml")
summary(quad_fitDEPTH4, fit.measures = T, standardized = T) ###some estimated lv variances are negative

anova(quad_fitDEPTH4,lg.DEPTH_model_fit)

##covariates

lg.DEPTH.lavaan_model <- '
# latent variable definitions
#intercept
eta_1 =~ 1*VOC_mean_SYNOPP1
eta_1 =~ 1*VOC_meanSYNOPP_2
eta_1 =~ 1*VOC_meanSYNOPP_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*VOC_mean_SYNOPP1
eta_2 =~ 1*VOC_meanSYNOPP_2
eta_2 =~ 2*VOC_meanSYNOPP_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2

# covariances among factors 
eta_1 ~~ eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
VOC_mean_SYNOPP1 ~~ theta*VOC_mean_SYNOPP1
VOC_meanSYNOPP_2 ~~ theta*VOC_meanSYNOPP_2
VOC_meanSYNOPP_3 ~~ theta*VOC_meanSYNOPP_3

# manifest means (fixed at zero)
VOC_mean_SYNOPP1 ~ 0*1
VOC_meanSYNOPP_2 ~ 0*1
VOC_meanSYNOPP_3 ~ 0*1


# modification_1#
VOC_meanSYNOPP_2 ~ VOC_mean_SYNOPP1
VOC_meanSYNOPP_3 ~ VOC_meanSYNOPP_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~Mother_EDU+SES+Age_T1+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~Mother_EDU+SES+Age_T1+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Mother_EDU ~~ Mother_EDU
SES ~~ SES
Age_T1 ~~ Age_T1

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ Mother_EDU
lenghEXPOSURE_ENG_months ~~ SES
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
CURRENT_ENGinputPERC ~~ Mother_EDU
CURRENT_ENGinputPERC ~~ SES
Mother_EDU ~~ SES
Age_T1 ~~ lenghEXPOSURE_ENG_months

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
Mother_EDU ~ 1
CURRENT_ENGinputPERC ~ 1
SES ~ 1
Age_T1 ~ 1

'


lg.DEPTH_model_fit <- sem(lg.DEPTH.lavaan_model, 
                          data = DataverySHORT_T3_NEW,
                          meanstructure = TRUE,
                          estimator = "ML",
                          missing = "fiml")

summary(lg.DEPTH_model_fit)
summary(lg.DEPTH_model_fit, fit.measures=TRUE, standardized = T)

lg2.DEPTH.lavaan_model <- '
# latent variable definitions
#intercept
eta_1 =~ 1*VOC_mean_SYNOPP1
eta_1 =~ 1*VOC_meanSYNOPP_2
eta_1 =~ 1*VOC_meanSYNOPP_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*VOC_mean_SYNOPP1
eta_2 =~ 1*VOC_meanSYNOPP_2
eta_2 =~ 2*VOC_meanSYNOPP_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2

# covariances among factors 
eta_1 ~~ eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
VOC_mean_SYNOPP1 ~~ theta*VOC_mean_SYNOPP1
VOC_meanSYNOPP_2 ~~ theta*VOC_meanSYNOPP_2
VOC_meanSYNOPP_3 ~~ theta*VOC_meanSYNOPP_3

# manifest means (fixed at zero)
VOC_mean_SYNOPP1 ~ 0*1
VOC_meanSYNOPP_2 ~ 0*1
VOC_meanSYNOPP_3 ~ 0*1


# modification_1#
VOC_meanSYNOPP_2 ~ VOC_mean_SYNOPP1
VOC_meanSYNOPP_3 ~ VOC_meanSYNOPP_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~0*Mother_EDU+0*SES+Age_T1+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~0*Mother_EDU+0*SES+Age_T1+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Mother_EDU ~~ Mother_EDU
SES ~~ SES
Age_T1 ~~ Age_T1

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ Mother_EDU
lenghEXPOSURE_ENG_months ~~ SES
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
CURRENT_ENGinputPERC ~~ Mother_EDU
CURRENT_ENGinputPERC ~~ SES
Mother_EDU ~~ SES
Age_T1 ~~ lenghEXPOSURE_ENG_months

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
Mother_EDU ~ 1
CURRENT_ENGinputPERC ~ 1
SES ~ 1
Age_T1 ~ 1

'


lg2.DEPTH_model_fit <- sem(lg2.DEPTH.lavaan_model, 
                          data = DataverySHORT_T3_NEW,
                          meanstructure = TRUE,
                          estimator = "ML",
                          missing = "fiml")

summary(lg2.DEPTH_model_fit)
summary(lg2.DEPTH_model_fit, fit.measures=TRUE, standardized = T)

anova(lg.DEPTH_model_fit,lg2.DEPTH_model_fit)

###final lg model for VOC DEPTH###

lgFINAL.DEPTH.lavaan_model <- '
# latent variable definitions
#intercept
eta_1 =~ 1*VOC_mean_SYNOPP1
eta_1 =~ 1*VOC_meanSYNOPP_2
eta_1 =~ 1*VOC_meanSYNOPP_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*VOC_mean_SYNOPP1
eta_2 =~ 1*VOC_meanSYNOPP_2
eta_2 =~ 2*VOC_meanSYNOPP_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2

# covariances among factors 
eta_1 ~~ eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
VOC_mean_SYNOPP1 ~~ theta*VOC_mean_SYNOPP1
VOC_meanSYNOPP_2 ~~ theta*VOC_meanSYNOPP_2
VOC_meanSYNOPP_3 ~~ theta*VOC_meanSYNOPP_3

# manifest means (fixed at zero)
VOC_mean_SYNOPP1 ~ 0*1
VOC_meanSYNOPP_2 ~ 0*1
VOC_meanSYNOPP_3 ~ 0*1


# modification_1#
VOC_meanSYNOPP_2 ~ VOC_mean_SYNOPP1
VOC_meanSYNOPP_3 ~ VOC_meanSYNOPP_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~Age_T1+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~Age_T1+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Age_T1 ~~ Age_T1

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
Age_T1 ~~ lenghEXPOSURE_ENG_months

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
Age_T1 ~ 1

'


lg.DEPTH_model_fitFINAL <- sem(lgFINAL.DEPTH.lavaan_model, 
                           data = DataverySHORT_T3_NEW,
                           meanstructure = TRUE,
                           estimator = "ML",
                           missing = "fiml")

summary(lg.DEPTH_model_fitFINAL)
summary(lg.DEPTH_model_fitFINAL, fit.measures=TRUE, standardized = T)

lgFINAL.DEPTH.lavaan_model_CONTROLage <- '
# latent variable definitions
#intercept
eta_1 =~ 1*VOC_mean_SYNOPP1
eta_1 =~ 1*VOC_meanSYNOPP_2
eta_1 =~ 1*VOC_meanSYNOPP_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*VOC_mean_SYNOPP1
eta_2 =~ 1*VOC_meanSYNOPP_2
eta_2 =~ 2*VOC_meanSYNOPP_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2

# covariances among factors 
eta_1 ~~ eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
VOC_mean_SYNOPP1 ~~ theta*VOC_mean_SYNOPP1
VOC_meanSYNOPP_2 ~~ theta*VOC_meanSYNOPP_2
VOC_meanSYNOPP_3 ~~ theta*VOC_meanSYNOPP_3

# manifest means (fixed at zero)
VOC_mean_SYNOPP1 ~ 0*1
VOC_meanSYNOPP_2 ~ 0*1
VOC_meanSYNOPP_3 ~ 0*1


# modification_1#
VOC_meanSYNOPP_2 ~ VOC_mean_SYNOPP1
VOC_meanSYNOPP_3 ~ VOC_meanSYNOPP_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~0*Age_T1+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~0*Age_T1+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Age_T1 ~~ Age_T1

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
Age_T1 ~~ lenghEXPOSURE_ENG_months

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
Age_T1 ~ 1

'


lg.DEPTH_model_fitFINAL_controlAGE <- sem(lgFINAL.DEPTH.lavaan_model_CONTROLage, 
                               data = DataverySHORT_T3_NEW,
                               meanstructure = TRUE,
                               estimator = "ML",
                               missing = "fiml")

summary(lg.DEPTH_model_fitFINAL_controlAGE)
summary(lg.DEPTH_model_fitFINAL_controlAGE, fit.measures=TRUE, standardized = T)
anova(lg.DEPTH_model_fitFINAL,lg.DEPTH_model_fitFINAL_controlAGE)

lgFINAL.DEPTH.lavaan_model_CONTROLinput1 <- '
# latent variable definitions
#intercept
eta_1 =~ 1*VOC_mean_SYNOPP1
eta_1 =~ 1*VOC_meanSYNOPP_2
eta_1 =~ 1*VOC_meanSYNOPP_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*VOC_mean_SYNOPP1
eta_2 =~ 1*VOC_meanSYNOPP_2
eta_2 =~ 2*VOC_meanSYNOPP_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2

# covariances among factors 
eta_1 ~~ eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
VOC_mean_SYNOPP1 ~~ theta*VOC_mean_SYNOPP1
VOC_meanSYNOPP_2 ~~ theta*VOC_meanSYNOPP_2
VOC_meanSYNOPP_3 ~~ theta*VOC_meanSYNOPP_3

# manifest means (fixed at zero)
VOC_mean_SYNOPP1 ~ 0*1
VOC_meanSYNOPP_2 ~ 0*1
VOC_meanSYNOPP_3 ~ 0*1


# modification_1#
VOC_meanSYNOPP_2 ~ VOC_mean_SYNOPP1
VOC_meanSYNOPP_3 ~ VOC_meanSYNOPP_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~Age_T1+0*lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~Age_T1+0*lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Age_T1 ~~ Age_T1

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
Age_T1 ~~ lenghEXPOSURE_ENG_months

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
Age_T1 ~ 1

'


lg.DEPTH_model_fitFINAL_controlINPUT1 <- sem(lgFINAL.DEPTH.lavaan_model_CONTROLinput1, 
                               data = DataverySHORT_T3_NEW,
                               meanstructure = TRUE,
                               estimator = "ML",
                               missing = "fiml")

summary(lg.DEPTH_model_fitFINAL_controlINPUT1)
summary(lg.DEPTH_model_fitFINAL_controlINPUT1, fit.measures=TRUE, standardized = T)
anova(lg.DEPTH_model_fitFINAL,lg.DEPTH_model_fitFINAL_controlINPUT1)

lgFINAL.DEPTH.lavaan_model_CONTROLinput2 <- '
# latent variable definitions
#intercept
eta_1 =~ 1*VOC_mean_SYNOPP1
eta_1 =~ 1*VOC_meanSYNOPP_2
eta_1 =~ 1*VOC_meanSYNOPP_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*VOC_mean_SYNOPP1
eta_2 =~ 1*VOC_meanSYNOPP_2
eta_2 =~ 2*VOC_meanSYNOPP_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2

# covariances among factors 
eta_1 ~~ eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
VOC_mean_SYNOPP1 ~~ theta*VOC_mean_SYNOPP1
VOC_meanSYNOPP_2 ~~ theta*VOC_meanSYNOPP_2
VOC_meanSYNOPP_3 ~~ theta*VOC_meanSYNOPP_3

# manifest means (fixed at zero)
VOC_mean_SYNOPP1 ~ 0*1
VOC_meanSYNOPP_2 ~ 0*1
VOC_meanSYNOPP_3 ~ 0*1


# modification_1#
VOC_meanSYNOPP_2 ~ VOC_mean_SYNOPP1
VOC_meanSYNOPP_3 ~ VOC_meanSYNOPP_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~Age_T1+lenghEXPOSURE_ENG_months+0*CURRENT_ENGinputPERC
eta_2~Age_T1+lenghEXPOSURE_ENG_months+0*CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Age_T1 ~~ Age_T1

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
Age_T1 ~~ lenghEXPOSURE_ENG_months

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
Age_T1 ~ 1

'


lg.DEPTH_model_fitFINAL_controlINPUT2 <- sem(lgFINAL.DEPTH.lavaan_model_CONTROLinput2, 
                                             data = DataverySHORT_T3_NEW,
                                             meanstructure = TRUE,
                                             estimator = "ML",
                                             missing = "fiml")

summary(lg.DEPTH_model_fitFINAL_controlINPUT2)
summary(lg.DEPTH_model_fitFINAL_controlINPUT2, fit.measures=TRUE, standardized = T)
anova(lg.DEPTH_model_fitFINAL,lg.DEPTH_model_fitFINAL_controlINPUT2)

###final plots###

lgFINAL.DEPTH.lavaan_model_NEW <- '
# latent variable definitions
#intercept
eta_1 =~ 1*VOC_mean_SYNOPP1
eta_1 =~ 1*VOC_meanSYNOPP_2
eta_1 =~ 1*VOC_meanSYNOPP_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*VOC_mean_SYNOPP1
eta_2 =~ 1*VOC_meanSYNOPP_2
eta_2 =~ 2*VOC_meanSYNOPP_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2

# covariances among factors 
eta_1 ~~ eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
VOC_mean_SYNOPP1 ~~ theta*VOC_mean_SYNOPP1
VOC_meanSYNOPP_2 ~~ theta*VOC_meanSYNOPP_2
VOC_meanSYNOPP_3 ~~ theta*VOC_meanSYNOPP_3

# manifest means (fixed at zero)
VOC_mean_SYNOPP1 ~ 0*1
VOC_meanSYNOPP_2 ~ 0*1
VOC_meanSYNOPP_3 ~ 0*1


# modification_1#
VOC_meanSYNOPP_2 ~ VOC_mean_SYNOPP1
VOC_meanSYNOPP_3 ~ VOC_meanSYNOPP_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~Age_T1+CURRENT_ENGinputPERC
eta_2~Age_T1+CURRENT_ENGinputPERC

#variance of TIV covariates
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Age_T1 ~~ Age_T1

#means of TIV covariates (freely estimated)
CURRENT_ENGinputPERC ~ 1
Age_T1 ~ 1

'


lg.DEPTH_model_fitFINAL_NEW <- sem(lgFINAL.DEPTH.lavaan_model_NEW, 
                               data = DataverySHORT_T3_NEW,
                               meanstructure = TRUE,
                               estimator = "ML",
                               missing = "fiml")

summary(lg.DEPTH_model_fitFINAL_NEW)
summary(lg.DEPTH_model_fitFINAL_NEW, fit.measures=TRUE, standardized = T)

###quad model for VOC DEPTH#
quad.DEPTH.lavaan_modelNEWNEW <- '
# latent variable definitions
#intercept
eta_1 =~ 1*VOC_mean_SYNOPP1
eta_1 =~ 1*VOC_meanSYNOPP_2
eta_1 =~ 1*VOC_meanSYNOPP_3

#linear slope (note intercept is a reserved term)
eta_2 =~ -0.5*VOC_mean_SYNOPP1
eta_2 =~ 0*VOC_meanSYNOPP_2
eta_2 =~ 0.5*VOC_meanSYNOPP_3

#quad slope (note intercept is a reserved term)
eta_3 =~ 0.5*VOC_mean_SYNOPP1
eta_3 =~ 0*VOC_meanSYNOPP_2
eta_3 =~ 0.25*VOC_meanSYNOPP_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2
eta_3 ~~ eta_3

# covariances among factors 
eta_1 ~~ eta_2
eta_1 ~~ eta_3
eta_3 ~~ eta_2

# factor means 
eta_1 ~ start(35)*1
eta_2 ~ start(4)*1

# manifest variances (made equivalent by naming theta)
VOC_mean_SYNOPP1 ~~ theta*VOC_mean_SYNOPP1
VOC_meanSYNOPP_2 ~~ theta*VOC_meanSYNOPP_2
VOC_meanSYNOPP_3 ~~ theta*VOC_meanSYNOPP_3

# manifest means (fixed at zero)
VOC_mean_SYNOPP1 ~ 0*1
VOC_meanSYNOPP_2 ~ 0*1
VOC_meanSYNOPP_3 ~ 0*1


# modification_1#
VOC_meanSYNOPP_2 ~ VOC_mean_SYNOPP1
VOC_meanSYNOPP_3 ~ VOC_meanSYNOPP_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~Age_T1+CURRENT_ENGinputPERC
eta_2~Age_T1+CURRENT_ENGinputPERC
#eta_3~Age_T1+CURRENT_ENGinputPERC

#variance of TIV covariates
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Age_T1 ~~ Age_T1

#covariance of TIV covaraites
CURRENT_ENGinputPERC ~~ Age_T1

#means of TIV covariates (freely estimated)
CURRENT_ENGinputPERC ~ 1
Age_T1 ~ 1

'


quad.DEPTH_model_fitNEWNEW <- sem(quad.DEPTH.lavaan_modelNEWNEW, 
                                  data = DataverySHORT_T3_NEW,
                                  meanstructure = TRUE,
                                  estimator = "ML",
                                  missing = "fiml")

summary(quad.DEPTH_model_fitNEWNEW, fit.measures=TRUE, standardized = T)


anova(lg.DEPTH_model_fitFINAL_NEW, quad.DEPTH_model_fitNEWNEW)

###final model VOC DEPTH###

summary(lg.DEPTH_model_fitFINAL_NEW, fit.measures=TRUE, standardized = T)


semPlot::semPaths(lg.DEPTH_model_fitFINAL_NEW, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 2, curvePivot = TRUE)
semPlot::semPaths(lg.DEPTH_model_fitFINAL_NEW, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram",nodeLabels = c("VOC DEPTH T1", "VOC DEPTH T2", "VOC DEPTH T3", "Age T1", "Current Input","Intercept","Slope"), label.cex=1.0, edge.label.cex=1, node.width = 2, curvePivot = TRUE)


###GRAMMAR WS######################################################################################################################


### no growth#####

summary(DataverySHORT_T3_NEW)



ng.WS.lavaan_model <- '
# latent variable definitions
#intercept
eta_1 =~ 1*GRAMM_WordStruct1
eta_1 =~ 1*GRAMM_WordStruct_2
eta_1 =~ 1*GRAMM_WordStruct_3

# factor variances
eta_1 ~~ eta_1

# covariances among factors 
#none (only 1 factor)

# factor means 
eta_1 ~ start(5)*1

# manifest variances (made equivalent by naming theta)
GRAMM_WordStruct1 ~~ theta*GRAMM_WordStruct1
GRAMM_WordStruct_2 ~~ theta*GRAMM_WordStruct_2
GRAMM_WordStruct_3 ~~ theta*GRAMM_WordStruct_3

# manifest means (fixed at zero)
GRAMM_WordStruct1 ~ 0*1
GRAMM_WordStruct_2 ~ 0*1
GRAMM_WordStruct_3 ~ 0*1

# modification_1#
GRAMM_WordStruct_2 ~ GRAMM_WordStruct1
GRAMM_WordStruct_3 ~ GRAMM_WordStruct_2

#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1~Mother_EDU+SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Mother_EDU ~~ Mother_EDU
SES ~~ SES

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ Mother_EDU
lenghEXPOSURE_ENG_months ~~ SES
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
CURRENT_ENGinputPERC ~~ Mother_EDU
CURRENT_ENGinputPERC ~~ SES
Mother_EDU ~~ SES

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
Mother_EDU ~ 1
CURRENT_ENGinputPERC ~ 1
SES ~ 1
'
ng.WS.lavaan_fit <- sem(ng.WS.lavaan_model, 
                          data = DataverySHORT_T3_NEW,
                          meanstructure = TRUE,
                          estimator = "ML",
                          missing = "fiml")


summary(ng.WS.lavaan_fit, fit.measures=TRUE, standardized = T)


summary(ng.WS.lavaan_fit)
fitMeasures(ng.BPVS.lavaan_fit)

semPlot::semPaths(ng.WS.lavaan_fit, title = FALSE,layout = "tree", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(ng.WS.lavaan_fit, title = FALSE,layout = "circle", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(ng.WS.lavaan_fit, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)

###linear growth model###

lg.WS.lavaan_model <- '
# latent variable definitions
#intercept
eta_1 =~ 1*GRAMM_WordStruct1
eta_1 =~ 1*GRAMM_WordStruct_2
eta_1 =~ 1*GRAMM_WordStruct_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*GRAMM_WordStruct1
eta_2 =~ 1*GRAMM_WordStruct_2
eta_2 =~ 2*GRAMM_WordStruct_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2

# covariances among factors 
eta_1 ~~ eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
GRAMM_WordStruct1 ~~ theta*GRAMM_WordStruct1
GRAMM_WordStruct_2 ~~ theta*GRAMM_WordStruct_2
GRAMM_WordStruct_3 ~~ theta*GRAMM_WordStruct_3

# manifest means (fixed at zero)
GRAMM_WordStruct1 ~ 0*1
GRAMM_WordStruct_2 ~ 0*1
GRAMM_WordStruct_3 ~ 0*1


# modification_1#
GRAMM_WordStruct_2 ~ GRAMM_WordStruct1
GRAMM_WordStruct_3 ~ GRAMM_WordStruct_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~Mother_EDU+SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~Mother_EDU+SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Mother_EDU ~~ Mother_EDU
SES ~~ SES

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ Mother_EDU
lenghEXPOSURE_ENG_months ~~ SES
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
CURRENT_ENGinputPERC ~~ Mother_EDU
CURRENT_ENGinputPERC ~~ SES
Mother_EDU ~~ SES

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
Mother_EDU ~ 1
CURRENT_ENGinputPERC ~ 1
SES ~ 1
'


lg.WS_model_fit <- sem(lg.WS.lavaan_model, 
                          data = DataverySHORT_T3_NEW,
                          meanstructure = TRUE,
                          estimator = "ML",
                          missing = "fiml")

summary(lg.WS_model_fit)
summary(lg.WS_model_fit, fit.measures=TRUE, standardized = T)
fitMeasures(lg.WS_model_fit)

semPlot::semPaths(lg.WS_model_fit, title = FALSE,layout = "tree", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(lg.WS_model_fit, title = FALSE,layout = "circle", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(lg.WS_model_fit, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)


anova(lg.WS_model_fit,ng.WS.lavaan_fit)


###Quadratic growth###

QUAD.WS.lavaan_model <- '
# latent variable definitions
#intercept
eta_1 =~ 1*GRAMM_WordStruct1
eta_1 =~ 1*GRAMM_WordStruct_2
eta_1 =~ 1*GRAMM_WordStruct_3

#linear slope (note intercept is a reserved term)
eta_2 =~ -0.5*GRAMM_WordStruct1
eta_2 =~ 0*GRAMM_WordStruct_2
eta_2 =~ 0.5*GRAMM_WordStruct_3

#quad slope (note intercept is a reserved term)
eta_3 =~ 0.5*GRAMM_WordStruct1
eta_3 =~ 0*GRAMM_WordStruct_2
eta_3 =~ 0.25*GRAMM_WordStruct_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2
eta_3 ~~ eta_3

# covariances among factors 
eta_1 ~~ eta_2
eta_1 ~~ eta_3
eta_3 ~~ eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
GRAMM_WordStruct1 ~~ theta*GRAMM_WordStruct1
GRAMM_WordStruct_2 ~~ theta*GRAMM_WordStruct_2
GRAMM_WordStruct_3 ~~ theta*GRAMM_WordStruct_3

# manifest means (fixed at zero)
GRAMM_WordStruct1 ~ 0*1
GRAMM_WordStruct_2 ~ 0*1
GRAMM_WordStruct_3 ~ 0*1


# modification_1#
GRAMM_WordStruct_2 ~ GRAMM_WordStruct1
GRAMM_WordStruct_3 ~ GRAMM_WordStruct_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~Mother_EDU+SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~Mother_EDU+SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Mother_EDU ~~ Mother_EDU
SES ~~ SES

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ Mother_EDU
lenghEXPOSURE_ENG_months ~~ SES
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
CURRENT_ENGinputPERC ~~ Mother_EDU
CURRENT_ENGinputPERC ~~ SES
Mother_EDU ~~ SES

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
Mother_EDU ~ 1
CURRENT_ENGinputPERC ~ 1
SES ~ 1
'
QUAD.WS_model_fit <- sem(QUAD.WS.lavaan_model, 
                       data = DataverySHORT_T3_NEW,
                       meanstructure = TRUE,
                       estimator = "ML",
                       missing = "fiml")
summary(QUAD.WS_model_fit, fit.measures=TRUE, standardized = T)

QUAD.WS.lavaan_model2 <- '
# latent variable definitions
#intercept
eta_1 =~ 1*GRAMM_WordStruct1
eta_1 =~ 1*GRAMM_WordStruct_2
eta_1 =~ 1*GRAMM_WordStruct_3

#linear slope (note intercept is a reserved term)
eta_2 =~ -0.5*GRAMM_WordStruct1
eta_2 =~ 0*GRAMM_WordStruct_2
eta_2 =~ 0.5*GRAMM_WordStruct_3

#quad slope (note intercept is a reserved term)
eta_3 =~ 0.5*GRAMM_WordStruct1
eta_3 =~ 0*GRAMM_WordStruct_2
eta_3 =~ 0.25*GRAMM_WordStruct_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2
eta_3 ~~ eta_3

# covariances among factors 
eta_1 ~~ eta_2
eta_1 ~~ 0*eta_3
eta_3 ~~ 0*eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
GRAMM_WordStruct1 ~~ theta*GRAMM_WordStruct1
GRAMM_WordStruct_2 ~~ theta*GRAMM_WordStruct_2
GRAMM_WordStruct_3 ~~ theta*GRAMM_WordStruct_3

# manifest means (fixed at zero)
GRAMM_WordStruct1 ~ 0*1
GRAMM_WordStruct_2 ~ 0*1
GRAMM_WordStruct_3 ~ 0*1


# modification_1#
GRAMM_WordStruct_2 ~ GRAMM_WordStruct1
GRAMM_WordStruct_3 ~ GRAMM_WordStruct_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~Mother_EDU+SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~Mother_EDU+SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Mother_EDU ~~ Mother_EDU
SES ~~ SES

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ Mother_EDU
lenghEXPOSURE_ENG_months ~~ SES
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
CURRENT_ENGinputPERC ~~ Mother_EDU
CURRENT_ENGinputPERC ~~ SES
Mother_EDU ~~ SES

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
Mother_EDU ~ 1
CURRENT_ENGinputPERC ~ 1
SES ~ 1
'
QUAD.WS_model_fit2 <- sem(QUAD.WS.lavaan_model2, 
                         data = DataverySHORT_T3_NEW,
                         meanstructure = TRUE,
                         estimator = "ML",
                         missing = "fiml")
summary(QUAD.WS_model_fit2, fit.measures=TRUE, standardized = T)
anova(lg.WS_model_fit,QUAD.WS_model_fit2)

###covariates###

lg2.WS.lavaan_model <- '
# latent variable definitions
#intercept
eta_1 =~ 1*GRAMM_WordStruct1
eta_1 =~ 1*GRAMM_WordStruct_2
eta_1 =~ 1*GRAMM_WordStruct_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*GRAMM_WordStruct1
eta_2 =~ 1*GRAMM_WordStruct_2
eta_2 =~ 2*GRAMM_WordStruct_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2

# covariances among factors 
eta_1 ~~ eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
GRAMM_WordStruct1 ~~ theta*GRAMM_WordStruct1
GRAMM_WordStruct_2 ~~ theta*GRAMM_WordStruct_2
GRAMM_WordStruct_3 ~~ theta*GRAMM_WordStruct_3

# manifest means (fixed at zero)
GRAMM_WordStruct1 ~ 0*1
GRAMM_WordStruct_2 ~ 0*1
GRAMM_WordStruct_3 ~ 0*1


# modification_1#
GRAMM_WordStruct_2 ~ GRAMM_WordStruct1
GRAMM_WordStruct_3 ~ GRAMM_WordStruct_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~0*Mother_EDU+0*SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~0*Mother_EDU+0*SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Mother_EDU ~~ Mother_EDU
SES ~~ SES

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ Mother_EDU
lenghEXPOSURE_ENG_months ~~ SES
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
CURRENT_ENGinputPERC ~~ Mother_EDU
CURRENT_ENGinputPERC ~~ SES
Mother_EDU ~~ SES

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
Mother_EDU ~ 1
CURRENT_ENGinputPERC ~ 1
SES ~ 1
'


lg2.WS_model_fit <- sem(lg2.WS.lavaan_model, 
                       data = DataverySHORT_T3_NEW,
                       meanstructure = TRUE,
                       estimator = "ML",
                       missing = "fiml")

summary(lg2.WS_model_fit)
summary(lg2.WS_model_fit, fit.measures=TRUE, standardized = T)
fitMeasures(lg2.WS_model_fit)

semPlot::semPaths(lg2.WS_model_fit, title = FALSE,layout = "tree", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(lg2.WS_model_fit, title = FALSE,layout = "circle", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(lg2.WS_model_fit, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)


anova(lg.WS_model_fit,lg2.WS_model_fit)

lg.WS.lavaan_modelINPUT1 <- '
# latent variable definitions
#intercept
eta_1 =~ 1*GRAMM_WordStruct1
eta_1 =~ 1*GRAMM_WordStruct_2
eta_1 =~ 1*GRAMM_WordStruct_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*GRAMM_WordStruct1
eta_2 =~ 1*GRAMM_WordStruct_2
eta_2 =~ 2*GRAMM_WordStruct_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2

# covariances among factors 
eta_1 ~~ eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
GRAMM_WordStruct1 ~~ theta*GRAMM_WordStruct1
GRAMM_WordStruct_2 ~~ theta*GRAMM_WordStruct_2
GRAMM_WordStruct_3 ~~ theta*GRAMM_WordStruct_3

# manifest means (fixed at zero)
GRAMM_WordStruct1 ~ 0*1
GRAMM_WordStruct_2 ~ 0*1
GRAMM_WordStruct_3 ~ 0*1


# modification_1#
GRAMM_WordStruct_2 ~ GRAMM_WordStruct1
GRAMM_WordStruct_3 ~ GRAMM_WordStruct_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC


#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
'


lg.WS.lavaan_modelINPUT1_fit <- sem(lg.WS.lavaan_modelINPUT1, 
                        data = DataverySHORT_T3_NEW,
                        meanstructure = TRUE,
                        estimator = "ML",
                        missing = "fiml")

summary(lg.WS.lavaan_modelINPUT1_fit)
summary(lg.WS.lavaan_modelINPUT1_fit, fit.measures=TRUE, standardized = T)
fitMeasures(lg.WS.lavaan_modelINPUT1_fit)

semPlot::semPaths(lg.WS.lavaan_modelINPUT1_fit, title = FALSE,layout = "tree", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(lg.WS.lavaan_modelINPUT1_fit, title = FALSE,layout = "circle", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(lg.WS.lavaan_modelINPUT1_fit, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)


lg.WS.lavaan_modelINPUT2 <- '
# latent variable definitions
#intercept
eta_1 =~ 1*GRAMM_WordStruct1
eta_1 =~ 1*GRAMM_WordStruct_2
eta_1 =~ 1*GRAMM_WordStruct_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*GRAMM_WordStruct1
eta_2 =~ 1*GRAMM_WordStruct_2
eta_2 =~ 2*GRAMM_WordStruct_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2

# covariances among factors 
eta_1 ~~ eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
GRAMM_WordStruct1 ~~ theta*GRAMM_WordStruct1
GRAMM_WordStruct_2 ~~ theta*GRAMM_WordStruct_2
GRAMM_WordStruct_3 ~~ theta*GRAMM_WordStruct_3

# manifest means (fixed at zero)
GRAMM_WordStruct1 ~ 0*1
GRAMM_WordStruct_2 ~ 0*1
GRAMM_WordStruct_3 ~ 0*1


# modification_1#
GRAMM_WordStruct_2 ~ GRAMM_WordStruct1
GRAMM_WordStruct_3 ~ GRAMM_WordStruct_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~lenghEXPOSURE_ENG_months+0*CURRENT_ENGinputPERC
eta_2~lenghEXPOSURE_ENG_months+0*CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC


#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
'


lg.WS.lavaan_modelINPUT2_fit <- sem(lg.WS.lavaan_modelINPUT2, 
                                    data = DataverySHORT_T3_NEW,
                                    meanstructure = TRUE,
                                    estimator = "ML",
                                    missing = "fiml")

summary(lg.WS.lavaan_modelINPUT2_fit)
summary(lg.WS.lavaan_modelINPUT2_fit, fit.measures=TRUE, standardized = T)
fitMeasures(lg.WS.lavaan_modelINPUT2_fit)


anova(lg.WS.lavaan_modelINPUT1_fit,lg.WS.lavaan_modelINPUT2_fit)

lg.WS.lavaan_modelINPUT3 <- '
# latent variable definitions
#intercept
eta_1 =~ 1*GRAMM_WordStruct1
eta_1 =~ 1*GRAMM_WordStruct_2
eta_1 =~ 1*GRAMM_WordStruct_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*GRAMM_WordStruct1
eta_2 =~ 1*GRAMM_WordStruct_2
eta_2 =~ 2*GRAMM_WordStruct_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2

# covariances among factors 
eta_1 ~~ eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
GRAMM_WordStruct1 ~~ theta*GRAMM_WordStruct1
GRAMM_WordStruct_2 ~~ theta*GRAMM_WordStruct_2
GRAMM_WordStruct_3 ~~ theta*GRAMM_WordStruct_3

# manifest means (fixed at zero)
GRAMM_WordStruct1 ~ 0*1
GRAMM_WordStruct_2 ~ 0*1
GRAMM_WordStruct_3 ~ 0*1


# modification_1#
GRAMM_WordStruct_2 ~ GRAMM_WordStruct1
GRAMM_WordStruct_3 ~ GRAMM_WordStruct_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~0*lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~0*lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC


#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
'


lg.WS.lavaan_modelINPUT3_fit <- sem(lg.WS.lavaan_modelINPUT3, 
                                    data = DataverySHORT_T3_NEW,
                                    meanstructure = TRUE,
                                    estimator = "ML",
                                    missing = "fiml")

summary(lg.WS.lavaan_modelINPUT3_fit)
summary(lg.WS.lavaan_modelINPUT3_fit, fit.measures=TRUE, standardized = T)
fitMeasures(lg.WS.lavaan_modelINPUT3_fit)


anova(lg.WS.lavaan_modelINPUT1_fit,lg.WS.lavaan_modelINPUT3_fit)


lg.WS.FINAL.lavaan_model <- '
# latent variable definitions
#intercept
eta_1 =~ 1*GRAMM_WordStruct1
eta_1 =~ 1*GRAMM_WordStruct_2
eta_1 =~ 1*GRAMM_WordStruct_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*GRAMM_WordStruct1
eta_2 =~ 1*GRAMM_WordStruct_2
eta_2 =~ 2*GRAMM_WordStruct_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2

# covariances among factors 
eta_1 ~~ eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
GRAMM_WordStruct1 ~~ theta*GRAMM_WordStruct1
GRAMM_WordStruct_2 ~~ theta*GRAMM_WordStruct_2
GRAMM_WordStruct_3 ~~ theta*GRAMM_WordStruct_3

# manifest means (fixed at zero)
GRAMM_WordStruct1 ~ 0*1
GRAMM_WordStruct_2 ~ 0*1
GRAMM_WordStruct_3 ~ 0*1


# modification_1#
GRAMM_WordStruct_2 ~ GRAMM_WordStruct1
GRAMM_WordStruct_3 ~ GRAMM_WordStruct_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~lenghEXPOSURE_ENG_months
eta_2~lenghEXPOSURE_ENG_months

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
'


lg.WS_FINAL_model_fit <- sem(lg.WS.FINAL.lavaan_model, 
                       data = DataverySHORT_T3_NEW,
                       meanstructure = TRUE,
                       estimator = "ML",
                       missing = "fiml")

summary(lg.WS_FINAL_model_fit)
summary(lg.WS_FINAL_model_fit, fit.measures=TRUE, standardized = T)  ###this model is great...but look at the parameters! they can't be higher than 1!!!
fitMeasures(lg.WS_FINAL_model_fit)

semPlot::semPaths(lg.WS_FINAL_model_fit, title = FALSE,layout = "tree", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(lg.WS_FINAL_model_fit, title = FALSE,layout = "circle", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(lg.WS_FINAL_model_fit, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)


###bivariate WS and BPVS##########################################################################################################################

##no links##

bivariateBPVS_WS <- '
# latent variable definitions
#intercept BPVS
eta_1BPVS =~ 1*VOC_BPVS1
eta_1BPVS =~ 1*VOC_BPVS_2
eta_1BPVS =~ 1*VOC_BPVS_3

#linear slope BPVS
eta_2BPVS =~ 0*VOC_BPVS1
eta_2BPVS =~ 1*VOC_BPVS_2
eta_2BPVS =~ 2*VOC_BPVS_3

#intercept
eta_1WS =~ 1*GRAMM_WordStruct1
eta_1WS =~ 1*GRAMM_WordStruct_2
eta_1WS =~ 1*GRAMM_WordStruct_3

#linear slope (note intercept is a reserved term)
eta_2WS =~ 0*GRAMM_WordStruct1
eta_2WS =~ 1*GRAMM_WordStruct_2
eta_2WS =~ 2*GRAMM_WordStruct_3

# factor variances
eta_1BPVS ~~ eta_1BPVS
eta_2BPVS ~~ eta_2BPVS
eta_1WS ~~ eta_1WS
eta_2WS ~~ eta_2WS

# covariances among factors 
eta_1BPVS ~~ eta_2BPVS
eta_1WS ~~ eta_2WS
eta_1BPVS ~~ eta_1WS
eta_2BPVS ~~ eta_2WS
eta_1BPVS ~~ eta_2WS
eta_1WS ~~ eta_2BPVS

# factor means 
eta_1BPVS ~ start(35)*1
eta_2BPVS ~ start(4)*1
eta_1WS ~ start(4)*1
eta_2WS ~ start(2)*1

# manifest variances (made equivalent by naming theta)
VOC_BPVS1 ~~ theta*VOC_BPVS1
VOC_BPVS_2 ~~ theta*VOC_BPVS_2
VOC_BPVS_3 ~~ theta*VOC_BPVS_3

# manifest variances (made equivalent by naming theta)
GRAMM_WordStruct1 ~~ theta2*GRAMM_WordStruct1
GRAMM_WordStruct_2 ~~ theta2*GRAMM_WordStruct_2
GRAMM_WordStruct_3 ~~ theta2*GRAMM_WordStruct_3

# manifest means (fixed at zero)
VOC_BPVS1 ~ 0*1
VOC_BPVS_2 ~ 0*1
VOC_BPVS_3 ~ 0*1

# manifest means (fixed at zero)
GRAMM_WordStruct1 ~ 0*1
GRAMM_WordStruct_2 ~ 0*1
GRAMM_WordStruct_3 ~ 0*1


# modification_1BPVS#
VOC_BPVS_2 ~ VOC_BPVS1
VOC_BPVS_3 ~ VOC_BPVS_2

# modification_1WS#
GRAMM_WordStruct_2 ~ GRAMM_WordStruct1
GRAMM_WordStruct_3 ~ GRAMM_WordStruct_2

###covariances###
GRAMM_WordStruct1 ~~ VOC_BPVS1
GRAMM_WordStruct_2 ~~ VOC_BPVS_2
GRAMM_WordStruct_3 ~~ VOC_BPVS_3

#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1BPVS~lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2BPVS~lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_1WS~lenghEXPOSURE_ENG_months
eta_2WS~lenghEXPOSURE_ENG_months

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC

#correlation
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC


#means of covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
'

bivariate_simple <- sem(bivariateBPVS_WS, 
                                  data = DataverySHORT_T3_NEW,
                                  meanstructure = TRUE,
                                  estimator = "ML",
                                  missing = "fiml")

summary(bivariate_simple, fit.measures=TRUE, standardized = T)

semPlot::semPaths(bivariate_simple, title = FALSE,layout = "tree", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(bivariate_simple, title = FALSE,layout = "circle", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(bivariate_simple, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)


##all links##

bivariateCOMPLEX_BPVS_WS <- '
# latent variable definitions
#intercept BPVS
eta_1BPVS =~ 1*VOC_BPVS1
eta_1BPVS =~ 1*VOC_BPVS_2
eta_1BPVS =~ 1*VOC_BPVS_3

#linear slope BPVS
eta_2BPVS =~ 0*VOC_BPVS1
eta_2BPVS =~ 1*VOC_BPVS_2
eta_2BPVS =~ 2*VOC_BPVS_3

#intercept
eta_1WS =~ 1*GRAMM_WordStruct1
eta_1WS =~ 1*GRAMM_WordStruct_2
eta_1WS =~ 1*GRAMM_WordStruct_3

#linear slope (note intercept is a reserved term)
eta_2WS =~ 0*GRAMM_WordStruct1
eta_2WS =~ 1*GRAMM_WordStruct_2
eta_2WS =~ 2*GRAMM_WordStruct_3

# factor variances
eta_1BPVS ~~ eta_1BPVS
eta_2BPVS ~~ eta_2BPVS
eta_1WS ~~ eta_1WS
eta_2WS ~~ eta_2WS

# covariances among factors 
eta_1BPVS ~~ eta_2BPVS
eta_1WS ~~ eta_2WS
eta_1BPVS ~~ eta_1WS
eta_2BPVS ~~ eta_2WS
eta_1BPVS ~~ eta_2WS
eta_1WS ~~ eta_2BPVS

# factor means 
eta_1BPVS ~ start(35)*1
eta_2BPVS ~ start(4)*1
eta_1WS ~ start(4)*1
eta_2WS ~ start(2)*1

# manifest variances (made equivalent by naming theta)
VOC_BPVS1 ~~ theta*VOC_BPVS1
VOC_BPVS_2 ~~ theta*VOC_BPVS_2
VOC_BPVS_3 ~~ theta*VOC_BPVS_3

# manifest variances (made equivalent by naming theta)
GRAMM_WordStruct1 ~~ theta2*GRAMM_WordStruct1
GRAMM_WordStruct_2 ~~ theta2*GRAMM_WordStruct_2
GRAMM_WordStruct_3 ~~ theta2*GRAMM_WordStruct_3

# manifest means (fixed at zero)
VOC_BPVS1 ~ 0*1
VOC_BPVS_2 ~ 0*1
VOC_BPVS_3 ~ 0*1

# manifest means (fixed at zero)
GRAMM_WordStruct1 ~ 0*1
GRAMM_WordStruct_2 ~ 0*1
GRAMM_WordStruct_3 ~ 0*1


# modification_1BPVS#
VOC_BPVS_2 ~ VOC_BPVS1
VOC_BPVS_3 ~ VOC_BPVS_2

# modification_1WS#
GRAMM_WordStruct_2 ~ GRAMM_WordStruct1
GRAMM_WordStruct_3 ~ GRAMM_WordStruct_2

###covariances###
GRAMM_WordStruct1 ~~ VOC_BPVS1
GRAMM_WordStruct_2 ~~ VOC_BPVS_2
GRAMM_WordStruct_3 ~~ VOC_BPVS_3

# BPVS to WS#
GRAMM_WordStruct_2 ~ VOC_BPVS1
GRAMM_WordStruct_3 ~ VOC_BPVS_2

# WS to BPVS#
VOC_BPVS_2 ~ GRAMM_WordStruct1
VOC_BPVS_3 ~ GRAMM_WordStruct_2

#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1BPVS~lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2BPVS~lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_1WS~lenghEXPOSURE_ENG_months
eta_2WS~lenghEXPOSURE_ENG_months

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC

#correlation
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC


#means of covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
'

bivariate_complex <- sem(bivariateCOMPLEX_BPVS_WS, 
                        data = DataverySHORT_T3_NEW,
                        meanstructure = TRUE,
                        estimator = "ML",
                        missing = "fiml")

summary(bivariate_complex, fit.measures=TRUE, standardized = T)

anova(bivariate_complex,bivariate_simple)


##voc to morph##

bivariateoneway_BPVS_WS <- '
# latent variable definitions
#intercept BPVS
eta_1BPVS =~ 1*VOC_BPVS1
eta_1BPVS =~ 1*VOC_BPVS_2
eta_1BPVS =~ 1*VOC_BPVS_3

#linear slope BPVS
eta_2BPVS =~ 0*VOC_BPVS1
eta_2BPVS =~ 1*VOC_BPVS_2
eta_2BPVS =~ 2*VOC_BPVS_3

#intercept
eta_1WS =~ 1*GRAMM_WordStruct1
eta_1WS =~ 1*GRAMM_WordStruct_2
eta_1WS =~ 1*GRAMM_WordStruct_3

#linear slope (note intercept is a reserved term)
eta_2WS =~ 0*GRAMM_WordStruct1
eta_2WS =~ 1*GRAMM_WordStruct_2
eta_2WS =~ 2*GRAMM_WordStruct_3

# factor variances
eta_1BPVS ~~ eta_1BPVS
eta_2BPVS ~~ eta_2BPVS
eta_1WS ~~ eta_1WS
eta_2WS ~~ eta_2WS

# covariances among factors 
eta_1BPVS ~~ eta_2BPVS
eta_1WS ~~ eta_2WS
eta_1BPVS ~~ eta_1WS
eta_2BPVS ~~ eta_2WS
eta_1BPVS ~~ eta_2WS
eta_1WS ~~ eta_2BPVS

# factor means 
eta_1BPVS ~ start(35)*1
eta_2BPVS ~ start(4)*1
eta_1WS ~ start(4)*1
eta_2WS ~ start(2)*1

# manifest variances (made equivalent by naming theta)
VOC_BPVS1 ~~ theta*VOC_BPVS1
VOC_BPVS_2 ~~ theta*VOC_BPVS_2
VOC_BPVS_3 ~~ theta*VOC_BPVS_3

# manifest variances (made equivalent by naming theta)
GRAMM_WordStruct1 ~~ theta2*GRAMM_WordStruct1
GRAMM_WordStruct_2 ~~ theta2*GRAMM_WordStruct_2
GRAMM_WordStruct_3 ~~ theta2*GRAMM_WordStruct_3

# manifest means (fixed at zero)
VOC_BPVS1 ~ 0*1
VOC_BPVS_2 ~ 0*1
VOC_BPVS_3 ~ 0*1

# manifest means (fixed at zero)
GRAMM_WordStruct1 ~ 0*1
GRAMM_WordStruct_2 ~ 0*1
GRAMM_WordStruct_3 ~ 0*1


# modification_1BPVS#
VOC_BPVS_2 ~ VOC_BPVS1
VOC_BPVS_3 ~ VOC_BPVS_2

# modification_1WS#
GRAMM_WordStruct_2 ~ GRAMM_WordStruct1
GRAMM_WordStruct_3 ~ GRAMM_WordStruct_2

###covariances###
GRAMM_WordStruct1 ~~ VOC_BPVS1
GRAMM_WordStruct_2 ~~ VOC_BPVS_2
GRAMM_WordStruct_3 ~~ VOC_BPVS_3

# BPVS to WS#
GRAMM_WordStruct_2 ~ VOC_BPVS1
GRAMM_WordStruct_3 ~ VOC_BPVS_2


#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1BPVS~lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2BPVS~lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_1WS~lenghEXPOSURE_ENG_months
eta_2WS~lenghEXPOSURE_ENG_months

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC

#correlation
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC


#means of covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
'

bivariate_vocTOmorph <- sem(bivariateoneway_BPVS_WS, 
                         data = DataverySHORT_T3_NEW,
                         meanstructure = TRUE,
                         estimator = "ML",
                         missing = "fiml")

summary(bivariate_vocTOmorph, fit.measures=TRUE, standardized = T)

anova(bivariate_vocTOmorph,bivariate_simple)



##morph to voc##

bivariateoneway_WS_BPVS <- '
# latent variable definitions
#intercept BPVS
eta_1BPVS =~ 1*VOC_BPVS1
eta_1BPVS =~ 1*VOC_BPVS_2
eta_1BPVS =~ 1*VOC_BPVS_3

#linear slope BPVS
eta_2BPVS =~ 0*VOC_BPVS1
eta_2BPVS =~ 1*VOC_BPVS_2
eta_2BPVS =~ 2*VOC_BPVS_3

#intercept
eta_1WS =~ 1*GRAMM_WordStruct1
eta_1WS =~ 1*GRAMM_WordStruct_2
eta_1WS =~ 1*GRAMM_WordStruct_3

#linear slope (note intercept is a reserved term)
eta_2WS =~ 0*GRAMM_WordStruct1
eta_2WS =~ 1*GRAMM_WordStruct_2
eta_2WS =~ 2*GRAMM_WordStruct_3

# factor variances
eta_1BPVS ~~ eta_1BPVS
eta_2BPVS ~~ eta_2BPVS
eta_1WS ~~ eta_1WS
eta_2WS ~~ eta_2WS

# covariances among factors 
eta_1BPVS ~~ eta_2BPVS
eta_1WS ~~ eta_2WS
eta_1BPVS ~~ eta_1WS
eta_2BPVS ~~ eta_2WS
eta_1BPVS ~~ eta_2WS
eta_1WS ~~ eta_2BPVS

# factor means 
eta_1BPVS ~ start(35)*1
eta_2BPVS ~ start(4)*1
eta_1WS ~ start(4)*1
eta_2WS ~ start(2)*1

# manifest variances (made equivalent by naming theta)
VOC_BPVS1 ~~ theta*VOC_BPVS1
VOC_BPVS_2 ~~ theta*VOC_BPVS_2
VOC_BPVS_3 ~~ theta*VOC_BPVS_3

# manifest variances (made equivalent by naming theta)
GRAMM_WordStruct1 ~~ theta2*GRAMM_WordStruct1
GRAMM_WordStruct_2 ~~ theta2*GRAMM_WordStruct_2
GRAMM_WordStruct_3 ~~ theta2*GRAMM_WordStruct_3

# manifest means (fixed at zero)
VOC_BPVS1 ~ 0*1
VOC_BPVS_2 ~ 0*1
VOC_BPVS_3 ~ 0*1

# manifest means (fixed at zero)
GRAMM_WordStruct1 ~ 0*1
GRAMM_WordStruct_2 ~ 0*1
GRAMM_WordStruct_3 ~ 0*1


# modification_1BPVS#
VOC_BPVS_2 ~ VOC_BPVS1
VOC_BPVS_3 ~ VOC_BPVS_2

# modification_1WS#
GRAMM_WordStruct_2 ~ GRAMM_WordStruct1
GRAMM_WordStruct_3 ~ GRAMM_WordStruct_2

###covariances###
GRAMM_WordStruct1 ~~ VOC_BPVS1
GRAMM_WordStruct_2 ~~ VOC_BPVS_2
GRAMM_WordStruct_3 ~~ VOC_BPVS_3

# WS to BPVS#
VOC_BPVS_2 ~ GRAMM_WordStruct1
VOC_BPVS_3 ~ GRAMM_WordStruct_2


#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1BPVS~lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2BPVS~lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_1WS~lenghEXPOSURE_ENG_months
eta_2WS~lenghEXPOSURE_ENG_months

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC

#correlation
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC


#means of covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
'

bivariate_morphTOvoc <- sem(bivariateoneway_WS_BPVS, 
                            data = DataverySHORT_T3_NEW,
                            meanstructure = TRUE,
                            estimator = "ML",
                            missing = "fiml")

summary(bivariate_morphTOvoc, fit.measures=TRUE, standardized = T)

semPlot::semPaths(bivariate_morphTOvoc, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)


anova(bivariate_morphTOvoc,bivariate_simple)


###final graph bivariate BPVS_WS###

summary(bivariate_simple, fit.measures=TRUE, standardized = T)
summary(bivariate_simple)


semPlot::semPaths(bivariate_simple, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(bivariate_simple, title = FALSE,layout = "spring", whatLabels = "std", nodeLabels = c("VOC B T1", "VOC B T2", "VOC B T3", "Current Input", "morph T2", "morph T3", "Lenght Exposure", "morph T1", "Intercept VOC", "Slope VOC", "Intercept morph","Slope morph"), intercepts = FALSE, style = "ram", label.cex=1.2, edge.label.cex=1, node.width = 1.5, curvePivot = TRUE)



###bivariate WS and VOC DEPTH##########################################################################################################################


##no links##

bivariateVOC_WS <- '
# latent variable definitions
#intercept BPVS
eta_1voc =~ 1*VOC_mean_SYNOPP1
eta_1voc =~ 1*VOC_meanSYNOPP_2
eta_1voc =~ 1*VOC_meanSYNOPP_3

#linear slope BPVS
eta_2voc =~ 0*VOC_mean_SYNOPP1
eta_2voc =~ 1*VOC_meanSYNOPP_2
eta_2voc =~ 2*VOC_meanSYNOPP_3

#intercept
eta_1WS =~ 1*GRAMM_WordStruct1
eta_1WS =~ 1*GRAMM_WordStruct_2
eta_1WS =~ 1*GRAMM_WordStruct_3

#linear slope (note intercept is a reserved term)
eta_2WS =~ 0*GRAMM_WordStruct1
eta_2WS =~ 1*GRAMM_WordStruct_2
eta_2WS =~ 2*GRAMM_WordStruct_3

# factor variances
eta_1voc ~~ eta_1voc
eta_2voc ~~ eta_2voc
eta_1WS ~~ eta_1WS
eta_2WS ~~ eta_2WS

# covariances among factors 
eta_1voc ~~ eta_2voc
eta_1WS ~~ eta_2WS
eta_1voc ~~ eta_1WS
eta_2voc ~~ eta_2WS
eta_1voc ~~ eta_2WS
eta_1WS ~~ eta_2voc

# factor means 
eta_1voc ~ start(4)*1
eta_2voc ~ start(2)*1
eta_1WS ~ start(4)*1
eta_2WS ~ start(2)*1

# manifest variances (made equivalent by naming theta)
VOC_mean_SYNOPP1 ~~ theta*VOC_mean_SYNOPP1
VOC_meanSYNOPP_2 ~~ theta*VOC_meanSYNOPP_2
VOC_meanSYNOPP_3 ~~ theta*VOC_meanSYNOPP_3

# manifest variances (made equivalent by naming theta)
GRAMM_WordStruct1 ~~ theta2*GRAMM_WordStruct1
GRAMM_WordStruct_2 ~~ theta2*GRAMM_WordStruct_2
GRAMM_WordStruct_3 ~~ theta2*GRAMM_WordStruct_3

# manifest means (fixed at zero)
VOC_mean_SYNOPP1 ~ 0*1
VOC_meanSYNOPP_2 ~ 0*1
VOC_meanSYNOPP_3 ~ 0*1

# manifest means (fixed at zero)
GRAMM_WordStruct1 ~ 0*1
GRAMM_WordStruct_2 ~ 0*1
GRAMM_WordStruct_3 ~ 0*1


# modification_1VOC#
VOC_meanSYNOPP_2 ~ VOC_mean_SYNOPP1
VOC_meanSYNOPP_3 ~ VOC_meanSYNOPP_2

# modification_1WS#
GRAMM_WordStruct_2 ~ GRAMM_WordStruct1
GRAMM_WordStruct_3 ~ GRAMM_WordStruct_2

###covariances###
GRAMM_WordStruct1 ~~ VOC_mean_SYNOPP1
GRAMM_WordStruct_2 ~~ VOC_meanSYNOPP_2
GRAMM_WordStruct_3 ~~ VOC_meanSYNOPP_3

#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1voc~CURRENT_ENGinputPERC+Age_T1
eta_2voc~CURRENT_ENGinputPERC+Age_T1
eta_1WS~lenghEXPOSURE_ENG_months
eta_2WS~lenghEXPOSURE_ENG_months

#correlation covariates
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
lenghEXPOSURE_ENG_months ~~ Age_T1

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Age_T1 ~~ Age_T1

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
Age_T1 ~ 1
'

bivariate_simple_VOC <- sem(bivariateVOC_WS, 
                        data = DataverySHORT_T3_NEW,
                        meanstructure = TRUE,
                        estimator = "ML",
                        missing = "fiml")

summary(bivariate_simple_VOC, fit.measures=TRUE, standardized = T)

semPlot::semPaths(bivariate_simple_VOC, title = FALSE,layout = "tree", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(bivariate_simple_VOC, title = FALSE,layout = "circle", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(bivariate_simple_VOC, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)



##full links##

bivariate_complex_VOC_WS <- '
# latent variable definitions
#intercept BPVS
eta_1voc =~ 1*VOC_mean_SYNOPP1
eta_1voc =~ 1*VOC_meanSYNOPP_2
eta_1voc =~ 1*VOC_meanSYNOPP_3

#linear slope BPVS
eta_2voc =~ 0*VOC_mean_SYNOPP1
eta_2voc =~ 1*VOC_meanSYNOPP_2
eta_2voc =~ 2*VOC_meanSYNOPP_3

#intercept
eta_1WS =~ 1*GRAMM_WordStruct1
eta_1WS =~ 1*GRAMM_WordStruct_2
eta_1WS =~ 1*GRAMM_WordStruct_3

#linear slope (note intercept is a reserved term)
eta_2WS =~ 0*GRAMM_WordStruct1
eta_2WS =~ 1*GRAMM_WordStruct_2
eta_2WS =~ 2*GRAMM_WordStruct_3

# factor variances
eta_1voc ~~ eta_1voc
eta_2voc ~~ eta_2voc
eta_1WS ~~ eta_1WS
eta_2WS ~~ eta_2WS

# covariances among factors 
eta_1voc ~~ eta_2voc
eta_1WS ~~ eta_2WS
eta_1voc ~~ eta_1WS
eta_2voc ~~ eta_2WS
eta_1voc ~~ eta_2WS
eta_1WS ~~ eta_2voc

# factor means 
eta_1voc ~ start(4)*1
eta_2voc ~ start(2)*1
eta_1WS ~ start(4)*1
eta_2WS ~ start(2)*1

# manifest variances (made equivalent by naming theta)
VOC_mean_SYNOPP1 ~~ theta*VOC_mean_SYNOPP1
VOC_meanSYNOPP_2 ~~ theta*VOC_meanSYNOPP_2
VOC_meanSYNOPP_3 ~~ theta*VOC_meanSYNOPP_3

# manifest variances (made equivalent by naming theta)
GRAMM_WordStruct1 ~~ theta2*GRAMM_WordStruct1
GRAMM_WordStruct_2 ~~ theta2*GRAMM_WordStruct_2
GRAMM_WordStruct_3 ~~ theta2*GRAMM_WordStruct_3

# manifest means (fixed at zero)
VOC_mean_SYNOPP1 ~ 0*1
VOC_meanSYNOPP_2 ~ 0*1
VOC_meanSYNOPP_3 ~ 0*1

# manifest means (fixed at zero)
GRAMM_WordStruct1 ~ 0*1
GRAMM_WordStruct_2 ~ 0*1
GRAMM_WordStruct_3 ~ 0*1


# modification_1VOC#
VOC_meanSYNOPP_2 ~ VOC_mean_SYNOPP1
VOC_meanSYNOPP_3 ~ VOC_meanSYNOPP_2

# modification_1WS#
GRAMM_WordStruct_2 ~ GRAMM_WordStruct1
GRAMM_WordStruct_3 ~ GRAMM_WordStruct_2

###covariances###
GRAMM_WordStruct1 ~~ VOC_mean_SYNOPP1
GRAMM_WordStruct_2 ~~ VOC_meanSYNOPP_2
GRAMM_WordStruct_3 ~~ VOC_meanSYNOPP_3

# VOC to WS#
GRAMM_WordStruct_2 ~ VOC_mean_SYNOPP1
GRAMM_WordStruct_3 ~ VOC_meanSYNOPP_2

# WS to VOC#
VOC_meanSYNOPP_2 ~ GRAMM_WordStruct1
VOC_meanSYNOPP_3 ~ GRAMM_WordStruct_2

#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1voc~CURRENT_ENGinputPERC+Age_T1
eta_2voc~CURRENT_ENGinputPERC+Age_T1
eta_1WS~lenghEXPOSURE_ENG_months
eta_2WS~lenghEXPOSURE_ENG_months

#correlation covariates
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
lenghEXPOSURE_ENG_months ~~ Age_T1

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Age_T1 ~~ Age_T1

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
Age_T1 ~ 1
'

bivariate_complex_VOC <- sem(bivariate_complex_VOC_WS, 
                            data = DataverySHORT_T3_NEW,
                            meanstructure = TRUE,
                            estimator = "ML",
                            missing = "fiml")

summary(bivariate_complex_VOC, fit.measures=TRUE, standardized = T)

anova(bivariate_simple_VOC,bivariate_complex_VOC)

##VOC to WS##

bivariate_VOCtoWS <- '
# latent variable definitions
#intercept BPVS
eta_1voc =~ 1*VOC_mean_SYNOPP1
eta_1voc =~ 1*VOC_meanSYNOPP_2
eta_1voc =~ 1*VOC_meanSYNOPP_3

#linear slope BPVS
eta_2voc =~ 0*VOC_mean_SYNOPP1
eta_2voc =~ 1*VOC_meanSYNOPP_2
eta_2voc =~ 2*VOC_meanSYNOPP_3

#intercept
eta_1WS =~ 1*GRAMM_WordStruct1
eta_1WS =~ 1*GRAMM_WordStruct_2
eta_1WS =~ 1*GRAMM_WordStruct_3

#linear slope (note intercept is a reserved term)
eta_2WS =~ 0*GRAMM_WordStruct1
eta_2WS =~ 1*GRAMM_WordStruct_2
eta_2WS =~ 2*GRAMM_WordStruct_3

# factor variances
eta_1voc ~~ eta_1voc
eta_2voc ~~ eta_2voc
eta_1WS ~~ eta_1WS
eta_2WS ~~ eta_2WS

# covariances among factors 
eta_1voc ~~ eta_2voc
eta_1WS ~~ eta_2WS
eta_1voc ~~ eta_1WS
eta_2voc ~~ eta_2WS
eta_1voc ~~ eta_2WS
eta_1WS ~~ eta_2voc

# factor means 
eta_1voc ~ start(4)*1
eta_2voc ~ start(2)*1
eta_1WS ~ start(4)*1
eta_2WS ~ start(2)*1

# manifest variances (made equivalent by naming theta)
VOC_mean_SYNOPP1 ~~ theta*VOC_mean_SYNOPP1
VOC_meanSYNOPP_2 ~~ theta*VOC_meanSYNOPP_2
VOC_meanSYNOPP_3 ~~ theta*VOC_meanSYNOPP_3

# manifest variances (made equivalent by naming theta)
GRAMM_WordStruct1 ~~ theta2*GRAMM_WordStruct1
GRAMM_WordStruct_2 ~~ theta2*GRAMM_WordStruct_2
GRAMM_WordStruct_3 ~~ theta2*GRAMM_WordStruct_3

# manifest means (fixed at zero)
VOC_mean_SYNOPP1 ~ 0*1
VOC_meanSYNOPP_2 ~ 0*1
VOC_meanSYNOPP_3 ~ 0*1

# manifest means (fixed at zero)
GRAMM_WordStruct1 ~ 0*1
GRAMM_WordStruct_2 ~ 0*1
GRAMM_WordStruct_3 ~ 0*1


# modification_1VOC#
VOC_meanSYNOPP_2 ~ VOC_mean_SYNOPP1
VOC_meanSYNOPP_3 ~ VOC_meanSYNOPP_2

# modification_1WS#
GRAMM_WordStruct_2 ~ GRAMM_WordStruct1
GRAMM_WordStruct_3 ~ GRAMM_WordStruct_2

###covariances###
GRAMM_WordStruct1 ~~ VOC_mean_SYNOPP1
GRAMM_WordStruct_2 ~~ VOC_meanSYNOPP_2
GRAMM_WordStruct_3 ~~ VOC_meanSYNOPP_3

# VOC to WS#
GRAMM_WordStruct_2 ~ VOC_mean_SYNOPP1
GRAMM_WordStruct_3 ~ VOC_meanSYNOPP_2

#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1voc~CURRENT_ENGinputPERC+Age_T1
eta_2voc~CURRENT_ENGinputPERC+Age_T1
eta_1WS~lenghEXPOSURE_ENG_months
eta_2WS~lenghEXPOSURE_ENG_months

#correlation covariates
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
lenghEXPOSURE_ENG_months ~~ Age_T1

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Age_T1 ~~ Age_T1

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
Age_T1 ~ 1
'

bivariate_VOCtoWS <- sem(bivariate_VOCtoWS, 
                             data = DataverySHORT_T3_NEW,
                             meanstructure = TRUE,
                             estimator = "ML",
                             missing = "fiml")

summary(bivariate_VOCtoWS, fit.measures=TRUE, standardized = T)

anova(bivariate_VOCtoWS,bivariate_simple_VOC)

##WS to VOC##

bivariate_WStoVOC <- '
# latent variable definitions
#intercept BPVS
eta_1voc =~ 1*VOC_mean_SYNOPP1
eta_1voc =~ 1*VOC_meanSYNOPP_2
eta_1voc =~ 1*VOC_meanSYNOPP_3

#linear slope BPVS
eta_2voc =~ 0*VOC_mean_SYNOPP1
eta_2voc =~ 1*VOC_meanSYNOPP_2
eta_2voc =~ 2*VOC_meanSYNOPP_3

#intercept
eta_1WS =~ 1*GRAMM_WordStruct1
eta_1WS =~ 1*GRAMM_WordStruct_2
eta_1WS =~ 1*GRAMM_WordStruct_3

#linear slope (note intercept is a reserved term)
eta_2WS =~ 0*GRAMM_WordStruct1
eta_2WS =~ 1*GRAMM_WordStruct_2
eta_2WS =~ 2*GRAMM_WordStruct_3

# factor variances
eta_1voc ~~ eta_1voc
eta_2voc ~~ eta_2voc
eta_1WS ~~ eta_1WS
eta_2WS ~~ eta_2WS

# covariances among factors 
eta_1voc ~~ eta_2voc
eta_1WS ~~ eta_2WS
eta_1voc ~~ eta_1WS
eta_2voc ~~ eta_2WS
eta_1voc ~~ eta_2WS
eta_1WS ~~ eta_2voc

# factor means 
eta_1voc ~ start(4)*1
eta_2voc ~ start(2)*1
eta_1WS ~ start(4)*1
eta_2WS ~ start(2)*1

# manifest variances (made equivalent by naming theta)
VOC_mean_SYNOPP1 ~~ theta*VOC_mean_SYNOPP1
VOC_meanSYNOPP_2 ~~ theta*VOC_meanSYNOPP_2
VOC_meanSYNOPP_3 ~~ theta*VOC_meanSYNOPP_3

# manifest variances (made equivalent by naming theta)
GRAMM_WordStruct1 ~~ theta2*GRAMM_WordStruct1
GRAMM_WordStruct_2 ~~ theta2*GRAMM_WordStruct_2
GRAMM_WordStruct_3 ~~ theta2*GRAMM_WordStruct_3

# manifest means (fixed at zero)
VOC_mean_SYNOPP1 ~ 0*1
VOC_meanSYNOPP_2 ~ 0*1
VOC_meanSYNOPP_3 ~ 0*1

# manifest means (fixed at zero)
GRAMM_WordStruct1 ~ 0*1
GRAMM_WordStruct_2 ~ 0*1
GRAMM_WordStruct_3 ~ 0*1


# modification_1VOC#
VOC_meanSYNOPP_2 ~ VOC_mean_SYNOPP1
VOC_meanSYNOPP_3 ~ VOC_meanSYNOPP_2

# modification_1WS#
GRAMM_WordStruct_2 ~ GRAMM_WordStruct1
GRAMM_WordStruct_3 ~ GRAMM_WordStruct_2

###covariances###
GRAMM_WordStruct1 ~~ VOC_mean_SYNOPP1
GRAMM_WordStruct_2 ~~ VOC_meanSYNOPP_2
GRAMM_WordStruct_3 ~~ VOC_meanSYNOPP_3

# WS to VOC#
VOC_meanSYNOPP_2 ~ GRAMM_WordStruct1
VOC_meanSYNOPP_3 ~ GRAMM_WordStruct_2

#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1voc~CURRENT_ENGinputPERC+Age_T1
eta_2voc~CURRENT_ENGinputPERC+Age_T1
eta_1WS~lenghEXPOSURE_ENG_months
eta_2WS~lenghEXPOSURE_ENG_months

#correlation covariates
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
lenghEXPOSURE_ENG_months ~~ Age_T1

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Age_T1 ~~ Age_T1

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
Age_T1 ~ 1
'

bivariate_WS_VOC <- sem(bivariate_WStoVOC, 
                             data = DataverySHORT_T3_NEW,
                             meanstructure = TRUE,
                             estimator = "ML",
                             missing = "fiml")

summary(bivariate_WS_VOC, fit.measures=TRUE, standardized = T)

anova(bivariate_simple_VOC,bivariate_WS_VOC)


### grammar TROG growth##########################################

### no growth#####

summary(DataverySHORT_T3_NEW)

summary(DataSHORTt1t2t3_NEW)

ggplot(DataSHORTt1t2t3_NEW, aes(x = TIME, y = GRAMM_TROGshort, color = as.factor(Subject), group = Subject)) + 
     geom_point() + 
     geom_line() + 
     theme_classic(base_size = 18) + 
    theme(legend.position = "none") + 
     labs(y = "syntax", x = "TIME") +
     scale_x_continuous(breaks=seq(1,3,by=1))

ng.TROG.lavaan_model <- '
# latent variable definitions
#intercept
eta_1 =~ 1*GRAMM_TROGshort1
eta_1 =~ 1*GRAMM_TROGshort_2
eta_1 =~ 1*GRAMM_TROGshort_3

# factor variances
eta_1 ~~ eta_1

# covariances among factors 
#none (only 1 factor)

# factor means 
eta_1 ~ start(5)*1

# manifest variances (made equivalent by naming theta)
GRAMM_TROGshort1 ~~ theta*GRAMM_TROGshort1
GRAMM_TROGshort_2 ~~ theta*GRAMM_TROGshort_2
GRAMM_TROGshort_3 ~~ theta*GRAMM_TROGshort_3

# manifest means (fixed at zero)
GRAMM_TROGshort1 ~ 0*1
GRAMM_TROGshort_2 ~ 0*1
GRAMM_TROGshort_3 ~ 0*1

# modification_1#
GRAMM_TROGshort_2 ~ GRAMM_TROGshort1
GRAMM_TROGshort_3 ~ GRAMM_TROGshort_2

#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1~Mother_EDU+SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Mother_EDU ~~ Mother_EDU
SES ~~ SES

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ Mother_EDU
lenghEXPOSURE_ENG_months ~~ SES
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
CURRENT_ENGinputPERC ~~ Mother_EDU
CURRENT_ENGinputPERC ~~ SES
Mother_EDU ~~ SES

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
Mother_EDU ~ 1
CURRENT_ENGinputPERC ~ 1
SES ~ 1
'
ng.TROG.lavaan_fit <- sem(ng.TROG.lavaan_model, 
                        data = DataverySHORT_T3_NEW,
                        meanstructure = TRUE,
                        estimator = "ML",
                        missing = "fiml")


summary(ng.TROG.lavaan_fit, fit.measures=TRUE, standardized = T)


summary(ng.TROG.lavaan_fit)
fitMeasures(ng.TROG.lavaan_fit)

semPlot::semPaths(ng.TROG.lavaan_fit, title = FALSE,layout = "tree", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(ng.TROG.lavaan_fit, title = FALSE,layout = "circle", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(ng.TROG.lavaan_fit, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)

###linear growth model###

lg.TROG.lavaan_model <- '
# latent variable definitions
#intercept
eta_1 =~ 1*GRAMM_TROGshort1
eta_1 =~ 1*GRAMM_TROGshort_2
eta_1 =~ 1*GRAMM_TROGshort_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*GRAMM_TROGshort1
eta_2 =~ 1*GRAMM_TROGshort_2
eta_2 =~ 2*GRAMM_TROGshort_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2

# covariances among factors 
eta_1 ~~ eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
GRAMM_TROGshort1 ~~ theta*GRAMM_TROGshort1
GRAMM_TROGshort_2 ~~ theta*GRAMM_TROGshort_2
GRAMM_TROGshort_3 ~~ theta*GRAMM_TROGshort_3

# manifest means (fixed at zero)
GRAMM_TROGshort1 ~ 0*1
GRAMM_TROGshort_2 ~ 0*1
GRAMM_TROGshort_3 ~ 0*1


# modification_1#
GRAMM_TROGshort_2 ~ GRAMM_TROGshort1
GRAMM_TROGshort_3 ~ GRAMM_TROGshort_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~Mother_EDU+SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~Mother_EDU+SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Mother_EDU ~~ Mother_EDU
SES ~~ SES

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ Mother_EDU
lenghEXPOSURE_ENG_months ~~ SES
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
CURRENT_ENGinputPERC ~~ Mother_EDU
CURRENT_ENGinputPERC ~~ SES
Mother_EDU ~~ SES

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
Mother_EDU ~ 1
CURRENT_ENGinputPERC ~ 1
SES ~ 1
'


lg.TROG_model_fit <- sem(lg.TROG.lavaan_model, 
                       data = DataverySHORT_T3_NEW,
                       meanstructure = TRUE,
                       estimator = "ML",
                       missing = "fiml")

summary(lg.TROG_model_fit)
summary(lg.TROG_model_fit, fit.measures=TRUE, standardized = T)
fitMeasures(lg.TROG_model_fit)

semPlot::semPaths(lg.TROG_model_fit, title = FALSE,layout = "tree", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(lg.TROG_model_fit, title = FALSE,layout = "circle", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(lg.TROG_model_fit, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)


anova(lg.TROG_model_fit,ng.TROG.lavaan_fit)

lg.TROG2.lavaan_model <- '
# latent variable definitions
#intercept
eta_1 =~ 1*GRAMM_TROGshort1
eta_1 =~ 1*GRAMM_TROGshort_2
eta_1 =~ 1*GRAMM_TROGshort_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*GRAMM_TROGshort1
eta_2 =~ 1*GRAMM_TROGshort_2
eta_2 =~ 2*GRAMM_TROGshort_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ 1*eta_2

# covariances among factors 
eta_1 ~~ 0.1*eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
GRAMM_TROGshort1 ~~ theta*GRAMM_TROGshort1
GRAMM_TROGshort_2 ~~ theta*GRAMM_TROGshort_2
GRAMM_TROGshort_3 ~~ theta*GRAMM_TROGshort_3

# manifest means (fixed at zero)
GRAMM_TROGshort1 ~ 0*1
GRAMM_TROGshort_2 ~ 0*1
GRAMM_TROGshort_3 ~ 0*1


# modification_1#
GRAMM_TROGshort_2 ~ GRAMM_TROGshort1
GRAMM_TROGshort_3 ~ GRAMM_TROGshort_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~Mother_EDU+SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~Mother_EDU+SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Mother_EDU ~~ Mother_EDU
SES ~~ SES

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ Mother_EDU
lenghEXPOSURE_ENG_months ~~ SES
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
CURRENT_ENGinputPERC ~~ Mother_EDU
CURRENT_ENGinputPERC ~~ SES
Mother_EDU ~~ SES

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
Mother_EDU ~ 1
CURRENT_ENGinputPERC ~ 1
SES ~ 1
'


lg.TROG2_model_fit <- sem(lg.TROG2.lavaan_model, 
                         data = DataverySHORT_T3_NEW,
                         meanstructure = TRUE,
                         estimator = "ML",
                         missing = "fiml")

summary(lg.TROG2_model_fit, fit.measures=TRUE, standardized = T)

lg.TROG3.lavaan_model <- '
# latent variable definitions
#intercept
eta_1 =~ 1*GRAMM_TROGshort1
eta_1 =~ 1*GRAMM_TROGshort_2
eta_1 =~ 1*GRAMM_TROGshort_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*GRAMM_TROGshort1
eta_2 =~ 1*GRAMM_TROGshort_2
eta_2 =~ 2*GRAMM_TROGshort_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ 1*eta_2

# covariances among factors 
eta_1 ~~ 0*eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
GRAMM_TROGshort1 ~~ theta*GRAMM_TROGshort1
GRAMM_TROGshort_2 ~~ theta*GRAMM_TROGshort_2
GRAMM_TROGshort_3 ~~ theta*GRAMM_TROGshort_3

# manifest means (fixed at zero)
GRAMM_TROGshort1 ~ 0*1
GRAMM_TROGshort_2 ~ 0*1
GRAMM_TROGshort_3 ~ 0*1


# modification_1#
GRAMM_TROGshort_2 ~ GRAMM_TROGshort1
GRAMM_TROGshort_3 ~ GRAMM_TROGshort_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~Mother_EDU+SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~Mother_EDU+SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Mother_EDU ~~ Mother_EDU
SES ~~ SES

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ Mother_EDU
lenghEXPOSURE_ENG_months ~~ SES
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
CURRENT_ENGinputPERC ~~ Mother_EDU
CURRENT_ENGinputPERC ~~ SES
Mother_EDU ~~ SES

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
Mother_EDU ~ 1
CURRENT_ENGinputPERC ~ 1
SES ~ 1
'


lg.TROG3_model_fit <- sem(lg.TROG3.lavaan_model, 
                          data = DataverySHORT_T3_NEW,
                          meanstructure = TRUE,
                          estimator = "ML",
                          missing = "fiml")

summary(lg.TROG3_model_fit, fit.measures=TRUE, standardized = T) ###correct lg model###

anova(lg.TROG3_model_fit,lg.TROG2_model_fit)

anova(lg.TROG3_model_fit, ng.TROG.lavaan_fit)

anova(lg.TROG3_model_fit,lg.TROG.lavaan_model)

###quad check###

QUAD.TROG3.lavaan_model <- '
# latent variable definitions
#intercept
eta_1 =~ 1*GRAMM_TROGshort1
eta_1 =~ 1*GRAMM_TROGshort_2
eta_1 =~ 1*GRAMM_TROGshort_3

#linear slope (note intercept is a reserved term)
eta_2 =~ -0.5*GRAMM_TROGshort1
eta_2 =~ 0*GRAMM_TROGshort_2
eta_2 =~ 0.5*GRAMM_TROGshort_3

#quadratic slope (note intercept is a reserved term)
eta_3 =~ 0.25*GRAMM_TROGshort1
eta_3 =~ 0*GRAMM_TROGshort_2
eta_3 =~ 0.25*GRAMM_TROGshort_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ 1*eta_2
eta_3 ~~ eta_3

# covariances among factors 
eta_1 ~~ eta_2
eta_1 ~~ eta_3
eta_3 ~~ eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
GRAMM_TROGshort1 ~~ theta*GRAMM_TROGshort1
GRAMM_TROGshort_2 ~~ theta*GRAMM_TROGshort_2
GRAMM_TROGshort_3 ~~ theta*GRAMM_TROGshort_3

# manifest means (fixed at zero)
GRAMM_TROGshort1 ~ 0*1
GRAMM_TROGshort_2 ~ 0*1
GRAMM_TROGshort_3 ~ 0*1


# modification_1#
GRAMM_TROGshort_2 ~ GRAMM_TROGshort1
GRAMM_TROGshort_3 ~ GRAMM_TROGshort_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~Mother_EDU+SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~Mother_EDU+SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Mother_EDU ~~ Mother_EDU
SES ~~ SES

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ Mother_EDU
lenghEXPOSURE_ENG_months ~~ SES
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
CURRENT_ENGinputPERC ~~ Mother_EDU
CURRENT_ENGinputPERC ~~ SES
Mother_EDU ~~ SES

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
Mother_EDU ~ 1
CURRENT_ENGinputPERC ~ 1
SES ~ 1
'


QUAD.TROG3_model_fit <- sem(QUAD.TROG3.lavaan_model, 
                          data = DataverySHORT_T3_NEW,
                          meanstructure = TRUE,
                          estimator = "ML",
                          missing = "fiml")

summary(QUAD.TROG3_model_fit, fit.measures=TRUE, standardized = T)
anova(QUAD.TROG3_model_fit,lg.TROG3_model_fit)

#negative variances, check with final model#

###covariates check###

lg.TROG.COVses.lavaan_model <- '
# latent variable definitions
#intercept
eta_1 =~ 1*GRAMM_TROGshort1
eta_1 =~ 1*GRAMM_TROGshort_2
eta_1 =~ 1*GRAMM_TROGshort_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*GRAMM_TROGshort1
eta_2 =~ 1*GRAMM_TROGshort_2
eta_2 =~ 2*GRAMM_TROGshort_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ 1*eta_2

# covariances among factors 
eta_1 ~~ 0*eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
GRAMM_TROGshort1 ~~ theta*GRAMM_TROGshort1
GRAMM_TROGshort_2 ~~ theta*GRAMM_TROGshort_2
GRAMM_TROGshort_3 ~~ theta*GRAMM_TROGshort_3

# manifest means (fixed at zero)
GRAMM_TROGshort1 ~ 0*1
GRAMM_TROGshort_2 ~ 0*1
GRAMM_TROGshort_3 ~ 0*1


# modification_1#
GRAMM_TROGshort_2 ~ GRAMM_TROGshort1
GRAMM_TROGshort_3 ~ GRAMM_TROGshort_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~Mother_EDU+0*SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~Mother_EDU+0*SES+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Mother_EDU ~~ Mother_EDU
SES ~~ SES

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ Mother_EDU
lenghEXPOSURE_ENG_months ~~ SES
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
CURRENT_ENGinputPERC ~~ Mother_EDU
CURRENT_ENGinputPERC ~~ SES
Mother_EDU ~~ SES

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
Mother_EDU ~ 1
CURRENT_ENGinputPERC ~ 1
SES ~ 1
'


lg.TROGnoSES_model_fit <- sem(lg.TROG.COVses.lavaan_model, 
                          data = DataverySHORT_T3_NEW,
                          meanstructure = TRUE,
                          estimator = "ML",
                          missing = "fiml")

summary(lg.TROGnoSES_model_fit, fit.measures=TRUE, standardized = T)

anova(lg.TROGnoSES_model_fit,lg.TROG3_model_fit)

lg.TROG.COVcorrect.lavaan_model <- '
# latent variable definitions
#intercept
eta_1 =~ 1*GRAMM_TROGshort1
eta_1 =~ 1*GRAMM_TROGshort_2
eta_1 =~ 1*GRAMM_TROGshort_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*GRAMM_TROGshort1
eta_2 =~ 1*GRAMM_TROGshort_2
eta_2 =~ 2*GRAMM_TROGshort_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ 1*eta_2

# covariances among factors 
eta_1 ~~ 0*eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
GRAMM_TROGshort1 ~~ theta*GRAMM_TROGshort1
GRAMM_TROGshort_2 ~~ theta*GRAMM_TROGshort_2
GRAMM_TROGshort_3 ~~ theta*GRAMM_TROGshort_3

# manifest means (fixed at zero)
GRAMM_TROGshort1 ~ 0*1
GRAMM_TROGshort_2 ~ 0*1
GRAMM_TROGshort_3 ~ 0*1


# modification_1#
GRAMM_TROGshort_2 ~ GRAMM_TROGshort1
GRAMM_TROGshort_3 ~ GRAMM_TROGshort_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~Mother_EDU+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~Mother_EDU+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Mother_EDU ~~ Mother_EDU

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ Mother_EDU
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
CURRENT_ENGinputPERC ~~ Mother_EDU

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
Mother_EDU ~ 1
CURRENT_ENGinputPERC ~ 1
'


lg.TROGcorrect_model_fit <- sem(lg.TROG.COVcorrect.lavaan_model, 
                              data = DataverySHORT_T3_NEW,
                              meanstructure = TRUE,
                              estimator = "ML",
                              missing = "fiml")

summary(lg.TROGcorrect_model_fit, fit.measures=TRUE, standardized = T)

lg.TROG.COVmum.lavaan_model <- '
# latent variable definitions
#intercept
eta_1 =~ 1*GRAMM_TROGshort1
eta_1 =~ 1*GRAMM_TROGshort_2
eta_1 =~ 1*GRAMM_TROGshort_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*GRAMM_TROGshort1
eta_2 =~ 1*GRAMM_TROGshort_2
eta_2 =~ 2*GRAMM_TROGshort_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ 1*eta_2

# covariances among factors 
eta_1 ~~ 0*eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
GRAMM_TROGshort1 ~~ theta*GRAMM_TROGshort1
GRAMM_TROGshort_2 ~~ theta*GRAMM_TROGshort_2
GRAMM_TROGshort_3 ~~ theta*GRAMM_TROGshort_3

# manifest means (fixed at zero)
GRAMM_TROGshort1 ~ 0*1
GRAMM_TROGshort_2 ~ 0*1
GRAMM_TROGshort_3 ~ 0*1


# modification_1#
GRAMM_TROGshort_2 ~ GRAMM_TROGshort1
GRAMM_TROGshort_3 ~ GRAMM_TROGshort_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~0*Mother_EDU+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~0*Mother_EDU+lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Mother_EDU ~~ Mother_EDU

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ Mother_EDU
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC
CURRENT_ENGinputPERC ~~ Mother_EDU

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
Mother_EDU ~ 1
CURRENT_ENGinputPERC ~ 1
'


lg.TROGmum_model_fit <- sem(lg.TROG.COVmum.lavaan_model, 
                                data = DataverySHORT_T3_NEW,
                                meanstructure = TRUE,
                                estimator = "ML",
                                missing = "fiml")

summary(lg.TROGmum_model_fit, fit.measures=TRUE, standardized = T)

anova(lg.TROGmum_model_fit,lg.TROGcorrect_model_fit)

lg.TROG.COVcorrectNEW.lavaan_model <- '
# latent variable definitions
#intercept
eta_1 =~ 1*GRAMM_TROGshort1
eta_1 =~ 1*GRAMM_TROGshort_2
eta_1 =~ 1*GRAMM_TROGshort_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*GRAMM_TROGshort1
eta_2 =~ 1*GRAMM_TROGshort_2
eta_2 =~ 2*GRAMM_TROGshort_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ 1*eta_2

# covariances among factors 
eta_1 ~~ 0*eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
GRAMM_TROGshort1 ~~ theta*GRAMM_TROGshort1
GRAMM_TROGshort_2 ~~ theta*GRAMM_TROGshort_2
GRAMM_TROGshort_3 ~~ theta*GRAMM_TROGshort_3

# manifest means (fixed at zero)
GRAMM_TROGshort1 ~ 0*1
GRAMM_TROGshort_2 ~ 0*1
GRAMM_TROGshort_3 ~ 0*1


# modification_1#
GRAMM_TROGshort_2 ~ GRAMM_TROGshort1
GRAMM_TROGshort_3 ~ GRAMM_TROGshort_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
'


lg.TROGcorrectNEW_model_fit <- sem(lg.TROG.COVcorrectNEW.lavaan_model, 
                            data = DataverySHORT_T3_NEW,
                            meanstructure = TRUE,
                            estimator = "ML",
                            missing = "fiml")

summary(lg.TROGcorrectNEW_model_fit, fit.measures=TRUE, standardized = T)

lg.TROG.COVinput1.lavaan_model <- '
# latent variable definitions
#intercept
eta_1 =~ 1*GRAMM_TROGshort1
eta_1 =~ 1*GRAMM_TROGshort_2
eta_1 =~ 1*GRAMM_TROGshort_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*GRAMM_TROGshort1
eta_2 =~ 1*GRAMM_TROGshort_2
eta_2 =~ 2*GRAMM_TROGshort_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ 1*eta_2

# covariances among factors 
eta_1 ~~ 0*eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
GRAMM_TROGshort1 ~~ theta*GRAMM_TROGshort1
GRAMM_TROGshort_2 ~~ theta*GRAMM_TROGshort_2
GRAMM_TROGshort_3 ~~ theta*GRAMM_TROGshort_3

# manifest means (fixed at zero)
GRAMM_TROGshort1 ~ 0*1
GRAMM_TROGshort_2 ~ 0*1
GRAMM_TROGshort_3 ~ 0*1


# modification_1#
GRAMM_TROGshort_2 ~ GRAMM_TROGshort1
GRAMM_TROGshort_3 ~ GRAMM_TROGshort_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~0*lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~0*lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
'


lg.TROGinput1_model_fit <- sem(lg.TROG.COVinput1.lavaan_model, 
                                   data = DataverySHORT_T3_NEW,
                                   meanstructure = TRUE,
                                   estimator = "ML",
                                   missing = "fiml")

summary(lg.TROGinput1_model_fit, fit.measures=TRUE, standardized = T)

anova(lg.TROGinput1_model_fit,lg.TROGcorrectNEW_model_fit)

lg.TROG.COVinput2.lavaan_model <- '
# latent variable definitions
#intercept
eta_1 =~ 1*GRAMM_TROGshort1
eta_1 =~ 1*GRAMM_TROGshort_2
eta_1 =~ 1*GRAMM_TROGshort_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*GRAMM_TROGshort1
eta_2 =~ 1*GRAMM_TROGshort_2
eta_2 =~ 2*GRAMM_TROGshort_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ 1*eta_2

# covariances among factors 
eta_1 ~~ 0*eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
GRAMM_TROGshort1 ~~ theta*GRAMM_TROGshort1
GRAMM_TROGshort_2 ~~ theta*GRAMM_TROGshort_2
GRAMM_TROGshort_3 ~~ theta*GRAMM_TROGshort_3

# manifest means (fixed at zero)
GRAMM_TROGshort1 ~ 0*1
GRAMM_TROGshort_2 ~ 0*1
GRAMM_TROGshort_3 ~ 0*1


# modification_1#
GRAMM_TROGshort_2 ~ GRAMM_TROGshort1
GRAMM_TROGshort_3 ~ GRAMM_TROGshort_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~lenghEXPOSURE_ENG_months+0*CURRENT_ENGinputPERC
eta_2~lenghEXPOSURE_ENG_months+0*CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC

#covariance of TIV covaraites
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
'


lg.TROGinput2_model_fit <- sem(lg.TROG.COVinput2.lavaan_model, 
                               data = DataverySHORT_T3_NEW,
                               meanstructure = TRUE,
                               estimator = "ML",
                               missing = "fiml")

summary(lg.TROGinput2_model_fit, fit.measures=TRUE, standardized = T)

anova(lg.TROGinput2_model_fit,lg.TROGcorrectNEW_model_fit)

###correct model###

summary(lg.TROGcorrectNEW_model_fit, fit.measures=TRUE, standardized = T)

lg.TROG.COVinputFINAL.lavaan_model <- '
# latent variable definitions
#intercept
eta_1 =~ 1*GRAMM_TROGshort1
eta_1 =~ 1*GRAMM_TROGshort_2
eta_1 =~ 1*GRAMM_TROGshort_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*GRAMM_TROGshort1
eta_2 =~ 1*GRAMM_TROGshort_2
eta_2 =~ 2*GRAMM_TROGshort_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ 1*eta_2

# covariances among factors 
eta_1 ~~ 0*eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
GRAMM_TROGshort1 ~~ theta*GRAMM_TROGshort1
GRAMM_TROGshort_2 ~~ theta*GRAMM_TROGshort_2
GRAMM_TROGshort_3 ~~ theta*GRAMM_TROGshort_3

# manifest means (fixed at zero)
GRAMM_TROGshort1 ~ 0*1
GRAMM_TROGshort_2 ~ 0*1
GRAMM_TROGshort_3 ~ 0*1


# modification_1#
GRAMM_TROGshort_2 ~ GRAMM_TROGshort1
GRAMM_TROGshort_3 ~ GRAMM_TROGshort_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~lenghEXPOSURE_ENG_months
eta_2~lenghEXPOSURE_ENG_months

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months

#covariance of TIV covaraites



#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
'


lg.TROG.COVinputFINAL_fit <- sem(lg.TROG.COVinputFINAL.lavaan_model, 
                               data = DataverySHORT_T3_NEW,
                               meanstructure = TRUE,
                               estimator = "ML",
                               missing = "fiml")

summary(lg.TROG.COVinputFINAL_fit, fit.measures=TRUE, standardized = T)

##check for quadratic growth##

quad.TROG.COVinputFINAL.lavaan_model <- '
# latent variable definitions
#intercept
eta_1 =~ 1*GRAMM_TROGshort1
eta_1 =~ 1*GRAMM_TROGshort_2
eta_1 =~ 1*GRAMM_TROGshort_3

#linear slope (note intercept is a reserved term)
eta_2 =~ -0.5*GRAMM_TROGshort1
eta_2 =~ 0*GRAMM_TROGshort_2
eta_2 =~ 0.5*GRAMM_TROGshort_3

#quadratic slope (note intercept is a reserved term)
eta_3 =~ 0.25*GRAMM_TROGshort1
eta_3 =~ 0*GRAMM_TROGshort_2
eta_3 =~ 0.25*GRAMM_TROGshort_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ 1*eta_2
eta_3 ~~ eta_3


# covariances among factors 
eta_1 ~~ 0*eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
GRAMM_TROGshort1 ~~ theta*GRAMM_TROGshort1
GRAMM_TROGshort_2 ~~ theta*GRAMM_TROGshort_2
GRAMM_TROGshort_3 ~~ theta*GRAMM_TROGshort_3

# manifest means (fixed at zero)
GRAMM_TROGshort1 ~ 0*1
GRAMM_TROGshort_2 ~ 0*1
GRAMM_TROGshort_3 ~ 0*1


# modification_1#
GRAMM_TROGshort_2 ~ GRAMM_TROGshort1
GRAMM_TROGshort_3 ~ GRAMM_TROGshort_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~lenghEXPOSURE_ENG_months
eta_2~lenghEXPOSURE_ENG_months

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months

#covariance of TIV covaraites



#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
'


QUAD.TROG.COVinputFINAL_fit <- sem(quad.TROG.COVinputFINAL.lavaan_model, 
                                 data = DataverySHORT_T3_NEW,
                                 meanstructure = TRUE,
                                 estimator = "ML",
                                 missing = "fiml")

summary(QUAD.TROG.COVinputFINAL_fit, fit.measures=TRUE, standardized = T)

###final model for TROG###

summary(lg.TROG.COVinputFINAL_fit, fit.measures=TRUE, standardized = T)


semPlot::semPaths(lg.TROG.COVinputFINAL_fit, title = FALSE,layout = "tree", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(lg.TROG.COVinputFINAL_fit, title = FALSE,layout = "circle", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(lg.TROG.COVinputFINAL_fit, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)


###bivariate TROG and BPVS##########################################################################################################################

##no links##

bivariateBPVS_TROG <- '
# latent variable definitions
#intercept BPVS
eta_1BPVS =~ 1*VOC_BPVS1
eta_1BPVS =~ 1*VOC_BPVS_2
eta_1BPVS =~ 1*VOC_BPVS_3

#linear slope BPVS
eta_2BPVS =~ 0*VOC_BPVS1
eta_2BPVS =~ 1*VOC_BPVS_2
eta_2BPVS =~ 2*VOC_BPVS_3

#intercept
eta_1TROG =~ 1*GRAMM_TROGshort1
eta_1TROG =~ 1*GRAMM_TROGshort_2
eta_1TROG =~ 1*GRAMM_TROGshort_3

#linear slope (note intercept is a reserved term)
eta_2TROG =~ 0*GRAMM_TROGshort1
eta_2TROG =~ 1*GRAMM_TROGshort_2
eta_2TROG =~ 2*GRAMM_TROGshort_3

# factor variances
eta_1BPVS ~~ eta_1BPVS
eta_2BPVS ~~ eta_2BPVS
eta_1TROG ~~ eta_1TROG
eta_2TROG ~~ 1*eta_2TROG

# covariances among factors 
eta_1BPVS ~~ eta_2BPVS
eta_1TROG ~~ 0*eta_2TROG
eta_1BPVS ~~ eta_1TROG
eta_2BPVS ~~ eta_2TROG
eta_1BPVS ~~ eta_2TROG
eta_1TROG ~~ eta_2BPVS

# factor means 
eta_1BPVS ~ start(35)*1
eta_2BPVS ~ start(4)*1
eta_1TROG ~ start(4)*1
eta_2TROG ~ start(2)*1

# manifest variances (made equivalent by naming theta)
VOC_BPVS1 ~~ theta*VOC_BPVS1
VOC_BPVS_2 ~~ theta*VOC_BPVS_2
VOC_BPVS_3 ~~ theta*VOC_BPVS_3

# manifest variances (made equivalent by naming theta)
GRAMM_TROGshort1 ~~ theta2*GRAMM_TROGshort1
GRAMM_TROGshort_2 ~~ theta2*GRAMM_TROGshort_2
GRAMM_TROGshort_3 ~~ theta2*GRAMM_TROGshort_3

# manifest means (fixed at zero)
VOC_BPVS1 ~ 0*1
VOC_BPVS_2 ~ 0*1
VOC_BPVS_3 ~ 0*1

# manifest means (fixed at zero)
GRAMM_TROGshort1 ~ 0*1
GRAMM_TROGshort_2 ~ 0*1
GRAMM_TROGshort_3 ~ 0*1


# modification_1BPVS#
VOC_BPVS_2 ~ VOC_BPVS1
VOC_BPVS_3 ~ VOC_BPVS_2

# modification_1WS#
GRAMM_TROGshort_2 ~ GRAMM_TROGshort1
GRAMM_TROGshort_3 ~ GRAMM_TROGshort_2

###covariances###
GRAMM_TROGshort1 ~~ VOC_BPVS1
GRAMM_TROGshort_2 ~~ VOC_BPVS_2
GRAMM_TROGshort_3 ~~ VOC_BPVS_3

#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1BPVS~lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2BPVS~lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_1TROG~lenghEXPOSURE_ENG_months
eta_2TROG~lenghEXPOSURE_ENG_months

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC

#correlation
CURRENT_ENGinputPERC ~~ lenghEXPOSURE_ENG_months

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
'

bivariate_simple_trog <- sem(bivariateBPVS_TROG, 
                        data = DataverySHORT_T3_NEW,
                        meanstructure = TRUE,
                        estimator = "ML",
                        missing = "fiml")

summary(bivariate_simple_trog, fit.measures=TRUE, standardized = T)

semPlot::semPaths(bivariate_simple_trog, title = FALSE,layout = "tree", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(bivariate_simple_trog, title = FALSE,layout = "circle", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(bivariate_simple_trog, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)


bivariateBPVS_TROG2 <- '
# latent variable definitions
#intercept BPVS
eta_1BPVS =~ 1*VOC_BPVS1
eta_1BPVS =~ 1*VOC_BPVS_2
eta_1BPVS =~ 1*VOC_BPVS_3

#linear slope BPVS
eta_2BPVS =~ 0*VOC_BPVS1
eta_2BPVS =~ 1*VOC_BPVS_2
eta_2BPVS =~ 2*VOC_BPVS_3

#intercept
eta_1TROG =~ 1*GRAMM_TROGshort1
eta_1TROG =~ 1*GRAMM_TROGshort_2
eta_1TROG =~ 1*GRAMM_TROGshort_3

#linear slope (note intercept is a reserved term)
eta_2TROG =~ 0*GRAMM_TROGshort1
eta_2TROG =~ 1*GRAMM_TROGshort_2
eta_2TROG =~ 2*GRAMM_TROGshort_3

# factor variances
eta_1BPVS ~~ eta_1BPVS
eta_2BPVS ~~ eta_2BPVS
eta_1TROG ~~ eta_1TROG
eta_2TROG ~~ 1*eta_2TROG

# covariances among factors 
eta_1BPVS ~~ eta_2BPVS
eta_1TROG ~~ 0*eta_2TROG
eta_1BPVS ~~ eta_1TROG
eta_2BPVS ~~ 0*eta_2TROG
eta_1BPVS ~~ 0*eta_2TROG
eta_1TROG ~~ eta_2BPVS

# factor means 
eta_1BPVS ~ start(35)*1
eta_2BPVS ~ start(4)*1
eta_1TROG ~ start(4)*1
eta_2TROG ~ start(2)*1

# manifest variances (made equivalent by naming theta)
VOC_BPVS1 ~~ theta*VOC_BPVS1
VOC_BPVS_2 ~~ theta*VOC_BPVS_2
VOC_BPVS_3 ~~ theta*VOC_BPVS_3

# manifest variances (made equivalent by naming theta)
GRAMM_TROGshort1 ~~ theta2*GRAMM_TROGshort1
GRAMM_TROGshort_2 ~~ theta2*GRAMM_TROGshort_2
GRAMM_TROGshort_3 ~~ theta2*GRAMM_TROGshort_3

# manifest means (fixed at zero)
VOC_BPVS1 ~ 0*1
VOC_BPVS_2 ~ 0*1
VOC_BPVS_3 ~ 0*1

# manifest means (fixed at zero)
GRAMM_TROGshort1 ~ 0*1
GRAMM_TROGshort_2 ~ 0*1
GRAMM_TROGshort_3 ~ 0*1


# modification_1BPVS#
VOC_BPVS_2 ~ VOC_BPVS1
VOC_BPVS_3 ~ VOC_BPVS_2

# modification_1WS#
GRAMM_TROGshort_2 ~ GRAMM_TROGshort1
GRAMM_TROGshort_3 ~ GRAMM_TROGshort_2

###covariances###
GRAMM_TROGshort1 ~~ VOC_BPVS1
GRAMM_TROGshort_2 ~~ VOC_BPVS_2
GRAMM_TROGshort_3 ~~ VOC_BPVS_3

#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1BPVS~lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2BPVS~lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_1TROG~lenghEXPOSURE_ENG_months
eta_2TROG~lenghEXPOSURE_ENG_months

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC

#correlation
CURRENT_ENGinputPERC ~~ lenghEXPOSURE_ENG_months

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
'

bivariate_simple_trog2 <- sem(bivariateBPVS_TROG2, 
                             data = DataverySHORT_T3_NEW,
                             meanstructure = TRUE,
                             estimator = "ML",
                             missing = "fiml")

summary(bivariate_simple_trog2, fit.measures=TRUE, standardized = T)

semPlot::semPaths(bivariate_simple_trog2, title = FALSE,layout = "tree", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)


##all links##

bivariateCOMPLEX_BPVS_TROG <- '
# latent variable definitions
#intercept BPVS
eta_1BPVS =~ 1*VOC_BPVS1
eta_1BPVS =~ 1*VOC_BPVS_2
eta_1BPVS =~ 1*VOC_BPVS_3

#linear slope BPVS
eta_2BPVS =~ 0*VOC_BPVS1
eta_2BPVS =~ 1*VOC_BPVS_2
eta_2BPVS =~ 2*VOC_BPVS_3

#intercept
eta_1TROG =~ 1*GRAMM_TROGshort1
eta_1TROG =~ 1*GRAMM_TROGshort_2
eta_1TROG =~ 1*GRAMM_TROGshort_3

#linear slope (note intercept is a reserved term)
eta_2TROG =~ 0*GRAMM_TROGshort1
eta_2TROG =~ 1*GRAMM_TROGshort_2
eta_2TROG =~ 2*GRAMM_TROGshort_3

# factor variances
eta_1BPVS ~~ eta_1BPVS
eta_2BPVS ~~ eta_2BPVS
eta_1TROG ~~ eta_1TROG
eta_2TROG ~~ 1*eta_2TROG

# covariances among factors 
eta_1BPVS ~~ eta_2BPVS
eta_1TROG ~~ 0*eta_2TROG
eta_1BPVS ~~ eta_1TROG
eta_2BPVS ~~ 0*eta_2TROG
eta_1BPVS ~~ 0*eta_2TROG
eta_1TROG ~~ eta_2BPVS

# factor means 
eta_1BPVS ~ start(35)*1
eta_2BPVS ~ start(4)*1
eta_1TROG ~ start(4)*1
eta_2TROG ~ start(2)*1

# manifest variances (made equivalent by naming theta)
VOC_BPVS1 ~~ theta*VOC_BPVS1
VOC_BPVS_2 ~~ theta*VOC_BPVS_2
VOC_BPVS_3 ~~ theta*VOC_BPVS_3

# manifest variances (made equivalent by naming theta)
GRAMM_TROGshort1 ~~ theta2*GRAMM_TROGshort1
GRAMM_TROGshort_2 ~~ theta2*GRAMM_TROGshort_2
GRAMM_TROGshort_3 ~~ theta2*GRAMM_TROGshort_3

# manifest means (fixed at zero)
VOC_BPVS1 ~ 0*1
VOC_BPVS_2 ~ 0*1
VOC_BPVS_3 ~ 0*1

# manifest means (fixed at zero)
GRAMM_TROGshort1 ~ 0*1
GRAMM_TROGshort_2 ~ 0*1
GRAMM_TROGshort_3 ~ 0*1


# modification_1BPVS#
VOC_BPVS_2 ~ VOC_BPVS1
VOC_BPVS_3 ~ VOC_BPVS_2

# modification_1WS#
GRAMM_TROGshort_2 ~ GRAMM_TROGshort1
GRAMM_TROGshort_3 ~ GRAMM_TROGshort_2

###covariances###
GRAMM_TROGshort1 ~~ VOC_BPVS1
GRAMM_TROGshort_2 ~~ VOC_BPVS_2
GRAMM_TROGshort_3 ~~ VOC_BPVS_3

# BPVS to trog#
GRAMM_TROGshort_2 ~ VOC_BPVS1
GRAMM_TROGshort_3 ~ VOC_BPVS_2

# trog to BPVS#
VOC_BPVS_2 ~ GRAMM_TROGshort1
VOC_BPVS_3 ~ GRAMM_TROGshort_2

#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1BPVS~lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2BPVS~lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_1TROG~lenghEXPOSURE_ENG_months
eta_2TROG~lenghEXPOSURE_ENG_months

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC

#correlation
CURRENT_ENGinputPERC ~~ lenghEXPOSURE_ENG_months

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
'

bivariate_complex_TROG <- sem(bivariateCOMPLEX_BPVS_TROG, 
                         data = DataverySHORT_T3_NEW,
                         meanstructure = TRUE,
                         estimator = "ML",
                         missing = "fiml")

summary(bivariate_complex_TROG, fit.measures=TRUE, standardized = T)

anova(bivariate_complex_TROG,bivariate_simple_trog2)


##voc to morph##

bivariateCOMPLEX_BPVStoTROG <- '
# latent variable definitions
#intercept BPVS
eta_1BPVS =~ 1*VOC_BPVS1
eta_1BPVS =~ 1*VOC_BPVS_2
eta_1BPVS =~ 1*VOC_BPVS_3

#linear slope BPVS
eta_2BPVS =~ 0*VOC_BPVS1
eta_2BPVS =~ 1*VOC_BPVS_2
eta_2BPVS =~ 2*VOC_BPVS_3

#intercept
eta_1TROG =~ 1*GRAMM_TROGshort1
eta_1TROG =~ 1*GRAMM_TROGshort_2
eta_1TROG =~ 1*GRAMM_TROGshort_3

#linear slope (note intercept is a reserved term)
eta_2TROG =~ 0*GRAMM_TROGshort1
eta_2TROG =~ 1*GRAMM_TROGshort_2
eta_2TROG =~ 2*GRAMM_TROGshort_3

# factor variances
eta_1BPVS ~~ eta_1BPVS
eta_2BPVS ~~ eta_2BPVS
eta_1TROG ~~ eta_1TROG
eta_2TROG ~~ 1*eta_2TROG

# covariances among factors 
eta_1BPVS ~~ eta_2BPVS
eta_1TROG ~~ 0*eta_2TROG
eta_1BPVS ~~ eta_1TROG
eta_2BPVS ~~ 0*eta_2TROG
eta_1BPVS ~~ 0*eta_2TROG
eta_1TROG ~~ eta_2BPVS

# factor means 
eta_1BPVS ~ start(35)*1
eta_2BPVS ~ start(4)*1
eta_1TROG ~ start(4)*1
eta_2TROG ~ start(2)*1

# manifest variances (made equivalent by naming theta)
VOC_BPVS1 ~~ theta*VOC_BPVS1
VOC_BPVS_2 ~~ theta*VOC_BPVS_2
VOC_BPVS_3 ~~ theta*VOC_BPVS_3

# manifest variances (made equivalent by naming theta)
GRAMM_TROGshort1 ~~ theta2*GRAMM_TROGshort1
GRAMM_TROGshort_2 ~~ theta2*GRAMM_TROGshort_2
GRAMM_TROGshort_3 ~~ theta2*GRAMM_TROGshort_3

# manifest means (fixed at zero)
VOC_BPVS1 ~ 0*1
VOC_BPVS_2 ~ 0*1
VOC_BPVS_3 ~ 0*1

# manifest means (fixed at zero)
GRAMM_TROGshort1 ~ 0*1
GRAMM_TROGshort_2 ~ 0*1
GRAMM_TROGshort_3 ~ 0*1


# modification_1BPVS#
VOC_BPVS_2 ~ VOC_BPVS1
VOC_BPVS_3 ~ VOC_BPVS_2

# modification_1WS#
GRAMM_TROGshort_2 ~ GRAMM_TROGshort1
GRAMM_TROGshort_3 ~ GRAMM_TROGshort_2

###covariances###
GRAMM_TROGshort1 ~~ VOC_BPVS1
GRAMM_TROGshort_2 ~~ VOC_BPVS_2
GRAMM_TROGshort_3 ~~ VOC_BPVS_3

# BPVS to trog#
GRAMM_TROGshort_2 ~ VOC_BPVS1
GRAMM_TROGshort_3 ~ VOC_BPVS_2

#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1BPVS~lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2BPVS~lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_1TROG~lenghEXPOSURE_ENG_months
eta_2TROG~lenghEXPOSURE_ENG_months

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC

#correlation
CURRENT_ENGinputPERC ~~ lenghEXPOSURE_ENG_months

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
'

bivariate_complex_BPVStoTROG <- sem(bivariateCOMPLEX_BPVStoTROG, 
                              data = DataverySHORT_T3_NEW,
                              meanstructure = TRUE,
                              estimator = "ML",
                              missing = "fiml")

summary(bivariate_complex_BPVStoTROG, fit.measures=TRUE, standardized = T)

anova(bivariate_complex_BPVStoTROG,bivariate_simple_trog2)

##morph to voc##

bivariateCOMPLEX_TROGtoBPVS <- '
# latent variable definitions
#intercept BPVS
eta_1BPVS =~ 1*VOC_BPVS1
eta_1BPVS =~ 1*VOC_BPVS_2
eta_1BPVS =~ 1*VOC_BPVS_3

#linear slope BPVS
eta_2BPVS =~ 0*VOC_BPVS1
eta_2BPVS =~ 1*VOC_BPVS_2
eta_2BPVS =~ 2*VOC_BPVS_3

#intercept
eta_1TROG =~ 1*GRAMM_TROGshort1
eta_1TROG =~ 1*GRAMM_TROGshort_2
eta_1TROG =~ 1*GRAMM_TROGshort_3

#linear slope (note intercept is a reserved term)
eta_2TROG =~ 0*GRAMM_TROGshort1
eta_2TROG =~ 1*GRAMM_TROGshort_2
eta_2TROG =~ 2*GRAMM_TROGshort_3

# factor variances
eta_1BPVS ~~ eta_1BPVS
eta_2BPVS ~~ eta_2BPVS
eta_1TROG ~~ eta_1TROG
eta_2TROG ~~ 1*eta_2TROG

# covariances among factors 
eta_1BPVS ~~ eta_2BPVS
eta_1TROG ~~ 0*eta_2TROG
eta_1BPVS ~~ eta_1TROG
eta_2BPVS ~~ 0*eta_2TROG
eta_1BPVS ~~ 0*eta_2TROG
eta_1TROG ~~ eta_2BPVS

# factor means 
eta_1BPVS ~ start(35)*1
eta_2BPVS ~ start(4)*1
eta_1TROG ~ start(4)*1
eta_2TROG ~ start(2)*1

# manifest variances (made equivalent by naming theta)
VOC_BPVS1 ~~ theta*VOC_BPVS1
VOC_BPVS_2 ~~ theta*VOC_BPVS_2
VOC_BPVS_3 ~~ theta*VOC_BPVS_3

# manifest variances (made equivalent by naming theta)
GRAMM_TROGshort1 ~~ theta2*GRAMM_TROGshort1
GRAMM_TROGshort_2 ~~ theta2*GRAMM_TROGshort_2
GRAMM_TROGshort_3 ~~ theta2*GRAMM_TROGshort_3

# manifest means (fixed at zero)
VOC_BPVS1 ~ 0*1
VOC_BPVS_2 ~ 0*1
VOC_BPVS_3 ~ 0*1

# manifest means (fixed at zero)
GRAMM_TROGshort1 ~ 0*1
GRAMM_TROGshort_2 ~ 0*1
GRAMM_TROGshort_3 ~ 0*1


# modification_1BPVS#
VOC_BPVS_2 ~ VOC_BPVS1
VOC_BPVS_3 ~ VOC_BPVS_2

# modification_1WS#
GRAMM_TROGshort_2 ~ GRAMM_TROGshort1
GRAMM_TROGshort_3 ~ GRAMM_TROGshort_2

###covariances###
GRAMM_TROGshort1 ~~ VOC_BPVS1
GRAMM_TROGshort_2 ~~ VOC_BPVS_2
GRAMM_TROGshort_3 ~~ VOC_BPVS_3

# trog to BPVS#
VOC_BPVS_2 ~ GRAMM_TROGshort1
VOC_BPVS_3 ~ GRAMM_TROGshort_2

#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1BPVS~lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2BPVS~lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_1TROG~lenghEXPOSURE_ENG_months
eta_2TROG~lenghEXPOSURE_ENG_months

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC

#correlation
CURRENT_ENGinputPERC ~~ lenghEXPOSURE_ENG_months

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
'

bivariate_complex_TROGtoBPVS <- sem(bivariateCOMPLEX_TROGtoBPVS, 
                         data = DataverySHORT_T3_NEW,
                         meanstructure = TRUE,
                         estimator = "ML",
                         missing = "fiml")

summary(bivariate_complex_TROGtoBPVS, fit.measures=TRUE, standardized = T)

anova(bivariate_complex_TROGtoBPVS,bivariate_simple_trog2)


###final graph bivariate BPVS_WS###

summary(bivariate_simple_trog2, fit.measures=TRUE, standardized = T)
summary(bivariate_simple_trog2)


semPlot::semPaths(bivariate_simple_trog2, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(bivariate_simple_trog2, title = FALSE,layout = "spring", whatLabels = "std", nodeLabels = c("VOC B T1", "VOC B T2", "VOC B T3", "English Input", "TROG T2", "TROG T3", "TROG T1", "Intercept VOC", "Slope VOC", "Intercept TROG","Slope TROG"), intercepts = FALSE, style = "ram", label.cex=1.2, edge.label.cex=1, node.width = 1.5, curvePivot = TRUE)





###bivariate TROG and VOC DEPTH##########################################################################################################################


##no links##

bivariateVOCdepth_TROG <- '
# latent variable definitions
#intercept BPVS
eta_1VOC =~ 1*VOC_mean_SYNOPP1
eta_1VOC =~ 1*VOC_meanSYNOPP_2
eta_1VOC =~ 1*VOC_meanSYNOPP_3

#linear slope BPVS
eta_2VOC =~ 0*VOC_mean_SYNOPP1
eta_2VOC =~ 1*VOC_meanSYNOPP_2
eta_2VOC =~ 2*VOC_meanSYNOPP_3

#intercept
eta_1TROG =~ 1*GRAMM_TROGshort1
eta_1TROG =~ 1*GRAMM_TROGshort_2
eta_1TROG =~ 1*GRAMM_TROGshort_3

#linear slope (note intercept is a reserved term)
eta_2TROG =~ 0*GRAMM_TROGshort1
eta_2TROG =~ 1*GRAMM_TROGshort_2
eta_2TROG =~ 2*GRAMM_TROGshort_3

# factor variances
eta_1VOC ~~ eta_1VOC
eta_2VOC ~~ eta_2VOC
eta_1TROG ~~ eta_1TROG
eta_2TROG ~~ 1*eta_2TROG

# covariances among factors 
eta_1VOC ~~ eta_2VOC
eta_1TROG ~~ 0*eta_2TROG
eta_1VOC ~~ eta_1TROG
eta_2VOC ~~ eta_2TROG
eta_1VOC ~~ eta_2TROG
eta_1TROG ~~ eta_2VOC

# factor means 
eta_1VOC ~ start(35)*1
eta_2VOC ~ start(4)*1
eta_1TROG ~ start(4)*1
eta_2TROG ~ start(2)*1

# manifest variances (made equivalent by naming theta)
VOC_mean_SYNOPP1 ~~ theta*VOC_mean_SYNOPP1
VOC_meanSYNOPP_2 ~~ theta*VOC_meanSYNOPP_2
VOC_meanSYNOPP_3 ~~ theta*VOC_meanSYNOPP_3

# manifest variances (made equivalent by naming theta)
GRAMM_TROGshort1 ~~ theta2*GRAMM_TROGshort1
GRAMM_TROGshort_2 ~~ theta2*GRAMM_TROGshort_2
GRAMM_TROGshort_3 ~~ theta2*GRAMM_TROGshort_3

# manifest means (fixed at zero)
GRAMM_TROGshort1 ~ 0*1
GRAMM_TROGshort_2 ~ 0*1
GRAMM_TROGshort_3 ~ 0*1

# manifest means (fixed at zero)
VOC_mean_SYNOPP1 ~ 0*1
VOC_meanSYNOPP_2 ~ 0*1
VOC_meanSYNOPP_3 ~ 0*1


# modification_1BPVS#
VOC_meanSYNOPP_2 ~ VOC_mean_SYNOPP1
VOC_meanSYNOPP_3 ~ VOC_meanSYNOPP_2

# modification_1WS#
GRAMM_TROGshort_2 ~ GRAMM_TROGshort1
GRAMM_TROGshort_3 ~ GRAMM_TROGshort_2

###covariances###
GRAMM_TROGshort1 ~~ VOC_mean_SYNOPP1
GRAMM_TROGshort_2 ~~ VOC_meanSYNOPP_2
GRAMM_TROGshort_3 ~~ VOC_meanSYNOPP_3

#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1VOC~CURRENT_ENGinputPERC+Age_T1
eta_2VOC~CURRENT_ENGinputPERC+Age_T1
eta_1TROG~lenghEXPOSURE_ENG_months
eta_2TROG~lenghEXPOSURE_ENG_months

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Age_T1 ~~ Age_T1

#correlation
CURRENT_ENGinputPERC ~~ lenghEXPOSURE_ENG_months
Age_T1 ~~ lenghEXPOSURE_ENG_months

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
Age_T1 ~ 1
'

bivariate_simple_DEPTHtrog <- sem(bivariateVOCdepth_TROG, 
                             data = DataverySHORT_T3_NEW,
                             meanstructure = TRUE,
                             estimator = "ML",
                             missing = "fiml")

summary(bivariate_simple_DEPTHtrog, fit.measures=TRUE, standardized = T)

semPlot::semPaths(bivariate_simple_DEPTHtrog, title = FALSE,layout = "tree", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(bivariate_simple_DEPTHtrog, title = FALSE,layout = "circle", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(bivariate_simple_DEPTHtrog, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)

##full links##

bivariateVOCdepth_TROG_full <- '
# latent variable definitions
#intercept BPVS
eta_1VOC =~ 1*VOC_mean_SYNOPP1
eta_1VOC =~ 1*VOC_meanSYNOPP_2
eta_1VOC =~ 1*VOC_meanSYNOPP_3

#linear slope BPVS
eta_2VOC =~ 0*VOC_mean_SYNOPP1
eta_2VOC =~ 1*VOC_meanSYNOPP_2
eta_2VOC =~ 2*VOC_meanSYNOPP_3

#intercept
eta_1TROG =~ 1*GRAMM_TROGshort1
eta_1TROG =~ 1*GRAMM_TROGshort_2
eta_1TROG =~ 1*GRAMM_TROGshort_3

#linear slope (note intercept is a reserved term)
eta_2TROG =~ 0*GRAMM_TROGshort1
eta_2TROG =~ 1*GRAMM_TROGshort_2
eta_2TROG =~ 2*GRAMM_TROGshort_3

# factor variances
eta_1VOC ~~ eta_1VOC
eta_2VOC ~~ eta_2VOC
eta_1TROG ~~ eta_1TROG
eta_2TROG ~~ 1*eta_2TROG

# covariances among factors 
eta_1VOC ~~ eta_2VOC
eta_1TROG ~~ 0*eta_2TROG
eta_1VOC ~~ eta_1TROG
eta_2VOC ~~ eta_2TROG
eta_1VOC ~~ eta_2TROG
eta_1TROG ~~ eta_2VOC

# factor means 
eta_1VOC ~ start(35)*1
eta_2VOC ~ start(4)*1
eta_1TROG ~ start(4)*1
eta_2TROG ~ start(2)*1

# manifest variances (made equivalent by naming theta)
VOC_mean_SYNOPP1 ~~ theta*VOC_mean_SYNOPP1
VOC_meanSYNOPP_2 ~~ theta*VOC_meanSYNOPP_2
VOC_meanSYNOPP_3 ~~ theta*VOC_meanSYNOPP_3

# manifest variances (made equivalent by naming theta)
GRAMM_TROGshort1 ~~ theta2*GRAMM_TROGshort1
GRAMM_TROGshort_2 ~~ theta2*GRAMM_TROGshort_2
GRAMM_TROGshort_3 ~~ theta2*GRAMM_TROGshort_3

# manifest means (fixed at zero)
GRAMM_TROGshort1 ~ 0*1
GRAMM_TROGshort_2 ~ 0*1
GRAMM_TROGshort_3 ~ 0*1

# manifest means (fixed at zero)
VOC_mean_SYNOPP1 ~ 0*1
VOC_meanSYNOPP_2 ~ 0*1
VOC_meanSYNOPP_3 ~ 0*1


# modification_1BPVS#
VOC_meanSYNOPP_2 ~ VOC_mean_SYNOPP1
VOC_meanSYNOPP_3 ~ VOC_meanSYNOPP_2

# modification_1WS#
GRAMM_TROGshort_2 ~ GRAMM_TROGshort1
GRAMM_TROGshort_3 ~ GRAMM_TROGshort_2

###covariances###
GRAMM_TROGshort1 ~~ VOC_mean_SYNOPP1
GRAMM_TROGshort_2 ~~ VOC_meanSYNOPP_2
GRAMM_TROGshort_3 ~~ VOC_meanSYNOPP_3

# voc to trog#
GRAMM_TROGshort_2 ~ VOC_mean_SYNOPP1
GRAMM_TROGshort_3 ~ VOC_meanSYNOPP_2

# trog to voc#
VOC_meanSYNOPP_2 ~ GRAMM_TROGshort1
VOC_meanSYNOPP_3 ~ GRAMM_TROGshort_2

#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1VOC~CURRENT_ENGinputPERC+Age_T1
eta_2VOC~CURRENT_ENGinputPERC+Age_T1
eta_1TROG~lenghEXPOSURE_ENG_months
eta_2TROG~lenghEXPOSURE_ENG_months

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Age_T1 ~~ Age_T1

#correlation
CURRENT_ENGinputPERC ~~ lenghEXPOSURE_ENG_months
Age_T1 ~~ lenghEXPOSURE_ENG_months

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
Age_T1 ~ 1
'

bivariate_complex_DEPTHtrog <- sem(bivariateVOCdepth_TROG_full, 
                                  data = DataverySHORT_T3_NEW,
                                  meanstructure = TRUE,
                                  estimator = "ML",
                                  missing = "fiml")

summary(bivariate_complex_DEPTHtrog, fit.measures=TRUE, standardized = T)

anova(bivariate_simple_DEPTHtrog,bivariate_complex_DEPTHtrog)



##VOC to TROG##

bivariateVOCdepthtoTROG <- '
# latent variable definitions
#intercept BPVS
eta_1VOC =~ 1*VOC_mean_SYNOPP1
eta_1VOC =~ 1*VOC_meanSYNOPP_2
eta_1VOC =~ 1*VOC_meanSYNOPP_3

#linear slope BPVS
eta_2VOC =~ 0*VOC_mean_SYNOPP1
eta_2VOC =~ 1*VOC_meanSYNOPP_2
eta_2VOC =~ 2*VOC_meanSYNOPP_3

#intercept
eta_1TROG =~ 1*GRAMM_TROGshort1
eta_1TROG =~ 1*GRAMM_TROGshort_2
eta_1TROG =~ 1*GRAMM_TROGshort_3

#linear slope (note intercept is a reserved term)
eta_2TROG =~ 0*GRAMM_TROGshort1
eta_2TROG =~ 1*GRAMM_TROGshort_2
eta_2TROG =~ 2*GRAMM_TROGshort_3

# factor variances
eta_1VOC ~~ eta_1VOC
eta_2VOC ~~ eta_2VOC
eta_1TROG ~~ eta_1TROG
eta_2TROG ~~ 1*eta_2TROG

# covariances among factors 
eta_1VOC ~~ eta_2VOC
eta_1TROG ~~ 0*eta_2TROG
eta_1VOC ~~ eta_1TROG
eta_2VOC ~~ eta_2TROG
eta_1VOC ~~ eta_2TROG
eta_1TROG ~~ eta_2VOC

# factor means 
eta_1VOC ~ start(35)*1
eta_2VOC ~ start(4)*1
eta_1TROG ~ start(4)*1
eta_2TROG ~ start(2)*1

# manifest variances (made equivalent by naming theta)
VOC_mean_SYNOPP1 ~~ theta*VOC_mean_SYNOPP1
VOC_meanSYNOPP_2 ~~ theta*VOC_meanSYNOPP_2
VOC_meanSYNOPP_3 ~~ theta*VOC_meanSYNOPP_3

# manifest variances (made equivalent by naming theta)
GRAMM_TROGshort1 ~~ theta2*GRAMM_TROGshort1
GRAMM_TROGshort_2 ~~ theta2*GRAMM_TROGshort_2
GRAMM_TROGshort_3 ~~ theta2*GRAMM_TROGshort_3

# manifest means (fixed at zero)
GRAMM_TROGshort1 ~ 0*1
GRAMM_TROGshort_2 ~ 0*1
GRAMM_TROGshort_3 ~ 0*1

# manifest means (fixed at zero)
VOC_mean_SYNOPP1 ~ 0*1
VOC_meanSYNOPP_2 ~ 0*1
VOC_meanSYNOPP_3 ~ 0*1


# modification_1BPVS#
VOC_meanSYNOPP_2 ~ VOC_mean_SYNOPP1
VOC_meanSYNOPP_3 ~ VOC_meanSYNOPP_2

# modification_1WS#
GRAMM_TROGshort_2 ~ GRAMM_TROGshort1
GRAMM_TROGshort_3 ~ GRAMM_TROGshort_2

###covariances###
GRAMM_TROGshort1 ~~ VOC_mean_SYNOPP1
GRAMM_TROGshort_2 ~~ VOC_meanSYNOPP_2
GRAMM_TROGshort_3 ~~ VOC_meanSYNOPP_3

# voc to trog#
GRAMM_TROGshort_2 ~ VOC_mean_SYNOPP1
GRAMM_TROGshort_3 ~ VOC_meanSYNOPP_2

#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1VOC~CURRENT_ENGinputPERC+Age_T1
eta_2VOC~CURRENT_ENGinputPERC+Age_T1
eta_1TROG~lenghEXPOSURE_ENG_months
eta_2TROG~lenghEXPOSURE_ENG_months

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Age_T1 ~~ Age_T1

#correlation
CURRENT_ENGinputPERC ~~ lenghEXPOSURE_ENG_months
Age_T1 ~~ lenghEXPOSURE_ENG_months

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
Age_T1 ~ 1
'

bivariate_complex_DEPTHtoTROG <- sem(bivariateVOCdepthtoTROG, 
                                   data = DataverySHORT_T3_NEW,
                                   meanstructure = TRUE,
                                   estimator = "ML",
                                   missing = "fiml")

summary(bivariate_complex_DEPTHtoTROG, fit.measures=TRUE, standardized = T)

anova(bivariate_simple_DEPTHtrog,bivariate_complex_DEPTHtoTROG)


##TROG to VOC##

bivariate_TROGtoVOCdepth <- '
# latent variable definitions
#intercept BPVS
eta_1VOC =~ 1*VOC_mean_SYNOPP1
eta_1VOC =~ 1*VOC_meanSYNOPP_2
eta_1VOC =~ 1*VOC_meanSYNOPP_3

#linear slope BPVS
eta_2VOC =~ 0*VOC_mean_SYNOPP1
eta_2VOC =~ 1*VOC_meanSYNOPP_2
eta_2VOC =~ 2*VOC_meanSYNOPP_3

#intercept
eta_1TROG =~ 1*GRAMM_TROGshort1
eta_1TROG =~ 1*GRAMM_TROGshort_2
eta_1TROG =~ 1*GRAMM_TROGshort_3

#linear slope (note intercept is a reserved term)
eta_2TROG =~ 0*GRAMM_TROGshort1
eta_2TROG =~ 1*GRAMM_TROGshort_2
eta_2TROG =~ 2*GRAMM_TROGshort_3

# factor variances
eta_1VOC ~~ eta_1VOC
eta_2VOC ~~ eta_2VOC
eta_1TROG ~~ eta_1TROG
eta_2TROG ~~ 1*eta_2TROG

# covariances among factors 
eta_1VOC ~~ eta_2VOC
eta_1TROG ~~ 0*eta_2TROG
eta_1VOC ~~ eta_1TROG
eta_2VOC ~~ eta_2TROG
eta_1VOC ~~ eta_2TROG
eta_1TROG ~~ eta_2VOC

# factor means 
eta_1VOC ~ start(35)*1
eta_2VOC ~ start(4)*1
eta_1TROG ~ start(4)*1
eta_2TROG ~ start(2)*1

# manifest variances (made equivalent by naming theta)
VOC_mean_SYNOPP1 ~~ theta*VOC_mean_SYNOPP1
VOC_meanSYNOPP_2 ~~ theta*VOC_meanSYNOPP_2
VOC_meanSYNOPP_3 ~~ theta*VOC_meanSYNOPP_3

# manifest variances (made equivalent by naming theta)
GRAMM_TROGshort1 ~~ theta2*GRAMM_TROGshort1
GRAMM_TROGshort_2 ~~ theta2*GRAMM_TROGshort_2
GRAMM_TROGshort_3 ~~ theta2*GRAMM_TROGshort_3

# manifest means (fixed at zero)
GRAMM_TROGshort1 ~ 0*1
GRAMM_TROGshort_2 ~ 0*1
GRAMM_TROGshort_3 ~ 0*1

# manifest means (fixed at zero)
VOC_mean_SYNOPP1 ~ 0*1
VOC_meanSYNOPP_2 ~ 0*1
VOC_meanSYNOPP_3 ~ 0*1


# modification_1BPVS#
VOC_meanSYNOPP_2 ~ VOC_mean_SYNOPP1
VOC_meanSYNOPP_3 ~ VOC_meanSYNOPP_2

# modification_1WS#
GRAMM_TROGshort_2 ~ GRAMM_TROGshort1
GRAMM_TROGshort_3 ~ GRAMM_TROGshort_2

###covariances###
GRAMM_TROGshort1 ~~ VOC_mean_SYNOPP1
GRAMM_TROGshort_2 ~~ VOC_meanSYNOPP_2
GRAMM_TROGshort_3 ~~ VOC_meanSYNOPP_3

# trog to voc#
VOC_meanSYNOPP_2 ~ GRAMM_TROGshort1
VOC_meanSYNOPP_3 ~ GRAMM_TROGshort_2

#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1VOC~CURRENT_ENGinputPERC+Age_T1
eta_2VOC~CURRENT_ENGinputPERC+Age_T1
eta_1TROG~lenghEXPOSURE_ENG_months
eta_2TROG~lenghEXPOSURE_ENG_months

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Age_T1 ~~ Age_T1

#correlation
CURRENT_ENGinputPERC ~~ lenghEXPOSURE_ENG_months
Age_T1 ~~ lenghEXPOSURE_ENG_months

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
Age_T1 ~ 1
'

bivariate_TROGtoDEPTH <- sem(bivariate_TROGtoVOCdepth, 
                                   data = DataverySHORT_T3_NEW,
                                   meanstructure = TRUE,
                                   estimator = "ML",
                                   missing = "fiml")

summary(bivariate_TROGtoDEPTH, fit.measures=TRUE, standardized = T)

anova(bivariate_simple_DEPTHtrog,bivariate_TROGtoDEPTH)


###final model###

summary(bivariate_simple_DEPTHtrog, fit.measures=TRUE, standardized = T)

semPlot::semPaths(bivariate_simple_DEPTHtrog, title = FALSE,layout = "tree", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(bivariate_simple_DEPTHtrog, title = FALSE,layout = "circle", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(bivariate_simple_DEPTHtrog, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)


###Control for modifications in simple growth models###

### final COV model BPVS####

lgCOV.BPVS.FINAL_model_NOmod <- '
# latent variable definitions
#intercept (note intercept is a reserved term)
eta_1 =~ 1*VOC_BPVS1
eta_1 =~ 1*VOC_BPVS_2
eta_1 =~ 1*VOC_BPVS_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*VOC_BPVS1
eta_2 =~ 1*VOC_BPVS_2
eta_2 =~ 2*VOC_BPVS_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2

# covariances among factors 
eta_1 ~~ eta_2

# factor means 
eta_1 ~ start(35)*1
eta_2 ~ start(4)*1

# manifest variances (made equivalent by naming theta)
VOC_BPVS1 ~~ theta*VOC_BPVS1
VOC_BPVS_2 ~~ theta*VOC_BPVS_2
VOC_BPVS_3 ~~ theta*VOC_BPVS_3

# manifest means (fixed at zero)
VOC_BPVS1 ~ 0*1
VOC_BPVS_2 ~ 0*1
VOC_BPVS_3 ~ 0*1


# modification_1#
#VOC_BPVS_2 ~ VOC_BPVS1
#VOC_BPVS_3 ~ VOC_BPVS_2

#Time invariant covariates
#regression of time-invariant covariate on intercept and slope factors
eta_1~lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC
eta_2~lenghEXPOSURE_ENG_months+CURRENT_ENGinputPERC

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC

#covariance of TIV covariates
lenghEXPOSURE_ENG_months ~~ CURRENT_ENGinputPERC

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
CURRENT_ENGinputPERC ~ 1
'

lg.BPVS.COVFINAL_fitNOmod <- sem(lgCOV.BPVS.FINAL_model_NOmod, 
                            data = DataverySHORT_T3_NEW,
                            meanstructure = TRUE,
                            estimator = "ML",
                            missing = "fiml")

summary(lg.BPVS.COVFINAL_fitNOmod, fit.measures=TRUE, standardized = T)
fitMeasures(lg.BPVS.COVFINAL_fit)
anova(lg.BPVS.COVFINAL_fitNOmod,lg.BPVS.COVFINAL_fit)

##use fixed.x =FALSE and fixed.x =TRUE to check

semPlot::semPaths(lg.BPVS.COVFINAL_fit, title = FALSE,layout = "tree", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(lg.BPVS.COVFINAL_fit, title = FALSE,layout = "circle", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(lg.BPVS.COVFINAL_fitNOmod, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)

lgFINAL.DEPTH.lavaan_model_NEW_NOmod <- '
# latent variable definitions
#intercept
eta_1 =~ 1*VOC_mean_SYNOPP1
eta_1 =~ 1*VOC_meanSYNOPP_2
eta_1 =~ 1*VOC_meanSYNOPP_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*VOC_mean_SYNOPP1
eta_2 =~ 1*VOC_meanSYNOPP_2
eta_2 =~ 2*VOC_meanSYNOPP_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2

# covariances among factors 
eta_1 ~~ eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
VOC_mean_SYNOPP1 ~~ theta*VOC_mean_SYNOPP1
VOC_meanSYNOPP_2 ~~ theta*VOC_meanSYNOPP_2
VOC_meanSYNOPP_3 ~~ theta*VOC_meanSYNOPP_3

# manifest means (fixed at zero)
VOC_mean_SYNOPP1 ~ 0*1
VOC_meanSYNOPP_2 ~ 0*1
VOC_meanSYNOPP_3 ~ 0*1


# modification_1#
#VOC_meanSYNOPP_2 ~ VOC_mean_SYNOPP1
#VOC_meanSYNOPP_3 ~ VOC_meanSYNOPP_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~Age_T1+CURRENT_ENGinputPERC
eta_2~Age_T1+CURRENT_ENGinputPERC

#variance of TIV covariates
CURRENT_ENGinputPERC ~~ CURRENT_ENGinputPERC
Age_T1 ~~ Age_T1

#means of TIV covariates (freely estimated)
CURRENT_ENGinputPERC ~ 1
Age_T1 ~ 1

'


lg.DEPTH_model_fitFINAL_NEW_NOmod <- sem(lgFINAL.DEPTH.lavaan_model_NEW_NOmod, 
                                   data = DataverySHORT_T3_NEW,
                                   meanstructure = TRUE,
                                   estimator = "ML",
                                   missing = "fiml")

summary(lg.DEPTH_model_fitFINAL_NEW_NOmod)
summary(lg.DEPTH_model_fitFINAL_NEW_NOmod, fit.measures=TRUE, standardized = T)
anova(lg.DEPTH_model_fitFINAL_NEW_NOmod,lg.DEPTH_model_fitFINAL_NEW)

summary(lg.DEPTH_model_fitFINAL_NEW, fit.measures=TRUE, standardized = T)


semPlot::semPaths(lg.DEPTH_model_fitFINAL_NEW, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 2, curvePivot = TRUE)
semPlot::semPaths(lg.DEPTH_model_fitFINAL_NEW, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram",nodeLabels = c("VOC DEPTH T1", "VOC DEPTH T2", "VOC DEPTH T3", "Age T1", "Current Input","Intercept","Slope"), label.cex=1.0, edge.label.cex=1, node.width = 2, curvePivot = TRUE)


lg.WS.FINAL.lavaan_model_NOmod <- '
# latent variable definitions
#intercept
eta_1 =~ 1*GRAMM_WordStruct1
eta_1 =~ 1*GRAMM_WordStruct_2
eta_1 =~ 1*GRAMM_WordStruct_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*GRAMM_WordStruct1
eta_2 =~ 1*GRAMM_WordStruct_2
eta_2 =~ 2*GRAMM_WordStruct_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ eta_2

# covariances among factors 
eta_1 ~~ eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
GRAMM_WordStruct1 ~~ theta*GRAMM_WordStruct1
GRAMM_WordStruct_2 ~~ theta*GRAMM_WordStruct_2
GRAMM_WordStruct_3 ~~ theta*GRAMM_WordStruct_3

# manifest means (fixed at zero)
GRAMM_WordStruct1 ~ 0*1
GRAMM_WordStruct_2 ~ 0*1
GRAMM_WordStruct_3 ~ 0*1


# modification_1#
#GRAMM_WordStruct_2 ~ GRAMM_WordStruct1
#GRAMM_WordStruct_3 ~ GRAMM_WordStruct_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~lenghEXPOSURE_ENG_months
eta_2~lenghEXPOSURE_ENG_months

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months

#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
'


lg.WS_FINAL_model_fit_NOmod <- sem(lg.WS.FINAL.lavaan_model_NOmod, 
                             data = DataverySHORT_T3_NEW,
                             meanstructure = TRUE,
                             estimator = "ML",
                             missing = "fiml")

summary(lg.WS_FINAL_model_fit_NOmod)
summary(lg.WS_FINAL_model_fit_NOmod, fit.measures=TRUE, standardized = T)  ###this model is great...but look at the parameters! they can't be higher than 1!!!
fitMeasures(lg.WS_FINAL_model_fit_NOmod)

semPlot::semPaths(lg.WS_FINAL_model_fit_NOmod, title = FALSE,layout = "tree", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(lg.WS_FINAL_model_fit_NOmod, title = FALSE,layout = "circle", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(lg.WS_FINAL_model_fit_NOmod, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)

anova(lg.WS_FINAL_model_fit_NOmod,lg.WS_FINAL_model_fit)

lg.TROG.COVinputFINAL.lavaan_model_NOmod <- '
# latent variable definitions
#intercept
eta_1 =~ 1*GRAMM_TROGshort1
eta_1 =~ 1*GRAMM_TROGshort_2
eta_1 =~ 1*GRAMM_TROGshort_3

#linear slope (note intercept is a reserved term)
eta_2 =~ 0*GRAMM_TROGshort1
eta_2 =~ 1*GRAMM_TROGshort_2
eta_2 =~ 2*GRAMM_TROGshort_3

# factor variances
eta_1 ~~ eta_1
eta_2 ~~ 1*eta_2

# covariances among factors 
eta_1 ~~ 0*eta_2

# factor means 
eta_1 ~ start(4)*1
eta_2 ~ start(2)*1

# manifest variances (made equivalent by naming theta)
GRAMM_TROGshort1 ~~ theta*GRAMM_TROGshort1
GRAMM_TROGshort_2 ~~ theta*GRAMM_TROGshort_2
GRAMM_TROGshort_3 ~~ theta*GRAMM_TROGshort_3

# manifest means (fixed at zero)
GRAMM_TROGshort1 ~ 0*1
GRAMM_TROGshort_2 ~ 0*1
GRAMM_TROGshort_3 ~ 0*1


# modification_1#
#GRAMM_TROGshort_2 ~ GRAMM_TROGshort1
#GRAMM_TROGshort_3 ~ GRAMM_TROGshort_2

#regression of time-invariant covariate on intercept and slope factors
eta_1~lenghEXPOSURE_ENG_months
eta_2~lenghEXPOSURE_ENG_months

#variance of TIV covariates
lenghEXPOSURE_ENG_months ~~ lenghEXPOSURE_ENG_months

#covariance of TIV covaraites



#means of TIV covariates (freely estimated)
lenghEXPOSURE_ENG_months ~ 1
'


lg.TROG.COVinputFINAL_fit_NOmod <- sem(lg.TROG.COVinputFINAL.lavaan_model_NOmod, 
                                 data = DataverySHORT_T3_NEW,
                                 meanstructure = TRUE,
                                 estimator = "ML",
                                 missing = "fiml")

summary(lg.TROG.COVinputFINAL_fit_NOmod, fit.measures=TRUE, standardized = T)
anova(lg.TROG.COVinputFINAL_fit_NOmod,lg.TROG.COVinputFINAL_fit)

semPlot::semPaths(lg.TROG.COVinputFINAL_fit_NOmod, title = FALSE,layout = "tree", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(lg.TROG.COVinputFINAL_fit_NOmod, title = FALSE,layout = "circle", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)
semPlot::semPaths(lg.TROG.COVinputFINAL_fit_NOmod, title = FALSE,layout = "spring", whatLabels = "std", intercepts = FALSE, style = "ram", label.cex=1.0, edge.label.cex=1, node.width = 1, curvePivot = TRUE)

