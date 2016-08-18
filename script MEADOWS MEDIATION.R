# Keith Lohse, PhD
# School of Kinesiology, Auburn University
# 2016 - 07 - 24

# Analysis script for Caroline Meadow's Reward Mediation Paper

## -----------------------------------------------------------------------
# First, make sure you have the following packages installed and then
# open them to your library:
library("ggplot2");library("lme4");library("dplyr");library("lmerTest");

# Make sure that the "data REWARD MASTER.txt" file is saved in your 
# working directory. In my case, I have a "data" folder inside of my 
# working directory. You will need to update your file path accordingly:
DATA<-read.table("./data/data REWARD MASTER.txt", header = TRUE, sep="\t") 
str(DATA)
# DATA should contain 3360 observations for 11 variables
head(DATA)
# subID = subject identification number
# Trial = trial number, 168 trials for each subject
# Reward = the incentive value, in USD, for each trial
# Reward.nom = categorical version of Reward variable
# Beta = natural log transform of Beta power averaged across electrodes. 
# Time = time of grip force onset in physio data
# RT = premotor reaction time, defined as time from go signal to EMG peak
# Max = RMSE recified value for EMG peak
# Mean = average RMSE-EMG value over the premotor window
# Stddev = standard deviation of RMSE-EMG value over the premotor window
# Plaus = Plausibility rating for each subject as a percentage (1.0 = believed incentives)

## The following variables all measure beta power at each electrode independently
# FC5 = natural log transform of Beta power at electrode FC5.
# FC3 = natural log transform of Beta power at electrode FC3.
# FC1 = natural log transform of Beta power at electrode FC1.
# C5 = natural log transform of Beta power at electrode C5.
# C3 = natural log transform of Beta power at electrode C3.
# C1 = natural log transform of Beta power at electrode C1.
# CP5 = natural log transform of Beta power at electrode CP5.
# CP3 = natural log transform of Beta power at electrode CP3.
# CP1 = natural log transform of Beta power at electrode CP1.


# Next, we want to remove missing cases for the natural log of beta power 
# (Beta) and reaction time (RT). We will store the result in a new data 
# frame called "FILTER":
FILTER<-subset(DATA,!Beta == 'NA')
FILTER<-subset(FILTER,!RT == 'NA')
head(FILTER)
# Values for monetary Rewards should range from 0 to 4.96$
summary(FILTER$Reward)

# Next, we can also condense our data to a single row for each subject.
# We will average across the dynamic variables to get a single observation 
# for each person.
IID<-summarize(group_by(FILTER, subID),
                      ID = subID[1], 
                      mReward = mean(Reward),
                      mBeta = mean(Beta),
                      mRT = mean(RT),
                      mPlaus = mean(Plaus))
head(IID)
# We can export this independent individual data to a new .csv file
write.csv(IID, file="Meadows_IID.csv")

# We can also use the IID file to look for between subject relationships
plot(IID$mPlaus, IID$mRT, pch = 19, col = 2, bty='n', cex.lab = 1)
cor.test(IID$mPlaus,IID$mRT)

plot(IID$mPlaus, IID$mBeta, pch = 19, col = 3, bty='n', cex.lab = 1)
cor.test(IID$mBeta,IID$mRT)


## -----------------------------------------------------------------------
# Descriptive Stats and Centering Variables
# Descriptive stats for the Beta variable
summary(FILTER$Beta)
FILTER$Beta.c<-FILTER$Beta+2.351

# Descriptive stats for the Reward variable
summary(FILTER$Reward)
FILTER$Reward.c<-FILTER$Reward-1.798

# Descriptive stats for the RT variable
summary(FILTER$RT)
plot(density(FILTER$RT))
FILTER$lnRT<-log(FILTER$RT)
plot(density(FILTER$lnRT))

## -----------------------------------------------------------------------
# Plots
## Plotting within subject relationships
set.seed(1)
randSub<-subset(FILTER, subID %in% sample(unique(FILTER$subID), size = 3))

## Plots of lnRT by Reward.
g1<-ggplot(data=randSub, aes(x = Reward, y = lnRT, group = subID))+geom_line()
g2<-g1+geom_point()+facet_wrap(~subID)+theme(text=element_text(size=20), panel.background=element_rect(fill="white"))
g3<-g2+geom_smooth(method = "lm", se=TRUE, color="blue", aes(group=subID))
g4<-g3 + scale_y_continuous(name = "ln(RT)") + scale_x_continuous(name = "Incentives ($)")
plot(g4)

## Plots of lnRT by Beta.
g1<-ggplot(data=randSub, aes(x = Beta, y = lnRT, group = subID))+geom_line()
g2<-g1+geom_point()+facet_wrap(~subID)+theme(text=element_text(size=20), panel.background=element_rect(fill="white"))
g3<-g2+geom_smooth(method = "lm", se=TRUE, color="green", aes(group=subID))
g4<-g3 + scale_y_continuous(name = "ln(RT)") + scale_x_continuous(name = "ln(Beta Power)")
plot(g4)

## Plots of Beta by Reward.
g1<-ggplot(data=randSub, aes(x = Reward, y = Beta, group = subID))+geom_line()
g2<-g1+geom_point()+facet_wrap(~subID)+theme(text=element_text(size=20), panel.background=element_rect(fill="white"))
g3<-g2+geom_smooth(method = "lm", se=TRUE, color="red", aes(group=subID))
g4<-g3 + scale_y_continuous(name = "ln(Beta Power)") + scale_x_continuous(name = "Incentives ($)")
plot(g4)
###########################################################################


## -----------------------------------------------------------------------
# Linear Mixed-Effect Regression Models
# Null Model controlling for the nested nature of the data
M0<-lmer(lnRT~1+(1|subID),data=FILTER, REML=FALSE)
summary(M0)
plot(resid(M0))
## Note that this random-intercepts model merely establishes a baseline.
## All critical model comparisons are between random-slopes models.


#Unmediated Path C
M1<-lmer(lnRT~1+Reward.c+(1+Reward.c|subID),data=FILTER, REML=FALSE)
summary(M1)
plot(density(resid(M1)))
plot(resid(M1)~fitted(M1))

anova(M0,M1)

# Unmediated Path A
M2<-lmer(Beta.c~1+Reward.c+(1+Reward.c|subID),data=FILTER, REML=FALSE)
summary(M2)

####################################################################
# Unmediated Path A by Electrode
## There is some concern about averaging across multiple electrodes
## for our calculation of beta power. Thus, we rerun the same model
## (M2) at each of the electrodes independently.
## There are no significant relationships between beta power and 
## incentive rewards at any of the electrodes.
FC5<-lmer(FC5~1+Reward.c+(1+Reward.c|subID),data=FILTER, REML=FALSE)
summary(FC5)

FC3<-lmer(FC3~1+Reward.c+(1+Reward.c|subID),data=FILTER, REML=FALSE)
summary(FC3)

FC1<-lmer(FC1~1+Reward.c+(1+Reward.c|subID),data=FILTER, REML=FALSE)
summary(FC1)

C5<-lmer(C5~1+Reward.c+(1+Reward.c|subID),data=FILTER, REML=FALSE)
summary(C5)

C3<-lmer(C3~1+Reward.c+(1+Reward.c|subID),data=FILTER, REML=FALSE)
summary(C3)

C1<-lmer(C1~1+Reward.c+(1+Reward.c|subID),data=FILTER, REML=FALSE)
summary(C1)

CP5<-lmer(CP5~1+Reward.c+(1+Reward.c|subID),data=FILTER, REML=FALSE)
summary(CP5)

CP3<-lmer(CP3~1+Reward.c+(1+Reward.c|subID),data=FILTER, REML=FALSE)
summary(CP3)

CP1<-lmer(CP1~1+Reward.c+(1+Reward.c|subID),data=FILTER, REML=FALSE)
summary(CP1)
#####################################################################



# Unmediated Path B
M3<-lmer(lnRT~1+Beta.c+(1+Beta.c|subID),data=FILTER, REML=FALSE)
summary(M3)

# Mediated Path C
M4<-lmer(lnRT~1+Reward.c+Beta.c+(1+Reward.c+Beta.c|subID),data=FILTER, REML=FALSE)
summary(M4)
plot(resid(M4))
plot(density(resid(M4)))
plot(resid(M4)~fitted(M4))

# Model Comparisons
anova(M1,M4)
anova(M3,M4) #adding both parameters as a set

# Approximating R-squared based on Verbeke and MolenBerghs (2000)
SST<-sum(residuals(M0)^2)
SST
SSR<-sum(residuals(M4)^2)
SSR
r2meta<-1-(SSR/SST)
r2meta


# It looks more like Reward and Beta Suppression make independent contributions to RT
# We will control for their interaction in the moderator model, M5

#Moderator Model
M5<-lmer(lnRT~1+Reward.c*Beta.c+(1+Reward.c+Beta.c|subID),data=FILTER, REML=FALSE)
summary(M5)

anova(M0,M1,M4,M5)

## -----------------------------------------------------------------------
# #Bootstrap analysis of the mediation path:
# #model.a <- lmer(M ~ X + (1|Subject), data=dat)
# model.a<-M2
# #model.b <- lmer(Y ~ M + (1|Subject), data=dat)
# model.b<-M3
# #model.y <- lmer(Y ~ X + (1|Subject), data=dat)
# model.y<-M1
# #model.m <- lmer(Y ~ X + M + (1|Subject), data=dat)
# model.m<-M4
# #saveToLocation <- "./mixedBoot.txt"
# nBoot <- 5000
# 
# bout <- matrix(nrow = nBoot,ncol = 7)
# for(i in 1:nBoot){
#   bout[i,1] <- fixef(refit(model.a, simulate(model.a, nsim = 1 , seed = NULL)))[2]
#   bout[i,2] <- fixef(refit(model.b, simulate(model.b, nsim = 1 , seed = NULL)))[2]
#   bout[i,3] <- fixef(refit(model.y, simulate(model.y, nsim = 1 , seed = NULL)))[2]
#   bout[i,4:5] <- fixef(refit(model.m, simulate(model.m, nsim = 1 , seed = NULL)))[2:3]
# }
# bout[,6] <- bout[,3]-bout[,4]
# bout[,7] <- bout[,1]*bout[,2]
# bout <- data.frame(bout)
# colnames(bout) <-c('a','b','c','cprime','m','c_cprime','ab')
# #save(bout, file=saveToLocation)
# write.csv(bout, file="mixedBoot.csv")
# 
# probs <- c(0.025, 0.975)
# quantile(bout$c_cprime, probs = probs)
# #######################################################################################


## -----------------------------------------------------------------------
# Effects of Plausability.
# Note that this removes participant 5 from the model (didn't answer question)
FILTER2<-subset(FILTER,!Plaus == 'NA')
#Descriptive stats for Plausability
summary(FILTER2$Plaus)
FILTER2$Plaus.c<-FILTER2$Plaus-0.6044

M6<-lmer(lnRT~1+Reward.c+Beta.c+Plaus.c+(1+Reward.c+Beta.c|subID),data=FILTER2, REML=FALSE)
summary(M6)

M7<-lmer(lnRT~1+Reward.c*Beta.c*Plaus.c+(1+Reward.c+Beta.c|subID),data=FILTER2, REML=FALSE)
summary(M7)



