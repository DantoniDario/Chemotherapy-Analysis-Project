library(nlme)
library(scales)
library(lme4)
library(lmerTest)
library(nlme)
library(ggplot2)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(MuMIn)
library(ggpubr)
library(ggcorrplot)
library(car)
library(robustlmm)

## A color function
gg_color <- function(n, alpha = 1) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100, alpha = alpha)[1:n]
}
Chemotherapy.version1 <- read.csv("~/Desktop/Mixed Linear Models/Documents-20201215/Chemotherapy-version1.csv")
df = Chemotherapy.version1

df$line[df$line >= 3] = 3
df[,c(1,3)] = lapply(df[,c(1,3)], as.factor)

# Visual description of the data
ggplot(df, aes(x=line, y=tumour)) + 
  geom_boxplot(fill="slateblue", alpha=0.5) +
  ggtitle("Tumour over line")
ggplot(df, aes(x=patient, y=tumour)) + 
  geom_boxplot(fill="slateblue", alpha=0.5) +
  ggtitle("Tumour over patient")
ggplot(df, aes(x=as.factor(month), y=tumour)) + # hint for month^2 --> banana shape
  geom_boxplot(fill="slateblue", alpha=0.5) +
  ggtitle("Tumour over month")
ggplot(df, aes(x=line, y=sensitivity)) + 
  geom_boxplot(fill="slateblue", alpha=0.5) +  # hint that over time/line, sensitivity doesn't
  ggtitle("Sensitivity over line")  #  explain well tumour --> negative slope while it is positive for tumour~line
ggplot(df, aes(x=as.factor(month), y=sensitivity)) + 
  geom_boxplot(fill="slateblue", alpha=0.5) +  # hint that over time/line, sensitivity doesn't
  ggtitle("Sensitivity over month") 

# Check for significant change of variance
summary(aov(tumour ~ patient, data = df))
summary(aov(tumour ~ line, data = df)) # We should maybe consider heteroskedastic model for each level of line
summary(aov(tumour ~ as.factor(month), data = df))
summary(aov(tumour ~ sensitivity, data = df))

# Average mean for different level of categoricals variable
plot.design(tumour~., main="plot.design",ylim=range(df$tumour),data=df)

# Interaction plot
cols = gg_color(nlevels(df$patient))
interaction.plot(df$line,df$patient,df$tumour,col=cols,main="interaction.plot",lwd=1.1,lty=1)
cols = gg_color(nlevels(df$line))
interaction.plot(df$month,df$line,df$tumour,col=cols,main="interaction.plot",lwd=1.1,lty=1)
cols = gg_color(nlevels(df$patient))
interaction.plot(df$month,df$patient,df$tumour,col=cols,main="interaction.plot",lwd=1.1,lty=1)

# Plot with grouped variables
# Same as ggplot (below) but I let the code here in case you prefer with that format
#Tumour.gD <- groupedData(tumour~line|patient,data=df, order.groups = T)
#plot(Tumour.gD)
#plot(Tumour.gD, outer = ~line) # easier to see here that maybe tumour is higher for the third line
#Tumour.gD <- groupedData(sensitivity~patient|tumour,data=df, order.groups = T)
#plot(Tumour.gD, outer = ~patient)

df %>%
  mutate(patient = fct_reorder(patient, tumour, .fun='max')) %>% # ordered by tumour (max fun function)
  ggplot(aes(x=tumour, y=patient, col = line)) +
  geom_point(fill="slateblue", alpha=0.7) +
  ggtitle("Patient over tumour (per line)") +
  xlab("Tumour") +
  ylab("Patient") +
  guides(col=guide_legend("Line"))

df %>%
  mutate(patient = fct_reorder(patient, tumour, .fun='max')) %>% # ordered by tumour (max fun function)
  ggplot(aes(x=tumour, y=patient, col = line)) +
  geom_point(fill="slateblue", alpha=0.7) +
  ggtitle("Patient over tumour (per line)") +
  xlab("Tumour") +
  ylab("Patient") +
  facet_wrap(line~., ncol = 2) +
  guides(col=guide_legend("Line"))

df %>%
  mutate(patient = fct_reorder(patient, tumour, .fun='max')) %>% # ordered by tumour (max fun function)
  ggplot(aes(x=tumour, y=patient, col = as.factor(month))) +
  geom_point(fill="slateblue", alpha=0.7) +
  ggtitle("Patient over Tumour (per month)") +
  xlab("Tumour") +
  ylab("Patient") +
  facet_wrap(month~., ncol = 2) +
  guides(col=guide_legend("Month"))

# Same as above but not ordered by tumour
ggplot(df, aes(x=tumour, y=patient, col = as.factor(month))) +
  geom_point(fill="slateblue", alpha=0.7) +
  ggtitle("Patient over sensitivity (per month)") +
  xlab("Tumour") +
  ylab("Patient") +
  facet_wrap(month~., ncol = 2) +
  guides(col=guide_legend("Month"))

ggplot(df, aes(x=sensitivity, y=tumour, col = as.factor(month))) +
  geom_point(fill="slateblue", alpha=0.7) +
  ggtitle("Tumour over sensitivity (per month)") +
  xlab("Sensitivity") +
  ylab("Tumour") +
  facet_wrap(month~., ncol = 2) +
  geom_smooth(method=lm, size=0.5, linetype="dashed"  , color="red", se=TRUE) +
  guides(col=guide_legend("Month"))

ggplot(df, aes(x=sensitivity, y=tumour, col = line)) +
  geom_point(fill="slateblue", alpha=0.7) +
  ggtitle("Tumour over sensitivity (per line)") +
  xlab("Sensitivity") +
  ylab("Tumour") +
  facet_wrap(line~., ncol = 2) +
  geom_smooth(method=lm, size=0.5, linetype="dashed"  , color="red", se=TRUE) +
  guides(col=guide_legend("Line"))

# For question 2, this graph (woth the code below) could answer to the question already because we can see that
# higher is line, smaller is the sensitivity slope which means that its explanatory power tumour
# become less significant. We need just to test statistically that the difference in slope is
# different

ggplot(df, aes(x=sensitivity, y=tumour, col = line)) +
  geom_point(fill="slateblue", alpha=0.7) +
  ggtitle("Tumour over sensitivity (per line)") +
  xlab("Sensitivity") +
  ylab("Tumour") +
  geom_smooth(aes(sensitivity,tumour, colour=line), size=0.5, linetype="dashed" , method=lm, se=TRUE) +
  guides(col=guide_legend("Line"))

ggplot(df, aes(x=sensitivity, y=tumour)) +
  geom_point(fill="slateblue", alpha=0.7) +
  ggtitle("Tumour over sensitivity (per line)") +
  xlab("Sensitivity") +
  ylab("Tumour") +
  geom_smooth( color='red', size=0.5, linetype="dashed" , method=lm, se=TRUE)

ggplot(df, aes(x=sensitivity, y=tumour)) +
  geom_point(fill="slateblue", alpha=0.5) +
  ggtitle("tumour over sensitivity (per patient)") +
  xlab("Sensitivity") +
  ylab("tumour") +
  ylim(-4, 4) +
  facet_wrap(patient~.) +
  geom_smooth(method=lm, size=0.5, linetype="dashed"  , color="red", se=TRUE)

ggplot(df, aes(x=month, y=tumour, color= line)) +
  geom_point(alpha=0.5) +
  ggtitle("tumour over month (per patient)") +
  xlab("month") +
  ylab("tumour") +
  facet_wrap(patient~.) +
  geom_smooth(method=lm, size=0.5, linetype="dashed"  , color="red", se=TRUE)

# Check for nomrality of the response (they have done this in the best group last year)

ggplot(df, aes(sample = tumour)) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("QQ-plot of tumour") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")

# Try to transform the data, we know they are log scale (the prof said it), let's try with exp or log again
ggplot(df, aes(sample = exp(tumour))) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("QQ-plot of tumour") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")

ggplot(df, aes(sample = log(tumour+4))) + # +4 because we have negative tumour -3.8,
  stat_qq() +                             # log doesn't accept negative values
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("QQ-plot of tumour") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")
# definitvely not a good idea

ggplot(df, aes(sample = tumour)) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("QQ-plot of tumour per line") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles") +
  facet_wrap(line~., ncol = 2)

ggplot(df, aes(sample = tumour)) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("QQ-plot of tumour per month") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles") +
  facet_wrap(month~., ncol = 2)

# I tried to build the most complex model while making sense as well (the goal here is not to add
# useless interactions without any sense)
mod1 = lme(fixed = tumour~sensitivity*line + month*line, random=list(patient=pdBlocked(list(
  pdIdent(~1), pdIdent(~sensitivity-1),pdIdent(~month-1),pdIdent(~line-1)))), data = df)

# Here another version of the first model but considering a correlation (≠ 0) between random intercept and slope effect for sensitivity-patient
mod1.ml = lme(fixed = tumour~sensitivity*line + month*line, random=list(patient=pdBlocked(list(
  pdSymm(~sensitivity),pdIdent(~month-1),pdIdent(~line-1)))), data = df, method = "ML")
mod1.ml = lme(fixed = tumour~sensitivity*line + month*line, random=list(patient=pdBlocked(list(
  pdDiag(~sensitivity),pdIdent(~month-1),pdIdent(~line-1)))), data = df, method = "ML")

# mod1 but I retrieve interaction between patient-line, patient-month, and month-line (this 
# interaction is not in our interest to estimate)
mod2.reml = lme(fixed = tumour~sensitivity*line + month, random=~sensitivity|patient, data = df, method = "REML")
mod2.ml = lme(fixed = tumour~sensitivity*line + month, random=~sensitivity|patient, data = df, method = "ML")

# Other version of mod2 (intercation between line and month)
mod2.reml = lme(fixed = tumour~sensitivity*line + month+line, random=~1|patient, data = df, method = "ML")
mod2.ml = lme(fixed = tumour~sensitivity*line + month*line, random=~1|patient, data = df, method = "ML")
mod2.ml.lme4 = lmer(tumour~sensitivity*line + month*line + (1|patient), REML = F, data = df)
mod2.reml.lme4 = lmer(tumour~sensitivity*line + month+line + (1|patient), REML = F, data = df)

####
n.r = 1000
mx.t.nr = simulate(mod2.reml.lme4, nsim=n.r)
pval.r = Chisq.r = lrt.r = rep(NA,n.r)
data.r <- df

for (rw in 1:n.r){
  data.r$tumour = mx.t.nr[,rw]
  fit.restr = lmer(tumour~sensitivity*line + month+line + (1|patient), REML = F, data = df)
  fit.full = lmer(tumour~sensitivity*line + month*line + (1|patient), REML = F, data = df)
  lrt.r[rw] = as.numeric(2*(logLik(fit.full) - logLik(fit.restr)))
  anova.r = anova(fit.restr, fit.full)
  pval.r[rw] = anova.r$Pr[2]
  Chisq.r[rw] = anova.r$Chisq[2]
}
AIC(mod2.reml, mod2.ml)
anova(mod2.ml.lme4, mod2.reml.lme4)
summary(mod2.reml.lme4)
mean(pval.r<0.05)

qbinom(c(0.025,0.975), prob = 0.05, size = 1000)/1000

LRT.real = as.numeric(2*(logLik(mod2.ml.lme4) - logLik(mod2.reml.lme4)))
mean(LRT.real<lrt.r)

pbkrtest::KRmodcomp(mod2.ml.lme4,mod2.reml.lme4)

Chisq.r.dens2 = as.data.frame(Chisq.r)



####


summary(mod2.reml)
summary(mod2.ml.lme4)
anova(mod2.reml, mod2.ml)
pbkrtest::KRmodcomp(mod2.reml,mod2.ml)


# mod2 but I retrieve correlation between random intercept and slope
mod3.reml = lme(fixed = tumour~sensitivity*line + month, random=list(patient=pdDiag(~sensitivity)), data = df, method = "REML")
mod3.ml = lme(fixed = tumour~sensitivity*line + month, random=list(patient=pdDiag(~sensitivity)), data = df, method = "ML")

# mod3 but I retrieve random slope
mod4.reml = lme(fixed = tumour~sensitivity*line + month, random=~1|patient, data = df, method = "REML")
mod4.ml = lme(fixed = tumour~sensitivity*line + month, random=~1|patient, data = df, method = "ML")

# summary of the model
# Clearly we can see that month and line are very small random effect, they could be neglected
summary(mod1)
summary(mod1.reml)
summary(mod2.reml)
summary(mod3.reml) # Random slope become very small, we can say it is not significant I think
summary(mod4.reml)

# Should I try another kind test? Are we allowed that kinf of test?
Anova(mod1.reml, type="III" )
Anova(mod2.reml, type="III" )
Anova(mod3.reml, type="III" )
Anova(mod4.reml, type="III" )

anova(mod1.reml, type="marginal" )
anova(mod2.reml, type="marginal" )
anova(mod3.reml, type="marginal" )
anova(mod4.reml, type="marginal" )

anova(mod1.reml, type="sequential" )
anova(mod2.reml, type="sequential" )
anova(mod3.reml, type="sequential" )
anova(mod4.reml, type="sequential" )

# Check consistency of estimator (with reml, the code runs well for only 2 models.. with ml, it runs for 3 models)
intervals(mod1.ml)
intervals(mod2.ml)
intervals(mod3.ml) # doesn't work
intervals(mod4.ml)

# Can we use this  check consistency? Except line and month, for mod1.reml the other variable
# doesn't seem to change a lot

mod1.reml$coefficients$fixed
mod2.reml$coefficients$fixed
mod3.reml$coefficients$fixed
mod4.reml$coefficients$fixed

mod1.reml$coefficients$random
mod2.reml$coefficients$random
mod3.reml$coefficients$random
mod4.reml$coefficients$random

# Residuals check for normality
# Same as ggplot (below) but I let the code here in case you prefer with that format
#qqnorm(resid(mod1),main="Residual\nQ-Q plot",pch=16)
#qqline(resid(mod1),col=2,lwd=2,lty=2)

qq1 = ggplot(df, aes(sample = resid(mod1.reml, type = "p"))) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")

qq2 = ggplot(df, aes(sample = resid(mod2.reml, type = "p"))) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")

qq3 = ggplot(df, aes(sample = resid(mod3.reml, type = "p"))) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")

qq11 = ggplot(df, aes(sample = resid(mod9.ml, type = "p"))) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")

ggarrange(qq1, qq2, qq3, qq4, 
          labels = c("QQ-plot (mod1)", "QQ-plot (mod2)", "QQ-plot (mod3)", "QQ-plot (mod4)"),
          ncol = 2, nrow = 2)

## SIMULATION of normal with n = 137, this step to be repeted 8 times in changint qqa, qqb, ..., qqh each time

qq_df = as.data.frame(rnorm(137, 0, mod9.ml$sigma))
qqa = ggplot(qq_df, aes(sample = `rnorm(137, 0, mod9.ml$sigma)`)) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")

qq_df = as.data.frame(rnorm(137, 0, mod9.ml$sigma))
qqb = ggplot(qq_df, aes(sample = `rnorm(137, 0, mod9.ml$sigma)`)) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")

qq_df = as.data.frame(rnorm(137, 0, mod9.ml$sigma))
qqc = ggplot(qq_df, aes(sample = `rnorm(137, 0, mod9.ml$sigma)`)) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")

qq_df = as.data.frame(rnorm(137, 0, mod9.ml$sigma))
qqd = ggplot(qq_df, aes(sample = `rnorm(137, 0, mod9.ml$sigma)`)) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")

qq_df = as.data.frame(rnorm(137, 0, mod9.ml$sigma))
qqe = ggplot(qq_df, aes(sample = `rnorm(137, 0, mod9.ml$sigma)`)) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")

qq_df = as.data.frame(rnorm(137, 0, mod9.ml$sigma))
qqf = ggplot(qq_df, aes(sample = `rnorm(137, 0, mod9.ml$sigma)`)) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")

qq_df = as.data.frame(rnorm(137, 0, mod9.ml$sigma))
qqg = ggplot(qq_df, aes(sample = `rnorm(137, 0, mod9.ml$sigma)`)) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")

qq_df = as.data.frame(rnorm(137, 0, mod9.ml$sigma))
qqh = ggplot(qq_df, aes(sample = `rnorm(137, 0, mod9.ml$sigma)`)) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")
summary(mod9.ml)
a = ggarrange(qqa, qqb, qqc, qqd, qq11, qqe, qqf, qqg, qqh,
          ncol = 3, nrow = 3)
annotate_figure(a,top = text_grob(expression("Simulated residuals"~ italic(N)(0, sigma[epsilon]^2) ~"versus fitted residuals"), color = "black", face = "bold", size = 15)) 


# Residuals check for heteroskeasticity or any pattern

# Same as ggplot (below) but I let the code here in case you prefer with that format
#plot(mod1,resid(.,type="p")~fitted(.)|line,abline=0,id=0.05)
#plot(mod1,resid(.,type="p")~fitted(.)|month,abline=0,id=0.05) 
#plot(mod1,resid(.,type="p")~fitted(.)|patient,abline=0,id=0.05)

gg1 = ggplot(df, aes(x=fitted(mod1.reml), y=resid(mod1.reml, type = "p"))) +
  geom_point(fill="slateblue", alpha=0.5) +
  ggtitle("") +
  xlab("Fitted tumour") + 
  ylab("Pearson residuals") +
  geom_smooth(method=lm, size=0.5, linetype="dashed" , color="red", se=TRUE)

gg2 = ggplot(df, aes(x=fitted(mod2.reml), y=resid(mod2.reml, type = "p"))) +
  geom_point(fill="slateblue", alpha=0.5) +
  ggtitle("") +
  xlab("Fitted tumour") + 
  ylab("Pearson residuals") +
  geom_smooth(method=lm, size=0.5, linetype="dashed" , color="red", se=TRUE)

gg3 = ggplot(df, aes(x=fitted(mod3.reml), y=resid(mod3.reml, type = "p"))) +
  geom_point(fill="slateblue", alpha=0.5) +
  ggtitle("") +
  xlab("Fitted tumour") + 
  ylab("Pearson residuals") +
  geom_smooth(method=lm, size=0.5, linetype="dashed" , color="red", se=TRUE)

gg4 = ggplot(df, aes(x=fitted(mod4.reml), y=resid(mod4.reml, type = "p"))) +
  geom_point(fill="slateblue", alpha=0.5) +
  ggtitle("") +
  xlab("Fitted tumour") + 
  ylab("Pearson residuals") +
  geom_smooth(method=lm, size=0.5, linetype="dashed" , color="red", se=TRUE)

ggarrange(gg1, gg2, gg3, gg4, 
          labels = c("mod1 (Residuals over month)", "mod2 (Residuals over month)", "mod3 (Residuals over month)", "mod4 (Residuals over month)"),
          ncol = 2, nrow = 2)

ggplot(df, aes(x=fitted(mod1.reml), y=resid(mod1.reml, type = "p"))) +
  geom_point(fill="slateblue", alpha=0.5) +
  ggtitle("Residuals over line") +
  xlab("Fitted tumour") +
  ylab("Pearson residuals") +
  facet_wrap(line~., ncol = 2) +
  geom_smooth(method=lm, size=0.5, linetype="dashed" , color="red", se=TRUE)

ggplot(df, aes(x=fitted(mod1.reml), y=resid(mod1.reml, type = "p"))) +
  geom_point(fill="slateblue", alpha=0.5) +
  ggtitle("Residuals over month") +
  xlab("Fitted tumour") +
  ylab("Pearson residuals") +
  facet_wrap(month~., ncol = 2) +
  geom_smooth(method=lm, size=0.5, linetype="dashed" , color="red", se=TRUE)

ggplot(df, aes(x=fitted(mod1.reml), y=resid(mod1.reml, type = "p"))) +
  geom_point(fill="slateblue", alpha=0.5) +
  ggtitle("Residuals over patient") +
  xlab("Fitted tumour") +
  ylab("Pearson residuals") +
  facet_wrap(patient~.) +
  geom_smooth(method=lm, size=0.5, linetype="dashed" , color="red", se=TRUE)

# Residuals check for pattern
ggplot(df, aes(x=line, y=resid(mod1.reml, type = "p"))) +
  geom_boxplot(fill="slateblue", alpha=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  ggtitle("Residuals over line") +
  ylab("Pearson residuals")

# Maybe heteroskedasticity and banana shape
ggplot(df, aes(x=as.factor(month), y=resid(mod1.reml, type = "p"))) +
  geom_boxplot(fill="slateblue", alpha=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  xlab("Month") +
  ggtitle("Residuals over month") +
  ylab("Pearson residuals")

ggplot(df, aes(x=patient, y=resid(mod1.reml, type = "p"))) +
  geom_boxplot(fill="slateblue", alpha=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  ggtitle("Residuals over patient") +
  ylab("Pearson residuals")

ggplot(df, aes(x=patient, y=resid(mod1.reml, type = "p"))) +
  geom_boxplot(fill="slateblue", alpha=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  ggtitle("Residuals over patient and month") +
  ylab("Pearson residuals") +
  facet_wrap(month~.)

ggplot(df, aes(x=patient, y=resid(mod1.reml, type = "p"))) +
  geom_boxplot(fill="slateblue", alpha=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  ggtitle("Residuals over patient and line") +
  ylab("Pearson residuals") +
  facet_wrap(line~.)

# Check for correlation between random effect
pairs(mod1.reml,~ranef(.),id=0.05,grid=T,abline=0)
pairs(mod2.reml,~ranef(.),id=0.05,grid=T,abline=0)
pairs(mod3.reml,~ranef(.),id=0.05,grid=T,abline=0) # Remember the random slope here is not significant
pairs(mod4.reml,~ranef(.),id=0.05,grid=T,abline=0) # doesn't work since we have only one random effect

# GGplot version of interraction between 2 random effects (mod2)
ggplot(data = ranef(mod2.reml),aes(y = ranef(mod2.reml)[,2], x = ranef(mod2.reml)[,1])) +
  geom_point() +
  geom_smooth(method=lm, size=0.5, linetype="dashed"  , color="red", se=TRUE) +
  theme_minimal() + 
  xlab("Random Intercept") + 
  ylab("Random Slope")

# Normality check for random effect
par(mfrow=c(2,3))
qqnorm(ranef(mod1.reml)[,1],pch=16,main="Intercept")
qqline(ranef(mod1.reml)[,1],col=2,lwd=2,lty=2)
qqnorm(ranef(mod1.reml)[,2],pch=16,main="sensitivity")
qqline(ranef(mod1.reml)[,2],col=2,lwd=2,lty=2)
qqnorm(ranef(mod1.reml)[,3],pch=16,main="month")
qqline(ranef(mod1.reml)[,3],col=2,lwd=2,lty=2)
qqnorm(ranef(mod1.reml)[,4],pch=16,main="line1")
qqline(ranef(mod1.reml)[,4],col=2,lwd=2,lty=2)
qqnorm(ranef(mod1.reml)[,5],pch=16,main="line2")
qqline(ranef(mod1.reml)[,5],col=2,lwd=2,lty=2)
qqnorm(ranef(mod1.reml)[,6],pch=16,main="line3")
qqline(ranef(mod1.reml)[,6],col=2,lwd=2,lty=2)
par(mfrow=c(1,1))

# I am not sure but I think I succeed to do the ggplot version of the quantile code above

qq5 = ggplot(ranef(mod5_4.ml), aes(sample = ranef(mod5_4.ml)[,1])) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")

qq6 = ggplot(ranef(mod1.reml), aes(sample = ranef(mod1.reml)[,2])) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")

qq7 = ggplot(ranef(mod1.reml), aes(sample = ranef(mod1.reml)[,3])) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")

qq8 = ggplot(ranef(mod9.ml), aes(sample = ranef(mod9.ml)[,1])) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")

qq9 = ggplot(ranef(mod1.reml), aes(sample = ranef(mod1.reml)[,5])) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")

qq10 = ggplot(ranef(mod1.reml), aes(sample = ranef(mod1.reml)[,6])) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")

ggarrange(qq5, qq6, qq7, qq8, qq9, qq10, 
          labels = c("Intercept", "Sensitivity", "month", "line 1", "line 2", "line 3"),
          ncol = 3, nrow = 2)


plot(nlme:::ACF.lme(mod5.ml,type="pearson"),alpha=0.05)
plot(nlme:::ACF.lme(mod2.reml,type="pearson"),alpha=0.05)
plot(nlme:::ACF.lme(mod3.reml,type="pearson"),alpha=0.05)
plot(nlme:::ACF.lme(mod4.reml,type="pearson"),alpha=0.05) # weird it shows only 11 lags but we should have
# at maximum 14 observation per patient
plot(nlme:::ACF.lme(mod1.reml,type="pearson", maxLag = 15),alpha=0.05)
plot(nlme:::ACF.lme(mod2.reml,type="pearson", maxLag = 15),alpha=0.05)
plot(nlme:::ACF.lme(mod3.reml,type="pearson", maxLag = 15),alpha=0.05)
plot(nlme:::ACF.lme(mod4.reml,type="pearson", maxLag = 15),alpha=0.05) # here a version with more lags

# Conclusion: some spikes are above the confidence interval but it doesn't seem alarming

round(r.squaredGLMM(mod1.reml),4)
round(r.squaredGLMM(mod2.reml),4)
round(r.squaredGLMM(mod3.reml),4)
round(r.squaredGLMM(mod4.reml),4)

# mod4 but with month^2
mod3.ml = lme(fixed = tumour~sensitivity*line + month, random=~1|patient, data = df, method = "ML")

mod5.ml = lme(fixed = tumour~sensitivity*line + month + I(month^2), random=~1|patient, data = df, method = "ML")

anova(mod3.ml, mod5.ml)
mod2.reml
# banana shape has disapeared
ggplot(df, aes(x=as.factor(month), y=resid(mod5.ml, type = "p"))) +
  geom_boxplot(fill="slateblue", alpha=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  xlab("Month") +
  ggtitle("Residuals over month") +
  ylab("Pearson residuals")

# Check that mod5.ml is the best compared to restricted other model

mod5_1.ml = lme(fixed = tumour~sensitivity*line, random=~1|patient, data = df, method = "ML")
mod5_2.ml = lme(fixed = tumour~sensitivity*line + month, random=~1|patient, data = df, method = "ML")
mod5_3.ml = lme(fixed = tumour~sensitivity*line + month + I(month^2), random=list(patient=pdDiag(~sensitivity)), data = df, method = "ML")
mod5_31.ml = lme(fixed = tumour~sensitivity*line + month + I(month^2), random=list(patient=pdSymm(~sensitivity)), data = df, method = "ML")

mod5_4.ml = lme(fixed = tumour~sensitivity*line + line:month + I(month^2), random=~1|patient, data = df, method = "ML")
modbest.ml = lme(fixed = tumour~sensitivity*line + line:month + I(month^2), random=~1|patient, data = df, method = "ML")


summary(mod5_3.ml)
anova(mod5.ml, mod5_1.ml)
anova(mod5.ml, mod5_2.ml)
anova(mod5.ml, mod5_3.ml)
anova(mod5.ml, mod5_31.ml)
anova(mod5.ml, mod5_4.ml)


summary(mod6.ml)

anova(mod5_3.ml, mod5_31.ml)
summary(modbest.ml)

# mod5 with hetero
mod6.ml = lme(fixed = tumour~sensitivity*line + line:month + 
                I(month^2), random=~1|patient, data = df, weights=varPower(form=~sensitivity|line), method = "ML")
mod7.ml = lme(fixed = tumour~sensitivity*line + line:month + 
                I(month^2), random=~1|patient, data = df, weights=varIdent(form=~1|line), method = "ML")
mod8.ml = lme(fixed = tumour~sensitivity*line + line:month + 
                I(month^2), random=~1|patient, data = df, weights=varPower(form=~month|line), method = "ML")
mod9.ml = lme(fixed = tumour~sensitivity*line + month + 
                I(month^2), random=~1|patient, data = df, weights=varExp(form=~sensitivity), method = "ML")

summary(mod9.ml)
# Check for visual improvement
ggplot(df, aes(x=sensitivity, y=resid(mod5.ml.hetero, type = "p"))) +
  geom_point(fill="slateblue", alpha=0.5) +
  ggtitle("Residuals over month") +
  xlab("Fitted tumour") +
  ylab("Pearson residuals") +
  facet_wrap(line~., ncol = 2) +
  geom_smooth(method=lm, size=0.5, linetype="dashed" , color="red", se=TRUE)

gg5 = ggplot(df, aes(x=sensitivity, y=resid(mod5.ml, type = "p"))) +
  geom_point(fill="slateblue", alpha=0.5) +
  ggtitle("Residuals over sensitivity") +
  xlab("Sensitivity") +
  ylab("Pearson residuals") +
  geom_smooth(method=lm, size=0.5, linetype="dashed" , color="red", se=TRUE) +
  theme(plot.title = element_text(size=10))

gg6 = ggplot(df, aes(x=sensitivity, y=resid(mod5_4.ml, type = "p"))) +
  geom_point(fill="slateblue", alpha=0.5) +
  ggtitle("Residuals over sensitivity \n(heteroskedasticity modelled)") +
  xlab("Sensitivity") +
  ylab("Pearson residuals") +
  geom_smooth(method=lm, size=0.5, linetype="dashed" , color="red", se=TRUE) +
  theme(plot.title = element_text(size=10))

gg7 = ggplot(df, aes(x=line, y=resid(mod5.ml, type = "p"))) +
  geom_boxplot(fill="slateblue", alpha=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  ggtitle("Residuals over line") +
  xlab("Line") +
  ylab("Pearson residuals") +
  theme(plot.title = element_text(size=10))

gg8 = ggplot(df, aes(x=line, y=resid(mod5_4.ml, type = "p"))) +
  geom_boxplot(fill="slateblue", alpha=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  ggtitle("Residuals over line \n(heteroskedasticity modelled)") +
  xlab("Line") +
  ylab("Pearson residuals") +
  theme(plot.title = element_text(size=10))

gg9 = ggplot(df, aes(x=as.factor(month), y=resid(mod5.ml, type = "p"))) +
  geom_boxplot(fill="slateblue", alpha=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  ggtitle("Residuals over month") +
  ylab("Pearson residuals") +
  xlab('Month') +
  theme(plot.title = element_text(size=10))

gg10 = ggplot(df, aes(x=as.factor(month), y=resid(mod5_4.ml, type = "p"))) +
  geom_boxplot(fill="slateblue", alpha=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  ggtitle("Residuals over month \n(heteroskedasticity modelled)") +
  ylab("Pearson residuals") +
  xlab('Month') +
  theme(plot.title = element_text(size=10))

ggarrange(gg5, gg7, gg9, gg6, gg8, gg10,
          ncol = 3, nrow = 2)
summary(mod9.ml)
ggplot(df, aes(x=fitted(mod5_4.ml), y=resid(mod5_4.ml, type = "p"))) +
  geom_point(fill="slateblue", alpha=0.5) +
  ggtitle("Residuals over sensitivity") +
  xlab("Sensitivity") +
  ylab("Pearson residuals") +
  geom_smooth(method=lm, size=0.5, linetype="dashed" , color="red", se=TRUE) +
  theme(plot.title = element_text(size=10))

# don't remember what we discussd but I think the exact text is not suitable here no?
anova(modbest.ml, mod6.ml)   # we should boostrap the distribution maybe
anova(modbest.ml, mod7.ml)
anova(modbest.ml, mod8.ml) 
anova(mod6.ml, mod9.ml) # This one is significantly better than mod5
# ARMA without hetero

mod9.ml = lme(fixed = tumour~sensitivity*line + line:month + 
                I(month^2), random=~1|patient, data = df, corr=corAR1(form=~1|patient), method = "ML")
mod10.ml = lme(fixed = tumour~sensitivity*line + line:month + 
                 I(month^2), random=~1|patient, data = df, corr=corARMA(form=~1|patient,p=2,q=0), method = "ML")
mod11.ml = lme(fixed = tumour~sensitivity*line + line:month + 
                 I(month^2), random=~1|patient, data = df, corr=corARMA(form=~1|patient,p=2,q=1), method = "ML")
mod12.ml = lme(fixed = tumour~sensitivity*line + line:month + 
                 I(month^2), random=~1|patient, data = df, corr=corARMA(form=~1|patient,p=2,q=2), method = "ML") # doesn't converge
mod13.ml = lme(fixed = tumour~sensitivity*line + line:month + 
                 I(month^2), random=~1|patient, data = df, corr=corARMA(form=~1|patient,p=0,q=2), method = "ML")
mod14.ml = lme(fixed = tumour~sensitivity*line + line:month + 
                 I(month^2), random=~1|patient, data = df, corr=corARMA(form=~1|patient,p=1,q=2), method = "ML")

# ARMA with hetero
mod9.ml = lme(fixed = tumour~sensitivity*line + month + 
                I(month^2), random=~1|patient, data = df, weights=varExp(form=~sensitivity), corr=corAR1(form=~1|patient), method = "ML")
mod10.ml = lme(fixed = tumour~sensitivity*line + line:month + 
                 I(month^2), random=~1|patient, data = df, weights=varExp(form=~sensitivity), corr=corARMA(form=~1|patient,p=2,q=0), method = "ML")
mod11.ml = lme(fixed = tumour~sensitivity*line + line:month + 
                 I(month^2), random=~1|patient, data = df, weights=varExp(form=~sensitivity), corr=corARMA(form=~1|patient,p=2,q=1), method = "ML")
mod12.ml = lme(fixed = tumour~sensitivity*line + line:month + 
                 I(month^2), random=~1|patient, data = df, weights=varExp(form=~sensitivity), corr=corARMA(form=~1|patient,p=2,q=2), method = "ML") # doesn't converge
mod13.ml = lme(fixed = tumour~sensitivity*line + line:month + 
                 I(month^2), random=~1|patient, data = df, weights=varExp(form=~sensitivity), corr=corARMA(form=~1|patient,p=0,q=2), method = "ML")
mod14.ml = lme(fixed = tumour~sensitivity*line + line:month + 
                 I(month^2), random=~1|patient, data = df, corr=corARMA(form=~1|patient,p=1,q=2), method = "ML")


anova(modbest.ml, mod9.ml)
anova(modbest.ml, mod10.ml)
anova(modbest.ml, mod11.ml) # all are significantly better than mod5 but mod11 a little better
anova(modbest.ml, mod13.ml)
anova(modbest.ml, mod14.ml)

Anova(mod13.ml, type = 3)
anova(mod13.ml, type = 'marginal')
summary(mod11.ml)
Anova(mod5.ml.hetero, type = 3)
anova(mod5.ml.hetero, type = 'marginal')
summary(mod5.ml.hetero)


# Modelling for Compsymm

#varExp Sensitivity 
mod5.ml = lme(fixed = tumour~sensitivity*line + month + 
                  I(month^2), random=~1|patient, data = df, method = "ML")
mod5.ml.Cor = lme(fixed = tumour ~ sensitivity*line + month + 
                    I(month^2), random = ~1|patient,correlation = corCompSymm(form = ~ 1|patient),data = df, method = "ML")

summary(mod5.ml.Cor)
summary(mod5.ml)
anova(mod5.ml,mod5.ml.Cor)
Anova(mod5.ml.rob, type = 'III')
# Robust model
df[df$patient == 'O',]
mod5.ml.rob <- rlmer(tumour ~ sensitivity * line + month + I(month^2) + (1|patient), data = df)
summary(mod5.ml.rob)

df$weights <- summary(mod5.ml.rob)$wgt.e

df %>%  mutate(patient = fct_reorder(patient, weights, .fun='min')) %>% 
  ggplot(aes(x = weights, y = patient)) + 
  geom_point()


gg11 = ggplot(df, aes(x=sensitivity, y=resid(mod5.ml.rob, type = "weighted"))) +
  geom_point(fill="slateblue", alpha=0.5) +
  ggtitle("Residuals over sensitivity \n(Robustness modelled)") +
  xlab("Sensitivity") +
  ylab("Weighted residuals") +
  geom_smooth(method=lm, size=0.5, linetype="dashed" , color="red", se=TRUE) +
  theme(plot.title = element_text(size=10))

gg12 = ggplot(df, aes(x=sensitivity, y=resid(mod5.ml, type = "r"))) +
  geom_point(fill="slateblue", alpha=0.5) +
  ggtitle("Residuals over sensitivity") +
  xlab("Sensitivity") +
  ylab("Residuals") +
  geom_smooth(method=lm, size=0.5, linetype="dashed" , color="red", se=TRUE) +
  theme(plot.title = element_text(size=10))

gg13 = ggplot(df, aes(x=line, y=resid(mod5.ml.rob, type = "weighted"))) +
  geom_boxplot(fill="slateblue", alpha=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  ggtitle("Residuals over line \n(Robustness modelled)") +
  xlab("Line") +
  ylab("Weighted residuals") +
  theme(plot.title = element_text(size=10))

gg14 = ggplot(df, aes(x=line, y=resid(mod5.ml, type = "r"))) +
  geom_boxplot(fill="slateblue", alpha=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  ggtitle("Residuals over line") +
  xlab("Line") +
  ylab("Residuals") +
  theme(plot.title = element_text(size=10))

gg15 = ggplot(df, aes(x=as.factor(month), y=resid(mod5.ml.rob, type = "weighted"))) +
  geom_boxplot(fill="slateblue", alpha=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  ggtitle("Residuals over month \n(Robustness modelled)") +
  ylab("Weighted residuals") +
  xlab('Month') +
  theme(plot.title = element_text(size=10))

gg16 = ggplot(df, aes(x=as.factor(month), y=resid(mod5.ml, type = "r"))) +
  geom_boxplot(fill="slateblue", alpha=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  ggtitle("Residuals over month") +
  ylab("Residuals") +
  xlab('Month') +
  theme(plot.title = element_text(size=10))

ggarrange(gg12, gg14, gg16, gg11, gg13, gg15, 
          ncol = 3, nrow = 2)

summary(mod5.ml.rob)
# final model is mod5.ml.hetero

summary(mod5.ml.hetero) # Here interaction between sens. and line 2 become signif. for one-sided test

# Check for normality of residuals
ggplot(df, aes(sample = resid(mod5.ml.hetero, type = "p"))) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("QQ-plot of residuals") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")
# It is the only issue I saw

# check for pattern, heteroskedasticity
ggplot(df, aes(x=fitted(mod5.ml.hetero), y=resid(mod5.ml.hetero, type = "p"))) +
  geom_point(fill="slateblue", alpha=0.5) +
  ggtitle("Residuals over fitted tumour") +
  xlab("Fitted tumour") + 
  ylab("Pearson residuals") +
  geom_smooth(method=lm, size=0.5, linetype="dashed" , color="red", se=TRUE)

ggplot(df, aes(x=fitted(mod5.ml.hetero), y=resid(mod5.ml.hetero, type = "p"))) +
  geom_point(fill="slateblue", alpha=0.5) +
  ggtitle("Residuals over line") +
  xlab("Fitted tumour") +
  ylab("Pearson residuals") +
  facet_wrap(line~., ncol = 2) +
  geom_smooth(method=lm, size=0.5, linetype="dashed" , color="red", se=TRUE)

ggplot(df, aes(x=fitted(mod5.ml.hetero), y=resid(mod5.ml.hetero, type = "p"))) +
  geom_point(fill="slateblue", alpha=0.5) +
  ggtitle("Residuals over month") +
  xlab("Fitted tumour") +
  ylab("Pearson residuals") +
  facet_wrap(month~., ncol = 2) +
  geom_smooth(method=lm, size=0.5, linetype="dashed" , color="red", se=TRUE)

ggplot(df, aes(x=fitted(mod5.ml.hetero), y=resid(mod5.ml.hetero, type = "p"))) +
  geom_point(fill="slateblue", alpha=0.5) +
  ggtitle("Residuals over patient") +
  xlab("Fitted tumour") +
  ylab("Pearson residuals") +
  facet_wrap(patient~.) +
  geom_smooth(method=lm, size=0.5, linetype="dashed" , color="red", se=TRUE)

# Residuals check for pattern
ggplot(df, aes(x=line, y=resid(mod5.ml.hetero, type = "p"))) +
  geom_boxplot(fill="slateblue", alpha=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  ggtitle("Residuals over line") +
  ylab("Pearson residuals")

ggplot(df, aes(x=as.factor(month), y=resid(mod5.ml.hetero, type = "p"))) +
  geom_boxplot(fill="slateblue", alpha=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  xlab("Month") +
  ggtitle("Residuals over month") +
  ylab("Pearson residuals")

ggplot(df, aes(x=patient, y=resid(mod5.ml.hetero, type = "p"))) +
  geom_boxplot(fill="slateblue", alpha=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  ggtitle("Residuals over patient") +
  ylab("Pearson residuals")

# Check for normality of random effect
ggplot(ranef(mod5.ml.hetero), aes(sample = ranef(mod5.ml.hetero)[,1])) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("QQ-plot random effect") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")

# Check for time correlation
plot(nlme:::ACF.lme(mod9.ml,type="pearson"),alpha=0.05) 
plot(nlme:::ACF.lme(mod5.ml.hetero,type="pearson", maxLag = 13),alpha=0.05)

round(r.squaredGLMM(mod5.ml.hetero),4)

# Compare with other model
summary(mod1.ml)
summary(mod11.ml)
summary(mod5.ml.hetero)

mod5.ml = lme(fixed = tumour~sensitivity*line + month + 
                I(month^2), random=~1|patient, data = df, method = "ML")
mod5.ml.hetero = lme(fixed = tumour~sensitivity*line + month + 
                I(month^2), random=~1|patient, data = df, weights=varExp(form=~sensitivity), method = "ML")

round(r.squaredGLMM(mod5.ml),4)
round(r.squaredGLMM(mod5.ml.hetero),4)
round(r.squaredGLMM(mod9.ml),4)

## Test for the interaction between line and sensitivity

mod5bis.ml = lme(fixed = tumour~sensitivity+line + month + 
                   I(month^2), random=~1|patient, data = df, method = "ML")
mod11bis.ml = lme(fixed = tumour~sensitivity+line + month + 
                 I(month^2), random=~1|patient, data = df, corr=corARMA(form=~1|patient,p=2,q=1), method = "ML")

Anova(mod11.ml, type = 3)
anova(mod11bis.ml, mod11.ml)
Anova(mod5.ml, type = 3)
anova(mod11bis.ml, mod5bis.ml)

## has to be forgotten ################################################################################################################
fit.restr.real <- lme(fixed = tumour~sensitivity+line + month + 
                                    I(month^2), random=~1|patient, data = df, corr=corARMA(form=~1|patient,p=2,q=1), method = "ML")
fit.full.real <- lme(fixed = tumour~sensitivity*line + month + 
                                   I(month^2), random=~1|patient, data = df, corr=corARMA(form=~1|patient,p=2,q=1), method = "ML")
############################################################################################################################################

mod5.ml.restr<-lmer(tumour~sensitivity+line + month + I(month^2) + 
                      (1|patient), data = df, REML = F)
mod5.ml = lmer(tumour~sensitivity*line + month + I(month^2) + 
                     (1|patient), data = df, REML = F)

n.r = 1000
mx.t.nr = simulate(mod5.ml.restr, nsim=n.r)
pval.r = Chisq.r = lrt.r = rep(NA,n.r)
data.r <- df

for (rw in 1:n.r){
  data.r$tumour = mx.t.nr[,rw]
  fit.restr = lmer(tumour~sensitivity+line + month + I(month^2) + (1|patient), data = data.r, REML = FALSE)
  fit.full = lmer(tumour~sensitivity*line + month + I(month^2) + (1|patient), data = data.r, REML = FALSE)
  lrt.r[rw] = as.numeric(2*(logLik(fit.full) - logLik(fit.restr)))
  anova.r = anova(fit.restr, fit.full)
  pval.r[rw] = anova.r$Pr[2]
  Chisq.r[rw] = anova.r$Chisq[2]
}
summary(mod11.ml)
anova(mod5.ml.full)

mean(pval.r<0.05)

qbinom(c(0.025,0.975), prob = 0.05, size = 1000)/1000

LRT.real = as.numeric(2*(logLik(mod5.ml.full) - logLik(mod5.ml.restr)))
mean(LRT.real<lrt.r) # This is the p-value of LRT corrected with Bootstrap

pbkrtest::KRmodcomp(mod5.ml.restr,mod5.ml.full)

Chisq.r.dens2 = as.data.frame(Chisq.r)

## chi-squared comparaison 1.0

chi_qq = data.frame(theor = qchisq(ppoints(1000), df = 2), sample = Chisq.r.dens2$Chisq.r[order(Chisq.r.dens2$Chisq.r)])

qq_chi = ggplot(chi_qq, aes(y = sample, x = theor)) +
  geom_point()+
  geom_abline(slope = 1, color = 'red') +
  ggtitle(expression('QQ-plot for'~chi[2]^2)) +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")

## chi-squared comparaison 2.0

simul.chi = dchisq(ppoints(1000), 2)*1000
xx = c(Chisq.r, simul.chi)
yy = c(rep('A', 1000), rep('B', 1000))
Chisq.r.dens = data_frame(yy = yy, xx = xx)

dchisq(x, df = 10)

x = data.frame(x = Chisq.r)
ggplot(Chisq.r.dens, aes(x=xx, fill = yy)) +
  geom_histogram(alpha=0.7, binwidth=0.2, position = "identity") + 
  geom_vline(xintercept = LRT.real, linetype="dashed", color = "red", size=0.5) +
  xlab('LRT') +
  ggtitle("Simulated vs real distribution") + 
  scale_fill_manual(name="Chi-squared",values=c("#eb0087","#554696"),labels=c("Simulated","Real"))

ggplot(x, aes(x=x)) +
  geom_histogram(aes(y = ..density.., fill = "Bootstrap"), alpha=0.7, binwidth=0.2, position = "identity") + 
  xlab('LRT') +
  ggtitle(expression(chi[2]^2)) + 
  scale_fill_manual(name="Chi-squared",values=c("#eb0087"),labels=c("Simulated")) + 
  stat_function(fun = dchisq, args = list(df = 2), aes(colour = "Chi-squared")) +
  scale_colour_manual("", values = c("red")) +
  scale_fill_manual("", values = c("darkgrey"))



gg_simul1 = ggplot(Chisq.r.dens, aes(x=xx, fill = yy)) +
  geom_histogram(alpha=0.7, binwidth=0.2, position = "identity") +
  xlab('LRT') +
  ggtitle(expression("Simulated vs", chi_{2})) + 
  scale_fill_manual(name="Chi-squared",values=c("#eb0087","#554696"),labels=c("Simulated","Real"))

gg_simul1 = ggplot(x, aes(x=x)) +
  geom_histogram(aes(y = ..density.., fill = "Bootstrap LRT"), alpha=0.5, binwidth=0.2, position = "identity") + 
  xlab('LRT') +
  ggtitle(expression('Bootstrap LRT vs'~ chi[2]^2~'distribution')) + 
  stat_function(fun = dchisq, args = list(df = 2), aes(colour = 'Chi-squared')) +
  scale_colour_manual("", values = c("red"), labels = expression(chi[2]^2~'distribution')) +
  scale_fill_manual("", values = c("black"))

gg_simul2 = ggplot(Chisq.r.dens2, aes(x=Chisq.r)) +
  geom_histogram(alpha=0.7, binwidth=0.2) + 
  geom_vline(aes(xintercept = LRT.real, color = "red"), linetype="dashed", size=0.5) +
  xlab('LRT') +
  ylim(c(0,100)) +
  ggtitle("LRT corrected with Bootstrap") +
  scale_colour_manual("", values = c("red"), labels = c('LRT'))


count(x)
x <- Chisq.r
x <- rchisq(1000,df = 2)
df <- data.frame(x = x)

ggplot(df, aes(x)) +
  geom_histogram(binwidth = 0.1,
                 alpha = 0.3, position = "identity") + 
  stat_theodensity(aes(y = stat(count) * 0.1, colour = "Chi_square"),
                   distri = "chisq", geom = "line") +
  coord_cartesian(xlim = c(0, 15))

ggarrange(gg_simul1, qq_chi, gg_simul2, ncol = 2, nrow = 2)

gg_simul1 = ggplot(x, aes(x=x)) +
  geom_histogram(aes(y = ..count.., fill = "Bootstrap LRT"), alpha=0.5, binwidth=0.2, position = "identity") + 
  xlab('LRT') +
  ggtitle(expression('Bootstrap LRT vs'~ chi[2]^2~'distribution')) + 
  stat_function(fun = function(x, df, n){n * dchisq(x = x, df = df)}, args = with(x, c(df = 2, n = 200)), aes(colour = 'Chi-squared')) +
  scale_colour_manual("", values = c("red"), labels = expression(chi[2]^2~'distribution')) +
  scale_fill_manual("", values = c("black")) +
  theme(legend.title=element_blank())


dchisq()
summary(mod5.ml.hetero)
# QQ.plot with transformed data

df1 = df
df2 = df
df1$tumour = exp(df$tumour)
df2$tumour = log(df$tumour+4)

qq12 = ggplot(df1, aes(sample = tumour)) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles") +
  facet_wrap(line~., ncol = 2)


qq13 = ggplot(df2, aes(sample = tumour)) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles") +
  facet_wrap(line~., ncol = 2)

qq14 = ggplot(df, aes(sample = tumour)) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("QQ-plot for Tumour") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")

qq15 = ggplot(df1, aes(sample = tumour)) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("QQ-plot exp(Tumour)") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")

qq16 = ggplot(df2, aes(sample = tumour)) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("QQ-plot log(Tumour)") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")

qq17 = ggplot(df[df$line == 1,], aes(sample = tumour)) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("QQ-plot Tumour line 1") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")

qq18 = ggplot(df[df$line == 2,], aes(sample = tumour)) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("QQ-plot Tumour line 2") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")

qq19 = ggplot(df[df$line == 3,], aes(sample = tumour)) +
  stat_qq() + 
  stat_qq_line(color = 'red', linetype="dashed") +
  ggtitle("QQ-plot Tumour line 3") +
  xlab("Theoritical Quantiles") +
  ylab("Sample Quantiles")

qq_all2 = ggarrange(qq14, qq15, qq16, qq17, qq18, qq19,
              ncol = 3, nrow = 2)
annotate_figure(qq_all2,top = text_grob("Simulated residuals versus fitted residuals", color = "black", face = "bold", size = 15))

mod5.ml.rob <- rlmer(tumour ~ sensitivity * line + month + I(month^2) + (1|patient), data = df)
mod5.ml.rob.bis <- rlmer(tumour ~ sensitivity + line + month + I(month^2) + (1|patient), data = df)

pbkrtest::KRmodcomp(mod5.ml,mod5.ml.hetero)
anova(mod5.ml,mod5.ml.hetero)
summary(mod5.ml.rob)
resid(mod5.ml.rob)
summary(mod9.ml)



