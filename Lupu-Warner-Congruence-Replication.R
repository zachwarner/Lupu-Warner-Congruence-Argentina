#########################################################
# Mass-Elite Congruence and Representation in Argentina
# Noam Lupu and Zach Warner
# Revised 11/2015
#########################################################

rm(list=ls()); gc()
require(foreign); require(ggplot2); require(arm); require(scales); require(MCMCglmm); require(boot)
require(plotMCMC); require(grid)
setwd("/mywd")


##### VERSION CONTROL #####
sessionInfo()
# R version 3.2.1 (2015-06-18)
# Platform: x86_64-apple-darwin13.4.0 (64-bit)
# Running under: OS X 10.10.5 (Yosemite)
# 
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
# [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] plotMCMC_2.0-0 MCMCglmm_2.22  ape_3.3        coda_0.18-1    scales_0.3.0   arm_1.8-6      lme4_1.1-10   
# [8] Matrix_1.2-2   MASS_7.3-44    ggplot2_1.0.1  foreign_0.8-66
# 
# loaded via a namespace (and not attached):
# [1] Rcpp_0.12.1        magrittr_1.5       splines_3.2.1      munsell_0.4.2      cubature_1.1-2    
# [6] colorspace_1.2-6   lattice_0.20-33    minqa_1.2.4        stringr_1.0.0      plyr_1.8.3        
# [11] caTools_1.17.1     tools_3.2.1        gtable_0.1.2       nlme_3.1-122       KernSmooth_2.23-15
# [16] corpcor_1.6.8      gtools_3.5.0       abind_1.4-3        digest_0.6.8       tensorA_0.36      
# [21] nloptr_1.0.4       reshape2_1.4.1     bitops_1.0-6       gdata_2.17.0       stringi_1.0-1     
# [26] gplots_2.17.0      proto_0.3-10    


##### ANALYSIS #####
### Some tools
# to rescale variables
rescalr <- function(x){
  x <- (x-min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T))
}
# scales for grid plotting
integer_breaks <- function(n = 3, ...) {
  breaker <- pretty_breaks(n, ...)
  function(x) {
    breaks <- breaker(x)
    breaks[breaks == breaks]
  }
}
# for R^2 from glmer/lmer model
r2.corr.mer <- function(m) {
  lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
  summary(lmfit)$r.squared
} 
# for subsetting dataframes
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

### Figure 1: Explaining the overlap measure
set.seed(1234567)
ARGcitizens <- read.csv("Lupu-Warner-ARGcitizens.csv")
ARGelites <- read.csv("Lupu-Warner-ARGelites.csv")
dens <- data.frame(dem_agg=c(ARGcitizens$P41,ARGelites$P67),type=c(rep("cit",1200),rep("el",140)))
fig1 <- ggplot(dens, aes(x=dem_agg,group=type)) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  geom_density(fill="gray20",alpha=.4, adjust=2) +
  annotate("text",x=7.5,y=.15,label="Citizens",family="Times") +
  annotate("text",x=1.5,y=.15,label="Elites",family="Times") +
  labs(title=NULL,x="Self-placement on a left-right scale",y="Density") +
  theme(text=element_text(family="Times"))
pdf("figure1.pdf",width=5,height=6)
fig1 # ignore the warning message, it's for NA values
dev.off()
rm(dens,fig1)


### Table 1: Mass-elite congruence
df <- read.csv("Lupu-Warner-PDFoverlap.csv")
tab <- data.frame(matrix(nrow=8,ncol=0))
tab$name <- c("Ideology","Democracy","Economic policy","Ideal society","Populism",
              "Order versus liberty","Decentralization","Most important problem")
tab$question <- as.character(unique(df$question)[c(14,16,29,1,32,19,25,10)])
tab$c.all_e.job.0 <- tab$c.all_e.job.1 <- tab$c.class.4_e.all <- tab$c.class.1_e.all <- tab$c.all_e.all <- rep(NA,nrow(tab)) 
for(i in 1:nrow(tab)){
  tab$c.all_e.all[i] <- df$congruence[which(df$elite.group == "e.party.all.job.all" & 
                                            df$citizen.group == "c.party.all.c.class.all.c.buenos.all" &
                                            df$question == tab$question[i])]
  tab$c.all_e.job.0[i] <- df$congruence[which(df$elite.group == "e.party.all.job.0" & 
                                              df$citizen.group == "c.party.all.c.class.all.c.buenos.all" &
                                              df$question == tab$question[i])]
  tab$c.all_e.job.1[i] <- df$congruence[which(df$elite.group == "e.party.all.job.1" & 
                                                      df$citizen.group == "c.party.all.c.class.all.c.buenos.all" &
                                                      df$question == tab$question[i])]
  tab$c.class.1_e.all[i] <- df$congruence[which(df$elite.group == "e.party.all.job.all" & 
                                                  df$citizen.group == "c.party.all.c.class.1.c.buenos.all" &
                                                  df$question == tab$question[i])]
  tab$c.class.4_e.all[i] <- df$congruence[which(df$elite.group == "e.party.all.job.all" & 
                                                df$citizen.group == "c.party.all.c.class.4.c.buenos.all" &
                                                df$question == tab$question[i])]
}
tab
rm(df,tab,i)


### DYADIC ANALYSIS
df <- read.csv("Lupu-Warner-Dyads.csv")

# rescale
df$cong_econ_state_rescale <- df$cong_econ_state - min(df$cong_econ_state, na.rm=T) # flip it to positive
df$cong_econ_state_rescale <- rescalr(df$cong_econ_state_rescale)
df$cong_pop_rescale <- df$cong_pop - min(df$cong_pop, na.rm=T) # flip it to positive
df$cong_pop_rescale <- rescalr(df$cong_pop_rescale)
df$cong_P41_rescale <- df$cong_P41 - min(df$cong_P41, na.rm=T) # flip it to positive
df$cong_P41_rescale <- rescalr(df$cong_P41_rescale)


### Figure 2/Table A1: Results from dyadic analysis
# estimate the models
mod1.P41 <- lmer(cong_P41_rescale ~ as.factor(c_buenosaires) + c_age + as.factor(c_female) + as.factor(c_P37) + 
                   as.factor(c_NSE_AGR)+ as.factor(e_job) + as.factor(e_Partido) + (1|c_ncue) + (1|e_n), data=df)
summary(mod1.P41)
r2.corr.mer(mod1.P41)

mod1.P50 <- glmer(cong_P50  ~ as.factor(c_buenosaires) + c_age + as.factor(c_female) + as.factor(c_P37) + 
                   as.factor(c_NSE_AGR)+ as.factor(e_job) + as.factor(e_Partido) + (1|c_ncue) + (1|e_n), 
                 family=binomial(link="logit"), data=df)
relgrad <- with(mod1.P50@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad)) # This is large so we should be careful with inferences
summary(mod1.P50)
r2.corr.mer(mod1.P50)

mod1.econ.state <- lmer(cong_econ_state_rescale ~ as.factor(c_buenosaires) + c_age + as.factor(c_female) + as.factor(c_P37) + 
                          as.factor(c_NSE_AGR)+ as.factor(e_job) + as.factor(e_Partido) + (1|c_ncue) + (1|e_n), data=df)
summary(mod1.econ.state)
r2.corr.mer(mod1.econ.state)

mod1.pop <- lmer(cong_pop_rescale ~ as.factor(c_buenosaires) + c_age + as.factor(c_female) + as.factor(c_P37) + 
                  as.factor(c_NSE_AGR)+ as.factor(e_job) + as.factor(e_Partido) + (1|c_ncue) + (1|e_n), data=df)
summary(mod1.pop)
r2.corr.mer(mod1.pop)

mod1.P55 <- glmer(cong_P55 ~ as.factor(c_buenosaires) + c_age + as.factor(c_female) + as.factor(c_P37) + 
                  as.factor(c_NSE_AGR)+ as.factor(e_job) + as.factor(e_Partido) + (1|c_ncue) + (1|e_n), 
                  family=binomial(link="logit"), data=df)
# Did it fail to converge? Well, on an absolute convergence criterion. Let's try a relative convergence criterion.
relgrad <- with(mod1.P55@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad)) # This is smaller than .001, the cutoff for the absolute criterion, so it's probably all right.
summary(mod1.P55)
r2.corr.mer(mod1.P55)

# Collect the coefficients and standard errors. Insert 0's for the baseline conditions in the plots
coef.mod1.econ.state <- c(0,fixef(mod1.econ.state)[2],0,fixef(mod1.econ.state)[c(5:6,8)],0,
                          fixef(mod1.econ.state)[c(9:11)],0,fixef(mod1.econ.state)[12],0,
                          fixef(mod1.econ.state)[c(13:15)])
se.mod1.econ.state <- c(0,se.fixef(mod1.econ.state)[2],0,se.fixef(mod1.econ.state)[c(5:6,8)],0,
                        se.fixef(mod1.econ.state)[c(9:11)],0,se.fixef(mod1.econ.state)[12],0,
                        se.fixef(mod1.econ.state)[c(13:15)])
coef.mod1.P50 <- c(0,fixef(mod1.P50)[2],0,fixef(mod1.P50)[c(5:6,8)],0,fixef(mod1.P50)[c(9:11)],0,
                   fixef(mod1.P50)[12],0,fixef(mod1.P50)[c(13:15)])
se.mod1.P50   <- c(0,se.fixef(mod1.P50)[2],0,se.fixef(mod1.P50)[c(5:6,8)],0,se.fixef(mod1.P50)[c(9:11)],0,
                   se.fixef(mod1.P50)[12],0,se.fixef(mod1.P50)[c(13:15)])
coef.mod1.pop <- c(0,fixef(mod1.pop)[2],0,fixef(mod1.pop)[c(5:6,8)],0,fixef(mod1.pop)[c(9:11)],0,
                   fixef(mod1.pop)[12],0,fixef(mod1.pop)[c(13:15)])
se.mod1.pop   <- c(0,se.fixef(mod1.pop)[2],0,se.fixef(mod1.pop)[c(5:6,8)],0,se.fixef(mod1.pop)[c(9:11)],0,
                   se.fixef(mod1.pop)[12],0,se.fixef(mod1.pop)[c(13:15)])
coef.mod1.P41 <- c(0,fixef(mod1.P41)[2],0,fixef(mod1.P41)[c(5:6,8)],0,fixef(mod1.P41)[c(9:11)],0,
                   fixef(mod1.P41)[12],0,fixef(mod1.P41)[c(13:15)])
se.mod1.P41   <- c(0,se.fixef(mod1.P41)[2],0,se.fixef(mod1.P41)[c(5:6,8)],0,se.fixef(mod1.P41)[c(9:11)],0,
                   se.fixef(mod1.P41)[12],0,se.fixef(mod1.P41)[c(13:15)])
coef.mod1.P55 <- c(0,fixef(mod1.P55)[2],0,fixef(mod1.P55)[c(5:6,8)],0,fixef(mod1.P55)[c(9:11)],0,
                   fixef(mod1.P55)[12],0,fixef(mod1.P55)[c(13:15)])
se.mod1.P55   <- c(0,se.fixef(mod1.P55)[2],0,se.fixef(mod1.P55)[c(5:6,8)],0,se.fixef(mod1.P55)[c(9:11)],0,
                 se.fixef(mod1.P55)[12],0,se.fixef(mod1.P55)[c(13:15)])

# Get ready to plot
plotmat <- data.frame(matrix(ncol=6,nrow=80))
colnames(plotmat) <- c("lo","est","hi","iter","issue","cov")
plotmat$est <- c(coef.mod1.econ.state,coef.mod1.P50,coef.mod1.pop,coef.mod1.P41,coef.mod1.P55)
plotmat$lo <- c((coef.mod1.econ.state - (qnorm(.975)*se.mod1.econ.state)),
                (coef.mod1.P50 - (qnorm(.975)*se.mod1.P50)),
                (coef.mod1.pop - (qnorm(.975)*se.mod1.pop)),
                (coef.mod1.P41 - (qnorm(.975)*se.mod1.P41)),
                (coef.mod1.P55 - (qnorm(.975)*se.mod1.P55)))
plotmat$hi <- c((coef.mod1.econ.state + (qnorm(.975)*se.mod1.econ.state)),
                (coef.mod1.P50 + (qnorm(.975)*se.mod1.P50)),
                (coef.mod1.pop + (qnorm(.975)*se.mod1.pop)),
                (coef.mod1.P41 + (qnorm(.975)*se.mod1.P41)),
                (coef.mod1.P55 + (qnorm(.975)*se.mod1.P55)))
plotmat$iter <- rev(rep(c(1:4,6:7,9:12,14:17,19:20),5))
plotmat$issue <- c(rep("Economic policy",16),rep("Democracy",16),rep("Populism",16),
                   rep("Ideology",16),rep("Order versus liberty",16))
plotmat$cov <- ifelse(plotmat$lo < 0 & plotmat$hi > 0,1,2)
plotmat$cov[which(plotmat$hi == plotmat$lo & plotmat$lo == plotmat$est)] <- 1
vars <- rev(c("GBA non-resident","GBA resident","Peronist (FPV)","Dissident Peronist",
              "Non-Peronist opposition","No party", "SES level 1 (A, B, C1)",
              "SES level 2 (C2)","SES level 3 (C3)","SES level 4 (D1, D2, E)","Elite legislative",
              "Elite executive","Elite Peronist (FPV)","Elite dissident Peronist",
              "Elite non-Peronist opposition","Elite regional party"))

# Add some fake data to get the x-axis into shape. This data won't be plotted; ggplot will throw a warning.
plotmat[c(81:85),2] <- rep(0,5)
plotmat[c(81:85),4] <- rep(-1,5)
plotmat[c(81:85),5] <- c("Economic policy","Populism","Ideology","Democracy","Order versus liberty")
plotmat[c(81:85),6] <- rep(1,5)
fakemin <- -.252
fakemax <- .252
plotmat[c(81:83),1] <- rep(fakemin,3)
plotmat[c(81:83),3] <- rep(fakemax,3)
fakemin <- -2.6770
fakemax <- 2.6770
plotmat[c(84:85),1] <- rep(fakemin,2)
plotmat[c(84:85),3] <- rep(fakemax,2)
plotmat$issue <- factor(plotmat$issue, levels=c("Ideology","Democracy","Economic policy",
                                                "Populism","Order versus liberty"))

# plot
fig2 <- ggplot(plotmat, aes(x=est,y=iter)) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  facet_grid(.~issue, scales="free_x") +
  scale_x_continuous(breaks = integer_breaks()) +
  geom_vline(xintercept=0,linetype="longdash",colour="darkgray") +
  geom_segment(data=plotmat,aes(x=lo,xend=hi,y=iter,yend=iter,colour=factor(cov)),size=.5,show_guide=F) +
  geom_point(show_guide=F,shape=21,aes(x=est,y=iter,fill=factor(cov)),size=2) +
  scale_fill_manual(values=c("white", "black")) +
  scale_colour_manual(values=c("black", "black")) +
  scale_y_continuous(limits = c(1,20), breaks = c(1:4,6:7,9:12,14:17,19:20),labels=vars) +
  labs(title=NULL,x="Difference in elite-mass congruence",y=NULL) +
  theme(text=element_text(family="Times"))
pdf("figure2.pdf", width=11, height=5)
fig2
dev.off() # warnings are the 5 fake data points not plotted, only used for x axis.

### Effect sizes/substantive significance
# GBA residency: Economic policy
testdat.GBA   <- c(1,1,round(44.51),1,0,0,0,1,0,0,1,0,0,1,0)
testdat.NoGBA <- c(1,0,round(44.51),1,0,0,0,1,0,0,1,0,0,1,0)
testdat <- cbind(testdat.GBA,testdat.NoGBA)
preds <- t(matrix(fixef(mod1.econ.state))) %*% testdat # predicted congruence for econ
(preds[1]-preds[2])/preds[2]
# Party ID: Ideology, Democracy, and Economic policy
testdat.P37_2 <- c(1,0,round(44.51),1,1,0,0,0,0,0,1,0,0,1,0)
testdat.P37_1 <- c(1,0,round(44.51),1,0,0,0,0,0,0,1,0,0,1,0)
testdat <- cbind(testdat.P37_2,testdat.P37_1)
preds <- t(matrix(fixef(mod1.P41))) %*% testdat # predicted congruence for econ
(preds[2]-preds[1])/preds[1]
preds <- inv.logit(t(matrix(fixef(mod1.P50))) %*% testdat) # predicted congruence for dem
(preds[2]-preds[1])/preds[1]
preds <- t(matrix(fixef(mod1.econ.state))) %*% testdat # predicted congruence for econ
(preds[2]-preds[1])/preds[1]
# SES: Populism and Order versus liberty
testdat.NSE4 <- c(1,0,round(44.51),1,0,0,0,1,0,0,1,0,0,1,0)
testdat.NSE1 <- c(1,0,round(44.51),1,0,0,0,1,0,0,0,0,0,1,0)
testdat <- cbind(testdat.NSE1,testdat.NSE4)
preds <- t(matrix(fixef(mod1.pop))) %*% testdat # predicted congruence for populism
(preds[1]-preds[2])/preds[1]
preds <- inv.logit(t(matrix(fixef(mod1.P55))) %*% testdat) # predicted congruence for order v liberty
preds
preds[1]/preds[2]

### Figure 3/Table A2: Differences by SES
# rescale variables to be on the unit interval
ARGcitizens$P50_rescale <- rescalr(ARGcitizens$P50)
ARGcitizens$econ_state_rescale <- rescalr(ARGcitizens$econ_state_agg)
ARGcitizens$pop_rescale <- rescalr(ARGcitizens$populist_agg)
ARGcitizens$P41_rescale <- rescalr(ARGcitizens$P41)
ARGcitizens$P55_2[which(ARGcitizens$P55 == 2)] <- 0
ARGcitizens$P55_2[which(ARGcitizens$P55 == 1)] <- 1

# estimate
options(na.action = na.omit)
mod2.P41 <-  glm(P41_rescale ~ as.factor(NSE_AGR), data=ARGcitizens)
summary(mod2.P41)
r2.corr.mer(mod2.P41)
mod2.P50 <- glm(P50_rescale ~ as.factor(NSE_AGR), data=ARGcitizens)
summary(mod2.P50)
r2.corr.mer(mod2.P50)
mod2.econ.state <- glm(econ_state_rescale ~ as.factor(NSE_AGR), data=ARGcitizens)
summary(mod2.econ.state)
r2.corr.mer(mod2.econ.state)
mod2.pop <-  glm(pop_rescale ~ as.factor(NSE_AGR), data=ARGcitizens)
summary(mod2.pop)
r2.corr.mer(mod2.pop)
#linear probability model for comparability of coefficients; results robust to using a logit (see below)
mod2.P55 <- glm(P55_2 ~ as.factor(NSE_AGR), data=ARGcitizens) 
summary(mod2.P55)
r2.corr.mer(mod2.P55)

# collect coefficients and SEs
coef.mod2.P50 <- c(0,coef(mod2.P50)[2:4])
se.mod2.P50 <- c(0,sqrt(diag(vcov(mod2.P50)))[2:4])
coef.mod2.econ.state <- c(0,coef(mod2.econ.state)[2:4])
se.mod2.econ.state <- c(0,sqrt(diag(vcov(mod2.econ.state)))[2:4])
coef.mod2.P41 <- c(0,coef(mod2.P41)[2:4])
se.mod2.P41 <- c(0,sqrt(diag(vcov(mod2.P41)))[2:4])
coef.mod2.P55 <- c(0,coef(mod2.P55)[2:4])
se.mod2.P55 <- c(0,sqrt(diag(vcov(mod2.P55)))[2:4])
coef.mod2.pop <- c(0,coef(mod2.pop)[2:4])
se.mod2.pop <- c(0,sqrt(diag(vcov(mod2.pop)))[2:4])

# get ready to plot
plotmat <- data.frame(matrix(ncol=6,nrow=20))
colnames(plotmat) <- c("lo","est","hi","iter","issue","cov")
plotmat$est <- c(coef.mod2.P50,coef.mod2.econ.state,coef.mod2.P41,coef.mod2.P55,coef.mod2.pop)
plotmat$lo <- c((coef.mod2.P50 - (qnorm(.975)*se.mod2.P50)),
                (coef.mod2.econ.state - (qnorm(.975)*se.mod2.econ.state)),
                (coef.mod2.P41 - (qnorm(.975)*se.mod2.P41)),
                (coef.mod2.P55 - (qnorm(.975)*se.mod2.P55)),
                (coef.mod2.pop - (qnorm(.975)*se.mod2.pop)))
plotmat$hi <- c((coef.mod2.P50 + (qnorm(.975)*se.mod2.P50)),
                (coef.mod2.econ.state + (qnorm(.975)*se.mod2.econ.state)),
                (coef.mod2.P41 + (qnorm(.975)*se.mod2.P41)),
                (coef.mod2.P55 + (qnorm(.975)*se.mod2.P55)),
                (coef.mod2.pop + (qnorm(.975)*se.mod2.pop)))
plotmat$iter <- rev(rep(c(1:4),5))
plotmat$issue <- c(rep("Democracy",4),rep("Economic policy",4),rep("Ideology",4),
                   rep("Order versus liberty",4),rep("Populism",4))
plotmat$cov <- ifelse(plotmat$lo < 0 & plotmat$hi > 0,1,2)
plotmat$cov[which(plotmat$hi == plotmat$lo & plotmat$lo == plotmat$est)] <- 1
vars <- rev(c("SES level 1 (A, B, C1)","SES level 2 (C2)","SES level 3 (C3)","SES level 4 (D1, D2, E)"))

# restrict attention to SES level 4
plotmat2 <- plotmat[which(plotmat$iter == 1),]
plotmat2 <- plotmat2[c(3,1,2,5,4),]
plotmat2$iter <- rev(1:5)
fig3 <- ggplot(plotmat2, aes(x=est,y=iter)) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  geom_vline(xintercept=0,linetype="longdash",colour="darkgray") +
  geom_segment(data=plotmat2,aes(x=lo,xend=hi,y=iter,yend=iter,colour=factor(cov)),size=.5,show_guide=F) +
  geom_point(show_guide=F,shape=21,aes(x=est,y=iter,fill=factor(cov)),size=2) +
  scale_fill_manual(values=c("white", "black")) +
  scale_colour_manual(values=c("black", "black")) +
  scale_y_continuous(breaks = seq(1:5),labels=rev(plotmat2$issue)) +
  labs(title=NULL,x="Divergence between high and low SES respondents",y=NULL) +
  theme(text=element_text(family="Times"))
pdf("figure3.pdf", width=5, height=5)
fig3
dev.off()


rm(list=ls(pattern=c("mod","coef","se","testdat"))); rm(fakemin,fakemax,preds,plotmat,plotmat2,fig2,fig3,relgrad,vars); gc() # ignore warning


### Table A3: Robustness to using only legislative elites for dyadic results ####
leg.df <- df[which(df$e_job == 0),]

mod3.P41 <- lmer(cong_P41_rescale ~ as.factor(c_buenosaires) + c_age + as.factor(c_female) + as.factor(c_P37) + 
                   as.factor(c_NSE_AGR)+ as.factor(e_Partido) + (1|c_ncue) + (1|e_n), data=leg.df)
summary(mod3.P41)
r2.corr.mer(mod3.P41)

# With citizen random effects, it does not evaluate
# mod3.P50 <- glmer(cong_P50 ~ as.factor(c_buenosaires) + c_age + as.factor(c_female) + as.factor(c_P37) + 
#                     as.factor(c_NSE_AGR)+ as.factor(e_Partido) + (1|c_ncue) + (1|e_n),
#                   family=binomial(link="logit"), data=leg.df)
# relgrad <- with(mod3.P50@optinfo$derivs,solve(Hessian,gradient))
# max(abs(relgrad)) # This is untolerably large
# summary(mod1.P50) # can't evaluate
# r2.corr.mer(mod1.P50) # can't evaluate
mod3.P50 <- glmer(cong_P50 ~ as.factor(c_buenosaires) + c_age + as.factor(c_female) + as.factor(c_P37) + 
                                      as.factor(c_NSE_AGR)+ as.factor(e_Partido) + (1|e_n),
                                      family=binomial(link="logit"), data=leg.df) # note no citizen REs
relgrad <- with(mod3.P50@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad)) # This is less than .001 so probably fine
summary(mod3.P50)
r2.corr.mer(mod3.P50)

mod3.econ.state <- lmer(cong_econ_state_rescale ~ as.factor(c_buenosaires) + c_age + as.factor(c_female) + 
                          as.factor(c_P37) + as.factor(c_NSE_AGR)+ as.factor(e_Partido) + (1|c_ncue) + 
                          (1|e_n), data=leg.df)
summary(mod3.econ.state)
r2.corr.mer(mod3.econ.state)

mod3.pop <- lmer(cong_pop_rescale ~ as.factor(c_buenosaires) + c_age + as.factor(c_female) + as.factor(c_P37) + 
                   as.factor(c_NSE_AGR)+ as.factor(e_Partido) + (1|c_ncue) + (1|e_n), data=leg.df)
summary(mod3.pop)
r2.corr.mer(mod3.pop)

mod3.P55 <- glmer(cong_P55 ~ as.factor(c_buenosaires) + c_age + as.factor(c_female) + as.factor(c_P37) + 
                    as.factor(c_NSE_AGR)+ as.factor(e_Partido) + (1|c_ncue) + (1|e_n),
                  family=binomial(link="logit"), data=leg.df)
# Again, convergence fails on absolute criterion, but not on relative criterion.
relgrad <- with(mod3.P55@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad)) # smaller than .001, so it's all right
summary(mod3.P55)
r2.corr.mer(mod3.P55)


### Table A4: Mass-elite congruence on economic policy by elite background
# Rescaling to make the tables more legible
df$e_parents_edu_rescale <- df$e_parents_edu/9 # original scale is 0-9 for each parent
df$e_gp_edu_rescale <- df$e_gp_edu/9 # same original scale

options(na.action = na.omit)
mod3.econ.state.1 <- lmer(cong_econ_state_rescale ~ as.factor(c_buenosaires) + c_age + as.factor(c_female) + 
                            as.factor(c_P37) + as.factor(c_NSE_AGR)+ as.factor(e_job) + as.factor(e_Partido) + 
                            as.numeric(e_parents_edu_rescale)  + (1|c_ncue) + (1|e_n), data=df)
summary(mod3.econ.state.1)
r2.corr.mer(mod3.econ.state.1)

mod3.econ.state.2 <- lmer(cong_econ_state_rescale ~ c_age + as.factor(c_female) + as.numeric(e_parents_edu_rescale) + 
                            (1|c_ncue) + (1|e_n), data=df)
summary(mod3.econ.state.2)
r2.corr.mer(mod3.econ.state.2)

mod3.econ.state.3 <- lmer(cong_econ_state_rescale ~ as.factor(c_buenosaires) + c_age + as.factor(c_female) + 
                            as.factor(c_P37) + as.factor(c_NSE_AGR)+ as.factor(e_job) + as.factor(e_Partido) + 
                            as.numeric(e_gp_edu_rescale) + (1|c_ncue) + (1|e_n), data=df)
summary(mod3.econ.state.3)
r2.corr.mer(mod3.econ.state.3)

mod3.econ.state.4 <- lmer(cong_econ_state_rescale ~ c_age + as.factor(c_female) + as.numeric(e_gp_edu_rescale) + 
                            (1|c_ncue) + (1|e_n), data=df)
summary(mod3.econ.state.4)
r2.corr.mer(mod3.econ.state.4)

# clean up
rm(list=ls(pattern="mod")); rm(relgrad, leg.df); gc() 


### Table A5: Bayesian replication of MLE results in Figure 2/Table A1
# This is memory-intensive and takes a while (approx 4 hours on my machine). I recommend storing the 
# coefficients and HPDs after each model, clearing the memory, then re-importing and plotting at the end. 

# MCMCglmm only uses 1 chain, so we'll run each model twice to generate two chains, then evaluate convergence.

### ideology
df.P41 <- completeFun(df, c("cong_P41_rescale","c_buenosaires","c_age","c_female","c_P37","c_NSE_AGR","e_job","e_Partido",
                            "c_ncue","e_n"))
BV.P41 <- gelman.prior(cong_P41_rescale ~ as.factor(c_buenosaires) + c_age + as.factor(c_female) + 
                         as.factor(c_P37) + as.factor(c_NSE_AGR)+ as.factor(e_job) + as.factor(e_Partido),data=df.P41)
my.pri.P41 <- list(B = list(mu = rep(0,15), V=BV.P41),
                   R=list(V=1, nu=0.05), G = list(G1=list(V=1,nu=0.05), G2=list(V=1, nu=0.05)))
bmod.P41 <- MCMCglmm(cong_P41_rescale ~ as.factor(c_buenosaires) + c_age + as.factor(c_female) + 
                       as.factor(c_P37) + as.factor(c_NSE_AGR)+ as.factor(e_job) + as.factor(e_Partido), 
                     random=~c_ncue + e_n, data=df.P41, prior=my.pri.P41)
bmod.P41.2 <- MCMCglmm(cong_P41_rescale ~ as.factor(c_buenosaires) + c_age + as.factor(c_female) +
                         as.factor(c_P37) + as.factor(c_NSE_AGR)+ as.factor(e_job) + as.factor(e_Partido), 
                       random=~c_ncue + e_n, data=df.P41, prior=my.pri.P41)

# convergence and autocorrelation
autocorr.diag(bmod.P41$Sol)        # should be less than .1
effectiveSize(bmod.P41$Sol)        # should be around 1,000
autocorr.diag(bmod.P41.2$Sol)      # should be less than .1
effectiveSize(bmod.P41.2$Sol)      # should be around 1,000
gelman.diag(mcmc.list(bmod.P41$Sol, bmod.P41.2$Sol)) # should be less than 1.1
# plotTrace(bmod.P41$Sol)
# plotTrace(bmod.P41.2$Sol)

# results
bmod.P41.chain <- as.mcmc(do.call(rbind,list(bmod.P41$Sol, bmod.P41.2$Sol)))
posterior.mode(bmod.P41.chain)
HPDinterval(bmod.P41.chain)
rm(df.P41,bmod.P41,bmod.P41.2,bmod.P41.chain); gc()


### democracy
df.P50 <- completeFun(df, c("cong_P50","c_buenosaires","c_age","c_female","c_P37","c_NSE_AGR","e_job","e_Partido",
                            "c_ncue","e_n"))
BV.P50 <- gelman.prior(cong_P50 ~ as.factor(c_buenosaires) + c_age + as.factor(c_female) + 
                         as.factor(c_P37) + as.factor(c_NSE_AGR)+ as.factor(e_job) + as.factor(e_Partido),data=df.P50)
my.pri.P50 <- list(B = list(mu = rep(0,15), V=BV.P50),
                   R=list(V=1, nu=0.05), G = list(G1=list(V=1,nu=0.05), G2=list(V=1, nu=0.05)))
bmod.P50 <- MCMCglmm(cong_P50 ~ as.factor(c_buenosaires) + c_age + as.factor(c_female) + 
                       as.factor(c_P37) + as.factor(c_NSE_AGR)+ as.factor(e_job) + as.factor(e_Partido), 
                     random=~c_ncue + e_n, data=df.P50,  family="categorical",burnin=5000,nitt=25000,thin=20,
                     prior=my.pri.P50)
bmod.P50.2 <- MCMCglmm(cong_P50 ~ as.factor(c_buenosaires) + c_age + as.factor(c_female) + 
                         as.factor(c_P37) + as.factor(c_NSE_AGR)+ as.factor(e_job) + as.factor(e_Partido), 
                       random=~c_ncue + e_n, data=df.P50,  family="categorical",burnin=5000,nitt=25000,thin=20,
                       prior=my.pri.P50)

# convergence and autocorrelation
autocorr.diag(bmod.P50$Sol)        # should be less than .1
effectiveSize(bmod.P50$Sol)        # should be around 1,000
autocorr.diag(bmod.P50.2$Sol)      # should be less than .1
effectiveSize(bmod.P50.2$Sol)      # should be around 1,000
gelman.diag(mcmc.list(bmod.P50$Sol, bmod.P50.2$Sol)) # should be less than 1.1
# plotTrace(bmod.P50$Sol)
# plotTrace(bmod.P50.2$Sol)

# results
bmod.P50.chain <- as.mcmc(do.call(rbind,list(bmod.P50$Sol, bmod.P50.2$Sol)))
posterior.mode(bmod.P50.chain)
HPDinterval(bmod.P50.chain)
rm(df.P50,bmod.P50,bmod.P50.2,bmod.P50.chain); gc()


### econ_state
df.econ_state <- completeFun(df, c("cong_econ_state_rescale","c_buenosaires","c_age","c_female","c_P37","c_NSE_AGR","e_job","e_Partido",
                                   "c_ncue","e_n"))
BV.econ_state <- gelman.prior(cong_econ_state_rescale ~ as.factor(c_buenosaires) + c_age + as.factor(c_female) + 
                                as.factor(c_P37) + as.factor(c_NSE_AGR)+ as.factor(e_job) + as.factor(e_Partido),data=df.econ_state)
my.pri.econ_state <- list(B = list(mu = rep(0,15), V=BV.econ_state),
                          R=list(V=1, nu=0.05), G = list(G1=list(V=1,nu=0.05), G2=list(V=1, nu=0.05)))
bmod.econ_state <- MCMCglmm(cong_econ_state_rescale ~ as.factor(c_buenosaires) + c_age + as.factor(c_female) + 
                              as.factor(c_P37) + as.factor(c_NSE_AGR)+ as.factor(e_job) + as.factor(e_Partido), 
                            random=~c_ncue + e_n, data=df.econ_state, prior=my.pri.econ_state)
bmod.econ_state.2 <- MCMCglmm(cong_econ_state_rescale ~ as.factor(c_buenosaires) + c_age + as.factor(c_female) +
                                as.factor(c_P37) + as.factor(c_NSE_AGR)+ as.factor(e_job) + as.factor(e_Partido), 
                              random=~c_ncue + e_n, data=df.econ_state, prior=my.pri.econ_state)

# convergence and autocorrelation
autocorr.diag(bmod.econ_state$Sol)        # should be less than .1
effectiveSize(bmod.econ_state$Sol)        # should be around 1,000
autocorr.diag(bmod.econ_state.2$Sol)      # should be less than .1
effectiveSize(bmod.econ_state.2$Sol)      # should be around 1,000
gelman.diag(mcmc.list(bmod.econ_state$Sol, bmod.econ_state.2$Sol)) # should be less than 1.1
# plotTrace(bmod.econ_state$Sol)
# plotTrace(bmod.econ_state.2$Sol)

# results
bmod.econ_state.chain <- as.mcmc(do.call(rbind,list(bmod.econ_state$Sol, bmod.econ_state.2$Sol)))
posterior.mode(bmod.econ_state.chain)
HPDinterval(bmod.econ_state.chain)
rm(df.econ_state,bmod.econ_state,bmod.econ_state.2,bmod.econ_state.chain); gc()


### populism
df.pop <- completeFun(df, c("cong_pop_rescale","c_buenosaires","c_age","c_female","c_P37","c_NSE_AGR","e_job","e_Partido",
                            "c_ncue","e_n"))
BV.pop <- gelman.prior(cong_pop_rescale ~ as.factor(c_buenosaires) + c_age + as.factor(c_female) + 
                         as.factor(c_P37) + as.factor(c_NSE_AGR)+ as.factor(e_job) + as.factor(e_Partido),data=df.pop)
my.pri.pop <- list(B = list(mu = rep(0,15), V=BV.pop),
                   R=list(V=1, nu=0.05), G = list(G1=list(V=1,nu=0.05), G2=list(V=1, nu=0.05)))
bmod.pop <- MCMCglmm(cong_pop_rescale ~ as.factor(c_buenosaires) + c_age + as.factor(c_female) + 
                       as.factor(c_P37) + as.factor(c_NSE_AGR)+ as.factor(e_job) + as.factor(e_Partido), 
                     random=~c_ncue + e_n, data=df.pop, prior=my.pri.pop)
bmod.pop.2 <- MCMCglmm(cong_pop_rescale ~ as.factor(c_buenosaires) + c_age + as.factor(c_female) +
                         as.factor(c_P37) + as.factor(c_NSE_AGR)+ as.factor(e_job) + as.factor(e_Partido), 
                       random=~c_ncue + e_n, data=df.pop, prior=my.pri.pop)

# convergence and autocorrelation
autocorr.diag(bmod.pop$Sol)        # should be less than .1
effectiveSize(bmod.pop$Sol)        # should be around 1,000
autocorr.diag(bmod.pop.2$Sol)      # should be less than .1
effectiveSize(bmod.pop.2$Sol)      # should be around 1,000
gelman.diag(mcmc.list(bmod.pop$Sol, bmod.pop.2$Sol)) # should be less than 1.1
# plotTrace(bmod.pop$Sol)
# plotTrace(bmod.pop.2$Sol)

# results
bmod.pop.chain <- as.mcmc(do.call(rbind,list(bmod.pop$Sol, bmod.pop.2$Sol)))
posterior.mode(bmod.pop.chain)
HPDinterval(bmod.pop.chain)
rm(df.pop,bmod.pop,bmod.pop.2,bmod.pop.chain); gc()

### Order vs liberty
df.P55 <- completeFun(df, c("cong_P55","c_buenosaires","c_age","c_female","c_P37","c_NSE_AGR","e_job","e_Partido",
                            "c_ncue","e_n"))  
BV.P55 <- gelman.prior(cong_P55 ~ as.factor(c_buenosaires) + c_age + as.factor(c_female) + 
                         as.factor(c_P37) + as.factor(c_NSE_AGR)+ as.factor(e_job) + as.factor(e_Partido),data=df.P55)
my.pri.P55 <- list(B = list(mu = rep(0,15), V=BV.P55),
                   R=list(V=1, nu=0.05), G = list(G1=list(V=1,nu=0.05), G2=list(V=1, nu=0.05)))
bmod.P55 <- MCMCglmm(cong_P55 ~ as.factor(c_buenosaires) + c_age + as.factor(c_female) + 
                       as.factor(c_P37) + as.factor(c_NSE_AGR)+ as.factor(e_job) + as.factor(e_Partido), 
                     random=~c_ncue + e_n, data=df.P55, family="categorical",burnin=5000,nitt=25000,thin=20,
                     prior=my.pri.P55)
bmod.P55.2 <- MCMCglmm(cong_P55 ~ as.factor(c_buenosaires) + c_age + as.factor(c_female) +
                         as.factor(c_P37) + as.factor(c_NSE_AGR)+ as.factor(e_job) + as.factor(e_Partido), 
                       random=~c_ncue + e_n, data=df.P55, family="categorical",burnin=5000,nitt=25000,thin=20,
                       prior=my.pri.P55)

# convergence and autocorrelation
autocorr.diag(bmod.P55$Sol)        # should be less than .1
effectiveSize(bmod.P55$Sol)        # should be around 1,000
autocorr.diag(bmod.P55.2$Sol)      # should be less than .1
effectiveSize(bmod.P55.2$Sol)      # should be around 1,000
gelman.diag(mcmc.list(bmod.P55$Sol, bmod.P55.2$Sol)) # should be less than 1.1
# plotTrace(bmod.P55$Sol)
# plotTrace(bmod.P55.2$Sol)

# store for plotting later and cleanup
bmod.P55.chain <- as.mcmc(do.call(rbind,list(bmod.P55$Sol, bmod.P55.2$Sol)))
posterior.mode(bmod.P55.chain)
HPDinterval(bmod.P55.chain)


# rm(list=ls())
