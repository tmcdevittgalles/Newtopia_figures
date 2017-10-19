### Newtopia figures #####
library(tidyr)
library(dplyr)
library(ggplot2)
library(lattice)
library(lme4)


## Imputting data
setwd("/Users/travismcdevitt-galles/Desktop/Current_Projects/Newtopia")

raw.df <- read.csv("newtopia.csv")

names(raw.df)

newt.df <- raw.df[complete.cases(raw.df$ttx_cm2),] ## removes all rows with NAs

## Figure 3, Violin plot of Newt TTX data between species


ttx.dist.log <- ggplot(newt.df, aes(y=log10(ttx_cm2), x = species, fill=species)) + 
    geom_violin(scale = "width") + 
    ylab(expression(TTX~concentration~log["10"](ng~cm^{"-2"}))) + theme_classic()

ttx.dist.log <- ttx.dist.log + theme(axis.title.x = element_text(size = rel(1.8)),
                      axis.text.x  = element_text(vjust=0.5, size=16),
                      axis.title.y = element_text(size = rel(1.8)),
                     axis.text.y  = element_text(vjust=0.5, size=16),
                      axis.line.x = element_line(color="black") ,
                      axis.line.y = element_line(color="black"),
                      axis.ticks.y = element_line(color="black"),
                      axis.ticks.x = element_line(color="black"),
                     legend.position = 'none') + xlab("Species")
    
    
ttx.dist.log


## Figure 1 but with adjustments based on random site effects

at1 <- lmer(log10(ttx_cm2+1) ~ 1 + (1|site) , data=newt.df)
at1
 
summary(at1)
site.effect <- ranef(at1)

site <- as.factor(row.names(site.effect$site))
effect <- as.numeric(site.effect$site$`(Intercept)`)
site.effect.df <- cbind.data.frame(site,effect)

vio.df <- newt.df[,c(7,11,52)]

vio.df$logTTX <- log10(vio.df$ttx_cm2 + 1)

vio.df <- plyr::join ( vio.df , site.effect.df, by="site" )

vio.df$adjTTX <- vio.df$logTTX + vio.df$effect

vio.df  <- vio.df[-304,]

ttx.dist.adj <- ggplot(vio.df, aes(y=adjTTX, x = species, fill=species)) + 
    geom_violin(scale = "width") + 
    ylab(expression(TTX~concentration~log["10"](ng~cm^{"-2"}))) + theme_classic()

ttx.dist.adj <- ttx.dist.adj + theme(axis.title.x = element_text(size = rel(1.8)),
                                     axis.text.x  = element_text(vjust=0.5, size=16),
                                     axis.title.y = element_text(size = rel(1.8)),
                                     axis.text.y  = element_text(vjust=0.5, size=16),
                                     axis.line.x = element_line(color="black") ,
                                     axis.line.y = element_line(color="black"),
                                     axis.ticks.y = element_line(color="black"),
                                     axis.ticks.x = element_line(color="black"),
                                     legend.position = 'none') + xlab("Species")


ttx.dist.log
ttx.dist.adj






newt.df$sSVL <- (newt.df$svl - mean(newt.df$svl))/sd(newt.df$svl)

newt.df <- raw.df[complete.cases(raw.df$ttx_cm2),]

newt.df$sTTX <- (log10(newt.df$ttx_cm2+1) - mean(log10(newt.df$ttx_cm2+1)))/sd(log10(newt.df$ttx_cm2+1))



newt.df <- newt.df[-304,]

newt.df$LogTTX <- log10(newt.df$ttx_cm2 )



newt.df$sTTX <- (newt.df$LogTTX - mean(newt.df$LogTTX)) / sd(newt.df$LogTTX)          
                  
at1 <- glmer(ALLPARASRICH~ species + sex + sSVL  + (1|site) +
                 (1|collection),  family='poisson', data=newt.df)


summary(at1)



## Figure 4 and 5 combinded 

newt.df <- newt.df[-304,] # removing indivudal with 0 ttx

# Figure 4 : Richness # 

# model to obtain the coeffients
rich <- glmer(ALLPARASRICH ~ sTTX + scale(svl) + sex +  (1|site) , data=newt.df,
              family='poisson') 

slop <- fixef(rich)

inter <- ranef(rich)$site

# model to obtain the residuals 
resid.m <- glmer(ALLPARASRICH ~  scale(svl) + sex +  (1|site) , data=newt.df,
              family='poisson') 


rich.resid <- resid(resid.m)
ttx <- seq(-4,4,length.out = 100 )

at2<-cbind.data.frame(inter=inter[,1])

dum.df <- expand.grid(inter[,1], ttx)

colnames(dum.df) <- c("site", "sTTX")

dum.df$pred <- rep(NA, nrow(dum.df))


for(i in 1:nrow(dum.df)){
    dum.df[i,3] <- exp(0.1857076+dum.df[i,1]) + dum.df[i,2] * -0.1001761
    
}


ttxrich <- ggplot( data=newt.df, aes(x=sTTX, y=ALLPARASRICH)) + 
     geom_point(alpha=.75,size=3, shape=15) +theme_classic() +
    geom_line(data=dum.df, aes(x=sTTX, y= pred,color=as.factor(site))) + 
    theme(legend.position = 'none') 

ttxrich <- ttxrich + theme(axis.title.x = element_text(size = rel(3)),
                axis.text.x  = element_text(vjust=0.5, size=18),
                axis.title.y = element_text(size = rel(3)),
                axis.text.y  = element_text(vjust=0.5, size=18),
                axis.line.x = element_line(color="black") ,
                axis.line.y = element_line(color="black"),
                axis.ticks.y = element_line(color="black"),
                axis.ticks.x = element_line(color="black"),
                legend.position = 'none',
                panel.border = element_rect(color = "black", fill=NA, size=2.5)) + ylab("Parasite richness") +
    xlab("TTX (z-score)")

ttxrich



ttxrich <- ttxrich +  annotate("text", x=-4 ,y=5, label="C", size=10)


bd <- glmer(BD_POZ ~ sTTX + (1|site) , data=newt.df,
              family='binomial', glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000000))) 

slop <- fixef(bd)

inter <- ranef(bd)$site


## New data frame

ttx <- seq(-4,4,length.out = 100 )

at2<-cbind.data.frame(inter=inter[,1])
   
             

dum.df <- expand.grid(inter[,1], ttx)

colnames(dum.df) <- c("site", "sTTX")

dum.df$pred <- rep(NA, nrow(dum.df))

for(i in 1:nrow(dum.df)){
    dum.df[i,3] <- (-1.5188245 + dum.df[i,1]) + (dum.df[i,2] * -0.4181254)
    
}

at4 <- boot::inv.logit(dum.df$pred)

bd.binom <-ggplot(dum.df,aes(x=sTTX, y = at4)) +
    geom_line(aes(color=as.factor(dum.df$site))) +theme(legend.position = 'none')
       

bd.binom <- bd.binom + geom_point(data=newt.df,aes(x=sTTX,y=BD_POZ, alpha=.05))

bd.binom <- bd.binom + theme_classic() + theme(axis.title.x = element_text(size = rel(2.5)),
                                  axis.text.x  = element_text(vjust=0.5, size=16),
                                  axis.title.y = element_text(size = rel(2.5)),
                                  axis.text.y  = element_text(vjust=0.5, size=16),
                                  axis.line.x = element_line(color="black") ,
                                  axis.line.y = element_line(color="black"),
                                  axis.ticks.y = element_line(color="black"),
                                  axis.ticks.x = element_line(color="black"),
                                  legend.position = 'none') + ylab("Bd infection status") +
    xlab("TTX concentration (scaled and centered)") 
bd.binom<-bd.binom + annotate("text", x=-4 ,y=1, label="A", size=10)

bd.binom

title1 <- gridExtra::textGrob("Scaled Log transformed TTX concentration",
                   gp=gpar(fontsize=20))
gridExtra::grid.arrange(ttxrich , 
                        bd.binom, nrow=1,
                        bottom = title1 )

## RV
rv1 <- glmer(RV_POZ ~ 1+   (1|site) , data=newt.df,
             family='binomial')

RV <- glmer(RV_POZ ~ sTTX + scale(svl) + sex +  (1|site) , data=newt.df,
            family='binomial')

slop <- fixef(RV)

sinter <- ranef(RV)$site
boot::inv.logit(-0.5056 )


## New data frame

ttx <- seq(-4,4,length.out = 327 )



dum.df <- expand.grid(inter[,1], ttx)

colnames(dum.df) <- c("site", "sTTX")


pred <-rep(NA,327)
for(i in 1:327){
   pred[i] <- (-2.6627162) + ttx[i]*-0.5056 
    
}

at4 <- boot::inv.logit(pred)

weird.df<- cbind.data.frame(ttx, at4)
rv.binom <-ggplot(data=newt.df,aes(x=sTTX,y=RV_POZ)) + geom_point(alpha=.5) +
   theme_classic() 


rv.binom <- rv.binom + geom_line(aes(x=weird.df$ttx, y = weird.df$at4))

rv.binom <- rv.binom + theme_classic() + ylab( "Ranavirus infection status" ) + 
    xlab( "TTX concentration (scaled and centered)"  )+theme(legend.position = 'none')


rv.binom <- rv.binom +annotate("text", x=-4 ,y=1, label="B", size=10)

rv.binom <- rv.binom + theme(axis.title.x = element_text(size = rel(2.5)),
                 axis.text.x  = element_text(vjust=0.5, size=16),
                 axis.title.y = element_text(size = rel(2.5)),
                 axis.text.y  = element_text(vjust=0.5, size=16),
                 axis.line.x = element_line(color="black") ,
                 axis.line.y = element_line(color="black"),
                 axis.ticks.y = element_line(color="black"),
                 axis.ticks.x = element_line(color="black"),
                 legend.position = 'none')
### exploring a figure for Cosmo
### 
### 

cosm.m <- glmer((P.OXYD) ~ sTTX + scale(svl) + sex + (1|site) , data=newt.df,
family='poisson') 

slop <- fixef(cosm.m)

inter <- ranef(cosm.m)$site

slop
ttx <- seq(-4,4,length.out = 100 )

at2<-cbind.data.frame(inter=inter[,1])

dum.df <- expand.grid(inter[,1], ttx)

colnames(dum.df) <- c("site", "sTTX")

dum.df$pred <- rep(NA, nrow(dum.df))


for(i in 1:nrow(dum.df)){
    dum.df[i,3] <- (-3.6618700 +dum.df[i,1]) + dum.df[i,2] * 0.4557495
    
}


ttxCosm <- ggplot( data=newt.df, aes(x=sTTX, y=log10(P.OXYD+1))) + 
    geom_point(alpha=.5) +theme_classic() +
    geom_line(data=dum.df, aes(x=sTTX, y= log10(exp(pred)+1),color=as.factor(site))) + 
    theme(legend.position = 'none') 

ttxCosm  <- ttxCosm + theme(axis.title.x = element_text(size = rel(1.8)),
                           axis.text.x  = element_text(vjust=0.5, size=16),
                           axis.title.y = element_text(size = rel(1.8)),
                           axis.text.y  = element_text(vjust=0.5, size=16),
                           axis.line.x = element_line(color="black") ,
                           axis.line.y = element_line(color="black"),
                           axis.ticks.y = element_line(color="black"),
                           axis.ticks.x = element_line(color="black"),
                           legend.position = 'none') +
    ylab(expression(COSM~abundance~~log["10"]~transformed~+1)) +
    xlab(expression(TTX~concentration~~log["10"]~transformed~scaled~and~centered))


ttxCosm  <- ttxCosm +  annotate("text", x=-4 ,y=3.5, label="D", size=10)


##### Combination of all ttx parasite relationships
bd.binom # Figure A BD
rv.binom # Figure B RV
ttxrich # Figure C Richness
ttxCosm # figure D Cosm

grid.arrange(bd.binom,rv.binom,nrow=1)

### building a model to explore patterns of parasite richness and ttx


Rfx.m <- glmer(ALLPARASRICH ~ (1|site), family = 'poisson', data= newt.df)


inter <- random.effects(Rfx.m)

inter <- (inter$site)

ttx <- seq(-3,3.4,length.out = 100 )



dum.df <- expand.grid(inter$site[,1], ttx)

colnames(dum.df) <- c("site", "sTTX")
head(dum.df)


dum.df$pred <- rep(NA, nrow(dum.df))


for( i in 1:nrow(dum.df)){
    
    dum.df[i,3] <- dum.df[i,1] + -0.116 * dum.df[i,2]
    
}


dum.df$tran <- exp(dum.df$pred)


richran <- ggplot(dum.df, aes(x= sTTX, y= tran, color= as.factor(site))) + geom_line()+
    theme_classic() + theme( legend.position = 'none')

richran <- richran + geom_point(data=newt.df, aes(x=sTTX, y = ALLPARASRICH,alpha=.05), color="black")

richran + xlab('TTX concentration (scaled and centered)') + ylab(" Host parasite richness")



spp.rich <- ggplot(newt.df, aes(x = species, y=ALLPARASRICH, fill = species))+
    geom_boxplot()

spp.rich <- spp.rich + theme_classic() + geom_jitter(width=.05, alpha=.1) + ylab("Host parasite richness")


### figure 2

## Making a new data frame to hold the transformed data

com.df <- data.frame(TAGR.m = as.numeric(),
                     TATO.m = as.numeric(),
                     TAGR.se = as.numeric(),
                     TATO.se = as.numeric(),
                     Parasite = factor()
                     )


for( i in 1:37 ) {
    
    at.df <- newt.df[ , c(11, 14+i) ]
    dum.m <- tapply( at.df[,2] , at.df[ ,1], mean, na.rm=T)
    dum.se <- tapply( at.df[,2] , at.df[ ,1], se)
    com.df[ i ,1 ] <- dum.m[1]
    com.df[ i ,2 ] <- dum.m[2]
    com.df[ i ,3 ] <- dum.se[1]
    com.df[ i ,4 ] <- dum.se[2]
}


par.df <- com.df[c(1,7,16,17,18,20,21,26:28,30,31, 34 ,36),]

par.names <- c("ACSP","CLSP", "MEMO","NEMA", "OXYD", "RHAB","RION", "TARP",
             "TITR", 'EITA', "NYSP","COSM", "BD", "RV"   )


par.df$parasites <- par.names

mic.df <- par.df[c(8:14),]

mic.df <- mic.df[-c(4,5),]

mac.df <- par.df[c(1:7),]

mac.df <- mac.df[-c(1,4,21),]

mac.df <- mac.df[-5,]


mac.par <- ggplot(mac.df,aes(x=log10(TATO.m+1) , y= log10(TAGR.m+1) ,color=parasites) )+
    geom_point() + theme_classic() + 
    geom_abline(intercept = 0, slope=1, linetype="dotted")

mac.par <- mac.par + geom_errorbar(aes(ymin = (log10(TAGR.m+1) - (log10(TAGR.se+1)))  , 
                        ymax = log10(TAGR.m+1) + (log10(TAGR.se+1)) ),width=0) +
     geom_errorbarh(aes(xmin = (log10(TATO.m+1) - (log10(TATO.se+1)))  , 
                        xmax = log10(TATO.m+1) + (log10(TATO.se+1)), height=0 ))

mac.par <- mac.par+ 
ylab(expression(italic("T. granulosa")~mean~infection~load~log["10"]~+1)) +
      xlim(c(0,2.7)) + ylim(c(0,2.7)) + 
    xlab(expression(italic("T. torosa")~mean~infection~load~log["10"]~+1))
mac.par<- mac.par + annotate("text", x=0.1 ,y=2.45, label="A", size=10)
mac.par <-mac.par + scale_color_discrete(name="Macroparasites")

### Microparasites

micro.par <- ggplot(mic.df,aes(x=log10(TATO.m+1) , y= log10(TAGR.m+1) ,color=parasites) )+
    geom_point() + theme_classic() + geom_abline(intercept = 0, slope=1,linetype="dotted")

micro.par<-micro.par + geom_errorbar(aes(ymin = (log10(TAGR.m+1) - (log10(TAGR.se+1)))  , 
                                     ymax = log10(TAGR.m+1) + (log10(TAGR.se+1))), width=0) +
    geom_errorbarh(aes(xmin = (log10(TATO.m+1) - (log10(TATO.se+1))) , 
                       xmax = log10(TATO.m+1) + (log10(TATO.se+1)), height=0) )

micro.par <- micro.par+ ylab(expression(italic("T. granulosa")~mean~infection~prevalence~log["10"]~+1)) +xlim(c(0,.12)) +
    ylim(c(0,.12)) +xlab(expression(italic("T. torosa")~mean~infection~prevalence~log["10"]~+1))
micro.par <- micro.par + annotate("text", x=0.01 ,y=.11, label="B", size=10)
micro.par<- micro.par +  scale_color_discrete(name="Microparasites")

gridExtra::grid.arrange(mac.par, micro.par, nrow=1)

at1 <- glmer(P.OXYD ~ scale(svl) + species + sex + 
          scale(log10(ttx_cm2+1)) + (1|site) ,
      data=newt.df, 
      family=poisson)

aggregate(RV.POZ ~ species * site, data = newt.df, FUN= n)

### New and improved figure 1

p_type <- rep(NA,9)
parasite <- rep(NA,9)
tagr.coef <- rep(NA,9)
tagr.se <- rep(NA,9)
tato.coef <- rep(NA,9)
tato.se <- rep(NA,9) 
fig1_df <- cbind.data.frame(p_type, parasite,tagr.coef,tagr.se,tato.coef,tato.se)
fig1_df[1:5,1] <- "Microparasites"
fig1_df[6:9,1] <- "Macroparasites"

# obtaining coefficents for BD#

dum.df <- newt.df[,c(7,11,48)]
fig1_df[1,2] <- "BD"
bd <- glmer(BD_POZ ~  0 +species + (1|site) , data=newt.df,
             family=binomial, glmerControl(optimizer="bobyqa",
                                            optCtrl = list(maxfun = 10000000)))


site.effect <- ranef(bd)
spp.effect <- fixef(bd)
site <- as.factor(row.names(site.effect$site))
effect <- as.numeric(site.effect$site$`(Intercept)`)

site.effect.df <- cbind.data.frame(site,effect)

species <- c('TATO', "TAGR")
Coef <- c(-1.359,-1.618)

fixed.eff.df <- cbind.data.frame(species,Coef)

dum1.df <- plyr::join(dum.df,site.effect.df,by="site")

dum2.df <- plyr::join(dum1.df, fixed.eff.df, by="species")

dum2.df$pred <- dum2.df$effect + dum2.df$Coef

dum2.df$Tpred <- boot::inv.logit(dum2.df$pred)

fig1_df[1,3] <- 0.2201948
fig1_df[1,4] <-  0.008142008
fig1_df[1,5] <- 0.2201334 
fig1_df[1,6] <- 0.008554723

# obtaining coefficents for RV
# 
fig1_df[2,2] <- "RV"
rv <- glmer(BD_POZ ~ 0 + species+(1|site) , data=newt.df,
            family='binomial', glmerControl(optimizer="bobyqa",
                                            optCtrl = list(maxfun = 10000000)))

summary(rv)

rv <- glm(BD_POZ ~ 0 + species , data=newt.df,
            family='binomial')

summary(rv)

site.effect <- ranef(rv)
spp.effect <- fixef(rv)
site <- as.factor(row.names(site.effect$site))
effect <- as.numeric(site.effect$site$`(Intercept)`)

site.effect.df <- cbind.data.frame(site,effect)

species <- c('TATO', "TAGR")
Coef <- c(-2.3238,-2.8717  )

fixed.eff.df <- cbind.data.frame(species,Coef)

dum1.df <- plyr::join(dum.df,site.effect.df,by="site")

dum2.df <- plyr::join(dum1.df, fixed.eff.df, by="species")

dum2.df$pred <- dum2.df$effect + dum2.df$Coef

dum2.df$Tpred <- boot::inv.logit(dum2.df$pred)

tapply(dum2.df$Tpred, dum2.df$species, mean, na.rm=T)

tapply(dum2.df$Tpred, dum2.df$species, se )
fig1_df[2,3] <-0.07938859
fig1_df[2,4] <- 0.004140802 
fig1_df[2,5] <- 0.10220451 
fig1_df[2,6] <- 0.005366071 

# obtaining coefficents for TARP
# 
fig1_df[3,2] <- "TAPR"
tapr <- glmer(P.TAPR ~ 0 + species +(1|site) , data=newt.df,
            family='binomial', glmerControl(optimizer="bobyqa",
                                            optCtrl = list(maxfun = 10000000)))

summary(tapr)
site.effect <- ranef(tapr)
spp.effect <- fixef(tapr)
site <- as.factor(row.names(site.effect$site))
effect <- as.numeric(site.effect$site$`(Intercept)`)

site.effect.df <- cbind.data.frame(site,effect)

species <- c('TATO', "TAGR")
Coef <- c(-8.864, -7.752  )

fixed.eff.df <- cbind.data.frame(species,Coef)

dum1.df <- plyr::join(dum.df,site.effect.df,by="site")

dum2.df <- plyr::join(dum1.df, fixed.eff.df, by="species")

dum2.df$pred <- inv.logit(dum2.df$effect) +  inv.logit(dum2.df$Coef)

dum2.df$Tpred <- boot::inv.logit(dum2.df$pred)

tapply(dum2.df$pred, dum2.df$species, mean, na.rm=T)

tapply(dum2.df$Tpred, dum2.df$species, se )

fig1_df[3,3] <- 0.0006872735
fig1_df[3,4] <-  boot::inv.logit(1.86231 )
fig1_df[3,5] <- 0.0001758270 
fig1_df[3,6] <- boot::inv.logit(1.89183 )


# obtaining coefficents for EITA
# 
fig1_df[4,2] <- "EITA"
eita <- glmer(P.EITA ~ 0 + species +scale(svl) + sex +(1|site) , data=newt.df,
              family='binomial', glmerControl(optimizer="bobyqa",
                                              optCtrl = list(maxfun = 10000000)))

summary(eita)

fig1_df[4,3] <- boot::inv.logit(-7.8611)
fig1_df[4,4] <- boot::inv.logit( 2.0900)
fig1_df[4,5] <- boot::inv.logit(-8.4610 )
fig1_df[4,6] <- boot::inv.logit(2.1115 )


# obtaining coefficents for TITR
# 
fig1_df[5,2] <- "TITR"
titr <- glmer(P.TITR ~ species +scale(svl) + sex +(1|site) , data=newt.df,
              family='binomial', glmerControl(optimizer="bobyqa",
                                              optCtrl = list(maxfun = 10000000)))

summary(titr)

fig1_df[5,3] <- boot::inv.logit(-11.5018 )
fig1_df[5,4] <- boot::inv.logit( 2.4751)
fig1_df[5,5] <- boot::inv.logit(-10.5846)
fig1_df[5,6] <- boot::inv.logit(2.1998 )


# obtaining coefficents for COSM
# 
fig1_df[6,2] <- "COSM"
cosm <- glm(P.OXYD ~ 0 + species, data=newt.df,
              family='poisson')

summary(cosm)

fig1_df[6,4] <- exp(-4.47703) 
fig1_df[6,5] <-  exp(0.73213 )
fig1_df[6,4] <- exp(-3.18555) 
fig1_df[6,5] <-  exp(0.72777 )


# obtaining coefficents for RHAB
# 
fig1_df[7,2] <- "RHAB"
rhab <- glmer(P.RHAB ~ 0 + species +scale(svl) + sex +(1|site) , data=newt.df,
              family='poisson', glmerControl(optimizer="bobyqa",
                                             optCtrl = list(maxfun = 10000000)))

summary(rhab)

fig1_df[7,3] <- exp(0.53755)
fig1_df[7,4] <-   exp(0.46376 )
fig1_df[7,5] <-  exp(-0.24229  )
fig1_df[7,6] <- exp(0.46548)

# obtaining coefficents for MEMO
# 
fig1_df[8,2] <- "MEMO"
memo<- glmer(P.MEMO ~ 0 + species  +(1|site) , data=newt.df,
              family='poisson', glmerControl(optimizer="bobyqa",
                                             optCtrl = list(maxfun = 10000000)))

summary(memo)

fig1_df[8,3] <- exp(-3.6606 ) 
fig1_df[8,4] <-    exp(0.8444)
fig1_df[8,5] <-  exp(-3.5031) 
fig1_df[8,6] <- exp(0.7995) 

# obtaining coefficents for CLSP
# 
fig1_df[9,2] <- "CLSP"
clsp<- glmer(P.CLSP ~ 0 + species +scale(svl) + sex +(1|site) , data=newt.df,
             family='poisson', glmerControl(optimizer="bobyqa",
                                            optCtrl = list(maxfun = 10000000)))

summary(clsp)
fig1_df[9,3] <- exp(-3.2743)
fig1_df[9,4] <-   exp( 0.7702) 
fig1_df[9,5] <-  exp(-4.1613 )
fig1_df[9,6] <- exp(0.7673) 


at1 <- boot::inv.logit(fig1_df$Coef[1:10])
at2 <- exp(fig1_df$Coef[11:18])
tran.Coef <- c(at1, at2)

fig1_df <- cbind.data.frame(fig1_df,tran.Coef)


at1 <- boot::inv.logit(fig1_df$SE[1:10])
at2 <- exp(fig1_df$SE[11:18])
tran.SE <- c(at1, at2)

fig1_df <- cbind.data.frame(fig1_df,tran.SE)



mac.par <- ggplot(subset(fig1_df,fig1_df$p_type=="Macroparasites"),aes(x=(tato.coef) , y= (tagr.coef) ,color=parasite))+
    geom_point() + theme_classic() + geom_abline(intercept = 0, slope=1)

mac.par 

mac.par<-mac.par + geom_errorbar(aes(ymin = (tagr.coef - tagr.se)  , 
                                     ymax = (tagr.coef + tagr.se) ))+
    geom_errorbarh(aes(xmin = (tato.coef - tato.se) , 
                       xmax = (tato.coef + tato.se)  ))

mac.par <- mac.par+ ylab("T. granulosa mean infection load log10(y+1)") +xlim(c(0,2.7)) +
    ylim(c(0,2.7)) +xlab("T. torosa mean infection load log10(x+1)")
mac.par<- mac.par + annotate("text", x=0.1 ,y=2.45, label="A", size=10)
mac.par <-mac.par + scale_color_discrete(name="Macroparasites")

micro.par <- ggplot(mic.df,aes(x=log10(TATO.m+1) , y= log10(TAGR.m+1) ,color=parasites) )+
    geom_point() + theme_classic() + geom_abline(intercept = 0, slope=1)

micro.par<-micro.par + geom_errorbar(aes(ymin = (log10(TAGR.m+1) - (log10(TAGR.se+1)))  , 
                                         ymax = log10(TAGR.m+1) + (log10(TAGR.se+1)) )) +
    geom_errorbarh(aes(xmin = (log10(TATO.m+1) - (log10(TATO.se+1)))  , 
                       xmax = log10(TATO.m+1) + (log10(TATO.se+1)) ))

micro.par <- micro.par+ ylab("T. granulosa mean infection load log10(y+1)") +xlim(c(0,.12)) +
    ylim(c(0,.12)) +xlab("T. torosa mean infection load log10(x+1)")
micro.par <- micro.par + annotate("text", x=0.01 ,y=.11, label="B", size=10)
micro.par<- micro.par +  scale_color_discrete(name="Microparasites")





## Figure 1 with each point being a site
## 

site.df <- subset(newt.df, newt.df$site == "PRPND015" | 
                      newt.df$site == "PRNTHOWL" |
                      newt.df$site == "PRNTHMIT" | 
                      newt.df$site == "PRNTHIDK" | 
                      newt.df$site == "PRNTHCHR" | 
                      newt.df$site == "PRNTH4" | 
                      newt.df$site == "PRNTH2" |
                      newt.df$site == "MUD46" |
                      newt.df$site == "MUD41" |
                      newt.df$site == "MUD40" | 
                      newt.df$site == "GDPND004" | 
                      newt.df$site == "CA-BN018" |
                      newt.df$site == "CA-BN017" |
                      newt.df$site == "CA-BN003" |
                      newt.df$site == "BNPND011" |
                      newt.df$site == "BNPND002" |
                      newt.df$site == "5CN1" |
                      newt.df$site == "5CN2")
# micro parasites

bd.df <- aggregate(BD_POZ ~ species * site, data = site.df, FUN= mean)
bd.df$bd_se <- aggregate(BD_POZ  ~ species * site, data = site.df, FUN= se)[,3]


RV.df <- aggregate(RV_POZ ~ species * site, data = site.df, FUN= mean)
RV.df$bd_se <- aggregate(RV_POZ  ~ species * site, data = site.df, FUN= se)[,3]

eita.df <- aggregate(P.EITA ~ species * site, data = site.df, FUN= mean)
eita.df$bd_se <- aggregate(P.EITA ~ species * site, data = site.df, FUN= se)[,3]

tarp.df <- aggregate(P.TAPR ~ species * site, data = site.df, FUN= mean)
tarp.df$bd_se <- aggregate(P.TAPR ~ species * site, data = site.df, FUN= se)[,3]


titr.df <- aggregate(P.TITR ~ species * site, data = site.df, FUN= mean)
titr.df$bd_se <- aggregate(P.TITR ~ species * site, data = site.df, FUN= se)[,3]


fig.1.df <- read.csv(file.choose(), header= T)


fig1.df <- subset(fig.1.df, fig.1.df$TAGR != 0 | fig.1.df$TATO != 0 )

site.fig1<- ggplot(fig1.df, aes( x= TATO, y= TAGR , color=Parasite, shape= Parasite)) +
    geom_point(size=5) + theme_classic() +  
    geom_abline(intercept = 0, slope=1, linetype="dotted")

site.fig1  + geom_errorbar(aes(ymin = TAGR - TAGR_SE  , 
                                ymax = TAGR + TAGR_SE ),width=0, alpha=.5 ) +
    geom_errorbarh (aes(xmin = TATO - TATO_SE  , 
                       xmax = TATO + TATO_SE  ),height=0, alpha=.5)


# macro parasites

OXYD.df <- aggregate(P.OXYD ~ species * site, data = site.df, FUN= mean)
OXYD.df$bd_se <- aggregate(P.OXYD  ~ species * site, data = site.df, FUN= se)[,3]


MEMO.df <- aggregate(P.MEMO ~ species * site, data = site.df, FUN= mean)
MEMO.df$bd_se <- aggregate(P.MEMO ~ species * site, data = site.df, FUN= se)[,3]

rhab.df <- aggregate(P.RHAB ~ species * site, data = site.df, FUN= mean)
rhab.df$bd_se <- aggregate(P.RHAB ~ species * site, data = site.df, FUN= se)[,3]

clsp.df <- aggregate(P.CLSP ~ species * site, data = site.df, FUN= mean)
clsp.df$bd_se <- aggregate(P.CLSP ~ species * site, data = site.df, FUN= se)[,3]


titr.df <- aggregate(P.TITR ~ species * site, data = site.df, FUN= mean)
titr.df$bd_se <- aggregate(P.TITR ~ species * site, data = site.df, FUN= se)[,3]


fig.1.df <- read.csv(file.choose(), header= T)


fig1.df <- subset(fig.1.df, fig.1.df$TAGR != 0 | fig.1.df$TATO != 0 )

site.fig1 <- ggplot(fig1.df, aes( x= TATO, y= TAGR , color=Parasite)) +
    geom_point(size=5) + theme_classic() +  
    geom_abline(intercept = 0, slope=1, linetype="dotted")

site.fig1 <- site.fig1  + geom_errorbar(aes(ymin = TAGR - TAGR_SE  , 
                               ymax = TAGR + TAGR_SE ),width=0, alpha=.5 ) +
    geom_errorbarh (aes(xmin = TATO - TATO_SE  , 
                        xmax = TATO + TATO_SE  ),height=0, alpha=.5)


site.fig1  <- site.fig1 + ylab(expression(italic("T. granulosa")~mean~site~level~infection~prevalence)) +
  xlab(expression(italic("T. torosa")~mean~site~level~infection~prevalence))
site.fig1 <- site.fig1 + annotate("text", x=-0 ,y= 1, label="B", size=10)
site.fig1  <- site.fig1 +  scale_color_discrete(name="Microparasites")

## Macro


fig.mac <- read.csv(file.choose(), header= T)

fig.mac.1 <- subset(fig.mac, fig.mac$TAGR != 0 | fig.mac$TATO != 0 )

mac.fig1<- ggplot(fig.mac.1, aes( x= log10(TATO+1), y= log10(TAGR+1) , color=Parasite)) +
    geom_point(size=5) + theme_classic() +  
    geom_abline(intercept = 0, slope=1, linetype="dotted")

mac.fig1 <- mac.fig1  + geom_errorbar(aes(ymin = log10(TAGR+1) - log10(TAGR_SE+1) , 
                               ymax = log10(TAGR+1) + log10(TAGR_SE+1) ),width=0, alpha=.5 ) +
    geom_errorbarh (aes(xmin = log10(TATO+1) - log10(TATO_SE+1)  , 
                        xmax = log10(TATO+1) + log10(TATO_SE+1)  ),height=0, alpha=.5)

mac.fig1 <- mac.fig1+ 
    ylab(expression(italic("T. granulosa")~mean~site~level~infection~load~(log["10"]~+1))) +
    xlim(c(0,2.7)) + ylim(c(0,2.7)) + 
    xlab(expression(italic("T. torosa")~mean~site~level~infection~load~(log["10"]~+1)))
mac.fig1 <- mac.fig1 + annotate("text", x=0 ,y=2.7, label="A", size=10)
mac.fig1 <- mac.fig1 + scale_color_discrete(name="Macroparasites")

mac.fig1 <- mac.fig1 +theme(axis.title.x = element_text(size = rel(1.8)),
               axis.text.x  = element_text(vjust=0.5, size=16),
               axis.title.y = element_text(size = rel(1.8)),
               axis.text.y  = element_text(vjust=0.5, size=16),
               axis.line.x = element_line(color="black") ,
               axis.line.y = element_line(color="black"),
               axis.ticks.y = element_line(color="black"),
               axis.ticks.x = element_line(color="black"))
site.fig1 <- site.fig1 +theme(axis.title.x = element_text(size = rel(1.8)),
                 axis.text.x  = element_text(vjust=0.5, size=16),
                 axis.title.y = element_text(size = rel(1.8)),
                 axis.text.y  = element_text(vjust=0.5, size=16),
                 axis.line.x = element_line(color="black") ,
                 axis.line.y = element_line(color="black"),
                 axis.ticks.y = element_line(color="black"),
                 axis.ticks.x = element_line(color="black"))

grid.arrange(mac.fig1, site.fig1, nrow=1)
