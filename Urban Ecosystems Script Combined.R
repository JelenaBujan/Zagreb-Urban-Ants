#########   ####	####	####	####
#############Final Script for publication############
####Jesovnik and Bujan 2021, Urban Ecosystems #######
############# Date : 24/02/2021 #####################
#########   ####	####	####	####

library(FSA) 
library(emmeans)
library(ggpubr)


####Fig S1####
dhmz.env=read.table("ZGmravDHMZ.txt",header=T)
gric=subset(dhmz.env,LOCATION=="GRIC",select=c( "MONTH","TEMPaver", "PERCIPITATIONmm", "CLOUDINESS"))
maxic=subset(dhmz.env,LOCATION=="MAKSIMIR",select=c( "MONTH","TEMPaver", "PERCIPITATIONmm", "CLOUDINESS"))

par(mfrow=c(1,3))
boxplot(gric$TEMPaver~gric$MONTH,col="#CA0020",boxwex=0.2,
        ylab="Monthly average temperature (?C )",xlab="",ylim=c(0,40), names =c("May","June","July","August"))
boxplot(maxic$TEMPaver~maxic$MONTH,col="#F4A582",add=T,boxwex=0.2,xaxt="n",at=1:4+0.35)
legend(c("Gric","Maksimir"),cex=0.8,pch=c(21,21),pt.bg=c("#CA0020","#F4A582"),bty="n",x="topright")
boxplot(gric$PERCIPITATIONmm~gric$MONTH,col="#0571B0",boxwex=0.2,
        ylab="Monthly rainfall average (mm)",xlab="",ylim=c(0,40), names =c("May","June","July","August"))
boxplot(maxic$PERCIPITATIONmm~maxic$MONTH,col="#92C5DE",add=T,boxwex=0.2,xaxt="n",at=1:4+0.35)
legend(c("Gric","Maksimir"),cex=0.8,pch=c(21,21),pt.bg=c("#0571B0","#92C5DE"),bty="n",x="topright")
boxplot(gric$CLOUDINESS~gric$MONTH,col="#0571B0",boxwex=0.2,ylab="Average cloudyness",xlab="",ylim=c(0,20),
        names =c("May","June","July","August"))
boxplot(maxic$CLOUDINESS~maxic$MONTH,col="#92C5DE",add=T,boxwex=0.2,xaxt="n",at=1:4+0.35)
legend(c("Gric","Maksimir"),cex=0.8,pch=c(21,21),pt.bg=c("#0571B0","#92C5DE"),bty="n",x="topright")


###Figure 2 NEW#####
pd <- position_dodge(0.2) # move them .05 to the left and right
V.plot<-
  ggplot(ZGrich, aes(x=WOODEDareaHA, y=RICHNESShp )) + 
  geom_point(size=3,position=position_jitter(h=0.1,w=0.1), alpha = 0.5)+
  theme(panel.background = element_rect(fill = "white"))+
  theme(panel.background = element_rect(fill = "white"))+
  theme(plot.background = element_blank())+
  theme(panel.border= element_blank())+
  theme(axis.line = element_line(color="black", size = 0.5))+
  theme(axis.text=element_text(size = 11))+
  coord_cartesian(ylim = c(5, 20))+
  labs(y="Species Richness")+
  labs(x="Wooded Area (ha)") +                      
  theme(legend.key = element_rect(fill = "white", color = NA)) 

ZGrich=read.table("ZGmrav per site RICHNESS summary.txt",header=T)
attach(ZGrich)

###Species Richness Table S2A####
SR.glm1=glm(RICHNESShp~MONTH+AREAha+NUMvisitors) 
SR.glm2=glm(RICHNESShp~MONTH+PATHSha+NUMvisitors ) 
SR.glm3=glm(RICHNESShp~MONTH+WOODEDareaHA+NUMvisitors ) 
SR.glm4=glm(RICHNESShp~NUMvisitors+AREAha) 
SR.glm5=glm(RICHNESShp~NUMvisitors+PATHSha) 
SR.glm6=glm(RICHNESShp~NUMvisitors+WOODEDareaHA) 
SR.glm7=glm(RICHNESShp~AREAha)  
SR.glm8=glm(RICHNESShp~PATHSha)
SR.glm9=glm(RICHNESShp~WOODEDareaHA)  
SR.glm10=update(SR.glm9,.~.-WOODEDareaHA)
AICs=AIC(SR.glm1,SR.glm2,SR.glm3,SR.glm4,SR.glm5,SR.glm6,SR.glm7,SR.glm8,
         SR.glm9,SR.glm10 )
###to compare models go to AICs comparisons at the bottom of the code

kruskal.test(ZGrich$RICHNESShp~ZGrich$LOCATION)
PH = dunnTest(ZGrich$RICHNESShp~ZGrich$LOCATION)
detach(ZGrich)

simpson=read.table("D-per-plot.txt",header=T)
kruskal.test(simpson$Simpson~simpson$LOCALITY)
PH = dunnTest(simpson$Simpson~simpson$LOCALITY)

sha=read.table("exponential-Shannon-per-plot.txt",header=T)
kruskal.test(sha$EXShannon~sha$LOCALITY)
PH = dunnTest(sha$EXShannon~sha$LOCALITY)  

###Figure 3 NEW####
boxplot(ZGrich$RICHNESShp~ZGrich$WOODED_cat, boxwex=0.4,col="#92C5DE",ylim=c(0,20),
        ylab="Ant Richness",xlab="",
        names =c("LE","SV","TO","VI","HR","ST","ZR", "BO"))
#legend(c("p = 0.01"),bty="n",x="bottomright")
legend(c("A)  "),cex=1, bty="n",x="topleft")
arrows(3,19,8,19, length=0.05, angle=90, code=0)
arrows(3,18.6,3,19, length=0.05, angle=180, code=0)
arrows(8,18.6,8,19, length=0.05, angle=180, code=0)
points(5.5,19.5,pch=8,cex=0.7)
boxplot(simpson$Simpson~simpson$WOODED_cat, boxwex=0.4,col="#92C5DE",ylim=c(0,1),
        ylab="Inverse Simpson Index",xlab="",
        names =c("LE","SV","TO","VI","HR","ST","ZR", "BO"))
#legend(c("p = 0.02"),bty="n",x="bottomright")
legend(c("B)  "),cex=1, bty="n",x="topleft")
arrows(2,0.85,8,0.85, length=0.05, angle=90, code=0)
arrows(2,0.83,2,0.85, length=0.05, angle=180, code=0)
arrows(8,0.83,8,0.85, length=0.05, angle=180, code=0)
points(5,0.9,pch=8,cex=0.7)
boxplot(sha$EXShannon~sha$WOODED_cat, boxwex=0.4,col="#92C5DE",ylim=c(0,8),
        ylab="Exponential Shannon Index",xlab="",
        names =c("LE","SV","TO","VI","HR","ST","ZR", "BO"))
#legend(c("p = 0.006"),bty="n",x="bottomright")
arrows(2,6.7,8,6.7, length=0.05, angle=90, code=0)
arrows(2,6.5,2,6.7, length=0.05, angle=180, code=0)
arrows(8,6.5,8,6.7, length=0.05, angle=180, code=0)
points(5,7,pch=8,cex=0.7)
legend(c("C)  "),cex=1, bty="n",x="topleft")


###Table S2B######
ZGmrav.baits=read.table("ZGmrav bait activity.txt",header=T)
attach(ZGmrav.baits)
y.pa<-cbind(HITS,TOTALbaits-HITS)

###Activity GLMs####
AB.glm1=glm(y.pa~MONTH+AREAha+NUMvisitors+TEMPaverIR,family=binomial)  
AB.glm2=glm(y.pa~MONTH+PATHSha+NUMvisitors+TEMPaverIR ,family=binomial)  
AB.glm3=glm(y.pa~MONTH+WOODEDareaHA+NUMvisitors+TEMPaverIR ,family=binomial) 
AB.glm4=glm(y.pa~NUMvisitors+AREAha+TEMPaverIR,family=binomial) 
AB.glm5=glm(y.pa~NUMvisitors+PATHSha+TEMPaverIR,family=binomial)  
AB.glm6=glm(y.pa~NUMvisitors+WOODEDareaHA+TEMPaverIR,family=binomial) 
AB.glm7=glm(y.pa~AREAha+TEMPaverIR,family=binomial) 
AB.glm8=glm(y.pa~PATHSha+TEMPaverIR,family=binomial) 
AB.glm9=glm(y.pa~WOODEDareaHA+TEMPaverIR,family=binomial) 
AB.glm10=glm(y.pa~AREAha,family=binomial) 
AB.glm11=glm(y.pa~PATHSha,family=binomial) 
AB.glm12=glm(y.pa~WOODEDareaHA,family=binomial) 
AB.glm13=glm(y.pa~TEMPaverIR,family=binomial)
AB.glm14=update(AB.glm13,.~.-TEMPaverIR)
AICs=AIC(AB.glm1,AB.glm2,AB.glm3,AB.glm4,AB.glm5,AB.glm6,AB.glm7,AB.glm8,
         AB.glm9,AB.glm10,AB.glm11,AB.glm12,AB.glm13,AB.glm14)

E1=resid(AB.glm9,type = "pearson")
sum(E1^2)/(AB.glm9$df.residual)

####Figure 5 NEW####
detach(ZGmrav.baits)
ZGmrav.baits$WOODEDareaHA=as.factor(ZGmrav.baits$WOODEDareaHA)
plot((ZGmrav.baits$PROPhits*100)~ZGmrav.baits$WOODEDareaHA, xlim=c(0.5,8.5),
     ylab="Ant Activity (%)",cex=1.5,
     boxwex=0.4,col="#92C5DE",ylim=c(0,100),xlab="",
     names =c("LE","SV","TO","VI","HR","ST","ZR", "BO"))

plot((ZGmrav.baits$PROPhits*100)~ZGmrav.baits$TEMPaverIR, xlim=c(20,50),xlab="Temperature (°C )",
     ylab="Ant Activity (%)",cex=1.5, bg="red",pch=21)


###AIC comaprisons#####
###compares AICS to find the optimal model (lowest AIC) and the AIC differences from it
myDf=AICs[,1]
AICsNum=AICs[,2]
minAW=min(AICsNum) ####detects the smallest AIC
Delta=AICsNum-minAW ###calculates the differences between min and all other AIC
RL=exp(-0.5*Delta)
wi=RL/sum(RL) ###converting the differences to Akaike weights
Z=data.frame(myDf,AICsNum,Delta,wi) ###adding columns and presenting the results
Z=round(Z,digits=3)
colnames(Z)=c("Df","AIC","SAIC differences","Akaike Weights")
Z



################################################################
####          Code for ZG MRAV data set analysis             ###
####  Written by multiple internet pages and Ana Ješovnik   ####
####         spring-summer 2020 + revision early 2021        ###
################################################################

setwd("~/Desktop/R-ZG-mrav")

install.packages("psych")
install.packages("vegan")
install.packages("ape")

library(readr)
library(lattice)
library(picante)
library(psych)

# RICHNESS ####

# read in data with number of species per plot 
dataset.div <- read_csv("div.csv")
dataset.div$group <- as.factor(dataset.div$group) 

# set an order of plots
levels(dataset.div$group)
dataset.div$group <- ordered(dataset.div$group,
                             levels = c("LE", "SV", "TO" , "VI" , "HR" , "ST" , "ZR" , "BO"))
#look at the data
boxplot(spnum ~ group,
        data = dataset.div,
        col=(c("#92C5DE")),boxwex=0.5,
        ylab="Ant Richness",
        xlab="Locality")

#create high res plot
tiff("Richness.tiff",width=2800,height=4000,units="px",compression="none",res=600)
boxplot(spnum ~ group,
        data = dataset.div,
        col=(c("#92C5DE")),boxwex=0.5,
        ylab="Ant Richness",
        xlab="Locality")
dev.off()


##   ABUNDANCE  ####

# this is total ant abundance (total # of individuals) per plot per month
dataset.abund <- read_csv("abund.csv")
dataset.abund$group <- as.factor(dataset.abund$group)

#set an order of treatments
levels(dataset.abund$group)
dataset.abund$group <- ordered(dataset.abund$group,
                               levels = c("LE", "SV", "TO" , "VI" , "HR" , "ST" , "ZR" , "BO"))
levels(dataset.abund$group)

boxplot(abund ~ group,
        data = dataset.abund,
        col=(c("#92C5DE")),boxwex=0.5,
        ylab="Ant Abundance",
        xlab="Locality")

### Diversity Indices ####


library(picante)

dataset.comm <- read.csv("comm.csv", header = TRUE, row.names = 1)

#attach the groups
landuse <- read_csv("comm-meta.csv")
landuse$group <- as.factor(landuse$group)
attach(landuse)

levels(landuse$group)
landuse$group <- ordered(landuse$group,
                         levels = c("LE", "SV", "TO" , "VI" , "HR" , "ST" , "ZR" , "BO"))


## Shannon Index
H <- diversity(dataset.comm, index = "shannon") 
summary(H) #gives summary statistics for the plots
View(H)

tiff("H-Diversity.tiff",width=2800,height=4000,units="px",compression="none",res=600)
boxplot(H ~ landuse$group,
        data = dataset.comm,
        col=(c("#92C5DE")),boxwex=0.5,main="Ant Diversity",
        ylab="Shannon Index",
        xlab="Locality")
dev.off()


## Exponential Shannon ####
eH <- exp (H)
View(eH)

tiff("expH-Diversity.tiff",width=2800,height=4000,units="px",compression="none",res=600)
boxplot(eH ~ landuse$group,
        data = dataset.comm,
        col=(c("#92C5DE")),boxwex=0.5,main="Ant Diversity",
        ylab="Exponential Shannon Index",
        xlab="Locality")
dev.off()

#### Simpson ####

D <- diversity(dataset.comm, index = "simpson") # Simpson Index
summary(D) #gives summary statistics for the plots
View(D)

tiff("SimpsonDiversity.tiff",width=2800,height=4000,units="px",compression="none",res=600)
boxplot(D ~ landuse$group,
        data = dataset.comm,
        col=(c("#92C5DE")),boxwex=0.5,
        ylab="Inverse Simpson Index",
        xlab="Locality")
dev.off()


# Richness-invSimpson-expShanon large graph ####
tiff("Richness-invSimpson-expShannon.tiff",width=5500,height=3000,units="px",compression="none",res=600)
par(mfrow=c(1,3))
boxplot(spnum ~ group,
        data = dataset.div,
        col=(c("#92C5DE")),boxwex=0.4,
        ylab="Ant Richness",
        xlab="Locality")
boxplot(D ~ landuse$group,
        data = dataset.comm,
        col=(c("#92C5DE")),boxwex=0.4,
        ylab="Inverse Simpson Index",
        xlab="Locality")
boxplot(eH ~ landuse$group,
        data = dataset.comm,
        col=(c("#92C5DE")),boxwex=0.4,
        ylab="Exponential Shannon Index",
        xlab="Locality")
dev.off()




### NMDS analysis ####

set.seed(3)
data.mds<-metaMDS(dataset.comm, distance = "bray", k = 3, maxit= 999, trymax = 500, autotransform =TRUE, noshare = 0.1, expand = TRUE, trace = 1, plot = FALSE)
data.mds

goodness(data.mds)

tiff("NMDS-stress-plot2.tiff",width=4500,height=3000,units="px",compression="none",res=600)
stressplot(data.mds)
dev.off()

plot(data.mds)
#in this plot the communities (“sites”) are open circles and species are red crosses

ordiplot(data.mds,type="n")
orditorp(data.mds,display="species",col="red",cex=0.75,air=0.01)
orditorp(data.mds,display="sites",cex=0.95,air=0.01)

orditorp(data.mds,display="species",col="darkgrey",air=0.01)
orditorp(data.mds,display="sites",col=c(rep("forestgreen",3),rep("turquoise3",3),rep("firebrick2",3), rep("slategrey",3), rep("orange1",3), rep("deepskyblue3",3), rep("darkslateblue",3), rep("maroon",3)),
         air=0.01,cex=1.5)
orditorp(data.mds,display="sites",
         air=0.01,cex=1.5)

# Anosim ####
# on 
# ratio between within-the-group and between-the-group variation
comm.ano <- anosim(dataset.comm, group, distance = "bray", permutations = 9999)
comm.ano
summary(comm.ano)
plot(comm.ano, notch=FALSE)


## final revised nmds plot ####
tiff("NMDS-final-rev.tiff",width=4000,height=2800,units="px",compression="none",res=600)
# make grey poligons that conect plots of same land use
treat=c(rep("BO",3),rep("SV",3), rep("TO",3), rep("ST",3), rep("ZR",3), rep("VI",3), rep("LE",3), rep("HR",3))
ordiplot(data.mds,type="n")
# plot just the plots, colour by locality, pch=19 means plot a circle
points(data.mds, "sites", pch = 19, cex=1.1, col = "#b10026", select = landuse$group == 
         "BO")
points(data.mds, "sites", pch = 19, cex=1.1, col = "#e31a1c", select = landuse$group == 
         "SV")
points(data.mds, "sites", pch = 19, cex=1.1, col = "#fc4e2a", select = landuse$group == 
         "TO")
points(data.mds, "sites", pch = 19, cex=1.1, col = "#fd8d3c", select = landuse$group == 
         "ST")
points(data.mds, "sites", pch = 19, cex=1.1, col = "#feb24c", select = landuse$group == 
         "ZR")
points(data.mds, "sites", pch = 19, cex=1.1, col = "#005a32", select = landuse$group == 
         "VI")
points(data.mds, "sites", pch = 19, cex=1.1, col = "#41ab5d", select = landuse$group == 
         "LE")
points(data.mds, "sites", pch = 19, cex=1.1, col = "#addd8e", select = landuse$group == 
         "HR")
legend(c("Stress: 0.087","ANOSIM: p=0.0001"),cex=0.7,
       bty="n",x="bottomright")
ordihull(data.mds,groups=treat,draw="polygon",col="grey",label=F)
dev.off()

