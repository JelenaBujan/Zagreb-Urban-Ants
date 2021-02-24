################################################################
####          Code for ZG MRAV data set analysis             ###
####  Written by multiple internet pages and Ana Ješovnik    ###
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
