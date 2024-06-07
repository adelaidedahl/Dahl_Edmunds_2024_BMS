#adelaide dahl
#processing April 2022 Moorea data
#made April 19, 2023

library(purrr)
library(BiodiversityR)
library(labdsv)
library(vegan)
library(readxl)
library(betapart)
library(tidyverse)
library(car)
library(ggpubr)
library(xlsx)
library(ggplot2)
library(dplyr) #plots
library(writexl)
library(lubridate)
library("plotrix")
library(nationalparkcolors)
library(lme4)
library(wesanderson)
library(plotly) #violin plots/interactive plot options
library(circular) #Watson's U2 test
library(reshape2) #stack and transform df to vertical (melt function)
library(RColorBrewer) #color setup
library(scales)
library(lmerTest) #add pvalue to lmer
library(here) #here fxn for outputs
library(ggdist) #rainclouds
library(gghalves) #rainclodus
library(patchwork)

setwd("~/Desktop/Mo'orea/2022")

#3 major hypotheses: 
##H1corals are distributed asymmetrically on bommies relative to flow direction
##H2 the growth rate of corals varies around perimeter of bommies
##H3 integrated flow speed differs around the perimeter of bommies

#H1#####
#testing H1 (corals are distributed asymmetrically on bommies relative to flow direction --> treatment)
#DV: corals per m^2
#IV: quadrant (relative to flow direction), site


#to get distribution relative to surface area of the bommie
#size of bommie corrected for SA of bommie as 1/2 cylinder
SAbeta = read.csv("Data/ALL_genera.xlsx - bommies.csv", header = TRUE)
#Data/ALL_genera.xlsx - bommies.csv
View(SAbeta)
SAbeta$bommie.ID <- as.factor(SAbeta$bommie.ID)

#test significant relationship between angle of the bommies primary axes####
#need to reduce df to just bommie angles and relative reef crest measurements
BommieANGLE <- SAbeta[-c(3, 10, 12:22)]
BommieANGLE <- na.omit(BommieANGLE)

#since there are 4 rows associated w every bommie right now we nee dto reduce this to one row each
#averages will let us do that since it will output the correct number
BommieANGLEsummary <- BommieANGLE %>%
  group_by(bommie, site) %>%
  summarize(meanangle=mean(angle, na.rm=TRUE))
BommieANGLEsummary$site <- as.factor(BommieANGLEsummary$site)

factor_to_value <- function(site) {
  # Define the mapping logic based on your specific requirements
  if (site == "West") {
    return(0)
  } else if (site == "MRB") {
    return(10)
  } else if (site == "East") {
    return(0)
  }
}

# Use mutate() to add a new column based on the factors column
BommieANGLEsummary <- BommieANGLEsummary %>%
  mutate(crestangle = factor_to_value(site))
view(BommieANGLEsummary)
#great now we have our df in the right shape

#meand and se of bommie angles####
mean(BommieANGLEsummary$meanangle)
sd(BommieANGLEsummary$meanangle, na.rm=TRUE)/sqrt(length(na.omit(BommieANGLEsummary$meanangle)))


#Watson's U2 test to see if angles differ relative to crest#####
meanvector <- c(BommieANGLEsummary$meanangle)
crestvector <- c(BommieANGLEsummary$crestangle)
meanvector_degree <-  circular(meanvector, units = "degrees", template = "geographics")
crestvector_degree <-  circular(crestvector, units = "degrees", template = "geographics")

watson.test(meanvector_degree)

####

#getting summary data per site and quadrant#### 
allsitemeanSURVEY2 <- SAbeta %>%
  group_by(site, quadrant) %>%
  summarize(meantotal=mean(Corals.per.0.25m2, na.rm=TRUE), se=sd(Corals.per.0.25m2, na.rm=TRUE)/sqrt(length(na.omit(Corals.per.0.25m2))))
view(allsitemeanSURVEY2)

#trying to stack and transform species columns####
SAbetaSIMPLE <- SAbeta[,-21:-22]

#remove rows with no values
SAbetaSIMPLEnull <- SAbetaSIMPLE[rowSums(SAbetaSIMPLE[, 15:20] == 0) < 6,] %>%
  mutate(treatment=as.factor(treatment)) %>%
  mutate(site=as.factor(site)) %>%
  mutate(bommie.ID=as.factor(bommie.ID))

view(SAbetaSIMPLEnull)

SAbetaSTACK<- melt(SAbetaSIMPLEnull, id=c("site", "date", "transect", "length", "width", "height", "SA.cm.2.", "SA.m.2.", "bommie", "bommie.ID", "angle", "quadrant", "treatment", "substrate"), variable.name = "species") %>%
  mutate(treatment=as.factor(treatment)) %>%
  mutate(site=as.factor(site)) %>%
  mutate(bommie.ID=as.factor(bommie.ID))
###LETS GO

#need to remove zero values again
SAbetaSTACK <- SAbetaSTACK[SAbetaSTACK$value != 0, ]

#SUMMARY INFO####
#want summary info for averages
speciesSUMMARY <- SAbetaSTACK %>%
  group_by(species, treatment) %>%
  summarize(meantotal=mean(value, na.rm=TRUE), se=sd(value, na.rm=TRUE)/sqrt(length(na.omit(value))), n = n())

byBommieSummary <- SAbetaSTACK %>%
  group_by(bommie.ID) %>%
  summarize(meantotal=mean(value, na.rm=TRUE), se=sd(value, na.rm=TRUE)/sqrt(length(na.omit(value))))

byQuadrantSummary <- SAbetaSTACK %>%
  group_by(treatment) %>%
  summarize(meantotal=mean(value, na.rm=TRUE), se=sd(value, na.rm=TRUE)/sqrt(length(na.omit(value))))

#want totals
speciesTOTALS <- SAbetaSTACK %>%
  group_by(species) %>%
  summarize(sum=sum(value, na.rm=TRUE))
#add percent column
speciesTOTALS$percent <- (speciesTOTALS$sum/(sum(speciesTOTALS$sum)))*100

bommieTOTALS <- SAbetaSTACK %>%
  group_by(bommie.ID) %>%
  summarize(sum=sum(value, na.rm=TRUE))
mean(bommieTOTALS$sum) #mean of total per bommie
sd(bommieTOTALS$sum, na.rm=TRUE)/sqrt(length(na.omit(bommieTOTALS$sum))) #se of total per bommie
#####

#removing juveniles
speciesadults<-speciesSUMMARY[!(speciesSUMMARY$species=="UJ"),]
speciesadults

speciesadults <- speciesadults %>%
  mutate(quadrant = case_when(treatment == "downstream" ~ 2,
                              treatment == "upstream" ~ 1,
                              treatment == "left" ~ 3,
                              treatment == "right" ~ 4))
                            


#####

#####VISUALIZE H1####
#bar plot with counts of species (y) by treatment (fill) and site (x) 
coralpal <- c("#0092B7","#78B7C5", "#EBCC2A","#E17300", "#B31D0B")
  #c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00")


show_col(coralpal)

#plotting coral species with a boxplotAP####
adultSPplot <-ggplot(speciesadults, aes(x=quadrant, y=meantotal, fill = species)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values=coralpal) +
  geom_errorbar(aes(ymax=speciesadults$meantotal+speciesadults$se, ymin=speciesadults$meantotal-speciesadults$se), position = position_dodge(0.90), width = 0.25) +
  #labs(title="Distribution of Corals Across Sites in the Backreef", x="location relative to reef crest",y = "mean coral abundance by m^2") +
  theme_classic()+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 5)),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 14, angle = 0,
                                   hjust = 0.5, vjust = 0.5,margin=margin(0,0,10,0)),
       axis.text.y = element_text(colour = "black", size = 14))

adultSPplot

ggsave(here("Outputs", "speciesSUM_1.pdf"), 
       dpi = 300,
       width = 15,
       height = 10)
#####
#violin plot by treatment and site####
wes_palette("Zissou1")
gp2 <- wes_palettes$Zissou1

fig <- SAbeta %>%
  plot_ly(type = 'violin') 
fig <- fig %>%
  add_trace(
    x = ~site[SAbeta$treatment == 'upstream'],
    y = ~Corals.per.0.25m2[SAbeta$treatment == 'upstream'],
    legendgroup = 'upstream',
    scalegroup = 'upstream',
    name = 'upstream',
    box = list(
      visible = T
    ),
    meanline = list(
      visible = T
    ),
    color = I("#005369")
  ) 
fig <- fig %>%
  add_trace(
    x = ~site[SAbeta$treatment == 'downstream'],
    y = ~Corals.per.0.25m2[SAbeta$treatment == 'downstream'],
    legendgroup = 'downstream',
    scalegroup = 'downstream',
    name = 'downstream',
    box = list(
      visible = T
    ),
    meanline = list(
      visible = T
    ),
    color = I("#78B7C5")
  )  
fig <- fig %>%
  add_trace(
    x = ~site[SAbeta$treatment == 'left'],
    y = ~Corals.per.0.25m2[SAbeta$treatment == 'left'],
    legendgroup = 'left',
    scalegroup = 'left',
    name = 'left',
    box = list(
      visible = T
    ),
    meanline = list(
      visible = T
    ),
    color = I("#C70000")
  )
fig <- fig %>%
  add_trace(
    x = ~site[SAbeta$treatment == 'right'],
    y = ~Corals.per.0.25m2[SAbeta$treatment == 'right'],
    legendgroup = 'right',
    scalegroup = 'right',
    name = 'right',
    box = list(
      visible = T
    ),
    meanline = list(
      visible = T
    ),
    color = I("#FF8B8B")
  )

fig <- fig %>%
  layout(
    yaxis = list(
      zeroline = F
    ),
    violinmode = 'group'
  )

fig
#wow this is so cool 
#####
#
#want to show species distribution across sites
#in raincloud

#####ANALYSIS H1####
#SAbetaSTACK is the df you want 

#BLOCK ANOVA looking at site and treatment as df and species as iv with bommie as the block
#remove zero counts
SAbetaSTACKnull <- SAbetaSTACK[SAbetaSTACK$value != 0, ]

#getting sample size for figure
nvaluesFIG3 <- 
  SAbetaSTACKnull %>%
  group_by(treatment, species) %>%
  summarize(num = n(),
            totalSP = sum(value))

write_xlsx(nvaluesFIG3,"Outputs/nvaluesFIG3.xlsx")

modelh1<-aov(value~bommie.ID + site*treatment, data=SAbetaSTACKnull)
summary(modelh1)
#Don't forget to check your model assumptions! --> not normal so were running a PERMANOVA

distributionPERM <-adonis2(SAbetaSIMPLEnull[,-1:-14]~treatment + site, 
                           data = SAbetaSIMPLEnull, 
                           strata = SAbetaSIMPLEnull$bommie.ID, 
                           permutations = 999, 
                           method="bray")

distributionPERM

#ratio of length to width and associated se/mean####
ratiodf <- SAbetaSIMPLEnull$length/SAbetaSIMPLEnull$width
ratiodf <- data.frame(ratiodf)
mean(ratiodf$ratiodf)
sd(ratiodf$ratiodf, na.rm=TRUE)/sqrt(length(na.omit(ratiodf$ratiodf))) #se of ratio
####

####testing H2######
#DV: change in BW
#IV: quadrant (relative to flow direction), site

#read GR DF####
GRbeta <- read.csv("Data/Buoyant Weights - Updated 2023.csv", header = TRUE)
bommieGR = GRbeta[which(GRbeta$location == "bommie"),] %>%
  na.omit(bommieGR) %>%
  mutate(quadrant = as.factor(quadrant)) %>%
  mutate(GR..g.day. = as.numeric(GR..g.day.)) %>%
  mutate(GR..mg.day. = as.numeric(GR..mg.day.)) %>%
  mutate(site = as.factor(site)) %>%
  mutate(bommie.ID = as.factor(bommie.ID)) 


bommieGRsubset <- subset(bommieGR, GR..mg.day.>-1) 


#summary data for GR df
GRsummary <- bommieGRsubset %>%
  #group_by(quadrant) %>%
  summarize(meantotal=mean(GR..mg.day., na.rm=TRUE), se=sd(GR..mg.day., na.rm=TRUE)/sqrt(length(na.omit(GR..mg.day.))))

range(bommieGRsubset$GR..mg.day.)
view(GRsummary)

#####

#####VISUALIZE H2####
treatmentPal <- c("#005369", "#78B7C5","#C70000", "#FF8B8B")
#use bommieGRsubset

###### 
GRfig <- bommieGRsubset %>%
  plot_ly(
    x = ~quadrant,
    y = ~GR..mg.day.,
    split = ~quadrant,
    type = 'violin',
    color = ~quadrant,
    colors = c("#005369", "#78B7C5","#C70000", "#FF8B8B"),
    box = list(
      visible = T
    ),
    meanline = list(
      visible = T
    )
  ) 

GRfig <- GRfig %>%
  layout(
    xaxis = list(
      title = "Treatment"
    ),
    yaxis = list(
      title = "nubbin growth mg/day",
      zeroline = F
    )
  )
GRfig

#raincloud plot
#####
GRrain <-
  ggplot(bommieGRsubset, aes(x = quadrant, y = GR..mg.day., fill = quadrant)) +
  theme_classic()+
  stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.2, 
    point_colour = NA) + 
  geom_boxplot(
    width = .1, 
    outlier.shape = NA) +
  stat_dots(aes(color = quadrant), 
            side = "left", 
            shape = 19,
            scale = .1,
            justification = 1.7) +
  scale_fill_manual(values = c("#006782", "#78B7C5","#CA3A3A", "#FFADAD"))+
  scale_color_manual(values = c("#006782", "#78B7C5","#CA3A3A", "#FFADAD"))

GRrain
#####


#####ANALYSIS H2#####
growrate <- bommieGRsubset$GR..mg.day.
hist(growrate,
     xlab = "Growth Rate of Coral Nubbins (mg cm^-2 day^-1)",
     main = "Distribution of Growth Rate in Pocillopora Coral Nubbins",
     xlim=c(-3,3),
     freq=FALSE)
qqp(bommieGRsubset$GR..mg.day., "norm")

bommieGRsubset$log <- log(bommieGRsubset$GR..mg.day. +1)

#test if treatment signifincatly impacted grwoth rate of coral nubbins during the treatment
GRmodel1<-lmer(log~site*quadrant + (1|bommie.ID), data=bommieGRsubset)
Anova(GRmodel1)
plot(GRmodel1)
#no signifincat effect of treatment or site
#assumptions look good

########testing H3 #####
#DV: dissolution of clod
#IV: quadrant (relative to flow direction), site 

#clod df reading####
clodbetafile = read.csv("Data/ALL_CLODS - Sheet1.csv", header = TRUE) %>%
  mutate(treatment = as.factor(treatment)) %>%
  mutate(site = as.factor(site)) %>%
  mutate(bommie.ID = as.factor(bommie.ID))

clodbetafile$DIlog <- log(clodbetafile$rate.dissolution+1)
qqp(clodbetafile$DIlog, "norm")

bommieclods = clodbetafile[which(clodbetafile$location == "bommie"),]
bommieclods$quadrant <- as.factor(bommieclods$quadrant)

#get sample sizes for each side
upstreamclod = bommieclods[which(bommieclods$quadrant == "1"),]
downstreamclod = bommieclods[which(bommieclods$quadrant == "2"),]
westclod = bommieclods[which(bommieclods$quadrant == "3"),]
eastclod = bommieclods[which(bommieclods$quadrant == "4"),]

CLODsummary <- bommieclods %>%
  group_by(site) %>%
  summarize(meantotal=mean(rate.dissolution, na.rm=TRUE), se=sd(rate.dissolution, na.rm=TRUE)/sqrt(length(na.omit(rate.dissolution))))
view(CLODsummary)

######VISUALIZE H3#####
DIfig <- bommieclods %>%
  plot_ly(
    x = ~treatment,
    y = ~rate.dissolution,
    split = ~treatment,
    type = 'violin',
    color = ~treatment,
    colors = c("#005369", "#78B7C5","#C70000", "#FF8B8B"),
    box = list(
      visible = T
    ),
    meanline = list(
      visible = T
    )
  ) 

DIfig <- DIfig %>%
  layout(
    xaxis = list(
      title = "Side of Bommie"
    ),
    yaxis = list(
      title = "Rate of Dissolution (g/hr)",
      zeroline = F
    )
  )
DIfig
#####
CLODrain <-
  ggplot(bommieclods, aes(x = quadrant, y = rate.dissolution, fill = quadrant)) +
  theme_classic()+
  stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.2, 
    point_colour = NA) + 
  geom_boxplot(
    width = .1, 
    outlier.shape = NA) +
  stat_dots(aes(color = quadrant), 
            side = "left", 
            shape = 19,
            scale = .1,
            justification = 1.7) +
  scale_fill_manual(values = c("#006782", "#78B7C5","#CA3A3A", "#FFADAD"))+
  scale_color_manual(values = c("#006782", "#78B7C5","#CA3A3A", "#FFADAD"))
CLODrain

GRrain / CLODrain

ggsave(here("Outputs", "Clod_GR.pdf"), 
       dpi = 300,
       width = 10,
       height = 15)

#####ANALYSIS H3#####

modelCLODS<-lmer(DIlog~treatment*site + (1|bommie.ID), data=bommieclods)

#Don't forget to check your model assumptions!
plot(modelCLODS)
#Looks good

#And then to get the ANOVA table:
anova(modelCLODS)
#no significance of site or treatment on clods

########EDITS
#need df with per bommie length and width 
library(dplyr)
summary_SA <- SAbeta %>%
  group_by(bommie.ID, length, width, height, site) %>%
  summarise_all(.funs = first) %>%
  mutate(length = as.numeric(length)) %>%
  mutate(width = as.numeric(width)) %>%
  mutate(height = as.numeric(height))

summary_SA2 <- summary_SA[, -5:-22]

frequency_SA <- melt(summary_SA2, id=c("bommie.ID"), variable.name = "dimension") %>%
  mutate(value = as.numeric(value)) %>%
  mutate_if(is.numeric , replace_na, replace = 0)

sizeplot <- ggplot(SAbeta, aes(x = length, y = width)) + 
  geom_point() + 
  theme_classic() + 
  labs(x = "bommie length (cm)",
    y = "bommie width (cm)")
sizeplot

#frequency
sizefrequencyPLOT <-
  ggplot(frequency_SA, aes(x = value)) +
  geom_histogram(aes(color = dimension, fill = dimension),
                 alpha = 0.4, position = "identity") +
  labs(x = "Value (cm)", y = "Frequency", title = "Frequency of dimension values on bommies") +
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#D8735B")) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#D8735B"))


sizefrequencyPLOT

#ggsave(ggsave(here("APCRS_fig3_frequency.png"), width = 10, height = 10, dpi = 300))


