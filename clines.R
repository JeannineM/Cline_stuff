# Load packages #######################################################################
library("fields")
library("geosphere")
library("ade4")
library("RColorBrewer")
library("gplots")
library("ggplot2")
library('ggbiplot')
library('MASS')
library('ape')
library('lme4')
library('scales')
library('maps')
library('raster')

# Functions ###########################################################################

# estimate clines per SNP
## input == genotype data format (0,1,2)
## regdat == x-axis of cline
## output is a table and a list of the sites which are significant

add_GLMclineLines=function(input, regdat){
  glmArray=c()
  #lowFreq=c()
  signSites = c()
  x <- regdat[order(regdat)]
  for (i in 11:ncol(input)){
    countA <- input[,i]
    countB <- 2 - countA
    
    GLM_func <- glm(cbind(countA, countB) ~ regdat, family = "quasibinomial", maxit=1000)
    
    # store non-NA GT-values
    GLM_dat <-data.frame(dist=regdat[!is.na(input[,i])],
                         SP=input$Sp2[!is.na(input[,i])], 
                         GT=(input[,i][!is.na(input[,i])]/2),
                         fiGT=fitted(GLM_func))
    slope = c()
    intercept = c()
    
    if (as.vector(GLM_func$coefficients[2]) < 0) {
      slope = -as.vector(GLM_func$coefficients[2])
      intercept = -as.vector(GLM_func$coefficients[1])
      } else {
      slope= as.vector(GLM_func$coefficients[2])
      intercept= as.vector(GLM_func$coefficients[1])
      }
    
    # save GLM results in array
    glmArray<-rbind(glmArray,
                    c(intercept=intercept,
                    slope=slope,
                    stErr=as.vector(sqrt(diag(vcov(GLM_func)))[2]),
                    centre=(-intercept / slope),
                    SLsign=as.logical(coef(summary(GLM_func))[2,4] < 0.05), 
                    Pmin=(1-(1/(1 + exp(intercept + slope*x))))[1],
                    Pmax=tail(1-(1/(1 + exp(intercept + slope*x))),1),
                    deltaP=tail(1-(1/(1 + exp(intercept + slope*x))),1) - (1-(1/(1 + exp(intercept + slope*x))))[1],
                    N=length(input$Sp2[!is.na(input[,i])]),
                    xAxis=regdat[!is.na(input[,i])],
                    SP=input$Sp2[!is.na(input[,i])], GT=(input[,i][!is.na(input[,i])]/2),
                    fiGT=fitted(GLM_func)))
    
    if (as.logical(coef(summary(GLM_func))[2,4] < 0.05) == TRUE & (-intercept / slope) > range(regdat)[1] & (-intercept / slope) < range(regdat)[2]) {
      signSites = c(signSites, i)
      if (mean(GLM_dat$fiGT[GLM_dat$SP == "ns"]) < 0.1 && mean(GLM_dat$fiGT[GLM_dat$SP == "hisp"]) < 0.1) {
          lowFreq=rbind(lowFreq, c(i, colnames(input[i])))}
    }
} #run the glm and store dat
  signSites <<- signSites -10
  #lowFreq <<- lowFreq
  return(glmArray)
}

## example call:
# glmArray2228_lat <- add_GLMclineLines(data[lat_order,], data$lat[lat_order]); signSites2228_lat <- signSites

# plot the estimated clines
## input == output of add_GLMclineLines (eg glmArray2228_lat)
## regdat == x-axis of cline
## colour == colour of the clines to be plotted

plot_GLMclines=function(input,regdat, color){
  
  for (i in 1:nrow(input)) {
  
  # only plot clines those centre are within the range of regdat
    if (input[i,5] == "TRUE" & as.numeric(input[i,4]) > range(regdat)[1] & as.numeric(input[i,4]) < range(regdat)[2]){
    
      b <- as.numeric(input[i, 1]) # intercept b
      a <- as.numeric(input[i, 2]) # slope m/a
      x <- regdat[order(regdat)]
      y <- 1-(1/(1 + exp(b + a*x)))
 
    points(x, y, col=color, t="l")
    }
  }
}
## example call: (eg all clines of a particular locus
# par(mfrow=c(1,1))
# plot(fitNSness304_new$dist, fitNSness304_new$admNS, xlab="distance in km", ylab="allele frequency", pch=1, cex=0.4, col="black")
# plot_GLMclines(input=glmArray2228HyS_dist[pos$CHROM == "LocusA" & signSitesDist$rowSum > 0,], regdat = distCentre304, color = "orange" )

# after GLM sift through the clines and plot
# input == output of add_GLMclineLines (eg glmArray2228_lat)
# regdat ==
# condition == any subset of glmArrayXXX
# outnameTab == name of file to store outout

selectSNPsbyCondition=function(input,regdat,condition,outnameTab){
  
  table=c()
  vector=c()
  locus=c()
  genes=c()
  
  # extract sites by condition
  for (i in which(condition==TRUE)){
    locus=c(locus,pos$CHROM_POS[i])
    genes=c(genes,pos$CHROM[i] )
      vector=c(gene=pos$CHROM[i], locus=pos$CHROM_POS[i],
             glmArray2228_dist[i,1:9], glmArray2228HyN_dist[i,1:9], glmArray2228HyS_dist[i,1:9],
             save_fis[,i], save_fis_sign[,i], sumTRUEFis=sum(save_fis_sign[,i], na.rm=TRUE),
             fst=save_fst[i], fst_sign=sign_fst[i])
    table=rbind(table, vector)
  }
  table<-as.data.frame(table)  
  rownames(table) <- locus
  table<<-table
  
  # output table to file
  write.csv(table, file = paste("~/",outnameTab,".csv", sep = ""))
  
  # summary gene table
  geneTab=c()
  for (i in levels(table$gene)){
    
    geneTab=rbind(geneTab, c(i, 
                             N_Snps = sum(pos$CHROM == i), 
                             N_signSnps = sum(pos$CHROM == i & condition == TRUE),
                             meanSlope=mean(as.numeric(input[pos$CHROM == i & condition == TRUE,2])),
                             sdSlope=sd(as.numeric(input[pos$CHROM == i & condition == TRUE,2])),
                             meanCentre=mean(as.numeric(input[pos$CHROM == i & condition == TRUE,4])),
                             sdCentre=sd(as.numeric(input[pos$CHROM == i & condition == TRUE,4])),
                             meanFst=mean(save_fst[pos$CHROM == i & condition == TRUE]),
                             sdFst=sd(save_fst[pos$CHROM == i & condition == TRUE])))
  }
  #
  geneTab <- as.data.frame(geneTab)
  colnames(geneTab) <- c("gene","N_Snps", "N_signSnps", "meanSlope", "sdSlope", "meanCentre", "sdCentre", "meanFst", "sdFst")
  
  geneTab <<- geneTab
  
  # plot the gene clines (all that are fitting the condition)
  par(mfrow=c(3,6))
  for (i in levels(table$gene)){
    plot(fitNSness304_new$dist, fitNSness304_new$admNS, ylab="non-scripta ness", pch=1, cex=0.4, col="gray31", xlab="")
    
    subtitle <- paste("N_Snps:",geneTab$N_signSnps[geneTab$gene == i],
                      "mS:",signif(as.numeric(as.vector(geneTab$meanSlope[geneTab$gene == i])),digits=4),
                      "mC:",signif(as.numeric(as.vector(geneTab$meanCentre[geneTab$gene == i])), digits=4))
    
    title(main=paste(i), xlab = subtitle)
    
    for (k in which(pos$CHROM == i & condition == TRUE)) {
      # according to input
      b <- as.numeric(input[k, 1]) # intercept b
      a <- as.numeric(input[k, 2]) # slope m/a
      x <- regdat[order(regdat)]
      y <- 1-(1/(1 + exp(b + a*x)))
 
    points(x, y, col="blue", t="l")
     }
    
    for (k in which(pos$CHROM == i & condition == TRUE & glmArray2228HyN_dist[,5] == TRUE)) {
      # HyN clines
      b <- as.numeric(glmArray2228HyN_dist[k, 1]) # intercept b
      a <- as.numeric(glmArray2228HyN_dist[k, 2]) # slope m/a
      x <- regdat[order(regdat)]
      y <- 1-(1/(1 + exp(b + a*x)))
 
    points(x, y, col="green4", t="l")
     }
    
    for (k in which(pos$CHROM == i & condition == TRUE & glmArray2228HyS_dist[,5] == TRUE)) {
      # HyS clines
      b <- as.numeric(glmArray2228HyS_dist[k, 1]) # intercept b
      a <- as.numeric(glmArray2228HyS_dist[k, 2]) # slope m/a
      x <- regdat[order(regdat)]
      y <- 1-(1/(1 + exp(b + a*x)))
 
    points(x, y, col="orange", t="l")
     }
    
  }
    
}
## example call:
# selectSNPsbyCondition(glmArray2228_dist, distCentre304, condition = signSitesDist$rowSum > 1 & as.numeric(input[,2]) < mean(as.numeric(input[signSites2228_dist,2])) & as.numeric(input[,4]) > -15 & as.numeric(input[,4]) < 15, "pups"); geneTab_cond3 <- geneTab; table_cond3 <- table



# Analyses ############################################################################
# Plot admixture bar plot
## with dynamic separator positions

## vector to order samples
lat_orderPop <- order(data_pop3$Sp2=="ns", data_pop3$Sp2=="hyN", data_pop3$Sp2=="hyS", data_pop3$Sp2=="hisp", data_pop3$lat, data_pop3$pop.ID, decreasing = TRUE)
# IDadm_orderPop <- order(data_pop3$Sp2=="ns", data_pop3$Sp2=="hyN", data_pop3$Sp2=="hyS", data_pop3$Sp2=="hisp",data_pop3$pop.ID, decreasing = TRUE)

xxxB=data_pop3$N[lat_orderPop] 
xxxA=data_pop3$pop.ID[lat_orderPop]

xxxx=0
sep=c()
for (i in 1:dim(data_pop3)[1]){
  xxxx= xxxx+ as.integer(xxxB)[i]
  sep=c(sep, xxxx)
}

## structure barplot for K=2 # sorted by pop and within as latitude
# take the co-ancestry table of structure like output

barplot(as.matrix(k2.meanQ), col=c("blue", "red"),
        border=TRUE, main="K = 2 run 7 * best",  
        names.arg = rep("",nrow(data)),
        space=0, width=1, las=2,
        axes=FALSE)

axis(side = 2)
at_sep = sep - 1.3
axis(side = 1, at= at_sep, labels = xxxA, tick = FALSE, line = FALSE, las=2, cex=0.1)

abline(v=c(sep), col="white")
abline(v=cumsum(c(sum(data$Sp2[lat_order] == "ns"), sum(data$Sp2[lat_order] == "hyN"), sum(data$Sp2[lat_order] == "hyS"), sum(data$Sp2[lat_order] == "hisp"))),col="black", lwd=3)

# Regress admixture (and organelle GT) proportions on latitude and extract info 
countA <- k2.meanQ$admNS 
countB <- 1- countA
NuclAdmGLM <- glm(cbind(countA, countB) ~ data$lat, family = "quasibinomial")
   
summary(NuclAdmGLM)
  I_Adm <- NuclAdmGLM$coefficients[1] # intercept b
  S_Adm <- NuclAdmGLM$coefficients[2] # slope m=a
  StdErr_Adm <- sqrt(diag(vcov(NuclAdmGLM)))[2]
  c_Adm <- -I_Adm/S_Adm
  CI_Adm <- confint(NuclAdmGLM)

  # -> values confidence intervals
  pr <- predict(NuclAdmGLM, se.fit = TRUE)
  family <- family(NuclAdmGLM)
  Adm_lower <- family$linkinv(pr$fit - qnorm(0.95) * pr$se.fit)
  Adm_upper <- family$linkinv(pr$fit + qnorm(0.95) * pr$se.fit)

# cp genotype per individual along latitude
countA <- 1-data$cp.GT
countB <- data$cp.GT
CPF_GLM <- glm(cbind(countA, countB) ~ data$lat, family = "binomial")

# Plot admixture cline along latitude

fitNSness309 <- data.frame(data[1:309, 1:10], glmFit=fitted(NuclAdmGLM), glmCPFFit=fitted(CPF_GLM))
lat_order <- order(data$lat)

colorFid <- popinfo$Sp2[order(popinfo$Sp2=="ns", popinfo$Sp2=="hyN", popinfo$Sp2=="hyS", popinfo$Sp2=="hisp", popinfo$lat, decreasing = TRUE)]
colorFid <- gsub("hisp", "red", colorFid)
colorFid <- gsub("ns", "blue", colorFid)
colorFid <- gsub("hyN", "green", colorFid)
colorFid <- gsub("hyS", "orange", colorFid)
 
plot(k2.meanQ.m[1,]~ popinfo$lat[order(popinfo$Sp2=="ns", popinfo$Sp2=="hyN", popinfo$Sp2=="hyS", popinfo$Sp2=="hisp", popinfo$lat, popinfo$sample.ID, decreasing = TRUE)], xlim=c(41.8,42.9), xlab="latitude", ylab="admixture proportions", col=colorFid, pch=20, cex=1.2)

title(main="Nuclear (violet, 95% CI) and Organelle (pink, 95% CI) cline and individual admixture proportion of Nonscripta-ness (dots)")
lines(fitNSness309$lat[lat_order], fitNSness309$glmCPFFit[lat_order], lty=2, lwd=3, col="deeppink2")
points(Adm_lower[lat_order] ~ fitNSness309$lat[lat_order], type="l", lty=2, lwd=1, col="blueviolet")
points(Adm_upper[lat_order] ~ fitNSness309$lat[lat_order], type="l", lty=2, lwd=1, col="blueviolet")

hist(k2.meanQ$admNS, breaks=11, ylim=c(0,110), main="Histogram of Non-scriptaness", xlab="admixture proportions across the whole hybrid zone")

# Geographic cline (obtain shortest straight distance between location and a 2D cline)

# prepare map
srtm_35_04 <- getData(name='SRTM', download = TRUE,lon=-7, lat=42)
srtm_36_04 <- getData(name='SRTM', download = TRUE,lon=-4, lat=42)
srtm <- mosaic(srtm_35_04, srtm_36_04, fun=mean)

# estimating shortest distance
cline05 <- read.csv("~/clineOFdot5.csv", h=FALSE) # required list of lon/lat coordiantes of a 2D cline
colnames(cline05) <- c("lon", "lat")

points304 <- cbind(data$lon, data$lat)

distinM304 <- as.data.frame(cbind((dist2Line(points304, cline05[47:265,], distfun=distHaversine)), points304))
rownames(distinM304) <- rownames(data304)
colnames(distinM304) <- c("Dist_m", "L_lon", "L_lat", "Lon", "Lat")

# plot distance and order

par(mfrow=c(1,2))
plot(srtm, xlim=c(-7.6, -5.6), ylim=c(41.75,43.1),col=rgb(colorRamp(c("chartreuse3","darkorange4", "darkgrey","whitesmoke"))((1:200)/200)/255)[20:200])
points(points304, xlab="Latitude", ylab="Longitude", bg=colorFid, pch=21) 
points(cline05[47:265,], type="l", col="gray31", lwd=3)

points(distinM304[,2], distinM304[,3], col='red', pch='x')
for (i in 1:nrow(distinM304)) lines(gcIntermediate(points304[i,], distinM304[i,2:3], 10), lwd=1, lty=2, cex=0.8)

# Line_dot
xL <- distinM304[,"L_lon"] # Longitude
yL <- distinM304[,"L_lat"] # Latitude

# Pop_dot
xP <- distinM304[,"Lon"] # Longitude
yP <- distinM304[,"Lat"] # Latitude

distinM_new <-c()
flipped <- c()
normal <- c()
for (i in 1:nrow(distinM304)) {
  if      (yP[i] < yL[i]) {
  distinM_new=rbind(distinM_new, -distinM304[i,1])
  flipped=rbind(flipped, i)}
  else if (yP[i] < yL[i] & xP[i] > xL[i]) {
    distinM_new=rbind(distinM_new, -distinM304[i,1])
    flipped=rbind(flipped, i)}
  else {
    distinM_new=rbind(distinM_new, distinM304[i,1])
    normal=rbind(normal, i)}
}

plot(distinM_new/1000, k2.meanQ[,1], col="black", xlab="distance in km from cline centre (adm = 0.5)", pch=21, bg=colorFid, ylab="non-scripta ness", cex=0.9)

distCentre304 <- c(distinM_new/1000)
distOrd304 <- order(distCentre304)

