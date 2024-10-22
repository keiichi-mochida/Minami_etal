#To cite R/qtl in publications, use
#Broman KW, Wu H, Sen S, Churchill GA (2003) R/qtl: QTL mapping in experimental crosses. Bioinformatics, 19: 889-890

setwd("C:/Users/Yoshi-RIKEN/Documents/R/Brachy analysis/2016_06_15 Bd21 x Koz4 F2 QTL 1st")

options(digits=4)
library(qtl)

mapthis_F2_ <- read.cross("csvr",
dir="C:/Users/yoshi/Documents/R/Brachy analysis/2016_06_15 Bd21 x Koz4 F2 QTL 1st",
file="Koz_1st_r.csv", estimate.map=FALSE)

summary(mapthis_F2_)
#individual 175, 175*0.9=157.5
#markers 442, 442*0.9=397.8
mapthis_F2_ <- est.rf(mapthis_F2_)
plotRF(mapthis_F2_, alternate.chrid=TRUE)

#Genetic map construction
#genotype data 欠損の確認
par(mfrow=c(1,1), las=1)
plotMissing(mapthis_F2_, main="missing genotype data")

par(mfrow=c(1,2), las=1)
plot(ntyped(mapthis_F2_, what="ind"), ylab="No. typed markers",
main="0 No. genotypes by individual", ylim=c(0,500))
plot(ntyped(mapthis_F2_, what="mar"), ylab="No. typed individuals",
main="0 No. genotypes by marker", ylim=c(0,200))

#makers genotyped: > 90% of population
mapthis_F2_1 <- subset(mapthis_F2_, ind=(ntyped(mapthis_F2_)>397.8))#397.8

nt.bymar <- ntyped(mapthis_F2_1, "mar")
todrop1 <- names(nt.bymar[nt.bymar < 157.5])#157.5
mapthis_F2_2 <- drop.markers(mapthis_F2_1, todrop1)

par(mfrow=c(3,2), las=1)

plot(ntyped(mapthis_F2_, what="ind"), ylab="No. typed markers"
, main="0 No. genotypes by individual", ylim=c(0,500))
plot(ntyped(mapthis_F2_, what="mar"), ylab="No. typed individuals"
, main="0 No. genotypes by marker", ylim=c(0,200))

plot(ntyped(mapthis_F2_1, what="ind"), ylab="No. typed markers"
, main="1 No. genotypes by individual", ylim=c(0,500))
plot(ntyped(mapthis_F2_1, what="mar"), ylab="No. typed individuals"
, main="1 No. genotypes by marker", ylim=c(0,200))

plot(ntyped(mapthis_F2_2, what="ind"), ylab="No. typed markers"
, main="2 No. genotypes by individual", ylim=c(0,500))
plot(ntyped(mapthis_F2_2, what="mar"), ylab="No. typed individuals"
, main="2 No. genotypes by marker", ylim=c(0,200))

summary(mapthis_F2_)
summary(mapthis_F2_1)
summary(mapthis_F2_2)

par(mfrow=c(1,1), las=1)
plotMissing(mapthis_F2_2, main="missing genotype data")

par(mfrow=c(3,1), las=1)
mapthis_F2_ <- est.rf(mapthis_F2_)
mapthis_F2_1 <- est.rf(mapthis_F2_1)
mapthis_F2_2 <- est.rf(mapthis_F2_2)
plotRF(mapthis_F2_, alternate.chrid=TRUE, main="F2_")
plotRF(mapthis_F2_, alternate.chrid=TRUE, main="F2_1")
plotRF(mapthis_F2_2, alternate.chrid=TRUE, main="F2_2")

#segregation
mapthis_F2_3 <- mapthis_F2_2
gt <- geno.table(mapthis_F2_3)
gt[gt$P.value < 0.05/totmar(mapthis_F2_3),]
todrop3 <- rownames(gt[gt$P.value < 1e-7,])
mapthis_F2_4 <- drop.markers(mapthis_F2_3, todrop3)
summary(mapthis_F2_4)

par(mfrow=c(1,1), las=1)
plotRF(mapthis_F2_4, alternate.chrid=TRUE)

#alelle freq.
g <- pull.geno(mapthis_F2_4)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq) / colSums(gfreq))
par(mfrow=c(1,3), las=1)
for(i in 1:3)
plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "AB", "BB")[i], ylim=c(0,1))

mapthis_F2_5 <- est.map(mapthis_F2_4)
par(mfrow=c(1,1), las=1)
plot.map(mapthis_F2_5)
summary(mapthis_F2_5)

mapthis_F2_6 <- replace.map(mapthis_F2_4, mapthis_F2_5)
plot.map(mapthis_F2_6)
summary(mapthis_F2_6)

rf6 <- pull.rf(mapthis_F2_6)
lod6 <- pull.rf(mapthis_F2_6, what="lod")
mn1 <- markernames(mapthis_F2_6, chr=1)
mn2 <- markernames(mapthis_F2_6, chr=2)
mn3 <- markernames(mapthis_F2_6, chr=3)
mn4 <- markernames(mapthis_F2_6, chr=4)
mn5 <- markernames(mapthis_F2_6, chr=5)

mn1
mn2
mn3
mn4
mn5

##marker ordering
pull.map(mapthis_F2_6, chr=5)
pull.map(mapthis_F2_6, chr=4)
pull.map(mapthis_F2_6, chr=3)
pull.map(mapthis_F2_6, chr=2)
pull.map(mapthis_F2_6, chr=1)
############################################################
summaryMap(mapthis_F2_6)

#plotMap(mapthis_F2_6, show.marker.names=TRUE)
plotMap(mapthis_F2_6, show.marker.names=FALSE)
plotRF(mapthis_F2_6)

#marker check
dropone <- droponemarker(mapthis_F2_6, error.prob=0.005)
par(mfrow=c(2,1))
plot(dropone, lod=1, ylim=c(-100,0))
plot(dropone, lod=2, ylab="Change in chromosome length")
summary(dropone, lod.column=2)

badmar <- rownames(summary(dropone, lod.column=2))[1:3]
mapthis_F2_7 <- drop.markers(mapthis_F2_6, badmar)

newmap1 <- est.map(mapthis_F2_7, error.prob=0.005)
mapthis_F2_8 <- replace.map(mapthis_F2_7, newmap1)

summaryMap(mapthis_F2_8)


#line check
par(mfrow=c(1,1), las=1)
plot(countXO(mapthis_F2_8), ylab="Number of crossovers")
mapthis_F2_9 <- subset(mapthis_F2_8, ind=(countXO(mapthis_F2_8) < 50))

newmap2 <- est.map(mapthis_F2_9, error.prob=0.005)
mapthis_F2_10 <- replace.map(mapthis_F2_9, newmap2)
summaryMap(mapthis_F2_10)

par(mfrow=c(2,1))
plot(countXO(mapthis_F2_8), ylab="Number of crossovers", main = "mapthis_F2_8", ylim=c(0,80))
plot(countXO(mapthis_F2_10), ylab="Number of crossovers", main = "mapthis_F2_10", ylim=c(0,80))

#Evidence for segregation distortion
gt <- geno.table(mapthis_F2_10, scanone.output=TRUE)
par(mfrow=c(2,1))
plot(gt, ylab=expression(paste(-log[10], " P-value")))
plot(gt, lod=3:5, ylab="Genotype frequency")
abline(h=c(0.25, 0.5), lty=2, col="gray")

mapthis_F2_10 <- est.rf(mapthis_F2_10)
par(mfrow=c(1,1))
plotRF(mapthis_F2_10, alternate.chrid=TRUE, main="F2_10")
summary(mapthis_F2_10)

pull.map(mapthis_F2_10, chr=1)
pull.map(mapthis_F2_10, chr=2)
pull.map(mapthis_F2_10, chr=3)
pull.map(mapthis_F2_10, chr=4)
pull.map(mapthis_F2_10, chr=5)

par(mfrow=c(1,1))
#plotMap(mapthis_F2_10, show.marker.names=TRUE)
plotMap(mapthis_F2_10, show.marker.names=FALSE)

###
#QTL analysis
###
qtl1 <- mapthis_F2_10
qtl1
geno.image(qtl1)

#imputation
fill <- fill.geno(qtl1)
geno.image(fill)

# autocofactor
autocofactors <- mqmautocofactors(fill, 100)
mqmplot.cofactors(fill, autocofactors, justdots=TRUE)

mqm_auto <- mqmscan(fill, autocofactors, pheno.col=c(1:3),　verbose = TRUE,　n.cluster=12)
write.csv(mqm_auto,  file = "mqm_auto.csv")

mqm_all <- mqmscanall(fill, cofactors = autocofactors, n.cluster=12)
pdf(file="mqmscanall(red=FW, green=nleaf blue=water).pdf")
mqmplot.multitrait(mqm_all, type="lines")
dev.off()
pdf(file="mqmscanall_chr=4(red=FW, green=nleaf blue=water).pdf")
mqmplot.multitrait(mqm_all, chr=4, type="lines")
dev.off()

mqm_auto.1 <- mqmscan(fill, autocofactors,
pheno.col= 1,　verbose = TRUE,　n.cluster=12)
mqm_auto.2 <- mqmscan(fill, autocofactors,
pheno.col= 2,　verbose = TRUE,　n.cluster=12)
mqm_auto.3 <- mqmscan(fill, autocofactors,
pheno.col= 3,　verbose = TRUE,　n.cluster=12)

# graph
plot(mqm_auto.1, ylim=range(c(0, 20)))
plot(mqm_auto.2, ylim=range(c(0, 20)))
plot(mqm_auto.3, ylim=range(c(0, 20)))

pdf(file="LOD.pdf")
plot(mqm_auto.1, ylim=range(c(0, 20)))
plot(mqm_auto.2, ylim=range(c(0, 20)))
plot(mqm_auto.3, ylim=range(c(0, 20)))
dev.off()

pdf(file="LOD chr=4.pdf")
plot(mqm_auto.1, chr=4, ylim=range(c(0, 20)))
plot(mqm_auto.2, chr=4, ylim=range(c(0, 20)))
plot(mqm_auto.3, chr=4, ylim=range(c(0, 20)))
dev.off()

# max LOD
mqmmax.1<-max(mqm_auto.1)
mqmmax.2<-max(mqm_auto.2)
mqmmax.3<-max(mqm_auto.3)

mqmmax <- data.frame(mqmmax.1, mqmmax.2, mqmmax.3)
write.csv(mqmmax,  file = "mqmmax.csv")

#genetic map
pdf(file="genetic map.pdf")
plotMap(qtl1)
dev.off()

pdf(file="genetic map2.pdf")
plotMap(qtl1, horizontal=TRUE)
plotMap(qtl1, horizontal=TRUE, chr=4)
dev.off()
save.image()

#missing genotype data
pdf(file="missing genotype data_F2_.pdf")
par(mfrow=c(1,1), las=1)
plotMissing(mapthis_F2_, main="missing genotype data")
dev.off()

pdf(file="missing genotype data_F2_10.pdf")
par(mfrow=c(1,1), las=1)
plotMissing(mapthis_F2_10, main="missing genotype data")
dev.off()

pdf(file="F2_ No. genotype.pdf")
par(mfrow=c(1,2), las=1)
plot(ntyped(mapthis_F2_, what="ind"), ylab="No. typed markers",
main="F2_ No. genotypes by individual", ylim=c(0,500))
plot(ntyped(mapthis_F2_, what="mar"), ylab="No. typed individuals",
main="F2_ No. genotypes by marker", ylim=c(0,200))
dev.off()

pdf(file="F2_ F2_2 No. genotype.pdf")
par(mfrow=c(2,2), las=1)
plot(ntyped(mapthis_F2_, what="ind"), ylab="No. typed markers"
, main="0 No. genotypes by individual", ylim=c(0,500))
plot(ntyped(mapthis_F2_, what="mar"), ylab="No. typed individuals"
, main="0 No. genotypes by marker", ylim=c(0,200))

plot(ntyped(mapthis_F2_2, what="ind"), ylab="No. typed markers"
, main="2 No. genotypes by individual", ylim=c(350,150))
plot(ntyped(mapthis_F2_2, what="mar"), ylab="No. typed individuals"
, main="2 No. genotypes by marker", ylim=c(100,200))
dev.off()

pdf(file="alle_freq.pdf")
par(mfrow=c(1,3), las=1)
for(i in 1:3)
plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "AB", "BB")[i], ylim=c(0,1))
dev.off()

#Evidence for segregation distortion
pdf(file="Evidence for segregation distortion.pdf")
par(mfrow=c(2,1), las=1)
plot(gt, ylab=expression(paste(-log[10], " P-value")))
plot(gt, lod=3:5, ylab="Genotype frequency")
abline(h=c(0.25, 0.5), lty=2, col="gray")
dev.off()

#estimated recombination fractions (upper-left triangle) and LOD scores (lower-right triangle)
pdf(file="recombinationfraction_LOD.pdf")
plotRF(mapthis_F2_10, alternate.chrid=TRUE, main="F2_10")
dev.off()

# graphical genotype
#The genotypes AA, AB, BB are displayed in the colors red, blue, and green
pdf(file="graphical_genotype.pdf")
geno.image(fill, main="FW_graphical_genotype", reorder=1)
geno.image(fill, main="nleaf_graphical_genotype", reorder=2)
geno.image(fill, main="water_graphical_genotype", reorder=3)
dev.off()

# graphical genotype
#The genotypes AA, AB, BB are displayed in the colors red, blue, and green
pdf(file="graphical_genotype_chr=4.pdf")
geno.image(fill, chr=4, main="FW_graphical_genotype", reorder=1)
geno.image(fill, chr=4, main="nleaf_graphical_genotype", reorder=2)
geno.image(fill, chr=4, main="water_graphical_genotype", reorder=3)
dev.off()

