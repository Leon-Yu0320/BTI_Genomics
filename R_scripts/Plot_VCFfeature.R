#!/usr/bin/Rscript

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(getopt)
library(ggplot2)

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
command = matrix(c(
    'DP', 'd', 1, "character",
    'QUAL', 'q', 1, "character",
    'QD', 'i', 1, "character"), byrow=TRUE, ncol=4)

args = getopt(command)

##help info-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if (!is.null(args$help) || is.null(args$DP) || is.null(args$QUAL) || is.null(args$QD)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Load depth
DP <- read.table(args$DP, header = F)

#Load quality
QUAL <- read.table(args$QUAL, header = F)

#Load quality score by depth (cutoff: 2)
QD <- read.table(args$QD,header = F)
QD$V1 <- as.numeric(QD$V1)

## -----------------------------------------------------------------------------------Plot depth-----------------------------------------------------------------------------------------


pdf("depth.pdf")

ggplot(data = DP, aes(x= V1)) + 
  geom_density(alpha = 0.3, color="darkblue", fill="lightblue") +
  geom_vline(aes(xintercept=400),color="blue", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=2000),color="blue", linetype="dashed", size=1) +
  xlab("Total mapping depth among F1 plants") +  ylab("Density") +
  xlim(0,3000) +
  theme_bw()

dev.off()

## -----------------------------------------------------------------------------------Plot quality score-------------------------------------------------------------------------------------
pdf("Quality_score.pdf")

ggplot(data = QUAL, aes(x= V1)) + 
  geom_density(alpha = 0.3, color="darkblue", fill="lightblue") +
  geom_vline(aes(xintercept=50),color="blue", linetype="dashed", size=1) +
  xlab("Quality score across SNPs among F1 plants") +  ylab("Density") +
  xlim(0,2000) +
  theme_bw()

dev.off()

## -----------------------------------------------------------------------------------Plot quality by depth-------------------------------------------------------------------------------------
pdf("Quality_by_depth.pdf")

ggplot(data = QD, aes(x= V1)) + 
  geom_density(alpha = 0.3, color="darkblue", fill="lightblue") +
  geom_vline(aes(xintercept=5),color="blue", linetype="dashed", size=1) +
  xlab("Quality by depth among F1 population") +  ylab("Density") +
  xlim(0,80) +
  theme_bw()

dev.off()
