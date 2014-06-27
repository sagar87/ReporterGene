library(ggplot2)
library(xlsx) 
library(plyr)


# READING AND RESHAPING DATA

## Read Beta-Gal and Reportergene Data
beta <- read.xlsx(file="/Volumes/Daten/Dropbox/Bachelor Results/0620_BETAGAL/Beta Gal 20.6.14.xlsx", sheetIndex=1)  
beta  <- beta[1:9, 1:13] # subset the real data
reporter  <- read.table("/Volumes/Daten/Dropbox/Bachelor Results/0618_REPORTER/Untitled 24.txt", header=TRUE)

## Conversion to a flatfile
val <- (c(t(beta[3, c(2:13)]), t(beta[4, c(2:13)]), t(beta[5, c(2:13)]))) # we have to use the transpose function in order to create flatfile
group <- c(rep("WT", 12), rep("E44A", 12), rep("N48A",12)) # name the Groups
treatment= c(rep("ctl", 6), rep("IFNg", 6), rep("ctl", 6), rep("IFNg", 6), rep("ctl", 6), rep("IFNg", 6)) # treatment

## Create a dataframe
df  <- data.frame(reporter=reporter, beta=val, group=group, treatment=treatment)
df <- transform(df, normalized=df$reporter/df$beta) # create a new column that contains beta-gal normalized data 
df$group  <- ordered(df$group, levels=c("WT", "E44A", "N48A")) # oder the levels

## Use dplyr to create statistics data.frame -> Counts N, Means,SD, SE 
stats <- ddply(df, .(group, treatment), summarise, N=sum(!is.na(normalized)), mean=mean(normalized, na.rm=TRUE), sd=sd(normalized, na.rm=TRUE), 
               se= sd(normalized, na.rm=TRUE) / sqrt(sum(!is.na(normalized))) )
stats$group  <- ordered(stats$group, levels=c("WT", "E44A", "N48A")) # again order the levels

## To create data in percent we have to select the unstimulated means of the statistics
percent = c(stats[1,4]/stats[1,4], stats[2,4]/stats[1,4], stats[3,4]/stats[3,4], stats[4,4]/stats[3,4], stats[5,4]/stats[5,4], stats[6,4]/stats[5,4])
percent <- round(100*percent,0)
stats  <- cbind(stats, percent)

# PLOTTING DATA

## Create percent plot
png(filename="/Volumes/Daten/Dropbox/Bachelor Results/percentage.png", width=600, height=500)
plot1 <- ggplot(stats, aes(x=factor(treatment), y=percent, fill=factor(treatment)))  + xlab("Treatment") + ylab("Percent Emission [%]") + ggtitle("Reportergen Assay") 
plot1 + geom_bar(stat="identity", position=position_dodge(), colour="black") + facet_grid(~ group) + theme(legend.position = "none") 
dev.off()

## Create a plot with normalized data 
png(filename="/Volumes/Daten/Dropbox/Bachelor Results/normalized.png", width=600, height=500)
plot2 <- ggplot(stats, aes(x=factor(treatment), y=mean, fill=factor(treatment)))  + xlab("Treatment") + ylab("Mean Light Emission (arbitrary units) ± SEM") + ggtitle("Mean Light Emission")
plot2 + geom_bar(stat="identity", position=position_dodge(), colour="black") + facet_grid(~ group) + geom_errorbar(aes(ymin=mean, ymax=mean+se), width=.5, colour="black") + theme(legend.position = "none") 
dev.off()

## Creates a boxplot
png(filename="/Volumes/Daten/Dropbox/Bachelor Results/boxplot.png", width=600, height=500)
plot3 <-qplot(factor(treatment), normalized, data = df)  + facet_grid(. ~ group) + geom_boxplot(aes(fill = factor(treatment))) + labs(title="Data", x="Treatment", y="Beta Gal Normalized Light Emission (arbitrary units)")
plot3  + theme(legend.position = "none") 
dev.off()
