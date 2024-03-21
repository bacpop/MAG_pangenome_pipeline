library(stringr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(tidyr)
library(ggsci)
library(ggpubr)

panaroo.stringency <- "strict"

# read in dfs generated from mmseqs2_performance.R and panaroo_performance.R
{
  cgt.rare.df <- read.csv("cgt_rare.csv")
  cgt.rare.df$analysis <- "CGT"
  
  cgt.core.df <- read.csv("cgt_core.csv")
  cgt.core.df$analysis <- "CGT"
  
  panaroo.rare.df <- read.csv("panaroo_rare.csv")
  panaroo.rare.df$analysis <- "PANAROO"
  
  panaroo.core.df <- read.csv("panaroo_core.csv")
  panaroo.core.df$analysis <- "PANAROO"
  
  rare.df <- rbind(cgt.rare.df, panaroo.rare.df)
  core.df <- rbind(cgt.core.df, panaroo.core.df)
}

# core genes
{
  core.df$variable <- NA
  core.df$split <- NA
  core.df$variable[core.df$type == "ori" & core.df$analysis == "CGT"] <- "True Total"
  core.df$split[core.df$type == "ori" & core.df$analysis == "CGT"] <- "mmseqs2"
  
  core.df$variable[core.df$type == "sim" & core.df$analysis == "CGT" & core.df$tool == "freq_only"] <- "Unadjusted"
  core.df$split[core.df$type == "sim" & core.df$analysis == "CGT" & core.df$tool == "freq_only"] <- "mmseqs2"
  
  core.df$variable[core.df$type == "sim" & core.df$analysis == "CGT" & core.df$tool == "cgt"] <- "Celebrimbor"
  core.df$split[core.df$type == "sim" & core.df$analysis == "CGT" & core.df$tool == "cgt"] <- "mmseqs2"
  
  core.df$variable[core.df$type == "ori" & core.df$analysis == "PANAROO" & core.df$tool == "panaroo" & core.df$stringency == panaroo.stringency] <- "True Total"
  core.df$split[core.df$type == "ori" & core.df$analysis == "PANAROO" & core.df$tool == "panaroo" & core.df$stringency == panaroo.stringency] <- "Panaroo"
  
  core.df$variable[core.df$type == "sim" & core.df$analysis == "PANAROO" & core.df$tool == "panaroo" & core.df$stringency == panaroo.stringency] <- "Unadjusted"
  core.df$split[core.df$type == "sim" & core.df$analysis == "PANAROO" & core.df$tool == "panaroo" & core.df$stringency == panaroo.stringency] <- "Panaroo"
  
  core.df$variable[core.df$type == "sim" & core.df$analysis == "PANAROO" & core.df$tool == "cgt" & core.df$stringency == panaroo.stringency] <- "Celebrimbor"
  core.df$split[core.df$type == "sim" & core.df$analysis == "PANAROO" & core.df$tool == "cgt" & core.df$stringency == panaroo.stringency] <- "Panaroo"
  
  core.df.subset <- subset(core.df, variable != "NA")
  core.df.subset <- subset(core.df.subset, error == 0.05 | error == 0.0)
  core.df.subset <- subset(core.df.subset, core_lim >= 0.7 & core_lim <= 0.99)
  
  core.df.subset$variable <- factor(core.df.subset$variable, levels = c("True Total", "Celebrimbor", "Unadjusted"))
  core.df.subset$split <- factor(core.df.subset$split, levels = c("mmseqs2", "Panaroo"))
  
  all_core_p <- ggplot(data=core.df.subset, aes(x=core_lim, y=core, colour=variable)) + facet_grid(.~split) + geom_point(size=2) + geom_line() + theme_light() + xlab("Core threshold") + ylab("Number of core genes") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 14), strip.text.y = element_text(size = 14, colour = "black", face = "bold"), legend.title=element_text(size=18,face="bold"), legend.text=element_text(size=16)) + guides(colour=guide_legend(title="Method")) + scale_colour_npg() + scale_y_continuous(breaks = seq(0, 100000, 1000)) + scale_x_continuous(breaks = seq(0, 1, 0.1))
  all_core_p
  
  ggsave(file="Celebrimbor_simulation_core_figure.svg", plot=all_core_p, width=11, height=5)
}


# rare genes
{
  rare.df$variable <- NA
  rare.df$split <- NA
  rare.df$variable[rare.df$type == "ori" & rare.df$analysis == "CGT"] <- "True Total"
  rare.df$split[rare.df$type == "ori" & rare.df$analysis == "CGT"] <- "mmseqs2"
  
  rare.df$variable[rare.df$type == "sim" & rare.df$analysis == "CGT" & rare.df$tool == "freq_only"] <- "Unadjusted"
  rare.df$split[rare.df$type == "sim" & rare.df$analysis == "CGT" & rare.df$tool == "freq_only"] <- "mmseqs2"
  
  rare.df$variable[rare.df$type == "sim" & rare.df$analysis == "CGT" & rare.df$tool == "cgt"] <- "Celebrimbor"
  rare.df$split[rare.df$type == "sim" & rare.df$analysis == "CGT" & rare.df$tool == "cgt"] <- "mmseqs2"
  
  rare.df$variable[rare.df$type == "ori" & rare.df$analysis == "PANAROO" & rare.df$tool == "panaroo" & rare.df$stringency == panaroo.stringency] <- "True Total"
  rare.df$split[rare.df$type == "ori" & rare.df$analysis == "PANAROO" & rare.df$tool == "panaroo" & rare.df$stringency == panaroo.stringency] <- "Panaroo"
  
  rare.df$variable[rare.df$type == "sim" & rare.df$analysis == "PANAROO" & rare.df$tool == "panaroo" & rare.df$stringency == panaroo.stringency] <- "Unadjusted"
  rare.df$split[rare.df$type == "sim" & rare.df$analysis == "PANAROO" & rare.df$tool == "panaroo" & rare.df$stringency == panaroo.stringency] <- "Panaroo"
  
  rare.df$variable[rare.df$type == "sim" & rare.df$analysis == "PANAROO" & rare.df$tool == "cgt" & rare.df$stringency == panaroo.stringency] <- "Celebrimbor"
  rare.df$split[rare.df$type == "sim" & rare.df$analysis == "PANAROO" & rare.df$tool == "cgt" & rare.df$stringency == panaroo.stringency] <- "Panaroo"
  
  rare.df.subset <- subset(rare.df, variable != "NA")
  rare.df.subset <- subset(rare.df.subset, error == 0.05 | error == 0.0)
  rare.df.subset <- subset(rare.df.subset, rare_lim <= 0.25)
  
  rare.df.subset$variable <- factor(rare.df.subset$variable, levels = c("True Total", "Celebrimbor", "Unadjusted"))
  rare.df.subset$split <- factor(rare.df.subset$split, levels = c("mmseqs2", "Panaroo"))
  
  all_rare_p <- ggplot(data=rare.df.subset, aes(x=rare_lim, y=rare, colour=variable)) + facet_grid(.~split) + geom_point(size=2) + geom_line() + theme_light() + xlab("Rare threshold") + ylab("Number of rare genes") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 14), strip.text.y = element_text(size = 14, colour = "black", face = "bold"), legend.title=element_text(size=18,face="bold"), legend.text=element_text(size=16)) + guides(colour=guide_legend(title="Method")) + scale_colour_npg() + scale_y_continuous(breaks = seq(0, 100000, 5000)) + scale_x_continuous(breaks = seq(0, 1, 0.1))
  all_rare_p
  
  ggsave(file="Celebrimbor_simulation_rare_figure.svg", plot=all_rare_p, width=11, height=5)
}

pub.table <- subset(core.df, core_lim == 0.95 & (error == 0.05 | error == 0) & variable != "NA")

