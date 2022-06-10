library("ggplot2")
library("scales")
library("RColorBrewer")
library("plyr")
library("magrittr")
library("dplyr")
library("grid")
library("cowplot")
library("here")
library(plyr)

update_geom_defaults("line", list(size = 1.0))
update_geom_defaults("step", list(size = 1.0))
update_geom_defaults("point", list(size = 1.5, stroke=0.8))
theme_set(theme_cowplot(
  font_size = 13,
  line_size = 0.75
) + theme(
  plot.margin = unit(c(0, 0, 0, 0), "cm")
))

group.colors <- c("DALI" = "#A495AF", "TMalign" = "#3B76AF", "Foldseek" = "#D1382C",
                  "TMalign-fast" ="#ADCDE1", "CE"="#559E3E", "CLE-SW"="#EE8632", 
                  "3D-BLAST" = "#643F95", "MMseqs2" = "#F4C27B" )

setwd("/Users/charlotte/plots-fs/SUPP_FIG")

##############################################################################################################################################
# SCOP TP/FP distribution
fp_fs = read.delim("FP_FS_AVERAGE.txt", header = FALSE)
fp_tma = read.delim("FP_TMA_AVERAGE.txt", header = FALSE)
fp_dali = read.delim("FP_DALI_AVERAGE.txt", header = FALSE)
fp_cle = read.delim("FP_CLESW_AVERAGE.txt", header = FALSE)
tp_fs_fold = read.delim("TP_FS_FOLD_AVERAGE.txt", header = FALSE)
tp_fs_sfam = read.delim("TP_FS_SFAM_AVERAGE.txt", header = FALSE)
tp_fs_fam = read.delim("TP_FS_FAM_AVERAGE.txt", header = FALSE)
tp_tma_fold = read.delim("TP_TMA_FOLD_AVERAGE.txt", header = FALSE)
tp_tma_sfam = read.delim("TP_TMA_SFAM_AVERAGE.txt", header = FALSE)
tp_tma_fam = read.delim("TP_TMA_FAM_AVERAGE.txt", header = FALSE)
tp_dali_fold = read.delim("TP_DALI_FOLD_AVERAGE.txt", header = FALSE)
tp_dali_sfam = read.delim("TP_DALI_SFAM_AVERAGE.txt", header = FALSE)
tp_dali_fam = read.delim("TP_DALI_FAM_AVERAGE.txt", header = FALSE)
tp_cle_fold = read.delim("TP_CLESW_FOLD_AVERAGE.txt", header = FALSE)
tp_cle_sfam = read.delim("TP_CLESW_SFAM_AVERAGE.txt", header = FALSE)
tp_cle_fam = read.delim("TP_CLESW_FAM_AVERAGE.txt", header = FALSE)

fp_all = rbind(fp_fs, fp_tma)
fp_all = rbind(fp_all, fp_dali)
fp_all = rbind(fp_all, fp_cle)
tp_fold_all = rbind(tp_fs_fold, tp_tma_fold)
tp_fold_all = rbind(tp_fold_all, tp_dali_fold)
tp_fold_all = rbind(tp_fold_all, tp_cle_fold)
tp_sfam_all = rbind(tp_fs_sfam, tp_tma_sfam)
tp_sfam_all = rbind(tp_sfam_all, tp_dali_sfam)
tp_sfam_all = rbind(tp_sfam_all, tp_cle_sfam)
tp_fam_all = rbind(tp_fs_fam, tp_tma_fam)
tp_fam_all = rbind(tp_fam_all, tp_dali_fam)
tp_fam_all = rbind(tp_fam_all, tp_cle_fam)

tp_all = rbind(tp_fold_all, tp_sfam_all)
tp_all = rbind(tp_all, tp_fam_all)

# all tools together now!
# Find the median of each group
fp_all$V2 <- 'FP'
tp_fold_all$V2 <- 'TP'

fold <- rbind(fp_all,tp_fold_all)
cdat_all_fold <- ddply(fold, "V2", summarise, V1.median=median(V1))

DENS_ALL_FOLD <- ggplot(fold, aes(x=V1, colour=V2)) +
  geom_density(size=1.5) +
  geom_vline(data=cdat_all_fold, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(), legend.position=c(0.7,0.8), axis.text.x = element_blank()) +  ylab("Frequency") + xlab("") +
  annotate("text",  x = 0.8, y = 4, vjust = -0.5, size=5, label = "TP Fold, FP Fold") +
  ylim(0,4.4) +
  xlim(0,1)

tp_sfam_all$V2 <- 'TP'

sfam <- rbind(fp_all,tp_sfam_all)
cdat_all_sfam <- ddply(sfam, "V2", summarise, V1.median=median(V1))

DENS_ALL_SFAM <-ggplot(sfam, aes(x=V1, colour=V2)) +
  geom_density(size=1.5) +
  geom_vline(data=cdat_all_sfam, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(),axis.text.x = element_blank(),   axis.text.y = element_blank(),legend.position="none") +  ylab("") + xlab("") +
  annotate("text",  x = 0.75, y = 4, vjust = -0.5, size=5, label = "TP Super Family, FP Fold") +
  ylim(0,4.4)+
  xlim(0,1)

tp_fam_all$V2 <- 'TP'

fam <- rbind(fp_all,tp_fam_all)
cdat_all_fam <- ddply(fam, "V2", summarise, V1.median=median(V1))

DENS_ALL_FAM <- ggplot(fam, aes(x=V1, colour=V2)) +
  geom_density(size=1.5) +
  geom_vline(data=cdat_all_fam, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(), legend.position="none") +  ylab("Frequency") + xlab("lddt") +
  theme(legend.title=element_blank(), axis.text.y = element_blank(), legend.position="none") +  ylab("") + xlab("lddt") +
  annotate("text",  x = 0.8, y = 4, vjust = -0.5, size=5, label = "TP Family, FP Fold") +
  ylim(0,4.4)+
  xlim(0,1)

tp_all$V2 <- 'TP'

all <- rbind(fp_all,tp_all)
cdat_all_all <- ddply(all, "V2", summarise, V1.median=median(V1))

DENS_ALL_ALL <- ggplot(all, aes(x=V1, colour=V2)) +
  geom_density(size=1.5) +
  geom_vline(data=cdat_all_all, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(), axis.text.y = element_blank(), legend.position="none") +  ylab("") + xlab("lddt") +
  annotate("text",  x = 0.75, y = 4, vjust = -0.5, size=5, label = "TP all categories, FP Fold") +
  ylim(0,4.4) +
  xlim(0,1)

DENS_ALL_ALL

grid <- plot_grid( DENS_ALL_FOLD, DENS_ALL_SFAM, DENS_ALL_FAM,DENS_ALL_ALL,
                    ncol = 2, nrow = 2, labels=c("a","b","c","d"))
grid
ggsave("SUPP_SCOP_ALL_LDDT_AVERAGE.pdf", grid, width=25, height=25, units="cm")


##############################################################################################################################################
# Alphafold 1/5/20 hits distribution
af_fs1 = read.delim("./AF_DIST1_PLDDT/AF_FS.txt", header = FALSE)
af_mm1 =read.delim("./AF_DIST1_PLDDT/AF_MM.txt", header = FALSE)
af_tm1 =read.delim("./AF_DIST1_PLDDT/AF_TMA.txt", header = FALSE)
af_dali1 =read.delim("./AF_DIST1_PLDDT/AF_DALI.txt", header = FALSE)
af_cle1 =read.delim("./AF_DIST1_PLDDT/AF_CLE.txt", header = FALSE)

af_fs1$V2 <- 'Foldseek'
af_mm1$V2 <- 'MMseqs2'
af_tm1$V2 <- 'TMalign'
af_dali1$V2 <- 'Dali'
af_cle1$V2 <- 'CLE-SW'

af_all1 <- rbind(af_fs1,af_mm1)
af_all1 <- rbind(af_all1,af_tm1)
af_all1 <- rbind(af_all1,af_dali1)
af_all1 <- rbind(af_all1,af_cle1)


# Find the median of each group
cdat_af_all1 <- ddply(af_all1, "V2", summarise, V1.median=median(V1))

# Density plots with medians
AF_PLOT_DENS1 <- ggplot(af_all1, aes(x=V1, colour=V2)) +
  geom_density(size=1.5) +
  geom_vline(data=cdat_af_all1, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(),axis.text.x = element_blank(),  legend.position=c(0.1,0.7)) +  ylab("Frequency") + xlab("") +
  scale_color_manual(values=group.colors, limits = force) + 
  annotate("text",  x = 0.9, y = 6.5, vjust = -0.5, size=5, label = "1st hit") +
  ylim(0,7.0)+
  xlim(0,1)
AF_PLOT_DENS1

af_fs5 = read.delim("./AF_DIST5_PLDDT/AF_FS.txt", header = FALSE)
af_mm5 =read.delim("./AF_DIST5_PLDDT/AF_MM.txt", header = FALSE)
af_tm5 =read.delim("./AF_DIST5_PLDDT/AF_TMA.txt", header = FALSE)
af_dali5 =read.delim("./AF_DIST5_PLDDT/AF_DALI.txt", header = FALSE)
af_cle5 =read.delim("./AF_DIST5_PLDDT/AF_CLE.txt", header = FALSE)

af_fs5$V2 <- 'Foldseek'
af_mm5$V2 <- 'MMseqs2'
af_tm5$V2 <- 'TMalign'
af_dali5$V2 <- 'Dali'
af_cle5$V2 <- 'CLE-SW'


af_all5 <- rbind(af_fs5,af_mm5)
af_all5 <- rbind(af_all5,af_tm5)
af_all5 <- rbind(af_all5,af_dali5)
af_all5 <- rbind(af_all5,af_cle5)


# Find the median of each group
cdat_af_all5 <- ddply(af_all5, "V2", summarise, V1.median=median(V1))

# Density plots with medians
AF_PLOT_DENS5 <- ggplot(af_all5, aes(x=V1, colour=V2)) +
  geom_density(size=1.5) +
  geom_vline(data=cdat_af_all5, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(),axis.text.x = element_blank(),  legend.position="none") +  ylab("Frequency") + xlab("") +
  scale_color_manual(values=group.colors, limits = force) + 
  annotate("text",  x = 0.9, y = 6.5, vjust = -0.5, size=5, label = "5 hits") +
  ylim(0,7.0)+
  xlim(0,1)
AF_PLOT_DENS5

af_fs20 = read.delim("./AF_DIST20_PLDDT/AF_FS.txt", header = FALSE)
af_mm20 =read.delim("./AF_DIST20_PLDDT/AF_MM.txt", header = FALSE)
af_tm20 =read.delim("./AF_DIST20_PLDDT/AF_TMA.txt", header = FALSE)
af_dali20 =read.delim("./AF_DIST20_PLDDT/AF_DALI.txt", header = FALSE)
af_cle20 =read.delim("./AF_DIST20_PLDDT/AF_CLE.txt", header = FALSE)

af_fs20$V2 <- 'Foldseek'
af_mm20$V2 <- 'MMseqs2'
af_tm20$V2 <- 'TMalign'
af_dali20$V2 <- 'Dali'
af_cle20$V2 <- 'CLE-SW'

af_all20 <- rbind(af_fs20,af_mm20)
af_all20 <- rbind(af_all20,af_tm20)
af_all20 <- rbind(af_all20,af_dali20)
af_all20 <- rbind(af_all20,af_cle20)
# Find the median of each group
cdat_af_all20 <- ddply(af_all20, "V2", summarise, V1.median=median(V1))

# Density plots with medians
AF_PLOT_DENS20 <- ggplot(af_all20, aes(x=V1, colour=V2)) +
  geom_density(size=1.5) +
  geom_vline(data=cdat_af_all20, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(), legend.position="none") +  ylab("Frequency") + xlab("lddt") +
  annotate("text",  x = 0.9, y = 6.5, vjust = -0.5, size=5, label = "20 hits") +
  scale_color_manual(values=group.colors, limits =  force) +  
  ylim(0,7.0)+
  xlim(0,1)
AF_PLOT_DENS20

grid2 <- plot_grid( AF_PLOT_DENS1, AF_PLOT_DENS5, AF_PLOT_DENS20, labels=c("a","b","c"), 
                    ncol = 1, nrow = 3)
grid2
ggsave("SUPP_AF_HITS_DIST.pdf", grid2, width=15, height=40, units="cm")

##############################################################################################################################################
# Coverage AF multidomain benchmark

### MULTI WITH PLDDT CHECK AND AVERAGE LDDT INSTEAD OF PER RESIDUE LDDT
### MULTI WITH PLDDT THRESHOLD
# Coverage TP 0.4
fs_countPA04 = read.delim("./04_final/foldseek_lddt_cutoff_final_04TP.txt", header = FALSE, sep=' ')
dali_countPA04 = read.delim("./04_final/dali_lddt_cutoff_final_04TP.txt", header = FALSE, sep=' ')
ms_countPA04 = read.delim("./04_final/mmseqs_lddt_cutoff_final_04TP.txt", header = FALSE, sep=' ')
tm_countPA04 = read.delim("./04_final/tmalign_lddt_cutoff_final_04TP.txt", header = FALSE, sep=' ')
clesw_countPA04 = read.delim("./04_final/clesw_lddt_cutoff_final_04TP.txt", header = FALSE, sep=' ')

colnames(ms_countPA04) <- c("V1", "V3")
colnames(tm_countPA04) <- c("V1", "V4")
colnames(dali_countPA04) <- c("V1", "V5")
colnames(clesw_countPA04) <- c("V1", "V6")

myfulldataPA04 <- merge(fs_countPA04, ms_countPA04, by.x='V1', all.x=TRUE)
myfulldataPA04 <- merge(myfulldataPA04, dali_countPA04, by.x='V1', all.x=TRUE)
myfulldataPA04 <- merge(myfulldataPA04, tm_countPA04, by.x='V1', all.x=TRUE)
myfulldataPA04 <- merge(myfulldataPA04, clesw_countPA04, by.x='V1', all.x=TRUE)

MULTI_PA04 <- ggplot(myfulldataPA04, aes(x=V1, y = value, color = group.colors)) + 
  scale_size_area() + theme(legend.title=element_blank(), legend.position="none",axis.text.y = element_blank()) +
  geom_point(aes( y = V2, col = "Foldseek"), size=2, shape=19, stroke=0.7) +
  geom_point(aes( y = V5, col = "DALI"), size=2, shape=19, stroke=0.7) +
  geom_point(aes(y = V3, col = "MMseqs2"), size=2, shape=19, stroke=0.7) +
  geom_point(aes(y = V4, col = "TMalign"), size=2,shape=19, stroke=0.7) +
  geom_point(aes(y = V6, col = "CLE"), size=2,shape=19, stroke=0.7) +
  geom_line(aes( y = V2, col = "Foldseek"),lwd=1)+
  geom_line(aes( y = V5, col = "DALI"),lwd=1)+
  geom_line(aes(y = V3, col = "MMseqs2"),lwd=1)+
  geom_line(aes(y = V4, col = "TMalign"),lwd=1)+
  geom_line(aes(y = V6, col = "CLE"),lwd=1) +
  ylim(0,1.0) +
  xlim(0,20) +
  scale_color_manual(values=group.colors, limits = force) + 
  annotate("text",  x = 15, y = .95, vjust = -0.5, size=5, label = "TP threshold 0.4") +
  xlab("TP hits up to 1st FP") + ylab("")

MULTI_PA04

# Coverage TP 0.5
fs_countPA05 = read.delim("./05_final/foldseek_lddt_cutoff_final_05TP.txt", header = FALSE, sep=' ')
dali_countPA05 = read.delim("./05_final/dali_lddt_cutoff_final_05TP.txt", header = FALSE, sep=' ')
ms_countPA05 = read.delim("./05_final/mmseqs_lddt_cutoff_final_05TP.txt", header = FALSE, sep=' ')
tm_countPA05 = read.delim("./05_final/tmalign_lddt_cutoff_final_05TP.txt", header = FALSE, sep=' ')
clesw_countPA05 = read.delim("./05_final/clesw_lddt_cutoff_final_05TP.txt", header = FALSE, sep=' ')

colnames(ms_countPA05) <- c("V1", "V3")
colnames(tm_countPA05) <- c("V1", "V4")
colnames(dali_countPA05) <- c("V1", "V5")
colnames(clesw_countPA05) <- c("V1", "V6")

myfulldataPA05 <- merge(fs_countPA05, ms_countPA05, by.x='V1', all.x=TRUE)
myfulldataPA05 <- merge(myfulldataPA05, dali_countPA05, by.x='V1', all.x=TRUE)
myfulldataPA05 <- merge(myfulldataPA05, tm_countPA05, by.x='V1', all.x=TRUE)
myfulldataPA05 <- merge(myfulldataPA05, clesw_countPA05, by.x='V1', all.x=TRUE)

MULTI_PA05 <- ggplot(myfulldataPA05, aes(x=V1, y = value, color = group.colors)) + 
  scale_size_area() + theme(legend.title=element_blank(), legend.position="none") +
  geom_point(aes( y = V2, col = "Foldseek"), size=2, shape=19, stroke=0.7) +
  geom_point(aes( y = V5, col = "DALI"), size=2, shape=19, stroke=0.7) +
  geom_point(aes(y = V3, col = "MMseqs2"), size=2, shape=19, stroke=0.7) +
  geom_point(aes(y = V4, col = "TMalign"), size=2,shape=19, stroke=0.7) +
  geom_point(aes(y = V6, col = "CLE"), size=2,shape=19, stroke=0.7) +
  geom_line(aes( y = V2, col = "Foldseek"),lwd=1)+
  geom_line(aes( y = V5, col = "DALI"),lwd=1)+
  geom_line(aes(y = V3, col = "MMseqs2"),lwd=1)+
  geom_line(aes(y = V4, col = "TMalign"),lwd=1)+
  geom_line(aes(y = V6, col = "CLE"),lwd=1) +
  ylim(0,1.0) +
  xlim(0,20) +
  scale_color_manual(values=group.colors, limits = force) + 
  annotate("text",  x = 15, y = .95, vjust = -0.5, size=5, label = "TP threshold 0.5") +
  xlab("TP hits up to 1st FP") + ylab("Query coverage")

MULTI_PA05

# Coverage TP 0.6
fs_countPA06 = read.delim("./06_final/foldseek_lddt_cutoff_final_06TP.txt", header = FALSE, sep=' ')
dali_countPA06 = read.delim("./06_final/dali_lddt_cutoff_final_06TP.txt", header = FALSE, sep=' ')
ms_countPA06 = read.delim("./06_final/mmseqs_lddt_cutoff_final_06TP.txt", header = FALSE, sep=' ')
tm_countPA06 = read.delim("./06_final/tmalign_lddt_cutoff_final_06TP.txt", header = FALSE, sep=' ')
clesw_countPA06 = read.delim("./06_final/clesw_lddt_cutoff_final_06TP.txt", header = FALSE, sep=' ')

colnames(ms_countPA06) <- c("V1", "V3")
colnames(tm_countPA06) <- c("V1", "V4")
colnames(dali_countPA06) <- c("V1", "V5")
colnames(clesw_countPA06) <- c("V1", "V6")

myfulldataPA06 <- merge(fs_countPA06, ms_countPA06, by.x='V1', all.x=TRUE)
myfulldataPA06 <- merge(myfulldataPA06, dali_countPA06, by.x='V1', all.x=TRUE)
myfulldataPA06 <- merge(myfulldataPA06, tm_countPA06, by.x='V1', all.x=TRUE)
myfulldataPA06 <- merge(myfulldataPA06, clesw_countPA06, by.x='V1', all.x=TRUE)

MULTI_PA06 <- ggplot(myfulldataPA06, aes(x=V1, y = value, color = group.colors)) + 
  scale_size_area() + theme(legend.title=element_blank(), legend.position="none",axis.text.x = element_blank(), ) +
  geom_point(aes( y = V2, col = "Foldseek"), size=2, shape=19, stroke=0.7) +
  geom_point(aes( y = V5, col = "DALI"), size=2, shape=19, stroke=0.7) +
  geom_point(aes(y = V3, col = "MMseqs2"), size=2, shape=19, stroke=0.7) +
  geom_point(aes(y = V4, col = "TMalign"), size=2,shape=19, stroke=0.7) +
  geom_point(aes(y = V6, col = "CLE"), size=2,shape=19, stroke=0.7) +
  geom_line(aes( y = V2, col = "Foldseek"),lwd=1)+
  geom_line(aes( y = V5, col = "DALI"),lwd=1)+
  geom_line(aes(y = V3, col = "MMseqs2"),lwd=1)+
  geom_line(aes(y = V4, col = "TMalign"),lwd=1)+
  geom_line(aes(y = V6, col = "CLE"),lwd=1) +
  scale_color_manual(values=group.colors, limits = force) + 
  ylim(0,1.0) +
  xlim(0,20) +
  ggtitle('TP 0.6, FP 0.25') +
  xlab("TP hits up to 1st FP") + ylab("Query coverage")

MULTI_PA06

# Coverage TP 0.7
fs_countPA07 = read.delim("./07_final/foldseek_lddt_cutoff_final_07TP.txt", header = FALSE, sep=' ')
dali_countPA07 = read.delim("./07_final/dali_lddt_cutoff_final_07TP.txt", header = FALSE, sep=' ')
ms_countPA07 = read.delim("./07_final/mmseqs_lddt_cutoff_final_07TP.txt", header = FALSE, sep=' ')
tm_countPA07 = read.delim("./07_final/tmalign_lddt_cutoff_final_07TP.txt", header = FALSE, sep=' ')
clesw_countPA07 = read.delim("./07_final/clesw_lddt_cutoff_final_07TP.txt", header = FALSE, sep=' ')

colnames(ms_countPA07) <- c("V1", "V3")
colnames(tm_countPA07) <- c("V1", "V4")
colnames(dali_countPA07) <- c("V1", "V5")
colnames(clesw_countPA07) <- c("V1", "V6")

myfulldataPA07 <- merge(fs_countPA07, ms_countPA07, by.x='V1', all.x=TRUE)
myfulldataPA07 <- merge(myfulldataPA07, dali_countPA07, by.x='V1', all.x=TRUE)
myfulldataPA07 <- merge(myfulldataPA07, tm_countPA07, by.x='V1', all.x=TRUE)
myfulldataPA07 <- merge(myfulldataPA07, clesw_countPA07, by.x='V1', all.x=TRUE)

MULTI_PA07 <- ggplot(myfulldataPA07, aes(x=V1, y = value, color = group.colors)) + 
  scale_size_area() + theme(legend.title=element_blank(), legend.position="none",axis.text.y = element_blank(), axis.text.x = element_blank(), ) +
  geom_point(aes( y = V2, col = "Foldseek"), size=2, shape=19, stroke=0.7) +
  geom_point(aes( y = V5, col = "DALI"), size=2, shape=19, stroke=0.7) +
  geom_point(aes(y = V3, col = "MMseqs2"), size=2, shape=19, stroke=0.7) +
  geom_point(aes(y = V4, col = "TMalign"), size=2,shape=19, stroke=0.7) +
  geom_point(aes(y = V6, col = "CLE"), size=2,shape=19, stroke=0.7) +
  geom_line(aes( y = V2, col = "Foldseek"),lwd=1)+
  geom_line(aes( y = V5, col = "DALI"),lwd=1)+
  geom_line(aes(y = V3, col = "MMseqs2"),lwd=1)+
  geom_line(aes(y = V4, col = "TMalign"),lwd=1)+
  geom_line(aes(y = V6, col = "CLE"),lwd=1) +
  scale_color_manual(values=group.colors, limits = force) + 
  ylim(0,1.0) +
  xlim(0,20) +
  annotate("text",  x = 15, y = .95, vjust = -0.5, size=5, label = "TP threshold 0.7") +
  xlab("") + ylab("") 

MULTI_PA07

# Coverage TP 0.8
fs_countPA08 = read.delim("./08_final/foldseek_lddt_cutoff_final_08TP.txt", header = FALSE, sep=' ')
dali_countPA08 = read.delim("./08_final/dali_lddt_cutoff_final_08TP.txt", header = FALSE, sep=' ')
ms_countPA08 = read.delim("./08_final/mmseqs_lddt_cutoff_final_08TP.txt", header = FALSE, sep=' ')
tm_countPA08 = read.delim("./08_final/tmalign_lddt_cutoff_final_08TP.txt", header = FALSE, sep=' ')
clesw_countPA08 = read.delim("./08_final/clesw_lddt_cutoff_final_08TP.txt", header = FALSE, sep=' ')

colnames(ms_countPA08) <- c("V1", "V3")
colnames(tm_countPA08) <- c("V1", "V4")
colnames(dali_countPA08) <- c("V1", "V5")
colnames(clesw_countPA08) <- c("V1", "V6")

myfulldataPA08 <- merge(fs_countPA08, ms_countPA08, by.x='V1', all.x=TRUE)
myfulldataPA08 <- merge(myfulldataPA08, dali_countPA08, by.x='V1', all.x=TRUE)
myfulldataPA08 <- merge(myfulldataPA08, tm_countPA08, by.x='V1', all.x=TRUE)
myfulldataPA08 <- merge(myfulldataPA08, clesw_countPA08, by.x='V1', all.x=TRUE)

MULTI_PA08 <- ggplot(myfulldataPA08, aes(x=V1, y = value, color = group.colors)) + 
  scale_size_area() + theme(legend.title=element_blank(), legend.position=c(0.6,0.7), axis.text.x = element_blank(), ) +
  geom_point(aes( y = V2, col = "Foldseek"), size=2, shape=19, stroke=0.7) +
  geom_point(aes( y = V5, col = "DALI"), size=2, shape=19, stroke=0.7) +
  geom_point(aes(y = V3, col = "MMseqs2"), size=2, shape=19, stroke=0.7) +
  geom_point(aes(y = V4, col = "TMalign"), size=2,shape=19, stroke=0.7) +
  geom_point(aes(y = V6, col = "CLE"), size=2,shape=19, stroke=0.7) +
  geom_line(aes( y = V2, col = "Foldseek"),lwd=1)+
  geom_line(aes( y = V5, col = "DALI"),lwd=1)+
  geom_line(aes(y = V3, col = "MMseqs2"),lwd=1)+
  geom_line(aes(y = V4, col = "TMalign"),lwd=1)+
  geom_line(aes(y = V6, col = "CLE"),lwd=1) +
  scale_color_manual(values=group.colors, limits = force) + 
  ylim(0,1.0) +
  xlim(0,20) +
  annotate("text",  x = 15, y = .95, vjust = -0.5, size=5, label = "TP threshold 0.8") +
  xlab("") + ylab("Query coverage")

MULTI_PA08


grid3 <- plot_grid(MULTI_PA08,MULTI_PA07, MULTI_PA05, MULTI_PA04,
                    ncol = 2, nrow = 2, labels=c("a","b","c","d"))
grid3
ggsave("SUPP_AF_COVERAGE.pdf", grid3, width=20, height=20, units="cm")




# EXTRA
#######################################################################################################################
# For distribution plots per tool
fp_fs_lddt = read.delim("FP_FS_AVERAGE.txt", header = FALSE)
fp_tma_lddt = read.delim("FP_TMA_AVERAGE.txt", header = FALSE)
fp_dali_lddt = read.delim("FP_DALI_AVERAGE.txt", header = FALSE)
fp_cle_lddt = read.delim("FP_CLESW_AVERAGE.txt", header = FALSE)
tp_fs_fold_lddt = read.delim("TP_FS_FOLD_AVERAGE.txt", header = FALSE)
tp_fs_sfam_lddt = read.delim("TP_FS_SFAM_AVERAGE.txt", header = FALSE)
tp_fs_fam_lddt = read.delim("TP_FS_FAM_AVERAGE.txt", header = FALSE)
tp_tma_fold_lddt = read.delim("TP_TMA_FOLD_AVERAGE.txt", header = FALSE)
tp_tma_sfam_lddt = read.delim("TP_TMA_SFAM_AVERAGE.txt", header = FALSE)
tp_tma_fam_lddt = read.delim("TP_TMA_FAM_AVERAGE.txt", header = FALSE)
tp_dali_fold_lddt = read.delim("TP_DALI_FOLD_AVERAGE.txt", header = FALSE)
tp_dali_sfam_lddt = read.delim("TP_DALI_SFAM_AVERAGE.txt", header = FALSE)
tp_dali_fam_lddt = read.delim("TP_DALI_FAM_AVERAGE.txt", header = FALSE)
tp_cle_fold_lddt = read.delim("TP_CLESW_FOLD_AVERAGE.txt", header = FALSE)
tp_cle_sfam_lddt = read.delim("TP_CLESW_SFAM_AVERAGE.txt", header = FALSE)
tp_cle_fam_lddt = read.delim("TP_CLESW_FAM_AVERAGE.txt", header = FALSE)

# FOLDSEEk
# FOLD
fp_fs_lddt$V2 <- 'FP_FS'
tp_fs_fold_lddt$V2 <- 'TP_FS_FOLD'

tpfp_fs_fold <- rbind(fp_fs_lddt,tp_fs_fold_lddt)

# Find the median of each group
cdat <- ddply(tpfp_fs_fold, "V2", summarise, V1.median=median(V1))
cdat

# Overlaid histograms with medians
DIST_HIST_FOLD <- ggplot(tpfp_fs_fold, aes(x=V1, fill=V2)) +
  geom_histogram(binwidth=.05, alpha=.5, position="identity") +
  geom_vline(data=cdat, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(), legend.position=c(0.6,0.8))  + ylab("Count") + xlab("") + ggtitle("Foldseek")

# Density plots with medians
DIST_DENS_FOLD <- ggplot(tpfp_fs_fold, aes(x=V1, colour=V2)) +
  geom_density() +
  geom_vline(data=cdat, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(), legend.position=c(0.6,0.8)) +  ylab("Density") + xlab("") + ggtitle("Foldseek") +
  ylim(0,7.0) +
  xlim(0,1)

# SFAMILY

fp_fs_lddt$V2 <- 'FP_FS'
tp_fs_sfam_lddt$V2 <- 'TP_FS_SFAM'

tpfp_fs_sfam <- rbind(fp_fs_lddt,tp_fs_sfam_lddt)

# Find the median of each group
cdat <- ddply(tpfp_fs_sfam, "V2", summarise, V1.median=median(V1))
cdat

# Overlaid histograms with medians
DIST_HIST_SFAM <- ggplot(tpfp_fs_sfam, aes(x=V1, fill=V2)) +
  geom_histogram(binwidth=.05, alpha=.5, position="identity") +
  geom_vline(data=cdat, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(), legend.position=c(0.6,0.8)) + ylab("Count") + xlab("")

# Density plots with medians
DIST_DENS_SFAM <- ggplot(tpfp_fs_sfam, aes(x=V1, colour=V2)) +
  geom_density() +
  geom_vline(data=cdat, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(), legend.position=c(0.6,0.8)) +  ylab("Density") + xlab("") +
  ylim(0,7.0)+
  xlim(0,1)

# FAMILY

fp_fs_lddt$V2 <- 'FP_FS'
tp_fs_fam_lddt$V2 <- 'TP_FS_FAM'

tpfp_fs_fam <- rbind(fp_fs_lddt,tp_fs_fam_lddt)

# Find the median of each group
cdat <- ddply(tpfp_fs_fam, "V2", summarise, V1.median=median(V1))
cdat

# Overlaid histograms with medians
DIST_HIST_FAM <- ggplot(tpfp_fs_fam, aes(x=V1, fill=V2)) +
  geom_histogram(binwidth=.05, alpha=.5, position="identity") +
  geom_vline(data=cdat, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(), legend.position=c(0.7,0.8)) +  ylab("Count") + xlab("lddt")
# Density plots with medians
DIST_DENS_FAM <- ggplot(tpfp_fs_fam, aes(x=V1, colour=V2)) +
  geom_density() +
  geom_vline(data=cdat, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(), legend.position=c(0.7,0.8)) +  ylab("Density") + xlab("lddt") +
  ylim(0,7.0)+
  xlim(0,1)

# TM-ALIGN
# FOLD
fp_tma_lddt$V2 <- 'FP_TMA'
tp_tma_fold_lddt$V2 <- 'TP_TMA_FOLD'

tpfp_tma_fold <- rbind(fp_tma_lddt,tp_tma_fold_lddt)

# Find the median of each group
cdat <- ddply(tpfp_tma_fold, "V2", summarise, V1.median=median(V1))
cdat

# Overlaid histograms with medians
DIST_HIST_FOLD2 <- ggplot(tpfp_tma_fold, aes(x=V1, fill=V2)) +
  geom_histogram(binwidth=.05, alpha=.5, position="identity") +
  geom_vline(data=cdat, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(), legend.position=c(0.6,0.8)) +  ylab("") + xlab("") + ggtitle("TM-Align")


# Density plots with medians
DIST_DENS_FOLD2 <- ggplot(tpfp_tma_fold, aes(x=V1, colour=V2)) +
  geom_density() +
  geom_vline(data=cdat, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(), legend.position=c(0.6,0.8)) +  ylab("") + xlab("") + ggtitle("TM-Align") +
  ylim(0,7.0)+
  xlim(0,1)
# SFAMILY

fp_tma_lddt$V2 <- 'FP_TMA'
tp_tma_sfam_lddt$V2 <- 'TP_TMA_SFAM'

tpfp_tma_sfam <- rbind(fp_tma_lddt,tp_tma_sfam_lddt)

# Find the median of each group
cdat <- ddply(tpfp_tma_sfam, "V2", summarise, V1.median=median(V1))
cdat

# Overlaid histograms with medians
DIST_HIST_SFAM2 <- ggplot(tpfp_tma_sfam, aes(x=V1, fill=V2)) +
  geom_histogram(binwidth=.05, alpha=.5, position="identity") +
  geom_vline(data=cdat, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(), legend.position=c(0.6,0.8)) +  ylab("") + xlab("")

# Density plots with medians
DIST_DENS_SFAM2 <- ggplot(tpfp_tma_sfam, aes(x=V1, colour=V2)) +
  geom_density() +
  geom_vline(data=cdat, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(), legend.position=c(0.6,0.8)) +  ylab("") + xlab("") +
  ylim(0,7.0)+
  xlim(0,1)

# FAMILY

fp_tma_lddt$V2 <- 'FP_TMA'
tp_tma_fam_lddt$V2 <- 'TP_TMA_FAM'

tpfp_tma_fam <- rbind(fp_tma_lddt,tp_tma_fam_lddt)

# Find the median of each group
cdat <- ddply(tpfp_tma_fam, "V2", summarise, V1.median=median(V1))
cdat

# Overlaid histograms with medians
DIST_HIST_FAM2 <- ggplot(tpfp_tma_fam, aes(x=V1, fill=V2)) +
  geom_histogram(binwidth=.05, alpha=.5, position="identity") +
  geom_vline(data=cdat, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(), legend.position=c(0.7,0.8)) +  ylab("") + xlab("lddt")

# Density plots with medians
DIST_DENS_FAM2 <- ggplot(tpfp_tma_fam, aes(x=V1, colour=V2)) +
  geom_density() +
  geom_vline(data=cdat, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(), legend.position=c(0.7,0.8)) +  ylab("") + xlab("lddt") +
  ylim(0,7.0)+
  xlim(0,1)


# DALI
# FOLD
fp_dali_lddt$V2 <- 'FP_DALI'
tp_dali_fold_lddt$V2 <- 'TP_DALI_FOLD'

tpfp_dali_fold <- rbind(fp_dali_lddt,tp_dali_fold_lddt)
tpfp_all_fold <- rbind(tp_dali_fold_lddt,tp_fs_fold_lddt)
tpfp_all_fold <- rbind(tpfp_all_fold, tp_tma_fold_lddt)
tpfp_all_fold <- rbind(tpfp_all_fold, fp_fs_lddt)
tpfp_all_fold <- rbind(tpfp_all_fold, fp_tma_lddt)
tpfp_all_fold <- rbind(tpfp_all_fold, fp_dali_lddt)
cdat_all_fold_one_plot <- ddply(tpfp_all_fold, "V2", summarise, V1.median=median(V1))


# Density plots with medians
DIST_DENS_ALL_FOLD <- ggplot(tpfp_all_fold, aes(x=V1, colour=V2)) +
  geom_density(aes(color=V2)) +
  # geom_vline(data=cdat_all_fold_one_plot, aes(xintercept=V1.median,  colour=V2),
  #            linetype="dashed", size=1) +
  theme(legend.title=element_blank(), legend.position=c(0.6,0.8)) +  ylab("") + xlab("")  +   ggtitle("Fold") +
  ylim(0,7.0)+
  xlim(0,1) +
  scale_color_manual(values=group.colors)
show(DIST_DENS_ALL_FOLD)


# Find the median of each group
cdat <- ddply(tpfp_dali_fold, "V2", summarise, V1.median=median(V1))
cdat

# Overlaid histograms with medians
DIST_HIST_FOLD3 <- ggplot(tpfp_dali_fold, aes(x=V1, fill=V2)) +
  geom_histogram(binwidth=.05, alpha=.5, position="identity") +
  geom_vline(data=cdat, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(), legend.position=c(0.6,0.8)) +  ylab("") + xlab("")  +   ggtitle("Dali")


# Density plots with medians
DIST_DENS_FOLD3 <- ggplot(tpfp_dali_fold, aes(x=V1, colour=V2)) +
  geom_density() +
  geom_vline(data=cdat, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(), legend.position=c(0.6,0.8)) +  ylab("") + xlab("")  +   ggtitle("Dali") +
  ylim(0,7.0)+
  xlim(0,1)
DIST_DENS_FOLD3

# SFAMILY

fp_dali_lddt$V2 <- 'FP_DALI'
tp_dali_sfam_lddt$V2 <- 'TP_DALI_SFAM'

tpfp_dali_sfam <- rbind(fp_dali_lddt,tp_dali_sfam_lddt)
tpfp_all_sfam <- rbind(tp_dali_sfam_lddt,tp_fs_sfam_lddt)
tpfp_all_sfam <- rbind(tpfp_all_sfam, tp_tma_sfam_lddt)
tpfp_all_sfam <- rbind(tpfp_all_sfam, fp_fs_lddt)
tpfp_all_sfam <- rbind(tpfp_all_sfam, fp_tma_lddt)
tpfp_all_sfam <- rbind(tpfp_all_sfam, fp_dali_lddt)
cdat_all_sfam_one_plot <- ddply(tpfp_all_sfam, "V2", summarise, V1.median=median(V1))

# Density plots with medians
DIST_DENS_ALL_SFAM <- ggplot(tpfp_all_sfam, aes(x=V1, colour=V2)) +
  geom_density(aes(color=V2)) +
  # geom_vline(data=cdat_all_sfam_one_plot, aes(xintercept=V1.median,  colour=V2),
  #            linetype="dashed", size=1) +
  theme(legend.title=element_blank(), legend.position=c(0.6,0.8)) +  ylab("") + xlab("")  +   ggtitle("SFam") +
  ylim(0,7.0)+
  xlim(0,1) +
  scale_color_manual(values=group.colors) 
show(DIST_DENS_ALL_SFAM)


# Find the median of each group
cdat <- ddply(tpfp_dali_sfam, "V2", summarise, V1.median=median(V1))
cdat

# Overlaid histograms with medians
DIST_HIST_SFAM3 <- ggplot(tpfp_dali_sfam, aes(x=V1, fill=V2)) +
  geom_histogram(binwidth=.05, alpha=.5, position="identity") +
  geom_vline(data=cdat, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(), legend.position=c(0.6,0.8)) +  ylab("") + xlab("")

# Density plots with medians
DIST_DENS_SFAM3 <- ggplot(tpfp_dali_sfam, aes(x=V1, colour=V2)) +
  geom_density() +
  geom_vline(data=cdat, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(), legend.position=c(0.6,0.8)) +  ylab("") + xlab("") +
  ylim(0,7.0)+
  xlim(0,1)

# FAMILY

fp_dali_lddt$V2 <- 'FP_DALI'
tp_dali_fam_lddt$V2 <- 'TP_DALI_FAM'

tpfp_dali_fam <- rbind(fp_dali_lddt,tp_dali_fam_lddt)


tpfp_all_fam <- rbind(tp_dali_fam_lddt,tp_fs_fam_lddt)
tpfp_all_fam <- rbind(tpfp_all_fam, tp_tma_fam_lddt)
tpfp_all_fam <- rbind(tpfp_all_fam, fp_fs_lddt)
tpfp_all_fam <- rbind(tpfp_all_fam, fp_tma_lddt)
tpfp_all_fam <- rbind(tpfp_all_fam, fp_dali_lddt)
cdat_all_fam_one_plot <- ddply(tpfp_all_fam, "V2", summarise, V1.median=median(V1))


# Find the median of each group
cdat <- ddply(tpfp_dali_fam, "V2", summarise, V1.median=median(V1))
cdat

# Overlaid histograms with medians
DIST_HIST_FAM3 <- ggplot(tpfp_dali_fam, aes(x=V1, fill=V2)) +
  geom_histogram(binwidth=.05, alpha=.5, position="identity") +
  geom_vline(data=cdat, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(), legend.position=c(0.7,0.8)) +  ylab("") + xlab("lddt")

# Density plots with medians
DIST_DENS_FAM3 <- ggplot(tpfp_dali_fam, aes(x=V1, colour=V2)) +
  geom_density() +
  geom_vline(data=cdat, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(), legend.position=c(0.7,0.8)) +  ylab("") + xlab("lddt") +
  ylim(0,7.0)+
  xlim(0,1)


# CLE-SW
# FOLD
fp_cle_lddt$V2 <- 'FP_CLE_SW'
tp_cle_fold_lddt$V2 <- 'TP_CLE_SW_FOLD'

tpfp_cle_fold <- rbind(fp_cle_lddt,tp_cle_fold_lddt)
# Find the median of each group
cdat_cle <- ddply(tpfp_cle_fold, "V2", summarise, V1.median=median(V1))
cdat_cle

# Overlaid histograms with medians
DIST_HIST_FOLD4 <- ggplot(tpfp_cle_fold, aes(x=V1, fill=V2)) +
  geom_histogram(binwidth=.05, alpha=.5, position="identity") +
  geom_vline(data=cdat_cle, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(), legend.position=c(0.6,0.8)) +  ylab("") + xlab("")  +   ggtitle("CLE-SW")


# Density plots with medians
DIST_DENS_FOLD4 <- ggplot(tpfp_cle_fold, aes(x=V1, colour=V2)) +
  geom_density() +
  geom_vline(data=cdat_cle, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(), legend.position=c(0.6,0.8)) +  ylab("") + xlab("")  +   ggtitle("CLE-SW") +
  ylim(0,7.0)+
  xlim(0,1)


# SFAMILY

fp_cle_lddt$V2 <- 'FP_CLE_SW'
tp_cle_sfam_lddt$V2 <- 'TP_CLE_SW_SFAM'

tpfp_cle_sfam <- rbind(fp_cle_lddt,tp_cle_sfam_lddt)

# Find the median of each group
cdat_cle_sfam <- ddply(tpfp_cle_sfam, "V2", summarise, V1.median=median(V1))
cdat_cle_sfam

# Overlaid histograms with medians
DIST_HIST_SFAM4 <- ggplot(tpfp_cle_sfam, aes(x=V1, fill=V2)) +
  geom_histogram(binwidth=.05, alpha=.5, position="identity") +
  geom_vline(data=cdat_cle_sfam, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(), legend.position=c(0.6,0.8)) +  ylab("") + xlab("")

# Density plots with medians
DIST_DENS_SFAM4 <- ggplot(tpfp_cle_sfam, aes(x=V1, colour=V2)) +
  geom_density() +
  geom_vline(data=cdat_cle_sfam, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(), legend.position=c(0.6,0.8)) +  ylab("") + xlab("") +
  ylim(0,7.0)+
  xlim(0,1)

# FAMILY

fp_cle_lddt$V2 <- 'FP_CLE_SW'
tp_cle_fam_lddt$V2 <- 'TP_CLE_SW_FAM'

tpfp_cle_fam <- rbind(fp_cle_lddt,tp_cle_fam_lddt)

# Find the median of each group
cdat_cle_fam <- ddply(tpfp_cle_fam, "V2", summarise, V1.median=median(V1))
cdat_cle_fam

# Overlaid histograms with medians
DIST_HIST_FAM4 <- ggplot(tpfp_cle_fam, aes(x=V1, fill=V2)) +
  geom_histogram(binwidth=.05, alpha=.5, position="identity") +
  geom_vline(data=cdat_cle_fam, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(), legend.position=c(0.7,0.8)) +  ylab("") + xlab("lddt")

# Density plots with medians
DIST_DENS_FAM4 <- ggplot(tpfp_cle_fam, aes(x=V1, colour=V2)) +
  geom_density() +
  geom_vline(data=cdat_cle_fam, aes(xintercept=V1.median,  colour=V2),
             linetype="dashed", size=1) +
  theme(legend.title=element_blank(), legend.position=c(0.7,0.8)) +  ylab("") + xlab("lddt") +
  ylim(0,7.0)+
  xlim(0,1)


grid4 <- plot_grid(DIST_DENS_FOLD, DIST_DENS_FOLD2,DIST_DENS_FOLD3, DIST_DENS_FOLD4,
                  DIST_DENS_SFAM,  DIST_DENS_SFAM2,DIST_DENS_SFAM3, DIST_DENS_SFAM4,
                  DIST_DENS_FAM,  DIST_DENS_FAM2, DIST_DENS_FAM3, DIST_DENS_FAM4,
                  vjust = 1, hjust = -1.5, axis = 'l',
                  ncol = 4, nrow = 3)
ggsave("FS_DALI_TMA_CLESW_LDDT_AVERAGE.pdf", grid4,width=47, height=25, units="cm")
show(grid4)