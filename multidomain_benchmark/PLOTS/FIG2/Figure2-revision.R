library("ggplot2")
library("scales")
library("RColorBrewer")
library("plyr")
library("magrittr")
library("dplyr")
library("grid")
library("cowplot")
library("here")
library("reshape2")
update_geom_defaults("line", list(size = 1.0))
update_geom_defaults("step", list(size = 1.0))
update_geom_defaults("point", list(size = 4.0, stroke=1.5))
theme_set(theme_cowplot(
  font_size = 13,
  line_size = 0.75
) + theme(
  plot.margin = unit(c(0, 0, 0, 0), "cm")
))

#setwd("/Users/mad/Downloads/foldseek(1)/Plots")
#setwd("~/Dropbox/Soeding/Master/HH_ss_roc/ss_evaluation/All_results_scop")
# setwd("/Users/stephaniekim/CLionProjects/prj-foldseek/lotti_plot/data")
setwd("/Users/charlotte/plots-fs/FIGURE_2_FINAL")

cbbPalette <- c( "#000000",  "#CC79A7", "#0072B2",  "#56B4E9", "#E69F00",   "#D55E00", "#5ecffc", "#00ac77")

#theme_set(theme_cowplot(font_size = 20))

scientific_10 <- function(x) {
  ifelse(x==0, "0", parse(text=gsub("1\\s\\%\\*\\%","", gsub("e\\+", " %*% 10^", scales::scientific_format()(x)))))
}

decimal_simple_zero <- function(x) {
  ifelse(x==1, "1",x)
}

# load data 
foldseek_df <- read.table("foldseek.rocx",header = TRUE)
foldseek_df$METHOD = "Foldseek"
blast3d_df <- read.table("3dblast.result.new.1stfp.rocx",header = TRUE)
blast3d_df$METHOD = "3D-BLAST"
cle_df <- read.table("clesw.rocx",header = TRUE)
cle_df$METHOD = "CLE-SW"
ce_df <- read.table("ce.rocx",header = TRUE)
ce_df$METHOD = "CE"
dali_df <- read.table("dali.rocx",header = TRUE)
dali_df$METHOD = "DALI"
tmalign_fast_df <- read.table("tmfast.rocx",header = TRUE)
tmalign_fast_df$METHOD = "TMalign-fast"
tmalign_df <- read.table("tm.rocx",header = TRUE)
tmalign_df$METHOD = "TMalign"
mmseqs2_df <- read.table("mmseqs.rocx",header = TRUE)
mmseqs2_df$METHOD = "MMseqs2"
evalue_data <- read.csv('SCOP_div_evalue_scatter_data.txt',header = FALSE)
appl_df<-read.table("rdrp.tsv", header=TRUE)

# Coverage
fs_countPA = read.delim("foldseek_lddt_cutoff_final_06TP.txt", header = FALSE, sep=' ')
dali_countPA = read.delim("dali_lddt_cutoff_final_06TP.txt", header = FALSE, sep=' ')
ms_countPA = read.delim("mmseqs_lddt_cutoff_final_06TP.txt", header = FALSE, sep=' ')
tm_countPA = read.delim("tmalign_lddt_cutoff_final_06TP.txt", header = FALSE, sep=' ')
clesw_countPA = read.delim("clesw_lddt_cutoff_final_06TP.txt", header = FALSE, sep=' ')

# alignment quality
df_ali_quali <- read.csv('sens_selec_5hit_06TP.csv',header = TRUE)

scientific_10_evalue <- function(x) {
  ifelse(x<1, parse(text=gsub("1\\s\\%\\*\\%","", gsub("e\\-", " %*% 10^-", scales::scientific_format()(x)))), parse(text=gsub("1\\s\\%\\*\\%","", gsub("e\\+", " %*% 10^", scales::scientific_format()(x)))))
}

# load lower panel data
df_speed_vs_sense <- read.table("foldseek-rawdata.csv", sep = ",", header=TRUE)
df_speed_vs_sense
# group.colors.Evalue <- c(" MMseqs2" = "#EE6677"," 3D-BLAST" = "#D55E00"," Foldseek" = "#AA3377")
# group.colors.APPL <- c("Foldseek" = "#FF0000","TMalign" = "#0072B2", "DALI" = "#000000")

group.colors <- c("DALI" = "#A495AF", "TMalign" = "#3B76AF", "Foldseek" = "#D1382C",
                  "TMalign-fast" ="#ADCDE1", "CE"="#559E3E", "CLE-SW"="#EE8632", 
                  "3D-BLAST" = "#643F95", "MMseqs2" = "#F4C27B" )

group.lines <- c("DALI" = "solid", "TMalign" = "dotdash", "Foldseek" = "solid", 
                 "TMalign-fast" ="dotted", "CE"="solid", "CLE-SW"="solid", 
                 "3D-BLAST" = "solid", "MMseqs2" = "solid")

group.alpha <- c("DALI" = 1, "TMalign" = 1, "Foldseek" = 1,
                 "TMalign-fast" = 1, "CE"= 1, "CLE-SW"= 1, 
                 "3D-BLAST" = 1, "MMseqs2" = 1)

group.shapes <- c("DALI" = 16, "TMalign" = 17, "Foldseek" = 15,
                  "TMalign-fast" = 17, "CE"= 18, "CLE-SW"= 4, 
                  "3D-BLAST" = 3, "MMseqs2" = 5)

group.size <- c("DALI" = 16, "TMalign" = 17, "Foldseek" = 15,
                "TMalign-fast" = 17, "CE"= 18, "CLE-SW"= 4, 
                "3D-BLAST" = 3, "MMseqs2" = 5)

groupSize = length(mmseqs2_df[,1])

# build data frame
df<- rbind(
  foldseek_df,
  blast3d_df,
  cle_df,
  tmalign_fast_df,
  ce_df,
  mmseqs2_df,
  tmalign_df,
  dali_df
)
df$METHOD <- factor(df$METHOD, levels = c("DALI", "TMalign", "Foldseek", "TMalign-fast", "CE", "CLE-SW", "3D-BLAST", "MMseqs2"))
levels(df$METHOD)

head(df)

FAM <- ggplot(df %>% 
                group_by(METHOD) %>% 
                arrange(FAM) %>% 
                dplyr::mutate(rn = row_number())) + 
  geom_step(aes(x=1-rn/groupSize, y=FAM,color=METHOD,linetype=METHOD,alpha=METHOD)) + #, alpha=0.8 for transparency or linetype="dashed"
  scale_color_manual(values=group.colors) +
  scale_linetype_manual(values=group.lines) +
  scale_alpha_manual(values=group.alpha) +
  annotate("text", x = 0.89, y = .95, vjust = -0.5, size=5, label = "Family") +
  ylab("Sensitivity up to the 1st FP") + xlab("Fraction of queries") + 
  scale_size_area() + theme(legend.title=element_blank(), legend.position=c(0.05,0.35))

show(FAM)

FAM_SPEED <- ggplot(df_speed_vs_sense, aes(x=TIME, y=FAM, group=METHOD)) +
  geom_point(aes(color=METHOD, shape=METHOD)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  scale_x_log10(breaks=c(1,10,100,1000,10000,100000,1000000,10000000,100000000), label=scientific_10) +
  scale_color_manual(values=group.colors) +  
  scale_shape_manual(values=group.shapes) + 
  ylab("Average Sensitivity up to the 1st FP") + xlab("Time (s)") + 
  scale_size_area() + theme(legend.title=element_blank(), legend.position=c(0.7,0.4))

SFAM <- ggplot(df %>% 
                 group_by(METHOD) %>% 
                 arrange(SFAM) %>% 
                 mutate(rn = row_number())) + 
  geom_step(aes(x=1-rn/groupSize, y=SFAM,color=METHOD,linetype=METHOD)) +
  scale_color_manual(values=group.colors) +
  scale_linetype_manual(values=group.lines) +
  scale_alpha_manual(values=group.alpha) +
  annotate("text", x = 0.6, y = .95, vjust = -0.5, size=5, label = "Super Family") +
  ylab("") + xlab("Fraction of queries") + 
  scale_size_area() + theme(legend.position = "none")

SFAM_SPEED <- ggplot(df_speed_vs_sense, aes(x=TIME, y=SFAM, group=METHOD)) +
  geom_point(aes(color=METHOD, shape=METHOD)) +
  scale_x_log10(breaks=c(1,10,100,1000,10000,100000,1000000,10000000,100000000), label=scientific_10) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  scale_color_manual(values=group.colors) +  
  scale_shape_manual(values=group.shapes) + 
  ylab("") + xlab("Time (s)") + 
  scale_size_area() + theme(axis.text.y = element_blank(), legend.position="none")


FOLD <- ggplot(df %>% 
                 group_by(METHOD) %>% 
                 arrange(FOLD) %>% 
                 mutate(rn = row_number())) + 
  geom_step(aes(x=1-rn/groupSize, y=FOLD, color=METHOD,linetype=METHOD)) +
  scale_color_manual(values=group.colors) +
  scale_linetype_manual(values=group.lines) +
  scale_alpha_manual(values=group.alpha) +
  annotate("text", x = 0.5, y = .95, vjust = -0.5, size=5, label = "Fold") +
  ylab("") + xlab("Fraction of queries") + 
  scale_size_area() + theme(legend.position = "none")

FOLD_SPEED <- ggplot(df_speed_vs_sense, aes(x=TIME, y=FOLD, group=METHOD)) +
  geom_point(aes(color=METHOD, shape=METHOD)) +
  scale_x_log10(breaks=c(1,10,100,1000,10000,100000,1000000,10000000,100000000), label=scientific_10) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  scale_color_manual(values=group.colors) +  
  scale_shape_manual(values=group.shapes) + 
  ylab("") + xlab("Time (s)") + 
  scale_size_area() + theme(axis.text.y = element_blank(), legend.position="none")


EVALUE <- ggplot(evalue_data[which(evalue_data$V2>0),], aes(x=V1, y=V2, group=V3, color=V3)) + 
  geom_point(aes(color=V3), size=1) +
  scale_x_log10(label=scientific_10_evalue) +
  scale_color_manual(values=group.colors, limits = force) + 
  scale_y_log10(label=scientific_10_evalue) +
  scale_size_area() + theme(legend.title=element_blank(), legend.position=c(0.25,0.8)) + 
  xlab("Expected E-value") + ylab("Observed E-value") +
  geom_abline(slope=1, intercept = 0, color = "#636363", linetype="dashed", size=1)
EVALUE

appl_df$TPCUMSUM <- ave(appl_df$TP, appl_df$METHOD, FUN=cumsum)/10
appl_df$FPCUMSUM <- ave(appl_df$FP, appl_df$METHOD, FUN=cumsum)/10
df_custom_before <- data.frame(METHOD=c("Foldseek", "TMalign", "DALI"), TP=c(0,0,0), FP=c(0,0,0), TPCUMSUM = c(0, 0, 0), FPCUMSUM = c(0, 0, 0))
df_custom_before
df_custom_after <- data.frame(METHOD=c("Foldseek", "TMalign", "DALI"), TP=c(0,0,0), FP=c(0,0,0), TPCUMSUM = c(1, 1, .5), FPCUMSUM = c(1.1, 1.1, 1.1))
appl_df<- rbind(df_custom_before, appl_df)
appl_df<- rbind(appl_df, df_custom_after)

APPL <- ggplot(appl_df, aes(x=FPCUMSUM,y=TPCUMSUM, color=METHOD, group=METHOD)) +
  geom_point(size=1, position=position_dodge(width=0.03)) + 
  geom_line(size=1, linetype="dashed", position=position_dodge(width=0.03)) + 
  scale_color_manual(values=group.colors, limits = force) + 
  xlab("False positive") + 
  ylab("True positive") +
  coord_cartesian(ylim=c(0,1), xlim=c(0,1)) +
  theme(legend.title=element_blank(), legend.position=c(0.25,0.8)) +
  scale_x_continuous(expand = c(0.03, 0.03), breaks=c(0.0,0.2,0.4,0.6,0.8,1.0)) + 
  scale_y_continuous(expand = c(0.03, 0.03), breaks=c(0.0,0.2,0.4,0.6,0.8,1.0)) 


APPL

df_speed_vs_sense 
# remove 3D-BLAST-SW from data frame
df_speed_vs_sense <- df_speed_vs_sense[ !(df_speed_vs_sense$METHOD %in% c("3D-BLAST-SW")), ]
df_speed_vs_sense

df_speed_vs_sense_long <- melt(df_speed_vs_sense, id.vars=c("METHOD", "TIME"))
df_speed_vs_sense_long
SPEED <- ggplot(df_speed_vs_sense_long, aes(y=value, x=TIME, shape=variable, size=variable, color=METHOD, group=c(METHOD))) +
  stat_summary(aes(group=variable, x=TIME, linetype=variable, color=variable), fun = identity, geom="line", size=0.5) +
  #geom_line(linetype="dashed", size=0.5) +
  geom_point(size=3) +
  scale_y_continuous(expand = c(0.05, 0.05), limits = c(0, 1)) +
  scale_x_log10(breaks=c(1,10,100,1000,10000,100000,1000000,10000000,100000000), label=scientific_10) +
  scale_color_manual(values=group.colors) +  
  scale_shape_manual(values=c(15,16,17), labels=c("Family", "Super Family", "Fold")) +
  scale_size_manual(values=c(5,5,6))+
  scale_linetype_manual(values=c("dotted", "longdash", "twodash")) +
  guides(
    size=guide_none(),
    linetype=guide_none(),
    color=guide_none(),  #guide_legend(ncol=2),
    shape=guide_legend(ncol=1)
  ) +
  ylab("Avg. sensitivity up to the 5th FP") + xlab("Time (s)") + 
  coord_flip() +
  theme(legend.title=element_blank(), legend.position=c(0.63,0.52))
SPEED

### MULTI WITH PLDDT CHECK AND AVERAGE LDDT INSTEAD OF PER RESIDUE LDDT
### MULTI WITH PLDDT THRESHOLD
colnames(ms_countPA) <- c("V1", "V3")
colnames(tm_countPA) <- c("V1", "V4")
colnames(dali_countPA) <- c("V1", "V5")
colnames(clesw_countPA) <- c("V1", "V6")
colnames(threedblast_countPA) <- c("V1", "V7")
colnames(ce_countPA) <- c("V1", "V8")

myfulldataPA <- merge(fs_countPA, ms_countPA, by.x='V1', all.x=TRUE)
myfulldataPA <- merge(myfulldataPA, dali_countPA, by.x='V1', all.x=TRUE)
myfulldataPA <- merge(myfulldataPA, tm_countPA, by.x='V1', all.x=TRUE)
myfulldataPA <- merge(myfulldataPA, clesw_countPA, by.x='V1', all.x=TRUE)
myfulldataPA <- merge(myfulldataPA, threedblast_countPA, by.x='V1', all.x=TRUE)
myfulldataPA <- merge(myfulldataPA, ce_countPA, by.x='V1', all.x=TRUE)
myfulldataPA
MULTI_PA <- ggplot(myfulldataPA, aes(x=V1, y = value, color = group.colors)) + 
  scale_size_area() + theme(legend.title=element_blank(), legend.position=c(0.55,0.65)) +
  # scale_size_area() + theme(legend.title=element_blank(),legend.position="none") +
  geom_point(aes( y = V2, col = "Foldseek"), size=2, shape=19, stroke=0.7) +
  geom_point(aes( y = V5, col = "DALI"), size=2, shape=19, stroke=0.7) +
  geom_point(aes(y = V3, col = "MMseqs2"), size=2, shape=19, stroke=0.7) +
  geom_point(aes(y = V4, col = "TMalign"), size=2,shape=19, stroke=0.7) +
  geom_point(aes(y = V6, col = "CLE-SW"), size=2,shape=19, stroke=0.7) +
  geom_line(aes( y = V2, col = "Foldseek"),lwd=1)+
  geom_line(aes( y = V5, col = "DALI"),lwd=1)+
  geom_line(aes(y = V3, col = "MMseqs2"),lwd=1)+
  geom_line(aes(y = V4, col = "TMalign"),lwd=1)+
  geom_line(aes(y = V6, col = "CLE-SW"),lwd=1) +
  scale_color_manual(values=group.colors, limits = force) + 
  ylim(0,1.0) +
  xlim(0,20) +
  xlab("TP hits up to 1st FP") + ylab("Query coverage")

MULTI_PA

df_ali_quali 
# remove 3D-BLAST-SW from data frame
df_ali_quali <- df_ali_quali[ !(df_ali_quali$TOOL %in% c("CE")), ]
df_ali_quali
ALI_Q <- ggplot(df_ali_quali, aes(x=SELECTIVITY, y=SENSITIVITY, group=TOOL)) +
  geom_point(aes(color=TOOL)) +
  scale_color_manual(values=group.colors, limits = force) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
  ylab("Sensitivity") + xlab("Precision") + 
  scale_size_area() +
  theme(legend.position="none")

FAM
SFAM
FOLD
EVALUE
SPEED
MULTI_PA
ALI_Q
APPL

grid <- plot_grid(SPEED, NULL,FAM, NULL,  SFAM, FOLD, EVALUE,NULL,  MULTI_PA, NULL, ALI_Q, APPL, labels=c("a","","b","","c","d","e","f","","g","","h"), ncol =6, nrow =2,
                  vjust = 1, hjust = -1.5, axis = 'l', rel_heights = c(1.0),  rel_widths = c(1, 0, 1, 0, 1,1))
ggsave("./Figure2-updated.pdf", grid, width=39, height=16.4, units="cm")
system("open Figure2-updated.pdf")









