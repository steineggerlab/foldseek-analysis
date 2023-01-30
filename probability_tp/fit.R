# fit.R <1> <2> <3>
# 1: csv file containing columns evalue and score
# 2: E-value cut-off
# 3: number of mixtures


library(mixtools)

args = commandArgs(trailingOnly=TRUE)

df <- read.csv(args[1])
df <- df[df$evalue <= as.numeric(args[2]),]

out <- gammamixEM(df$score, k=strtoi(args[3]), verb = FALSE)

write(out$lambda, stdout())
write(out$gamma.pars, stdout())

