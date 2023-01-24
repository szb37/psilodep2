rm(list = ls())
library(infotheo)

data_dir = 'C:/My Drive/Projects/psilodep2/codebase/data/'
scores <- read.csv(paste(data_dir, 'scores.csv', sep='/'))

scores$patient_id <- as.character(scores$patient_id)
scores$trt = as.factor(scores$trt)
scores <- within(scores, trt <- relevel(trt, ref='E'))
scores$tp = as.factor(scores$tp)
scores <- within(scores, tp <- relevel(tp, ref='wk0'))
scores <- scores[c("patient_id", "tp", "scale", "trt", "score_chg", "rtx")]

for (outcome in c('hamd', 'bdi', 'madrs', 'qids', 'stait', 'wemwbs')) {
  for (arm in c('E', 'P')) {
    tmp <- subset(scores, (tp=='wk6') & (scale==outcome) & (trt==arm))
    tmp <- tmp[complete.cases(tmp), ]
    rtx <- tmp$rtx
    delta <- tmp$score_chg
    mi <- as.character(round(natstobits(mutinformation(rtx, delta)), 3))
    print(sprintf("MI on the % s scale in the % s arm is % s", outcome, arm, mi))
  }
}