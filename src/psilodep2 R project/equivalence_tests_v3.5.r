rm(list = ls())
library(lme4)
library(lmerTest)
library(TOSTER)
library(tidyr)
library(parameters)
library(effsize)

load_data <- function(){
  data_dir = 'C:/Users/szb37/My Drive/Projects/psilodep2/codebase/data/'
  scores <- read.csv(paste(data_dir, 'scores.csv', sep='/'))
  scores$patient_id <- as.character(scores$patient_id)
  scores$trt = as.factor(scores$trt)
  scores <- within(scores, trt <- relevel(trt, ref='E'))
  scores$tp = as.factor(scores$tp)
  scores <- within(scores, tp <- relevel(tp, ref='wk0'))
  return(scores)
}

### WITHIN-ARM MODELS
### Get 0.5 SMD equivalence bound value
if (TRUE) {
  
  rm(list = ls()[ls() != "load_data"])
  scores <- load_data()
  #scores$rtx <- scale(scores$rtx)
  
  ### Validation
  df <- scores %>% subset(scale=='hamd' & trt=='P')
  wk0 <- subset(df,tp=='wk0')
  wk6 <- subset(df,tp=='wk6')
  
  cohen.d(wk0$score, wk6$score, na.rm=TRUE)$estimate #2.6
  
  mean(wk0$score, na.rm=TRUE) #  8.32
  mean(wk6$score, na.rm=TRUE) # 19.2
  mean_diff <-  mean(wk0$score, na.rm=TRUE) - mean(wk6$score, na.rm=TRUE)
  pooled_sd <- sqrt((sd(wk0$score, na.rm=TRUE)^2+sd(wk6$score, na.rm=TRUE)^2)/2)
  mean_diff/pooled_sd # should also get 2.6, same as cohen.d function
  
  ### Calculate 0.5 SMD for each scale
  for (measure in c('hamd', 'bdi', 'madrs', 'qids', 'stait', 'wemwbs')) {
    df <- scores %>% subset(scale==measure)
    wk0 <- subset(df,tp=='wk0')
    wk6 <- subset(df,tp=='wk6')
    print(paste(measure,':',
                round(0.5* #How many SMD?
                        sqrt((sd(wk0$score, na.rm=TRUE)^2+sd(wk6$score, na.rm=TRUE)^2)/2),2)))
  }
}

### Calculate SDs
if (TRUE) {
  df <- unique(subset(scores, trt=='P', select = c("patient_id", "rtx")))
  sd(df$rtx, na.rm=TRUE)
  
  df <- unique(subset(scores, trt=='E', select = c("patient_id", "sss")))
  sd(df$sss, na.rm=TRUE)
}

### Eq testing - tp*rtx in P arm
if (TRUE) { 
  
  rm(list = ls()[ls() != "load_data"])
  scores <- load_data()
  scores$rtx <- c((scores$rtx - median(scores$rtx, na.rm=TRUE))/sd(scores$rtx, na.rm=TRUE))
  scores$sss <- c((scores$sss - median(scores$sss, na.rm=TRUE))/sd(scores$sss, na.rm=TRUE))
  
  subscores = subset(scores, scale=='hamd' & trt=='P' & ((tp=='wk0')|(tp=='wk6')))
  model = lmer(score~(1|patient_id)+tp*rtx, subscores) 
  eqb=2.34
  equivalence_test(range = c(-eqb,eqb), model, verbose=TRUE)
  
  subscores = subset(scores, scale=='bdi' & trt=='P' & ((tp=='wk0')|(tp=='wk6')))
  model = lmer(score~(1|patient_id)+tp*rtx, subscores) 
  eqb=4.74
  equivalence_test(range = c(-eqb,eqb), model)
  
  subscores = subset(scores, scale=='madrs' & trt=='P' & ((tp=='wk0')|(tp=='wk6')))
  model = lmer(score~(1|patient_id)+tp*rtx, subscores) 
  eqb=3.85
  equivalence_test(range = c(-eqb,eqb), model)
  
  subscores = subset(scores, scale=='qids' & trt=='P' & ((tp=='wk0')|(tp=='wk6')))
  model = lmer(score~(1|patient_id)+tp*rtx, subscores) 
  eqb=2.57
  equivalence_test(range = c(-eqb,eqb), model)
  
  subscores = subset(scores, scale=='stait' & trt=='P' & ((tp=='wk0')|(tp=='wk6')))
  model = lmer(score~(1|patient_id)+tp*rtx, subscores) 
  eqb=5.05
  equivalence_test(range = c(-eqb,eqb), model)
  
  subscores = subset(scores, scale=='wemwbs' & trt=='P' & ((tp=='wk0')|(tp=='wk6')))
  model = lmer(score~(1|patient_id)+tp*rtx, subscores) 
  eqb=4.49
  equivalence_test(range = c(-eqb,eqb), model)
}

### Eq testing - tp*sss in E arm
if (TRUE) { 
  
  rm(list = ls()[ls() != "load_data"])
  scores <- load_data()
  scores$rtx <- c((scores$rtx - median(scores$rtx, na.rm=TRUE))/sd(scores$rtx, na.rm=TRUE))
  scores$sss <- c((scores$sss - median(scores$sss, na.rm=TRUE))/sd(scores$sss, na.rm=TRUE))
  
  subscores = subset(scores, scale=='hamd' & trt=='E' & ((tp=='wk0')|(tp=='wk6')))
  model = lmer(score~(1|patient_id)+tp*sss, subscores) 
  eqb=2.34
  equivalence_test(range = c(-eqb,eqb), model)
  
  subscores = subset(scores, scale=='bdi' & trt=='E' & ((tp=='wk0')|(tp=='wk6')))
  model = lmer(score~(1|patient_id)+tp*sss, subscores) 
  eqb=4.74
  equivalence_test(range = c(-eqb,eqb), model)
  
  subscores = subset(scores, scale=='madrs' & trt=='E' & ((tp=='wk0')|(tp=='wk6')))
  model = lmer(score~(1|patient_id)+tp*sss, subscores) 
  eqb=3.85
  equivalence_test(range = c(-eqb,eqb), model)
  
  subscores = subset(scores, scale=='qids' & trt=='E' & ((tp=='wk0')|(tp=='wk6')))
  model = lmer(score~(1|patient_id)+tp*sss, subscores) 
  eqb=2.57
  equivalence_test(range = c(-eqb,eqb), model)
  
  subscores = subset(scores, scale=='stait' & trt=='E' & ((tp=='wk0')|(tp=='wk6')))
  model = lmer(score~(1|patient_id)+tp*sss, subscores) 
  eqb=5.05
  equivalence_test(range = c(-eqb,eqb), model)
  
  subscores = subset(scores, scale=='wemwbs' & trt=='E' & ((tp=='wk0')|(tp=='wk6')))
  model = lmer(score~(1|patient_id)+tp*sss, subscores) 
  eqb=4.49
  equivalence_test(range = c(-eqb,eqb), model)
  
}


### BETWEEN-ARMS MODELS
### Get 0.5 SMD equivalence bound value
if (TRUE) {
  
  rm(list = ls()[ls() != "load_data"])
  scores <- load_data()
  scores$rtx <- c((scores$rtx - median(scores$rtx, na.rm=TRUE))/sd(scores$rtx, na.rm=TRUE))
  scores$sss <- c((scores$sss - median(scores$sss, na.rm=TRUE))/sd(scores$sss, na.rm=TRUE)
  
  ### Validation
  df <- scores %>% subset(scale=='hamd' & tp=='wk6')
  trtE <- subset(df,trt=='E')
  trtP <- subset(df,trt=='P')
  
  cohen.d(trtE$score_chg, trtP$score_chg, na.rm=TRUE)$estimate #1
  
  mean(trtE$score_chg, na.rm=TRUE) #  -4.86
  mean(trtP$score_chg, na.rm=TRUE) # -10.75
  
  mean_diff <- mean(trtE$score_chg, na.rm=TRUE) - mean(trtP$score_chg, na.rm=TRUE)
  pooled_sd <- sqrt((sd(trtE$score_chg, na.rm=TRUE)^2+sd(trtP$score_chg, na.rm=TRUE)^2)/2)
  
  mean_diff/pooled_sd
  
  ### Calculate 0.5 SMD for each scale
  for (measure in c('hamd', 'bdi', 'madrs', 'qids', 'stait', 'wemwbs')) {
    df <- scores %>% subset(scale==measure & tp=='wk6')
    trtE <- subset(df,trt=='E')
    trtP <- subset(df,trt=='P')
    print(paste(measure,':',
      round(0.5* #How many SMDs?
        sqrt((sd(trtE$score_chg, na.rm=TRUE)^2+sd(trtP$score_chg, na.rm=TRUE)^2)/2),2)))
  }
}

### Eq testing for between arm models
if (TRUE) { 
  
  rm(list = ls()[ls() != "load_data"])
  scores <- load_data()
  scores$rtx <- c((scores$rtx - median(scores$rtx, na.rm=TRUE))/sd(scores$rtx, na.rm=TRUE))
  scores$sss <- c((scores$sss - median(scores$sss, na.rm=TRUE))/sd(scores$sss, na.rm=TRUE))
                  
  
  subscores <- subset(scores, scale=='hamd' & ((tp=='wk0')|(tp=='wk6')))
  model = lmer(score~(1|patient_id)+tp*trt*rtx, subscores)
  eqb=2.96
  equivalence_test(range = c(-eqb,eqb), model)
  
  subscores <- subset(scores, scale=='bdi' & ((tp=='wk0')|(tp=='wk6')))
  model = lmer(score~(1|patient_id)+tp*trt*rtx, subscores)
  eqb=5.56
  equivalence_test(range = c(-eqb,eqb), model)
  
  subscores <- subset(scores, scale=='madrs' & ((tp=='wk0')|(tp=='wk6')))
  model = lmer(score~(1|patient_id)+tp*trt*rtx, subscores)
  eqb=4.64
  equivalence_test(range = c(-eqb,eqb), model)
  
  subscores <- subset(scores, scale=='qids' & ((tp=='wk0')|(tp=='wk6')))
  model = lmer(score~(1|patient_id)+tp*trt*rtx, subscores)
  eqb=2.93
  equivalence_test(range = c(-eqb,eqb), model)
  
  subscores <- subset(scores, scale=='stait' & ((tp=='wk0')|(tp=='wk6')))
  model = lmer(score~(1|patient_id)+tp*trt*rtx, subscores)
  eqb=6.22
  equivalence_test(range = c(-eqb,eqb), model)
  
  subscores <- subset(scores, scale=='wemwbs' & ((tp=='wk0')|(tp=='wk6')))
  model = lmer(score~(1|patient_id)+tp*trt*rtx, subscores)
  eqb=5.41
  equivalence_test(range = c(-eqb,eqb), model)
}


