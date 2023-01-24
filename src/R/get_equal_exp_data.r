rm(list = ls())
library(lme4)
library(lmerTest)
library(stringr)
library(dplyr)

data_dir = 'C:/My Drive/Projects/psilodep2/codebase/data/'
out_dir = 'C:/My Drive/Projects/psilodep2/codebase/data/equal expectancy analysis/'

### Load data & set types and reference levels
scores <- read.csv(paste(data_dir, 'scores.csv', sep='/'))
master_predict <- read.csv(paste(data_dir, 'to_predict_template.csv', sep='/')) 

scores$patient_id <- as.character(scores$patient_id)
scores$trt = as.factor(scores$trt)
scores$tp = as.factor(scores$tp)
scores <- within(scores, trt <- relevel(trt, ref='E'))
scores <- within(scores, tp <- relevel(tp, ref='wk0'))

master_predict$patient_id <- as.character(master_predict$patient_id)
master_predict$trt = as.factor(master_predict$trt)
master_predict$tp = as.factor(master_predict$tp)
master_predict <- within(master_predict, trt <- relevel(trt, ref='E'))

generate_eqexp_data <- function() {

  questionnaires = c('hamd', 'bdi', 'madrs', 'qids', 'stait', 'wemwbs')

  predictions_dir = paste(out_dir,'predicted values',sep='/') 
  results_dir = paste(out_dir,'models from predictied',sep='/') 
    
  esc_values = distinct(subset(scores, select = c(patient_id, esc)))$esc
  psi_values = distinct(subset(scores, select = c(patient_id, psi)))$psi
  stdev_exp = sd(c(esc_values, psi_values), na.rm=TRUE)

  for(mean_exp in seq(0,100,3)) {
    print(mean_exp)

      for(iter in 0:32) {
      predict = add_expectancy_values(mean_exp, stdev_exp)

      for (questionnaire in questionnaires) {

        ### Get prediction from RTX model in E arm
        # Get model and relevant data subset
        subscores_E = subset(scores, trt=='E' & scale==questionnaire & ((tp=='wk0')|(tp=='wk6')))
        model_rtx_E = lmer(score~(1|patient_id)+tp*rtx, subscores_E)
        predicts_rtx_E <- subset(predict, trt=='E' & scale==questionnaire)

        # Get rid of unneeded columns
        predicts_rtx_E$avg <- NULL
        predicts_rtx_E$score_chg <- NULL

        # Make predictions
        predicts_rtx_E_wk0 <- subset(predicts_rtx_E, tp=='wk0')
        predicts_rtx_E_wk6 <- subset(predicts_rtx_E, tp=='wk6')
        predicts_rtx_E_wk6$score <- round(predict(model_rtx_E, newdata=predicts_rtx_E_wk6), digits=0)
        rm(subscores_E, model_rtx_E, predicts_rtx_E)

        # Get prediction from RTX model in P arm
        # Get model and relevant data subset
        subscores_P = subset(scores, trt=='P' & scale==questionnaire & ((tp=='wk0')|(tp=='wk6')))
        model_rtx_P = lmer(score~(1|patient_id)+tp, subscores_P) # delete RTX
        predicts_rtx_P <- subset(predict, trt=='P' & scale==questionnaire)

        # Get rid of unneeded columns
        predicts_rtx_P$avg <- NULL
        predicts_rtx_P$score_chg <- NULL

        # Make predictions
        predicts_rtx_P_wk0 <- subset(predicts_rtx_P, tp=='wk0')
        predicts_rtx_P_wk6 <- subset(predicts_rtx_P, tp=='wk6')
        predicts_rtx_P_wk6$score <- round(predict(model_rtx_P, newdata=predicts_rtx_P_wk6), digits=0)
        rm(subscores_P, model_rtx_P, predicts_rtx_P)

        # Save Predictions
        rtx_predict <- rbind(predicts_rtx_E_wk0, predicts_rtx_P_wk0, predicts_rtx_E_wk6, predicts_rtx_P_wk6)
        fname = paste('predictions_scale',toupper(questionnaire),'_rtx',str_pad(mean_exp, width=3, side=c('left'), pad='0'),'_iter',str_pad(iter, width=3, side=c('left'), pad='0'),'.csv', sep='')
        write.csv(rtx_predict, paste(predictions_dir, fname, sep='/'), row.names=FALSE, quote=FALSE)
        rm(predicts_rtx_E_wk0, predicts_rtx_P_wk0, predicts_rtx_E_wk6, predicts_rtx_P_wk6)

        # Generate model results from predictions
        model_fromRTXpreds = lmer(score~(1|patient_id)+tp*trt, rtx_predict)
        fname = paste('results_scale',toupper(questionnaire),'_rtx',str_pad(mean_exp, width=3, side=c('left'), pad='0'),'_iter',str_pad(iter, width=3, side=c('left'), pad='0'),'.csv', sep='')
        write.csv(coef(summary(model_fromRTXpreds)), paste(results_dir, fname, sep='/'), quote=FALSE)
        rm(model_fromRTXpreds)


        ### Get prediction from AVG model in E arm
        # Get model and relevant data subset
        subscores_E = subset(scores, trt=='E' & scale==questionnaire & ((tp=='wk0')|(tp=='wk6')))
        model_avg_E = lmer(score~(1|patient_id)+tp*avg, subscores_E)
        predicts_avg_E <- subset(predict, trt=='E' & scale==questionnaire)

        # Get rid of unneeded columns
        predicts_avg_E$rtx <- NULL
        predicts_avg_E$score_chg <- NULL

        # Make predictions
        predicts_avg_E_wk0 <- subset(predicts_avg_E, tp=='wk0')
        predicts_avg_E_wk6 <- subset(predicts_avg_E, tp=='wk6')
        predicts_avg_E_wk6$score <- round(predict(model_avg_E, newdata=predicts_avg_E_wk6), digits=0)
        rm(subscores_E, model_avg_E, predicts_avg_E)

        # Get prediction from AVG model in P arm
        # Get model and relevant data subset
        subscores_P = subset(scores, trt=='P' & scale==questionnaire & ((tp=='wk0')|(tp=='wk6')))
        model_avg_P = lmer(score~(1|patient_id)+tp, subscores_P) # delete AVG
        predicts_avg_P <- subset(predict, trt=='P' & scale==questionnaire)

        # Get rid of unneeded columns
        predicts_avg_P$rtx <- NULL
        predicts_avg_P$score_chg <- NULL

        # Make predictions
        predicts_avg_P_wk0 <- subset(predicts_avg_P, tp=='wk0')
        predicts_avg_P_wk6 <- subset(predicts_avg_P, tp=='wk6')
        predicts_avg_P_wk6$score <- round(predict(model_avg_P, newdata=predicts_avg_P_wk6), digits=0)
        rm(subscores_P, model_avg_P, predicts_avg_P)

        # Save Predictions
        avg_predict <- rbind(predicts_avg_E_wk0, predicts_avg_P_wk0, predicts_avg_E_wk6, predicts_avg_P_wk6)
        fname = paste('predictions_scale',toupper(questionnaire),'_avg',str_pad(mean_exp, width=3, side=c('left'), pad='0'),'_iter',str_pad(iter, width=3, side=c('left'), pad='0'),'.csv', sep='')
        write.csv(avg_predict, paste(predictions_dir, fname, sep='/'), row.names=FALSE, quote=FALSE)
        rm(predicts_avg_E_wk0, predicts_avg_P_wk0, predicts_avg_E_wk6, predicts_avg_P_wk6)

        # Generate model results from predictions
        model_fromAVGpreds = lmer(score~(1|patient_id)+tp*trt, avg_predict)
        fname = paste('results_scale',toupper(questionnaire),'_avg',str_pad(mean_exp, width=3, side=c('left'), pad='0'),'_iter',str_pad(iter, width=3, side=c('left'), pad='0'),'.csv', sep='')
        write.csv(coef(summary(model_fromAVGpreds)), paste(results_dir, fname, sep='/'), quote=FALSE)
        rm(model_fromAVGpreds)

      }
      rm(predict)
    }
  }
}
add_expectancy_values <- function(mean_exp, stdev_exp) {
  
  to_predict <- master_predict
  
  for (patient_id in unique(master_predict$patient_id)) {
    esc = round(rnorm(1, mean_exp, stdev_exp), digits=0)
    psi = round(rnorm(1, mean_exp, stdev_exp), digits=0)
    to_predict[to_predict$patient_id==patient_id,]$esc <- esc
    to_predict[to_predict$patient_id==patient_id,]$psi <- psi
  }
  
  to_predict$avg <- (to_predict$esc+to_predict$psi)/2
  to_predict$rtx[to_predict$trt=='E'] <- to_predict$esc[to_predict$trt=='E']
  to_predict$rtx[to_predict$trt=='P'] <- to_predict$psi[to_predict$trt=='P']
  
  return(to_predict)
}

generate_eqexp_data ()