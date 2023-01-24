rm(list = ls())
library(lme4)
library(lmerTest)
library(tidyr)

out_dir = 'C:/My Drive/Projects/psilodep2/codebase/export models/tmp/'

### Load data
data_dir = 'C:/MY Drive/Projects/psilodep2/codebase/data/'
scores <- read.csv(paste(data_dir, 'scores.csv', sep='/'))
scores$patient_id <- as.character(scores$patient_id)
scores$trt = as.factor(scores$trt)
scores <- within(scores, trt <- relevel(trt, ref='E'))
scores$tp = as.factor(scores$tp)
scores <- within(scores, tp <- relevel(tp, ref='wk0'))


get_withinarm <- function(df_name='within_arm_exp_models_bonferroni') {
  
  df <- data.frame()
  
  for (expmes in c('avg', 'rtx')) {
    for (measure in c('hamd', 'bdi', 'madrs', 'qids', 'stait', 'wemwbs')) {
      for (arm in c('E', 'P')) {
        
        if (expmes=='rtx'){
          colnames(scores)[colnames(scores) == 'rtx'] <- 'exp'
        }
        if (expmes=='avg'){
          colnames(scores)[colnames(scores) == 'avg'] <- 'exp'
        }
        
        subscores = subset(scores, scale==measure & trt==arm & ((tp=='wk0')|(tp=='wk6')))
        subscores$exp <- subscores$exp - mean(subscores$exp, na.rm=TRUE)
        subscores$sss <- subscores$sss - mean(subscores$sss, na.rm=TRUE)
        subscores$modtas <- subscores$modtas - mean(subscores$modtas, na.rm=TRUE)
        
        # Results from the 3 models below were manually merged for a single table presented in paper
        model = lmer(score~(1|patient_id)+tp*exp, subscores) 
        #model = lmer(score~(1|patient_id)+tp*sss, subscores) 
        #model = lmer(score~(1|patient_id)+tp*modtas, subscores) 
        rm(subscores)
        
        # Format df
        tmp_df <- data.frame(coef(summary(model)))
        tmp_df$exp <- toupper(expmes)
        tmp_df$scale <- toupper(measure)
        tmp_df$trt <- toupper(arm)
        tmp_df$component <- row.names(tmp_df)
        df <- rbind(df, tmp_df)
      }
    }
  }
  
  #fname = paste('raw_stats_',df_name,'.csv', sep='')
  #write.csv(df, paste(out_dir, fname, sep='/'), row.names = FALSE)
  
  ### Format output
  df <- format_df(df, type='within')
  fname = paste(df_name,'csv', sep='.')
  write.csv(df, paste(out_dir, fname, sep='/'), row.names = FALSE)
  
}
get_betweenarms <- function(df_name='tmp_between_arm_bonferroni_v2.4') {
  
  df <- data.frame()
  
  for (expmes in c('avg', 'rtx')) {
    for (measure in c('hamd', 'bdi', 'madrs', 'qids', 'stait', 'wemwbs')) {
      
      if (expmes=='rtx'){
        colnames(scores)[colnames(scores) == 'rtx'] <- 'exp'
      }
      if (expmes=='avg'){
        colnames(scores)[colnames(scores) == 'avg'] <- 'exp'
      }
      
      subdf = subset(scores, scale==measure & ((tp=='wk0')|(tp=='wk6')))
      subdf$exp <- subdf$exp - mean(subdf$exp, na.rm=TRUE) # zero-center exp measure
      model = lmer(score~(1|patient_id)+tp*trt*exp, subdf)
      rm(subdf)
      
      # Format dataframe
      tmp_df <- data.frame(coef(summary(model)))
      tmp_df$exp <- toupper(expmes)
      tmp_df$scale <- toupper(measure)
      tmp_df$component <- row.names(tmp_df)
      df <- rbind(df, tmp_df)
      rm(model)  
      
    }
  }
  
  #fname = paste('raw_stats_',df_name,'.csv', sep='')
  #write.csv(df, paste(out_dir, fname, sep='/'), row.names = FALSE)
  
  ### Format output
  df <- format_df(df, type='between')
  fname = paste(df_name,'csv', sep='.')
  write.csv(df, paste(out_dir, fname, sep='/'), row.names = FALSE)
}

adjust_within_ps <- function(df){
  # https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust
  
  tmp <- df[c('exp', 'trt', 'component', 'scale', 'p')]
  ps <- as.data.frame(pivot_wider(tmp, names_from='scale', values_from='p'))
  #write.csv(ps, paste(out_dir,'within_raw_ps.csv', sep='/'), row.names = FALSE)
  adj_ps <- ps
  
  for(k in 1:nrow(ps)) {      
    
    if (ps[k,2]=='Intercept') {
      next
    }
    
    adj_ps[k, 4:9] <- p.adjust(ps[k, 4:9], method='bonferroni')
  }
  #write.csv(adj_ps, paste(out_dir,'within_adjusted_ps.csv', sep='/'), row.names = FALSE)  
  
  ### Add adjusted p-values to original dataframe
  ncol_trt = which(colnames(adj_ps)=='trt')
  ncol_exp = which(colnames(adj_ps)=='exp')
  ncol_comp = which(colnames(adj_ps)=='component')
  
  for(k in 1:nrow(adj_ps)) {
    for(scale in c('HAMD', 'BDI', 'MADRS', 'QIDS', 'STAIT', 'WEMWBS')) {
      
      trt = adj_ps[k, ncol_trt]
      exp = adj_ps[k, ncol_exp]
      comp = adj_ps[k, ncol_comp]
      adj_p = adj_ps[k, which(colnames(adj_ps)==scale)]
      
      df[df$trt==trt & df$exp==exp & df$comp==comp & df$scale==scale, ]$adj_p <- adj_p  
      
    }
  }
  return(df)
}
adjust_between_ps <- function(df){
  # https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust
  
  tmp <- df[c('exp', 'component', 'scale', 'p')]
  ps <- as.data.frame(pivot_wider(tmp, names_from='scale', values_from='p'))
  #write.csv(ps, paste(out_dir,'between_raw_ps.csv', sep='/'), row.names = FALSE)
  adj_ps <- ps
  
  for(k in 1:nrow(ps)) {      
    
    if (ps[k,2]=='Intercept') {
      next
    }
    
    adj_ps[k, 3:8] <- p.adjust(ps[k, 3:8], method='bonferroni')
  }
  #write.csv(adj_ps, paste(out_dir,'between_adjusted_ps.csv', sep='/'), row.names = FALSE)  
  
  ### Add adjusted p-values to original dataframe
  ncol_exp = which(colnames(adj_ps)=='exp')
  ncol_comp = which(colnames(adj_ps)=='component')
  
  for(k in 1:nrow(adj_ps)) {
    for(scale in c('HAMD', 'BDI', 'MADRS', 'QIDS', 'STAIT', 'WEMWBS')) {
      
      exp = adj_ps[k, ncol_exp]
      comp = adj_ps[k, ncol_comp]
      adj_p = adj_ps[k, which(colnames(adj_ps)==scale)]
      
      df[df$exp==exp & df$comp==comp & df$scale==scale, ]$adj_p <- adj_p  
      
    }
  }
  return(df)
}
format_signifacance <- function(df) {
  ncol_sig <- which(colnames(df)=="sig")
  ncol_p <- which(colnames(df)=="p")
  ncol_adj_sig <- which(colnames(df)=="adj_sig")
  ncol_adjp <- which(colnames(df)=="adj_p")
  
  for(i in 1:nrow(df)) {
    
    # Unadjusted p-values
    if (df[i, ncol_p] >= 0.05) {
      df[i, ncol_sig] <- ''
    }
    
    if (0.01 <= df[i, ncol_p] & df[i, ncol_p] < 0.05) {
      df[i, ncol_sig] <- '*'
    }
    
    if (0.001 <= df[i, ncol_p] & df[i, ncol_p] < 0.01) {
      df[i, ncol_sig] <- '**'
    }
    
    if (0.001 > df[i, ncol_p]) {
      df[i, ncol_sig] <- '***'
    }
    
    # Adjusted p-values
    if (df[i, ncol_adjp] >= 0.05) {
      df[i, ncol_adj_sig] <- ''
    }
    
    if (0.01 <= df[i, ncol_adjp] & df[i, ncol_adjp] < 0.05) {
      df[i, ncol_adj_sig] <- '*'
    }
    
    if (0.001 <= df[i, ncol_adjp] & df[i, ncol_adjp] < 0.01) {
      df[i, ncol_adj_sig] <- '**'
    }
    
    if (0.001 > df[i, ncol_adjp]) {
      df[i, ncol_adj_sig] <- '***'
    }
  }
  return(df)
}
format_df <- function(df, type){
  
  df$estimate <- round(df$Estimate, 2)
  df$Estimate <- NULL
  df$SE <- round(df$Std..Error, 2)
  df$Std..Error <- NULL
  df$df <- round(df$df, 2)
  df$t <- round(df$t.value, 2)
  df$t.value <- NULL
  df$p <- df$Pr...t..
  df$Pr...t.. <- NULL
  
  df$sig <- ''
  df$adj_p <- 999 # placeholder
  df$adj_sig <- ''
  
  if(type=='within'){
    df <- adjust_within_ps(df)
  }
  
  if(type=='between'){
    df <- adjust_between_ps(df)
  }
  
  
  df <- format_signifacance(df)
  
  df$p <- round(df$p, 3)
  df['p'][df['p'] == '0'] <- '<0.001'  
  df$adj_p <- round(df$adj_p, 3)
  df['adj_p'][df['adj_p'] == '0'] <- '<0.001'  
  
  if(type=='within'){
    df <- df[, c('exp', 'scale', 'trt', 'component', 'estimate', 'SE', 't', 'df', 'p', 'sig', 'adj_p', 'adj_sig')]  
  }
  
  if(type=='between'){
    df <- df[, c('exp', 'scale', 'component', 'estimate', 'SE', 't', 'df', 'p', 'sig', 'adj_p', 'adj_sig')]  
  }
  
  return(df)  
}

get_withinarm()
get_betweenarms()