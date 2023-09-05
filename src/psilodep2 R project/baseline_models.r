rm(list = ls())
library(lme4)
library(lmerTest)

### Load data
data_dir = 'C:/Users/szb37/My Drive/Projects/psilodep2/codebase/data/'
exp <- read.csv(paste(data_dir, 'expectancy.csv', sep='/'))
bslvars <- read.csv(paste(data_dir, 'baseline_variables.csv', sep='/'))

### Set types and reference levels
exp$patient_id <- as.character(exp$patient_id)
exp$trt <- as.factor(exp$trt)
exp <- within(exp, trt <- relevel(trt, ref='E'))
exp$type <- as.character(exp$type)

bslvars$patient_id <- as.character(bslvars$patient_id)
bslvars$trt <- as.factor(bslvars$trt)
bslvars <- within(bslvars, trt <- relevel(trt, ref='E'))
bslvars$sex <- as.factor(bslvars$sex)
bslvars <- within(bslvars, sex <- relevel(sex, ref='M'))
bslvars$edu <- as.factor(bslvars$edu)
bslvars <- within(bslvars, edu <- relevel(edu, ref='0'))
bslvars$edu <- as.character(bslvars$edu)

### Baseline differences models
summary(lmer(exp ~ (1|patient_id)+type, exp))
summary(lmer(exp ~ (1|patient_id)+type*trt, exp))
summary(lm(sss ~ trt, bslvars))
summary(lm(modtas ~ trt, bslvars))

# MODTAS subscales
summary(lm(modtas_ii ~ trt, bslvars))  # Imaginative involvement
summary(lm(modtas_asc ~ trt, bslvars)) # Altered States of Consciousness
summary(lm(modtas_ain ~ trt, bslvars)) # Aesthetic involvement in nature
summary(lm(modtas_esp ~ trt, bslvars)) # Extra Sensory Perception
summary(lm(modtas_syn ~ trt, bslvars)) # Synesthesia