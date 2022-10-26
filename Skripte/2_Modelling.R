# install.packages("Matrix")
library(Matrix)
# install.packages("tidyverse")
require(tidyverse)
# install.packages("lubridate")
library(lubridate)
# install.packages("ergm")
library(ergm)
# install.packages("network")
library(network)
# install.packages("data.table")
library(data.table)
# install.packages("sf")
library(sf)
# IMPORTANT: Before being able to make this code work you need to install the package ergm.patent 
# (the respective tar.gz file is saved in the replication material)
library(ergm.patent)
rm(list=ls()) 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Functions/functions.R")
Rcpp::sourceCpp('Functions/helper_fun.cpp')
# cor(patent_fillings$docdb_family_size,patent_fillings$cit_db_US_5yrs)
gc(full = T,reset = T)
# df includes all filed patents 
df = fread("../Data/df.csv")
# Decide on the stepsize in months -> in this application we analyse yearly data hence the stepsize is 12 months
stepsize_months = 12
# How long is an actor included in the analysis after not filling a patent? -> in our analysis 2 steps (which is equivalent to 2 years)
n_steps_inclusion =2
df$steps = ceiling(df$month_consec/stepsize_months)
# do the pre-processing for each month separately 
time_steps = unique(df$steps)
# as we condition on the past leave out the first step
time_steps_eff = unique(df$steps[df$year > 1999]) # 
tmp_number_pat = table(df$docdb_family_id)
# Which patents are singletons?
ind_singelton = names(tmp_number_pat)[tmp_number_pat==1]
# Which patents are collaborative patents?
tmp_number_pat = names(tmp_number_pat)[tmp_number_pat>1]
df_singelton = df[docdb_family_id %in% ind_singelton]
df = df[docdb_family_id %in% tmp_number_pat]
equivalent_step_year = df[,.(year = min(year)), by = steps]
write_rds(equivalent_step_year,file = "../Data/equivalent_step_year.rds")
# In the following analyses we concentrate on the subfield of electrical engineering
df_sub2 = df[mainarea34 == "Electrical Engineering"]
# What how far back into the past do we go for calculating the covariates? -> in our paper we go back 5 years -> 5 steps 
n_steps_inclusion_past = 5
n_steps_inclusion_actor = n_steps_inclusion
n_steps_inclusion_actors = n_steps_inclusion
continues = T
# We incorporate two offset terms in our analysis
# 1. An offset term to guarantee that no patent has degree zero -> then the patent would not have been filled in the first place and could not be in our data
#    For this offset we include a statistic counting the patents with degree larger or equal to two (b2mindegree) and set it to Inf
# 2. An offset term to control for the changing size of networks according to Krivitsky et al. (detailed in our manuscript), which we set to 1
offset_coefs = c(Inf, 1)

model_formula = net1 ~ b1factor("onset") + 
  b1factor("continuation") + 
  b1star(2)+
  b1factor("gender") +
  b1nodematch("gender",diff=T) +
  b1factor("number_patents_fac") + 
  b2star(2)+
  two_star_avr(repetition_matrix, "repetition_mean", nrow(repetition_matrix)) +
  two_star_avr(common_partner_matrix, "common_mean", nrow(common_partner_matrix)) +
  two_star_avr(distance_matrix_bin, "dist_mean", nrow(distance_matrix_bin)) +
  offset(b2mindegree(2)) +
  offset(dyadcov(matrix(log(1/(actor_count_time1 + patient_count_time1)),
                        nrow = actor_count_time1,ncol = patient_count_time1)))
number_changes = 15000

# Note that the function calculates more covariates than detailed in the paper  
# Explanation of each provided parameter: 
# df: includes the data.table of all patents 
# df_singelton: includes a data.table of all singleton patents (we use it to derive covariates that were not used in the final model)
# n_steps_inclusion_actor, n_steps_inclusion_past: detailed above
# time_steps_eff: time steps that we effectively model 
# stepsize_months: detailed above
# mainarea: name of the folder to be generated
# continues: logical value to indicate whether the model should continue where a prior run left of or not
# offset_coefs, model_formula: detailed above 
# numer_changes: the number of proposals in the MCMC algorithm 
# seed: numerical value to seed the random number generators (for exact reproducability)
# Needed time around 3-4 hours
full_estimation(df = df,df_singelton = df_singelton, 
           n_steps_inclusion_actor = n_steps_inclusion,
           n_steps_inclusion_past = 5,
           time_steps_eff = time_steps_eff,
           stepsize_months = stepsize_months, 
           mainarea = "Final", 
           continues = F, bs = F, 
           offset_coefs = c(Inf, 1),
           model_formula = model_formula, 
           number_changes = number_changes, seed = 123)

corr_data = numeric(length = length(time_steps_eff))
n = 1
for(i in time_steps_eff){
  corr_data[n] = prepare_data(i = i,df = df,actor_data = actor_data) 
  n = n+1
}
