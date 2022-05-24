# install.packages("Matrix")
library(Matrix)
# install.packages("tidyverse")
require(tidyverse)
# install.packages("lubridate")
library(lubridate)
# install.packages("tibbletime")
library(tibbletime)
# install.packages("mgcv")
library(mgcv)
# install.packages("network")
library(network)
# install.packages("ergm")
library(ergm)
# install.packages("network")
library(network)
# install.packages("data.table")
library(data.table)
# install.packages("parallel")
library(parallel)
# install.packages("parallelMap")
library(parallelMap)
# install.packages("sf")
library(sf)
# install.packages("ggpubr")
library(ggpubr)
# install.packages("giscoR")
library(giscoR)
library(haven)
library(sp)
# install.packages("rworldmap")
library(rworldmap)
# install.packages("ergm.userterms")
library(ergm.userterms)
library(ergm.patent)
# install.packages("geosphere")
library(geosphere)
# install.packages("Rfast")
library(Rfast)
rm(list=ls()) 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Functions/functions.R")
Rcpp::sourceCpp('Functions/helper_fun.cpp')
# cor(patent_fillings$docdb_family_size,patent_fillings$cit_db_US_5yrs)
gc(full = T,reset = T)
df = fread("../Data/df.csv")
stepsize_months = 12
n_steps_inclusion =2
df$steps = ceiling(df$month_consec/stepsize_months)
# do the pre-processing for each month separately 
time_steps = unique(df$steps)
# as we condition on the past leave out the first step
# only look at the steps after 1997 -> but we condition on the first two years 
time_steps_eff = unique(df$steps[df$year > 1997]) # 
tmp_number_pat = table(df$docdb_family_id)
# Which patents are singletons?
ind_singelton = names(tmp_number_pat)[tmp_number_pat==1]
# Which patents are collaborative patents?
tmp_number_pat = names(tmp_number_pat)[tmp_number_pat>1]
df_singelton = df[docdb_family_id %in% ind_singelton]
df = df[docdb_family_id %in% tmp_number_pat]

equivalent_step_year = df[,.(year = min(year)), by = steps]
write_rds(equivalent_step_year,file = "../Data/equivalent_step_year.rds")

df_sub2 = df[mainarea34 == "Electrical Engineering"]
n_steps_inclusion_past = 5
n_steps_inclusion_actor = n_steps_inclusion
n_steps_inclusion_actors = n_steps_inclusion
mainarea = "Constraint"
continues = T
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
