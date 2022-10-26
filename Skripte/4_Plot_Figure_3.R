library(Matrix)
library(tidyverse)
library(lubridate)
library(tibbletime)
library(mgcv)
library(network)
library(ergm)
library(network)
library(data.table)
library(parallel)
library(sf)
library(ggpubr)
library(giscoR)
library(sp)
library(rworldmap)
rm(list=ls()) 
gc(full = T,reset = T)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Functions/functions.R")

df = fread("../Data/df.csv")
df = df[mainarea34 == "Electrical Engineering" & date > ymd("2000-01-01")]
# df = df[mainarea34 == "Instruments"]
length(unique(df$docdb_family_id))
length(unique(df$person_id))

# we could also only look at patents with a specific number of inventors 
tmp_number_pat = table(df$docdb_family_id)
tmp_number_pat = names(tmp_number_pat)[tmp_number_pat>1]
df = df[docdb_family_id %in% tmp_number_pat]
step_data = read_rds("../Data/equivalent_step_year.rds")
df$year = year(df$date)
df$month = month(df$date)
df$month_consec = df$month + 12* (df$year - min(df$year))
df$date_consec = df$date -min(df$date) +1
df$steps = step_data$steps[match(df$year,step_data$year)]

# Get aggregate data per actor -> from when to when was he/she active, what number of patents were submittet ect.
actor_data = df[,.(start = min(month_consec), 
                   end = max(month_consec), 
                   start_date = min(date), 
                   end_date = max(date), 
                   month_list = list(month_consec),
                   date_list = list(date_consec),
                   number_pat = length(date),
                   length= max(month_consec) - min(month_consec)), 
                by = person_id]
# If the actor was only active at one time point exlucde it
actor_data[start == end]$end = NA
stepsize_months = 12*5
actor_data$steps = ceiling(actor_data$start/stepsize_months)
actor_data$duration = actor_data$end - actor_data$start
actor_data$duration_cens = actor_data$duration
actor_data$duration_cens = actor_data$duration_cens  +1
actor_data$cens = actor_data$duration > 12*5
actor_data$ind = !actor_data$cens
actor_data$duration_cens[actor_data$cens] = 12*5

patent_data = df[,.(date = min(date), 
                    people_included = length(month_consec), 
                    area = area34[1]), 
                 by = docdb_family_id]

patent_data$docdb_family_id_num = 1:length(patent_data$docdb_family_id)
actor_data$person_id_num = 1:length(actor_data$person_id)
df$person_id_num = match(df$person_id,actor_data$person_id)
df$docdb_family_id_num = match(df$docdb_family_id,patent_data$docdb_family_id)
first_patent_per_year = df[,.(first_patent = min(docdb_family_id_num)),by = year]
first_patent_per_year$y_val = max(df$person_id_num) + 1000
first_patent_per_year$first_patent_average = first_patent_per_year$first_patent+0.5*diff(first_patent_per_year$first_patent)
# Figure 3 from the manuscript
png("../Plots/patent_descriptive_full.png",width = 1000,height = 1000)
ggplot() + 
  geom_vline(xintercept = first_patent_per_year$first_patent, color = "grey") +
  geom_point(data=df, mapping=aes(y=person_id_num, x = docdb_family_id_num),
             color = "black", cex = 0.0001, alpha = 0.1) +
  geom_text(data = first_patent_per_year, aes(x =first_patent_average,y =y_val, label = year ), size = 6.5)  +
  ylab("Inventor ID") +
  xlab("Patent ID") +
  theme_pubr(base_size = 20)
dev.off()
