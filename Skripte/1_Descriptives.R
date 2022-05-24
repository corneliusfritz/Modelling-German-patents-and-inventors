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
library(survival)
library(colourvalues)
library(giscoR)
library(sp)
library(rworldmap)
rm(list=ls()) 
gc(full = T,reset = T)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Functions/functions.R")

# This is a script of the descriptives 
# df.sub2 <- df.sub[date %within% interval(ymd("1990-01-01"),ymd("1990-06-01"))]

df = fread("../Data/df.csv")
df = df[mainarea34 == "Electrical Engineering" & date > ymd("2000-01-01")]
# df = df[mainarea34 == "Instruments"]

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

gender_df = df[,.(mean(gender == "M")), by = year]

a = ggplot(data = gender_df)+
  geom_line(aes(y = V1, x = year), size =1) +
  theme_pubr(base_size = 20) +
  xlab("Year") +
  ylab("Percentage of male inventors")
ggsave(a,filename = "../Plots/gender.pdf",width = 5,height = 5)




actor_data = df[,.(start = min(month_consec), 
                   end = max(month_consec), 
                   start_date = min(date), 
                   end_date = max(date), 
                   month_list = list(month_consec),
                   date_list = list(date_consec),
                   number_pat = length(date),
                   length= max(month_consec) - min(month_consec)), 
                by = person_id]

actor_data[start == end]$end = NA
stepsize_months = 12*5
actor_data$steps = ceiling(actor_data$start/stepsize_months)
actor_data$duration = actor_data$end - actor_data$start
actor_data$duration_cens = actor_data$duration
actor_data$duration_cens = actor_data$duration_cens  +1
actor_data$cens = actor_data$duration > 12*5
actor_data$ind = !actor_data$cens
actor_data$duration_cens[actor_data$cens] = 12*5

tmp_function = function(x){
  tmp = x$date_list
  from = tmp[-length(tmp)]
  to = tmp[-1]
  return(data.table(from = from, to = to, status = 1, time_beg = x$start))
}
# debugonce(tmp_function)
tmp = apply(actor_data, MARGIN = 1, FUN = tmp_function)
between_data = rbindlist(tmp)
between_data = between_data[!is.na(from)]
# between_data = between_data[from != to]
between_data$time = between_data$to - between_data$from
between_data$status[between_data$time == 0] = 0
between_data$time[between_data$time == 0] = between_data$time[between_data$time == 0] +1

between_data$steps = ceiling(between_data$time_beg/stepsize_months)
between_data$time_beg = 0

# this takes pretty long -> 10-15 min
surv = survfit(Surv(time = time,event = status, type = "left") ~ 1,  type="kaplan-meier",
               data=between_data)

time_patent = data.table(time = c(0,surv$time/365), 
                         surv= c(1,surv$surv))

# Figure 4
pdf("../Plots/km_estimator_1.pdf",width = 7,height = 5)
ggplot(data = time_patent)+
  geom_hline(yintercept = time_patent$surv[time_patent$time ==2],lty = 2, color = "grey")+ 
  geom_vline(xintercept = 2,lty = 2, color = "grey") + 
  geom_step(aes(y = surv, x = time), size =1) +
  theme_pubr(base_size = 20) +
  xlab("Duration in years between consecutive patents") +
  ylab("Survival probability") 

dev.off()
