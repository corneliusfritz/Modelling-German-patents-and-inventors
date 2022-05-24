library(Matrix)
library(tidyverse)
library(lubridate)
library(tibbletime)
library(mgcv)
library(network)
library(ergm)
library(network)
library(cowplot)
library(scales)
library(data.table)
library(parallel)
library(parallelMap)
library(sf)
library(ggpubr)
library(giscoR)
library(readr)
library(sna)
library(sp)
library(gridExtra)
library(grid)
library(rworldmap)
library(ergm.userterms)
library(ergm.patent)
library(geosphere)
library(Rfast)
library(statnet)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Functions/functions.R")
Rcpp::sourceCpp('Functions/helper_fun.cpp')
library(haven)
library(ggmcmc)


# Minor Descriptive Plots (Figure 1 in the SM)----

dir_path = paste0("../Data/preprocessing_Final_2_5_12/networks/")
data_files = paste0(dir_path,list.files(dir_path))
network_files = data_files
result = list()
steps = sub('.+_(.+)', '\\1', network_files)
steps = unique(as.numeric(substr(steps,start = 1, stop = nchar(steps) - 4)))
step_data = read_rds("../Data/equivalent_step_year.rds")

for(i in 1:length(network_files)){
  year = step_data$year[match(steps[i],step_data$steps)]
  result[[i]] = readRDS(file = network_files[i])  
  names(result)[i] = year
}

density_bp = function(x){
  dat <- as.edgelist.sna(x)
  return((dim(dat)[1]/2)/((length(attr(dat, "vnames"))-attr(dat, "bipartite")) *
                     attr(dat, "bipartite")))
  # return(sum(mat)/prod(dim(mat)))
}
patents_bp = function(x){
  dat <- as.edgelist.sna(x)
  return(attr(dat, "bipartite"))
  # return(sum(mat)/prod(dim(mat)))
}
inventors_bp = function(x){
  dat <- as.edgelist.sna(x)
  return(length(attr(dat, "vnames"))-attr(dat, "bipartite"))
  # return(sum(mat)/prod(dim(mat)))
}
components_per_step = unlist(lapply(result, components, connected = "recursive"))
connectedness_per_step = unlist(lapply(result, connectedness))
density_per_step = unlist(lapply(result, density_bp))
inventors_per_step = unlist(lapply(result, inventors_bp))
patents_per_step = unlist(lapply(result, patents_bp))

plot_data = data.table(steps = as.numeric(names(components_per_step)), 
                       components_per_step,
                       connectedness_per_step,
                       density_per_step = density_per_step*30000000, 
                       inventors_per_step, 
                       patents_per_step)
plot_data = plot_data[steps>= 10]
plot_data = melt.data.table(plot_data,id.vars = "steps")

options(scipen = 10000)
ggplot(plot_data[variable %in% c("patents_per_step","inventors_per_step","density_per_step")]) +
  geom_line(aes(x = steps , y = value, lty = variable)) + 
  scale_y_continuous(
    name = "Count of inventors and patents", breaks = c(3000,5000,10000,15000),limits = c(3000,15000),
    sec.axis = sec_axis( trans=~./30000000, name="Density")
  ) + 
  xlab("Years") +
  theme_pubr(base_size = 20) +
  scale_linetype_discrete("",label = c("Density", "Active inventors", "Patents")) 

ggsave(filename = "../Plots/descriptive.pdf",width = 9,height = 6)

a = ggplot(plot_data[variable %in% c("patents_per_step","inventors_per_step")]) +
  geom_line(aes(x = steps, y = (value), lty = (variable)), size = 1) + 
  scale_y_continuous(
    name = "Count of inventors and patents", breaks = c(3000,5000,10000,15000),limits = c(3000,15000)
  ) + 
  xlab("Years") +
  theme_pubr(base_size = 20) +
  scale_linetype_discrete("",label = c("Active inventors", "Patents")) 
b = ggplot(plot_data[variable %in% c("density_per_step")]) +
  geom_line(aes(x = steps , y =value/30000000, lty = (variable)), size = 1) + 
  scale_y_continuous(
    name = "Density"
    ) + 
  xlab("Years") +
  theme_pubr(base_size = 20)  +
  scale_linetype_discrete("",label = c(" ")) 
ggsave(a,filename = "../Plots/descriptive_tmpa.pdf",width = 5,height = 5)
ggsave(b,filename = "../Plots/descriptive_tmpb.pdf",width = 5,height = 5)

ggarrange(a,b,ncol = 2)
ggsave(filename = "../Plots/descriptive_tmp.pdf",width = 10,height = 5)


# Model Plots----

dir_path = paste0("../Data/preprocessing_Final_2_5_12/covariates/")
model_files = paste0(dir_path,list.files(dir_path))
dir_path = paste0("../Data/preprocessing_Final_2_5_12/estimates/")
data_files = paste0(dir_path,list.files(dir_path))
dir_path = paste0("../Data/preprocessing_Final_2_5_12/sampled_networks/")
sampled_data_files = paste0(dir_path,list.files(dir_path))
steps = sub('.+_(.+)', '\\1', data_files)
steps = unique(as.numeric(substr(steps,start = 1, stop = nchar(steps) - 4)))
estimators = list()
ses = list()

# Do the simulations for the gof assessment ----
for(i in 1:length(model_files)){
  model_tmp = readRDS(model_files[i])
  # tmp_networks = simulate(model_tmp,nsim = 100,sequential =T,seed = 123,control = control.simulate.formula(MCMC.burnin = 1000,MCMC.interval = 5,MCMC.maxedges = 12000,MCMC.runtime.traceplot = T))
  simulated_networks = simulate(model_tmp,nsim = 200,sequential =F,seed = 123)
  gc(reset = T,full = T)
  write_rds(x = simulated_networks,
            file = paste0("../Data/preprocessing_Final_2_5_12/sampled_networks_for_gof/networks_",steps[i],".rds"))

  # write_rds(x = simulated_networks,
  #           file = paste0("../Data/preprocessing_Electrical_Engineering_2_5_12/sampled_networks_for_gof/networks_",steps[i],".rds"))
  rm(simulated_networks)
  gc(reset = T,full = T)
  cat(i," Step finished \n")
}
dir_path = paste0("../Data/preprocessing_Final_2_5_12/covariates/")
models = paste0(dir_path,list.files(dir_path))
dir_path = paste0("../Data/preprocessing_Final_2_5_12/sampled_networks_for_gof/")
sampled_networks = paste0(dir_path,list.files(dir_path))
# Plot the MCMC diagnostics ----
for(i in 1:length(sampled_networks)){
  year = step_data$year[match(steps[i],step_data$steps)]
  cat("Starting year",year,"\n")
  tmp_model = readr::read_rds(models[[i]])
  
  # tmp_model$sample[[1]] = tmp_model$sample[[1]][,-c(14)]
  tmp_model$sample[[1]] = tmp_model$sample[[1]][,-c(12,13)]
  
  attr(tmp_model$sample[[1]], "dimnames")[[2]] = c("Intercept onset", "Intercept continuation", 
                                                   "Inventor two-star", 
                                                   "Gender male","Female homophily","Male homophily",  
                                                 "Seniority", 
                                                 "Patent two-star", 
                                                 "Team persistence", "Collaboration interlocking", 
                                                 "Spatial proximity"
  )

  
  S <- ggs(tmp_model$sample)
  # Plot the densities of the samples 
  density_1 = ggs_density(S[S$Parameter %in% c("Intercept onset", "Intercept continuation", "Inventor two-star", 
                                   "Gender male","Patent two-star"),]) + theme_pubr() + ylab("Density")
  density_2 = ggs_density(S[S$Parameter %in% c("Seniority", 
                                               "Team persistence", "Collaboration interlocking", 
                                               "Spatial proximity"),]) + theme_pubr() + ylab("Density")
  # Plot the chains of the samples 
  trace_1 = ggs_traceplot(S[S$Parameter %in% c("Intercept onset", "Intercept continuation", "Inventor two-star", 
                                               "Gender male","Patent two-star"),]) + theme_pubr() + ylab("Trace plot")
  trace_2 = ggs_traceplot(S[S$Parameter %in% c("Seniority", 
                                               "Team persistence", "Collaboration interlocking", 
                                               "Spatial proximity"),]) + theme_pubr()+ ylab("Trace plot")

  
  res_1 = grid.arrange(trace_1,density_1, ncol = 2, 
                     top=textGrob( paste("Year ",year), gp=gpar(fontsize=20,font=8)))
  res_2 = grid.arrange(trace_2,density_2, ncol = 2, 
                       top=textGrob( paste("Year ",year), gp=gpar(fontsize=20,font=8)))
  
  ggsave(res_1,filename = paste0("../Plots/mcmc_diag/fin_diag_1_",year,".pdf"),width = 7,height = 10)
  ggsave(res_2,filename = paste0("../Plots/mcmc_diag/fin_diag_2_",year,".pdf"),width = 7,height = 10)
  cat("Finished year",year,"\n")
  rm(tmp_model)
  gc(reset = T,full = T)
  }

# Generate the GOF models (Figure 2 in the SM) ----
for(i in 1:length(sampled_networks)){
  year = step_data$year[match(steps[i],step_data$steps)]
  cat("Starting year",year,"\n")
  tmp_networks = readr::read_rds(sampled_networks[[i]])
  tmp_model = readr::read_rds(model_files[[i]])

  degreedist_bipartite = function(network,vers = 1) {
    n1 = network$gal$bipartite
    n2 = network$gal$n -n1
    if(vers == 1) {
      degb1 = table(unlist(lapply(network$oel[1:n1],length)))
      degb2 = table(unlist(lapply(network$iel[(n1 + 1):(n1 + n2)],length)))
      # degb2[2] = degb2[2] + degb2[1]
      # degb2[1] = 0
    } else  {
      degb1 = table(unlist(lapply(network$iel[1:n1],length)))
      degb2 = table(unlist(lapply(network$oel[(n1 + 1):(n1 + n2)],length)))
    }

    return(list(degb1 = degb1, degb2 = degb2))
  }
  # Do it for the degree of the inventors
  
  trying = lapply(tmp_networks,degreedist_bipartite)
  degb1_data = rbindlist(lapply(trying, function(x) {
    tmp = data.table(x$degb1)
    max =  1:60
    missing = max[!max %in% as.numeric(tmp$V1)]
    tmp = rbind(tmp, data.table(V1 = missing,N = 0))
    return(tmp)
  }))
  
  min = degb1_data[,.(min = quantile(N,probs = c(0))), by = V1]
  
  max = degb1_data[,.(max = quantile(N,probs = c(1))), by = V1]
  mean = degb1_data[,.(mean = mean(N)), by = V1]
  
  data_plot = cbind(min,max$max,mean$mean)
  
  names(data_plot) = c("degree", "min", "max", "mean")
  obs_degb1 =   degreedist_bipartite(tmp_model$network,vers = 2)$degb1
  max_sim = max(as.numeric(data_plot$degree))
  max_obs = max(names(obs_degb1))
  unique_obs = unique(names(obs_degb1))
  
  data_plot$obs =  obs_degb1[match(data_plot$degree, names(obs_degb1))]
  data_plot$obs[is.na(data_plot$obs)] = 0
  data_plot$degree = as.numeric(data_plot$degree)
  data_plot = data_plot[degree %in% 1:min(data_plot$degree[data_plot$max == 0])]
  breaks = c(2,5,10,20,50,100,500,seq(0,max(data_plot$obs),by = 1000))
  
  degreeb1 = ggplot(data_plot, aes(x = degree, y = (obs))) +
    geom_ribbon(aes(ymin = (min), ymax = (max)), fill = "grey") +
    geom_line(aes(y = (mean))) +
    geom_point(color = "red") +
    theme_pubr(base_size = 17) +
    xlab("Degree") +
    ylab("Number of patents") +
    scale_y_continuous(trans = log1p_trans(), breaks = breaks) 
  
  # Do it for the degree of the patents
  degb2_data = rbindlist(lapply(trying, function(x) {
    tmp = data.table(x$degb2)
    max =  1:40
    missing = max[!max %in% as.numeric(tmp$V1)]
    tmp = rbind(tmp, data.table(V1 = missing,N = 0))
    return(tmp)
  }))
  
  min = degb2_data[,.(min = quantile(N,probs = c(0))), by = V1]
  max = degb2_data[,.(max = quantile(N,probs = c(1))), by = V1]
  mean = degb2_data[,.(mean = mean(N)), by = V1]
  
  data_plot = cbind(min,max$max,mean$mean)
  
  names(data_plot) = c("degree", "min", "max", "mean")
  obs_degb2 =   degreedist_bipartite(tmp_model$network,vers = 2)$degb2
  max_sim = max(as.numeric(data_plot$degree))
  max_obs = max(names(obs_degb2))
  unique_obs = unique(names(obs_degb2))
  
  data_plot$obs =  obs_degb2[match(data_plot$degree, names(obs_degb2))]
  data_plot$obs[is.na(data_plot$obs)] = 0
  data_plot$degree = as.numeric(data_plot$degree)
  data_plot = data_plot[degree >1]
  data_plot = data_plot[degree %in% 2:min(data_plot$degree[data_plot$max == 0])]
  breaks = c(2,5,10,20,50,100,500,seq(0,max(data_plot$obs),by = 1000))
  
  degreeb2 = ggplot(data_plot, aes(x = degree, y = (obs))) +
    geom_ribbon(aes(ymin = (min), ymax = (max)), fill = "grey") +
    geom_line(aes(y = (mean))) +
    geom_point(color = "red") +
    theme_pubr(base_size = 17) +
    xlab("Degree") +
    ylab("Number of inventors") +
    scale_y_continuous(trans = log1p_trans(), breaks = breaks) 
  
  res = grid.arrange(degreeb1,degreeb2, ncol = 2, 
               top=textGrob( paste("Year ",steps[i] +1990), gp=gpar(fontsize=20,font=8)))
  ggsave(res,filename = paste0("../Plots/gof/fin_degree_",steps[i] +1990,".pdf"),width = 10,height = 7)
  cat("Finished year",steps[i] +1990,"\n")
  gc(reset = T,full = T)
}

# Plot the coefficients (Figure 5 and 7)----

for(i in 1:length(data_files)){
  readRDS(sampled_data_files[i])
  ses[[i]] = data.table(t(sqrt(diag(readRDS(sampled_data_files[i])))),steps[i])
  estimators[[i]] = data.table(t(readRDS(data_files[i])),steps[i])
  }


estimators = rbindlist(estimators)
ses = rbindlist(ses)
estimators$`offset(b2mindeg2)` = NULL
estimators$`offset(dyadcov.matrix(log(1/(actor_count_time1 + patient_count_time1)), nrow = actor_count_time1, ncol = patient_count_time1))` = NULL
estimators = melt.data.table(estimators,id.vars = "V2")
ses = melt.data.table(ses,id.vars = "V2")

plot_data_comp = data.table(steps = estimators$V2, coefficients = estimators$value, 
                            upper = estimators$value + qnorm(p = 0.975)*ses$value, 
                            lower = estimators$value - qnorm(p = 0.975)*ses$value, 
                            names = estimators$variable)

plot_data_comp$steps = step_data$year[match(plot_data_comp$steps,step_data$steps)]
plot_data_comp = plot_data_comp[steps>=2000]
unique_names = unique(plot_data_comp$names)
plot_list = list()
max_time = max(plot_data_comp$steps)


names = c("Onset",
          "Continuation" ,
          "Inventor star", 
          "Gender male",
          "Female homophily",
          "Male homophily",
          # "Inventor max fillings",
          "Inventor colab patents",
          "Patent star",
          "Patent-inventor repetition",
          "Patent-inventor common partner", 
          "Patent-inventor distance")

plot_list = list()
max_time = max(plot_data_comp$steps)
comb = list(c(1,2), c(3,8),4,c(5,6),7,9,10,11)
names_comb = c("Propensity to invent",
               "Inventor and patent two-star",
               "Gender male",
               "Gender homophily",
               # "Patent family size", 
               "Seniority", 
               "Team persistence", 
               "Collaboration interlocking",
               "Spatial proximity")

size = 1
i = length(comb)
for(i in 1:length(comb)){
  cat(i)
  plot_data = data.frame("coefficients" = 0, "steps"= 0,
                         "lower"= 0,
                         "upper"= 0, 
                         "names"= 0, "zero_included" = F)
  for(j in comb[[i]]){
    plot_data_tmp = plot_data_comp[names == unique_names[j]]
    plot_data_tmp$zero_included = sum(sign(range(range(plot_data_tmp$lower,na.rm = T),range(plot_data_tmp$upper,na.rm = T)))) == 0
    plot_data = rbind(plot_data, plot_data_tmp)
  }
  plot_data = plot_data[-1,]
  zero_included = min(plot_data_tmp$zero_included) == 1
  
  if(zero_included){
    plot_list[[i]] =  ggplot(data = plot_data, aes(x = as.numeric(steps), 
                                                   y = (coefficients),ymin =(lower), 
                                                   ymax =(upper),
                                                   color = names,
                                                   group = names))+
      geom_hline(yintercept = 0,lty = 2) +
      geom_errorbar(width = 0.2,position =position_dodge(width=0.3),size = size) +
      geom_point(size = size+ 1.5,position = position_dodge(width=0.3)) +
      theme_pubr(base_size =  22) +
      theme(axis.text.x = element_text(hjust = 1, angle = 50)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(names_comb[i]) +
      xlab("Year") +
      ylab(expression(paste(hat(theta)[t]))) +
      scale_x_continuous(breaks = 2:max_time)
    
  } else {
    plot_list[[i]] =   ggplot(data = plot_data, aes(x = as.numeric(steps), 
                                                    y = (coefficients),ymin =(lower), 
                                                    ymax =(upper),
                                                    color = names,
                                                    group = names))+
      geom_errorbar(width = 0.2,position =position_dodge(width=0.3),size = size) +
      geom_point(size = size +1.5,position = position_dodge(width=0.3)) +
      theme_pubr(base_size =  22) +
      theme(axis.text.x = element_text(hjust = 1, angle = 50)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(names_comb[i]) +
      xlab("Year") +
      ylab(expression(paste(hat(theta)[t]))) +
      scale_x_continuous(breaks = 2:max_time)
  }
  if(i == 1){
    plot_list[[i]] = plot_list[[i]]+
      scale_color_manual("",labels = c("Experienced", "New"),values = c("black", "grey55"))
  }
  if(i == 2){
    plot_list[[i]] = plot_list[[i]]+
      ggtitle("Two star") +
      scale_color_manual("",labels = c("Inventor", "Patent"),values = c("black", "grey55"))
  }
  if(i == 3){
    plot_list[[i]] = plot_list[[i]]+
      ggtitle("Gender Male") +
      scale_color_manual("",labels = c(""),values = c("black", "grey55"))
  }
  if(i == 4){
    plot_list[[i]] = plot_list[[i]]+
      ggtitle("Gender homophily") +
      scale_color_manual("",labels = c("Female", "Male"),values = c("grey55","black"))
  }
  if(!i %in% c(1,2,3,4)){
    plot_list[[i]] = plot_list[[i]]+
      scale_color_manual("",labels = c("Inventor", "Patent"),values = c("black", "grey55"))+
      guides(color = "none")
  }
  ggsave(paste0("../Plots/smooth_effect/",names_comb[i],"_smooth.pdf"),
         plot = plot_list[[i]],width = 7,height = 7)
}
gc(reset = T,full = T)

# Function to make two coordinate systems equal to one another 
equal_coords = function(plot1, plot2){
  plot_list_tmp_cut = list(plot1, plot2)
  plot_info1_cut = ggplot_build(plot_list_tmp_cut[[1]])
  range_1_cut = c(min(plot_info1_cut$data[[2]]$ymin,na.rm = T), max(plot_info1_cut$data[[2]]$ymax,na.rm = T))
  plot_info2_cut = ggplot_build(plot_list_tmp_cut[[2]])
  range_2_cut = c(min(plot_info2_cut$data[[2]]$ymin,na.rm = T), max(plot_info2_cut$data[[2]]$ymax,na.rm = T))
  ranges_cut = c(min(range_1_cut, range_2_cut,0),max(range_1_cut, range_2_cut,0))
  plot_list_tmp_cut[[1]] = plot_list_tmp_cut[[1]] + coord_cartesian(ylim = ranges_cut)
  plot_list_tmp_cut[[2]] = plot_list_tmp_cut[[2]] + coord_cartesian(ylim = ranges_cut)
  if(sum(sign(ranges_cut)) == 0) {
    plot_list_tmp_cut[[1]] = plot_list_tmp_cut[[1]] + geom_hline(yintercept = 0,lty = 2) 
    plot_list_tmp_cut[[2]] = plot_list_tmp_cut[[2]] + geom_hline(yintercept = 0,lty = 2) 
  }
  return(plot_list_tmp_cut)
}

plot_a = grid.arrange(plot_list[[1]],plot_list[[2]],ncol=2, nrow=1)
plots = equal_coords(plot_list[[6]],plot_list[[7]])
plot_b = grid.arrange(plots[[1]], plots[[2]],ncol=2, nrow=1)
pdf( paste("../Plots/fin_plots_1.pdf", sep = ""),width = 16,height =15)
grid.arrange(plot_a,
             plot_b, ncol = 1)
dev.off()
plots = equal_coords(plot_list[[8]],plots[[1]])
plot_c = grid.arrange(plot_list[[8]], plot_list[[5]],ncol=2, nrow=1)
plots = equal_coords(plot_list[[3]],plot_list[[4]])
plot_d =  grid.arrange(plots[[1]], plots[[2]],ncol=2, nrow=1)
pdf( paste("../Plots/fin_plots_2.pdf", sep = ""),width = 16,height =15)
grid.arrange(plot_c,
             plot_d, ncol = 1)
dev.off()
plots = equal_coords(plot_list[[9]],plot_list[[10]])
plot_e = grid.arrange(plots[[1]], plots[[2]],ncol=2, nrow=1)
# plots = equal_coords(plot_list[[10]],plot_list[[11]])
plots = list(plot_list[[11]],plot_list[[12]])
plot_f = grid.arrange(plots[[1]], plots[[2]],ncol=2, nrow=1)
pdf( paste("../Plots/fin_plots_3.pdf", sep = ""),width = 16,height =15)
grid.arrange(plot_e,
             plot_f, ncol = 1)
dev.off()

gc(full = T)
