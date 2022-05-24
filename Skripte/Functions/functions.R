from_el_to_bimatrix = function(el, mode1, mode2) {
  
  bimatrix = spMatrix(nrow = nrow(as.vector(unique(el[,mode1, with = F]))),
                       ncol = nrow(unique(el[,mode2, with = F])),
                       i = as.numeric(factor(t(el[,mode1, with = F]))),
                       j = as.numeric(factor(t(el[,mode2, with = F]))),
                       x = rep(1, length(as.numeric(t(el[,mode1, with = F])))))
  
  rownames(bimatrix) = levels(factor(t(el[,mode1, with = F])))
  colnames(bimatrix) = levels(factor(t(el[,mode2, with = F])))
  bimatrix = as.matrix(bimatrix)
  return(bimatrix)
}

# The single argument to this function, points, is a data.frame in which:
#   - column 1 contains the longitude in degrees
#   - column 2 contains the latitude in degrees
coords2country = function(points)
{  
  countriesSP = getMap(resolution='low')
  #countriesSP = getMap(resolution='high') #you could use high res map from rworldxtra if you were concerned about detail
  
  # convert our list of points to a SpatialPoints object
  
  # pointsSP = SpatialPoints(points, proj4string=CRS(" +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
  
  #setting CRS directly to that from rworldmap
  pointsSP = SpatialPoints(points, proj4string=CRS(proj4string(countriesSP)))  
  
  
  # use 'over' to get indices of the Polygons object containing each point 
  indices = over(pointsSP, countriesSP)
  
  # return the ADMIN names of each country
  indices$ADMIN  
  #indices$ISO3 # returns the ISO3 code 
  #indices$continent   # returns the continent (6 continent model)
  #indices$REGION   # returns the continent (7 continent model)
}

#version of from_el_to_bimatrix, where the actor set includes a specific given set of actors 
from_el_to_bimatrix_alt = function(el, mode1, mode2, actor_set) {
  # match the names with the actors in the actor_set
  i = match(as.character(t(el[,mode1, with = F])), actor_set)
  
  bimatrix = spMatrix(nrow = length(actor_set),
                       ncol = nrow(unique(el[,mode2, with = F])),
                       i = i,
                       j = as.numeric(factor(t(el[,mode2, with = F]))),
                       x = rep(1, length(as.numeric(t(el[,mode1, with = F])))))
  
  rownames(bimatrix) = actor_set
  colnames(bimatrix) = levels(factor(t(el[,mode2, with = F])))
  bimatrix = as.matrix(bimatrix)
  return(bimatrix)
}

prepare_data = function(i,df_sub2,actor_data, n_steps_inclusion = 5, dir_path = NULL) {
  now = Sys.time()
  # define the bimatrix at time point t 
  df_time1 = df_sub2[steps == i]
  # df_time1$area34 = factor(df_time1$area34,levels = unique(df_sub2$area34))
  df_past_present = df_sub2[steps %in% ((i-(n_steps_inclusion -1)):i)]
  # df_past_present$area34 = factor(df_past_present$area34,levels = unique(df_sub2$area34))
  actor_data_tmp = actor_data[steps < i]
  actor_data_tmp = actor_data_tmp[order(steps,decreasing = T)]
  actor_data_tmp = actor_data_tmp[!duplicated(actor_data_tmp$person_id)]
  
  df_past = df_past_present[steps != i]
  # only look at the single patents
  cov1 = df_time1[!duplicated(df_time1[, "docdb_family_id"]), ]
  actors_time1 = unique(df_past_present$person_id)
  df_time1_mat = from_el_to_bimatrix_alt(el = df_time1,mode1 =  "person_id",
                                         mode2 = "docdb_family_id",actor_set = actors_time1)
  patient_count_time1 = length(colnames(df_time1_mat))
  patients_time1 = colnames(df_time1_mat)
  
  # df.time1 = df.time1 %>% arrange(docdb_family_id)
  # cov1 = df.time1[!duplicated(df.time1[, "docdb_family_id"]), ]
  # 
  net1 = network(df_time1_mat, matrix.type = "bipartite", ignore.eval = FALSE)
  # data_tmp = data.table(
  #   from = as.character(df_time1$person_id),
  #   to =as.character(df_time1$docdb_family_id),
  #   stringsAsFactors = FALSE
  # )
  # vertices = data.table(c(patients_time1,actors_time1), 
  #            is_actor = c(rep(FALSE,times = length(patients_time1)), 
  #                         rep(TRUE,times = length(actors_time1))))
  # debugonce(as.network)
  # trying = as.network(data_tmp, bipartite = T, directed = F,vertices = vertices)
  
  # define the bimatrix at time point t  -1 
  df_time0 = df_past[person_id %in% actors_time1]
  df_time0_mat = from_el_to_bimatrix_alt(el = df_time0, mode1 = "person_id",
                                          mode2 =  "docdb_family_id", actor_set = actors_time1)
  net1%v%"area_patent" = c(rep(as.character(NA), nrow(df_time1_mat)), cov1$area34)
  net1%v%"lon" =   c(df_past_present$lon[match(actors_time1, df_past_present$person_id)],rep(NA,patient_count_time1))
  net1%v%"lat" =   c(df_past_present$lat[match(actors_time1, df_past_present$person_id)],rep(NA,patient_count_time1))
  net1%v%"number_patent_lag" = c(rowSums(df_time0_mat),rep(NA,patient_count_time1),rep(NA,patient_count_time1))
  net1%v%"one_star_patent_lag" = c(rowSums(df_time0_mat[,colSums(df_time0_mat) == 1]),rep(NA,patient_count_time1))
  net1%v%"two_star_patent_lag" = c(rowSums(df_time0_mat[,colSums(df_time0_mat) == 2]),rep(NA,patient_count_time1))
  net1%v%"three_star_patent_lag" = c(rowSums(df_time0_mat[,colSums(df_time0_mat) == 3]),rep(NA,patient_count_time1))
  net1%v%"rep_patent_lag" = net1%v%"number_patent_lag" > 0
  tmp_actor_area = actor_data_tmp$area34[match(actors_time1, actor_data_tmp$person_id)]
  tmp_actor_area[is.na(tmp_actor_area)] = "No prior patent"
  net1%v%"area_actor" = c(tmp_actor_area,rep(NA,patient_count_time1))
  
  number_steps_per_actors = df_time0[,.(number_steps = length(unique(steps))), by = person_id]
  tmp_number = c(rep(0,nrow(df_time1_mat)),rep(NA,patient_count_time1))
  tmp_match = match(number_steps_per_actors$person_id,net1%v%"vertex.names")
  tmp_number[tmp_match] = number_steps_per_actors$number_steps
  net1%v%"number_steps" = tmp_number
  
  # get the change statistics
  dta.array = ergmMPLE(net1 ~ gwb1dsp(decay=0, fixed=T) +
                          b1star(2) + 
                          b1star(3) + 
                          # b2nodematch("area_patent", beta = 0, diff = TRUE) +
                          b2factor("area_patent") +
                          b2mindegree(1)+
                          b1cov("number_patent_lag") + 
                          # b1cov("one_star_patent_lag")  + 
                          b1cov("two_star_patent_lag") + 
                          b1cov("three_star_patent_lag") + 
                          # b1nodematch("area_actor", beta = 0,byb2attr = "area_patent") +
                          b1factor("area_actor") +
                          b1cov("lon") + 
                          b1cov("lat") + 
                          b1cov("rep_patent_lag")+ 
                          b1cov("number_steps"), 
                        output = "array",
                        control = control.ergm(MPLE.max.dyad.types = network.dyadcount(net1)*2))
  
  ncoef = length(dta.array$predictor[1,2,])

  # data_full = data.table(matrix(nrow = length(Y),ncol = ncoef))
  # names(data_full) = attr(dta.array$predictor,"dimnames")$term
  data_full = apply(dta.array$predictor,MARGIN = 3, FUN = function(x){
    as.numeric(x)
  })
  data_full = as.data.frame(data_full)
  # setDT(data)
  # head(data)
  # for(p in 1:ncoef) {
  #   data_full[,p] = as.numeric(dta.array$predictor[,,p])
  # }
  # 
  data_full$Y = as.numeric(dta.array$response)
  data_full$inventor = rep(actors_time1,times = patient_count_time1 )
  setDT(data_full)
  
  # data_full$patent = rep(cov1$docdb_family_id ,each = length(actors_time1))
  
  # There is a problem with the fact that some subgroups are not present in some steps  
  # tmp_matrix = data.table(model.matrix(docdb_family_id~ area34, data = cov1))
  # 
  # data_full = cbind(data_full,tmp_matrix[match(data_full$patent,cov1$docdb_family_id),])
  # 
  # data_full$patent =data_full$`(Intercept)` = NULL
  
  # carry out 'weight trick'  
  data_full = data_full[,.(weight = .N),by = names(data_full)]
  data_full$times = i
  cat("Step ",i," finished\n")
  cat("Needed time was:",Sys.time() -now," Seconds \n")
  if(!is.null(dir_path)) {
    save_path = paste0(dir_path,"/Step_",i,".rds")
    saveRDS(data_full,file = save_path)
  }
  
  return(data_full)
}


estimate_list_function = function(x, repetition_matrix,common_partner_matrix, 
                                  distance_matrix_bin, homophily_area_mat, 
                                  actor_count_time1,patient_count_time1 ) {
  mod_tmp =   ergm(x ~ b1factor("onset") + 
                     b1factor("continuation") + 
                     b1star(2)+
                     b1factor("gender") +
                     b1cov("internationalisation_bin") + 
                     b1factor("max_fillings_fac") +
                     b1factor("number_patents_fac") + 
                     b1factor("singeltons_bin")  +
                     b2star(2)+
                     two_star_avr(repetition_matrix, "repetition_mean", nrow(repetition_matrix)) +
                     two_star_avr(common_partner_matrix, "common_mean", nrow(common_partner_matrix)) +
                     # two_star_avr(distance_matrix_bin, "dist_mean", nrow(distance_matrix_bin)) +
                     dyadcov(homophily_area_mat) + 
                     offset(dyadcov(matrix(log(1/(actor_count_time1 + patient_count_time1)),
                                           nrow = actor_count_time1,ncol = patient_count_time1))) + 
                     offset(b2mindegree(2)),estimate = "MPLE",
                   control = control.ergm(MPLE.max.dyad.types = network.dyadcount(net1) *4),
                   offset.coef = c(1,Inf))
  # 
  # 
  # 
  # mod_tmp = ergm(x ~ b1factor("onset") + 
  #                  b1factor("continuation") + 
  #                  b1star(2)+
  #                  b1factor("gender") +
  #                  b1cov("internationalisation_bin") + 
  #                  b1factor("max_fillings_fac") +
  #                  b1factor("number_patents_fac") + 
  #                  b1factor("singeltons_bin")  +
  #                  b2star(2)+
  #                  two_star_avr(repetition_matrix, "repetition_mean", nrow(repetition_matrix)) +
  #                  two_star_avr(common_partner_matrix, "common_mean", nrow(common_partner_matrix)) +
  #                  two_star_avr(distance_matrix_bin, "dist_mean", nrow(distance_matrix_bin)) +
  #                  dyadcov(homophily_area_mat),estimate = "MPLE",verbose = F)
  gc(full = T)
  return(coef(mod_tmp))
}


prepare_data_no_re = function(i,df_sub2,df_singelton = df_singelton,actor_data, bs, offset_coefs = c(1),
                              model_formula, number_changes = 1000, seed, 
                              n_steps_inclusion_actors = 2,n_steps_inclusion_past = 5, dir_path = NULL) {
  now = Sys.time()
  # Step 1: Define networks ----
  # define the bimatrix at time point t 
  df_time1 = df_sub2[steps == i]
  # df_time1$area34 = factor(df_time1$area34,levels = unique(df_sub2$area34))
  df_past_present = df_sub2[steps %in% ((i-(n_steps_inclusion_actors -1)):i)]
  # df_past_present$area34 = factor(df_past_present$area34,levels = unique(df_sub2$area34))
  actor_data_tmp = actor_data[steps < i]
  actor_data_tmp = actor_data_tmp[order(steps,decreasing = T)]
  
  # actor_data_tmp = actor_data_tmp[!duplicated(actor_data_tmp$person_id)]
  
  df_past = df_sub2[steps %in% (((i-1)-(n_steps_inclusion_past -1)):(i-1))]
  df_singelton_past = df_singelton[steps %in% (((i-1)-(n_steps_inclusion_past -1)):(i-1))]
  
  number_of_complete_patents = length(unique(df_past$docdb_family_id)) + length(unique(df_singelton_past$docdb_family_id))
  
  number_actors_time0 = length(unique(df_past$person_id))
  number_patents_time0 = length(unique(df_past$docdb_family_id))
  
  # only look at the single patents
  cov1 = df_time1[!duplicated(df_time1[, "docdb_family_id"]), ]
  actors_time1 = unique(df_past_present$person_id)
  df_time1_mat = from_el_to_bimatrix_alt(el = df_time1,mode1 =  "person_id",
                                         mode2 = "docdb_family_id",actor_set = actors_time1)
  patient_count_time1 = length(colnames(df_time1_mat))
  actor_count_time1 = length(actors_time1)
  patients_time1 = colnames(df_time1_mat)
  
  # Step 2: Define network object and add attributed----
  # df.time1 = df.time1 %>% arrange(docdb_family_id)
  # cov1 = df.time1[!duplicated(df.time1[, "docdb_family_id"]), ]
  # 
  net1 = network(df_time1_mat, matrix.type = "bipartite", ignore.eval = FALSE)
  # data_tmp = data.table(
  #   from = as.character(df_time1$person_id),
  #   to =as.character(df_time1$docdb_family_id),
  #   stringsAsFactors = FALSE
  # )
  # vertices = data.table(c(patients_time1,actors_time1), 
  #            is_actor = c(rep(FALSE,times = length(patients_time1)), 
  #                         rep(TRUE,times = length(actors_time1))))
  # debugonce(as.network)
  # trying = as.network(data_tmp, bipartite = T, directed = F,vertices = vertices)
  
  # Step 2 a) Patent info ----
  df_time0 = df_past[person_id %in% actors_time1]
  df_time0_mat = from_el_to_bimatrix_alt(el = df_time0, mode1 = "person_id",
                                          mode2 =  "docdb_family_id", actor_set = actors_time1)
  cov1$area34 = factor(cov1$area34,levels = c("Audiovisual", "Electr/Energy", "Telecom", "ComputerTech", 
                                                               "BasicCommProcess", "DigitalComm", "Semiconductors", "IT_Methods", 
                                                               "No prior patent"))
  cov1 = cov1[match(net1 %v% "vertex.names",cov1$docdb_family_id)]
  cov1 = cov1[!is.na(steps)]
  
  # order cov1 just like net1 %v% "vertex.names"
  net1%v%"area_patent" = c(rep(NA, nrow(df_time1_mat)), cov1$area34)
  net1%v%"lon" =   c(round(df_past_present$lon[match(actors_time1, df_past_present$person_id)],digits = 2),rep(NA,patient_count_time1))
  net1%v%"lat" =   c(round(df_past_present$lat[match(actors_time1, df_past_present$person_id)],digits = 2),rep(NA,patient_count_time1))
  # plot(match(cov1$docdb_family_id,net1 %v% "vertex.names"))
  
  
  # Step 2 b) Inventor info ----
  # actor info for all actors in the network 
  gender_data = actor_data[person_id %in% actors_time1]
  gender_dt = gender_data[,.(gender = gender[1]),by = person_id]
  gender_info = gender_dt$gender[match(actors_time1, gender_dt$person_id)]
  
  net1%v%"gender" =   c(gender_info,rep(NA,patient_count_time1))
  
  # How many number of steps was each inventor active?
  number_steps_per_actors = df_time0[,.(number_steps = i- min(steps)), by = person_id]
  number_steps = number_steps_per_actors$number_steps[match(actors_time1, number_steps_per_actors$person_id)]
  # everyone who was not active in the past but is now is active 0 steps
  number_steps[is.na(number_steps)] = 0
  
  # tmp_number = c(rep(0,nrow(df_time1_mat)),rep(NA,patient_count_time1))
  # tmp_match = number_steps_per_actors$number_steps[match(actors_time1,number_steps_per_actors$person_id)]
  
  # match(number_steps_per_actors$person_id,net1%v%"vertex.names")
  # tmp_number[tmp_match] = number_steps_per_actors$number_steps
  # match(number_steps_per_actors$person_id,net1%v%"vertex.names")
  # net1%v%"number_steps" = tmp_number
  
  
  # Max and mean number of fillings in the past 
  filling_dt = df_time0[,.(max_filling = max(fillings_number), mean_filling = mean(fillings_number)),by = person_id]
  max_fillings = filling_dt$max_filling[match(actors_time1, filling_dt$person_id)]
  # where there are no fillings in the past the number of fillings is zero
  max_fillings[is.na(max_fillings)] = 0
  mean_fillings = filling_dt$mean_filling[match(actors_time1, filling_dt$person_id)]
  # where there are no fillings in the past the number of fillings is zero
  mean_fillings[is.na(mean_fillings)] = 0
  net1%v%"max_fillings" =   c(max_fillings,rep(NA,patient_count_time1))
  # 
  # max_fillings_1 = max_fillings == 1
  # max_fillings_2 = max_fillings %in% c(2:3)
  # max_fillings_3 = max_fillings > 4
  # net1%v%"max_fillings_1" =   c(max_fillings_1,rep(NA,patient_count_time1))
  # net1%v%"max_fillings_2" =   c(max_fillings_2,rep(NA,patient_count_time1))
  # net1%v%"max_fillings_3" =   c(max_fillings_3,rep(NA,patient_count_time1))
  #
  quantiles_fillings = quantile(max_fillings[max_fillings>0],probs = c(0,0.5))
  # if 1 = 2 then we set 2 equal to the next unique value 
  tmp1 = quantiles_fillings[1] == quantiles_fillings[2]
  if(tmp1){
    unique_vals = sort(unique(max_fillings[max_fillings>0]))
    quantiles_fillings[2] = unique_vals[min(which(unique_vals>quantiles_fillings[2]))]
  }
 
  # tmp = which(quantiles_fillings[-length(quantiles_fillings)] == quantiles_fillings[-1])
  # if(length(tmp)>0){
  #   quantiles_fillings[(tmp+1):(length(quantiles_fillings))] = quantiles_fillings[(tmp+1):(length(quantiles_fillings))] +1
  # }

  max_fillings_fac = max_fillings
  max_fillings_fac[ (max_fillings < quantiles_fillings[2])] = "Level 1"
  max_fillings_fac[(max_fillings >= quantiles_fillings[2]) ] = "Level 2"
  
  net1%v%"max_fillings_fac" =  c(max_fillings_fac,rep(NA,patient_count_time1))
  net1%v%"mean_fillings" =   c(mean_fillings,rep(NA,patient_count_time1))
  # Where inventors part of international collaborations in the past? 
  internationalisation_dt = df_time0[,.(internationalisation = mean(percentage_internationals)),by = person_id]
  internationalisation_dt$internationalisation_bin = internationalisation_dt$internationalisation== 1
  internationalisation = internationalisation_dt$internationalisation[match(actors_time1,
                                                                            internationalisation_dt$person_id)]
  # where there are no fillings in the past the percentage of international collaboration is 0 -> national = 1
  internationalisation[is.na(internationalisation)] = 1
  internationalisation_bin = internationalisation_dt$internationalisation_bin[match(actors_time1, 
                                                                                    internationalisation_dt$person_id)]
  internationalisation_bin[is.na(internationalisation_bin)] = T
  internationalisation = 1-internationalisation
  internationalisation_bin = 1-internationalisation_bin
  net1%v%"internationalisation" =   c(internationalisation,rep(NA,patient_count_time1))
  net1%v%"internationalisation_bin" =   c(internationalisation_bin,rep(NA,patient_count_time1))
  
  df_singelton_tmp = df_singelton_past[,.(number_of_singeltons = length(area34)),by = person_id]
  singeltons = df_singelton_tmp$number_of_singeltons[match(actors_time1,
                                                           df_singelton_tmp$person_id)]
  # if the actors could not be found in the singelton dt then the respective actors did not yet have any singeltons -> 0
  singeltons[is.na(singeltons)] = 0
  net1%v%"number_singeltons_lag" =   c(singeltons/number_of_complete_patents,rep(NA,patient_count_time1))
  net1%v%"singeltons_bin" =   c(singeltons>0,rep(NA,patient_count_time1))
  # How many patents did the inventors have in the previous year?
  number_patents = rowSums(df_time0_mat)
  net1%v%"number_patent_lag" = c(number_patents/number_of_complete_patents,rep(NA,patient_count_time1))
  mean_patents = number_patents/number_steps
  mean_patents[is.nan(mean_patents)] = 0
  net1%v%"mean_patents" = c(mean_patents,rep(NA,patient_count_time1))
  # Average number patents per step per person
  quantiles_number_patent = quantile(number_patents[number_patents>0]/number_steps[number_patents>0],probs = c(0,0.5))
  # is any number simultaneously two quantiles?
  # tmp = which(quantiles_number_patent[-length(quantiles_number_patent)] == quantiles_number_patent[-1])
  # # if 1 = 3 then we set 2 equal to the next unique value and 3 equal to the value after  
  # tmp0 = quantiles_number_patent[1] == quantiles_number_patent[3]
  # if(tmp0){
  #   unique_vals = sort(unique(number_patents[number_patents>0]/number_steps[number_patents>0]))
  #   quantiles_number_patent[2] = unique_vals[min(which(unique_vals>quantiles_number_patent[2]))]
  #   quantiles_number_patent[3] = unique_vals[min(which(unique_vals>quantiles_number_patent[2]))]
  # }
  # if 1 = 2 then we set 2 equal to the next unique value 
  tmp1 = quantiles_number_patent[1] == quantiles_number_patent[2]
  if(tmp1){
    unique_vals = sort(unique(number_patents[number_patents>0]/number_steps[number_patents>0]))
    quantiles_number_patent[2] = unique_vals[min(which(unique_vals>quantiles_number_patent[2]))]
  }
  # # if then 2 = 3 we also set 3 equal to the next unique value 
  # # if 2 = 3 we also set 3 equal to the next unique value 
  # tmp2 = quantiles_number_patent[3] == quantiles_number_patent[2]
  # if(tmp2){
  #   unique_vals = sort(unique(number_patents[number_patents>0]/number_steps[number_patents>0]))
  #   quantiles_number_patent[3] = unique_vals[min(which(unique_vals>quantiles_number_patent[3]))]
  # }
  
  # if(length(tmp)>0){
  #   quantiles_number_patent[(tmp+1):(length(quantiles_number_patent))] = quantiles_number_patent[(tmp+1):(length(quantiles_number_patent))] +1
  # }
  
  number_patents_fac = number_patents/number_steps
  number_patents_fac[is.nan(number_patents_fac)] = "Level 1"
  number_patents_fac[(number_patents/number_steps < quantiles_number_patent[2])] = "Level 1"
  number_patents_fac[(number_patents/number_steps >= quantiles_number_patent[2])] = "Level 2"
  # number_patents_fac[number_patents/number_steps >= quantiles_number_patent[3]] = "Level 3"

  
  # 
  # number_patents_fac[number_patents %in% c(1:(quantiles_number_patent[2]-1))] = "Level 1"
  # number_patents_fac[number_patents %in% c(quantiles_number_patent[2]:(quantiles_number_patent[3]-1))] = "Level 2"
  # number_patents_fac[number_patents >= quantiles_number_patent[3]] = "Level 3"
  
  
  # number_patents_fac = number_patents
  # number_patents_fac[number_patents == 1] = "Level 1"
  # number_patents_fac[number_patents %in% c(2:3)] = "Level 2"
  # number_patents_fac[number_patents >= 4] = "Level 3"
  net1%v%"number_patents_fac" = c(number_patents_fac,rep(NA,patient_count_time1))
  net1%v%"continuation" = c(number_patents>0,rep(NA,patient_count_time1))
  net1%v%"onset" = c(number_patents==0,rep(NA,patient_count_time1))
  
  # number_patents_1 = number_patents == 1
  # number_patents_2 = number_patents %in% c(2:3)
  # number_patents_3 = number_patents > 4
  # 
  # net1%v%"number_patents_1" =   c(number_patents_1,rep(NA,patient_count_time1))
  # net1%v%"number_patents_2" =   c(number_patents_2,rep(NA,patient_count_time1))
  # net1%v%"number_patents_3" =   c(number_patents_3,rep(NA,patient_count_time1))
  # 
  
  
  # In what percentage of patents of the previous year was the inventor involved?
  # net1%v%"number_patent_scaled_lag" = c(rowSums(df_time0_mat)/number_patents_time0,rep(NA,patient_count_time1))
  # How many single-inventor patents did the inventor have in the previous year?
  # colsums_0 = colSums(df_time0_mat)
  # net1%v%"one_star_patent_lag" = c(rowSums(df_time0_mat[,colsums_0 == 1]),rep(NA,patient_count_time1))
  # # What percentage of single-inventor patents in the previous year are attributed to the inventor?
  # net1%v%"one_star_patent_scaled_lag" = c(rowSums(df_time0_mat[,colsums_0 == 1])/sum(colsums_0 == 1),rep(NA,patient_count_time1))
  # # How many two-inventor patents did the inventor have in the previous year?
  # net1%v%"two_star_patent_lag" = c(rowSums(df_time0_mat[,colsums_0 == 2]),rep(NA,patient_count_time1))
  # # What percentage of two-inventor patents in the previous year are attributed to the inventor?
  # net1%v%"two_star_patent_scaled_lag" = c(rowSums(df_time0_mat[,colsums_0 == 2])/sum(colsums_0 == 2),rep(NA,patient_count_time1))
  # # How many three-inventor patents did the inventor have in the previous year?
  # net1%v%"three_star_patent_lag" = c(rowSums(df_time0_mat[,colsums_0== 3]),rep(NA,patient_count_time1))
  # # What percentage of two-inventor patents in the previous year are attributed to the inventor?
  # net1%v%"three_star_patent_scaled_lag" = c(rowSums(df_time0_mat[,colsums_0== 3])/sum(colsums_0 == 3),rep(NA,patient_count_time1))
  # 
  # 
  # net1%v%"rep_patent_lag" = net1%v%"number_patent_lag" > 0
  
  # actor_data_tmp[person_id == "171657"]
  # tmp = actor_data_tmp[,.(names(table(area34))[which.max(table(area34))]), by = person_id]
  # tmp[person_id == "171657"]
  # tmp = actor_data_tmp[,.(paste(unique(area34),collapse = "_")), by = person_id]
  # 
  tmp_actor_area = actor_data_tmp$area34[match(actors_time1, actor_data_tmp$person_id)]
  
  tmp_actor_area[is.na(tmp_actor_area)] = "No prior patent"
  tmp_actor_area = factor(tmp_actor_area,levels = c("Audiovisual", "Electr/Energy", "Telecom", "ComputerTech", 
                                   "BasicCommProcess", "DigitalComm", "Semiconductors", "IT_Methods", 
                                   "No prior patent"))
  
  net1%v%"area_actor" = c(tmp_actor_area,rep(NA,patient_count_time1))
  # # Include a dyadic covariate that only checks for i-k if inventor i has already had a patent in the same area as patent k is  
  actor_areas = data.table(person_id = actors_time1,"ComputerTech" = 0, "Electr/Energy" = 0,"Telecom" = 0, "Audiovisual" = 0,
                           "DigitalComm" = 0, "Semiconductors" = 0, "BasicCommProcess" = 0, "IT_Methods" = 0)
  actor_data_needed = actor_data_tmp[person_id %in% actor_areas$person_id]
  actor_data_needed$person_id_num = match(actor_data_needed$person_id, actor_areas$person_id)
  actor_data_needed$area34_id = match(actor_data_needed$area34,names(actor_areas))
  # actor_data_needed[is.na(area34_id)]$area34_id = 3
  actor_areas = as.matrix(actor_areas)
  actor_areas[cbind(actor_data_needed$person_id_num,actor_data_needed$area34_id)] = 
    rep(1,times = length(actor_data_needed$person_id_num))
  # What are the patent-specific areas?
  # We here transform the character to a numeric one
  cov1$area34_id = match(cov1$area34,colnames(actor_areas))
  # cov1[is.na(area34_id)]$area34_id = 3
  
  homophily_area_mat = actor_areas[,cov1$area34_id]
  colnames(homophily_area_mat) = cov1$docdb_family_id
  rownames(homophily_area_mat) = actors_time1
  # actor_data_needed
  # cov1

  
  # Step 2 c) Inventor pairwise information - Distance ----
  # Calculate the pairwise distances
  # centroids = data.table(person_id = 1:length(actors_time1),
  #   lon = df_past_present$lon[match(actors_time1, df_past_present$person_id)],
  #   lat = df_past_present$lat[match(actors_time1, df_past_present$person_id)])
  tmp_lat = df_past_present$lat[match(actors_time1, df_past_present$person_id)]
  tmp_lon = df_past_present$lon[match(actors_time1, df_past_present$person_id)]
  tmp = data.table(matrix(combinations(x = length(actors_time1)),ncol = 2,byrow = T))
  names(tmp) = c("from", "to")
  tmp = tmp[from <to]
  # tmp$from_to = paste(tmp$from, tmp$to)
  
  tmp$lon_from = tmp_lon[tmp$from]
  tmp$lat_from = tmp_lat[tmp$from]
  tmp$lon_to = tmp_lon[tmp$to]
  tmp$lat_to = tmp_lat[tmp$to]
  tmp$dist_km2 = distHaversine(cbind(tmp$lon_from, tmp$lat_from),   
                               cbind(tmp$lon_to, tmp$lat_to), r=6378137)/1000
  
  distance_matrix = matrix(data = 0,nrow = nrow(df_time0_mat),ncol = nrow(df_time0_mat))
  rownames(distance_matrix) = colnames(distance_matrix) =actors_time1
  distance_matrix[cbind(tmp$from,tmp$to)] = tmp$dist_km2
  distance_matrix[cbind(tmp$to,tmp$from)] = tmp$dist_km2
  distance_matrix = floor(distance_matrix/10)*10
  
  distance_matrix_bin = distance_matrix
  distance_matrix_bin[distance_matrix>=50] = 0
  distance_matrix_bin[distance_matrix<50] = 1
  distance_matrix_log = log1p(distance_matrix)
  gc(full = T)
  
  # Step 2 c) Inventor pairwise information - Common Partner ----
  # calculate the number of common partners and some other endog. statistics from the past 
  # Generate a dt of all active actots in step t
  # df_past$person_id = factor(df_past$person_id, levels = unique(c(df_past$person_id,actors_time1)))
  patent_data_time0 = df_past[,.(inventor = list(person_id)), by = docdb_family_id]
  inventor_data_time0 = df_past[,.(patents = list(docdb_family_id)), by = person_id]
  # delete all actors that are not active in t and t -1 but in the last n_steps_inclusion_past steps
  inventor_data_time0 = inventor_data_time0[person_id %in% actors_time1]
  
  # add inventors that were not active in the last three years 
  # inventor_data_time0 = rbind(inventor_data_time0, 
  #                             data.table(person_id = actors_time1[!actors_time1%in% inventor_data_time0$person_id],
  #                                        patents = list()))
  # inventor_data_time0 = inventor_data_time0[order(person_id)]
  # patent_data_time0 = patent_data_time0[order(docdb_family_id)]
  tmp_fun = function(patents, inventor, patent_data_time0){
    # whether or not the unique stays depends on how we count repeated joint patents 
    res = unique(unlist(patent_data_time0[match(patents[[1]], patent_data_time0$docdb_family_id)]$inventor))
    # res = (unlist(patent_data_time0[match(patents[[1]], patent_data_time0$docdb_family_id)]$inventor))
    # return(res[res != inventor])
    return(paste(res[res != inventor],sep = "_"))
  }
  # in tmp_save we save each person <-> person that we observed (this interaction is from a joint patent in the past)
  tmp_save = inventor_data_time0[,.(involved_actots = (tmp_fun(patents = patents,
                                                               inventor = person_id, 
                                                               patent_data_time0 = patent_data_time0))), 
                                 by = person_id]
  
  number_part_per_id = data.table(number_partners = table(tmp_save$person_id))
  names(number_part_per_id) = c("person_id", "number_partners")
  number_partners = number_part_per_id$number_partners[match(actors_time1, number_part_per_id$person_id)]
  number_partners[is.na(number_partners)] = 0
  
  quantiles_number_partners = quantile(number_partners[number_partners>0],type =8,probs = c(0,0.33,0.66))
  # is any number simultaneously two quantiles?
  tmp = which(quantiles_number_partners[-length(quantiles_number_partners)] == quantiles_number_partners[-1])
  if(length(tmp)>0){
    quantiles_number_partners[(tmp+1):(length(quantiles_number_partners))] = quantiles_number_partners[(tmp+1):(length(quantiles_number_partners))] +1
  }
  
  number_partners_fac = number_partners
  number_partners_fac[(number_partners >= quantiles_number_partners[1]) & (number_partners < quantiles_number_partners[2])] = "Level 1"
  number_partners_fac[(number_partners >= quantiles_number_partners[2]) & (number_partners< quantiles_number_partners[3])] = "Level 2"
  number_partners_fac[number_partners >= quantiles_number_partners[3]] = "Level 3"
  
  
  # 
  # 
  # number_partners_fac = number_partners
  # number_partners_fac[number_partners == 1] = "Level 1"
  # number_partners_fac[number_partners %in% c(2:3)] = "Level 2"
  # number_partners_fac[number_partners >= 4] = "Level 3"
  # 
  net1%v%"number_partners_fac" = c(number_partners_fac,rep(NA,patient_count_time1))
  
  # number_partners_2 = number_partners %in% c(2:3)
  # number_partners_3 = number_partners > 4
  # net1%v%"number_partners_lag" =   c(number_partners,rep(NA,patient_count_time1))
  # net1%v%"number_partners_1" =   c(number_partners_1,rep(NA,patient_count_time1))
  # net1%v%"number_partners_2" =   c(number_partners_2,rep(NA,patient_count_time1))
  # net1%v%"number_partners_3" =   c(number_partners_3,rep(NA,patient_count_time1))

  # data.table(t(comb_n(tmp_save[involved_actots == 12605401]$person_id,k = 2)))
  fun_tmp = function(x){
    if(length(x)>1){
      return(list(data.table(t(comb_n(x,k = 2)))))
    } else {
      list()
    }
  }
  # tmp_save = data.table(person_id = c(1,2),involved_actots = c(3,3))
  # Now for each actor on the right side (the invovled actor),
  # we look at the set of actors she/he/it is involved with and 
  # return a dt of all possible pairwise combinations of them
  # e.g. if 1 and 2 have had a patent with 3, 1 and 2 have at 
  # least one common partner (being 3)
  tmp_save_common_partner = tmp_save[,.(tmp = fun_tmp(person_id)), by = involved_actots]
  tmp_save_common_partner = rbindlist(tmp_save_common_partner$tmp)
  tmp_save_common_partner$ind = paste(tmp_save_common_partner$V1,tmp_save_common_partner$V2,sep = "_")
  tmp_save_common_partner = tmp_save_common_partner[,.(from = V1[1], to = V2[1], times = length(V1)), by = .(ind)]
  
  # tmp$from_id = actors_time1[tmp$from]
  # tmp$to_id = actors_time1[tmp$to]
  # 
  tmp_save_common_partner$from_number = match(tmp_save_common_partner$from,actors_time1)
  tmp_save_common_partner$to_number = match(tmp_save_common_partner$to,actors_time1)
  common_partner_matrix = matrix(data = 0,nrow = nrow(df_time0_mat),ncol = nrow(df_time0_mat))
  rownames(common_partner_matrix) = colnames(common_partner_matrix) =actors_time1
  # common_partner_matrix[cbind(tmp_save_common_partner$from_number, tmp_save_common_partner$to_number)] = tmp_save_common_partner$times
  # common_partner_matrix[cbind(tmp_save_common_partner$to_number, tmp_save_common_partner$from_number)] = tmp_save_common_partner$times
  common_partner_matrix[cbind(tmp_save_common_partner$from_number, tmp_save_common_partner$to_number)] = 1
  common_partner_matrix[cbind(tmp_save_common_partner$to_number, tmp_save_common_partner$from_number)] = 1
  
  # Step 2 c) Inventor pairwise information - Repetition ----
  # Now derive which actors already interacted 
  repetition_matrix = matrix(data = 0,nrow = nrow(df_time0_mat),ncol = nrow(df_time0_mat))
  rownames(repetition_matrix) = colnames(repetition_matrix) =actors_time1
  tmp_save$from_number = match(tmp_save$person_id,actors_time1)
  tmp_save$to_number = match(tmp_save$involved_actots,actors_time1)
  tmp_save$ind = paste(tmp_save$from_number,tmp_save$to_number,sep = "_")
  
  # there are some involved inventors that were only active in the past but are now not anymore in the network 
  # -> we delete them as they are not anymore in the network  
  tmp_save = tmp_save[!is.na(to_number)]
  # table(tmp_save[,.(from_number = from_number[1], to_number = to_number[1], times = length(to_number)),by = ind]$times)
  
  repetition_matrix[cbind(tmp_save$from_number, tmp_save$to_number)] = 1
  repetition_matrix[cbind(tmp_save$to_number, tmp_save$from_number)] = 1
 # browser()
  # two_star_matrix = matrix(data = 0,nrow = nrow(df_time0_mat),ncol = nrow(df_time0_mat))
  # rownames(two_star_matrix) = colnames(two_star_matrix) =actors_time1
  # two_star_matrix[,] = 1
  # diag(two_star_matrix) = 0
  # 
  # Add a list variable indicating the patent number ids each actor was involved in 
  # inventor_data_time0$involved_actors = tmp_save$involved_actots
  # inventor_data_time0
  # 
  # trying = outerList(inventor_data_time0$involved_actors, inventor_data_time0$involved_actors, FUN = function(x,y){length(intersect(x,y))})
  
  # unlist(c(inventor_data_time0$involved_actors[c(1,2)],inventor_data_time0$involved_actors[c(3,4)]))
  
  # tmp$from_involved_actors =  inventor_data_time0$involved_actors[tmp$from]
  # tmp$to_involved_actors =  inventor_data_time0$involved_actors[tmp$to]
  # 
  # parallelStart(mode = "socket",cpus = 20)
  # tmp$common = parallelMap(intersect, tmp$from_involved_actors,tmp$to_involved_actors)
  # parallelStop()

  # Derive common partners of actors 
  
  # get the change statistics
  # trying = ergm(net1 ~ two_start_attr(distance_matrix_bin, "dist_coef", nrow(distance_matrix)) +
  #            b2factor("area_patent") +
  #            b1cov("number_patent_lag") + 
  #            # b1cov("number_patent_scaled_lag") + 
  #            # b1cov("one_star_patent_lag")  + 
  #            b1cov("two_star_patent_lag") +
  #            # b1cov("two_star_patent_scaled_lag") + 
  #            b1cov("three_star_patent_lag") + 
  #            # b1cov("three_star_patent_scaled_lag") + 
  #            # b1nodematch("area_actor", beta = 0,byb2attr = "area_patent") +
  #            b1factor("area_actor") +
  #            b1cov("lon") + 
  #            b1cov("lat") + 
  #            # b1cov("number_steps") + 
  #            b1cov("rep_patent_lag"),estimate = "MPLE")
  # 
  # g.sim = simulate(net1 ~  two_star_avr(distance_matrix_bin, "dist_mean", nrow(distance_matrix_bin)) +
  #                     two_star_avr(repetition_matrix, "repetition_mean", nrow(repetition_matrix)) +
  #                     two_star_avr(common_partner_matrix, "common_mean", nrow(common_partner_matrix)) +
  #                     b1star(2)+ b1cov("number_steps"), coef=ergm_tmp$coefficients,basis = net1)
  # 
  # ergm_tmp = ergm(net1 ~ two_star_avr(distance_matrix_bin, "dist_mean", nrow(distance_matrix_bin)) +
  #                    two_star_avr(repetition_matrix, "repetition_mean", nrow(repetition_matrix)) +
  #                    two_star_avr(common_partner_matrix, "common_mean", nrow(common_partner_matrix)) +
  #                    b1star(2)+ b1cov("number_steps"),estimate = "MPLE")
 
 
  
 # Step 3 Calculate and return change stats
  # data_full = ergmMPLE(net1 ~ b1factor("onset") + 
  #                        b1factor("continuation") + 
  #                        b1star(2)+
  #                        b1factor("gender") +
  #                        b1cov("internationalisation_bin") + 
  #                        b1factor("max_fillings_fac") +
  #                        b1factor("number_patents_fac") + 
  #                        b1factor("singeltons_bin")  +
  #                        b2star(2)+
  #                        two_star_avr(repetition_matrix, "repetition_mean", nrow(repetition_matrix)) +
  #                        two_star_avr(common_partner_matrix, "common_mean", nrow(common_partner_matrix)) +
  #                        two_star_avr(distance_matrix_bin, "dist_mean", nrow(distance_matrix_bin)) +
  #                        dyadcov(homophily_area_mat),
  #                      # b1cov("mean_fillings") + 
  #                      # b1cov("number_steps") +
  #                      # b1cov("rep_patent_lag") +, 
  #                      output = "matrix",
  #                      control = control.ergm(MPLE.max.dyad.types = network.dyadcount(net1) *4))
  # # data_full = ergmMPLE(net1 ~ 
  # #                        # two_star_coef(distance_matrix_bin, "dist_sum", nrow(distance_matrix_bin)) +
  # #                        # two_star_coef(repetition_matrix, "repetition_sum", nrow(repetition_matrix)) +
  # #                        # two_star_coef(common_partner_matrix, "common_sum", nrow(common_partner_matrix)) +
  # #                        two_star_avr(distance_matrix_bin, "dist_mean", nrow(distance_matrix_bin)) +
  # #                        two_star_avr(repetition_matrix, "repetition_mean", nrow(repetition_matrix)) +
  # #                        two_star_avr(common_partner_matrix, "common_mean", nrow(common_partner_matrix)) +
  # #                        # two_star_avr(two_star_matrix, "twostar_mean", nrow(common_partner_matrix)) +
  # #                        # two_star_indep(two_star_matrix, "twostar_indep", nrow(common_partner_matrix),patient_count_time1) +
  # #                        # two_star_indep(repetition_matrix, "repetition_indep", nrow(repetition_matrix),patient_count_time1) +
  # #                        # two_star_indep(common_partner_matrix, "common_indep", nrow(common_partner_matrix),patient_count_time1) +
  # #                        b2star(2)+
  # #                        # b2cov("area_patent") +
  # #                        # b1cov("area_actor") +
  # #                        b1factor("gender") +
  # #                        b1star(2)+
  # #                        dyadcov(homophily_area_mat)+
  # #                        b1factor("singeltons_bin") + 
  # #                        b1factor("number_patents_fac") +
  # #                        # b1cov("number_patent_scaled_lag") + 
  # #                        # b1cov("one_star_patent_lag")  + 
  # #                        # b1cov("two_star_patent_lag") +
  # #                        # b1cov("two_star_patent_scaled_lag") +
  # #                        # b1cov("three_star_patent_lag") + 
  # #                        # b1cov("three_star_patent_scaled_lag") +
  # #                        # b1factor("area_actor") +
  # #                        # b1cov("lon") +
  # #                        # b1cov("lat") +
  # #                        # b1cov("internationalisation") + 
  # #                        b1cov("internationalisation_bin") + 
  # #                        b1factor("max_fillings_fac") + 
  # #                        b1factor("onset") + 
  # #                        b1factor("continuation"),
  # #                        # b1cov("mean_fillings") + 
  # #                        # b1cov("number_steps") +
  # #                        # b1cov("rep_patent_lag") +, 
  # #                       output = "matrix",
  # #                       control = control.ergm(MPLE.max.dyad.types = network.dyadcount(net1) *4))
  # 
  # 
  # gc(full = T)
  # data_full = data.table(cbind(data_full$response, data_full$predictor, weight = data_full$weights))
  # names(data_full)[1] = "Y"
  # # Offset 
  # data_full$offset = log(1/sum(dim(df_time0_mat)))
  # names(data_full) = gsub(pattern = " ",replacement = ".",x = names(data_full))
  # Step 3: Estimation -----
  gc(reset = T,full = T)
  # table(net1%v% "singeltons_bin")
  
  # model_alt= ergm(net1 ~ b1factor("onset") + 
  #                   b1factor("continuation") + 
  #                   b1star(2)+
  #                   b1factor("gender") +
  #                   b1cov("internationalisation_bin") + 
  #                   b1factor("max_fillings_fac") +
  #                   b1factor("number_patents_fac") + 
  #                   b1factor("singeltons_bin")  +
  #                   b2star(2)+
  #                   two_star_avr(repetition_matrix, "repetition_mean", nrow(repetition_matrix)) +
  #                   two_star_avr(common_partner_matrix, "common_mean", nrow(common_partner_matrix)) +
  #                   two_star_avr(distance_matrix_bin, "dist_mean", nrow(distance_matrix_bin)) +
  #                   dyadcov(homophily_area_mat) + 
  #                   offset(dyadcov(matrix(log(1/(actor_count_time1 + patient_count_time1)),
  #                                         nrow = actor_count_time1,ncol = patient_count_time1))),
  #                 offset.coef = c(1), control = control.ergm(main.method = "Stochastic-Approximation"))
  # model_sa = model_alt
  # mcmc.diagnostics(model_tmp)
  # 
  # model_alt= ergm(net1 ~ b1factor("onset") + 
  #                   b1factor("continuation") + 
  #                   b1star(2)+
  #                   b1factor("gender") +
  #                   b1cov("internationalisation_bin") + 
  #                   b1factor("max_fillings_fac") +
  #                   b1factor("number_patents_fac") + 
  #                   b1factor("singeltons_bin")  +
  #                   b2star(2)+
  #                   two_star_avr(repetition_matrix, "repetition_mean", nrow(repetition_matrix)) +
  #                   two_star_avr(common_partner_matrix, "common_mean", nrow(common_partner_matrix)) +
  #                   two_star_avr(distance_matrix_bin, "dist_mean", nrow(distance_matrix_bin)) +
  #                   dyadcov(homophily_area_mat) + 
  #                   isolates +
  #                   offset(b2mindegree(2)) +
  #                   offset(dyadcov(matrix(log(1/(actor_count_time1 + patient_count_time1)),
  #                                         nrow = actor_count_time1,ncol = patient_count_time1))),estimate = "MPLE",
  #                 offset.coef = c(Inf,1))
  # 
  # model_tmp = ergm(net1 ~ b1factor("onset") + 
  #                    b1factor("continuation") + 
  #                    b1star(2)+
  #                    b1factor("gender") +
  #                    b1cov("internationalisation_bin") +
  #                    b1factor("max_fillings_fac") +
  #                    b1factor("number_patents_fac") + 
  #                    b1factor("singeltons_bin")  +
  #                    b2star(2)+
  #                    two_star_avr(repetition_matrix, "repetition_mean", nrow(repetition_matrix)) +
  #                    two_star_avr(common_partner_matrix, "common_mean", nrow(common_partner_matrix)) +
  #                    two_star_avr(distance_matrix_bin, "dist_mean", nrow(distance_matrix_bin)) +
  #                    dyadcov(homophily_area_mat) + 
  #                    isolates +
  #                    offset(b2mindegree(2)) +
  #                    offset(dyadcov(matrix(log(1/(actor_count_time1 + patient_count_time1)),
  #                                          nrow = actor_count_time1,ncol = patient_count_time1))),estimate = "CD",
  #                  offset.coef = c(Inf,1),control = control.ergm(CD.maxit = 100,CD.nsteps = 1000,
  #                                                            init =coefs_beg))
  # 
  # 

  # browser()
  attr(model_formula, ".Environment") = environment()  

  model_alt= ergm(model_formula,estimate = "MPLE",
                  offset.coef = offset_coefs)
  coefs_beg = coef(model_alt)
  rm(model_alt)
  gc(reset = T,full = T)
  model_tmp = ergm(model_formula,estimate = "CD",
                   offset.coef = offset_coefs,control = control.ergm(CD.maxit = 100,CD.nsteps = number_changes,
                                                             init =coefs_beg,seed = seed))
  # model_tmp2 = ergm(model_formula,estimate = "CD",
  #                  offset.coef = offset_coefs,control = control.ergm(CD.maxit = 100,CD.nsteps = 1000,
  #                                                                    init =coefs_beg))
  # model_tmp3 = ergm(model_formula,estimate = "CD",
  #                   offset.coef = offset_coefs,control = control.ergm(CD.maxit = 100,CD.nsteps = 15000,
  #                                                                     init =coefs_beg))
  # as.vector(coef(model_tmp))
  # as.vector(coef(model_tmp2))
  # as.vector(coef(model_tmp3))
  # model_stepping = ergm(net1 ~ b1factor("onset") + 
  #                    b1factor("continuation") + 
  #                    b1star(2)+
  #                    b1factor("gender") +
  #                    b1cov("internationalisation_bin") +
  #                    b1factor("max_fillings_fac") +
  #                    b1factor("number_patents_fac") + 
  #                    b1factor("singeltons_bin")  +
  #                    b2star(2)+
  #                    two_star_avr(repetition_matrix, "repetition_mean", nrow(repetition_matrix)) +
  #                    two_star_avr(common_partner_matrix, "common_mean", nrow(common_partner_matrix)) +
  #                    two_star_avr(distance_matrix_bin, "dist_mean", nrow(distance_matrix_bin)) +
  #                    dyadcov(homophily_area_mat) + 
  #                    offset(dyadcov(matrix(log(1/(actor_count_time1 + patient_count_time1)),
  #                                          nrow = actor_count_time1,ncol = patient_count_time1))),
  #                  offset.coef = c(1),control = control.ergm(main.method = "Stepping",
  #                                                            init =model_tmp$coefficients))
  
  # trying = gof.formula(model_tmp$formula, coef = model_tmp$coefficients)
  
  # as.vector(model_stepping$coefficients)
  # as.vector(model_tmp$coefficients)
  # x = repetition_matrix
  # x[,] = 1
  # diag(x) = 0
  # 
  # model_alt = ergm(net1 ~ b1factor("onset") + 
  #                    b1factor("continuation") + 
  #                    b1star(2)+
  #                    b1factor("gender") +
  #                    b1cov("internationalisation_bin") +
  #                    b1factor("max_fillings_fac") +
  #                    b1factor("number_patents_fac") +
  #                    b1factor("singeltons_bin")  +
  #                    b2star(2)+
  #                    dyadcov(homophily_area_mat)+
  #                    offset(dyadcov(matrix(log(1/(actor_count_time1 + patient_count_time1)),
  #                                          nrow = actor_count_time1,ncol = patient_count_time1))),
  #                  estimate = "MPLE",
  #                  offset.coef = c(1))
  # 
  # 
  # summary(model_alt)
  # summary(model_tmp)
  # table(model_tmp$glm$offset)
  
  # model_alt = ergm(net1 ~ b1factor("onset") + 
  #                    b1factor("continuation") + 
  #                    b1star(2)+
  #                    b1factor("gender") +
  #                    b1cov("internationalisation_bin") + 
  #                    b1factor("max_fillings_fac") +
  #                    b1factor("number_patents_fac") + 
  #                    b1factor("singeltons_bin")  +
  #                    b2star(2)+
  #                    two_star_avr(repetition_matrix, "repetition_mean", nrow(repetition_matrix)) +
  #                    two_star_avr(common_partner_matrix, "common_mean", nrow(common_partner_matrix)) +
  #                    # two_star_avr(distance_matrix_bin, "dist_mean", nrow(distance_matrix_bin)) +
  #                    dyadcov(homophily_area_mat) + 
  #                    offset(dyadcov(matrix(log(1/(actor_count_time1 + patient_count_time1)),
  #                                          nrow = actor_count_time1,ncol = patient_count_time1))),estimate = "MPLE",
  #                  control = control.ergm(MPLE.max.dyad.types = network.dyadcount(net1) *4),
  #                  offset.coef = c(1))
  
  
  
  
  cat("Starting Simulations\n")
  
  # browser()
  rm(df_singelton, df_sub2,df_time0,df_time1,df_time1_mat,df_time0_mat,df_past,distance_matrix,distance_matrix_log,filling_dt, 
     gender_data, gender_dt, internationalisation,internationalisation_bin,inventor_data_time0,number_part_per_id,
     tmp_save,tmp_save_common_partner, actors_time1, number_partners,number_partners_fac)
  gc(reset = T,full = T)
  # 
  # simulations2 = simulate(model_tmp2,nsim = 100,sequential =F,seed = 123, output = "stats")
  # simulations3 = simulate(model_tmp3,nsim = 100,sequential =F,seed = 123, output = "stats")
  simulations = simulate(model_tmp,nsim = 200,sequential =F,seed = 123, output = "stats")
  # colmeans(simulations2)
  # colmeans(simulations3)
  # colmeans(simulations)
  # as.vector(model_tmp$nw.stats)
  # colmeans(simulations)
  # as.vector(model_tmp$nw.stats)
  fisher_inv = solve(var(simulations)[1:(length(coefs_beg)- length(offset_coefs)),1:(length(coefs_beg)- length(offset_coefs))])
  cat(dim(fisher_inv), "Dimension of the Fisher info \n")
  cat(length(coefs_beg), "Length of the coefficient vector \n")
                                      
  # 
  # simulations_alt = simulate(model_tmp,nsim = 100,sequential =F,seed = 123,output = "stats")
  # 
  # fisher_inv = solve(var(simulations_alt)[1:13,1:13])
  # diag(fisher_inv)
  # diag(model_alt$covar)
  # as.vector(apply(simulations_alt,MARGIN = 2,var))
  # 
  # colmeans(simulations_alt)
  # 
  # as.vector(model_tmp$nw.stats)
  # 
  # colmeans(simulations)
  # as.vector(model_tmp$nw.stats)
  # as.vector(tmp)
  # colnames(simulations2)
  # colnames(simulations)
  # 
  # colmeans(simulations2)
  # as.vector(model_alt$nw.stats)
  # 
  # plot(simulations[,9],simulations2[,9])
  # mean(simulations[,12])
  # 
  # net1$gal$mnext
  # 
  # mean(unlist(lapply(simulations, function(x){x$gal$mnext})))
  # 
  # model_tmp = glm(Y ~ -1 + b1factor.onset.TRUE + 
  #                   b1factor.continuation.TRUE + 
  #                   b1star2+
  #                   b1factor.gender.M + 
  #                   b1cov.internationalisation_bin + 
  #                   b1factor.max_fillings_fac.Level.2 +
  #                   b1factor.number_patents_fac.Level.2+
  #                   b1factor.singeltons_bin.TRUE+
  #                   b2star2+
  #                   repetition_mean+ 
  #                   common_mean+ 
  #                   dist_mean + 
  #                   dyadcov.homophily_area_mat, data = data_full, 
  #                 weights = weight, offset = offset, family = binomial)
  
  # 
  #     
  #     simulations = simulate(net1 ~ 
  #                              b1factor("onset") + 
  #                              b1factor("continuation") + 
  #                              b1star(2)+
  #                              b1factor("gender") +
  #                              b1cov("internationalisation_bin") + 
  #                              b1factor("max_fillings_fac") +
  #                              b1factor("number_patents_fac") + 
  #                              b1factor("singeltons_bin")  +
  #                              b2star(2)+
  #                              two_star_avr(repetition_matrix, "repetition_mean", nrow(repetition_matrix)) +
  #                              two_star_avr(common_partner_matrix, "common_mean", nrow(common_partner_matrix)) +
  #                              two_star_avr(distance_matrix_bin, "dist_mean", nrow(distance_matrix_bin)) +
  #                              dyadcov(homophily_area_mat), coef = model_tmp$coefficients,
  #                            nsim = 1,
  #                            sequential = T
  #     )
  gc(reset = T,full = T)
  
  if(bs){
    cat("Starting estimating the model on the BS\n")
    trying = lapply(X = simulations,FUN = estimate_list_function,
                    repetition_matrix = repetition_matrix,
                    common_partner_matrix = common_partner_matrix,
                    distance_matrix_bin = distance_matrix_bin, 
                    homophily_area_mat = homophily_area_mat, 
                    patient_count_time1 = patient_count_time1, 
                    actor_count_time1 = actor_count_time1)
    bs_results = rbindlist(lapply(trying, function(x){data.table(t(x))}))
    
  }
  
  
 
  cat("Step ",i," finished\n")
  cat("Needed time was:",Sys.time() -now," Seconds \n")
  # cat("Number of covariates:",ncol(data_full),"\n")
  if(!is.null(dir_path)) {
    if(bs){
      save_path = paste0(dir_path,"/bs_results/bs_res_",i,".rds")
      saveRDS(bs_results,file = save_path)
    }
    
    # Save estimates
    if(dir.exists(paste0(dir_path,"/estimates/"))){
      save_path = paste0(dir_path,"/estimates/coef_",i,".rds")
      saveRDS(model_tmp$coefficients,file = save_path)
    } else {
      dir.create(paste0(dir_path,"/estimates/"))
      save_path = paste0(dir_path,"/estimates/coef_",i,".rds")
      saveRDS(model_tmp$coefficients,file = save_path)
    }
    # Save fisher info
    if(dir.exists(paste0(dir_path,"/sampled_networks/"))){
      save_path = paste0(dir_path,"/sampled_networks/bs_res_",i,".rds")
      saveRDS(fisher_inv,file = save_path)
    } else {
      dir.create(paste0(dir_path,"/sampled_networks/"))
      save_path = paste0(dir_path,"/sampled_networks/bs_res_",i,".rds")
      saveRDS(fisher_inv,file = save_path)
    }
    # Save model 
    if(dir.exists(paste0(dir_path,"/covariates/"))){
      save_path = paste0(dir_path,"/covariates/Model_",i,".rds")
      saveRDS(model_tmp,file = save_path)
    } else {
      dir.create(paste0(dir_path,"/covariates/"))
      save_path = paste0(dir_path,"/covariates/Model_",i,".rds")
      saveRDS(model_tmp,file = save_path)
    }
    # Save network 
    if(dir.exists(paste0(dir_path,"/networks/"))){
      save_path = paste0(dir_path,"/networks/Network_",i,".rds")
      saveRDS(net1,file = save_path)
    } else {
      dir.create(paste0(dir_path,"/networks/"))
      save_path = paste0(dir_path,"/networks/Network_",i,".rds")
      saveRDS(net1,file = save_path)
    }
   
  
    } else {
    res = list(model = model_tmp, network = net1)
    return(res)
  }
  
}



prepare_network = function(i,df_sub2, n_steps_inclusion = 5,mainarea = "Electrical Engineering", include_actor_data= F) {
  now = Sys.time()
  df_sub2 = df_sub2[mainarea34 == mainarea]
  
  # define the bimatrix at time point t 
  df_time1 = df_sub2[steps == i]
  # df_time1$area34 = factor(df_time1$area34,levels = unique(df_sub2$area34))
  df_past_present = df_sub2[steps %in% ((i-(n_steps_inclusion -1)):i)]
  # df_past_present$area34 = factor(df_past_present$area34,levels = unique(df_sub2$area34))
  df_past = df_past_present[steps != i]
  # only look at the single patents
  cov1 = df_time1[!duplicated(df_time1[, "docdb_family_id"]), ]
  actors_time1 = unique(df_past_present$person_id)
  df_time1_mat = from_el_to_bimatrix_alt(el = df_time1,mode1 =  "person_id",
                                         mode2 = "docdb_family_id",actor_set = actors_time1)
  
  patient_count_time1 = length(colnames(df_time1_mat))
  patients_time1 = colnames(df_time1_mat)
  
  # df.time1 = df.time1 %>% arrange(docdb_family_id)
  # cov1 = df.time1[!duplicated(df.time1[, "docdb_family_id"]), ]
  # 
  net1 = network(df_time1_mat, matrix.type = "bipartite", ignore.eval = FALSE)
  # data_tmp = data.table(
  #   from = as.character(df_time1$person_id),
  #   to =as.character(df_time1$docdb_family_id),
  #   stringsAsFactors = FALSE
  # )
  # vertices = data.table(c(patients_time1,actors_time1), 
  #            is_actor = c(rep(FALSE,times = length(patients_time1)), 
  #                         rep(TRUE,times = length(actors_time1))))
  # debugonce(as.network)
  # trying = as.network(data_tmp, bipartite = T, directed = F,vertices = vertices)
  
  # define the bimatrix at time point t  -1 
  df_time0 = df_past[person_id %in% actors_time1]
  df_time0_mat = from_el_to_bimatrix_alt(el = df_time0, mode1 = "person_id",
                                          mode2 =  "docdb_family_id", actor_set = actors_time1)
  net1%v%"area_patent" = c(rep(as.character(NA), nrow(df_time1_mat)), cov1$area34)
  net1%v%"lon" =   c(df_past_present$lon[match(actors_time1, df_past_present$person_id)],rep(NA,patient_count_time1))
  net1%v%"lat" =   c(df_past_present$lat[match(actors_time1, df_past_present$person_id)],rep(NA,patient_count_time1))
  net1%v%"number_patent_lag" =   c(rowSums(df_time0_mat),rep(NA,patient_count_time1))
  net1%v%"one_star_patent_lag" =    c(rowSums(df_time0_mat[,colSums(df_time0_mat) == 1]),rep(NA,patient_count_time1))
  net1%v%"two_star_patent_lag" =    c(rowSums(df_time0_mat[,colSums(df_time0_mat) == 2]),rep(NA,patient_count_time1))
  net1%v%"three_star_patent_lag" =    c(rowSums(df_time0_mat[,colSums(df_time0_mat) == 3]),rep(NA,patient_count_time1))
  net1%v%"rep_patent_lag" = net1%v%"number_patent_lag" > 0
  net1%v%"area_actor" = c(df_past_present$area34[match(actors_time1, df_past_present$person_id)],rep(NA,patient_count_time1))
  
  number_steps_per_actors = df_time0[,.(number_steps = length(unique(steps))), by = person_id]
  tmp_number = c(rep(0,nrow(df_time1_mat)),rep(NA,patient_count_time1))
  tmp_match = match(number_steps_per_actors$person_id,net1%v%"vertex.names")
  tmp_number[tmp_match] = number_steps_per_actors$number_steps
  net1%v%"number_steps" = tmp_number
  net1%v%"is_inventor" = is.na(net1%v%"area_patent")
  if(include_actor_data){
    actor_data = data.table(names =  net1%v%"vertex.names", 
                            area = net1%v%"area_patent",
                            lon =  net1%v%"lon", 
                            lat =  net1%v%"lat",   
                            number_patent_lag =  net1%v%"number_patent_lag", 
                            one_star_patent_lag =  net1%v%"one_star_patent_lag", 
                            two_star_patent_lag =  net1%v%"two_star_patent_lag", 
                            three_star_patent_lag =  net1%v%"three_star_patent_lag", 
                            number_steps =  tmp_number, 
                            is_inventor = is.na(net1%v%"area_patent")
    )
    return(list(network =net1, actor_data = actor_data))
  } else {
    return(net1)
  }
  
  
}

prepare_network_year = function(df_sub2, year_from, year_to,mainarea = "Electrical Engineering",subarea = NULL, include_actor_data= F) {
  now = Sys.time()
  if(is.null(subarea)){
    df_sub2 = df_sub2[mainarea34 == mainarea & year %in% year_from:year_to]  
  } else {
    df_sub2 = df_sub2[area34 == subarea & year %in% year_from:year_to]
  }
  
  
  # define the bimatrix at time point t 
  df_time1 = df_sub2
  df_time1_mat = from_el_to_bimatrix(el = df_time1,mode1 =  "person_id",
                                         mode2 = "docdb_family_id")
  patient_count_time1 = length(colnames(df_time1_mat))
  patients_time1 = colnames(df_time1_mat)
  
  # df.time1 = df.time1 %>% arrange(docdb_family_id)
  # cov1 = df.time1[!duplicated(df.time1[, "docdb_family_id"]), ]
  # 
  net1 = network(df_time1_mat, matrix.type = "bipartite", ignore.eval = FALSE)
  net1%v%"is_inventor" = net1%v%"vertex.names" %in% rownames(df_time1_mat)
  return(net1)
}


full_estimation = function(df,df_singelton, n_steps_inclusion_actor, n_steps_inclusion_past,
                      time_steps_eff,stepsize_months, seed,
                      mainarea = "Electrical Engineering", cl = NULL, continues = F, bs= T, 
                      offset_coefs ,
                      model_formula,number_changes = 1000 ) {
  # df = df[mainarea34 == mainarea]
  actor_data = df
  actor_data = actor_data[!duplicated(cbind(actor_data$person_id,actor_data$area34,actor_data$steps))]
  actor_data = actor_data[,.(person_id, area34, lat, lon, steps,gender)]
  
  i = time_steps_eff[2]
  gc(full = T)
  
  dir_path = paste0("../Data/preprocessing_", gsub(" ",replacement = "_",x = mainarea), "_",n_steps_inclusion_actor, "_",n_steps_inclusion_past,"_",stepsize_months)
  if(continues){
    if(!dir.exists(dir_path)){
      dir.create(dir_path)
    } else{
      # If we continue and some steps are already carried out we need to check which steps are already completed 
      data_files = list.files(dir_path,full.names = T,all.files = T,recursive = T)
      data_files = gsub(".rds",replacement = "", data_files)
      completed_steps = as.numeric(sub(".*_", "", data_files))
      # Those completed steps are then excluded from the preprocessing loop
      time_steps_eff = time_steps_eff[!time_steps_eff%in% completed_steps]
      cat("Since the estimation is to be continued the excecuted steps are: ",time_steps_eff, "\n")
      }
  } else {
    if(!dir.exists(dir_path)){
      dir.create(dir_path)
    } else{
      # if the directory exists delete all data that was there before
      data_files = list.files(dir_path,full.names = T,all.files = T,recursive = T)
      file.remove(data_files)
      dir.create(dir_path)
    }
  }
  
  
  if(is.null(cl)){
    for(i in time_steps_eff){
      prepare_data_no_re(i = i,df_sub2 = df,df_singelton = df_singelton, actor_data = actor_data, 
                         n_steps_inclusion_actor = n_steps_inclusion_actor,
                         n_steps_inclusion_past = n_steps_inclusion_past,
                         dir_path = dir_path, bs = bs, seed = seed, 
                         offset_coefs = offset_coefs,
                         model_formula = model_formula, number_changes = number_changes)
      
      # if(re){
      #   result_tmp = prepare_data(i = i,df_sub2 = df,actor_data = actor_data, 
      #                             n_steps_inclusion = n_steps_inclusion_actor,dir_path = dir_path)
      #   
      # } else {
      #   result_tmp = prepare_data_no_re(i = i,df_sub2 = df,actor_data = actor_data, 
      #                                   n_steps_inclusion_actor = n_steps_inclusion_actor,
      #                                   n_steps_inclusion_past = n_steps_inclusion_past,
      #                                   dir_path = dir_path)
      #   
      # }
      # saveRDS(paste0("../Data/preprocessing/Step_",i,".rds"),object =  result_tmp)
      # n = n +1
      gc(full = T)
      gc()
    }
  } else {
    result = parLapply(cl, time_steps_eff,prepare_data_no_re,df_sub2 = df,dir_path = dir_path)
  }
}


