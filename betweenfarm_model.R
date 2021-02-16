betweenfarm_model <- function(
  # beta transmission probabilities
  beta_net,
  beta_local,
  beta_reinfected,
  beta_truck,
  beta_formula,
  
  # initial state nodes and infectius period
  I,
  tI_time,
  R,
  tR_time,
  
  # detection
  farm_detection_ratio,
  efetive_surveillance,
  
  # force farms to be infected or suceptible in determinated time
  force_events_infected,
  force_events_suceptible,
  
  # additonal information
  tSim, # number of simulated time steps
  time_step, # unit of time from the time steps
  time_start, # start date model
  seasonality, 
  
  # nodes and networks 
  nodes,
  gravity_matrix, # matrix gravity model between nodes
  distance_matrix, # matrix distances between nodes
  movement_dataframe, # database with netwotk
  movement_matrix, # list of movement matrix by each time step
  truck_matrix,
  formula
  
){
  
  # conditional to know the month we are modelling to fit with the seasonality
  if(time_step == "1 week"){
    time_start <- time_start + seq(0, (tSim-1)*7, 7)
  } else if(time_step == "2 week"){
    time_start <- time_start + seq(0, (tSim-1)*14, 14)
  } else if(time_step == "day"){
    time_start <- time_start + seq(0, (tSim-1), 1)
  } 
  time_start <- as.numeric(format(time_start, "%m"))
  
  # matrix with the status of the farms in each day
  M_Sim_I = matrix(data = NA, nrow = nrow(nodes), ncol = (tSim+1)) #
  M_Sim_newI = matrix(data = NA, nrow = nrow(nodes), ncol = (tSim+1))
  M_Sim_R = matrix(data = NA, nrow = nrow(nodes), ncol = (tSim+1))
  M_Sim_O = matrix(data = NA, nrow = nrow(nodes), ncol = (tSim+1))
  
  
  # forced events with time equal to 0 - this mean we force farms to be infected
  I[ force_events_infected$network.id[force_events_infected$t1 == 0] ] <- 1
  
  # define which farms have initial infectius state time step 1
  tI <- I
  tR <- R
  
  # time infected: define a number of time steps to be infected
  # simulated time the intial infected farms are going to be infeted - the elapsed infected time
  time_first_infected <- (VGAM::rpospois(length(which(I != 0)), tI_time)) - tI[tI>0] 
  # sometimes the previus operation can result in 0,
  # to avoid conflict we say those farms have still 1 time step to be infected
  time_first_infected[time_first_infected <= 0] <- 1
  tI[tI>0] <- time_first_infected 
  
  # we add that information to the matrix at day 1 (which us t0)
  M_Sim_I[, 1] = tI
  M_Sim_R[, 1] = tR # here is 0 for all because there are not farms recovery yet
  M_Sim_O[, 1] = 0
  M_Sim_O[tI > 0, 1] = 1
  M_Sim_newI[, 1] = 0
  M_Sim_newI[tI > 0, 1] = 1
  
  
  new_detected_neighbords <- c()
  for (nl in  which(quarenten_vector > 0)) {
    new_detected_neighbords <- c(new_detected_neighbords, which(distance_matrix[,nl] <= distance_q))
  }
  new_detected_neighbords <- new_detected_neighbords[!new_detected_neighbords %in% which(quarenten_vector > 0)]
  new_detected_neighbords <- unique(new_detected_neighbords)
  
  # quareten initial days
  quarenten_vector[c( which(quarenten_vector > 0), new_detected_neighbords)] <- days_q
  
  trans_routes <- data.frame()
  time_infected <- as.numeric(tI>0)
  
  for (i in 1:tSim) {
    
    # network movements in time step i
    movement_matrix_timestep <- movement_matrix[[i]]
    
    # network truck  in time step i
    truck_matrix_timestep <- truck_matrix[[i]] 
    
    # formula in time step i
    formula_timestep <- formula[,i]
    
    
    I = as.numeric(tI > 0) # identified infected farms as 1
    R = as.numeric(tR > 0) # identified recovery farms as 1
    S = as.numeric(!I & !R) # identified Suceptible farms as 1
    
    # quarent
    # quarenten_vector_acumulates <- (quarenten_vector + (quarenten_vector > 0)) * (quarenten_vector != days_q)
    # I <- I * as.numeric(quarenten_vector == 0) # if 0 meas there is not quarentaine and it can infect
    quarenten_vector[which(quarenten_vector>0)] <- quarenten_vector[which(quarenten_vector>0)] - 1
    
    # plot( 0.95/(1 + exp(-0.5*(x_plogis - 15))) )
    # 0.95 maximun value, -0.5 = curve, time infected vector with the days, farm_detection_ratio mid point
    # detetect just for infeted farms
    ratedetection <- I * ((0.95/(1 + exp(-0.5*(time_infected - farm_detection_ratio))) ) * efetive_surveillance )
    
    time_infected <- time_infected + as.numeric(tI>0)
    
    # force suceptible by time have been infected
    time_infected_nursery <- time_infected[ nodes$network.id[nodes$farm_type == "Nursery"] ]
    names(time_infected_nursery) <- nodes$network.id[nodes$farm_type == "Nursery"]
    time_infected_nursery <- as.numeric(names(which(time_infected_nursery > 8)))
    
    time_infected_finisher <- time_infected[ nodes$network.id[nodes$farm_type == "Finisher"] ]
    names(time_infected_finisher) <- nodes$network.id[nodes$farm_type == "Finisher"]
    time_infected_finisher <- as.numeric(names(which(time_infected_finisher > 26)))
    
    time_infected_gilt <- time_infected[ nodes$network.id[nodes$farm_type == "Gilt"] ]
    names(time_infected_gilt) <- nodes$network.id[nodes$farm_type == "Gilt"]
    time_infected_gilt <- as.numeric(names(which(time_infected_gilt > 23)))
    
    
    time_infected <- time_infected * I # if the farm is not infected the time will retunr to 0
    
    # detect just for farms have not been detected previously
    ratedetection <- ratedetection*(M_Sim_O[, i] == 0)
    # here decided threshold to detect the farm
    threshold_detection <- runif(length(ratedetection), 0, 1)
    
    # identified new farm with outbreak
    new_outbreak <- as.numeric(ratedetection>threshold_detection)
    # force some farms not to be an outbreak
    new_outbreak[new_outbreak %in%
                   force_events_infected$network.id[force_events_infected$t2 >= i &
                                                      force_events_infected$t2 <= (tI_time - farm_detection_ratio)]
    ] <- 0
    
    # we acumulated the number of days since the farm was detected
    # if the farm is not infected anymore, it will become 0 again
    M_Sim_O[, i + 1] <- (M_Sim_O[, i] + (M_Sim_O[, i] > 0) + new_outbreak)*I
    
    # new outbreaks
    new_detected_nodes <- which(M_Sim_O[, i + 1] == 1)
    # exclude from new_detected_nodes the nodes which were identified previosly as detected
    
    if(auto_q == T){
      
      new_detected_neighbords <- c()
      for (nl in  new_detected_nodes) {
        new_detected_neighbords <- c(new_detected_neighbords, which(distance_matrix[,nl] <= distance_q))
      }
      new_detected_neighbords <- new_detected_neighbords[!new_detected_neighbords %in% new_detected_nodes]
      new_detected_neighbords <- unique(new_detected_neighbords)
      
      quarenten_vector[c(new_detected_nodes, new_detected_neighbords)] <- days_q
      # quarenten_vector[!quarenten_vector %in% allowed_quarenten] <- 0 #----- check
    }
    
    # auto vaccine outbreaks
    if(auto_v == T){
      time_vac_ratio_increase <- 1
      # time_vac_ratio_decrease <- 8
      # plot( 0.1/(1 + exp(-1.5*(1:10 - time_vac_ratio_increase))) )
      nodes_outbreaks <- which(M_Sim_O[, i + 1] > 0 & M_Sim_O[, i + 1] <= inmunity_days) 
      # nodes_outbreaks <- nodes_outbreaks[!nodes_outbreaks %in% allowed_vac] #----- check
      nodes_vac_out <- M_Sim_O[nodes_outbreaks, i + 1]
      ratevacine <-  (effectivity_control/(1 + exp(-1.5*(nodes_vac_out - time_vac_ratio_increase))) ) 
      
      # sweep(gravity_matrix[nodes_outbreaks,124:125] , MARGIN=1, ratevacine, `*`)
      # if(length( nodes_outbreaks ) > 0 ){
      # I replaced here effectivity_control by ratevacine
      movement_matrix_timestep <- as.matrix(movement_matrix_timestep)
      movement_matrix_timestep[nodes_outbreaks,] <- movement_matrix_timestep[nodes_outbreaks,] -
        sweep(movement_matrix_timestep[nodes_outbreaks,] , MARGIN=1, ratevacine, `*`)
      
      gravity_matrix[nodes_outbreaks,] <- gravity_matrix[nodes_outbreaks,] -
        sweep(gravity_matrix[nodes_outbreaks,] , MARGIN=1, ratevacine, `*`)
      
      # movement_matrix_timestep[nodes_outbreaks,] <- movement_matrix_timestep[nodes_outbreaks,] -
      #   (movement_matrix_timestep[nodes_outbreaks,] * effectivity_control )
      
      # gravity_matrix[nodes_outbreaks,] <- gravity_matrix[nodes_outbreaks,] -
      #   (gravity_matrix[nodes_outbreaks,] * effectivity_control )
      # }
    }
    
    # matrix multiplication infected farms and movements 
    # vector x columns
    # with this we know if a farm receibed animals from an infected farm
    # if column is not 0 mean that the farm recibed
    # a movement from an infected farm
    # Pvector <- I %*% movement_matrix_timestep
    # quarenten_vector 0 = T or 1 and diff 0 = F or 0, if we multiply by the status
    # 1 infected by 0 quarenten = 0 transmission
    Pvector <- (I * as.numeric(quarenten_vector == 0) ) %*% movement_matrix_timestep
    
    truck_vector <- (I * as.numeric(quarenten_vector == 0) ) %*% truck_matrix_timestep
    formula_vector <- (I < 1) * formula_timestep
    
    Pvector <- Pvector[1,] # new line
    truck_vector <- truck_vector[1,] # new line
    # example
    # matrix(c(F,T,T,F), nrow = 2, ncol = 2)
    #       [,1]  [,2]
    # [1,] FALSE  TRUE
    # [2,]  TRUE FALSE
    # c(0,1) %*%  matrix(c(F,T,T,F), nrow = 2, ncol = 2)
    #       [,1] [,2]
    # [1,]    1    0
    # 0*F = 0, 1*T = 1 --- 0 + 1 = 1
    # 0*T = 0, 1*F = 0 --- 0 + 0 = 0
    
    # here we calculate a transmission value
    # if a farm did not recebe movements from infected farms
    # the Pvertex = 0 and will not have infection there
    # but if the farms is 1 will have the same value as the
    # beta defined at the beginnig of the model
    # if a farm recebe from more than infected farm that value
    # will increase - that is defined by the Pvector
    # Pvertex = 1 - (1 - beta)^Pvector # farm infectiuos threshold
    
    
    rateTransmission_net <- (beta_net * seasonality[names(seasonality) == time_start[i]]) * Pvector # new line
    
    rateTransmission_truck <- (beta_truck * seasonality[names(seasonality) == time_start[i]]) * truck_vector
    rateTransmission_formula <- (beta_formula * seasonality[names(seasonality) == time_start[i]]) * formula_vector
    
    spatial_vector <- (I * as.numeric(quarenten_vector == 0)) %*% gravity_matrix
    spatial_vector <- spatial_vector[1,]
    spatial_vector <- spatial_vector * evi_matrix[,i]
    rateTransmission_sp<- (beta_local*seasonality[names(seasonality) == time_start[i]] ) * spatial_vector
    
    # high risk reinfected farms
    # the infected farms cannot reinfected again 
    ratereinfected <- (reinfected * (beta_reinfected* seasonality[names(seasonality) == time_start[i]] ) ) * as.numeric(!I)
    
    # farms were not infected a long time ago
    # they are different from the intial farms infected
    # they are different from the farm initial reinfected
    ratehardlyinfected <- hadrlyinfected * as.numeric(!I)
    
    
    trans_routes_aux <- data.frame(id = 1:length(movement_matrix_timestep[1, ]),
                                   t = i,
                                   net = 1 - exp(-rateTransmission_net),
                                   sp = 1 - exp(-rateTransmission_sp),
                                   rebreak = 1 - exp(-ratereinfected),
                                   truck = 1 - exp(-rateTransmission_truck),
                                   formula = 1 - exp(-rateTransmission_formula)
    )
    
    
    # sum transmission
    rateTransmission <- rateTransmission_net +
      rateTransmission_sp +
      ratereinfected +
      rateTransmission_truck+
      rateTransmission_formula
    
    
    rateTransmission <- 1 - exp(-rateTransmission) # new line
    
    # rate vacination
    ratevacination <- effectivity_control * Control # proportion animal immune
    rateTransmission <- rateTransmission - (rateTransmission * ratevacination)
    rateTransmission <- rateTransmission - (rateTransmission * ratehardlyinfected)
    
    
    
    
    # random probability a farm get infected
    Prand <- runif(n = length(movement_matrix_timestep[1, ]), min = 0, max = 1)
    
    
    # if the farm have a Pvertex higher than the prand it will be infected
    # Pinfection = as.numeric(Prand <= Pvertex)
    
    Pinfection <- as.numeric(rateTransmission >= Prand)
    
    # force farm to be infected when time > t0
    Pinfection[force_events_infected$network.id[force_events_infected$t1 == i]] <- 1
    
    # force farm to be outbreak
    # new_outbreak[force_events_infected$network.id[force_events_infected$days == i]] <- 1
    M_Sim_O[force_events_infected$network.id[force_events_infected$t2 == i], i+1] <- 1
    
    # if there is a control in the farm that will not be infected
    # Pinfection = Pinfection * !(Control)
    # new infected farm 
    Inew = as.numeric(S & (Pinfection))
    M_Sim_newI[, i + 1] = Inew
    trans_routes_aux$infected <- Inew
    trans_routes <- rbind(trans_routes, trans_routes_aux)
    
    # new recovery farms
    Rnew = rep(0, length(R))
    # if the farm infected had 1 day left that will be recovered in the next step
    if(tR_time > 0){Rnew[tI == 1] = 1} 
    ### aqui se resta um dia para de contagio
    # reduce the number days that the infected farms are infected
    temp_tI = tI - 1
    tInew = as.numeric(temp_tI < 0) + temp_tI # here farm become suceptible again if the infected time is less than 0
    # reduce the number days that recovery farms are recovery
    temp_tR = tR - 1
    tRnew = as.numeric(temp_tR < 0) + temp_tR
    
    # days infected for new infected farms
    aux_tI = rep(0, length(Inew))
    if(length(aux_tI[Inew > 0]) > 0){aux_tI[Inew > 0] <- VGAM::rpospois(length(which(Inew != 0)), tI_time) }
    
    
    # days recovered for new recovred farms
    aux_tR = rep(0, length(Rnew))
    if(tR_time > 0){
      if(length(aux_tR[Rnew > 0]) > 0){
        aux_tR[Rnew > 0] <- VGAM::rpospois(length(which(Rnew != 0)), tR_time) 
      }
    }
    
    
    tI = tInew + aux_tI # we add the new information of infected farms (days)
    tR = tRnew + aux_tR # + Control # we add information about recovery farms
    # (days), also, we can maybe add a control as days to increase the recovey
    # days
    
    # force infected to suceptible by type of farm and movements
    # if finisher or nursery and move 60% of the farm capacity it will become suceptible
    tI[force_events_suceptible[[i]]] <- 0
    tI[time_infected_nursery] <- 0
    tI[time_infected_finisher] <- 0
    tI[time_infected_gilt] <- 0
    
    # add the new results to the infected and recovery matrix time
    M_Sim_I[, i + 1] = tI
    M_Sim_R[, i + 1] = tR
    
    
    # print(paste(i, "i"))
  }
  
  answer = list(M_Sim_I, M_Sim_R, M_Sim_O, M_Sim_newI, trans_routes)
  return(answer)
}
