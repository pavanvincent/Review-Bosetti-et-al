require(tidyverse)
vect_name_region <- c('ARA', 'BFC', 'BRE', 'CVL', 'COR', 'GES', 'HDF', 'IDF', 'NAQ', 'NOR', 'OCC', 'PAC', 'PDL')

### Functions used to generate the vaccination schedule based on hypotheses regarding
### the delay between doses, the vaccine roll-out pace, and the number of doses being delivered every month.

### The script is suited to generate matrices for the following age-groups:
### 0-9 ; 10-17 ; 18-29 ; 30-39 ; 40-44 ; 45-49 ; 50-54 ; 55-59 ; 60-64 ; 65-69 ; 70-74 ; 75-79 ; 80+
age_groups_to_target <- c('0-17y', '18-49y', '50-64y', '65-74y', '75y+')
corresponding_ages <- c(rep(1, 2),
                        rep(2, 4),
                        rep(3, 3),
                        rep(4, 2),
                        rep(5, 2))

vaccine_priority_Age <- list(list(list('75y+', 0:3)),
                             list(list('65-74y', 0:3)),
                             list(list('50-64y', 0:3)),
                             list(list('18-49y', 0:3)))

vaccine_priority_Age_From6574 <- list(list(list('65-74y', 0:3)),
                                      list(list('50-64y', 0:3)),
                                      list(list('18-49y', 0:3)))

vaccine_priority_Age_From5064 <- list(list(list('50-64y', 0:3)),
                                      list(list('18-49y', 0:3)))

vaccine_priority_Age_Just50p <- list(list(list('75y+', 0:3)),
                                     list(list('65-74y', 0:3)),
                                     list(list('50-64y', 0:3)))


### Functions used assuming no change in the roll-out pace
get_nbFirstDosesPerDay_nVaccines <- function(df_Doses,
                                             seq_n_maxDosesPerDay,
                                             date_beginning_vaccine = as.Date('2021-02-01'),
                                             date_end_vaccine = as.Date('2021-12-31'),
                                             seq_DelayBetweenDoses){
  
  # df_Doses is a dataframe with the following columns:
  ### DateDelivery : the date at which vaccines doses will be delivered
  ### Doses_i : the number of doses of vaccine i for by date of delivery
  
  seq_Dates <- seq.Date(from = date_beginning_vaccine,
                        to = date_end_vaccine,
                        by = 'day')
  
  nVaccines <- length(seq_n_maxDosesPerDay)
  
  list_seq_FirstDoses <- lapply(1:nVaccines, FUN = function(iVaccine){
    rep(0, length(seq_Dates))
  })
  list_seq_SecondDoses <- lapply(1:nVaccines, FUN = function(iVaccine){
    rep(0, length(seq_Dates) + seq_DelayBetweenDoses[iVaccine])
  })
  
  list_seq_DosesAvailable <- lapply(1:nVaccines, FUN = function(iVaccine){
    sapply(seq_Dates, FUN = function(tmp_date){
      DosesAlreadyDistributed <- sum(df_Doses[df_Doses$DateDelivery <= tmp_date, paste0('Doses_', iVaccine)])
    })
  })  # Vector where you store the available doses on a given day (remove second doses)
 
  list_seq_DosesAvailable <- lapply(1:nVaccines, FUN = function(iVaccine){
    c(list_seq_DosesAvailable[[iVaccine]], rep(1e8, seq_DelayBetweenDoses[iVaccine]))
  })
  
  
  for(i in 1:length(seq_Dates)){
    tmp_date <- seq_Dates[i]
    list_tmpSecondDosesToDistribute <- lapply(1:nVaccines, FUN = function(iVaccine){
      list_seq_SecondDoses[[iVaccine]][i]
    })
    ## Get a number of doses that can be distributed (and checking that there will be enough)
    list_tmpDosesAvailable <- lapply(1:nVaccines, FUN = function(iVaccine){
      tmp_seqDosesAvailable <- list_seq_DosesAvailable[[iVaccine]]
      min(tmp_seqDosesAvailable[i:length(tmp_seqDosesAvailable)])
    })
    list_tmpSecondDosesAvailable <- lapply(1:nVaccines, FUN = function(iVaccine){
      tmp_seqDosesAvailable <- list_seq_DosesAvailable[[iVaccine]]
      min(tmp_seqDosesAvailable[(i + seq_DelayBetweenDoses[iVaccine]):length(tmp_seqDosesAvailable)]/2)
    })
    list_tmpDosesThatCanBeDistributed <- lapply(1:nVaccines, FUN = function(iVaccine){
      seq_n_maxDosesPerDay[iVaccine] - list_tmpSecondDosesToDistribute[[iVaccine]]
    })
    list_tmp_FirstDosesToDistribute <- lapply(1:nVaccines, FUN = function(iVaccine){
      min(list_tmpDosesAvailable[[iVaccine]],
          list_tmpSecondDosesAvailable[[iVaccine]],
          list_tmpDosesThatCanBeDistributed[[iVaccine]])
    })
    for(iVaccine in 1:nVaccines){
      if(seq_DelayBetweenDoses[iVaccine] == 0){ # Single-dose vaccine
        list_seq_FirstDoses[[iVaccine]][[i]] <- list_tmp_FirstDosesToDistribute[[iVaccine]]
        
        list_seq_DosesAvailable[[iVaccine]][i:length(list_seq_DosesAvailable[[iVaccine]])] <- 
          list_seq_DosesAvailable[[iVaccine]][i:length(list_seq_DosesAvailable[[iVaccine]])] - 
          list_tmp_FirstDosesToDistribute[[iVaccine]]
        
        
      } else{
        list_seq_FirstDoses[[iVaccine]][[i]] <- list_tmp_FirstDosesToDistribute[[iVaccine]]
        list_seq_SecondDoses[[iVaccine]][i + seq_DelayBetweenDoses[iVaccine]] <- list_tmp_FirstDosesToDistribute[[iVaccine]]
        
        list_seq_DosesAvailable[[iVaccine]][i:length(list_seq_DosesAvailable[[iVaccine]])] <- 
          list_seq_DosesAvailable[[iVaccine]][i:length(list_seq_DosesAvailable[[iVaccine]])] - 
          list_tmp_FirstDosesToDistribute[[iVaccine]]
        list_seq_DosesAvailable[[iVaccine]][(i + seq_DelayBetweenDoses[iVaccine]):(length(list_seq_DosesAvailable[[iVaccine]]))] <- 
          list_seq_DosesAvailable[[iVaccine]][(i + seq_DelayBetweenDoses[iVaccine]):(length(list_seq_DosesAvailable[[iVaccine]]))] - 
          list_tmp_FirstDosesToDistribute[[iVaccine]]
      }
      
    }
  }
  
  list_schema <- list(seq_Dates = seq_Dates,
                      list_seq_FirstDoses = list_seq_FirstDoses,
                      list_seq_SecondDoses = list_seq_SecondDoses,
                      list_seq_DosesAvailable = list_seq_DosesAvailable)
                    
  return(list_schema)
}

get_Mat_VaccineStrategy_nVaccines_withComorbidities <- function(date_beginning_vaccine,
                                                                date_end_vaccine,
                                                                vect_maxVC,
                                                                seq_n_maxDosesPerDay,
                                                                df_Doses,
                                                                seq_DelayBetweenDoses,
                                                                list_vaccine_priority,
                                                                N_0, N_1, N_2, N_3){
  ## Defining the matrix
  time_beginning_vaccine <- as.numeric(date_beginning_vaccine)
  time_end_vaccine <- as.numeric(date_end_vaccine)
  
  nVaccines <- length(seq_n_maxDosesPerDay)
  n.ageG <- length(vect_maxVC)
  
  
  list_FirstDoses <- get_nbFirstDosesPerDay_nVaccines(df_Doses,
                                                      seq_n_maxDosesPerDay, 
                                                      date_beginning_vaccine, date_end_vaccine,
                                                      seq_DelayBetweenDoses)
  
  list_vect_nMaxVaccinePerDay <- list_FirstDoses[["list_seq_FirstDoses"]]
  
  list_Mat_VaccineStrategy <- lapply(1:nVaccines, FUN = function(iVaccine){
    Mat_VaccineStrategy <- matrix(0, nrow = time_end_vaccine - time_beginning_vaccine + 1,
                                  ncol = 4*n.ageG)
    
  })
  ## Maximum number of individuals that can be vaccinated in the different age-comorbidities groups
  max_IndivToVaccinate <- c(N_0*vect_maxVC, N_1*vect_maxVC, N_2*vect_maxVC, N_3*vect_maxVC)
  
  list_nGroupsToVaccine <- lapply(1:nVaccines, FUN = function(iVaccine){
    length(list_vaccine_priority[[iVaccine]])
  })
  list_ColumnIndicesMat_VaccineStrategy <- lapply(1:nVaccines, FUN = function(iVaccine){
    lapply(list_vaccine_priority[[iVaccine]], FUN = function(l1){
      sort(Reduce('c', lapply(l1, FUN = function(l){
        tmp_comorbidities <- l[[2]]
        tmp_age <- which(age_groups_to_target ==  l[1])
        tmp_indices_age <- which(corresponding_ages == tmp_age)
        Reduce('c',lapply(tmp_comorbidities, FUN = function(tmp_comorbidity){
          tmp_indices_age + (tmp_comorbidity)*n.ageG
        }))
      })))
    })
  })
  list_RepartitionVaccineWithinGroups <- lapply(1:nVaccines, FUN = function(iVaccine){
    lapply(list_ColumnIndicesMat_VaccineStrategy[[iVaccine]], FUN = function(tmp_indices){
      tmp_num <- max_IndivToVaccinate[tmp_indices]
      tmp_den <- sum(max_IndivToVaccinate[tmp_indices])
      if(tmp_den == 0){
        tmp_num
      } else{
        tmp_num/tmp_den
      }
    })
  })
  list_MaxIndivToVaccinateWithinGroups <- lapply(1:nVaccines, FUN = function(iVaccine){
    lapply(list_ColumnIndicesMat_VaccineStrategy[[iVaccine]], FUN = function(tmp_indices){
      max_IndivToVaccinate[tmp_indices]
    })
  })
  seq_N.AllocatedVaccines <- rep(0, nVaccines)
  seq_nVaccinePerDay <- sapply(1:nVaccines, FUN = function(iVaccine){
    list_vect_nMaxVaccinePerDay[[iVaccine]][1]
  })
  
  t <- 1
  seq_index_groupToVaccinate <- sapply(1:nVaccines, FUN = function(iVaccine){
    1
  })
  
  
  seq_runvacc <- sapply(1:nVaccines, FUN = function(iVaccine){
    TRUE
  })# Boolean storing whether you finished or not the vaccination of the given group on the given day
  seq_change <- rep(T, nVaccines)
  
  
  while(((t <= nrow(list_Mat_VaccineStrategy[[1]])) &
         (sum(seq_runvacc) >=1)
  )){
    ## Number of doses already allocated to the age group we consider
    if(t == 1){
      list_nDosesAlreadyAllocated_All_PerVacc <- lapply(1:nVaccines, FUN = function(iVaccine){
        list_Mat_VaccineStrategy[[iVaccine]][t, ]
      })
      
    } else{
      list_nDosesAlreadyAllocated_All_PerVacc <- lapply(1:nVaccines, FUN = function(iVaccine){
        apply(list_Mat_VaccineStrategy[[iVaccine]][1:t, ], 2, sum)
      })
      
    }
    nDosesAlreadyAllocated_All <- Reduce('+', list_nDosesAlreadyAllocated_All_PerVacc)
    
    
    for(iVaccine in 1:nVaccines){
      tmp_run_vacc <- seq_runvacc[iVaccine]
      if(tmp_run_vacc){
        ## Number of doses to distribute
        tmp_nDosesToAllocate <- seq_nVaccinePerDay[iVaccine]
        ## Repartition of the doses to allocate
        tmp_index_groupToVaccinate <- seq_index_groupToVaccinate[iVaccine]
        tmp_Current_RepartitionAllocation <- list_RepartitionVaccineWithinGroups[[iVaccine]][[tmp_index_groupToVaccinate]]
        tmp_Current_MaxIndivToVaccinate <- list_MaxIndivToVaccinateWithinGroups[[iVaccine]][[tmp_index_groupToVaccinate]]
        
        tmp_nDosesAlreadyAllocated <- nDosesAlreadyAllocated_All[(list_ColumnIndicesMat_VaccineStrategy[[iVaccine]])[[tmp_index_groupToVaccinate]]]
        
        tmp_nDosesToAllocateThisGroupThisDay <- ifelse(tmp_Current_MaxIndivToVaccinate - tmp_nDosesAlreadyAllocated - tmp_nDosesToAllocate*tmp_Current_RepartitionAllocation < 0,
                                                       tmp_Current_MaxIndivToVaccinate - tmp_nDosesAlreadyAllocated,
                                                       tmp_nDosesToAllocate*tmp_Current_RepartitionAllocation)
        list_Mat_VaccineStrategy[[iVaccine]][t, list_ColumnIndicesMat_VaccineStrategy[[iVaccine]][[tmp_index_groupToVaccinate]]] <- 
          list_Mat_VaccineStrategy[[iVaccine]][t, list_ColumnIndicesMat_VaccineStrategy[[iVaccine]][[tmp_index_groupToVaccinate]]] +
          tmp_nDosesToAllocateThisGroupThisDay
        
        seq_N.AllocatedVaccines[iVaccine] <- seq_N.AllocatedVaccines[iVaccine] + sum(tmp_nDosesToAllocateThisGroupThisDay)
        
        nDosesAlreadyAllocated_All[list_ColumnIndicesMat_VaccineStrategy[[iVaccine]][[tmp_index_groupToVaccinate]]] <-  
          nDosesAlreadyAllocated_All[list_ColumnIndicesMat_VaccineStrategy[[iVaccine]][[tmp_index_groupToVaccinate]]] +
          tmp_nDosesToAllocateThisGroupThisDay
        
        
        if(abs(sum(tmp_nDosesToAllocateThisGroupThisDay) -  seq_nVaccinePerDay[iVaccine]) < 1e-8){
          seq_nVaccinePerDay[iVaccine] <- list_vect_nMaxVaccinePerDay[[iVaccine]][t + 1]
          seq_change[iVaccine] <- F
        } else{
          seq_nVaccinePerDay[iVaccine] <- seq_nVaccinePerDay[iVaccine] - sum(tmp_nDosesToAllocateThisGroupThisDay)
          seq_index_groupToVaccinate[iVaccine] <- tmp_index_groupToVaccinate + 1
          seq_change[iVaccine] <- T
        }
      }
    }
    if(sum(seq_change) == 0){
      t <- t + 1
    }
    seq_runvacc <- sapply(1:nVaccines, FUN = function(iVaccine){
      (seq_change[iVaccine] || sum(seq_change) == 0)& # Do you need to run it for another iteration if t did not increase
        (t <= (time_end_vaccine - time_beginning_vaccine + 1) & # Are you by the end of the simulation ?
           (seq_index_groupToVaccinate[iVaccine] <= list_nGroupsToVaccine[[iVaccine]]))  # Is there more groups to vaccinate
    })
    
  }
  
  list_Mat <- list_Mat_VaccineStrategy
  names(list_Mat) <- paste0('Mat_VaccineStrategy_', 1:nVaccines)
  return(list_Mat)
}

get_Mat_VaccineStrategy_nVaccines_justAge <- function(date_beginning_vaccine,
                                                      date_end_vaccine,
                                                      vect_maxVC,
                                                      seq_n_maxDosesPerDay,
                                                      df_Doses,
                                                      seq_DelayBetweenDoses,
                                                      list_vaccine_priority,
                                                      popSize.ageG){
  
  n.ageG <- length(popSize.ageG)
  N_0 <- popSize.ageG
  N_1 <- N_2 <- N_3 <- rep(0, n.ageG)
  
  list_Mat <- get_Mat_VaccineStrategy_nVaccines_withComorbidities(date_beginning_vaccine,
                                                                  date_end_vaccine,
                                                                  vect_maxVC,
                                                                  seq_n_maxDosesPerDay,
                                                                  df_Doses,
                                                                  seq_DelayBetweenDoses,
                                                                  list_vaccine_priority,
                                                                  N_0, N_1, N_2, N_3)
  
  list_Mat_AgeOnly <- lapply(list_Mat, FUN = function(tmp_Mat){
    sapply(1:n.ageG, FUN = function(tmp_age){
      round(tmp_Mat[, tmp_age] + tmp_Mat[, tmp_age + n.ageG] +
              tmp_Mat[, tmp_age + 2*n.ageG] + tmp_Mat[, tmp_age + 3*n.ageG],
            4)
    })
  })
  
  
  return(list_Mat_AgeOnly)
}


get_Mat_Rates_nVaccines_justAge <- function(date_beginning_vaccine,
                                            date_end_vaccine,
                                            vect_maxVC,
                                            seq_n_maxDosesPerDay,
                                            df_Doses,
                                            seq_DelayBetweenDoses,
                                            list_vaccine_priority,
                                            popSize.ageG){
  
  n.ageG <- length(popSize.ageG)
  
  list_Mat_AgeOnly <- get_Mat_VaccineStrategy_nVaccines_justAge(date_beginning_vaccine,
                                                                date_end_vaccine,
                                                                vect_maxVC,
                                                                seq_n_maxDosesPerDay,
                                                                df_Doses,
                                                                seq_DelayBetweenDoses,
                                                                list_vaccine_priority,
                                                                popSize.ageG)
  
  nVaccines <- length(list_Mat_AgeOnly)
  nDays <- nrow(list_Mat_AgeOnly$Mat_VaccineStrategy_1)
  
  
  
  list_Mat_Rates <- lapply(1:nVaccines, FUN = function(i){
    matrix(0, ncol = n.ageG, nrow = nDays)
  })
  names(list_Mat_Rates) <- names(list_Mat_AgeOnly)
  nNonVaccinated.ageG <- popSize.ageG
  for(iDay in 1:nDays){
    list_VectToVaccinate <- lapply(list_Mat_AgeOnly, FUN = function(l){
      l[iDay,]
    })
    list_Vectrates <- lapply(list_VectToVaccinate, FUN = function(l){
      nNewNonVaccinated.ageG <- nNonVaccinated.ageG - l
      log(nNonVaccinated.ageG) - log(nNewNonVaccinated.ageG)
    })
    ## Update the number of non vaccinated individuals in the different age groups
    nNonVaccinated.ageG <- nNonVaccinated.ageG - Reduce('+', list_VectToVaccinate)
    ## Update the list with the rates
    for(iVaccine in 1:nVaccines){
      list_Mat_Rates[[paste0('Mat_VaccineStrategy_', iVaccine)]][iDay, ] <- list_Vectrates[[paste0('Mat_VaccineStrategy_', iVaccine)]] 
    }
  }
  return(list_Mat_Rates)
}



### Functions used when using a change in the roll-out pace for the different vaccines
get_nbFirstDosesPerDay_nVaccines_ChancePace <- function(df_Doses,
                                                        seq_n_maxDosesPerDay_Before,
                                                        seq_n_maxDosesPerDay_After,
                                                        date_beginning_vaccine = as.Date('2021-02-01'),
                                                        date_end_vaccine = as.Date('2021-12-31'),
                                                        seq_DelayBetweenDoses,
                                                        seq_DateChangePace
                                                        ){
  
  # df_Doses is a dataframe with the following columns:
  ### DateDelivery : the date at which vaccines doses will be delivered
  ### Doses_i : the number of doses of vaccine i for by date of delivery
  
  seq_Dates <- seq.Date(from = date_beginning_vaccine,
                        to = date_end_vaccine,
                        by = 'day')
  
  nVaccines <- length(seq_n_maxDosesPerDay_Before)
  
  list_seq_FirstDoses <- lapply(1:nVaccines, FUN = function(iVaccine){
    rep(0, length(seq_Dates))
  })
  list_seq_SecondDoses <- lapply(1:nVaccines, FUN = function(iVaccine){
    rep(0, length(seq_Dates) + seq_DelayBetweenDoses[iVaccine])
  })
  
  list_seq_DosesAvailable <- lapply(1:nVaccines, FUN = function(iVaccine){
    sapply(seq_Dates, FUN = function(tmp_date){
      DosesAlreadyDistributed <- sum(df_Doses[df_Doses$DateDelivery <= tmp_date, paste0('Doses_', iVaccine)])
    })
  })  # Vector where you store the available doses on a given day (remove second doses)
  
  list_seq_DosesAvailable <- lapply(1:nVaccines, FUN = function(iVaccine){
    c(list_seq_DosesAvailable[[iVaccine]], rep(1e8, seq_DelayBetweenDoses[iVaccine]))
  })
  list_vect_n_maxDosesPerDay <- lapply(1:nVaccines, FUN = function(iVaccine){
    tmp_date_change <- seq_DateChangePace[iVaccine]
    tmp_nDoses_Before <- seq_n_maxDosesPerDay_Before[iVaccine]
    tmp_nDoses_After <- seq_n_maxDosesPerDay_After[iVaccine]
    
    ifelse(seq_Dates < tmp_date_change, tmp_nDoses_Before, tmp_nDoses_After)
  })
  
  
  for(i in 1:length(seq_Dates)){
    seq_n_maxDosesPerDay <- sapply(1:nVaccines, FUN = function(iVaccine){
      list_vect_n_maxDosesPerDay[[iVaccine]][i]
    })
    tmp_date <- seq_Dates[i]
    list_tmpSecondDosesToDistribute <- lapply(1:nVaccines, FUN = function(iVaccine){
      list_seq_SecondDoses[[iVaccine]][i]
    })
    ## Get a number of doses that can be distributed (and checking that there will be enough)
    list_tmpDosesAvailable <- lapply(1:nVaccines, FUN = function(iVaccine){
      tmp_seqDosesAvailable <- list_seq_DosesAvailable[[iVaccine]]
      min(tmp_seqDosesAvailable[i:length(tmp_seqDosesAvailable)])
    })
    list_tmpSecondDosesAvailable <- lapply(1:nVaccines, FUN = function(iVaccine){
      tmp_seqDosesAvailable <- list_seq_DosesAvailable[[iVaccine]]
      min(tmp_seqDosesAvailable[(i + seq_DelayBetweenDoses[iVaccine]):length(tmp_seqDosesAvailable)]/2)
    })
    list_tmpDosesThatCanBeDistributed <- lapply(1:nVaccines, FUN = function(iVaccine){
      seq_n_maxDosesPerDay[iVaccine] - list_tmpSecondDosesToDistribute[[iVaccine]]
    })
    list_tmp_FirstDosesToDistribute <- lapply(1:nVaccines, FUN = function(iVaccine){
      min(list_tmpDosesAvailable[[iVaccine]],
          list_tmpSecondDosesAvailable[[iVaccine]],
          list_tmpDosesThatCanBeDistributed[[iVaccine]])
    })
    for(iVaccine in 1:nVaccines){
      list_seq_FirstDoses[[iVaccine]][[i]] <- list_tmp_FirstDosesToDistribute[[iVaccine]]
      list_seq_SecondDoses[[iVaccine]][i + seq_DelayBetweenDoses[iVaccine]] <- list_tmp_FirstDosesToDistribute[[iVaccine]]
      
      list_seq_DosesAvailable[[iVaccine]][i:length(list_seq_DosesAvailable[[iVaccine]])] <- 
        list_seq_DosesAvailable[[iVaccine]][i:length(list_seq_DosesAvailable[[iVaccine]])] - 
        list_tmp_FirstDosesToDistribute[[iVaccine]]
      list_seq_DosesAvailable[[iVaccine]][(i + seq_DelayBetweenDoses[iVaccine]):(length(list_seq_DosesAvailable[[iVaccine]]))] <- 
        list_seq_DosesAvailable[[iVaccine]][(i + seq_DelayBetweenDoses[iVaccine]):(length(list_seq_DosesAvailable[[iVaccine]]))] - 
        list_tmp_FirstDosesToDistribute[[iVaccine]]
    }
  }
  
  list_schema <- list(seq_Dates = seq_Dates,
                      list_seq_FirstDoses = list_seq_FirstDoses,
                      list_seq_SecondDoses = list_seq_SecondDoses,
                      list_seq_DosesAvailable = list_seq_DosesAvailable)
  
  return(list_schema)
}

get_Mat_VaccineStrategy_nVaccines_withComorbidities_ChangePace <- function(date_beginning_vaccine,
                                                                          date_end_vaccine,
                                                                          vect_maxVC,
                                                                          seq_n_maxDosesPerDay_Before,
                                                                          seq_n_maxDosesPerDay_After,
                                                                          df_Doses,
                                                                          seq_DelayBetweenDoses,
                                                                          seq_DateChangePace,
                                                                          list_vaccine_priority,
                                                                          N_0, N_1, N_2, N_3){
  ## Defining the matrix
  time_beginning_vaccine <- as.numeric(date_beginning_vaccine)
  time_end_vaccine <- as.numeric(date_end_vaccine)
  
  nVaccines <- length(seq_n_maxDosesPerDay_Before)
  n.ageG <- length(vect_maxVC)
  
  list_FirstDoses <- get_nbFirstDosesPerDay_nVaccines_ChancePace(df_Doses,
                                                                 seq_n_maxDosesPerDay_Before,
                                                                 seq_n_maxDosesPerDay_After,
                                                                 date_beginning_vaccine, date_end_vaccine,
                                                                 seq_DelayBetweenDoses,
                                                                 seq_DateChangePace)
  
  list_vect_nMaxVaccinePerDay <- list_FirstDoses[["list_seq_FirstDoses"]]
  
  list_Mat_VaccineStrategy <- lapply(1:nVaccines, FUN = function(iVaccine){
    Mat_VaccineStrategy <- matrix(0, nrow = time_end_vaccine - time_beginning_vaccine + 1,
                                  ncol = 4*n.ageG)
    
  })
  ## Maximum number of individuals that can be vaccinated in the different age-comorbidities groups
  max_IndivToVaccinate <- c(N_0*vect_maxVC, N_1*vect_maxVC, N_2*vect_maxVC, N_3*vect_maxVC)
  
  list_nGroupsToVaccine <- lapply(1:nVaccines, FUN = function(iVaccine){
    length(list_vaccine_priority[[iVaccine]])
  })
  list_ColumnIndicesMat_VaccineStrategy <- lapply(1:nVaccines, FUN = function(iVaccine){
    lapply(list_vaccine_priority[[iVaccine]], FUN = function(l1){
      sort(Reduce('c', lapply(l1, FUN = function(l){
        tmp_comorbidities <- l[[2]]
        tmp_age <- which(age_groups_to_target ==  l[1])
        tmp_indices_age <- which(corresponding_ages == tmp_age)
        Reduce('c',lapply(tmp_comorbidities, FUN = function(tmp_comorbidity){
          tmp_indices_age + (tmp_comorbidity)*n.ageG
        }))
      })))
    })
  })
  list_RepartitionVaccineWithinGroups <- lapply(1:nVaccines, FUN = function(iVaccine){
    lapply(list_ColumnIndicesMat_VaccineStrategy[[iVaccine]], FUN = function(tmp_indices){
      tmp_num <- max_IndivToVaccinate[tmp_indices]
      tmp_den <- sum(max_IndivToVaccinate[tmp_indices])
      if(tmp_den == 0){
        tmp_num
      } else{
        tmp_num/tmp_den
      }
    })
  })
  list_MaxIndivToVaccinateWithinGroups <- lapply(1:nVaccines, FUN = function(iVaccine){
    lapply(list_ColumnIndicesMat_VaccineStrategy[[iVaccine]], FUN = function(tmp_indices){
      max_IndivToVaccinate[tmp_indices]
    })
  })
  seq_N.AllocatedVaccines <- rep(0, nVaccines)
  seq_nVaccinePerDay <- sapply(1:nVaccines, FUN = function(iVaccine){
    list_vect_nMaxVaccinePerDay[[iVaccine]][1]
  })
  
  t <- 1
  seq_index_groupToVaccinate <- sapply(1:nVaccines, FUN = function(iVaccine){
    1
  })
  
  
  seq_runvacc <- sapply(1:nVaccines, FUN = function(iVaccine){
    TRUE
  })# Boolean storing whether you finished or not the vaccination of the given group on the given day
  seq_change <- rep(T, nVaccines)
  
  
  while(((t <= nrow(list_Mat_VaccineStrategy[[1]])) &
         (sum(seq_runvacc) >=1)
  )){
    ## Number of doses already allocated to the age group we consider
    if(t == 1){
      list_nDosesAlreadyAllocated_All_PerVacc <- lapply(1:nVaccines, FUN = function(iVaccine){
        list_Mat_VaccineStrategy[[iVaccine]][t, ]
      })
      
    } else{
      list_nDosesAlreadyAllocated_All_PerVacc <- lapply(1:nVaccines, FUN = function(iVaccine){
        apply(list_Mat_VaccineStrategy[[iVaccine]][1:t, ], 2, sum)
      })
      
    }
    nDosesAlreadyAllocated_All <- Reduce('+', list_nDosesAlreadyAllocated_All_PerVacc)
    
    
    for(iVaccine in 1:nVaccines){
      tmp_run_vacc <- seq_runvacc[iVaccine]
      if(tmp_run_vacc){
        ## Number of doses to distribute
        tmp_nDosesToAllocate <- seq_nVaccinePerDay[iVaccine]
        ## Repartition of the doses to allocate
        tmp_index_groupToVaccinate <- seq_index_groupToVaccinate[iVaccine]
        tmp_Current_RepartitionAllocation <- list_RepartitionVaccineWithinGroups[[iVaccine]][[tmp_index_groupToVaccinate]]
        tmp_Current_MaxIndivToVaccinate <- list_MaxIndivToVaccinateWithinGroups[[iVaccine]][[tmp_index_groupToVaccinate]]
        
        tmp_nDosesAlreadyAllocated <- nDosesAlreadyAllocated_All[(list_ColumnIndicesMat_VaccineStrategy[[iVaccine]])[[tmp_index_groupToVaccinate]]]
        
        tmp_nDosesToAllocateThisGroupThisDay <- ifelse(tmp_Current_MaxIndivToVaccinate - tmp_nDosesAlreadyAllocated - tmp_nDosesToAllocate*tmp_Current_RepartitionAllocation < 0,
                                                       tmp_Current_MaxIndivToVaccinate - tmp_nDosesAlreadyAllocated,
                                                       tmp_nDosesToAllocate*tmp_Current_RepartitionAllocation)
        list_Mat_VaccineStrategy[[iVaccine]][t, list_ColumnIndicesMat_VaccineStrategy[[iVaccine]][[tmp_index_groupToVaccinate]]] <- 
          list_Mat_VaccineStrategy[[iVaccine]][t, list_ColumnIndicesMat_VaccineStrategy[[iVaccine]][[tmp_index_groupToVaccinate]]] +
          tmp_nDosesToAllocateThisGroupThisDay
        
        seq_N.AllocatedVaccines[iVaccine] <- seq_N.AllocatedVaccines[iVaccine] + sum(tmp_nDosesToAllocateThisGroupThisDay)
        
        nDosesAlreadyAllocated_All[list_ColumnIndicesMat_VaccineStrategy[[iVaccine]][[tmp_index_groupToVaccinate]]] <-  
          nDosesAlreadyAllocated_All[list_ColumnIndicesMat_VaccineStrategy[[iVaccine]][[tmp_index_groupToVaccinate]]] +
          tmp_nDosesToAllocateThisGroupThisDay
        
        
        if(abs(sum(tmp_nDosesToAllocateThisGroupThisDay) -  seq_nVaccinePerDay[iVaccine]) < 1e-8){
          seq_nVaccinePerDay[iVaccine] <- list_vect_nMaxVaccinePerDay[[iVaccine]][t + 1]
          seq_change[iVaccine] <- F
        } else{
          seq_nVaccinePerDay[iVaccine] <- seq_nVaccinePerDay[iVaccine] - sum(tmp_nDosesToAllocateThisGroupThisDay)
          seq_index_groupToVaccinate[iVaccine] <- tmp_index_groupToVaccinate + 1
          seq_change[iVaccine] <- T
        }
      }
    }
    seq_SimulationOver <- sapply(1:nVaccines, FUN = function(iVaccine){
      (seq_index_groupToVaccinate[iVaccine] > list_nGroupsToVaccine[[iVaccine]])
    }) 
    
    if(sum(seq_change[! seq_SimulationOver]) == 0){
      t <- t + 1
    }
    seq_runvacc <- sapply(1:nVaccines, FUN = function(iVaccine){
      (seq_change[iVaccine] || sum(seq_change[! seq_SimulationOver]) == 0) & # Do you need to run it for another iteration if t did not increase
        (t <= (time_end_vaccine - time_beginning_vaccine + 1) & # Are you by the end of the simulation ?
           (seq_index_groupToVaccinate[iVaccine] <= list_nGroupsToVaccine[[iVaccine]]))  # Is there more groups to vaccinate
    })
    
  }
  
  list_Mat <- list_Mat_VaccineStrategy
  names(list_Mat) <- paste0('Mat_VaccineStrategy_', 1:nVaccines)
  return(list_Mat)
}

get_Mat_VaccineStrategy_nVaccines_justAge_ChangePace <- function(date_beginning_vaccine,
                                                                 date_end_vaccine,
                                                                 vect_maxVC,
                                                                 seq_n_maxDosesPerDay_Before,
                                                                 seq_n_maxDosesPerDay_After,
                                                                 df_Doses,
                                                                 seq_DelayBetweenDoses,
                                                                 seq_DateChangePace,
                                                                 list_vaccine_priority,
                                                                 popSize.ageG){
  
  n.ageG <- length(popSize.ageG)
  N_0 <- popSize.ageG
  N_1 <- N_2 <- N_3 <- rep(0, n.ageG)
  
  list_Mat <- get_Mat_VaccineStrategy_nVaccines_withComorbidities_ChangePace(date_beginning_vaccine,
                                                                             date_end_vaccine,
                                                                             vect_maxVC,
                                                                             seq_n_maxDosesPerDay_Before,
                                                                             seq_n_maxDosesPerDay_After,
                                                                             df_Doses,
                                                                             seq_DelayBetweenDoses,
                                                                             seq_DateChangePace,
                                                                             list_vaccine_priority,
                                                                             N_0, N_1, N_2, N_3)
  
  list_Mat_AgeOnly <- lapply(list_Mat, FUN = function(tmp_Mat){
    sapply(1:n.ageG, FUN = function(tmp_age){
      round(tmp_Mat[, tmp_age] + tmp_Mat[, tmp_age + n.ageG] +
              tmp_Mat[, tmp_age + 2*n.ageG] + tmp_Mat[, tmp_age + 3*n.ageG],
            4)
    })
  })
  
  
  return(list_Mat_AgeOnly)
}

get_Mat_Rates_nVaccines_justAge_ChangePace <- function(date_beginning_vaccine,
                                                       date_end_vaccine,
                                                       vect_maxVC,
                                                       seq_n_maxDosesPerDay_Before,
                                                       seq_n_maxDosesPerDay_After,
                                                       df_Doses,
                                                       seq_DelayBetweenDoses,
                                                       seq_DateChangePace,
                                                       list_vaccine_priority,
                                                       popSize.ageG){
  
  n.ageG <- length(popSize.ageG)
  
  list_Mat_AgeOnly <- get_Mat_VaccineStrategy_nVaccines_justAge_ChangePace(date_beginning_vaccine,
                                                                           date_end_vaccine,
                                                                           vect_maxVC,
                                                                           seq_n_maxDosesPerDay_Before,
                                                                           seq_n_maxDosesPerDay_After,
                                                                           df_Doses,
                                                                           seq_DelayBetweenDoses,
                                                                           seq_DateChangePace,
                                                                           list_vaccine_priority,
                                                                           popSize.ageG)
  
  nVaccines <- length(list_Mat_AgeOnly)
  nDays <- nrow(list_Mat_AgeOnly$Mat_VaccineStrategy_1)
  
  
  
  list_Mat_Rates <- lapply(1:nVaccines, FUN = function(i){
    matrix(0, ncol = n.ageG, nrow = nDays)
  })
  names(list_Mat_Rates) <- names(list_Mat_AgeOnly)
  nNonVaccinated.ageG <- popSize.ageG
  for(iDay in 1:nDays){
    list_VectToVaccinate <- lapply(list_Mat_AgeOnly, FUN = function(l){
      l[iDay,]
    })
    list_Vectrates <- lapply(list_VectToVaccinate, FUN = function(l){
      nNewNonVaccinated.ageG <- nNonVaccinated.ageG - l
      log(nNonVaccinated.ageG) - log(nNewNonVaccinated.ageG)
    })
    ## Update the number of non vaccinated individuals in the different age groups
    nNonVaccinated.ageG <- nNonVaccinated.ageG - Reduce('+', list_VectToVaccinate)
    ## Update the list with the rates
    for(iVaccine in 1:nVaccines){
      list_Mat_Rates[[paste0('Mat_VaccineStrategy_', iVaccine)]][iDay, ] <- list_Vectrates[[paste0('Mat_VaccineStrategy_', iVaccine)]] 
    }
  }
  return(list_Mat_Rates)
}




################## Scripts to get the real vaccination schedule
get_Retrospective_Mat_VaccinationStrategy <- function(my_region_name){
  
  ## Loading the VACSI data
  date_file_VACSI <- as.Date('2021-04-17')
  df_VACSI <- read.csv2(paste0('Data/VACSI/', format(date_file_VACSI, '%Y%m%d'), '/vacsi-a-reg-', date_file_VACSI, '-19h10.csv'))
  
  df_VACSI$jour <- as.Date(df_VACSI$jour)
  INSEE_codes <- read.csv2('Data/population_data/INSEE_Code_RegionNames.csv')
  df_VACSI$region_name <- sapply(df_VACSI$reg, FUN = function(tmp_reg){
    tmp_x <- INSEE_codes[INSEE_codes$insee_code == tmp_reg, 'region_name']
  })
  
  
  
  if(my_region_name == 'Metro'){
    
    tmp_df <- df_VACSI %>% filter(region_name %in% vect_name_region, clage_vacsi !=0) %>%
      group_by(jour, clage_vacsi) %>%
      summarise(n_dose1_sum = sum(n_dose1)) %>%
      select(jour, n_dose1_sum, clage_vacsi) %>%
      ungroup() %>% group_by(clage_vacsi) %>%
      select(jour, n_dose1_sum) %>%
      group_split() %>%
      reduce(left_join, by = 'jour')
    
    
  } else{
    tmp_df <- df_VACSI %>% filter(region_name == my_region_name, clage_vacsi !=0) %>%
      group_by(clage_vacsi) %>%
      select(jour, n_dose1) %>%
      group_split() %>%
      reduce(left_join, by = 'jour')
    
  }
  
  name_ages <- tmp_df  %>%  select(starts_with('clage_vacsi')) %>% filter(row_number() == 1) %>% as.numeric()
  tmp_mat <- tmp_df  %>%  select(starts_with('n_dose1')) %>% as.matrix()
  colnames(tmp_mat) <- name_ages
  rownames(tmp_mat) <- (tmp_df  %>%  select(jour))[[1]] %>% as.character()
  
  tmp_mat_2 <- cbind(rep(0, nrow(tmp_mat)),
                     rep(0, nrow(tmp_mat)),
                     tmp_mat[, '24'] + tmp_mat[, '29'],
                     tmp_mat[, '39'],
                     tmp_mat[, '49']/2, tmp_mat[, '49']/2,
                     tmp_mat[, '59']/2, tmp_mat[, '59']/2,
                     tmp_mat[, '64'], tmp_mat[, '69'],
                     tmp_mat[, '74'], tmp_mat[, '79'],
                     tmp_mat[, '80'])
  return(tmp_mat_2)
}

get_Mat_Rates_from_Mat_VaccinationStrategy <- function(RealMatVaccinationStrategy,
                                                       popSize.ageG){
 
  list_Mat_AgeOnly <- list(RealMatVaccinationStrategy)
  names(list_Mat_AgeOnly) <- 'Mat_VaccineStrategy_1'
  nVaccines <- length(list_Mat_AgeOnly)
  nDays <- nrow(list_Mat_AgeOnly[[1]])
  n.ageG <- length(popSize.ageG)
  
  list_Mat_Rates <- lapply(1:nVaccines, FUN = function(i){
    matrix(0, ncol = n.ageG, nrow = nDays)
  })
  names(list_Mat_Rates) <- names(list_Mat_AgeOnly)
  nNonVaccinated.ageG <- popSize.ageG
  for(iDay in 1:nDays){
    list_VectToVaccinate <- lapply(list_Mat_AgeOnly, FUN = function(l){
      l[iDay,]
    })
    list_Vectrates <- lapply(list_VectToVaccinate, FUN = function(l){
      nNewNonVaccinated.ageG <- nNonVaccinated.ageG - l
      log(nNonVaccinated.ageG) - log(nNewNonVaccinated.ageG)
    })
    ## Update the number of non vaccinated individuals in the different age groups
    nNonVaccinated.ageG <- nNonVaccinated.ageG - Reduce('+', list_VectToVaccinate)
    ## Update the list with the rates
    for(iVaccine in 1:nVaccines){
      list_Mat_Rates[[paste0('Mat_VaccineStrategy_', iVaccine)]][iDay, ] <- list_Vectrates[[paste0('Mat_VaccineStrategy_', iVaccine)]] 
    }
  }
  return(list_Mat_Rates)
}

get_Mat_Rates_RetrospectiveAndProspective <- function(my_region_name,
                                                      date_end_vaccine,
                                                      vect_maxVC,
                                                      seq_n_maxDosesPerDay,
                                                      df_Doses,
                                                      seq_DelayBetweenDoses,
                                                      seq_DateChangePace,
                                                      list_vaccine_priority){
  
  if(! (my_region_name %in% c('Metro', vect_name_region))){
    stop('Schedule not implemented for region: ', my_region_name)
  }
  
  ## Get the approximate matrix describing the number of doses received by the different individuals since the beginning of the 
  RealMatVaccinationStrategy <- get_Retrospective_Mat_VaccinationStrategy(my_region_name)
  RealMatRates <- get_Mat_Rates_from_Mat_VaccinationStrategy(RealMatVaccinationStrategy, get_populationVector(my_region_name))
  RealMatVaccinationStrategy_Metro <- get_Retrospective_Mat_VaccinationStrategy('Metro')
  
  DatesRetrospective <- as.Date(rownames(RealMatVaccinationStrategy))
  date_beginning_vaccine <- min(DatesRetrospective)
  date_beginning_projections_vaccination <- max(DatesRetrospective) + 1
  
  PropDosesToRegion <- sum(RealMatVaccinationStrategy)/sum(RealMatVaccinationStrategy_Metro)
  
  ## Correcting the population vector and the vector of maxVC to account for already vaccinated individuals
  popSize.ageG <- get_populationVector(my_region_name) - apply(RealMatVaccinationStrategy, 2, sum)
  maxVaccinated <- get_populationVector(my_region_name)*vect_maxVC
  maxVaccinated_RemoveAlreadyVaccinated <- maxVaccinated - apply(RealMatVaccinationStrategy, 2, sum)
  vect_maxVC_corr <- maxVaccinated_RemoveAlreadyVaccinated/popSize.ageG
  
  ## Correct the dataframe with the number of doses available
  df_Doses_RemoveRetrospective <- df_Doses[df_Doses$DateDelivery >= lubridate::floor_date(date_beginning_projections_vaccination, '%m'), ]
  
  ## Get the number of first doses that were distributed at the national level
  df_Doses_Retrospective <- (apply(df_Doses[df_Doses$DateDelivery < lubridate::floor_date(date_beginning_projections_vaccination, '%m'), 
                                                             grep(colnames(df_Doses), pattern = 'Doses_')],
                                     2, sum))
  DosesDistributedMetropolitanFrance <- sum(RealMatVaccinationStrategy_Metro)*2
  RemainingRetrospective <- df_Doses_Retrospective - df_Doses_Retrospective/sum(df_Doses_Retrospective)*DosesDistributedMetropolitanFrance
  
  df_Doses_RemoveRetrospective[, grep(names(df_Doses_RemoveRetrospective), pattern = 'Doses_')][1,] <- 
    df_Doses_RemoveRetrospective[, grep(names(df_Doses_RemoveRetrospective), pattern = 'Doses_')][1,] +
    RemainingRetrospective
  
  ## Distribute remaining doses across regions
  df_Doses_RemoveRetrospective[, grep(names(df_Doses_RemoveRetrospective), pattern = 'Doses_')] <- 
    PropDosesToRegion*df_Doses_RemoveRetrospective[, grep(names(df_Doses_RemoveRetrospective), pattern = 'Doses_')]
  
  
  ## Determine the change in the vaccination pace to mimic the distribution of the second doses
  LengthTimeWindowDistribution <- 28
  seq_DateChangePace <- rep(date_beginning_projections_vaccination + LengthTimeWindowDistribution - 1, 3)
  NSecondDosesToRemovePerDay <- sum(RealMatVaccinationStrategy[(nrow(RealMatVaccinationStrategy) - LengthTimeWindowDistribution + 1):nrow(RealMatVaccinationStrategy), ])/LengthTimeWindowDistribution

  tmp_seq_n_maxDosesPerDay_Before <- PropDosesToRegion*seq_n_maxDosesPerDay - seq_n_maxDosesPerDay/sum(seq_n_maxDosesPerDay)*NSecondDosesToRemovePerDay
  tmp_seq_n_maxDosesPerDay_After <- PropDosesToRegion*seq_n_maxDosesPerDay
  
  if(sum(tmp_seq_n_maxDosesPerDay_Before <0) >0){
    stop('Negative number of doses distributed in region ', my_region_name)
  }
  
  MatRatesProspective <- get_Mat_Rates_nVaccines_justAge_ChangePace(date_beginning_vaccine = date_beginning_projections_vaccination,
                                                                    date_end_vaccine = date_end_vaccine,
                                                                    vect_maxVC = vect_maxVC_corr,
                                                                    seq_n_maxDosesPerDay_Before = tmp_seq_n_maxDosesPerDay_Before,
                                                                    seq_n_maxDosesPerDay_After = tmp_seq_n_maxDosesPerDay_After,
                                                                    df_Doses = df_Doses_RemoveRetrospective,
                                                                    seq_DelayBetweenDoses = seq_DelayBetweenDoses,
                                                                    seq_DateChangePace = seq_DateChangePace,
                                                                    list_vaccine_priority = list_vaccine_priority,
                                                                    popSize.ageG = popSize.ageG)
  
  
  
  MatRatesAll <- rbind(RealMatRates[[1]], Reduce('+', MatRatesProspective))
  return(MatRatesAll)
}
