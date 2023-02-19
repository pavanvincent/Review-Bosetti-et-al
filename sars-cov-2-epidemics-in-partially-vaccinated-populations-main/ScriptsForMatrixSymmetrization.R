## Getting age distribution of the population for a defined geographical area
## (metropolitan France or specific region)
get_populationVector <- function(name_area){
  # region name is the name of the area of interest
  # National for metropolitan France or the short abbreviations of the regions otherwise
  
  if(name_area == 'Metro'){
    population_data <- read.csv('Data/population_data/PopByRegionAndAge_MetropolianFrance.csv')
    vect_name_region <- c('ARA', 'BFC', 'BRE', 'CVL', 'COR', 'GES', 'HDF', 'IDF', 'NAQ', 'NOR', 'OCC', 'PAC', 'PDL')
    PopGeographicalArea <- apply(population_data[, vect_name_region], 1, sum)
  } else{
    population_data <- read.csv('Data/population_data/PopByRegionAndAge_EntireFrance.csv')
    PopGeographicalArea <- population_data[, name_area]
  }
  names(PopGeographicalArea) <- population_data$Age
  return(PopGeographicalArea)
}


## Function used to adjust the matrices for different infectivity/susceptibility
adjust_matrix_infectivity_susceptibility <- function(contactMatrix, susceptibility, infectivity){
  # contactMatrix is of size nAge*nAge
  # susceptibility and infectivity are vectors of length nAge
  
  contactMatrix_corr<- matrix(0, ncol = ncol(contactMatrix), nrow = nrow(contactMatrix))
  for(i in 1:nrow(contactMatrix_corr)){
    for(j in 1:ncol(contactMatrix_corr)){
      contactMatrix_corr[i,j] <- contactMatrix[i,j]*infectivity[j]*susceptibility[i]
    }
  }
  return(contactMatrix_corr)
}


## Function used to symmetrize the matrices
## Symmetrization well explained in Funk et al., 2019 (BMC Medicine)
symmetrize_contact_matrix <- function(contactMatrix,
                                      popSize.ageG){
  
  n.ageG <- length(popSize.ageG)
  tmp_mat_pop <- matrix(rep(popSize.ageG, n.ageG), ncol = n.ageG, nrow = n.ageG)
  MatPop <- contactMatrix*tmp_mat_pop
  
  NormalizedMat <- (MatPop + t(MatPop))/(2*tmp_mat_pop)
  return(NormalizedMat)
}


## Getting the maximum eigenvalue of a matrix
get_max_eigenval <- function(M){
  eigenvalues <- eigen(M)$values
  max_eigenval <- max(Re(eigenvalues[abs(Im(eigenvalues)) < 1e-6]))
  return(max_eigenval)
}

## Compute the reduced matrix from a matrix of reduction of contacts for all the age groups
compute_matrix_reduction_all_contacts <- function(CM,
                                                  vect_alpha){
  mat_min_alpha <- sapply(vect_alpha, FUN = function(tmp_alpha){
    sapply(vect_alpha, FUN = function(tmp_alpha2){
      min(tmp_alpha, tmp_alpha2)
    })
  })
  mat_res <- CM*mat_min_alpha
  return(mat_res)
}

compute_normalized_matrix_reduction_all_contacts <- function(CM,
                                                             vect_alpha){
  mat_res <- compute_matrix_reduction_all_contacts(CM, vect_alpha)
  mat_res <- mat_res/get_max_eigenval(mat_res)
  return(mat_res)
}

get_vect_alpha_from8 <- function(vect_8 # vector of size 8
                           ){
  return(c(vect_8[1:2], #0-9 ; 10-17
    1, #18-29
    vect_8[3], #30-39
    rep(vect_8[4], 2), #40-49
    rep(vect_8[5], 2), #50-59
    rep(vect_8[6], 2), #60-69
    rep(vect_8[7], 2), #70-79
    vect_8[8] #80p
    ))
}
