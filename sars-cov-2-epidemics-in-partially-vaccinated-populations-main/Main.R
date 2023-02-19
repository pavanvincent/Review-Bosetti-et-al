source("models.R")
source('ScriptsForMatrixSymmetrization.R')
library(odin)

# Load the vector of population size for metropolitan France
chosen_region <- "Metro"
pop_by_region_newAge <- read.csv("Data/population_data/PopByRegionAndAge_MetropolianFrance.csv")
popSize.ageG <- get_populationVector(chosen_region)
nAge <- length(popSize.ageG)

# Load the matrices and adjust for the role of the children
contactMatrix <- readRDS('Data/contact_matrix/mat_regular.rds')
contactMatrix <- symmetrize_contact_matrix(contactMatrix, popSize.ageG)
contactMatrix <- adjust_matrix_infectivity_susceptibility(contactMatrix, susceptibility = c(0.5, 0.75, rep(1, nAge - 2)), infectivity = rep(1, nAge))
c_mat1_regular <- contactMatrix/get_max_eigenval(contactMatrix)

# Proportion of 12-17y.o. among 10-17y.o.
prop_12_17_to_10_17 <- 0.74823467

# Load data on severity
df_severity <- read.csv2('Data/severity_estimates/SeverityEstimatesOurAgeGroups.csv')
p.hosp.inf <- df_severity$IHR*0.85

# Load initial conditions from calibration at September 1st
init_cond_Sept1st <- readRDS("Data/Initial_Conditions_Sept1st.rds")
cnames <- names(init_cond_Sept1st)

# Extract values from  compartments of unvaccinated individuals for the two variants
S_ini_tot <- init_cond_Sept1st[grep("S\\[", cnames)]
E1_ini_1 <- init_cond_Sept1st[grep("E1\\[1", cnames)]
E1_ini_2 <- init_cond_Sept1st[grep("E1\\[2", cnames)]
E2_ini_1 <- init_cond_Sept1st[grep("E2\\[1", cnames)]
E2_ini_2 <- init_cond_Sept1st[grep("E2\\[2", cnames)]
IMild_ini_1 <- init_cond_Sept1st[grep("IMild\\[1", cnames)]
IMild_ini_2 <- init_cond_Sept1st[grep("IMild\\[2", cnames)]
R_ini_1 <- init_cond_Sept1st[grep("R\\[1", cnames)]
R_ini_2 <- init_cond_Sept1st[grep("R\\[2", cnames)]
IHosp_ini_1 <- init_cond_Sept1st[grep("IHosp\\[1", cnames)]
IHosp_ini_2 <- init_cond_Sept1st[grep("IHosp\\[2", cnames)]
IBarHosp_ini_1 <- init_cond_Sept1st[grep("IBarHosp\\[1", cnames)]
IBarHosp_ini_2 <- init_cond_Sept1st[grep("IBarHosp\\[2", cnames)]
H_ini_1 <- init_cond_Sept1st[grep("H\\[1", cnames)]
H_ini_2 <- init_cond_Sept1st[grep("H\\[2", cnames)]

#Extract values from  compartments  of vaccinated individuals for the two variants
Sv_ini_tot <- init_cond_Sept1st[grep("Sv\\[", cnames)]
E1v_ini_1 <- init_cond_Sept1st[grep("E1v\\[1", cnames)]
E1v_ini_2 <- init_cond_Sept1st[grep("E1v\\[2", cnames)]
E2v_ini_1 <- init_cond_Sept1st[grep("E2v\\[1", cnames)]
E2v_ini_2 <- init_cond_Sept1st[grep("E2v\\[2", cnames)]
IMildv_ini_1 <- init_cond_Sept1st[grep("IMildv\\[1", cnames)]
IMildv_ini_2 <- init_cond_Sept1st[grep("IMildv\\[2", cnames)]
Rv_ini_1 <- init_cond_Sept1st[grep("Rv\\[1", cnames)]
Rv_ini_2 <- init_cond_Sept1st[grep("Rv\\[2", cnames)]
IHospv_ini_1 <- init_cond_Sept1st[grep("IHospv\\[1", cnames)]
IHospv_ini_2 <- init_cond_Sept1st[grep("IHospv\\[2", cnames)]
IBarHospv_ini_1 <- init_cond_Sept1st[grep("IBarHospv\\[1", cnames)]
IBarHospv_ini_2 <- init_cond_Sept1st[grep("IBarHospv\\[2", cnames)]
Hv_ini_1 <- init_cond_Sept1st[grep("Hv\\[1", cnames)]
Hv_ini_2 <- init_cond_Sept1st[grep("Hv\\[2", cnames)]

#Merge the two variants 
E1_ini_tot <- E1_ini_1 + E1_ini_2
E2_ini_tot <- E2_ini_1 + E2_ini_2
IMild_ini_tot <- IMild_ini_1 + IMild_ini_2
R_ini_tot <- R_ini_1 + R_ini_2
IHosp_ini_tot <- IHosp_ini_1 + IHosp_ini_2
IBarHosp_ini_tot <- IBarHosp_ini_1 + IBarHosp_ini_2
H_ini_tot <- H_ini_1 + H_ini_2
ini_tot_no_v <- S_ini_tot + E1_ini_tot + E2_ini_tot + IMild_ini_tot + R_ini_tot + IHosp_ini_tot + IBarHosp_ini_tot + H_ini_tot 

E1v_ini_tot <- E1v_ini_1 + E1v_ini_2
E2v_ini_tot <- E2v_ini_1 + E2v_ini_2
IMildv_ini_tot <- IMildv_ini_1 + IMildv_ini_2
Rv_ini_tot <- Rv_ini_1 + Rv_ini_2
IHospv_ini_tot <- IHospv_ini_1 + IHospv_ini_2
IBarHospv_ini_tot <- IBarHospv_ini_1 + IBarHospv_ini_2
Hv_ini_tot <- Hv_ini_1 + Hv_ini_2
ini_tot_v <-Sv_ini_tot+ E1v_ini_tot + E2v_ini_tot + IMildv_ini_tot + Rv_ini_tot + IHospv_ini_tot + IBarHospv_ini_tot + Hv_ini_tot 

# Define increase transmissibility of one variant with respect to the other (in this case no increase in transmissibility)
alpha_variants <- c(1.0, 1.0) 

# Change points for R0 an contact matrix 
Nday <- 200
interp_ts_R0 <- c(0, Nday)
interp_ts_c_mat <- c(0, Nday)
interp_c_mats <- array(NA, dim = c(2, length(popSize.ageG), length(popSize.ageG)))
interp_c_mats[1,,] <- c_mat1_regular
interp_c_mats[2,,] <- c_mat1_regular

# Define matrix for the vaccination rates 
#(In this case is set to 0 because vaccination is stopped during the simulation period)
interp_ts_v_rate <- 0:Nday
interp_v_rate <-array(0, dim = c(length(interp_ts_v_rate), length(popSize.ageG), 1))
 
#Increase in the probability of hospitalization upon infection due to Delta VOC.
#We assume that Delta VOC is 50% more severe than Alpha VOC (Twohig et al 2021), 
#while Alpha VOC is 40% more severe than previously circulating strains (Bager et al 2021)
pHosp_mat <- rbind(1.4*1.5*p.hosp.inf, 1.4*1.5*p.hosp.inf) 

# Vaccine scenarios
vaccine_scenarios <- expand.grid(coverage_60_plus = c(0.90, 0.95),
                                 coverage_18_60  = c(0.60,0.80,0.90),
                                 coverage_12_18 = c(0.0,0.3,0.7))
vaccine_scenarios <- vaccine_scenarios[vaccine_scenarios$coverage_18_60 >= vaccine_scenarios$coverage_12_18,]
vaccine_scenarios <- vaccine_scenarios[vaccine_scenarios$coverage_60_plus >= vaccine_scenarios$coverage_18_60,]
rownames(vaccine_scenarios) <- NULL

# Define parameters to run the simulation
index_scenario_vacc <- 13 #Row nummber in the vaccine_scenarios matrix (i.e. 13 for Baseline)

R0 <- 5.0 # R0 
reduction_R0 <- 0.0 # Reduction in the transmission rate
reduction_onlyunvacc <- 0 # Reduction in the transmission rate applied to all (0) or to unvaccinated (1)
prop_I <- 0.25 # Proportion of infected individuals at September 1st 

prop_tested <- 0.0 # Proportion of people tested 
delay_betw_tests <- 7 # Delay between tests (e.g. 7 days)
tests_random <- 0 # Testing unvaccinated individuals (0) or randmly (1)  
sensitivity <- 0.75 # Sensitivity of the test (0.75 for autotests and 0.9 for antigenic tests)
beta_reduction <- (1 - 0.75) # Reduction in transmission after testing



# Define Vaccine efficacy against infection (VE_susc), transmission (VE_inf), hospitalisation (VE_hosp).
# VE_sev is parametrized as in Tran Kiem et al 2021 (EclinicalMedicine) from VE_hosp to account for 
# the Vaccine efficacy against infection (VE_susc).
# rows = efficacy to variants 1 or 2, columns = vaccines 1 or 2 (in this case a single vaccine is considered)
# e.g. position [1,2] efficacy of vaccines 2 to variant 1 
VE_susc <- cbind(c(0.6,0.6))
VE_inf <- cbind(c(0.5,0.5))
VE_hosp <- c(0.95, 0.95)
VE_sev <- 1-((1-VE_hosp)/(1-VE_susc))

# Define the testing rates (see Additional information).  
if (tests_random == 0){
  testing_rate <- c(0.0,
                    prop_12_17_to_10_17*sensitivity*prop_tested*(1.0/delay_betw_tests),
                    rep(sensitivity*prop_tested*(1.0/delay_betw_tests),11))
  testing_rate_v <- rep(0.0,13)
}else if(tests_random == 1){
  
  testing_rate <- c(0.0,
                    prop_12_17_to_10_17*sensitivity*prop_tested*(1 - vaccine_scenarios[index_scenario_vacc,"coverage_12_18"])*(1.0/delay_betw_tests),
                    rep(sensitivity*prop_tested*(1 - vaccine_scenarios[index_scenario_vacc,"coverage_18_60"])*(1.0/delay_betw_tests),6),
                    rep(sensitivity*prop_tested*(1 - vaccine_scenarios[index_scenario_vacc,"coverage_60_plus"])*(1.0/delay_betw_tests),5))
  testing_rate_v <- c(0.0,
                      prop_12_17_to_10_17*sensitivity*prop_tested*(1 - vaccine_scenarios[index_scenario_vacc,"coverage_12_18"])*(1.0/delay_betw_tests),
                      rep(sensitivity*prop_tested*(1 - vaccine_scenarios[index_scenario_vacc,"coverage_18_60"])*(1.0/delay_betw_tests),6),
                      rep(sensitivity*prop_tested*(1 - vaccine_scenarios[index_scenario_vacc,"coverage_60_plus"])*(1.0/delay_betw_tests),5))
}else{print("test_random take values 0 or 1")}

# Reduction in transmissibility applied to all/unvaccinated 
if (reduction_onlyunvacc == 0) { # All
  red_unvacc <- 0.0 
  interp_R0s <- c(R0*(1-reduction_R0),R0*(1-reduction_R0))
}else if(reduction_onlyunvacc == 1){ # Unvaccinated
  red_unvacc <- reduction_R0
  interp_R0s <- c(R0, R0)
}else{print("reduction_onlyunvacc take values 0 or 1")}

# Recompute immunity acquired through natural infection on September 1st (see dditional information)
tmp_tomove <- (sum(S_ini_tot + Sv_ini_tot)/sum(popSize.ageG) - (1-prop_I))*sum(popSize.ageG)
tmp_S_ini <- ((Rv_ini_tot + R_ini_tot)/sum(Rv_ini_tot + R_ini_tot))*tmp_tomove*(1-ini_tot_v/popSize.ageG)
tmp_Sv_ini <- ((Rv_ini_tot + R_ini_tot)/sum(Rv_ini_tot + R_ini_tot))*tmp_tomove*(ini_tot_v/popSize.ageG)
S_ini_tot_addInf <- S_ini_tot - tmp_S_ini
R_ini_tot_addInf <- R_ini_tot + tmp_S_ini
Sv_ini_tot_addInf <- Sv_ini_tot - tmp_Sv_ini 
Rv_ini_tot_addInf <- Rv_ini_tot + tmp_Sv_ini

#Recompute vaccination coverage on September 1st according to the selected Scenario  
vaccination_start <- c(0.0,
                       prop_12_17_to_10_17*vaccine_scenarios[index_scenario_vacc,"coverage_12_18"],
                       rep(vaccine_scenarios[index_scenario_vacc,"coverage_18_60"],6),
                       rep(vaccine_scenarios[index_scenario_vacc,"coverage_60_plus"],5))

add_vaccinate <- (vaccination_start -(1 - (ini_tot_no_v/popSize.ageG)) ) * (popSize.ageG/ini_tot_no_v )
add_vaccinate[add_vaccinate <0] <- 0

S_ini_tot_new <- S_ini_tot_addInf - S_ini_tot_addInf*add_vaccinate
E1_ini_tot_new <-  E1_ini_tot - E1_ini_tot*add_vaccinate
E2_ini_tot_new <- E2_ini_tot - E2_ini_tot*add_vaccinate
IMild_ini_tot_new <- IMild_ini_tot - IMild_ini_tot*add_vaccinate 
R_ini_tot_new <- R_ini_tot_addInf - R_ini_tot_addInf*add_vaccinate
IHosp_ini_tot_new <- IHosp_ini_tot - IHosp_ini_tot*add_vaccinate
IBarHosp_ini_tot_new <-  IBarHosp_ini_tot - IBarHosp_ini_tot*add_vaccinate
H_ini_tot_new <- H_ini_tot - H_ini_tot*add_vaccinate

Sv_ini_tot_new <- Sv_ini_tot_addInf + S_ini_tot_addInf*add_vaccinate
E1v_ini_tot_new <-  E1v_ini_tot + E1_ini_tot*add_vaccinate
E2v_ini_tot_new <- E2v_ini_tot + E2_ini_tot*add_vaccinate
IMildv_ini_tot_new <- IMildv_ini_tot + IMild_ini_tot*add_vaccinate 
Rv_ini_tot_new <- Rv_ini_tot_addInf + R_ini_tot_addInf*add_vaccinate
IHospv_ini_tot_new <- IHospv_ini_tot + IHosp_ini_tot*add_vaccinate
IBarHospv_ini_tot_new <-  IBarHospv_ini_tot + IBarHosp_ini_tot*add_vaccinate
Hv_ini_tot_new <- Hv_ini_tot + H_ini_tot*add_vaccinate
  
#Initial conditions
#Unvaccinated compartments
start_perc_newV <- (1.0-0.05)
S_ini <- as.double(S_ini_tot_new) 
E1_ini <- rbind(as.double(E1_ini_tot_new)*start_perc_newV, as.double(E1_ini_tot_new)*(1.0 - start_perc_newV))
E2_ini <- rbind(as.double(E2_ini_tot_new)*start_perc_newV, as.double(E2_ini_tot_new)*(1.0 - start_perc_newV))
IMild_ini <- rbind(as.double(IMild_ini_tot_new)*start_perc_newV, as.double(IMild_ini_tot_new)*(1.0 - start_perc_newV))
IHosp_ini <- rbind(as.double(IHosp_ini_tot_new)*start_perc_newV, as.double(IHosp_ini_tot_new)*(1.0 - start_perc_newV))
R_ini <-  rbind(as.double(R_ini_tot_new + H_ini_tot_new),rep(0, length(popSize.ageG)))
H_ini <-  rbind(rep(0, length(popSize.ageG)) ,rep(0, length(popSize.ageG)))
IBarHosp_ini <- rbind(as.double(IBarHosp_ini_tot_new),rep(0, length(popSize.ageG)))
E2_iso_ini <- IMild_iso_ini <- IHosp_iso_ini <- rbind(rep(0, length(popSize.ageG)), rep(0, length(popSize.ageG)))

#Vaccinated compartments
Sv_ini <-  cbind( as.double(Sv_ini_tot_new))
E1v_ini <- array(0, dim = c(2,  length(popSize.ageG), 1))
E2v_ini <- array(0, dim = c(2,  length(popSize.ageG), 1))
IMildv_ini <-  array(0, dim = c(2,  length(popSize.ageG), 1))
IHospv_ini <-  array(0, dim = c(2,  length(popSize.ageG), 1))
Rv_ini <- array(0, dim = c(2,  length(popSize.ageG), 1))
Hv_ini <-  array(0, dim = c(2,  length(popSize.ageG), 1))
IBarHospv_ini <-  array(0, dim = c(2,  length(popSize.ageG), 1))

E1v_ini[1,,] <- cbind(as.double(E1v_ini_tot_new)*start_perc_newV)
E1v_ini[2,,] <- cbind(as.double(E1v_ini_tot_new)*(1.0 - start_perc_newV))
E2v_ini[1,,] <- cbind(as.double(E2v_ini_tot_new)*start_perc_newV)
E2v_ini[2,,] <- cbind(as.double(E2v_ini_tot_new)*(1.0 - start_perc_newV))
IMildv_ini[1,,] <- cbind(as.double(IMildv_ini_tot_new)*start_perc_newV)
IMildv_ini[2,,] <- cbind(as.double(IMildv_ini_tot_new)*(1.0 - start_perc_newV))
IHospv_ini[1,,] <- cbind(as.double(IHospv_ini_tot_new)*start_perc_newV)
IHospv_ini[2,,] <- cbind(as.double(IHospv_ini_tot_new)*(1.0 - start_perc_newV))
IBarHospv_ini[1,,] <- cbind(as.double(IBarHospv_ini_tot_new))
Rv_ini[1,,] <- cbind(as.double( Rv_ini_tot_new + Hv_ini_tot_new))
E2v_iso_ini <- IMildv_iso_ini <- IHospv_iso_ini <-  array(0, dim = c(2,  length(popSize.ageG), 1))

#Compile the model  
mod_var <- seeir_age_variants_nvacc_testing$new(pop = popSize.ageG,
                                            pHosp = pHosp_mat,
                                            alpha_variants = alpha_variants,
                                            beta_reduction = beta_reduction,
                                            red_unvacc = red_unvacc,
                                            testing_rate = testing_rate,
                                            testing_rate_v = testing_rate_v,
                                            interp_ts_R0 = interp_ts_R0,
                                            interp_R0s = interp_R0s,
                                            interp_ts_c_mat = interp_ts_c_mat,
                                            interp_c_mats = interp_c_mats,
                                            interp_ts_v_rate = interp_ts_v_rate,
                                            interp_v_rate = interp_v_rate,
                                            S_ini = S_ini,
                                            E1_ini = E1_ini,
                                            E2_ini = E2_ini,
                                            E2_iso_ini = E2_iso_ini,
                                            IMild_ini = IMild_ini,
                                            IMild_iso_ini = IMild_iso_ini,
                                            R_ini = R_ini,
                                            IHosp_ini = IHosp_ini,
                                            IHosp_iso_ini = IHosp_iso_ini,
                                            IBarHosp_ini = IBarHosp_ini,
                                            H_ini = H_ini,
                                            Sv_ini = Sv_ini,
                                            E1v_ini = E1v_ini,
                                            E2v_ini = E2v_ini,
                                            E2v_iso_ini = E2v_iso_ini,
                                            IMildv_ini = IMildv_ini,
                                            IMildv_iso_ini = IMildv_iso_ini,
                                            Rv_ini = Rv_ini,
                                            IHospv_ini = IHospv_ini,
                                            IHospv_iso_ini = IHospv_iso_ini,
                                            IBarHospv_ini = IBarHospv_ini,
                                            Hv_ini = Hv_ini,
                                            VE_susc = VE_susc,
                                            VE_inf = VE_inf,
                                            VE_sev = VE_sev)

# Solve ODE system
ts <- seq(0, Nday, length.out = Nday+1)
y <- mod_var$run(ts)


# RESULTS

# Number of hospital admissions at peak
max(y[,"iHosp"])
  
# Function to merge the 13 age groups into 3 age groups (0-17,18-59,60+)  
group_age <- function(v){
  age_groups <- rep(c(1,2,3),c(2,6,5))
  return( round(unname(tapply(as.double(v), age_groups, sum)),6))
}

# Number of unvaccinated/vaccinated for each age group
cnames_y <- colnames(y)
nb_notvacc <- group_age(popSize.ageG - y[1,grep("Vacc_age\\[", cnames_y)])
nb_vacc <- group_age(y[1,grep("Vacc_age\\[", cnames_y)])
# Number of hospitalisations for unvaccinated/vaccinated for each age group
nb_hosp_vacc <- group_age(colSums( y[,grep("iHosp_vacc_age\\[", cnames_y)]))
nb_hosp_notvacc <- group_age(colSums( y[,grep("iHosp_notvacc_age\\[", cnames_y)]))
# Number of infections for unvaccinated/vaccinated for each age group
nb_iInf_vacc <- group_age(colSums( y[,grep("iInf_vacc_age\\[", cnames_y)]))
nb_iInf_notvacc <- group_age(colSums( y[,grep("iInf_notvacc_age\\[", cnames_y)]))
# Number of infections due to unvaccinated/vaccinated for each age group
nb_iInf_byvacc <- group_age(colSums( y[,grep("iInf_byvacc_age\\[", cnames_y)]))
nb_iInf_bynotvacc <- group_age(colSums( y[,grep("iInf_bynotvacc_age\\[", cnames_y)]))

# Risk Ratios
RR_hospitalized_age <- (nb_hosp_notvacc/nb_notvacc)/(nb_hosp_vacc/nb_vacc)
RR_hospitalized_all <-(sum(nb_hosp_notvacc)/sum(nb_notvacc))/(sum(nb_hosp_vacc)/sum(nb_vacc))

RR_iInf_age <- (nb_iInf_notvacc/nb_notvacc)/(nb_iInf_vacc /nb_vacc)
RR_iInf_all <- (sum(nb_iInf_notvacc)/sum(nb_notvacc))/(sum(nb_iInf_vacc)/sum(nb_vacc))

RR_trasm_age <- (nb_iInf_bynotvacc/nb_notvacc)/(nb_iInf_byvacc/nb_vacc)
RR_trasm_all <- (sum(nb_iInf_bynotvacc)/sum(nb_notvacc))/(sum(nb_iInf_byvacc)/sum(nb_vacc))


# Contribution of groups defined by their age and vaccination status to infections
((nb_iInf_vacc + nb_iInf_notvacc)/sum(nb_iInf_vacc + nb_iInf_notvacc))*100 # All
(nb_iInf_vacc /sum(nb_iInf_vacc + nb_iInf_notvacc))*100 # Vaccinated individuals
(nb_iInf_notvacc/sum(nb_iInf_vacc + nb_iInf_notvacc))*100 # Unvaccinated individuals

# Contribution of groups defined by their age and vaccination status to disease spread 
((nb_iInf_byvacc + nb_iInf_bynotvacc)/sum(nb_iInf_vacc + nb_iInf_notvacc))*100 # All
(nb_iInf_byvacc /sum(nb_iInf_vacc + nb_iInf_notvacc))*100 # Vaccinated individuals
(nb_iInf_bynotvacc/sum(nb_iInf_vacc + nb_iInf_notvacc))*100 # Unvaccinated individuals

# Contribution of groups defined by their age and vaccination status to hospital burden
((nb_hosp_vacc + nb_hosp_notvacc)/sum(nb_hosp_vacc + nb_hosp_notvacc))*100 # All
(nb_hosp_vacc/sum(nb_hosp_vacc + nb_hosp_notvacc))*100 # Vaccinated individuals
(nb_hosp_notvacc/sum(nb_hosp_vacc + nb_hosp_notvacc))*100 # Unvaccinated individuals

# Age distribution of the different groups in the population
population_age <- group_age(popSize.ageG)
(population_age/sum(population_age))*100 # All
(nb_vacc/sum(population_age))*100# Vaccinated individuals
(nb_notvacc/sum(population_age))*100 # Unvaccinated individuals
  
