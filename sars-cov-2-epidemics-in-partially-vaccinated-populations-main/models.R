require(odin)


################################################################################
# Age-structured model with 2 variants and n vaccines and testing
################################################################################
seeir_age_variants_nvacc_testing <- odin({
  # NOTE: First index = variant, second (and third for matrices) = age, third = vaccines
  R0 <- interpolate(interp_ts_R0, interp_R0s, "constant")
  c_mat[, ] <- interpolate(interp_ts_c_mat, interp_c_mats, "constant")
  v_rate[,] <- interpolate(interp_ts_v_rate, interp_v_rate, "constant")
  
  D <- 1. / g2 + 1. / g3 # Infectious period
  beta[] <- alpha_variants[i] * R0 / D # Transmission rate
  prop_Iv[,,] <-  (1-VE_inf[i,k])*(E2v[i,j,k] + IMildv[i,j,k] + IHospv[i,j,k])
  prop_Iv_iso[,,] <-  (1-VE_inf[i,k])*(beta_reduction)*(E2v_iso[i,j,k] + IMildv_iso[i,j,k] + IHospv_iso[i,j,k])
  prop_I[,] <- ( (1.0 - red_unvacc )*(E2[i, j] + IMild[i,j] + IHosp[i,j]) +
                   (1.0 - red_unvacc )*beta_reduction*(E2_iso[i, j] + IMild_iso[i,j] + IHosp_iso[i,j]) +
                   sum(prop_Iv[i,j,]) +  sum(prop_Iv_iso[i,j,]) ) / pop[j] # Proportion infectious
  c_prop_I[,,] <- c_mat[j, k] * prop_I[i,k] # Contact matrix * prop infectious
  
  
  prop_I_notvacc[,]  <- ((1.0 - red_unvacc )*(E2[i, j] + IMild[i,j] + IHosp[i,j]) +
                           (1.0 - red_unvacc )*beta_reduction*(E2_iso[i, j] + IMild_iso[i,j] + IHosp_iso[i,j])) / pop[j]
  c_prop_I_notvacc[,,] <- c_mat[j, k] * prop_I_notvacc[i,k]
  
  
  iInf_bynotvacc_tmp[,,] <- beta[i] * S[j] * c_prop_I_notvacc[i, j,k] 
  iInf_bynotvaccv_tmp[,,,] <- (1 - VE_susc[i,k])* beta[i] * Sv[j,k] * c_prop_I_notvacc[i, j, l] 
  
  prop_I_vacc[,] <- (sum(prop_Iv[i,j,]) +  sum(prop_Iv_iso[i,j,])) / pop[j] 
  c_prop_I_vacc[,,] <- c_mat[j, k] * prop_I_vacc[i,k]
  
  iInf_byvacc_tmp[,,] <- beta[i] * S[j] * c_prop_I_vacc[i, j,k] 
  iInf_byvaccv_tmp[,,,] <- (1 - VE_susc[i,k])* beta[i] * Sv[j,k] * c_prop_I_vacc[i, j, l] 
  
  
  ##############################################################################
  # Transitions
  ##############################################################################
  # Infections
  S_to_E1[,] <- beta[i] * S[j] * sum(c_prop_I[i, j,])
  E1_to_E2[,] <- g1 * E1[i, j]
  E2_to_I[,] <- g2 * E2[i, j]
  E2_to_E2_iso[,] <- testing_rate[j] * E2[i, j]
  E2_iso_to_I_iso[,] <- g2 * E2_iso[i, j]
  # Mild cases
  E2_to_IMild[,] <- (1. - pHosp[i,j]) * E2_to_I[i, j]
  E2_iso_to_IMild_iso[,] <- (1. - pHosp[i,j]) * E2_iso_to_I_iso[i, j]
  IMild_to_R[,] <- g3 * IMild[i,j]
  IMild_to_IMild_iso[,] <- testing_rate[j] * IMild[i, j]
  IMild_iso_to_R[,] <- g3 * IMild_iso[i,j]
  # Hospitalized cases
  E2_to_IHosp[,] <- pHosp[i,j] * E2_to_I[i,j]
  E2_iso_to_IHosp_iso[,] <- pHosp[i,j] * E2_iso_to_I_iso[i, j]
  IHosp_to_IBarHosp[,] <- g3 * IHosp[i,j]
  IHosp_to_IHosp_iso[,] <- testing_rate[j] * IHosp[i, j]
  IHosp_iso_to_IBarHosp[,] <- g3 * IHosp_iso[i,j]
  IBarHosp_to_H[,] <- gToHosp * IBarHosp[i,j]
  H_to_R[,] <- gHtoR * H[i,j]
  
  
  #vaccination
  S_to_Sv[,] <- v_rate[i,j]*S[i]
  E1_to_E1v[,,] <- v_rate[j,k]*E1[i,j]
  E2_to_E2v[,,] <- v_rate[j,k]*E2[i,j]
  E2_iso_to_E2v_iso[,,] <- v_rate[j,k]*E2_iso[i,j]
  IMild_to_IMildv[,,] <- v_rate[j,k]*IMild[i,j]
  IMild_iso_to_IMildv_iso[,,] <- v_rate[j,k]*IMild_iso[i,j]
  IHosp_to_IHospv[,,] <- v_rate[j,k]*IHosp[i,j]
  IHosp_iso_to_IHospv_iso[,,] <- v_rate[j,k]*IHosp_iso[i,j]
  IBarHosp_to_IBarHospv[,,] <- v_rate[j,k]*IBarHosp[i,j]
  H_to_Hv[,,] <- v_rate[j,k]*H[i,j]
  R_to_Rv[,,] <- v_rate[j,k]*R[i,j]
  
  
  # Infections vaccination 
  Sv_to_E1v[,,] <- (1 - VE_susc[i,k])* beta[i] * Sv[j,k] * sum(c_prop_I[i, j, ])
  E1v_to_E2v[,,] <- g1 * E1v[i,j,k]
  E2v_to_Iv[,,] <- g2 * E2v[i,j,k]
  E2v_to_E2v_iso[,,] <- testing_rate_v[j] * E2v[i,j,k]
  E2v_iso_to_Iv_iso[,,] <- g2 * E2v_iso[i,j,k]
  # Mild cases vaccination
  E2v_to_IMildv[,,] <- (1. - pHosp[i, j]*(1 - VE_sev[i,k])) * E2v_to_Iv[i,j,k]
  E2v_iso_to_IMildv_iso[,,] <- (1. - pHosp[i, j]*(1 - VE_sev[i,k]))* E2v_iso_to_Iv_iso[i,j,k]
  IMildv_to_Rv[,,] <- g3 * IMildv[i,j,k]
  IMildv_to_IMildv_iso[,,] <- testing_rate_v[j] * IMildv[i, j,k]
  IMildv_iso_to_Rv[,,] <- g3 * IMildv_iso[i,j,k]
  # Hospitalized cases vaccination
  E2v_to_IHospv[,,] <- pHosp[i, j]*(1 - VE_sev[i,k]) * E2v_to_Iv[i,j,k]
  E2v_iso_to_IHospv_iso[,,] <-  pHosp[i, j]*(1 - VE_sev[i,k])* E2v_iso_to_Iv_iso[i,j,k]
  IHospv_to_IBarHospv[,,] <- g3 * IHospv[i,j,k]
  IHospv_to_IHospv_iso[,,] <- testing_rate_v[j] * IHospv[i, j,k]
  IHospv_iso_to_IBarHospv[,,] <- g3 * IHospv_iso[i,j,k]
  IBarHospv_to_Hv[,,] <- gToHosp * IBarHospv[i,j,k]
  Hv_to_Rv[,,] <- gHtoR * Hv[i,j,k]
  
  
  
  ##############################################################################
  # Derivatives
  ##############################################################################
  deriv(S[]) <- -sum(S_to_E1[,i]) - sum(S_to_Sv[i,])
  deriv(E1[,]) <- S_to_E1[i,j] - E1_to_E2[i, j] - sum(E1_to_E1v[i,j,])
  deriv(E2[,]) <- E1_to_E2[i,j] - E2_to_I[i, j] - E2_to_E2_iso[i,j] -sum(E2_to_E2v[i,j,])
  deriv(E2_iso[,]) <- E2_to_E2_iso[i,j] - E2_iso_to_I_iso[i, j] -sum(E2_iso_to_E2v_iso[i,j,])
  deriv(IMild[,]) <- E2_to_IMild[i,j] - IMild_to_IMild_iso[i,j]- IMild_to_R[i,j] - sum(IMild_to_IMildv[i,j,])
  deriv(IMild_iso[,]) <- E2_iso_to_IMild_iso[i,j] + IMild_to_IMild_iso[i,j] - IMild_iso_to_R[i,j] - sum(IMild_iso_to_IMildv_iso[i,j,])
  deriv(R[,]) <- IMild_to_R[i,j] + IMild_iso_to_R[i,j] + H_to_R[i,j]  - sum(R_to_Rv[i,j,])
  deriv(IHosp[,]) <- E2_to_IHosp[i, j] - IHosp_to_IHosp_iso[i, j] - IHosp_to_IBarHosp[i, j] - sum(IHosp_to_IHospv[i,j,])
  deriv(IHosp_iso[,]) <- E2_iso_to_IHosp_iso[i, j] + IHosp_to_IHosp_iso[i, j] - IHosp_iso_to_IBarHosp[i, j] - sum(IHosp_iso_to_IHospv_iso[i,j,])
  deriv(IBarHosp[,]) <- IHosp_to_IBarHosp[i,j] + IHosp_iso_to_IBarHosp[i, j] - IBarHosp_to_H[i,j] -sum(IBarHosp_to_IBarHospv[i,j,])
  deriv(H[,]) <- IBarHosp_to_H[i,j] - H_to_R[i,j]  - sum(H_to_Hv[i,j,])
  
  
  deriv(Sv[,]) <- -sum(Sv_to_E1v[,i,j]) + S_to_Sv[i,j]
  deriv(E1v[,,]) <- Sv_to_E1v[i,j,k] - E1v_to_E2v[i,j,k] + E1_to_E1v[i,j,k]
  deriv(E2v[,,]) <- E1v_to_E2v[i,j,k] - E2v_to_Iv[i,j,k] - E2v_to_E2v_iso[i,j,k] +E2_to_E2v[i,j,k] 
  deriv(E2v_iso[,,]) <- E2v_to_E2v_iso[i,j,k] - E2v_iso_to_Iv_iso[i,j,k] +E2_iso_to_E2v_iso[i,j,k] 
  deriv(IMildv[,,]) <- E2v_to_IMildv[i,j,k] - IMildv_to_IMildv_iso[i,j,k] - IMildv_to_Rv[i,j,k] + IMild_to_IMildv[i,j,k] 
  deriv(IMildv_iso[,,]) <- E2v_iso_to_IMildv_iso[i,j,k] + IMildv_to_IMildv_iso[i,j,k] - IMildv_iso_to_Rv[i,j,k] + IMild_iso_to_IMildv_iso[i,j,k]
  deriv(Rv[,,]) <- IMildv_to_Rv[i,j,k]  + IMildv_iso_to_Rv[i,j,k] + Hv_to_Rv[i,j,k] + R_to_Rv[i,j,k]
  deriv(IHospv[,,]) <- E2v_to_IHospv[i,j,k] - IHospv_to_IHospv_iso[i,j,k] - IHospv_to_IBarHospv[i,j,k] + IHosp_to_IHospv[i,j,k] 
  deriv(IHospv_iso[,,]) <- E2v_iso_to_IHospv_iso[i, j,k] + IHospv_to_IHospv_iso[i,j,k] - IHospv_iso_to_IBarHospv[i,j,k] + IHosp_iso_to_IHospv_iso[i,j,k]
  deriv(IBarHospv[,,]) <- IHospv_to_IBarHospv[i,j,k] + IHospv_iso_to_IBarHospv[i,j,k] - IBarHospv_to_Hv[i,j,k] +IBarHosp_to_IBarHospv[i,j,k]
  deriv(Hv[,,]) <- IBarHospv_to_Hv[i,j,k] - Hv_to_Rv[i,j,k] + H_to_Hv[i,j,k]
  
  
  ##############################################################################
  # Initial conditions
  ##############################################################################
  initial(S[]) <- S_ini[i]
  initial(E1[,]) <- E1_ini[i,j]
  initial(E2[,]) <- E2_ini[i,j]
  initial(E2_iso[,]) <- E2_iso_ini[i,j]
  initial(IMild[,]) <- IMild_ini[i,j]
  initial(IMild_iso[,]) <- IMild_iso_ini[i,j]
  initial(R[,]) <- R_ini[i,j]
  initial(IHosp[,]) <- IHosp_ini[i,j]
  initial(IHosp_iso[,]) <- IHosp_iso_ini[i,j]
  initial(IBarHosp[,]) <- IBarHosp_ini[i,j]
  initial(H[,]) <- H_ini[i,j]
  
  initial(Sv[,]) <- Sv_ini[i,j]
  initial(E1v[,,]) <- E1v_ini[i,j,k]
  initial(E2v[,,]) <- E2v_ini[i,j,k]
  initial(E2v_iso[,,]) <- E2v_iso_ini[i,j,k]
  initial(IMildv[,,]) <- IMildv_ini[i,j,k]
  initial(IMildv_iso[,,]) <- IMildv_iso_ini[i,j,k]
  initial(Rv[,,]) <- Rv_ini[i,j,k]
  initial(IHospv[,,]) <- IHospv_ini[i,j,k]
  initial(IHospv_iso[,,]) <- IHospv_iso_ini[i,j,k]
  initial(IBarHospv[,,]) <- IBarHospv_ini[i,j,k]
  initial(Hv[,,]) <- Hv_ini[i,j,k]
  
  
  ##############################################################################
  # Variables to add to output 
  ##############################################################################
  output(iInf_v[]) <- sum(S_to_E1[i, ]) + sum(Sv_to_E1v[i,,]) # + sum(R_to_E1_sec[i,]) + sum(Rv_to_E1v_sec[i,,]) # Incident n. of infections per variant
  output(iHosp_v[]) <- sum(IBarHosp_to_H[i,]) + sum(IBarHospv_to_Hv[i,,] ) # Incident n. of hospitalizations per variant
  output(iInf) <- sum(S_to_E1[,]) + sum(Sv_to_E1v[,,]) # Incident n. of infections
  output(iHosp) <- sum(IBarHosp_to_H[, ]) + sum(IBarHospv_to_Hv[,,]) # Incident n. of hospitalizations
  output(iHosp_age[]) <- sum(IBarHosp_to_H[,i]) + sum(IBarHospv_to_Hv[,i,])
  
  output(iHosp_notvacc) <- sum(IBarHosp_to_H[, ]) 
  output(iHosp_vacc) <- sum(IBarHospv_to_Hv[,,])
  output(iHosp_notvacc_age[]) <- sum(IBarHosp_to_H[,i])
  output(iHosp_vacc_age[]) <- sum(IBarHospv_to_Hv[,i,])
  
  
  output(iInf_bynotvacc) <- sum(iInf_bynotvacc_tmp[,,]) + sum(iInf_bynotvaccv_tmp[,,,])
  output(iInf_byvacc) <- sum(iInf_byvacc_tmp[,,]) + sum(iInf_byvaccv_tmp[,,,])
  output(iInf_bynotvacc_age[]) <- sum(iInf_bynotvacc_tmp[,,i]) + sum(iInf_bynotvaccv_tmp[,,,i])
  output(iInf_byvacc_age[]) <- sum(iInf_byvacc_tmp[,,i]) + sum(iInf_byvaccv_tmp[,,,i])
  
  output(iInf_vacc) <- sum(Sv_to_E1v[,,])
  output(iInf_notvacc) <- sum(S_to_E1[,]) 
  output(iInf_vacc_age[]) <- sum(Sv_to_E1v[,i,])
  output(iInf_notvacc_age[]) <- sum(S_to_E1[,i])
  
  output(Vacc) <- sum(Sv[,]) + sum(E1v[,,]) + sum(E2v[,,]) + sum(E2v_iso[,,])+ sum(IMildv[,,]) +sum(IMildv_iso[,,]) + sum(Rv[,,])+
    sum(IHospv[,,]) + sum(IHospv_iso[,,]) + sum(IBarHospv[,,]) +  sum(Hv[,,])
  
  output(Vacc_age[]) <- sum(Sv[i,]) + sum(E1v[,i,]) + sum(E2v[,i,]) + sum(E2v_iso[,i,])+ sum(IMildv[,i,]) +sum(IMildv_iso[,i,]) + sum(Rv[,i,])+
    sum(IHospv[,i,]) + sum(IHospv_iso[,i,]) + sum(IBarHospv[,i,]) +  sum(Hv[,i,])
  
  ##############################################################################
  # User defined parameters (default in parentheses)
  ##############################################################################
  g1 <- user(1. / 4.)
  g2 <- user(1. / 1.)
  g3 <- user(1. / 3.)
  gToHosp <- user(1. / 4.)
  gHtoR <- user(1. / 13.)
  pop[] <- user()
  alpha_variants[] <- user() # Increase in transmissibility
  VE_sev[,] <- user()
  VE_susc[,] <- user()
  VE_inf[,] <- user()
  pHosp[, ] <- user()
  beta_reduction <- user()
  red_unvacc <- user()
  testing_rate[] <- user()
  testing_rate_v[] <- user()
  
  S_ini[] <- user()
  E1_ini[,] <- user()
  E2_ini[,] <- user()
  E2_iso_ini[,] <- user()
  IMild_ini[,] <- user()
  IMild_iso_ini[,] <- user()
  R_ini[,] <- user()
  IHosp_ini[,] <- user()
  IHosp_iso_ini[,] <- user()
  IBarHosp_ini[,] <- user()
  H_ini[,] <- user()
  
  Sv_ini[,] <- user()
  E1v_ini[,,] <- user()
  E2v_ini[,,] <- user()
  E2v_iso_ini[,,] <- user()
  IMildv_ini[,,] <- user()
  IMildv_iso_ini[,,] <- user()
  Rv_ini[,,] <- user()
  IHospv_ini[,,] <- user()
  IHospv_iso_ini[,,] <- user()
  IBarHospv_ini[,,] <- user()
  Hv_ini[,,] <- user()
  
  interp_ts_R0[] <- user()
  interp_R0s[] <- user()
  interp_ts_c_mat[] <- user()
  interp_c_mats[,,] <- user()
  interp_ts_v_rate[] <- user()
  interp_v_rate[,,] <- user()
  
  ##############################################################################
  # Dimensions
  ##############################################################################
  dim(pop) <- user()
  n_age <- length(pop)
  dim(alpha_variants) <- user()
  n_variants <- length(alpha_variants)
  
  dim(interp_ts_R0) <- user()
  dim(interp_R0s) <- length(interp_ts_R0)
  dim(interp_ts_c_mat) <- user()
  dim(interp_c_mats) <- c(length(interp_ts_c_mat), n_age, n_age)
  dim(interp_ts_v_rate) <- user()
  dim(interp_v_rate) <- user()
  n_vaccines <- dim(interp_v_rate,3) 
  dim(testing_rate) <- n_age
  dim(testing_rate_v) <- n_age
  
  dim(VE_sev) <- c(n_variants, n_vaccines)
  dim(VE_susc) <-c(n_variants, n_vaccines)
  dim(VE_inf) <- c(n_variants, n_vaccines)
  
  dim(beta) <- c(n_variants)
  dim(c_mat) <- c(n_age, n_age)
  dim(v_rate) <- c(n_age, n_vaccines)
  dim(pHosp) <- c(n_variants, n_age)
  
  dim(prop_Iv) <- c(n_variants, n_age, n_vaccines)
  dim(prop_Iv_iso) <- c(n_variants, n_age, n_vaccines)
  dim(prop_I) <- c(n_variants, n_age)
  dim(c_prop_I) <- c(n_variants, n_age, n_age)
  
  dim(prop_I_notvacc) <- c(n_variants, n_age)
  dim(c_prop_I_notvacc) <- c(n_variants, n_age, n_age)
  dim(iInf_bynotvacc_tmp) <- c(n_variants, n_age,n_age)
  dim(iInf_bynotvaccv_tmp) <- c(n_variants, n_age, n_vaccines, n_age)
  dim(prop_I_vacc) <- c(n_variants, n_age)
  dim(c_prop_I_vacc) <- c(n_variants, n_age, n_age)
  dim(iInf_byvacc_tmp) <- c(n_variants, n_age, n_age)
  dim(iInf_byvaccv_tmp) <- c(n_variants, n_age, n_vaccines, n_age)
  
  dim(S_to_E1) <- c(n_variants, n_age)
  dim(E1_to_E2) <- c(n_variants, n_age)
  dim(E2_to_I) <- c(n_variants, n_age)
  dim(E2_to_E2_iso) <- c(n_variants, n_age)
  dim(E2_iso_to_I_iso) <- c(n_variants, n_age)
  dim(E2_to_IMild) <- c(n_variants, n_age)
  dim(E2_iso_to_IMild_iso) <- c(n_variants, n_age)
  dim(IMild_to_R) <- c(n_variants, n_age)
  dim(IMild_to_IMild_iso) <- c(n_variants, n_age)
  dim(IMild_iso_to_R) <- c(n_variants, n_age)
  dim(E2_to_IHosp) <- c(n_variants, n_age)
  dim(E2_iso_to_IHosp_iso) <- c(n_variants, n_age)
  dim(IHosp_to_IBarHosp) <- c(n_variants, n_age)
  dim(IHosp_to_IHosp_iso) <- c(n_variants, n_age)
  dim(IHosp_iso_to_IBarHosp) <- c(n_variants, n_age)
  dim(IBarHosp_to_H) <- c(n_variants, n_age)
  dim(H_to_R) <- c(n_variants, n_age)
  
  dim(S) <- n_age
  dim(E1) <- c(n_variants, n_age)
  dim(E2) <- c(n_variants, n_age)
  dim(E2_iso) <- c(n_variants, n_age)
  dim(IMild) <- c(n_variants, n_age)
  dim(IMild_iso) <- c(n_variants, n_age)
  dim(R) <- c(n_variants, n_age)
  dim(IHosp) <- c(n_variants, n_age)
  dim(IHosp_iso) <- c(n_variants, n_age)
  dim(IBarHosp) <- c(n_variants, n_age)
  dim(H) <- c(n_variants, n_age)
  dim(S_ini) <- n_age
  dim(E1_ini) <- c(n_variants, n_age)
  dim(E2_ini) <- c(n_variants, n_age)
  dim(E2_iso_ini) <- c(n_variants, n_age)
  dim(IMild_ini) <- c(n_variants, n_age)
  dim(IMild_iso_ini) <- c(n_variants, n_age)
  dim(R_ini) <- c(n_variants, n_age)
  dim(IHosp_ini) <- c(n_variants, n_age)
  dim(IHosp_iso_ini) <- c(n_variants, n_age)
  dim(IBarHosp_ini) <- c(n_variants, n_age)
  dim(H_ini) <- c(n_variants, n_age)
  
  
  dim(S_to_Sv) <- c(n_age, n_vaccines)
  dim(E1_to_E1v) <- c(n_variants, n_age, n_vaccines)
  dim(E2_to_E2v) <- c(n_variants, n_age, n_vaccines)
  dim(E2_iso_to_E2v_iso) <- c(n_variants, n_age, n_vaccines)
  dim(IMild_to_IMildv) <- c(n_variants, n_age, n_vaccines)
  dim(IMild_iso_to_IMildv_iso) <- c(n_variants, n_age, n_vaccines)
  dim(IHosp_to_IHospv) <- c(n_variants, n_age, n_vaccines)
  dim(IHosp_iso_to_IHospv_iso) <- c(n_variants, n_age, n_vaccines)
  dim(IBarHosp_to_IBarHospv) <- c(n_variants, n_age, n_vaccines)
  dim(H_to_Hv) <- c(n_variants, n_age, n_vaccines)
  dim(R_to_Rv) <- c(n_variants, n_age, n_vaccines)
  
  
  dim(Sv_to_E1v) <- c(n_variants, n_age, n_vaccines)
  dim(E1v_to_E2v) <- c(n_variants, n_age, n_vaccines)
  dim(E2v_to_Iv) <- c(n_variants, n_age, n_vaccines)
  dim(E2v_to_E2v_iso) <- c(n_variants, n_age, n_vaccines)
  dim(E2v_iso_to_Iv_iso) <- c(n_variants, n_age, n_vaccines)
  dim(E2v_to_IMildv) <- c(n_variants, n_age, n_vaccines)
  dim(E2v_iso_to_IMildv_iso) <- c(n_variants, n_age, n_vaccines)
  dim(IMildv_to_Rv) <- c(n_variants, n_age, n_vaccines)
  dim(IMildv_to_IMildv_iso) <- c(n_variants, n_age, n_vaccines)
  dim(IMildv_iso_to_Rv) <- c(n_variants, n_age, n_vaccines)
  dim(E2v_to_IHospv) <- c(n_variants, n_age, n_vaccines)
  dim(E2v_iso_to_IHospv_iso) <- c(n_variants, n_age, n_vaccines)
  dim(IHospv_to_IBarHospv) <- c(n_variants, n_age, n_vaccines)
  dim(IHospv_to_IHospv_iso) <- c(n_variants, n_age, n_vaccines)
  dim(IHospv_iso_to_IBarHospv) <- c(n_variants, n_age, n_vaccines)
  dim(IBarHospv_to_Hv) <- c(n_variants, n_age, n_vaccines)
  dim(Hv_to_Rv) <- c(n_variants, n_age, n_vaccines)
  
  
  dim(Sv) <- c(n_age, n_vaccines)
  dim(E1v) <- c(n_variants, n_age, n_vaccines)
  dim(E2v) <- c(n_variants, n_age, n_vaccines)
  dim(E2v_iso) <- c(n_variants, n_age, n_vaccines)
  dim(IMildv) <- c(n_variants, n_age, n_vaccines)
  dim(IMildv_iso) <- c(n_variants, n_age, n_vaccines)
  dim(Rv) <- c(n_variants, n_age, n_vaccines)
  dim(IHospv) <- c(n_variants, n_age, n_vaccines)
  dim(IHospv_iso) <- c(n_variants, n_age, n_vaccines)
  dim(IBarHospv) <- c(n_variants, n_age, n_vaccines)
  dim(Hv) <- c(n_variants, n_age, n_vaccines)
  dim(Sv_ini) <- c(n_age, n_vaccines)
  dim(E1v_ini) <- c(n_variants, n_age, n_vaccines)
  dim(E2v_ini) <- c(n_variants, n_age, n_vaccines)
  dim(E2v_iso_ini) <- c(n_variants, n_age, n_vaccines)
  dim(IMildv_ini) <- c(n_variants, n_age, n_vaccines)
  dim(IMildv_iso_ini) <- c(n_variants, n_age, n_vaccines)
  dim(Rv_ini) <- c(n_variants, n_age, n_vaccines)
  dim(IHospv_ini) <- c(n_variants, n_age, n_vaccines)
  dim(IHospv_iso_ini) <- c(n_variants, n_age, n_vaccines)
  dim(IBarHospv_ini) <- c(n_variants, n_age, n_vaccines)
  dim(Hv_ini) <- c(n_variants, n_age, n_vaccines)
  
 
  dim(iInf_v) <- c(n_variants)
  dim(iHosp_v) <- c(n_variants)
  dim(iHosp_age) <- c(n_age)
  dim(Vacc_age) <- c(n_age)
  dim(iHosp_notvacc_age) <- c(n_age)
  dim(iHosp_vacc_age) <- c(n_age)
  dim(iInf_vacc_age) <- c(n_age)
  dim(iInf_notvacc_age) <- c(n_age)
  dim(iInf_bynotvacc_age) <- c(n_age)
  dim(iInf_byvacc_age) <- c(n_age)
})
