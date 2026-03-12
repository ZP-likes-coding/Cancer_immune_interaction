library(rstan)
library(dplyr)

# 1. Load and clean data
df <- read.csv("T11_ICI.csv")

# Map text timepoints to numbers (assuming 'end' is day 21)
df <- df %>%
  mutate(time_num = case_when(
    Timepoint == "day3" ~ 3,
    Timepoint == "day7" ~ 7,
    Timepoint == "end"  ~ 14,
    TRUE ~ 0
  )) %>%
  arrange(time_num)

# 2. Extract observed vectors (D, N, Tr, T8)
y_obs_matrix <- df %>%
  select(Dendritic_quanTIseq, NK_quanTIseq, Tregs_quanTIseq, T.CD8_quanTIseq) %>%
  as.matrix()

# 3. Define the 14 Fixed Parameters (Order must match the Stan block)
# Replace these with your actual values
fixed_params <- c(
  lambda_C = 1.5, C_M = 0.8, eta_8 = 328.55, eta_N = 300,
  d_C = 0.17, K_C = 0.4, d_D = 0.1, gamma_N = 150,
  a_C = 0.5, beta_1 = 0.4, beta_2 = 0.0002, d_Tr = 0.2,
  lambda_Tr_comb = 0.0002, K_D = 0.0004
)

# 4. Define Initial Conditions (C, D, N, Tr, T8 at t=0)
y0 <- c(C=1000, D=0.001, N=0.001, Tr=0.0001, T8=0.0001)

# Define centers for: lambda_T8_comb, lambda_DC_comb, d_T8, sigma_N, K_Tr, d_N
mu_priors <- c(0.00108, 0.00008, 0.18, 0.00005, 0.00025,0.1) 
sigma_priors <- c(
  0.00108,  # for lambda_T8_comb
  0.00008,  # for lambda_DC_comb
  0.18,     # for d_T8
  0.00005,  # for sigma_N
  0.00025,  # for K_Tr
  0.1       # for d_N
) # priors are weakly informative 

# 5. Prepare data list for Stan
stan_data <- list(
  N_obs = nrow(df),
  ts = df$time_num,
  y_obs = y_obs_matrix,
  y0 = y0,
  fixed_p = fixed_params,
  mu_priors = mu_priors,    # Pass them here
  sigma_priors = sigma_priors
)

# 6. Run MCMC
fit <- stan(
  file = "cancer_model.stan",
  data = stan_data,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = parallel::detectCores()
)

# 7. Check results
print(fit, pars = c("theta", "sigma"))
plot(fit, pars = "theta")