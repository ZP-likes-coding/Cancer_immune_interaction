functions {
  vector cancer_ode(real t, vector y, vector theta, data vector fixed_p) {
    vector[5] dydt;
    
    // Unpack variables for clarity
    real C  = y[1]; real D  = y[2]; real N  = y[3];
    real Tr = y[4]; real T8 = y[5];
    
    // Unpack parameters to fit
    real lambda_T8_comb = theta[1];
    real lambda_DC_comb = theta[2];
    real d_T8           = theta[3];
    real sigma_N        = theta[4];
    real K_Tr           = theta[5]; // k_treg
    real d_N            = theta[6];
    
    // Unpack fixed parameters (passed from R)
    real lambda_C = fixed_p[1];  real C_M    = fixed_p[2];
    real eta_8    = fixed_p[3];  real eta_N  = fixed_p[4];
    real d_C      = fixed_p[5];  real K_C    = fixed_p[6];
    real d_D      = fixed_p[7];  real gamma_N = fixed_p[8];
    real a_C      = fixed_p[9];  real beta_1 = fixed_p[10];
    real beta_2   = fixed_p[11]; real d_Tr   = fixed_p[12];
    real lambda_Tr_comb = fixed_p[13]; real K_D = fixed_p[14];

    // The ODE system
    dydt[1] = lambda_C * C * (1 - C / C_M) - eta_8 * T8 * C - eta_8 * N * C - d_C * C;
    dydt[2] = lambda_DC_comb * C / (C + K_C) - d_D * D;
    dydt[3] = sigma_N - d_N * N - gamma_N * Tr * N + a_C * N * C / (1 + C/beta_1 + N/beta_2);
    dydt[4] = -d_Tr * Tr + lambda_Tr_comb * C / (K_C + C);
    dydt[5] = -d_T8 * T8 + lambda_T8_comb * D / (K_D + D) / (1 + Tr / K_Tr);
    
    return dydt;
  }
}

data {
  int<lower=1> N_obs;               // Number of data points
  array[N_obs] real ts;             // Observed time points
  array[N_obs] vector[4] y_obs;     // Measured: D, N, Tr, T8 (Order from CSV)
  vector[5] y0;                     // Initial conditions at t=0
  vector[14] fixed_p;               // The 14 parameters to keep fixed
  vector[6] mu_priors;      // Centers for your 6 parameters
  vector[6] sigma_priors;   // Widths (uncertainty) for your 6 parameters
}


parameters {
  vector<lower=0>[6] theta;         // The 6 parameters to fit
  real<lower=0> sigma;              // Measurement error noise
}

transformed parameters {
  // Solve ODE for all time points
  // We solve starting from t=0.0 to the max time in ts
  array[N_obs] vector[5] mu = ode_rk45(cancer_ode, y0, 0.0, ts, theta, fixed_p);
}

model {
  // Each parameter gets its own specific normal distribution
  theta ~ normal(mu_priors, sigma_priors); 
  sigma ~ exponential(1);
  
  // Likelihood
  for (n in 1:N_obs) {
    // We map mu[n] states (2,3,4,5) to y_obs (D, N, Tr, T8)
    y_obs[n] ~ normal([mu[n, 2], mu[n, 3], mu[n, 4], mu[n, 5]]', sigma);
  }
}


