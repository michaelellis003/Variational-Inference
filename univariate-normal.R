vi <- function(
    X,
    mu_0 = 0,
    sigmasq_0 = 10^8,
    alpha_0 = 0.01,
    beta_0 = 0.01,
    max_iter = 10
) {
    
    N <- length(X) # sample size
    sum_X <- sum(X)
    
    ## initialize B_{q(\sigma^2)}
    sigmasq_q_mu <- 0
    mu_q_mu <- 0
    A_q_sigmasq <- 1
    B_q_sigmasq <- 1
    
    ## save parameters
    mu_q_mu_save <- rep(NA, max_iter)
    sigmasq_q_mu_save <- rep(NA, max_iter)
    B_q_sigmasq_save <- rep(NA, max_iter)
    elbos <- rep(NA, max_iter)
    
    ## Calculate shape parameter for sigma^2 distribution
    A_q_sigmasq <- alpha_0 + N/2
    
    for(i in 1:max_iter) {
        
        E_sigmasq <- B_q_sigmasq / A_q_sigmasq
        
        ## update 
        sigmasq_q_mu <- 1 / (N/E_sigmasq + 1/sigmasq_0)
        mu_q_mu <- sigmasq_q_mu * (sum_X/E_sigmasq + mu_0/sigmasq_0)
        
        sum_sq <- sum((X - mu_q_mu)^2)
        B_q_sigmasq <- beta_0 + (1/2) * (sum_sq + N*sigmasq_q_mu)
        
        elbos[i] <- ELBO(N, sigmasq_q_mu, mu_q_mu, B_q_sigmasq,
                         mu_0, sigmasq_0, alpha_0, beta_0)
        
    }
    
    return(
        list(
            mu_q_mu = mu_q_mu,
            sigmasq_q_mu = sigmasq_q_mu,
            A_q_sigmasq = A_q_sigmasq,
            B_q_sigmasq = B_q_sigmasq,
            elbos = elbos
        )
    )
    
}

ELBO <- function(
    N,
    sigmasq_q_mu,
    mu_q_mu,
    B_q_sigmasq,
    mu_0,
    sigmasq_0,
    alpha_0,
    beta_0
) {
    
    ELBO <- 1/2 - (N/2)*log(2*pi) + (1/2)*log(sigmasq_q_mu/sigmasq_0) - 
        (((mu_q_mu - mu_0)^2 + sigmasq_q_mu)/(2*sigmasq_0)) + 
        alpha_0*log(beta_0) - (alpha_0 + N/2)*log(B_q_sigmasq) + 
        log(gamma(alpha_0 + N/2)) - log(gamma(alpha_0))
    
    return(ELBO)
    
}