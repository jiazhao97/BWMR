###### VEM and Inference for BWMR ######
##
### INPUT
## gammahat:                        SNP-exposure effect;
## Gammahat:                        SNP-outcome effect;
## sigmaX:                          standard error of SNP-exposure effect;
## sigmaY:                          standard error of SNP-outcome effect;
##
### OUTPUT
## beta                             the estimate of beta;
## se_beta                          the estimate the standard error of beta;
## P-value                          P-value
## plot1                            Plot of Data with Standard Error Bar;
## plot2                            Plot of Evidence Lower Bound (ELBO);
## plot3                            Posterior Mean of Weight of Each Observation;
## plot4                            Plot of Weighted Data and Its Regression Result.


# packages
library("ggplot2")


# known parameters for the prior distributions
sqsigma0 <- (1e+6)^2
alpha <- 100


## define function to calculate ELBO and E[Lc] (the approximate log-likelihood)
ELBO_func <- function(N, gammahat, Gammahat, sqsigmaX, sqsigmaY, mu_beta, sqsigma_beta, mu_gamma, sqsigma_gamma, a, b, pi_w, sqsigma, sqtau) {
  # + E[log p(gammahat, sqsigmaX | gamma)]
  l <- - 0.5*sum(log(sqsigmaX)) - 0.5*sum(((gammahat-mu_gamma)^2+sqsigma_gamma)/sqsigmaX)    
  # + E[log p(Gammahat, sqsigmaY| beta, gamma, w, sqtau)]
  l <- l - 0.5*log(2*pi)*sum(pi_w) - 0.5*sum(pi_w*log(sqsigmaY+sqtau)) 
  l <- l - 0.5*sum(pi_w*((mu_beta^2+sqsigma_beta)*(mu_gamma^2+sqsigma_gamma)-2*mu_beta*mu_gamma*Gammahat+Gammahat^2)/(sqsigmaY+sqtau))
  # + E[log p(beta | sqsigma0)]
  l <- l - 0.5*(mu_beta^2+sqsigma_beta)/sqsigma0
  # + E[log p(gamma | sqsigma)]
  l <- l - 0.5*N*log(sqsigma) - 0.5/sqsigma*sum(mu_gamma^2+sqsigma_gamma)
  # + E[log p(w | pi1)]
  l <- l + (digamma(a)-digamma(a+b))*sum(pi_w) + (digamma(b)-digamma(a+b))*(N-sum(pi_w))
  # + E[log p(pi1)]
  l <- l + (alpha-1)*(digamma(a)-digamma(a+b))
  
  # - E[log q(beta)]
  l <- l + 0.5*log(sqsigma_beta)
  # - E[log q(gamma)]
  l <- l + 0.5*sum(log(sqsigma_gamma))
  # - E[log q(pi1)]
  l <- l - (a-1)*(digamma(a)-digamma(a+b)) - (b-1)*(digamma(b)-digamma(a+b)) + lbeta(a, b)
  # - E[log q(w)]
  # Need to check if pi_w = 0 or pi_w = 1, since there are log terms of pi_w and 1 - pi_w.
  # e1 <- pi_w*log(pi_w)
  # e2 <- (1-pi_w)*log(1-pi_w)
  # e1[which(pi_w == 0)] <- 0
  # e2[which(pi_w == 1)] <- 0
  # l <- l - sum(e1+e2)
  ## A TRICK TO SIMPLIFY THE CODE
  l <- l - sum(pi_w*log(pi_w+(pi_w==0)) + (1-pi_w)*log(1-pi_w+(pi_w==1)))
  # l: ELBO
}


BWMR <- function(gammahat, Gammahat, sigmaX, sigmaY) {
  ## data
  N <- length(gammahat)
  sqsigmaX <- sigmaX^2
  sqsigmaY <- sigmaY^2
  
  ### Variational EM algorithm ###
  # initialize
  # initial parameters of BWMR
  beta <- 0
  sqtau <- 1^2
  sqsigma <- 1^2
  # initial parameters of variational distribution
  mu_gamma <- gammahat
  sqsigma_gamma <- rep(0.1, N)
  pi_w <- rep(0.5, N)
  # declare sets of ELBO and approximate log-likelihood
  ELBO_set <- numeric(0)
  
  for (iter in 1:5000) {
    ## Variational E-Step
    # beta
    sqsigma_beta <- 1/(1/sqsigma0 + sum(pi_w*(mu_gamma^2+sqsigma_gamma)/(sqsigmaY+sqtau)))
    mu_beta <- sum(pi_w*mu_gamma*Gammahat/(sqsigmaY+sqtau))*sqsigma_beta
    # gamma
    sqsigma_gamma <- 1/(1/sqsigmaX + pi_w*(mu_beta^2+sqsigma_beta)/(sqsigmaY+sqtau) + 1/sqsigma)
    mu_gamma <- (gammahat/sqsigmaX + pi_w*Gammahat*mu_beta/(sqsigmaY+sqtau))*sqsigma_gamma
    # pi1
    a <- alpha + sum(pi_w)
    b <- N + 1 - sum(pi_w)
    # w
    q0 <- exp(digamma(b) - digamma(a+b))
    q1 <- exp(- 0.5*log(2*pi) - 0.5*log(sqsigmaY+sqtau) - 0.5*((mu_beta^2+sqsigma_beta)*(mu_gamma^2+sqsigma_gamma)-2*mu_beta*mu_gamma*Gammahat+Gammahat^2)/(sqsigmaY+sqtau) + digamma(a)-digamma(a+b))
    pi_w <- q1/(q0+q1)
    
    if (sum(pi_w) == 0){
      message("Invalid IVs!")
      mu_beta = NA
      se_beta = NA
      P_value = NA
      return(list(beta=NA, se_beta=NA, P_value=NA))
    }
    
    ## M-Step
    sqsigma <- sum(mu_gamma^2 + sqsigma_gamma)/N
    sqtau <- sum(pi_w*((mu_beta^2+sqsigma_beta)*(mu_gamma^2+sqsigma_gamma)-2*mu_beta*mu_gamma*Gammahat+Gammahat^2)*sqtau^2/((sqsigmaY+sqtau)^2)) / sum(pi_w/(sqsigmaY+sqtau))
    sqtau <- sqrt(sqtau)
    
    ## check ELBO
    ELBO <- ELBO_func(N, gammahat, Gammahat, sqsigmaX, sqsigmaY, mu_beta, sqsigma_beta, mu_gamma, sqsigma_gamma, a, b, pi_w, sqsigma, sqtau)
    ELBO_set <- c(ELBO_set, ELBO)
    if (iter > 1 && (abs((ELBO_set[iter]-ELBO_set[iter-1])/ELBO_set[iter-1]) < 1e-6)) {
      break
    }
  }
  # message("Iteration=", iter, ", beta=", mu_beta, ", tau=", sqrt(sqtau), ", sigma=", sqrt(sqsigma), ".")
  
  
  
  ### visualize the result of VEM algorithm
  # Plot1: Plot of Data with Standard Error Bar
  df1 <- data.frame(
    gammahat = gammahat,
    Gammahat = Gammahat,
    sigmaX = sigmaX,
    sigmaY = sigmaY
  )
  plot1 <- ggplot(data = df1, aes(x = gammahat, y = Gammahat)) +  
    geom_pointrange(aes(ymin = Gammahat - sigmaY, ymax = Gammahat + sigmaY), color="gray59", size = 0.3) +
    geom_errorbarh(aes(xmin = gammahat - sigmaX, xmax = gammahat + sigmaX, height = 0), color="gray59") +
    labs(x = "SNP-exposure effect", y = "SNP-outcome effect", title = "Plot1: Plot of data with standard error bar")
  
  # Plot2: Plot of Evidence Lower Bound (ELBO)
  iteration <- seq(1, (length(ELBO_set)), by = 1)
  df2 <- data.frame(
    iteration = iteration,
    ELBO_iter = ELBO_set
  )
  plot2 <- ggplot(df2, aes(x=iteration, y=ELBO_iter)) + geom_line(size = 0.5, color = "tomato1") + geom_point(size=0.5, color = "tomato1") +
    labs(x = "iteration", y="elbo", title = "Plot2: Plot of evidence lower bound (elbo)")
  
  # Plot3: Posterior Mean of Weight of Each Observation
  serial_number <- seq(1, N, by = 1)
  df3 <- data.frame(
    weight = pi_w,
    serial_number = serial_number
  )
  plot3 <- ggplot(data = df3, mapping = aes(x = factor(serial_number), y = weight, fill = weight)) + geom_bar(stat = 'identity', position = 'dodge') +
    labs(x = "observation No.", y = "weight", title = "Plot3: Posterior mean of weight of each observation") +
    ylim(0, 1) +
    theme(axis.text.x = element_text(size = 5))
  # scale_x_discrete(breaks = seq(10, N, 20)) +
  
  # Plot4: Plot of Weighted Data and Its Regression Result
  df4 <- data.frame(
    gammahat = gammahat,
    Gammahat = Gammahat,
    sqsigmaX = sqsigmaX,
    sqsigmaY = sqsigmaY,
    w = pi_w
  )
  plot4 <- ggplot(df4, aes(x=gammahat, y=Gammahat, color=w)) + geom_point(size = 0.3) +
    geom_pointrange(aes(ymin = Gammahat - sigmaY, ymax = Gammahat + sigmaY), size = 0.3) +
    geom_errorbarh(aes(xmin = gammahat - sigmaX, xmax = gammahat + sigmaX, height = 0)) +
    geom_abline(intercept=0, slope=mu_beta, color="#990000", linetype="dashed", size=0.5) +
    labs(x = "SNP-exposure effect", y = "SNP-outcome effect", title = "Plot4: Plot of weighted data and its regression result")
  
  
  ### LRVB and Standard Error ###
  ## matrix V
  forV <- matrix(nrow = N, ncol = 4)
  forV[ ,1] <- sqsigma_gamma
  forV[ ,2] <- 2*mu_gamma*sqsigma_gamma
  forV[ ,3] <- forV[ ,2]
  forV[ ,4] <- 2*sqsigma_gamma^2 + 4*mu_gamma^2*sqsigma_gamma
  V <- matrix(rep(0, (3*N+4)*(3*N+4)), nrow = 3*N+4, ncol = 3*N+4)
  for (j in 1:N) {
    V[(3*j):(3*j+1), (3*j):(3*j+1)] <- matrix(forV[j, ], 2, 2)
    V[3*j+2, 3*j+2] <- pi_w[j] - (pi_w[j]^2)
  }
  V[1:2, 1:2] <- matrix(c(sqsigma_beta, 2*mu_beta*sqsigma_beta, 2*mu_beta*sqsigma_beta, 2*sqsigma_beta^2+4*mu_beta^2*sqsigma_beta), 2, 2)
  V[(3*N+3):(3*N+4), (3*N+3):(3*N+4)] <- matrix(c(trigamma(a)-trigamma(a+b), -trigamma(a+b), -trigamma(a+b), trigamma(b)-trigamma(a+b)), 2, 2)
  
  ## matrix H
  H <- matrix(rep(0, (3*N+4)*(3*N+4)), nrow = 3*N+4, ncol = 3*N+4)
  forH <- matrix(nrow = N, ncol = 6)
  forH[ ,1] <- pi_w*Gammahat/(sqsigmaY+sqtau)
  forH[ ,2] <- mu_gamma*Gammahat/(sqsigmaY+sqtau)
  forH[ ,3] <- -0.5*pi_w/(sqsigmaY+sqtau)
  forH[ ,4] <- -0.5*(mu_gamma^2+sqsigma_gamma)/(sqsigmaY+sqtau)
  forH[ ,5] <- mu_beta*Gammahat/(sqsigmaY+sqtau)
  forH[ ,6] <- -0.5*(mu_beta^2+sqsigma_beta)/(sqsigmaY+sqtau)
  for (j in 1:N) {
    H[1, 3*j] <- forH[j, 1]
    H[1, 3*j+2] <- forH[j, 2]
    H[2, 3*j+1] <- forH[j, 3]
    H[2, 3*j+2] <- forH[j, 4]
    H[(3*N+3):(3*N+4), 3*j+2] <- c(1, -1)
    H[3*j+2, (3*N+3):(3*N+4)] <- c(1, -1)
    H[(3*j):(3*j+1), 3*j+2] <- forH[j, 5:6]
    H[3*j+2, (3*j):(3*j+1)] <- forH[j, 5:6]
  }
  H[ ,1] <- H[1, ]
  H[ ,2] <- H[2, ]
  
  
  ## accurate covariance estimate and standard error
  I <- diag(3*N+4)
  Sigma_hat <- try(solve(I-V%*%H)%*%V)
  
  if (inherits(Sigma_hat, "try-error")){
    message("Invalid IVs!")
    return(list(beta=NA, se_beta=NA, P_value=NA))
  } else{
    se_beta <- sqrt(Sigma_hat[1, 1])
  }
  
  
  ## test
  W <- (mu_beta/se_beta)^2
  # P_value <- 1 - pchisq(W, 1)
  P_value <- pchisq(W, 1, lower.tail=F)
  
  message("Estimate of beta=", mu_beta, ", se of beta=", se_beta, ", P-value=", P_value, ".")
  
  ## output
  output <- list(beta=mu_beta, se_beta=se_beta, P_value=P_value, 
                 weights=pi_w, tau=sqrt(sqtau), sigma=sqrt(sqsigma), mu_pi=a/(a+b),
                 plot1=plot1, plot2=plot2, plot3=plot3, plot4=plot4)
}
