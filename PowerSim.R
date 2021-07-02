
# Load packages -----------------------------------------------------------
library(lme4)
library(lmerTest)
set.seed(10408) # for reproducibility 

# Set constants -----------------------------------------------------------

# Experimental/simulation parameters
N <- seq(10, 100, by=10)   # sample sizes to test
nTrial <- 100              # number of trials per subject
nSim <- 100                # number of simulations per sample size per subject
alpha <- .05               # alpha thereshold

# Model parameters
G00   <- 350.05            # gamma00 -- grand mean RT
G10   <- 35.55             # gamma11 -- congruency effect for Stroop (assuming effects coding)
tau20 <- 35000             # tau20 -- cluster-level variance about G00
tau21 <- 2500              # tau21 -- cluster-level variance about G10 
sigma2 <- 50000            # sigma2 -- within-subject variance

# Create design matrix ----------------------------------------------------

options <- c('green', 'blue', 'red', 'purple')
words <- sample(rep(options, each=nTrial/4))
cols  <- sapply(1:length(words), function(x) ifelse(x<=50, words[x], sample(words[words!=words[x]], size = 1)))
congruency <- ifelse(words==cols, -1, 1) # effects coding

stroop = data.frame(words, cols, congruency)

# Simluate ----------------------------------------------------

PowerSim <- data.frame(stringsAsFactors = F)

for(sim in 1:nSim) {

  for(thisN in N) {
    
    # Keep track of progress
    cat('simulation', sim, '/', nSim, 'for sample size:', thisN, '\n') 
    
    # Simulate!
    thisb0 <- G00 + rnorm(thisN, 0, sqrt(tau20))
    thisb1 <- G10 + rnorm(thisN, 0, sqrt(tau21))
    
    RTsim = c()
    for(sub in 1:thisN) {
      thisRT <- thisb0[sub] + thisb1[sub]*stroop$congruency + rnorm(nTrial, 0, sqrt(sigma2))
      RTsim <- c(RTsim, thisRT)
    }
    
    # Model! 
    thisStroop <- do.call("rbind", replicate(thisN, stroop, simplify = FALSE))
    thisStroop$RT <- RTsim
    thisStroop$id <- rep(1:thisN, each=nTrial)
    
    thisMLM <- lmer(RTsim ~ congruency + (congruency|id), data=thisStroop)
    
    # store estimates 
    b0 <- fixef(thisMLM)[['(Intercept)']]
    b1 <- fixef(thisMLM)[['congruency']]
    tau20_sim <- as.data.frame(VarCorr(thisMLM))[1,'vcov']
    tau21_sim <- as.data.frame(VarCorr(thisMLM))[2,'vcov']
    sigma2_sim <- sigma(thisMLM)^2
    p <- summary(thisMLM)$coefficients[,'Pr(>|t|)'][['congruency']]
    
    thisSim <- data.frame(thisN, sim, b0, b1, tau20_sim, tau21_sim, sigma2_sim, p)
    PowerSim <- rbind(PowerSim, thisSim)
    
  }
}


# Visualize ---------------------------------------------------------------

PowerSim$sig <- ifelse(PowerSim$p < alpha, 1, 0)
pwr <- tapply(PowerSim$sig, PowerSim$thisN, mean)

plot(names(pwr), pwr, type='b',
     xlab='Sample Size', ylab='p(Significant)', 
     main=paste0(nTrial, ' Trials\n alpha=', alpha))
legend('bottomright', legend=paste0(nSim, ' simulations'), bty='n')
abline(h=0.80, lty='dashed')
