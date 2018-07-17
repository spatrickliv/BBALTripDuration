###### pacakges ######
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

###### data ######
data <- read.table("bugsdata5.txt")
## no_of_cores <- 3

data <- subset(data, Tripdur < 600) # remove trips longer than 600 hours
### removing NAs
I <- is.na(data$ids) |
  is.na(data$Sexe) |
  is.na(data$Cycle) |
  is.na(data$CycleBand) |
  is.na(data$Stade) |
  is.na(data$Tripdur) |
  is.na(data$Age) |
  is.na(data$IODAnnual) 
 
data <- data[!I,]

### scaling and centering
data$Tripdur <- as.numeric(scale(sqrt(data$Tripdur)))
data$Age <- as.numeric(scale(data$Age))
data$IODAnnual <- as.numeric(scale(data$IODAnnual))
data$Stade <- relevel(data$Stade, ref = "INCUBATION")

db <- data.frame(Cycle = sort(unique(data$Cycle)), Cycle.n = 1:length(unique(data$Cycle)))
data$Cycle.n <- db$Cycle.n[match(data$Cycle,db$Cycle)]

##### creating BUGS datalist #####
mk.stan.data <- function(Xm, Xv, dat){
    out <- list(
        Km = ncol(Xm),
        Kv = ncol(Xv),
        N = nrow(Xm),
        nid = length(unique(dat$ids)),
        ## npr = length(unique(data$PairID)),
        nyr = length(unique(dat$Cycle)),
        id = as.numeric(droplevels(dat$ids)),
        ## pr = as.numeric(droplevels(as.factor(data$PairID))),
        yr = dat$Cycle.n,
        Xm = unname(Xm),
        Xv = unname(Xv),
        Y = dat$Tripdur
    )
    out
}

Xm.full <- model.matrix(~ 1 + Age + (Sexe + Stade) * IODAnnual, data)
Xv.full <- model.matrix(~ 1 + Age + (Sexe + Stade) * IODAnnual, data)
stan.data.full <- mk.stan.data(Xm.full,Xv.full,data)

dh_full <- stan(file="Stan_Model_nopair.stan",
                  data = stan.data.full,
                  iter = 5000,
                  warmup = 3000,
                  thin = 5,
                  chains = 5,
                control = list(adapt_delta = 0.96))
save(dh_full,file="dhglm_full.rda")

load("dhglm_full.rda")

params <- c("beta","gamma","sigma_id","sigma_yr","L_O_id","L_O_yr", "Sigma_id","Sigma_yr","Omega_id[1,2]","Omega_yr[1,2]")
params2 <- c("beta","gamma", "Sigma_id","Sigma_yr","Omega_id[1,2]","Omega_yr[1,2]")

results <- summary(dh_full, pars = params, probs = c(0.025, 0.5, 0.975))$summary
rownames(results)[1:(ncol(Xm.full)+ncol(Xv.full))] <- c(paste0("b_", colnames(Xm.full)),
                            paste0("g_", colnames(Xv.full)))
print(results, digits = 2)

stan_plot(dh_full, show_density = TRUE, pars=params2) + ggtitle("Posterior distributions model QG")
stan_dens(dh_full, pars=params2) + ggtitle("Posterior distributions model QG")


Xm.sig <- model.matrix(~ 1 +  Age + Sexe + Stade * IODAnnual, data)
Xv.sig <- model.matrix(~ 1 + Age + Sexe + Stade + IODAnnual, data)
stan.data.sig <- mk.stan.data(Xm.sig, Xv.sig, data)

dh_sig <- stan(file="Stan_Model_nopair.stan",
                  data = stan.data.sig,
                  iter = 5000,
                  warmup = 3000,
                  thin = 5,
                  chains = 5,
                control = list(adapt_delta = 0.96))
save(dh_sig,file="dhglm_sig.rda")


Xm.Asig <- model.matrix(~ 1 +  Sexe + (Age + Stade) * IODAnnual, data)
Xv.Asig <- model.matrix(~ 1 + Sexe + Stade + Age * IODAnnual, data)
stan.data.Asig <- mk.stan.data(Xm.Asig, Xv.Asig, data)

dh_Asig <- stan(file="Stan_Model_nopair.stan",
                  data = stan.data.Asig,
                  iter = 5000,
                  warmup = 3000,
                  thin = 5,
                chains = 5,
                fit=NA,
                control = list(adapt_delta = 0.96))
save(dh_Asig,file="dhglm_Asig.rda")

Xm.Afull <- model.matrix(~ 1 + (Age + Sexe + Stade) * IODAnnual, data)
Xv.Afull <- model.matrix(~ 1 + (Age + Sexe + Stade) * IODAnnual, data)
stan.data.Afull <- mk.stan.data(Xm.Afull,Xv.Afull,data)

dh_Afull <- stan(file="Stan_Model_nopair.stan",
                  data = stan.data.Afull,
                  iter = 5000,
                  warmup = 3000,
                  thin = 5,
                  chains = 5,
                control = list(adapt_delta = 0.96))
save(dh_Afull,file="dhglm_Afull.rda")
