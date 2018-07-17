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

Xm <- unname( model.matrix(~ 1 + Age + (Sexe + Stade) * IODAnnual, data))
Xv <- unname( model.matrix(~ 1 + Age + (Sexe + Stade) * IODAnnual, data))

##### creating BUGS datalist #####
stan.data <- list(
    Km = ncol(Xm),
    Kv = ncol(Xv),
    N = nrow(Xm),
    nid = length(unique(data$ids)),
    ## npr = length(unique(data$PairID)),
    nyr = length(unique(data$Cycle)),
    id = as.numeric(droplevels(data$ids)),
    ## pr = as.numeric(droplevels(as.factor(data$PairID))),
    yr = data$Cycle.n,
    Xm = Xm,
    Xv = Xv,
    Y = data$Tripdur
)


dh_nopair <- stan(file="Stan_Model_nopair.stan", data = stan.data, iter = 4000, warmup = 2000, thin = 1, chains = 4, control = list(adapt_delta = 0.95))
save(dh_nopair,file="dhglm_nopair.rda")
