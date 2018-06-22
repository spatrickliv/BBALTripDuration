#!/usr/bin/env Rscript

#### script arguments
args <- commandArgs(trailingOnly = TRUE)
if(!(length(args) %in% c(0,3))) {
    cat("\n When used a script, usage is:\n");
    cat("   <script_name> input_file output_file no_of_cores\n\n");
    q()
}

input_file  <- args[1]
output_file <- args[2]
no_of_cores <- as.integer(args[3])

if(file.access(input_file, mode=0)) {
   cat("input file is not readable\n\n")
}

data <- read.table(input_file, header=TRUE)

###### packages ######
library(runjags)
library(coda)

###### if on computer ######
if(exists("input.file") == FALSE) {
    data <- read.table("bugsdata5.txt")
    no_of_cores <- 3
}
###### data processing ######
### removing NAs
I <-is.na(data$ids)|
  is.na(data$Sexe)|
  is.na(data$Cycle)|
  is.na(data$CycleBand)|
  is.na(data$Rsnum)|
  is.na(data$RSplus1num)|
##  is.na(data$CyclePairBand)|
  is.na(data$Stade)|
  is.na(data$ntrajet)|
  is.na(data$Tripdur)|
  is.na(data$Age)|
  is.na(data$IODAnnual)| 
  is.na(data$SOIAnnual)|  
  is.na(data$AAOAnnual) ##|  
##  is.na(data$ONIAnnual)
 
data <- data[!I,]

###### Not needed ######
## birdID <- as.data.frame(data$ids)
## birdIDnum <- as.data.frame(as.numeric(droplevels(data$ids)))
## unique(birdIDnum)
## write.table(cbind(birdID,birdIDnum), "BirdIDs.txt")
## write.table(data, "croppedbugsdata2.txt")
######

### scaling and centering
## data$ntrajet<- as.numeric(scale(data$ntrajet))
data$Tripdur <- as.numeric(scale(sqrt(data$Tripdur)))
data$Age<- as.numeric(scale(data$Age))
data$IODAnnual<- as.numeric(scale(data$IODAnnual))
data$SOIAnnual<- as.numeric(scale(data$SOIAnnual))
data$AAOAnnual<- as.numeric(scale(data$AAOAnnual))
##data$ONIAnnual<- as.numeric(scale(data$ONIAnnual))

db <- data.frame(Cycle=sort(unique(data$Cycle)),Cycle.n=1:length(unique(data$Cycle)))
data$Cycle.n <- db$Cycle.n[match(data$Cycle,db$Cycle)]

### obtaining the Pairs (i.e. removing the year, first four characters, from CyclePairBand)
data$PairBand <- substring(data$CyclePairBand,5)

##### creating BUGS datalist #####
BUGS6.data <- list(n = nrow(data), #121
                   n.Pairs = length(unique(data$PairID)),  #n.pair
                   n.Cycles = length(unique(data$Cycle)), #n.cycle
                   n.ids = length(unique(data$ids)),
                   y = data$Tripdur,
                   ###### Sexe is coded as female 1 and male 2 now to code for the interaction more easily
                   Sexe = as.numeric(data$Sexe) - 1, ##female coded as 0 and Males as 1
                   Rsnum = data$Rsnum,
                   IODAnnual = data$IODAnnual,
                   SOIAnnual = data$SOIAnnual,
                   AAOAnnual = data$AAOAnnual,
                   ###### Same thing reproductive status is coded differently. instead of having one variable for each status now, we have only 1 variable varying between 1-3
                   Repro = as.numeric(data$Stade) , ##"BROODING"=1, "INCUBATION"=2, "PRE-PONTE"=3,
                   Brooding =ifelse(data$Stade=="BROODING",1,0), 
                   Preponte = ifelse(data$Stade=="PRE-PONTE",1,0),
                   Incubation = ifelse(data$Stade=="INCUBATION",1,0),
                   id = as.numeric(droplevels(data$ids)),
                   PairBand = as.numeric(droplevels(as.factor(data$PairID))),
                   Cycle = data$Cycle.n, 
                   Age = data$Age,
                   mat.id = diag(2),
                   zero = c(0,0),
                   mat.pair = diag(2),
                   zerop = c(0,0),
                   mat.yr = diag(2),
                   zeroyr = c(0,0)
)


##### number of fixed effects #####
nbeta <- 6
ngamma <- 6


##### parameters so save #####
BUGS6.parameters <- c("beta","gamma", "betaSi", "betaRi", "gammaSi", "gammaRi", "sig2.pair", "sig2.id", "sig2.yr", "tau.yr", "tau.pair", "tau.id") #,"uid")



##### intial values #####
BUGS6.inits.jags <- function() {list(beta = rnorm(nbeta), #fixed
                                     gamma = rnorm(ngamma), 
                                     betaSi = rnorm(1), # interaction with Sexe
                                     betaRi = rnorm(2), # interaction with reproduction status
                                     gammaSi = rnorm(1),
                                     gammaRi = rnorm(2),
                                    # sig.cycle=runif(1), #year
                                     tau.pair = structure(.Data= c(1, 0, 0, 1), .Dim=c(2, 2)), #runif(1),
                                     tau.id = structure(.Data= c(1.2, 0, 0, 1.2), .Dim=c(2, 2)),
                                     tau.yr = structure(.Data= c(0.9, 0, 0, 0.9), .Dim=c(2, 2)),
                                    .RNG.name=sample(c("base::Super-Duper","base::Wichmann-Hill"),1), .RNG.seed=sample(1:100,1))}

##### chain parameters #####
n.c.jags <- no_of_cores # ICS set number of chains = number of cores
n.a.jags <- 5000 # number of adaptive iterations
n.b.jags <- 20000 # number of burnin
n.t.jags <- 100 # thining interval
n.s.jags <- 1000  # number of samples (i.e. total number of iterations is nb.samples * thinning interval)


##### model #####
jagsIODnoRS <- run.jags("DHGLM_BBAL_IODnoRRSexEnvtInt.R", BUGS6.parameters, BUGS6.data, n.chains=n.c.jags, inits = BUGS6.inits.jags, 
                        sample = n.s.jags, thin = n.t.jags, burnin = n.b.jags, adapt = n.a.jags, method="parallel")


if(exists("output.file") == TRUE) {
    save(jagsIODnoRS, file = output_file)
}


  
