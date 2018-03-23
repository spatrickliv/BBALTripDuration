library(runjags)
library(coda)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 3) {
    cat("usage is:\n");
    cat("   <script_name> input_file output_file no_of_cores\n\n");
    q()
}

input_file  <- args[1]
output_file <- args[2]
no_of_cores <- as.integer(args[3])

if(file.access(input_file, mode=0)) {
   cat("input file is not readable\n\n")
}

data=read.table(input_file, header=TRUE)

###### if on computer ########

setwd("C:\\Users\\spatrick\\Desktop\\LivAntarctica Changed\\Stats\\jags var")
data<-read.table("bugsdata2.txt")
library(runjags)
library(coda)
no_of_cores <- 3

#################################

I <-is.na(data$ids)|
  is.na(data$Sexe)|
  is.na(data$Cycle)|
  is.na(data$CycleBand)|
  is.na(data$Rsnum)|
  is.na(data$RSplus1num)|
  is.na(data$CyclePairBand)|
  is.na(data$Stade)|
  is.na(data$ntrajet)|
  is.na(data$Tripdur)|
  is.na(data$Age)|
  is.na(data$IODAnnual)| 
  is.na(data$SOIAnnual)|  
  is.na(data$AAOAnnual)|  
  is.na(data$ONIAnnual)
  
data<-data[!I,]

## remove preponte

birdID<-as.data.frame(data$ids)
birdIDnum<-as.data.frame(as.numeric(droplevels(data$ids)))
unique(birdIDnum)
write.table(cbind(birdID,birdIDnum), "BirdIDs.txt")
write.table(data, "croppedbugsdata2.txt")

data$ntrajet<- as.numeric(scale(data$ntrajet))
data$Tripdur<- as.numeric(scale(sqrt(data$Tripdur)))
data$Age<- as.numeric(scale(data$Age))
data$IODAnnual<- as.numeric(scale(data$IODAnnual))
data$SOIAnnual<- as.numeric(scale(data$SOIAnnual))
data$AAOAnnual<- as.numeric(scale(data$AAOAnnual))
data$ONIAnnual<- as.numeric(scale(data$ONIAnnual))

db <- data.frame(Cycle=sort(unique(data$Cycle)),Cycle.n=1:length(unique(data$Cycle)))
data$Cycle.n <- db$Cycle.n[match(data$Cycle,db$Cycle)]
  
BUGS6.data <- list(n=nrow(data), #121
                   n.Pairs = length(unique(data$CyclePairBand)),  #n.pair
                   n.Cycles = length(unique(data$Cycle)), #n.cycle
                   n.ids=length(unique(data$ids)),
                   y = data$Tripdur,
                   ###### Sexe is coded as female 1 and male 2 now to code for the interaction more easily
                   Sexe = as.numeric(data$Sexe), ##female coded as 1 and Males as 2
                   Rsnum = data$Rsnum,
                   IODAnnual = data$IODAnnual,
                   SOIAnnual = data$SOIAnnual,
                   AAOAnnual = data$AAOAnnual,
                   ###### Same thing reproductive status is coded differently. instead of having one variable for each status now, we have only 1 variable varying between 1-3
                   Repro = as.numeric(data$Stade) , ##"BROODING"=1, "INCUBATION"=2, "PRE-PONTE"=3,
                   Brooding =ifelse(data$Stade=="BROODING",1,0), 
                   Preponte = ifelse(data$Stade=="PRE-PONTE",1,0),
                   #Incubation = ifelse(data$Stade=="INCUBATION",1,0),
                   id = as.numeric(droplevels(data$ids)),
                   CyclePairBand = as.numeric(droplevels(data$CyclePairBand)),
                   Cycle = data$Cycle.n, 
                   Age = data$Age,
                   mat.id = diag(2),
                   zero = c(0,0)
)

nbeta <- 5
ngamma <- 5
#you have to provide inital values for the priors since OpenBUGS can crash when reasonable initial values are not specied

BUGS6.parameters <- c("beta","gamma", "betaSi", "betaRi", "gammaSi", "gammaRi", "sig.pair","sig2.pair","tau.id","sig2.id","SexEnv","ReproEnv")#,"uid")

#using JAGS

###### multiple paramters added to be able to code for the interaction easily
BUGS6.inits.jags <- function() {list(beta = rnorm(nbeta), #fixed
                                     gamma = rnorm(ngamma), 
                                     betaSi = rnorm(2), # interaction with Sexe
                                     betaRi = rnorm(3), # interaction with reproduction status
                                     gammaSi = rnorm(2),
                                     gammaRi = rnorm(3),
                                    # sig.cycle=runif(1), #year
                                     sig.pair=runif(1), #pair
                                     tau.id = structure(.Data= c(1.2, 0, 0, 1.2), .Dim=c(2, 2)),
                                    .RNG.name=sample(c("base::Super-Duper","base::Wichmann-Hill"),1), .RNG.seed=sample(1:100,1))}

n.c.jags <- 3 # number of chains
n.c.jags <- no_of_cores # ICS set number of chains = number of cores
n.a.jags <- 1000 # number of adaptive iterations
n.b.jags <- 2000 # number of burnin
n.s.jags <- 1000 # number of samples (i.e. total number of iterations is nb.samples * thinning interval)
n.t.jags <- 10 # thining interval
  
BUGS.6.jagsIODnoRS <- run.jags("DHGLM_BBAL_07-11IODnoRRSexEnvtInt_contrasts.txt", BUGS6.parameters, BUGS6.data, n.chains=n.c.jags, inits = BUGS6.inits.jags, 
                        sample = n.s.jags, thin = n.t.jags, burnin = n.b.jags, adapt = n.a.jags, n.sims=no_of_cores, method="parallel")


##save(BUGS.6.jags, file=output_file)

  
