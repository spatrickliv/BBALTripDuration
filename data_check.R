library(lme4)

data.old <- read.table("bugsdata2.txt")
data.new <- read.table("bugsdata5.txt")


###### data processing ######
### removing NAs
data.proc <- function(data){
    I <-is.na(data$ids)|
        is.na(data$Sexe)|
        is.na(data$Cycle)|
        is.na(data$CycleBand)|
        ##is.na(data$Rsnum)|
        ##is.na(data$RSplus1num)|
        ##  is.na(data$CyclePairBand)|
        is.na(data$Stade)|
        ##is.na(data$ntrajet)|
        is.na(data$Tripdur)|
        is.na(data$Age)|
        is.na(data$IODAnnual) 
        ##is.na(data$SOIAnnual)|  
        ##is.na(data$AAOAnnual) ##|  
##  is.na(data$ONIAnnual)
    data <- data[!I,]
    data <-subset(data, Tripdur < 600)

    data$sTripdur <- as.numeric(scale(sqrt(data$Tripdur)))
    data$Age<- as.numeric(scale(data$Age))
    data$IODAnnual<- as.numeric(scale(data$IODAnnual))

    db <- data.frame(Cycle=sort(unique(data$Cycle)),Cycle.n=1:length(unique(data$Cycle)))
    data$Cycle.n <- db$Cycle.n[match(data$Cycle,db$Cycle)]

    data
}

data.old <- data.proc(data.old)
data.new <- data.proc(data.new)

m.old <- lmer(sTripdur ~ Age + (Sexe + Stade) * IODAnnual + (1|ids) + (1|Cycle), data.old)
m.new <- update(m.old, data = data.new)

summary(m.old, correlation = FALSE)
summary(m.new, correlation = FALSE)

sum.data <- function(data){
    summ <- matrix(c(nrow(data),
                     length(unique(data$ids)),
                     sum(data$Tripdur<7),
                     max(data$Tripdur),
                     sum(data$Tripdur>1000)
                     ),
                   ncol=1)
              
    rownames(summ) <- c("N",
                     "N_individuals",
                     "N_trip_<7",
                     "max_tripdur",
                     "N_trip_>1000"
                     )
    print(summ, digits=2)
}

    
sum.data(data.old)
sum.data(data.new)


