library(rsubgroup)
library(discretization)

path <- "/home/antonio/Escritorio/CARySD/SD/BRCA"

setwd(path)

dataset <- read.csv("discretizaen5all.csv")

# creating a target binary variable
target <- as.target("Class", "T")

#creating a config object for passing to subgroup discovery algorithms
namesAtt <- colnames(dataset) 
namesAtt <- namesAtt[-length(namesAtt)]
config <- new("SDTaskConfig", attributes= namesAtt, qf = "wracc", minsize = 33, k = 500)
task <- CreateSDTask(source= dataset, target = target, config = config)

patterns <-  DiscoverSubgroupsByTask(task, as.df= TRUE)
head(patterns)
filtrado <- strsplit(as.character(patterns$description), ",")

# Filtramos por aquellas reglas que tengan al menos 2 items

mequedo <- c()
for(i in 1:length(filtrado)){
  if(length(filtrado[[i]]) >= 2){
    mequedo <- c(mequedo, i)
  }
}
patterns <- patterns[mequedo,]

# Calculamos el límite superior e inferior para calcular el wraccn

ub <- 0.5*0.5
lb <- -0.25

# Calculamos la métrica para cada una de las muestras

wraccn <- c()
for(i in 1:nrow(patterns)){
  x = (patterns$quality[i] - (lb)) / (ub - lb)
  wraccn = c(wraccn, x)
}
patterns <- cbind(wraccn, patterns)

write.csv(x = patterns, file = "sd.csv", row.names = FALSE)
