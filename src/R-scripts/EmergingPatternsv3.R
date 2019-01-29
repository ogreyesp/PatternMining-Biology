library(arules)

path <- "/home/oscar/Dropbox/Publicaciones/Paper patter mining-datos biologicos/1-Diabetes/"

datasetT <- read.csv(paste(path,"T-withoutMissingValues-norm-disc.csv",sep = ""))

datasetN <- read.csv(paste(path,"N-withoutMissingValues-norm-disc.csv",sep = ""))

numberAtt <- ncol(datasetT)

patterns <- apriori(datasetT,
                    parameter = list(
                      support = 0.95,
                      target = "frequent itemsets",
                      maxlen = numberAtt
                    ))

# Arules uses its own objects
patterns <- as(patterns, "data.frame")

if(nrow(patterns) > 0){
  
  newData <- data.frame(
    "itemset" = character(0),
    "supportT" = numeric(0),
    "supportN" = numeric(0),
    "growthRate" = numeric(0))
  
  for (t in 1:nrow(patterns)) {
    
    transaction <- patterns[t, "items"]
    
    #convert the factor in a character vector
    transaction <- as.character(transaction)
    
    itemset <- unlist(strsplit(transaction, ","))
    
    supportN <- 0
    
    #for each row in datasetN. In this part is computed the support of the pattern in the other dataset
    for (f in 1:nrow(datasetN)) {
      
      flag = TRUE
      
      for (item in itemset) {
        
        itemValues <- unlist(strsplit(item, "="))
        
        itemName <- itemValues[1]
        
        if (startsWith(itemName, "{")) {
          itemName <- unlist(strsplit(itemName, '\\{'))[2]
        }
        
        itemValue <- itemValues[2]
        
        if (endsWith(itemValue, "}")) {
          itemValue <- unlist(strsplit(itemValue, '\\}'))[1]
        }
        
        currentValue <- as.character(datasetN[f, itemName])
        
        if (currentValue != itemValue)
        {
          flag <- FALSE
          break()
        }
      }
      
      # The pattern matchs
      if (flag) {
        supportN <- supportN + 1
      }
    }
    
    supportN <- supportN / nrow(datasetN)
    supportP <- patterns[t, "support"]
    growthRate <- supportP / supportN
    
    #register the frequent itemset
    newData <-
      rbind(
        newData,
        data.frame(
          "itemset" = transaction,
          "supportP" = supportP,
          "supportN" = supportN,
          "growthRate" = growthRate
        )
      )
  }
  
  # Sort by growthRate
  newData <-
    newData[order(-newData$growthRate,-newData$supportP),]
  
  write.table(
    newData,
    file = paste(path,"patterns-growthRate.csv",sep = ""),
    quote = TRUE,
    sep = "," ,
    row.names = FALSE,
    col.names = TRUE)
}