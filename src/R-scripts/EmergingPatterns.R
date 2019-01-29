library(arules)

dataset <-
  read.csv("/home/oscar/workspace/IMIBIC/datasets/cancer/multi-label-norm-disc.csv")

#filter by label
positiveTarget = 'S'
negativeTarget = 'N'

numLabels = 13

numberAtt = ncol(dataset) - numLabels

labelNames <- colnames(dataset[(numberAtt+1):(numberAtt+13)])

for (label in labelNames) {
  
  datasetP <- dataset[dataset[label] == positiveTarget, 1:numberAtt]
  
  datasetN <- dataset[dataset[label] == negativeTarget, 1:numberAtt]
  
  patterns <- apriori(datasetP,
                      parameter = list(
                        support = 0.9,
                        target = "maximally frequent itemsets",
                        maxlen = numberAtt
                      ))
  
  # Arules uses its own objects
  patterns <- as(patterns, "data.frame")
  
  newData <- data.frame(
      "itemset" = character(0),
      "supportP" = numeric(0),
      "supportN" = numeric(0),
      "growthRate" = numeric(0))
  
  for (t in 1:nrow(patterns)) {
    
    transaction <- patterns[t, "items"]
    
    #conver the factor in a character vector
    transaction <- as.character(transaction)
    
    itemset <- unlist(strsplit(transaction, ","))
    
    supportN <- 0
    
    #for each row in datasetN
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
    newData[order(-newData$growthRate,-newData$supportP), ]
  
  write.table(
    newData,
    file = paste("/home/oscar/workspace/IMIBIC/datasets/cancer/patterns-growthRate-",label,".csv",sep = ""),
    quote = TRUE,
    sep = "," ,
    row.names = FALSE,
    col.names = TRUE)
}