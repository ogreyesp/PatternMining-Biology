library(arules)

labels <- c("dlpN","dlpS", "dmN", "dmS", "extraprostN", "extraprostS", "gleason7N","gleason7S","gleason8N","gleason8S",
            "invperinueralN","invperinueralS","mtxaltacarN","mtxN","mtxS","mtxoseaN","mtxoseaS",
            "metforminaN","metforminaS","resistenciaN","resistenciaS","respuesta7mesesN","respuesta7mesesS")

path <- "/home/oscar/workspace/IMIBIC/datasets/cancer/"

for(label in labels){
  
  dataset <-
    read.csv(paste(path,"clinico-",label,"-norm-disc.csv",sep = ""))
  
  #filter by target variable
  positiveTarget = 'S'
  negativeTarget = 'N'
  
  numberAtt = ncol(dataset) - 1
  targetIndex= ncol(dataset) 
  
  datasetP <- dataset[dataset[targetIndex] == positiveTarget, 1:numberAtt]
  
  datasetN <- dataset[dataset[targetIndex] == negativeTarget, 1:numberAtt]
  
  patterns <- apriori(datasetP,
                      parameter = list(
                        support = 1,
                        target = "maximally frequent itemsets",
                        maxlen = numberAtt
                      ))
  
  # Arules uses its own objects
  patterns <- as(patterns, "data.frame")
  
  if(nrow(patterns) > 0){
  
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
    newData[order(-newData$growthRate,-newData$supportP), ]
  
  write.table(
    newData,
    file = paste(path,"patterns-growthRate-",label,".csv",sep = ""),
    quote = TRUE,
    sep = "," ,
    row.names = FALSE,
    col.names = TRUE)
  }
  }