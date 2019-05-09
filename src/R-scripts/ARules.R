rm(list = ls())
library(arules)

dataset <- read.csv("/home/antonio/Escritorio/CARySD/CAR/BRCA/discretizaen5all.csv")

numAtt <- ncol(dataset) - 1 

dataset$Class <- as.factor(dataset$Class)

rules <- apriori(dataset,
                 parameter = list(support = 0.16,
                                  confidence = 0.5,
                                  target = "rules",
                                  minlen = 3,
                                  maxlen = numAtt),
                 appearance = list(rhs = c("Class=T"),
                                   default="lhs"))

# Arules uses its own objects
df <- as(rules, "data.frame")

# Sort by lift and confidence
df <- df[order(-df$support, df$lift, df$confidence), ]

write.table(df, file="/home/antonio/Escritorio/CARySD/CAR/BRCA/rules-test.csv", 
            quote = TRUE, sep="," , row.names = FALSE, col.names = TRUE)
