library(edgeR)

rm(list = ls())
setwd("/home/antonio/Escritorio/Riñoncete")
condicion <- factor(c(rep("kich", 66), rep("kirc", 529), rep("kirp", 289)))

# Cargamos los datasets con los subtipos de la enfermedad

load("KICH.rda")
load("KIRC.rda")
load("KIRP.rda")

# Construimos un único dataset con todos los subtipos y borramos los demás
# para ahorrar espacio

dfinal <- data.frame(kichfinal, kircfinal, kirpfinal)
rm(kichfinal)
rm(kirpfinal)
rm(kircfinal)
dfinal <- dfinal[-1,]

# Comprobación de que todos los datos estén sean del tipo "numeric"

indx <- sapply(dfinal, is.factor)
dfinal[indx] <- lapply(dfinal[indx], function(x) as.numeric(as.character(x)))

# Análisis de expresión diferencial

dge <- DGEList(dfinal, group = condicion)
keep <- rowSums(cpm(dge)>1) >= 15 # Filtro para aquellos genes que tengan al menos cpm >= 1 en la mitad de las muestras
dge <- dge[keep,]
dge$samples$lib.size <- colSums(dge$counts) # Recalculo el tamaño de libreria

dge.n <- calcNormFactors(dge)

design <- model.matrix( ~ condicion)
colnames(design) <- c("kich", "kirc", "kirp")


design
dge.c <- estimateCommonDisp(dge.n)
dge.t <- estimateTagwiseDisp(dge.c)

fit <- glmQLFit(dge.t, design)
qlf <- glmQLFTest(fit, coef = 2:3) # ANOVA edgeR

topTags(qlf, n = 100)

genes <- topTags(qlf, n = 100)
genes <- rownames(genes$table)

df <- dfinal[genes,]


save(df, file = "df100genesnormal.rda")
write.csv(df, file = "subtype_riñonkichkirckirp.csv", quote = F)