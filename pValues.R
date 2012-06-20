# OJO!!!! Pon el MINE.jar en el directorio donde tengas esto que si no te fallara. 
# Tambien puedes cambiar una funcion del MINE.r para indicar donde lo tienes

install.packages("rJava") # 1-time initialization step
source("MINE.r")



getMIC <- function(orig,perm, output="tmp"){
  rMINE(rbind(as.numeric(orig),as.numeric(perm)),output,style="all.pairs") # hago el as.numeric porque orig era un data.frame y protestaba el MINE
  tmpF <- read.csv(paste(output,"allpairs,cv=0.0,B=n^0.6,Results.csv",sep=","))
  tmpF[1,3]
}

readMINEResults <- function(file){
  tmpF <- read.csv(file)
  mynames <- t(apply(tmpF[,1:2], 1, sort))
  mynames <- apply(mynames, 1, paste, collapse = ";")
  tmpF <- tmpF[,3]
  names(tmpF) <- mynames
  tmpF
}



factorial(7) # 7! = 5040 -> number of different permutations of variable Y

misdatos <- read.csv("data/datos.originales.txt",head=F,row.names=1) # misdatos : 1611 rows, 7 columns




### MIC from permutated data

nper <- 1000  # number of permutations to be done

MICper <- matrix(NA, nrow = nrow(misdatos)*(nrow(misdatos)-1)/2, ncol = nper)


# generating gene comparisons for MIC and using them as rownames of MICper
mygenes <- rownames(misdatos)

mycombi <- combn(length(mygenes), 2)

mycompar <- apply(mycombi, 2, function (x) { paste(sort(mygenes[x]), collapse = ";") }) 

rownames(MICper) <- mycompar



# computing MICper
for (i in 1:nper) {

 # permutation of original data (no need of permuting row 1)
 miper <- t(apply(misdatos, 1, sample))

 # now we compute MIC between each row of original data
 # and the rest of permutated rows

 for (j in 1:length(mycompar)) {

   myMIC <- getMIC(misdatos[mycombi[1,j],], miper[mycombi[2,j],])

   MICper[j,i] <- myMIC
 }
}




### True MIC

MICbuenos <- readMINEResults("data/MINE.output.csv") # FG: Esta funcion devuelve un vector con los MICs del fichero csv


# same gene/comparison order in MICper and MICbuenos
MICper <- MICper[names(MICbuenos),]




### p-values

mispvalores <- NULL

for (i in 1:length(MICbuenos)) {

 elpval <- length(which(MICper[i,] > MICbuenos[i])) / nper

 mispvalores <- c(mispvalores, elpval)

}



### adjusted p-values
pval.ajustado <- p.adjust(mispvalores)



### SUMMARY

mysummary <- cbind(mycompar, MICbuenos, mispvalores, pval.ajustado)
colnames(mysummary) <- c("comparison", "MIC", "pvalue", "adjusted_pvalue")

write.table(mysummary, "results.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
