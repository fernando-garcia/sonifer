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
  tmpF[,3]
}



factorial(7) # 7! = 5040 -> number of different permutations of variable Y

misdatos <- read.csv("data/datos.originales.txt",head=F,row.names=1) # misdatos : 1611 rows, 7 columns


### MIC from permutated data

nper <- 1000  # number of permutations to be done

MICper <- matrix(NA, nrow = nrow(misdatos)*(nrow(misdatos)-1)/2, ncol = nper)

for (i in 1:nper) {

 # permutation of original data (no need of permuting row 1)
 miper <- t(apply(misdatos[-1,], 1, sample))

 # now we compute MIC between each row of original data
 # and the rest of permutated rows
 fila <- 1

 for (j in 1:(nrow(misdatos)-1)) {

   for (k in j:nrow(miper)) {

     myMIC <- getMIC(misdatos[j,], miper[k,])
       
       
     # El myMIC esta patatero, lo tienes que arreglar para que lo calcule.
     # rMINE te coge el objeto R (no se si como yo lo he escrito)
     # pero no se si el resultado de lo devuelve en txt o en objeto R
     # Como devuelve MIC y mas cosas tienes que coger el que quieras
     # para que myMIC sea un solo valor

     MICper[fila,i] <- myMIC

     fila <- fila + 1

   }
 }

}




### True MIC

MICbuenos <- readMINEResults("MINE.output.csv") # FG: Esta funcion devuelve un vector con los MICs del fichero csv
# hay que tomar solo el valor MIC
# tienes que asegurarte de que las comparaciones de parejas de genes estan
# en el mismo orden en MICper y en MICbuenos
# si estan en el mismo orden, me lo dices y modificamos codigo



### p-values

mispvalores <- NULL

for (i in 1:length(MICbuenos)) {

 elpval <- length(which(MICper[i,] > MICbuenos[i])) / nper

 mispvalores <- c(mispvalores, elpval)

}


### adjusted p-values
pval.ajustado <- p.adjust(mispvalores)
