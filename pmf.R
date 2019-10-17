masses <- as.numeric(c("4478", "1600", "1489", "1344")) # sample masses
proteome <- readFasta("il1403.fasta") # sample fasta
tolerance <- 1 # sample tolerance, typically: 0.1

pmf <- function(proteome, tolerance, masses){
# Proteome: uniprot fasta. 
# Tolerance: number (+- Da).
# Masses: as.numeric(c("1", "2", "3"))

library(OrgMassSpecR)
library(microseq)
library(data.table)

new <- c()
finds <- c()

for (x in 1:length(proteome$Sequence)){
  a <- Digest(proteome$Sequence[x], enzyme = "trypsin", missed = 0, IAA = TRUE, N15 = FALSE)
  for (y in masses){
     new <- append(between(a$mz1, y-tolerance, y+tolerance, incbounds=TRUE), new)
  }
  if (sum(new, na.rm = TRUE) >= length(masses)){
    finds <- append(proteome$Header[x], finds)
  }
  new <- c()
}
return(finds)
}