library(microseq)
library(OrgMassSpecR)
library(data.table)
library(protViz)
library(Peptides)

new_2 <- data.frame(Protein=NA, Mass=NA, Score=NA)
new <- c()
## Empty vectors
proteins <- c()
mass <- c()
scores <- c()


# print(new[with(new, order(Protein)),])


proteome <- readFasta("il1403.fasta")
masses <- as.numeric(c("796.044986046", "818.0249"))
tolerance <- 1.0


for (x in 1:length(proteome$Sequence)){
  a <- Digest(proteome$Sequence[x], enzyme = "trypsin", missed = 1, IAA = FALSE, N15 = FALSE)
  for (y in masses){
    new <- append(between(a$mz1, y-tolerance, y+tolerance, incbounds=TRUE), new)
  }
  proteins[x] <- substr(proteome$Header[x], 1, gregexpr('OS=', proteome$Header[x])[[1]][1]-2)
  mass[x] <- paste(signif(mw(proteome$Sequence[x], monoisotopic = FALSE)/1000, digits = 3), "kDa")
  scores[x] <- sum(new, na.rm = TRUE)
  new <- c()
}

results <- cbind(mass, proteins, scores)