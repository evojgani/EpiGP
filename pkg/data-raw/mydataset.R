## code to prepare `mydataset` dataset goes here

usethis::use_data("mydataset")
file.create("R/data.R")

library(BGLR)
data(wheat)
geno <- wheat.X
t1 <- sample(1:ncol(geno), 20)
t2 <- sample(1:ncol(geno), 20)
y1 <- rowSums((geno[,t1]==2) * (geno[,t2]==2))
t1 <- sample(1:ncol(geno), 20)
t2 <- sample(1:ncol(geno), 20)
y2 <- rowSums((geno[,t1]==2) * (geno[,t2]==0))
t1 <- sample(1:ncol(geno), 20)
t2 <- sample(1:ncol(geno), 20)
y3 <- rowSums((geno[,t1]==0) * (geno[,t2]==2))
t1 <- sample(1:ncol(geno), 20)
t2 <- sample(1:ncol(geno), 20)
y4 <- rowSums((geno[,t1]==0) * (geno[,t2]==0))
y <- y1+y2+y3+y4
Phenotype <- scale(y)
names(Phenotype) <- names(wheat.Y[,1])
save(Phenotype, file=paste('data-raw/Phenotype_wheat', '.Rdata', sep = ''))

usethis::use_data(Phenotype_wheat, overwrite = TRUE)


install.packages("sinew")
devtools::install_github("mdlincoln/docthis")

