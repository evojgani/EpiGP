

usethis::use_data("mydataset_add")
file.create("R/data_add.R")

library(BGLR)
data(wheat)
geno <- wheat.X

haplo <- t(geno)
haplo <- haplo[,sort(rep(1:599,2))]

library(MoBPS)

t1 <- cbind(sample(1:nrow(haplo), 80), 1, sample(1:nrow(haplo), 80), 1,
            c(rep(1,20), rep(0,60)), 0, c(rep(0,20), rep(1,20), rep(0,40)), 0,0,0, c(rep(0,40), rep(1,20), rep(0,20)), 0, c(rep(0,60), rep(1,20)))
t2 <- cbind(sample(1:nrow(haplo), 80), 1, sample(1:nrow(haplo), 80), 1,
            c(rep(1,20), rep(0,60)), 0, c(rep(0,20), rep(1,20), rep(0,40)), 0,0,0, c(rep(0,40), rep(1,20), rep(0,20)), 0, c(rep(0,60), rep(1,20)))

# shuffle.cor is correlation between genetic component
# new.phenotype.correlation is residual correlation
# heritability to control heritablity
population <- creating.diploid(haplo, real.bv.mult = list(t1,t2),
                               shuffle.traits = TRUE,
                               shuffle.cor = matrix(c(1,0.9,
                                                      0.9,1), nrow=2))


population <- breeding.diploid(population, new.phenotype.correlation = matrix(c(1,0.9,
                                                                                0.9,1), nrow=2),
                               phenotyping.gen = 1, heritability = c(0.8,0.8))

y <- t(get.pheno(population, gen=1))
Phenotype_BV <- y
rownames(Phenotype_BV) <- names(wheat.Y[,1])

save(Phenotype_BV, file=paste('data-raw/Phenotype_Bivar_wheat', '.Rdata', sep = ''))

usethis::use_data(Phenotype_BV, overwrite = TRUE)



install.packages("sinew")
devtools::install_github("mdlincoln/docthis")

