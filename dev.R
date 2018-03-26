source("functions.R")
library(stringr);library(plyr)

sim <- coalesce(nsamples = 8,n=10000,mu=2e-7,nsites = 1000,max_gen = 1e8,sims =1,cores=8)

e <- sim$tables[[1]]
e$n_desc_1 <- str_count(e$l1,"\\(")+1
e$n_desc_2 <- str_count(e$l2,"\\(")+1
sfs <- c(e$mut1,e$mut2,e$n_desc_1,e$n_desc_2) %>% 
  matrix(ncol=2,byrow=F) %>% 
  data.frame() %>% 
  ddply(.(X2),summarize,count=sum(X1))
names(sfs) <- c("n_individuals","n_sites")
return(sfs)