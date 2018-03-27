source("functions.R")
library(stringr);library(plyr);library(ggplot2);library(magrittr);library(ggridges)

#tracking generation of mutations for plots of sfs class frequency ~ generation

sim_coal_discrete_mut <- function(nsamples=8,n=10000,max_gen=1e7,mu=2e-8,nsites=1000){
  lineages <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))[1:nsamples]
  #lineages <- sample(1:n,nsamples)
  mu <- mu*nsites
  gen <- 0
  j <- 1
  node_ages <- rep(0,nsamples)
  node_tbl <- data.frame(id=1:nsamples,age=rep(0,nsamples))
  coal_table <- data.frame(gen=integer(),l1=character(),l2=character(),
                           bl1=integer(),bl2=integer(),mut1=integer(),
                           mut2=integer())
  mutations <- data.frame(gen=NA,ndesc=NA)[-1,]
  while(length(lineages)>1 & gen<max_gen){
    if(length(n)>1){
      n <- n[j]
    }
    coal <- rbinom(length(lineages),1,(length(lineages)-1)/n)           #each of k lineages has a (k-1)/n probability of coalescence per generation
    if(any(coal)){ 
      coal_pair <- sample(1:length(lineages),2)                         #choose two lineages to coalesce
      l1 <- lineages[coal_pair[1]]
      l2 <- lineages[coal_pair[2]]
      bl1 <- gen-node_ages[coal_pair[1]]                                #record branch lengths
      bl2 <- gen-node_ages[coal_pair[2]]
      mut1 <- sum(rbinom(bl1,1,mu))                                     #get n mutations for each branch by binomial sampling with prob=mu each generation
      if(mut1>0){
        mut_gens <-  sample(node_ages[coal_pair[1]]:gen,mut1)
        ndesc1 <- str_count(l1,"\\(")+1
        mutations <- rbind(mutations,data.frame(gen=mut_gens,ndesc=ndesc1))
      }
      mut2 <- sum(rbinom(bl2,1,mu))
      if(mut2>0){
        mut_gens <-  sample(node_ages[coal_pair[2]]:gen,mut2)
        ndesc2 <- str_count(l2,"\\(")+1
        mutations <- rbind(mutations,data.frame(gen=mut_gens,ndesc=ndesc1))
      }
      
      coal_table <- rbind(coal_table,data.frame(gen,l1,l2,bl1,bl2,mut1,mut2))
      
      lineages[coal_pair[1]] <- paste0("(",l1,":",                      #update lineage and node age vectors
                                       bl1,",",
                                       l2,":",    
                                       bl2,")")     
      lineages <- lineages[-coal_pair[2]]                                         
      node_ages[coal_pair[1]] <- gen
      node_ages <- node_ages[-coal_pair[2]] 
    }
    gen <- gen+1
    j <- j+1
  }
  newick <- paste0(lineages,";")
  tree <- read.tree(text=newick)
  return(list(tree,coal_table,mutations))
}

coalesce_mut <- function(nsamples=8,n=10000,max_gen=1e7,mu=2e-8,nsites=1000,sims=1,cores=1){
  require(ape);require(magrittr);require(ape);require(ggtree);require(plyr);require(foreach);require(doMC)
  registerDoMC(cores=cores)
  out <- foreach(i=1:sims) %dopar% sim_coal_discrete_mut(nsamples = nsamples,
                                                     n = n,max_gen = max_gen,
                                                     mu = mu,nsites = nsites)
  trees <- lapply(out,function(e) e[[1]])
  tables <- lapply(out,function(e) e[[2]])
  mutations <- lapply(out,function(e) e[[3]])
  out <- list(trees,tables,mutations)
  names(out) <- c("trees","tables","mutations")
  return(out)
}

sim <- coalesce_mut(nsamples = 30,nsites=16000,sims=1000,cores=8,mu=2e-7)

muts <- do.call(rbind.data.frame,sim$mutations) 
muts <- subset(muts,ndesc<=20)
ggplot(data=muts,aes(x=gen,y=factor(ndesc)))+
  ylab("SFS Class")+xlab("Generation")+
  theme_bw()+theme(panel.grid=element_blank())+
  xlim(0,5000)+geom_density_ridges()


