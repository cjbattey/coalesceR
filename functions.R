# discrete coalescent tree simulation 

#sim_coal_discrete() simulates coalescence of a set of k lineages within a randomly-mating population of size n. 
#Each generation, every lineage has a (k-1)/n probability of coalescing with another sampled lineage, which
#is approximated by running n_lineages binomial trials with prob=(k-1)/n. If any coalescence events occur in a 
#given generation, two lineages are randomly chosen to coalesce and the age of the new node is recorded. 
#The number of mutations along the branch is then estimated by running n_gen binomial trials with prob=mu. 

#Mechanistically, this simulator tracks two vectors: one listing the sampled lineages and one the corresponding node ages. 
#As we move back in the tree, coalescent events cause one of the coalescing lineages to be replaced by a 
#new node written in newick tree format while the other is dropped [for example, if lineages A and B coalesce in a generation
#the "A" in the lineages vector is replaced with "(A,B)" and "B" is dropped]. 
#The same operation is conducted for the corresponding node ages. When all lineages have coalesced, 
#the lineages vector should be a complete newick string of the resulting tree. Mutations, coalescence events,
#and branch lengths are also stored in a table. 

#simulates coalescence of a set of individuals sampled from a population and returns a 
#phylogenetic tree (as a phylo object - see ape manual for info), a vector of mutation 
#counts for each branch in the order of the branches in the phylo object, and a table
#summarizing each coalescent event. 

sim_coal_discrete <- function(nsamples=8,n=10000,max_gen=1e7,mu=2e-8,nsites=1000){
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
      mut2 <- sum(rbinom(bl2,1,mu))
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
  return(list(tree,coal_table))
}

coalesce <- function(nsamples=8,n=10000,max_gen=1e7,mu=2e-8,nsites=1000,sims=1,cores=1){
  require(ape);require(magrittr);require(ape);require(ggtree);require(plyr);require(foreach);require(doMC)
  registerDoMC(cores=cores)
  out <- foreach(i=1:sims) %dopar% sim_coal_discrete(nsamples = nsamples,
                                                     n = n,max_gen = max_gen,
                                                     mu = mu,nsites = nsites)
  trees <- lapply(out,function(e) e[[1]])
  tables <- lapply(out,function(e) e[[2]])
  out <- list(trees,tables)
  names(out) <- c("trees","tables")
  return(out)
}

get_sfs <- function(sim,plot=F,sum=T){
  require(stringr);require(plyr)
  sfs_list <- lapply(sim$tables,function(e){
    e$n_desc_1 <- str_count(e$l1,"\\(")+1
    e$n_desc_2 <- str_count(e$l2,"\\(")+1
    sfs <- c(e$mut1,e$mut2,e$n_desc_1,e$n_desc_2) %>% 
      matrix(ncol=2,byrow=F) %>% 
      data.frame() %>% 
      ddply(.(X2),summarize,count=sum(X1))
    names(sfs) <- c("n_individuals","n_sites")
    return(sfs)
  })
  if(sum==T){
    sfs <- do.call(rbind.data.frame,sfs_list)
  } else {
    sfs <- sfs_list()
  }
  if(plot==T){
    print(ggplot(data=sfs,aes(x=factor(n_individuals),y=n_sites))+
            theme_bw()+theme(panel.grid = element_blank())+
            xlab("derived allele count")+
            ylab("sites")+
            geom_bar(stat="identity"))
  }
  return(sfs)
}

plot_coal_tree <- function(sim,n=NULL,n_col=1){
  if(is.null(n)){
    n <- length(sim$trees)
  }
  par(mfrow=c(ceiling(min(n,length(sim$trees))/n_col),n_col),
      mar=rep(2,4))
  for(j in 1:min(n,length(sim$trees))){
    tree <- sim$trees[[j]]
    coal_table <- sim$tables[[j]]
    mutations <- c()
    for(i in tree$edge.length){
      if(i %in% coal_table$bl1){
        mutations <- append(mutations,subset(coal_table,bl1==i)$mut1)
      } else if(i %in% coal_table$bl2){
        mutations <- append(mutations,subset(coal_table,bl2==i)$mut2)  }
    }
    mutations[mutations==0] <- NA
    plot(tree,cex=0.65,xpd=F)
    axisPhylo(cex.axis=0.65)
    edgelabels(mutations,frame="none",cex=0.65,adj=c(0.5,-0.3))
  }
  par(mfrow=c(1,1))
}







