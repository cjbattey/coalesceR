# discrete coalescent tree simulation 

#sim_coal_discrete() simulates coalescence of a set of k lineages within a randomly-mating population of size n. 
#Each generation, every lineage has a (k-1)/n probability of coalescing with another sampled lineage, which
#is approximated by running n_lineages binomial trials with prob=(k-1)/n. If any coalescence events occur in a 
#given generation, two lineages are randomly chosen to coalesce and the age of the new node is recorded. 
#The number of mutations along the branch is then estimated by running n_gen binomial trials with prob=mu. 

#Mechanistically, this simulator amounts to a randomized reordering of two vectors: 
#one listing extant lineages, and one listing the corresponding node ages. 
#Each generation there is a (k-1)/n probability of a coalescence event (where k is the number of lineages
#and n is the population size). Whenever a coalescence occurs we randomly choose two indexes from the lineage 
#vector (bc mating is random), replace the first index slot with a new node written in newick
#tree format, and drop the second index. For example, if lineages A and B coalesce in a generation
#the "A" in the lineages vector is replaced with "(A,B)" and "B" is dropped. 
#The same operation is conducted for the corresponding node ages. When all lineages have coalesced, 
#the lineages vector should be a complete newick string of the resulting tree. Mutations, coalescence events,
#and branch lengths are also stored in a table. 

#TLDR:
#simulates coalescence of a set of individuals sampled from a population and returns a 
#phylogenetic tree (as a phylo object - see ape manual for info), a vector of mutation 
#counts for each branch in the order of the branches in the phylo object, and a table
#summarizing each coalescent event. 
sim_coal_discrete <- function(nsamples=8,n=10000,max_gen=1e7,mu=2e-8,nsites=1000,track_mut_gen=T){
  require(ape);require(stringr)
  lineages <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))[1:nsamples]
  #lineages <- sample(1:n,nsamples)
  if(length(n)>1 & length(n)<max_gen){
    n <- append(n,rep(n[length(n)],max_gen-length(n))) #if n vector is shorter than max gen, repeat the last value to avoid infinite loop
  }
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
    if(length(n)>1){                                                    #update pop size if different across generations
      n <- n[j]
    }
    coal <- rbinom(length(lineages),1,(length(lineages)-1)/n)           #each of k lineages has a (k-1)/n probability of coalescence per generation
    if(any(coal)){ 
      coal_pair <- sample(1:length(lineages),2)                         #choose two lineages to coalesce
      l1 <- lineages[coal_pair[1]]
      l2 <- lineages[coal_pair[2]]
      bl1 <- gen-node_ages[coal_pair[1]]                                #record branch lengths
      bl2 <- gen-node_ages[coal_pair[2]]
      ndesc1 <- str_count(l1,"\\(")+1
      ndesc2 <- str_count(l2,"\\(")+1
      mut1 <- sum(rbinom(bl1,1,mu))                                     #get n mutations for each branch by binomial sampling with prob=mu each generation
      
      if(track_mut_gen==T){                                             #optionally track the generation of mutations in different SFS classes
        if(mut1>0){
          mut_gens <-  sample(node_ages[coal_pair[1]]:gen,mut1)
          mutations <- rbind(mutations,data.frame(gen=mut_gens,ndesc=ndesc1))
        }
        mut2 <- sum(rbinom(bl2,1,mu))
        if(mut2>0){
          mut_gens <-  sample(node_ages[coal_pair[2]]:gen,mut2)
          mutations <- rbind(mutations,data.frame(gen=mut_gens,ndesc=ndesc1))
        }
      }
      
      coal_table <- rbind(coal_table,data.frame(gen,l1,l2,bl1,bl2,mut1,mut2,ndesc1,ndesc2))
      
      lineages[coal_pair[1]] <- paste0("(",l1,":",                      #update lineage and node age vectors
                                       bl1,",",
                                       l2,":",    
                                       bl2,")")     
      lineages <- lineages[-coal_pair[2]]                                         
      node_ages[coal_pair[1]] <- gen
      node_ages <- node_ages[-coal_pair[2]] 
    } #end updates for one coalescence
    gen <- gen+1
    j <- j+1
  } #end coalescence loop
  newick <- paste0(lineages,";")
  tree <- read.tree(text=newick)
  if(track_mut_gen==T){
    return(list(tree,coal_table,mutations))
  } else {
    return(list(tree,coal_table))
  }
}

#run multiple simulations in parallel
coalesce <- function(nsamples=8,n=10000,max_gen=1e7,mu=2e-8,nsites=1000,sims=1,cores=1){
  require(foreach);require(doMC)
  registerDoMC(cores=cores)
  out <- foreach(i=1:sims) %dopar% sim_coal_discrete(nsamples = nsamples,
                                                     n = n,max_gen = max_gen,
                                                     mu = mu,nsites = nsites)
  trees <- lapply(out,function(e) e[[1]])
  tables <- lapply(out,function(e) e[[2]])
  mutations <- lapply(out,function(e) e[[3]])
  out <- list(trees,tables,mutations)
  names(out) <- c("trees","tables","mutations")
  return(out)
}

#plot a sample of coalescent trees generated by coalesce()
plot_coal_tree <- function(sim,n=NULL,n_col=1,axis=F,tip_labels=F,mar=0.5,
                           cex=0.75,edge_colors=T,pal="viridis",edge_width=1.75){
  require(RColorBrewer);require(ggplot2);require(ape);require(plyr);require(viridis)
  if(is.null(n)){
    n <- length(sim$trees)
  }
  par(mfrow=c(ceiling(min(n,length(sim$trees))/n_col),n_col),
      mar=rep(mar,4))
  for(j in 1:min(n,length(sim$trees))){
    tree <- sim$trees[[j]]
    coal_table <- sim$tables[[j]]
    
    sfs_class <- c() #build vector of n_descendents in order of branches in phylo object
    for(i in tree$edge.length){
      if(i %in% coal_table$bl1){
        sfs_class <- append(sfs_class,subset(coal_table,bl1==i)$ndesc1)
      } else if(i %in% coal_table$bl2){
        sfs_class <- append(sfs_class,subset(coal_table,bl2==i)$ndesc2)  }
    }
    
    mutations <- c() #build vector of mutations in order of branches in phylo object
    for(i in tree$edge.length){
      if(i %in% coal_table$bl1){
        mutations <- append(mutations,subset(coal_table,bl1==i)$mut1)
      } else if(i %in% coal_table$bl2){
        mutations <- append(mutations,subset(coal_table,bl2==i)$mut2)  }
    }
    mutations[mutations==0] <- NA
    
    if(edge_colors==T){
      if(pal %in% c("greys","bw")){
        col <- mapvalues(sfs_class,from=1:(length(tree$tip.label)-1),
                         to=grey.colors(length(tree$tip.label)-1,start=0))
        plot(tree,xpd=F,edge.color = col,show.tip.label = tip_labels,edge.width = edge_width)
      } else if(pal %in% c("viridis","magma","inferno","plasma","cividis")){
        col <- mapvalues(sfs_class,from=1:(length(tree$tip.label)-1),
                         to=viridis(length(tree$tip.label)-1,option = pal))
        plot(tree,xpd=F,edge.color = col,show.tip.label = tip_labels,edge.width = edge_width)
      } else if(pal %in% rownames(brewer.pal.info)){
        col <- mapvalues(sfs_class,from=1:(length(tree$tip.label)-1),
                         to=rev(brewer.pal(length(tree$tip.label)-1,pal)))
        plot(tree,xpd=F,edge.color = col,show.tip.label = tip_labels,edge.width = edge_width)
      }
    } else if(edge_colors==F){
      mutations[mutations==0] <- NA
      plot(tree,cex=cex,xpd=F,edge.width=edge_width,show.tip.label=tip_labels)
      edgelabels(mutations,frame="none",cex=cex,adj=c(0.5,-0.3))
    }
    if(axis==T){
      axisPhylo(cex.axis=0.65)
    }
  }
  par(mfrow=c(1,1))
}

#helper function for reformatting mutation/descendent data 
getsfsdata <- function(i){
  matrix(c(i$ndesc1,i$ndesc2,i$mut1,i$mut2),ncol=2,byrow=F)
}

#returns site-frequency spectrum of coalescent simulations generated by coalesce()
get_sfs <- function(sim,plot=F,sum=T,cores=8){
  require(foreach);require(ggplot2)
  registerDoMC(cores=cores)
  sfs <- foreach(i=sim$tables,.combine = rbind) %dopar% getsfsdata(i)
  sfs <- tapply(sfs[,2],sfs[,1],function(e) sum(e))
  sfs <- data.frame(n_ind=1:9,n_sites=sfs)
  if(plot==T){
    print(ggplot(data=sfs,aes(x=factor(n_ind),y=n_sites))+
            theme_bw()+theme(panel.grid = element_blank())+
            xlab("derived allele count")+
            ylab("sites")+
            geom_bar(stat="identity"))
  }
  return(sfs)
}





