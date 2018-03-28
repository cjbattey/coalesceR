source("functions.R")

sim <- coalesce(nsamples = 10,nsites=16000,sims=10000,cores=8,mu=2e-8)

pdf(width=3,height=3,"sfs_ridge_plot_color.pdf")
#plot coalescent trees with edges colored by the number of descendent lineages
plot_coal_tree(sim,n_col=3,n=9,edge_colors = T,mar=1,cex=1,bw=F,brewerpal="YlOrRd")

#get mutation generations
muts <- do.call(rbind.data.frame,sim$mutations) 
muts <- subset(muts,ndesc<=10)
med_mut_times <- ddply(muts,.(ndesc),summarize,med_gen=median(gen))

#ridge plot for time distributions of mutations in each SFS class
ridges <- ggplot(data=muts,aes(x=gen,y=factor(ndesc),fill=factor(ndesc)))+
  #scale_fill_manual(values=grey.colors(max(muts$ndesc),start=0.1),guide=F)+
  scale_fill_manual(values=rev(brewer.pal(max(muts$ndesc),"YlOrRd")),guide=F)+
  ylab("Derived Allele Count")+xlab("Mutation Time\n (generations before present)")+
  theme_bw()+theme(panel.grid=element_blank(),
                   axis.text = element_text(size=8),
                   axis.title = element_text(size=8))+
  scale_x_continuous(breaks=c(0,6000,12000),limits=c(0,12000))+
  geom_density_ridges(lwd=0.25)+
  geom_segment(data=med_mut_times,aes(x=med_gen,y=ndesc,xend=med_gen,yend=ndesc-0.25),col="black",lwd=0.65)

vp.right <- grid::viewport(just="right",width=0.65,x = 0.65)
print(ridges,vp=vp.right)
dev.off()


ggplot(data=muts,aes(x=gen,y=ndesc))+
  ylab("SFS Class")+xlab("Generation")+
  theme_bw()+theme(panel.grid = element_blank())+
  scale_y_continuous(breaks=1:29)+
  xlim(0,10000)+
  geom_point(position=position_jitter(height=0.3),alpha=0.03,size=0.5)+
  geom_segment(data=mean_mut_times,aes(x=mean_gen,y=ndesc-0.3,xend=mean_gen,yend=ndesc+0.3),col="red",lwd=1)



