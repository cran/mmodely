
get.phylo.stats <- function(phylo, data, trait.clmn, gs.clmn='gn_sp', ace.method='REML',ace.scaled=TRUE,ace.kappa=1){
       
        rownames(data) <- data[,gs.clmn] 
        data <- data[!is.na(data[,trait.clmn]),]
        phylo <- drop.tip(phylo, phylo$tip.label[!phylo$tip.label %in% data$gn_sp])

        trait.v <- with(data[phylo$tip.label,], round(nv(get(trait.clmn),gn_sp),2)); #112
        names(trait.v) <- phylo$tip.label 

        ancest.label <- names(which.max(branching.times(phylo)));ancest.label  
        .get.ace <- function(a,pa=ancest.label){ a$ace[as.character(as.numeric(pa)+1)]}

        a <- ace(x=trait.v, phy=phylo, method=ace.method,scaled=ace.scaled, kappa=ace.kappa)
        #print(.get.ace(a))
        
        # uses v 0.0-8 phytools fucntion 'phylosig()' in phylosig.R 
        sig.lambda <- .phylosig(tree=phylo, method='lambda',x=trait.v,test=TRUE);print(sig.lambda)
        sig.K      <- .phylosig(tree=phylo, method='K',     x=trait.v,test=TRUE);print(sig.K)
        return(list(lambda=sig.lambda, K=sig.K, ace=a))
}



plot.transformed.phylo <- function(x, delta=1,kappa=1,...){
  dd <- VCV.array(x) ^ delta    # x is a tree object
  ddh <- as.phylo(hclust(dist(dd)))
  ddhh <- ddh
  ddhh$edge.length <- ddh$edge.length^kappa
  plot(ddhh,...)
}   


trim.phylo <- function(phylo, gs.vect){ 
  if(!ape::is.binary.tree(phylo)){
   warning('not a fully dicotomous tree, running the recommended "multi2di" to fix')
   phylo <- multi2di(phylo)
  } 
 phylo <- drop.tip(phylo, phylo$tip.label[!phylo$tip.label %in% gs.vect])
 return(phylo)
}
