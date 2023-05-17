# the formula for calculating AIC via the "small sample size" correction 

get.mod.outcome <- function(model) {as.character(as.formula(model))[2]}
get.mod.vars <- function(model) {strsplit(x=as.character(as.formula(model))[3],split=' \\+ ')[[1]]}
get.mod.clmns <- function(model, gs.clmn='gn_sp') { c(get.mod.outcome(model),get.mod.vars(model),gs.clmn) }

# simple model varibale count function (works by breaking a formula into it's three parts (via as.character) then counting the split on "+") 
count.mod.vars <- function(model) if(inherits(model,'formula')) length(get.mod.vars(model)) else NA


## create a list of all possible models based on combinations of all candiate variables
get.model.combos <- function(outcome.var, predictor.vars, min.q=1){
  all.models <- paste(predictor.vars, collapse="+")
  for(mi in (length(predictor.vars)-1):min.q){ 
    new.combos <-  t(combn(predictor.vars,  m=mi))
    new.models <- apply(new.combos, 1, paste, collapse="+")
    all.models <- c(all.models, new.models)
  }
 ## prepend the new composite "loud.calls" variable to the beginning
 full.models <- paste(outcome.var,all.models,sep='~')
}


ct.possible.models  <- function(q){
 sum(unlist(sapply(1:q,function(i) ncol(combn(q,i)))))
}
