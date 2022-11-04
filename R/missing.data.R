


compare.data.gs.vs.tree.tips <- function(data, phylo, match.on=c('gn_sp','rownames')[1]){
 if(match.on=='rownames'){
  gnsp.vect <- rownames(data)
 }else{
  gnsp.vect  <-  data$gn_sp
 }
 mismatch.TF <- !gnsp.vect %in% phylo$tip.label
 mismatches  <-  gnsp.vect[mismatch.TF]

 mismatch.count <- sum(mismatch.TF)
 if(mismatch.count==0){
  message('no mismatches found!')
 }else{
  message(paste('found',mismatch.count,'mismatches between tree and data')) 
  if(any(is.na(mismatches)))
   message('NAs in mismatch list indicate the alias table may need updating')
  return(mismatches)
 }
}


gs.names.mismatch.check <- function(df, alias.table.path, gs.clmn='gn_sp'){
 mismatches <- nv(df[which(rownames(df) != df$gn_sp),],gs.clmn)
 alias.tab <- read.delim(alias.table.path)
 if(length(mismatches) > nrow(alias.tab)){
  warning("more mismatches than rows in alias table!")
  non.alias.mm <- mismatches[!names(mismatches)%in% unlist(alias.tab)]
  print(paste('found',length(non.alias.mm),'(non alias) mismatches'))
  print(non.alias.mm)
 }
 print(paste(length(mismatches),'total mismatches'))
 invisible(mismatches)
}

missing.data <- function(x, cols=NULL, rows=NULL){ 
 msng.by.col <- apply(x, 2, function(x) sum(is.na(x)))
 msng.by.row <- apply(x, 1, function(x) sum(is.na(x)))
 message('NA counts by columns:')
 print(sort(msng.by.col[msng.by.col!=0], decreasing=TRUE))
 message('NA counts by rows:')
 print(sort(msng.by.row[msng.by.row!=0], decreasing=TRUE))
 ret <- list()
 if(!is.null(cols))
  ret[['by.col']] <- sapply(cols, function(v) { print(v);rownames(x)[which(is.na(x[,v]))]})
 if(!is.null(rows))
  ret[['by.row']] <- sapply(rows, function(r) { print(r);colnames(x)[which(is.na(x[r,]))]})
 if(!is.null(rows) | !is.null(cols)){
  message('Missingness by specific rows and/or columns:')
  return(ret)
 }
}

missing.fill.in <- function(x, var.from, var.to){
 missing.var.to <- is.na(x[,var.to])
 na.ct <- sum(missing.var.to)
 if(na.ct==0) stop("found zero missing values in 'var.to'")
 x[missing.var.to,var.to] <- x[missing.var.to,var.from]
 return(x)
}

drop.na.data <- function(df, vars=names(df)){ 
 has.missing <- apply(df[,vars], 1, function(x) any(is.na(x)))
 df[!has.missing,]
}


## accessory funciton to help with pgls "comparative.data funciton which requires non-missing data
interpolate <- function(df, taxa=c('genus','family'), clmns=1:length(df)){
   clmn.nms <- names(df)[clmns]
   for(tx in taxa){ 
    ints <- groupBy(df, by=tx, clmns=clmn.nms, aggregation=rep('mean',length(clmn.nms)),na.rm=TRUE) 
    for(gn in unique(df[,tx])){ 
     for(clmn.nm in clmn.nms){
     df[df[[tx]]==gn & is.na(df[, clmn.nm]), clmn.nm] <- ints[gn, clmn.nm] 
   }}}
   return(df)
}

